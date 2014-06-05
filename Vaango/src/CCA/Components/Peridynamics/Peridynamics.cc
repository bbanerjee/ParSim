#include <CCA/Components/Peridynamics/Peridynamics.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Components/Peridynamics/PeridynamicsDomainBoundCond.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/LinearInterpolator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
#include <Core/Util/DebugStream.h>
#include <Core/Thread/Mutex.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

using namespace Vaango;

using Uintah::UintahParallelComponent;
using Uintah::SchedulerP;
using Uintah::ProcessorGroup;
using Uintah::Task;
using Uintah::DataWarehouse;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::Ghost;
using Uintah::MaterialSet;
using Uintah::MaterialSubset;
using Uintah::ParticleSubset;
using Uintah::ProblemSpecP;

using Uintah::ParticleVariable;
using Uintah::constParticleVariable;
using Uintah::NCVariable;
using Uintah::constNCVariable;
using Uintah::NodeIterator;
using Uintah::long64;
using Uintah::Matrix3;
using Uintah::NeighborList;
using SCIRun::Point;
using SCIRun::Vector;
using SCIRun::IntVector;

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern SCIRun::Mutex cerrLock;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDDoing:+,PDDebug:+".....
//  bash     : export SCI_DEBUG="PDDoing:+,PDDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDDoing", false);
static DebugStream cout_dbg("PDDebug", false);

/*! Construct */
Peridynamics::Peridynamics(const ProcessorGroup* myworld) :
  UintahParallelComponent(myworld)
{
  d_labels = scinew PeridynamicsLabel();
  d_flags = scinew PeridynamicsFlags(myworld);
  d_interpolator = scinew Uintah::LinearInterpolator();

  d_dataArchiver = 0;
  d_familyComputer = 0;
  d_defGradComputer = 0;
  d_bondIntForceComputer = 0;
  d_intForceComputer = 0;

  d_numGhostNodes = 1;
  d_recompile = false;
}

/*! Delete */
Peridynamics::~Peridynamics()
{
  delete d_labels;
  delete d_flags;
  delete d_interpolator;

  if (d_familyComputer) {
    delete d_familyComputer;
  }

  if (d_defGradComputer) {
    delete d_defGradComputer;
  }

  if (d_bondIntForceComputer) {
    delete d_bondIntForceComputer;
  }

  if (d_intForceComputer) {
    delete d_intForceComputer;
  }
}

/*! Read input file and set up problem */
void 
Peridynamics::problemSetup(const ProblemSpecP& prob_spec, 
                           const ProblemSpecP& restart_prob_spec,
                           Uintah::GridP& grid,
                           Uintah::SimulationStateP& sharedState)
{
  cout_doing << "Doing problemSetup: Peridynamics " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  d_sharedState = sharedState;
  dynamic_cast<Uintah::Scheduler*>(getPort("scheduler"))->setPositionVar(d_labels->pPositionLabel);
  
  d_dataArchiver = dynamic_cast<Uintah::Output*>(getPort("output"));
  if(!d_dataArchiver){
    throw Uintah::InternalError("Peridynamics:couldn't get output port", __FILE__, __LINE__);
  }

  // Set up a pointer to the problem spec that will also work with restarts
  ProblemSpecP restart_mat_ps = prob_spec;
  if (restart_prob_spec) restart_mat_ps = restart_prob_spec;

  // Read all Peridynamics <SimulationFlags> d_flags (look in PeridynamicsFlags.cc)
  // Also set up <PhysicalConstants>
  d_flags->readPeridynamicsFlags(restart_mat_ps, d_dataArchiver);
  d_sharedState->setParticleGhostLayer(Ghost::AroundNodes, d_flags->d_numCellsInHorizon);

  // Creates Peridynamics material w/ constitutive models and damage models
  // Looks for <MaterialProperties> and then <Peridynamics>
  materialProblemSetup(restart_mat_ps);

  // Create the family computer
  d_familyComputer = scinew FamilyComputer(d_flags, d_labels);

  // Create the deformation gradient computer object
  d_defGradComputer = scinew PeridynamicsDefGradComputer(d_flags, d_labels);

  // Create the internal force computer objects
  d_bondIntForceComputer = scinew BondInternalForceComputer(d_flags, d_labels);
  d_intForceComputer = scinew ParticleInternalForceComputer(d_flags, d_labels);

  // Set up contact model (TODO: Create Peridynamics version)
  //d_contactModel = 
  //  Uintah::ContactFactory::create(UintahParallelComponent::d_myworld, restart_mat_ps, sharedState, 
  //                         d_labels, d_flags);
}

void 
Peridynamics::materialProblemSetup(const ProblemSpecP& prob_spec) 
{
  cout_doing << "Doing material problem set up: Peridynamics " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  //Search for the MaterialProperties block and then get the Peridynamics section
  ProblemSpecP mat_ps = prob_spec->findBlockWithOutAttribute("MaterialProperties");
  ProblemSpecP peridynamics_mat_ps = mat_ps->findBlock("Peridynamics");

  for (ProblemSpecP ps = peridynamics_mat_ps->findBlock("material"); ps != 0;
       ps = ps->findNextBlock("material") ) {

    std::string index("");
    ps->getAttribute("index",index);
    std::stringstream id(index);
    const int DEFAULT_VALUE = -1;
    int index_val = DEFAULT_VALUE;

    id >> index_val;

    if( !id ) {
      // stringstream parsing failed... on many (most) systems, the
      // original value assigned to index_val would be left
      // intact... but on some systems (redstorm) it inserts garbage,
      // so we have to manually restore the value.
      index_val = DEFAULT_VALUE;
    }
    // cout << "Material attribute = " << index_val << ", " << index << ", " << id << "\n";

    //Create and register as an Peridynamics material
    PeridynamicsMaterial *mat = scinew PeridynamicsMaterial(ps, d_sharedState, d_flags);

    // When doing restart, we need to make sure that we load the materials
    // in the same order that they were initially created.  Restarts will
    // ALWAYS have an index number as in <material index = "0">.
    // Index_val = -1 means that we don't register the material by its 
    // index number.
    if (index_val > -1){
      d_sharedState->registerPeridynamicsMaterial(mat, index_val);
    }
    else{
      d_sharedState->registerPeridynamicsMaterial(mat);
    }
  }
}

void 
Peridynamics::outputProblemSpec(ProblemSpecP& root_ps)
{
  cout_doing << "Doing output problem spec: Peridynamics" 
             << __FILE__ << ":" << __LINE__ << std::endl;

  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP d_flags_ps = root->appendChild("Peridynamics");
  d_flags->outputProblemSpec(d_flags_ps);

  ProblemSpecP mat_ps = 0;
  mat_ps = root->findBlockWithOutAttribute("MaterialProperties");

  if (mat_ps == 0)
    mat_ps = root->appendChild("MaterialProperties");
    
  ProblemSpecP peridynamic_ps = mat_ps->appendChild("Peridynamics");
  for (int i = 0; i < d_sharedState->getNumPeridynamicsMatls();i++) {
    PeridynamicsMaterial* mat = d_sharedState->getPeridynamicsMaterial(i);
    ProblemSpecP cm_ps = mat->outputProblemSpec(peridynamic_ps);
  }

  d_contactModel->outputProblemSpec(peridynamic_ps);
}

/*----------------------------------------------------------------------------------------
 *  Method:  scheduleInitialize
 *  The initialization task is actually composed of two tasks
 *    - actuallyInitialize (for initializing data for which no previous information
 *                          is required)
 *    - findNeighborsInHorizon (finding the list of neighbor particles inside the
 *                              horizon)
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleInitialize(const Uintah::LevelP& level,
                                 SchedulerP& sched)
{
  cout_doing << "Doing schedule initialize: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  //--------------------------------------------------------------------------
  // Task 1: actuallyInitialize
  //--------------------------------------------------------------------------
  Task* t = scinew Task("Peridynamics::actuallyInitialize",
                        this, &Peridynamics::actuallyInitialize);

  const PatchSet* patches = level->eachPatch();
  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_labels->partCountLabel);
  t->computes(d_labels->pPositionLabel);
  t->computes(d_labels->pParticleIDLabel);
  t->computes(d_labels->pHorizonLabel);
  t->computes(d_labels->pDisplacementLabel);
  t->computes(d_labels->pStressLabel);

  t->computes(d_labels->pMassLabel);
  t->computes(d_labels->pVolumeLabel);
  t->computes(d_labels->pSizeLabel);
  t->computes(d_labels->pVelocityLabel);
  t->computes(d_labels->pExternalForceLabel);
  t->computes(d_sharedState->get_delt_label(),level.get_rep());
  t->computes(d_labels->pCellNAPIDLabel,zeroth_matl);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for(int body = 0; body < numBodies; body++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

    // Add deformation gradient initial computes
    d_defGradComputer->addInitialComputesAndRequires(t, peridynamic_matl, patches);

    // Add constitutive model computes
    PeridynamicsMaterialModel* cm = peridynamic_matl->getMaterialModel();
    cm->addInitialComputesAndRequires(t, peridynamic_matl, patches);

    // Add damage model computes
    PeridynamicsDamageModel* fm = peridynamic_matl->getDamageModel();
    fm->addInitialComputesAndRequires(t, peridynamic_matl, patches);
  }

  sched->addTask(t, patches, d_sharedState->allPeridynamicsMaterials());

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) delete zeroth_matl; // shouln't happen, but...

  //--------------------------------------------------------------------------
  // Task 2: findNeighborsInHorizon
  //--------------------------------------------------------------------------
  // After the basic quantities have been initialized, find the neighbor information
  // for the particles
  cout_doing << "Doing schedule compute neighbors/family: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* neighborFinderTask = scinew Task("Peridynamics::findNeighborsInHorizon",
                                      this, &Peridynamics::findNeighborsInHorizon);

  for(int body = 0; body < numBodies; body++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

    // Add family computer requires and computes
    d_familyComputer->addInitialComputesAndRequires(neighborFinderTask, peridynamic_matl, patches);
  }

  sched->addTask(neighborFinderTask, patches, d_sharedState->allPeridynamicsMaterials());
}

/*----------------------------------------------------------------------------------------
 *  Method:  actuallyInitialize
 *  Purpose: initializing data for which no previous information is required
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::actuallyInitialize(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls,
                                 DataWarehouse*,
                                 DataWarehouse* new_dw)
{
  cout_doing << "Doing actually initialize: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;
  particleIndex totalParticles=0;
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    Uintah::CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_labels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    for(int body=0 ; body < matls->size(); body++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      //int indx = peridynamic_matl->getDWIndex();

      // Create particles
      particleIndex numParticles = peridynamic_matl->countParticles(patch);
      totalParticles+=numParticles;
      peridynamic_matl->createParticles(numParticles, cellNAPID, patch, new_dw);

      // Initialize deformation gradient and shape
      d_defGradComputer->initialize(patch, peridynamic_matl, new_dw);

      // Initialize constitutive model
      peridynamic_matl->getMaterialModel()->initialize(patch, peridynamic_matl, new_dw);

      // Initialize damage model
      peridynamic_matl->getDamageModel()->initialize(patch, peridynamic_matl, new_dw); 
    }
  }

  // Initialize some per patch variables
  new_dw->put(Uintah::sumlong_vartype(totalParticles), d_labels->partCountLabel);
}

/*----------------------------------------------------------------------------------------
 *  Method:  findNeighborsInHorizon
 *  Purpose: Finding the list of neighbor particles inside the horizon of a particle
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::findNeighborsInHorizon(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* matls,
                                     DataWarehouse*,
                                     DataWarehouse* new_dw)
{
  cout_doing << "Doing find neighbors in horizon: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    for(int body=0; body < matls->size(); body++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

      // Create neighbor list
      d_familyComputer->createNeighborList(peridynamic_matl, patch, new_dw);

    }
  }
}

/*------------------------------------------------------------------------------------------------
 *  Method:  restartInitialize
 *  Purpose: Set variables that are normally set during the initialization
 *           phase, but get wiped clean when you restart
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::restartInitialize()
{
  cout_doing << "Doing restart initialize: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeStableTimestep
 *  Purpose: A stable timestep is compute for explicit (forward Euler) time integration 
 *           based on the CFL condition and a timestep multiplier
 *  Note:    The timestep size depends on the bulk wave speed in the material and
 *           is therefore actually computed in the MaterialModel
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeStableTimestep(const Uintah::LevelP& level,
                                            SchedulerP& sched)
{
  cout_doing << "Doing schedule compute time step: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  // However, this task needs to do something in the case that Peridynamics
  // is being run on more than one level. (**NOT IMPLEMENTED** BB 28 May 2014)
  Task* t = 0;
  t = scinew   Task("Peridynamics::actuallyComputeStableTimestep",
                    this, &Peridynamics::actuallyComputeStableTimestep);

  const MaterialSet* peridynamic_matls = d_sharedState->allPeridynamicsMaterials();
  t->computes(d_sharedState->get_delt_label(),level.get_rep());
  sched->addTask(t,level->eachPatch(), peridynamic_matls);
}

/*----------------------------------------------------------------------------------------
 *  Method:  actuallyComputeStableTimestep
 *  Purpose: Does nothing here unless there are multiple levels (NOT IMPLEMENTED)
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::actuallyComputeStableTimestep(const ProcessorGroup*,
                                            const PatchSubset* patches,
                                            const MaterialSubset* ,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw)
{
  cout_doing << "Doing actually compute stable time step: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Put something here to satisfy the need for a reduction operation in
  // the case that there are multiple levels present
  const Uintah::Level* level = getLevel(patches);
  new_dw->put(Uintah::delt_vartype(999.0), d_labels->delTLabel, level);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleTimeAdvance
 *  Purpose: These are the tasks that are performed sequentially in each timestep 
 *           of the simulation
 * ------------------------------------------------------------------------------------------------
 */
void
Peridynamics::scheduleTimeAdvance(const Uintah::LevelP & level,
                                  SchedulerP & sched)
{
  cout_doing << "Doing schedule time advance: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  const PatchSet* patches = level->eachPatch();
  const MaterialSet* matls = d_sharedState->allPeridynamicsMaterials();

  // Apply any external loads that may have been prescribed
  scheduleApplyExternalLoads(sched, patches, matls);

  // Interpolate data from particles to the grid so that peridynamics can use
  // the MPM algorithm for contact detection
  scheduleInterpolateParticlesToGrid(sched, patches, matls);

  // Apply any loads that may have been generated due to contact
  //scheduleApplyContactLoads(sched, patches, matls);

  // Compute the peridynamics deformation gradient
  // ** NOTE ** the accuracy of the gradient depends on the density of particles
  //            used to represent a volume of material
  scheduleComputeDeformationGradient(sched, patches, matls);

  // Compute the Cauchy stress and PK1 stress using continuum mechanics models
  scheduleComputeStressTensor(sched, patches, matls);

  // Compute the internal force for the bonds and then the sum of the volume
  // weighted bond internal forces
  scheduleComputeInternalForce(sched, patches, matls);

  // Compute the accleration and integrate to find velocity and displacement
  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);

  // Project the acceleration and velocity to the grid
  scheduleProjectParticleAccelerationToGrid(sched, patches, matls);

  // Correct the contact forces
  //scheduleCorrectContactLoads(sched, patches, matls);

  // Set the grid boundary conditions
  scheduleSetGridBoundaryConditions(sched, patches, matls);

  // Update the particle velocities and displacements
  scheduleUpdateParticleKinematics(sched, patches, matls);

  // Compute the damage tensor
  scheduleComputeDamage(sched, patches, matls);

  // Finalize the particle state for this time step
  scheduleFinalizeParticleState(sched, patches, matls);

  // Now that everything has been computed create a task that
  // takes the updated information and relocates particles if needed
  sched->scheduleParticleRelocation(level, 
                                    d_labels->pPositionLabel_preReloc,
                                    d_sharedState->d_particleState_preReloc,
                                    d_labels->pPositionLabel,
                                    d_sharedState->d_particleState,
                                    d_labels->pParticleIDLabel,
                                    matls, 1);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleApplyExternalLoads
 *  Purpose: This tasks sets up the external load application to each particle
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleApplyExternalLoads(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  cout_doing << "Doing schedule apply external loads: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::applyExternalLoads",
                        this, &Peridynamics::applyExternalLoads);
                  
  t->requires(Task::OldDW, d_labels->pPositionLabel,          Ghost::None);

  t->computes(d_labels->pExternalForceLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  applyExternalLoads
 *  Purpose: Actually apply the external loads to the particle in a consistent manner.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::applyExternalLoads(const ProcessorGroup* ,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  cout_doing << "Doing **NOTHING TODO** apply external loads: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Loop thru patches to update external force vector
  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);
    
    int numBodies=d_sharedState->getNumPeridynamicsMatls();
    
    for(int body = 0; body < numBodies; body++){

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

      int matlIndex = peridynamic_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      // Get the particle data
      constParticleVariable<Point>   pPosition;
      old_dw->get(pPosition, d_labels->pPositionLabel, pset);

      // Allocate space for the external force vector
      ParticleVariable<Vector>       pExternalForce_new;
      new_dw->allocateAndPut(pExternalForce_new, d_labels->pExternalForceLabel_preReloc,  pset);

      for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        pExternalForce_new[*iter] = 0.;
      }
    } // matl loop
  }  // patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleInterpolateParticlesToGrid
 *  Purpose: This task used the grdi shape functions to inetrpolate particle data
 *           to the grid.  MPM uses this extensively and we need this to do MPM-like
 *           contact in peridynamics.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  cout_doing << "Doing schedule interpolate particle to grid: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew   Task("Peridynamics::interpolateParticlesToGrid",
                          this, &Peridynamics::interpolateParticlesToGrid);

  Ghost::GhostType aroundNodes = Ghost::AroundNodes;
  int numGhostNodes = 1;
  t->requires(Task::OldDW, d_labels->pPositionLabel, aroundNodes, numGhostNodes);
  t->requires(Task::OldDW, d_labels->pMassLabel, aroundNodes, numGhostNodes);
  t->requires(Task::OldDW, d_labels->pVolumeLabel, aroundNodes, numGhostNodes);
  t->requires(Task::OldDW, d_labels->pSizeLabel, aroundNodes, numGhostNodes);
  t->requires(Task::OldDW, d_labels->pVelocityLabel, aroundNodes, numGhostNodes);
  t->requires(Task::NewDW, d_labels->pExternalForceLabel_preReloc, aroundNodes, numGhostNodes);
  t->requires(Task::OldDW, d_labels->pDefGradLabel, aroundNodes, numGhostNodes);

  t->computes(d_labels->gMassLabel);
  t->computes(d_labels->gVolumeLabel);
  t->computes(d_labels->gVelocityLabel);
  t->computes(d_labels->gExternalForceLabel);
  t->computes(d_labels->gMassLabel,     d_sharedState->getAllInOneMatl(), Task::OutOfDomain);
  t->computes(d_labels->gVolumeLabel,   d_sharedState->getAllInOneMatl(), Task::OutOfDomain);
  t->computes(d_labels->gVelocityLabel, d_sharedState->getAllInOneMatl(), Task::OutOfDomain);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  interpolateParticlesToGrid
 *  Purpose: Use the grid cell shape functions to take particle information to the
 *           background grid.  This is needed for contact detection a la MPM.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::interpolateParticlesToGrid(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset* ,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{
  cout_doing << "Doing interpolate particles to grid: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Matrix3 dummyMatrix(0.0);

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    // Create a copy of the LinearInterpolator
    Uintah::ParticleInterpolator* interpolator = d_interpolator->clone(patch);
    std::vector<IntVector> nodeIndices(interpolator->size());
    std::vector<double> shapeFunction(interpolator->size());

    // Allocate and initialize global data
    NCVariable<double> gMassGlobal;
    new_dw->allocateAndPut(gMassGlobal, d_labels->gMassLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gMassGlobal.initialize(std::numeric_limits<double>::epsilon());

    NCVariable<double> gVolGlobal;
    new_dw->allocateAndPut(gVolGlobal, d_labels->gVolumeLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gVolGlobal.initialize(std::numeric_limits<double>::epsilon());

    NCVariable<Vector> gVelGlobal;
    new_dw->allocateAndPut(gVelGlobal, d_labels->gVelocityLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gVelGlobal.initialize(Vector(0.0));

    // Loop through peridynamics objects
    int numBodies = d_sharedState->getNumPeridynamicsMatls();
    for(int body = 0; body < numBodies; body++){

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get the particle subset
      int numGhostNodes = 1;  // For LinearInterpolator
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch, 
                                                       Ghost::AroundNodes, numGhostNodes, 
                                                       d_labels->pPositionLabel);

      // Get the particle data needed for projection of the relevant quantities
      // to the grid
      constParticleVariable<Point>  pPosition;
      old_dw->get(pPosition, d_labels->pPositionLabel, pset);

      constParticleVariable<double> pMass; 
      old_dw->get(pMass, d_labels->pMassLabel, pset);

      constParticleVariable<double> pVolume;
      old_dw->get(pVolume, d_labels->pVolumeLabel, pset);

      constParticleVariable<Matrix3> pSize;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);

      constParticleVariable<Vector> pVelocity;
      old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);

      constParticleVariable<Vector> pExtForce;
      new_dw->get(pExtForce, d_labels->pExternalForceLabel_preReloc, pset);

      constParticleVariable<Matrix3> pDefGrad_old;
      old_dw->get(pDefGrad_old, d_labels->pDefGradLabel, pset);

      // Create and initialize arrays for the grid data that will be computed
      NCVariable<double> gMass;
      new_dw->allocateAndPut(gMass, d_labels->gMassLabel, matlIndex, patch);
      gMass.initialize(std::numeric_limits<double>::epsilon());

      NCVariable<double> gVolume;
      new_dw->allocateAndPut(gVolume, d_labels->gVolumeLabel, matlIndex, patch);
      gVolume.initialize(std::numeric_limits<double>::epsilon());

      NCVariable<Vector> gVelocity;
      new_dw->allocateAndPut(gVelocity, d_labels->gVelocityLabel, matlIndex, patch);
      gVelocity.initialize(Vector(0,0,0));

      NCVariable<Vector> gExtForce;
      new_dw->allocateAndPut(gExtForce, d_labels->gExternalForceLabel, matlIndex, patch);
      gExtForce.initialize(Vector(0,0,0));

      // loop over all particles in the patch and project particle data to the grid
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {

        particleIndex idx = *iter;
        interpolator->findCellAndWeights(pPosition[idx], nodeIndices, shapeFunction, dummyMatrix, dummyMatrix);

        // Compute particle momentum 
        // The momentum is projected to the grid and then divided by the grid mass to get the grid velocity
        Vector pMomentum = pVelocity[idx]*pMass[idx];

        // Add each particle's contribution to the local mass & velocity 
        // Must use the node indices
        // Iterates through the nodes which receive information from the current particle
        IntVector node;
        for (unsigned int k = 0; k < nodeIndices.size(); k++) {
          node = nodeIndices[k];
          if (patch->containsNode(node)) {
            gMass[node] += pMass[idx]*shapeFunction[k];
            gVelocity[node] += pMomentum*shapeFunction[k];
            gVolume[node] += pVolume[idx]*shapeFunction[k];
          }
        }

      } // End of particle loop

      for (NodeIterator iter=patch->getExtraNodeIterator(); !iter.done();iter++) {
        IntVector node = *iter; 
        gMassGlobal[node] += gMass[node];
        gVolGlobal[node]  += gVolume[node];
        gVelGlobal[node]  += gVelocity[node];
        gVelocity[node] /= gMass[node];
      }

    }  // End loop over materials

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done();iter++) {
      IntVector node = *iter;
      gVelGlobal[node] /= gMassGlobal[node];
    }

    // remove the interpolator clone
    delete interpolator;

  }  // End loop over patches
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleApplyContactLoads
 *  Purpose: This is a **TODO** task.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleApplyContactLoads(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  cout_doing << "Doing schedule apply contact loads: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;
  d_contactModel->addComputesAndRequiresInterpolated(sched, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeDeformationGradient
 *  Purpose: This task sets up the quantities requires for the deformation gradient 
 *           to be computed.  Needs neighbor information.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeDeformationGradient(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute deformation gradient: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::computeDeformationGradient",
                        this, &Peridynamics::computeDeformationGradient);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_defGradComputer->addComputesAndRequires(t, matl, patches);
  }

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeDeformationGradient
 *  Purpose: Compute the peridynamics deformation gradient and the inverse of the
 *           shape tensor.  Needs neighbor information.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeDeformationGradient(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset* ,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{
  cout_doing << "Doing compute deformation gradient: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int numBodies = d_sharedState->getNumPeridynamicsMatls();
    for (int body = 0; body < numBodies; body++) {
      PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
      d_defGradComputer->computeDeformationGradient(patch, matl, old_dw, new_dw);
    } // end matl loop
  } // end patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeStressTensor
 *  Purpose: This task sets up the quantities required for the Cauchy stress 
 *           to be computed.  
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeStressTensor(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute stress tensor: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::computeStressTensor",
                        this, &Peridynamics::computeStressTensor);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);

    // Add computes and requires specific to the material model
    PeridynamicsMaterialModel* cm = matl->getMaterialModel();
    cm->addComputesAndRequires(t, matl, patches);
  }

  t->computes(d_sharedState->get_delt_label(),getLevel(patches));

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeStressTensor
 *  Purpose: Use the continuum mechanics constitutive models to compute the
 *           Cauchy stress and the PK1 stress.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeStressTensor(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* ,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  cout_doing << "Doing compute stress tensor: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {

    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);

    PeridynamicsMaterialModel* cm = matl->getMaterialModel();
    cm->computeStressTensor(patches, matl, old_dw, new_dw);

  } // end matl loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeInternalForce
 *  Purpose: This method sets up two tasks:
 *           1) compute the bond internal forces for individual bonds
 *           2) a volume weighted sum of the bond internal forces to compute a particle
 *              internal force.task sets up the quantities required for the Cauchy stress 
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeInternalForce(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute internal force: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  //-------------------------------------------
  // Task 1: Compute bond internal forces
  //-------------------------------------------
  Task* task1 = scinew Task("Peridynamics::computeBondInternalForce",
                            this, &Peridynamics::computeBondInternalForce);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_bondIntForceComputer->addComputesAndRequires(task1, matl, patches);
  }

  sched->addTask(task1, patches, matls);
  
  //-------------------------------------------
  // Task 2: Compute particle internal forces
  //-------------------------------------------
  Task* task2 = scinew Task("Peridynamics::computeInternalForce",
                            this, &Peridynamics::computeInternalForce);

  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_intForceComputer->addComputesAndRequires(task2, matl, patches);
  }

  sched->addTask(task2, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeBondInternalForce
 *  Purpose: Compute the peridynamics bond internal force.
 *           Needs neighbor information.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeBondInternalForce(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* ,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  cout_doing << "Doing compute bond internal force: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_bondIntForceComputer->computeInternalForce(patches, matl, old_dw, new_dw);
  } // end matl loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeInternalForce
 *  Purpose: Compute the internal force at each particle.
 *           Needs neighbor information.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeInternalForce(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  cout_doing << "Doing compute particle internal force: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_intForceComputer->computeInternalForce(patches, matl, old_dw, new_dw);
  } // end matl loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeAndIntegrateAcceleration
 *  Purpose: This method sets up the required variables and the computed
 *           variables for solving the momentum equation using
 *           Forward Euler.  The acceleration and an intermediate 
 *           velocity are computed.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
                                                      const PatchSet* patches,
                                                      const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute and integrate acceleration: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::computeAndIntegrateAcceleration",
                        this, &Peridynamics::computeAndIntegrateAcceleration);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );
  t->requires(Task::OldDW, d_labels->pMassLabel,                   Ghost::None);
  t->requires(Task::OldDW, d_labels->pVolumeLabel,                 Ghost::None);
  //t->requires(Task::OldDW, d_labels->pDisplacementLabel,           Ghost::None);
  t->requires(Task::OldDW, d_labels->pVelocityLabel,               Ghost::None);
  t->requires(Task::NewDW, d_labels->pInternalForceLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_labels->pExternalForceLabel_preReloc, Ghost::None);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for(int body = 0; body < numBodies; body++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
    const MaterialSubset* matlset = peridynamic_matl->thisMaterial();

    //t->computes(d_labels->pDisplacementLabel_preReloc, matlset);
    t->computes(d_labels->pVelocityStarLabel, matlset);
    t->computes(d_labels->pAccelerationLabel, matlset);
  }
  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeAndIntegrateAcceleration
 *  Purpose: Use the Peridynamics balance of momentum equation to solve for the acceleration
 *           at each particle in the patch.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeAndIntegrateAcceleration(const ProcessorGroup*,
                                              const PatchSubset* patches,
                                              const MaterialSubset*,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw)
{
  cout_doing << "Doing compute and integrate acceleration: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Uintah::delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );
 
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    Vector gravity = d_flags->d_gravity;

    for(int body = 0; body < d_sharedState->getNumPeridynamicsMatls(); body++){

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get required variables for this patch
      //constParticleVariable<Vector> pDisp;
      constParticleVariable<Vector> pInternalForce, pExternalForce;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<double> pMass, pVolume;

      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      old_dw->get(pMass,          d_labels->pMassLabel,          pset);
      old_dw->get(pVolume,        d_labels->pVolumeLabel,        pset);
      //old_dw->get(pDisp,          d_labels->pDisplacementLabel,  pset);
      old_dw->get(pVelocity,      d_labels->pVelocityLabel,      pset);

      new_dw->get(pInternalForce, d_labels->pInternalForceLabel_preReloc, pset);
      new_dw->get(pExternalForce, d_labels->pExternalForceLabel_preReloc, pset);

      // Create variables for the results
      //ParticleVariable<Vector> pDisp_star;
      //new_dw->allocateAndPut(pDisp_star,     d_labels->pDisplacementLabel_preReloc, pset);
      ParticleVariable<Vector> pVelocity_star, pAcceleration;
      new_dw->allocateAndPut(pVelocity_star, d_labels->pVelocityStarLabel, pset);
      new_dw->allocateAndPut(pAcceleration,  d_labels->pAccelerationLabel, pset);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;

        if (pMass[idx] > 0.0) {

          // compute the specific volume (volume/mass = 1/density)
	  double specific_volume = pVolume[idx]/pMass[idx];

          // compute the acceleration due to the internal force
          Vector internal_force_acc = pInternalForce[idx]*specific_volume;

          // add the accleration due to the body force to that due to the internal force
          internal_force_acc += gravity;

          // compute the accleration due to the external force
          Vector external_force_acc = pExternalForce[idx]*specific_volume;

          // substract the external_force_acc from the internal_force_acc
          // to get the acceleration
          Vector acc = internal_force_acc - external_force_acc;
          pAcceleration[idx] = acc;

          // Integrate the acceleration to get the velocity
          // Forward Euler
          Vector vel = pVelocity[idx] + acc*delT;
          pVelocity_star[idx] = vel;

          // Integrate the velocity to get the displacement
          // Forward Euler
          // Vector disp = pDisp[idx] + vel*delT;
          // pDisp_star[idx] = disp;

        } else {
          pAcceleration[idx] = Vector(0.0);
          pVelocity_star[idx] = Vector(0.0);
          // pDisp_star[idx] = Vector(0.0);
        }
      } // end particle loop
    } // end matls loop
  } // end patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleProjectAccelerationToGrid
 *  Purpose: Sets up task for projecting the computed particle velocities and
 *           acclerations to the background grid
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleProjectParticleAccelerationToGrid(SchedulerP& sched,
                                                       const PatchSet* patches,
                                                       const MaterialSet* matls)

{
  cout_doing << "Doing schedule project particle acceleration to grid: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::projectParticleAccelerationToGrid",
                        this, &Peridynamics::projectParticleAccelerationToGrid);
                  
  int numGhostNodes = 1;  // Linear interpolation only
  t->requires(Task::OldDW, d_labels->pPositionLabel,     Ghost::AroundNodes, numGhostNodes);
  t->requires(Task::NewDW, d_labels->pVelocityStarLabel, Ghost::AroundNodes, numGhostNodes);
  t->requires(Task::NewDW, d_labels->pAccelerationLabel, Ghost::AroundNodes, numGhostNodes);

  t->computes(d_labels->gAccelerationLabel);
  t->computes(d_labels->gVelocityStarLabel);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  projectParticleAccelerationToGrid
 *  Purpose: Grid quantities are used to determine what happens when particles 
 *           reach a domain boundary.  This requires that the particle quantities
 *           be projected on to the grid.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::projectParticleAccelerationToGrid(const ProcessorGroup*,
                                                const PatchSubset* patches,
                                                const MaterialSubset* ,
                                                DataWarehouse* old_dw,
                                                DataWarehouse* new_dw)
{
  cout_doing << "Doing project particle acceleration to grid: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Matrix3 dummyMatrix(0.0);

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    int numBodies = d_sharedState->getNumPeridynamicsMatls();

    // Create a copy of the LinearInterpolator
    Uintah::ParticleInterpolator* interpolator = d_interpolator->clone(patch);
    std::vector<IntVector> nodeIndices(interpolator->size());
    std::vector<double> shapeFunction(interpolator->size());

    for(int body = 0; body < numBodies; body++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get the particle subset needed for projection to the grid
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch, Ghost::AroundNodes,
                                                       d_flags->d_numCellsInHorizon,
                                                       d_labels->pPositionLabel);

      // Get the required particle data
      constParticleVariable<Point> pPosition;
      old_dw->get(pPosition, d_labels->pPositionLabel, pset);
      
      constParticleVariable<Vector> pVelocityStar;
      new_dw->get(pVelocityStar, d_labels->pVelocityStarLabel, pset);
      
      constParticleVariable<Vector> pAcceleration;
      new_dw->get(pAcceleration, d_labels->pVelocityStarLabel, pset);

      // Allocate the computed grid variables
      NCVariable<Vector> gVelocity_star;
      new_dw->allocateAndPut(gVelocity_star, d_labels->gVelocityStarLabel,  matlIndex, patch);
      gVelocity_star.initialize(Vector(0.0, 0.0, 0.0));

      NCVariable<Vector> gAcceleration;
      new_dw->allocateAndPut(gAcceleration,  d_labels->gAccelerationLabel,  matlIndex, patch);
      gAcceleration.initialize(Vector(0.0, 0.0, 0.0));

      // Interpolate the particle velocity and acceleration to the grid
      // Loop through the particles in the set
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {

        particleIndex pidx = *iter;

        interpolator->findCellAndWeights(pPosition[pidx], nodeIndices, shapeFunction, 
                                         dummyMatrix, dummyMatrix);

        // Loop through the nodes
        IntVector node;
        for (unsigned int nidx = 0; nidx < nodeIndices.size(); nidx++) {

          node = nodeIndices[nidx]; 

          if (patch->containsNode(node)) {
            gVelocity_star[node] += pVelocityStar[pidx]*shapeFunction[nidx];
            gAcceleration[node] += pAcceleration[pidx]*shapeFunction[nidx];
          }

        } // end node loop

      } // end particle loop

    } // end matl loop

    // remove the interpolator clone
    delete interpolator;

  }  //end  patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleCorrectContactLoads
 *  Purpose: Use the contact model to correct any contact loads (grid based contact)
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleCorrectContactLoads(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  cout_doing << "Doing schedule correct contact loads: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  d_contactModel->addComputesAndRequiresIntegrated(sched, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleSetGridBoundaryConditions
 *  Purpose: Sets up task for applying grid (domain) boundary conditions for the
 *           simulation
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)

{
  cout_doing << "Doing schedule set grid boundary conditions: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::setGridBoundaryConditions",
                        this, &Peridynamics::setGridBoundaryConditions);
                  
  const MaterialSubset* all_materials = matls->getUnion();
  t->requires(Task::OldDW, d_sharedState->get_delt_label() );
  
  t->requires(Task::NewDW, d_labels->gVelocityLabel, Ghost::None);

  t->modifies(d_labels->gAccelerationLabel, all_materials);
  t->modifies(d_labels->gVelocityStarLabel, all_materials);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  setGridBoundaryConditions
 *  Purpose: Grid quantities are used to determine what happens when particles 
 *           reach a domain boundary.  This requires that the particle quantities
 *           be projected on to the grid.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::setGridBoundaryConditions(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* ,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  cout_doing << "Doing set grid boundary conditions: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Uintah::delt_vartype delT;            
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    int numBodies = d_sharedState->getNumPeridynamicsMatls();

    for(int body = 0; body < numBodies; body++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      int numGhostCells = 0;
      constNCVariable<Vector> gVelocity_old;  
      new_dw->get(gVelocity_old, d_labels->gVelocityLabel, matlIndex, patch, Ghost::None, numGhostCells);

      NCVariable<Vector> gVelocity_star, gAcceleration;
      new_dw->getModifiable(gVelocity_star, d_labels->gVelocityStarLabel,  matlIndex,patch);
      new_dw->getModifiable(gAcceleration,  d_labels->gAccelerationLabel,  matlIndex,patch);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      PeridynamicsDomainBoundCond bc;
      bc.setBoundaryCondition(patch, matlIndex, "Velocity", gVelocity_star); 
      bc.setBoundaryCondition(patch, matlIndex, "Symmetric", gVelocity_star);

      // Now back calculate acceleration as the difference between the velocity
      // interpolated to the grid (no bcs applied) and the new velocity_star
      for(NodeIterator iter=patch->getExtraNodeIterator();!iter.done(); iter++){
        IntVector c = *iter;
        gAcceleration[c] = (gVelocity_star[c] - gVelocity_old[c])/delT;
      } // node loop
    } // matl loop
  }  // patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleUpdateParticleKinematics
 *  Purpose: Set up the task for updating particle kinematics based on the updated
 *           grid velocities and accelerations.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleUpdateParticleKinematics(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)

{
  cout_doing << "Doing schedule update particle kinematics: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::updateParticleKinematics",
                        this, &Peridynamics::updateParticleKinematics);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );

  t->requires(Task::OldDW, d_labels->pPositionLabel,     Ghost::None);
  t->requires(Task::OldDW, d_labels->pVelocityLabel,     Ghost::None);
  t->requires(Task::OldDW, d_labels->pDisplacementLabel, Ghost::None);

  int numGhostCells = 1;  // Linear interpolation only
  t->requires(Task::NewDW, d_labels->gAccelerationLabel, Ghost::AroundCells, numGhostCells);
  t->requires(Task::NewDW, d_labels->gVelocityStarLabel, Ghost::AroundCells, numGhostCells);

  t->computes(d_labels->pPositionLabel_preReloc);
  t->computes(d_labels->pDisplacementLabel_preReloc);
  t->computes(d_labels->pVelocityLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  updateParticleKinematics
 *  Purpose: After the boundary velocities and accelerations have been corrected 
 *           on the grid by application of boundary conditions, update
 *           particle positions, velocities, and displacements
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::updateParticleKinematics(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* ,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  cout_doing << "Doing update particle state: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Matrix3 dummyMatrix(0.0);

  Uintah::delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for (int p=0; p<patches->size(); p++) {

    const Patch* patch = patches->get(p);

    // Create a copy of the LinearInterpolator
    Uintah::ParticleInterpolator* interpolator = d_interpolator->clone(patch);
    std::vector<IntVector> nodeIndices(interpolator->size());
    std::vector<double> shapeFunctions(interpolator->size());

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively
    int numBodies = d_sharedState->getNumPeridynamicsMatls();
    for(int body = 0; body < numBodies; body++) {

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get the particle set for this patch
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      // Get the arrays of particle values that are needed to update the kinematics
      constParticleVariable<Point> pPosition;
      old_dw->get(pPosition, d_labels->pPositionLabel, pset);

      constParticleVariable<Vector> pDisp;
      old_dw->get(pDisp, d_labels->pDisplacementLabel, pset);

      constParticleVariable<Vector> pVelocity;
      old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);

      // Get the arrays of grid data on which the new particle values depend
      int numGhostCells = 1;
      constNCVariable<Vector> gVelocity_star;
      new_dw->get(gVelocity_star, d_labels->gVelocityStarLabel, matlIndex, patch,
                  Ghost::AroundCells, numGhostCells);

      constNCVariable<Vector> gAcceleration;
      new_dw->get(gAcceleration, d_labels->gAccelerationLabel, matlIndex, patch,
                  Ghost::AroundCells, numGhostCells);

      // Create the arrays into which the updated results will be put
      ParticleVariable<Point> pPosition_new;
      new_dw->allocateAndPut(pPosition_new, d_labels->pPositionLabel_preReloc, pset);

      ParticleVariable<Vector> pDisp_new;
      new_dw->allocateAndPut(pDisp_new, d_labels->pDisplacementLabel_preReloc, pset);

      ParticleVariable<Vector> pVelocity_new;
      new_dw->allocateAndPut(pVelocity_new, d_labels->pVelocityLabel_preReloc, pset);

      // Loop over particles
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {

        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pPosition[idx], nodeIndices, shapeFunctions, dummyMatrix, dummyMatrix);

        // Accumulate the contribution from each surrounding vertex
        Vector vel(0.0,0.0,0.0);
        Vector acc(0.0,0.0,0.0);
        for (unsigned int kk = 0; kk < nodeIndices.size(); kk++) {
          IntVector node = nodeIndices[kk];
          vel += gVelocity_star[node]  * shapeFunctions[kk];
          acc += gAcceleration[node]   * shapeFunctions[kk];
        }

        // Update the particle's position and velocity
	pPosition_new[idx] = pPosition[idx] + vel*delT;
        pDisp_new[idx] = pDisp[idx] + vel*delT;
        pVelocity_new[idx] = pVelocity[idx] + acc*delT;

      } // end particle loop

    }  // end of matl loop

    // remove the interpolator clone
    delete interpolator;

  } // end of patch loop
  
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeDamage
 *  Purpose: Set up the task for computing the damage tensor based on the updated kinematics
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeDamage(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute damage: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::computeDamage",
                        this, &Peridynamics::computeDamage);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for(int body = 0; body < numBodies; body++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

    PeridynamicsDamageModel* d_damageModel = peridynamic_matl->getDamageModel();
    d_damageModel->addComputesAndRequires(t, peridynamic_matl, patches);
  }

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeDamage
 *  Purpose: Use the damage model to compute the damage tensor based on the updated kinematics
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeDamage(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* ,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  cout_doing << "Doing compute damage: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  for(int body = 0; body < d_sharedState->getNumPeridynamicsMatls(); body++){

    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);

    PeridynamicsDamageModel* damageModel = peridynamic_matl->getDamageModel();
    damageModel->computeDamageTensor(patches, peridynamic_matl, old_dw, new_dw);

  }
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleFinalizeParticleState
 *  Purpose: Set up the task for copying the data that will be used in the 
 *           next timestep.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleFinalizeParticleState(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls)

{
  cout_doing << "Doing schedule finalize particle state: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::finalizeParticleState",
                        this, &Peridynamics::finalizeParticleState);

  t->requires(Task::OldDW, d_labels->pParticleIDLabel,        Ghost::None);
  t->requires(Task::OldDW, d_labels->pMassLabel,              Ghost::None);
  t->requires(Task::OldDW, d_labels->pVolumeLabel,            Ghost::None); // TODO: pVolume should be 
                                                                            // updated in the material model
  t->requires(Task::OldDW, d_labels->pSizeLabel,              Ghost::None);

  t->requires(Task::OldDW, d_labels->pHorizonLabel,           Ghost::None);
  t->requires(Task::OldDW, d_labels->pNeighborListLabel,      Ghost::None);
  t->requires(Task::OldDW, d_labels->pNeighborCountLabel,     Ghost::None);

  t->requires(Task::NewDW, d_labels->pPositionLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_labels->pVelocityLabel_preReloc, Ghost::None);

  t->computes(d_labels->pParticleIDLabel_preReloc);
  t->computes(d_labels->pMassLabel_preReloc);
  t->computes(d_labels->pVolumeLabel_preReloc);
  t->computes(d_labels->pSizeLabel_preReloc);

  t->computes(d_labels->pHorizonLabel_preReloc);
  t->computes(d_labels->pNeighborCountLabel_preReloc);
  t->computes(d_labels->pNeighborListLabel_preReloc);

  //__________________________________
  //  reduction variables
  t->computes(d_labels->TotalMomentumLabel);
  t->computes(d_labels->KineticEnergyLabel);
  t->computes(d_labels->CenterOfMassPositionLabel);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  updateParticleKinematics
 *  Purpose: Copy data that will be needed for the next timestep to the new datawarehouse
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::finalizeParticleState(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset* ,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  cout_doing << "Doing finalize particle state: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Uintah::delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for (int p=0; p<patches->size(); p++) {

    const Patch* patch = patches->get(p);

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively
    double kineticEnergy = 0.0;
    Vector centerOfMassPosition(0.0,0.0,0.0);
    Vector totalMomentum(0.0,0.0,0.0);

    int numBodies = d_sharedState->getNumPeridynamicsMatls();
    for (int body = 0; body < numBodies; body++) {

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get the particle set in this patch
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      // Get the arrays of particle values needed for this task
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);

      constParticleVariable<double> pMass;
      old_dw->get(pMass, d_labels->pMassLabel, pset);

      constParticleVariable<double> pVolume;
      old_dw->get(pVolume, d_labels->pVolumeLabel, pset);

      constParticleVariable<Matrix3> pSize;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);

      constParticleVariable<double> pHorizon;
      old_dw->get(pHorizon, d_labels->pHorizonLabel, pset);

      constParticleVariable<double> pNeighborCount;
      old_dw->get(pNeighborCount, d_labels->pNeighborCountLabel, pset);

      constParticleVariable<NeighborList> pNeighborList;
      old_dw->get(pNeighborList, d_labels->pNeighborListLabel, pset);

      constParticleVariable<Point> pPosition_new;
      new_dw->get(pPosition_new, d_labels->pPositionLabel_preReloc, pset);

      constParticleVariable<Vector> pVelocity_new;
      new_dw->get(pVelocity_new, d_labels->pVelocityLabel_preReloc, pset);

      // Allocate the arrays of particle data that will be updated
      ParticleVariable<long64> pParticleID_new;
      new_dw->allocateAndPut(pParticleID_new, d_labels->pParticleIDLabel_preReloc, pset);

      ParticleVariable<double> pMass_new;
      new_dw->allocateAndPut(pMass_new, d_labels->pMassLabel_preReloc, pset);

      ParticleVariable<double> pVolume_new;
      new_dw->allocateAndPut(pVolume_new, d_labels->pVolumeLabel_preReloc, pset);

      ParticleVariable<Matrix3> pSize_new;
      new_dw->allocateAndPut(pSize_new, d_labels->pSizeLabel_preReloc, pset);

      ParticleVariable<double> pHorizon_new;
      new_dw->allocateAndPut(pHorizon_new, d_labels->pHorizonLabel_preReloc, pset);

      ParticleVariable<double> pNeighborCount_new;
      new_dw->allocateAndPut(pNeighborCount_new, d_labels->pNeighborCountLabel_preReloc, pset);

      ParticleVariable<NeighborList> pNeighborList_new;
      new_dw->allocateAndPut(pNeighborList_new, d_labels->pNeighborListLabel_preReloc, pset);

      // Copy old data to new arrays
      pParticleID_new.copyData(pParticleID);
      pMass_new.copyData(pMass);
      pVolume_new.copyData(pVolume);
      pSize_new.copyData(pSize);
      pHorizon_new.copyData(pHorizon);
      pNeighborCount_new.copyData(pNeighborCount);
      pNeighborList_new.copyData(pNeighborList);

      // Compute energy and momentum for particles in this patch
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {
        particleIndex idx = *iter;

        kineticEnergy += 0.5*pMass[idx]*pVelocity_new[idx].length2();
        centerOfMassPosition += (pPosition_new[idx]*pMass[idx]).asVector();
        totalMomentum  += pVelocity_new[idx]*pMass[idx];

      } // end particle loop

      //  reduction variables
      new_dw->put(Uintah::sumvec_vartype(totalMomentum),        d_labels->TotalMomentumLabel);
      new_dw->put(Uintah::sum_vartype(kineticEnergy),           d_labels->KineticEnergyLabel);
      new_dw->put(Uintah::sumvec_vartype(centerOfMassPosition), d_labels->CenterOfMassPositionLabel);

    }  // end of matl loop

  } // end of patch loop
  
}

bool 
Peridynamics::needRecompile(double , double , const Uintah::GridP& )
{
  cout_doing << "Doing need recompile: Peridynamics " << __FILE__ << ":" << __LINE__ << std::endl;
  if(d_recompile){
    d_recompile = false;
    return true;
  }
  else{
    return false;
  }
}

