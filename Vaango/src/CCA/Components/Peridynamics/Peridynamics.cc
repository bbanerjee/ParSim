/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <CCA/Components/Peridynamics/Peridynamics.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>

#include <CCA/Components/Peridynamics/PeridynamicsDomainBoundCond.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCFactory.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>
#include <CCA/Components/Peridynamics/ContactModels/ContactModelFactory.h>

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
using Uintah::LevelP;
using Uintah::GridP;
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
using Uintah::Point;
using Uintah::Vector;
using Uintah::IntVector;

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern Uintah::Mutex cerrLock;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDDoing:+,PDDebug:+".....
//  bash     : export SCI_DEBUG="PDDoing:+,PDDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDDoing", false);
static DebugStream cout_dbg("PDDebug", false);
static DebugStream dbg_extra("PDDebugExtra", false);

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
  delete d_contactModel;

  ParticleLoadBCFactory::clean();

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

  // Set up load bcs on particles, if any
  ParticleLoadBCFactory::create(restart_mat_ps);
 
  // Creates Peridynamics material w/ constitutive models and damage models
  // Looks for <MaterialProperties> and then <Peridynamics>
  materialProblemSetup(restart_mat_ps, grid);

  // Create the family computer
  d_familyComputer = scinew FamilyComputer(d_flags, d_labels);

  // Create the deformation gradient computer object
  d_defGradComputer = scinew PeridynamicsDefGradComputer(d_flags, d_labels);

  // Create the internal force computer objects
  d_bondIntForceComputer = scinew BondInternalForceComputer(d_flags, d_labels);
  d_intForceComputer = scinew ParticleInternalForceComputer(d_flags, d_labels);

  // Set up contact model 
  d_contactModel = 
    Vaango::ContactModelFactory::create(UintahParallelComponent::d_myworld, restart_mat_ps, sharedState, 
                                        d_labels, d_flags);
}

void 
Peridynamics::materialProblemSetup(const ProblemSpecP& prob_spec,
                                   const GridP grid) 
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
    PeridynamicsMaterial *mat = scinew PeridynamicsMaterial(ps, grid, d_sharedState, d_flags);

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

  ProblemSpecP part_bc_ps = root->appendChild("ParticleBC");
  ProblemSpecP load_bc_ps = part_bc_ps->appendChild("Load");
  for (auto iter = ParticleLoadBCFactory::particleLoadBCs.begin();
            iter != ParticleLoadBCFactory::particleLoadBCs.end();
            iter++) {
    (*iter)->outputProblemSpec(load_bc_ps);
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
Peridynamics::scheduleInitialize(const LevelP& level,
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
  t->computes(d_labels->pLoadCurveIDLabel);

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

  // Add another task that will initialize particle load bcs
  scheduleInitializeParticleLoadBCs(level, sched);

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
 *  Method:  scheduleInitializeParticleLoadBCs
 *  Purpose: Initialize load bcs for particles
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleInitializeParticleLoadBCs(const LevelP& level,
                                                SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();
  
  d_loadCurveIndex = scinew MaterialSubset();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int curLoadBCID = 0;
  for (auto iter = ParticleLoadBCFactory::particleLoadBCs.begin();
            iter != ParticleLoadBCFactory::particleLoadBCs.end();
            iter++) {
    d_loadCurveIndex->add(curLoadBCID);
    curLoadBCID++;
  }

  if (ParticleLoadBCFactory::particleLoadBCs.size() > 0) {

    // Create a task that calculates the total number of particles
    // associated with each load curve.  
    Task* t = scinew Task("Peridynamics::countSurfaceParticlesPerLoadCurve",
                          this, &Peridynamics::countSurfaceParticlesPerLoadCurve);
    t->requires(Task::NewDW, d_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_labels->surfaceParticlesPerLoadCurveLabel, d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_sharedState->allPeridynamicsMaterials());

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("Peridynamics::initializeParticleLoadBC",
                    this, &Peridynamics::initializeParticleLoadBC);
    t->requires(Task::NewDW, d_labels->pPositionLabel,                 Ghost::None);
    t->requires(Task::NewDW, d_labels->pSizeLabel,                     Ghost::None);
    t->requires(Task::NewDW, d_labels->pDefGradLabel,                  Ghost::None);
    t->requires(Task::NewDW, d_labels->pLoadCurveIDLabel,              Ghost::None);
    t->requires(Task::NewDW, d_labels->surfaceParticlesPerLoadCurveLabel,
                             d_loadCurveIndex, Task::OutOfDomain,      Ghost::None);
    t->modifies(d_labels->pExternalForceLabel);

    sched->addTask(t, patches, d_sharedState->allPeridynamicsMaterials());
  }

  if (d_loadCurveIndex->removeReference()) {
    delete d_loadCurveIndex;
  }
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
 *  Method:  countSurfaceParticlesPerLoadCurve
 *  Purpose: Calculate the number of material points per load curve
 *           Assumes that each surface particle has been assigned a load curve index.
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::countSurfaceParticlesPerLoadCurve(const ProcessorGroup*,
                                                const PatchSubset* patches,
                                                const MaterialSubset*,
                                                DataWarehouse* ,
                                                DataWarehouse* new_dw)
{
  // Find the number of pressure BCs in the problem
  int curLoadBCID = 0;
  for (auto iter = ParticleLoadBCFactory::particleLoadBCs.begin();
            iter != ParticleLoadBCFactory::particleLoadBCs.end(); iter++) {

    curLoadBCID++;

    // Loop through the patches and count
    for (int p=0; p < patches->size(); p++) {

      const Patch* patch = patches->get(p);
      int numPeridynamicsMatls = d_sharedState->getNumPeridynamicsMatls();

      int numPts = 0;
      for (int m = 0; m < numPeridynamicsMatls; m++) {

        PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial( m );
        int matlIndex = matl->getDWIndex();
        ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);

        constParticleVariable<int> pLoadCurveID;
        new_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);

        for (auto piter = pset->begin(); piter != pset->end(); piter++) {
          particleIndex idx = *piter;
          if (pLoadCurveID[idx] == curLoadBCID) ++numPts;
        }
      } // matl loop

      new_dw->put(Uintah::sumlong_vartype(numPts), 
                  d_labels->surfaceParticlesPerLoadCurveLabel, 0, curLoadBCID-1);

    }  // patch loop
  } // bc loop
}

/*----------------------------------------------------------------------------------------
 *  Method:  initializeParticleLoadBC
 *  Purpose: Use the type of LoadBC to find the initial external force at each particle
 *----------------------------------------------------------------------------------------
 */
void 
Peridynamics::initializeParticleLoadBC(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* ,
                                       DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;
  cout_dbg << "Current Time (Initialize Pressure BC) = " << time << std::endl;


  // Calculate the force vector at each particle
  int curLoadBCID  = 0;
  for (auto iter = ParticleLoadBCFactory::particleLoadBCs.begin();
            iter != ParticleLoadBCFactory::particleLoadBCs.end(); iter++) {

    // Get the type of the BC
    std::string bc_type = (*iter)->getType();

    // Get the surface particles per load curve 
    Uintah::sumlong_vartype numPart = 0;
    new_dw->get(numPart, d_labels->surfaceParticlesPerLoadCurveLabel, 0, curLoadBCID);
    (*iter)->numParticlesOnLoadSurface(numPart);
    curLoadBCID = (*iter)->loadCurveID();
    cout_dbg << "    Load Curve = " << curLoadBCID << " Num Particles = " << numPart << std::endl;

    // Loop through the patches and calculate the force vector
    // at each particle
    for(int p=0;p<patches->size();p++){
      const Patch* patch = patches->get(p);
      int numPeridynamicsMatls=d_sharedState->getNumPeridynamicsMatls();
      for(int m = 0; m < numPeridynamicsMatls; m++){

        // Get the particle set for this material
        PeridynamicsMaterial* mpm_matl = d_sharedState->getPeridynamicsMaterial( m );
        int matlIndex = mpm_matl->getDWIndex();
        ParticleSubset* pset = new_dw->getParticleSubset(matlIndex, patch);

        constParticleVariable<Point> px;
        new_dw->get(px, d_labels->pPositionLabel, pset);

        //constParticleVariable<Matrix3> psize;
        //new_dw->get(psize, d_labels->pSizeLabel, pset);

        //constParticleVariable<Matrix3> pDefGrad;
        //new_dw->get(pDefGrad, d_labels->pDefGradLabel, pset);

        constParticleVariable<int> pLoadCurveID;
        new_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);

        ParticleVariable<Vector> pExternalForce;
        new_dw->getModifiable(pExternalForce, d_labels->pExternalForceLabel, pset);

        if (bc_type == "Pressure") {

          ParticlePressureBC* bc = dynamic_cast<ParticlePressureBC*>(*iter);

          // Calculate the force per particle at t = 0.0
          double forcePerPart = bc->forcePerParticle(time);

          // Assign initial external force
          for (auto piter = pset->begin(); piter != pset->end(); piter++) {
            particleIndex idx = *piter;
            if (pLoadCurveID[idx] == curLoadBCID) {
              pExternalForce[idx] = bc->getForceVector(px[idx], forcePerPart, time);
            }
          }
        } else if (bc_type == "NormalForce") {
          ParticleNormalForceBC* bc = dynamic_cast<ParticleNormalForceBC*>(*iter);

          // Get load from load curve
          double load = bc->getLoad(time);

          // Assign initial external force
          for (auto piter = pset->begin(); piter != pset->end(); piter++) {
            particleIndex idx = *piter;
            if (pLoadCurveID[idx] == curLoadBCID) {
              pExternalForce[idx] = bc->getForceVector(px[idx], load, time);
            }
          }
        } else if (bc_type == "Force") {
          ParticleForceBC* bc = dynamic_cast<ParticleForceBC*>(*iter);

          // Get load from load curve
          Uintah::Vector load = bc->getLoad(time);
          std::cout << "Particle BC: at t = 0: Force = " << load << std::endl;

          // Assign initial external force
          for (auto piter = pset->begin(); piter != pset->end(); piter++) {
            particleIndex idx = *piter;
            std::cout << "\t Particle = " << idx << " loadCurveID = " << pLoadCurveID[idx] 
                     << " BCid = " << curLoadBCID << std::endl;
            if (pLoadCurveID[idx] == curLoadBCID) {
              pExternalForce[idx] = load;
            }
          }
        } else {
          // Assign initial external force
          for (auto piter = pset->begin(); piter != pset->end(); piter++) {
            particleIndex idx = *piter;
            if (pLoadCurveID[idx] == curLoadBCID) {
              pExternalForce[idx] = Uintah::Vector(0, 0, 0);
            }
          }
        }
      } // matl loop
    }  // patch loop
  } // bc loop
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
Peridynamics::scheduleComputeStableTimestep(const LevelP& level,
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
Peridynamics::scheduleTimeAdvance(const LevelP & level,
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
  scheduleContactMomentumExchangeAfterInterpolate(sched, patches, matls);

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
  // on the grid
  scheduleComputeAndIntegrateGridAcceleration(sched, patches, matls);

  // Compute the accleration and integrate to find velocity and displacement
  // directly on the particles
  scheduleComputeAndIntegrateParticleAcceleration(sched, patches, matls);

  // Project the acceleration and velocity to the grid
  scheduleProjectParticleAccelerationToGrid(sched, patches, matls);

  // Correct the contact forces
  scheduleContactMomentumExchangeAfterIntegration(sched, patches, matls);

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
  // **NOTE** The <>Label_preReloc data are copied into <>Label variables
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
  t->requires(Task::OldDW, d_labels->pSizeLabel,              Ghost::None);
  t->requires(Task::OldDW, d_labels->pMassLabel,              Ghost::None);
  t->requires(Task::OldDW, d_labels->pExternalForceLabel,     Ghost::None);
  t->requires(Task::OldDW, d_labels->pLoadCurveIDLabel,     Ghost::None);
  t->computes(d_labels->pLoadCurveIDLabel_preReloc);
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
  cout_doing << "Doing apply external loads: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  // Get the current time
  double time = d_sharedState->getElapsedTime();
  cout_doing << "Current Time (applyExternalLoads) = " << time << std::endl;

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

      // Get the external force data and allocate new space for
      // external force
      constParticleVariable<Vector> pExternalForce;
      old_dw->get(pExternalForce, d_labels->pExternalForceLabel, pset);

      // Allocate space for the external force vector
      ParticleVariable<Vector>       pExternalForce_new;
      new_dw->allocateAndPut(pExternalForce_new, d_labels->pExternalForceLabel_preReloc,  pset);

      // If there are no particle load BCs
      for (auto piter = pset->begin(); piter != pset->end(); piter++){
        pExternalForce_new[*piter] = Uintah::Vector(0.0, 0.0, 0.0);
      }

      // If there are particle load BCs
      // Get the load curve data
      constParticleVariable<int> pLoadCurveID;
      old_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);

      // Recycle the loadCurveIDs
      ParticleVariable<int> pLoadCurveID_new;
      new_dw->allocateAndPut(pLoadCurveID_new, d_labels->pLoadCurveIDLabel_preReloc, pset);
      pLoadCurveID_new.copyData(pLoadCurveID);

      // Iterate over the particles
      for (auto piter = pset->begin(); piter != pset->end(); piter++) {
        particleIndex idx = *piter;
        int loadCurveID = pLoadCurveID[idx];

        // If there is no load curve associated with this particle, copy and go to
        // the next particle.  The deafult value is -1 for particles without associated
        // load curves.  See ParticleCreator.cc->getLoadCurveID.
        if (loadCurveID < 0) {
          pExternalForce_new[idx] = pExternalForce[idx];
          continue;
        }

        // If there is a load curve associated with the particle
        // Loop through particle load bc list
        int curLoadBCID = 0;
        for (auto bciter = ParticleLoadBCFactory::particleLoadBCs.begin();
                  bciter != ParticleLoadBCFactory::particleLoadBCs.end();
                  bciter++) {
          curLoadBCID = (*bciter)->loadCurveID();

          // Check that the load curve associated with the current particle
          // is equal to the id from the BC list we are iterating over
          std::cout << "\t Particle = " << idx << " loadCurveID = " << pLoadCurveID[idx] 
                     << " BCid = " << curLoadBCID << std::endl;
          if (loadCurveID != curLoadBCID) continue;

          std::string bc_type = (*bciter)->getType();
          if (bc_type == "Pressure") {
            ParticlePressureBC* pbc = dynamic_cast<ParticlePressureBC*>(*bciter);

            // Calculate the force per particle at current time
            double forcePerPart = pbc->forcePerParticle(time);

            // Update applied external force
            pExternalForce_new[idx] = pbc->getForceVector(pPosition[idx], forcePerPart, time);

          } else if (bc_type == "NormalForce") {
            ParticleNormalForceBC* pbc = dynamic_cast<ParticleNormalForceBC*>(*bciter);

            // Get load from load curve
            double load = pbc->getLoad(time);

            // Assign applied external force
            pExternalForce_new[idx] = pbc->getForceVector(pPosition[idx], load, time);

          } else if (bc_type == "Force") {
            ParticleForceBC* pbc = dynamic_cast<ParticleForceBC*>(*bciter);

            // Get load from load curve
            Uintah::Vector load = pbc->getLoad(time);
            std::cout << "Particle BC: at t = " << time << " : Force = " << load << std::endl;

            // Assign applied external force
            pExternalForce_new[idx] = load;
          } else { // Unknown BC type

            // Copy old force into new
            pExternalForce_new[idx] = pExternalForce[idx];
          }

          // A particle cannot have more than one BC.  Therefore, break out of the bc loop.
          break;

        } // end loop through BC list
      } // end loop through particles
    } // body loop
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
    auto interpolator = d_interpolator->clone(patch);
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
            gExtForce[node] += pExtForce[idx]*shapeFunction[k];
          }
          //if (gVelocity[node].length() > 0.0) {
          //  cout_dbg << " node = " << node << " gvel = " << gVelocity[node]
          //            << " pidx = " << idx << " pvel = " << pVelocity[idx]
          //            << std::endl;
          //}
        }

      } // End of particle loop

      cout_dbg << "After interpolate:" << std::endl;
      for (NodeIterator iter=patch->getExtraNodeIterator(); !iter.done();iter++) {
        IntVector node = *iter; 
        gMassGlobal[node] += gMass[node];
        gVolGlobal[node]  += gVolume[node];
        gVelGlobal[node]  += gVelocity[node];
        gVelocity[node] /= gMass[node];
        //if (gVelocity[node].length() > 0.0) {
        //  cout_dbg << " node = " << node << " gvel = " << gVelocity[node]
        //           << std::endl;
        //}
      }

    }  // End loop over materials

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done();iter++) {
      IntVector node = *iter;
      gVelGlobal[node] /= gMassGlobal[node];
    }

    // remove the interpolator clone
    //delete interpolator;

  }  // End loop over patches
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleContactMomentumExchangeAfterInterpolate
 *  Purpose: Exchange momentum between contacting objects on the grid
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleContactMomentumExchangeAfterInterpolate(SchedulerP& sched,
                                                              const PatchSet* patches,
                                                              const MaterialSet* matls)
{
  cout_doing << "Doing schedule contact momentum exchage: Peridynamics " 
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
 *  Purpose: This method sets up three tasks:
 *           1) compute the internal force on the background grid
 *           2) compute the bond internal forces for individual bonds
 *           3) a volume weighted sum of the bond internal forces to compute a particle
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
  // Task 1: Compute internal forces on the grid
  //-------------------------------------------
  Task* task1 = scinew Task("Peridynamics::computeGridInternalForce",
                            this, &Peridynamics::computeGridInternalForce);

  int numGhostCells = 1;  // Linear interpolation
  task1->requires(Task::NewDW, d_labels->gVolumeLabel, Ghost::None);
  task1->requires(Task::NewDW, d_labels->gVolumeLabel, d_sharedState->getAllInOneMatl(), 
                  Task::OutOfDomain, Ghost::None);
  task1->requires(Task::OldDW, d_labels->pStressLabel, Ghost::AroundNodes, numGhostCells);
  task1->requires(Task::OldDW, d_labels->pVolumeLabel, Ghost::AroundNodes, numGhostCells);
  task1->requires(Task::OldDW, d_labels->pPositionLabel, Ghost::AroundNodes, numGhostCells);
  task1->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::AroundNodes, numGhostCells);
  task1->requires(Task::OldDW, d_labels->pDefGradLabel, Ghost::AroundNodes, numGhostCells);

  task1->computes(d_labels->gInternalForceLabel);
  task1->computes(d_labels->gStressLabel);
  task1->computes(d_labels->gStressLabel, d_sharedState->getAllInOneMatl(), Task::OutOfDomain);
  
  sched->addTask(task1, patches, matls);

  //-------------------------------------------
  // Task 2: Compute bond internal forces
  //-------------------------------------------
  Task* task2 = scinew Task("Peridynamics::computeBondInternalForce",
                            this, &Peridynamics::computeBondInternalForce);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_bondIntForceComputer->addComputesAndRequires(task2, matl, patches);
  }

  sched->addTask(task2, patches, matls);
  
  //-------------------------------------------
  // Task 3: Compute particle internal forces
  //-------------------------------------------
  Task* task3 = scinew Task("Peridynamics::computeParticleInternalForce",
                            this, &Peridynamics::computeParticleInternalForce);

  for (int body = 0; body < numBodies; body++) {
    PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial(body);
    d_intForceComputer->addComputesAndRequires(task3, matl, patches);
  }

  sched->addTask(task3, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeGridInternalForce
 *  Purpose: Compute the internal force on the grid using linear interpolation.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeGridInternalForce(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* ,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0/dx.x();
    oodx[1] = 1.0/dx.y();
    oodx[2] = 1.0/dx.z();
    Matrix3 Id;
    Id.Identity();

    auto interpolator = d_interpolator->clone(patch); 
    std::vector<IntVector> nodeIndices(interpolator->size());
    std::vector<double> shapeFunction(interpolator->size());
    std::vector<Vector> shapeGradient(interpolator->size());

    int numPeridynamicsMatls = d_sharedState->getNumPeridynamicsMatls();

    constNCVariable<double>   gVolumeGlobal;
    new_dw->get(gVolumeGlobal,  d_labels->gVolumeLabel,
                d_sharedState->getAllInOneMatl()->get(0), patch, Ghost::None,0);

    NCVariable<Matrix3> gStressGlobal;
    new_dw->allocateAndPut(gStressGlobal, d_labels->gStressLabel, 
                           d_sharedState->getAllInOneMatl()->get(0), patch);

    int numGhostCells = 1;
    for(int m = 0; m < numPeridynamicsMatls; m++){

      PeridynamicsMaterial* mpm_matl = d_sharedState->getPeridynamicsMaterial( m );
      int matlIndex = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch,
                                                       Ghost::AroundNodes, numGhostCells,
                                                       d_labels->pPositionLabel);

      constParticleVariable<Point>   pPosition;
      old_dw->get(pPosition,  d_labels->pPositionLabel, pset);

      constParticleVariable<double>  pVolume;
      old_dw->get(pVolume, d_labels->pVolumeLabel, pset);

      constParticleVariable<Matrix3> pDefGrad;
      old_dw->get(pDefGrad, d_labels->pDefGradLabel,  pset);

      constParticleVariable<Matrix3> pStress;
      old_dw->get(pStress, d_labels->pStressLabel,  pset);

      constParticleVariable<Matrix3> pSize;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);

      constNCVariable<double>        gVolume;
      new_dw->get(gVolume, d_labels->gVolumeLabel, matlIndex, patch, Ghost::None, 0);

      NCVariable<Matrix3>            gStress;
      new_dw->allocateAndPut(gStress, d_labels->gStressLabel, matlIndex, patch);

      NCVariable<Vector>             gInternalForce;
      new_dw->allocateAndPut(gInternalForce, d_labels->gInternalForceLabel, matlIndex, patch);
      gInternalForce.initialize(Vector(0,0,0));

      for (auto iter = pset->begin(); iter != pset->end(); iter++) {

        particleIndex idx = *iter;
  
        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(pPosition[idx], nodeIndices, shapeFunction, 
                                                            shapeGradient, pSize[idx], pDefGrad[idx]);

        Matrix3 stressVol  = pStress[idx]*pVolume[idx];

        IntVector node;
        for (unsigned int k = 0; k < nodeIndices.size(); k++) {
          node = nodeIndices[k];
          if (patch->containsNode(node)) {
            Vector div(shapeGradient[k].x()*oodx[0],shapeGradient[k].y()*oodx[1],
                         shapeGradient[k].z()*oodx[2]);
            gInternalForce[node] -= (div * pStress[idx])  * pVolume[idx];
            gStress[node] += stressVol * shapeFunction[k];
          }
        }
      }

      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        gStressGlobal[c] += gStress[c];
        gStress[c] /= gVolume[c];
      }

    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gStressGlobal[c] /= gVolumeGlobal[c];
    }

    //delete interpolator;
  }
  
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
Peridynamics::computeParticleInternalForce(const ProcessorGroup*,
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
 *  Method:  scheduleComputeAndIntegrateGridAcceleration
 *  Purpose: This method sets up the required variables and the computed
 *           variables for solving the momentum equation using
 *           Forward Euler on the background grid.  The acceleration and an intermediate 
 *           velocity are computed.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeAndIntegrateGridAcceleration(SchedulerP& sched,
                                                          const PatchSet* patches,
                                                          const MaterialSet* matls)
{
  Task* t = scinew Task("Peridynamics::computeAndIntegrateGridAcceleration",
                        this, &Peridynamics::computeAndIntegrateGridAcceleration);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );

  t->requires(Task::NewDW, d_labels->gMassLabel,          Ghost::None);
  t->requires(Task::NewDW, d_labels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gExternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gVelocityLabel,      Ghost::None);

  t->computes(d_labels->gVelocityStarLabel);
  t->computes(d_labels->gAccelerationLabel);

  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeAndIntegrateGridAcceleration
 *  Purpose: Solve the momentum equations on the background grid
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeAndIntegrateGridAcceleration(const ProcessorGroup*,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset*,
                                                  DataWarehouse* old_dw,
                                                  DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    Ghost::GhostType gnone = Ghost::None;
    Vector gravity = d_flags->d_gravity;

    for(int m = 0; m < d_sharedState->getNumPeridynamicsMatls(); m++){
      PeridynamicsMaterial* matl = d_sharedState->getPeridynamicsMaterial( m );
      int matlIndex = matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<Vector> internalforce, externalforce, velocity;
      constNCVariable<double> mass;

      Uintah::delt_vartype delT;
      old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );
 
      new_dw->get(internalforce,d_labels->gInternalForceLabel, matlIndex, patch, gnone, 0);
      new_dw->get(externalforce,d_labels->gExternalForceLabel, matlIndex, patch, gnone, 0);
      new_dw->get(mass,         d_labels->gMassLabel,          matlIndex, patch, gnone, 0);
      new_dw->get(velocity,     d_labels->gVelocityLabel,      matlIndex, patch, gnone, 0);

      // Create variables for the results
      NCVariable<Vector> velocity_star,acceleration;
      new_dw->allocateAndPut(velocity_star, d_labels->gVelocityStarLabel, matlIndex, patch);
      new_dw->allocateAndPut(acceleration,  d_labels->gAccelerationLabel, matlIndex, patch);
      acceleration.initialize(Vector(0.,0.,0.));

      for(NodeIterator iter=patch->getExtraNodeIterator(); !iter.done();iter++){

        IntVector c = *iter;

        Vector acc(0.,0.,0.);
        if (mass[c] > std::numeric_limits<double>::epsilon()){
          acc  = (internalforce[c] + externalforce[c])/mass[c];
        }
        acceleration[c] = acc +  gravity;
        velocity_star[c] = velocity[c] + acceleration[c] * delT;
      }
    }    // matls
  }
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleComputeAndIntegrateParticleAcceleration
 *  Purpose: This method sets up the required variables and the computed
 *           variables for solving the momentum equation using
 *           Forward Euler directly on the particles.  The acceleration and an intermediate 
 *           velocity are computed.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleComputeAndIntegrateParticleAcceleration(SchedulerP& sched,
                                                              const PatchSet* patches,
                                                              const MaterialSet* matls)
{
  cout_doing << "Doing schedule compute and integrate acceleration: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Task* t = scinew Task("Peridynamics::computeAndIntegrateParticleAcceleration",
                        this, &Peridynamics::computeAndIntegrateParticleAcceleration);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );
  t->requires(Task::OldDW, d_labels->pMassLabel,                   Ghost::None);
  t->requires(Task::NewDW, d_labels->pVolumeLabel_preReloc,        Ghost::None);
  t->requires(Task::OldDW, d_labels->pPositionLabel,               Ghost::None);
  t->requires(Task::OldDW, d_labels->pDisplacementLabel,           Ghost::None);
  t->requires(Task::OldDW, d_labels->pVelocityLabel,               Ghost::None);
  t->requires(Task::NewDW, d_labels->pInternalForceLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_labels->pExternalForceLabel_preReloc, Ghost::None);

  int numBodies = d_sharedState->getNumPeridynamicsMatls();
  for(int body = 0; body < numBodies; body++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
    const MaterialSubset* matlset = peridynamic_matl->thisMaterial();

    t->computes(d_labels->pPositionStarLabel, matlset);
    t->computes(d_labels->pDisplacementStarLabel, matlset);
    t->computes(d_labels->pVelocityStarLabel, matlset);
    t->computes(d_labels->pAccelerationLabel, matlset);
  }
  sched->addTask(t, patches, matls);
}

/*------------------------------------------------------------------------------------------------
 *  Method:  computeAndIntegrateParticleAcceleration
 *  Purpose: Use the Peridynamics balance of momentum equation to solve for the acceleration
 *           at each particle in the patch.
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::computeAndIntegrateParticleAcceleration(const ProcessorGroup*,
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

      // Get particle set
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      // Get required variables for this patch
      constParticleVariable<double> pMass;
      old_dw->get(pMass, d_labels->pMassLabel, pset);

      constParticleVariable<double> pVolume;
      new_dw->get(pVolume, d_labels->pVolumeLabel_preReloc, pset);

      constParticleVariable<Point> pPosition;
      old_dw->get(pPosition, d_labels->pPositionLabel, pset);

      constParticleVariable<Vector> pDisp;
      old_dw->get(pDisp, d_labels->pDisplacementLabel, pset);

      constParticleVariable<Vector> pVelocity;
      old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);

      constParticleVariable<Vector> pInternalForce;
      new_dw->get(pInternalForce, d_labels->pInternalForceLabel_preReloc, pset);

      constParticleVariable<Vector> pExternalForce;
      new_dw->get(pExternalForce, d_labels->pExternalForceLabel_preReloc, pset);

      // Create variables for the results
      ParticleVariable<Point> pPosition_star;
      new_dw->allocateAndPut(pPosition_star, d_labels->pPositionStarLabel, pset);

      ParticleVariable<Vector> pDisp_star;
      new_dw->allocateAndPut(pDisp_star, d_labels->pDisplacementStarLabel, pset);

      ParticleVariable<Vector> pVelocity_star;
      new_dw->allocateAndPut(pVelocity_star, d_labels->pVelocityStarLabel, pset);

      ParticleVariable<Vector> pAcceleration;
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
          pVelocity_star[idx] = pVelocity[idx] + acc*delT;

          // Integrate the velocity to get the displacement
          // Forward Euler
          Vector disp = pVelocity[idx]*delT + 0.5*acc*delT*delT;
          pDisp_star[idx] = pDisp[idx] + disp;

          // Update position
          pPosition_star[idx] = pPosition[idx] + disp;

	  cout_dbg << "Compute particle acc: " << std::endl;
	  cout_dbg << " Particle " << idx << " "
                   << " f_int = " << internal_force_acc << " "
                   << " f_ext = " << external_force_acc << " " 
                   << " acc = " << acc << std::endl
                   << " v* = " << pVelocity_star[idx]
                   << " u* = " << pDisp_star[idx]
                   << " x* = " << pPosition_star[idx] << std::endl;
        } else {
          pAcceleration[idx] = Vector(0.0);
          pVelocity_star[idx] = Vector(0.0);
          pDisp_star[idx] = Vector(0.0);
          pPosition_star[idx] = pPosition[idx];
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

  t->computes(d_labels->gpAccelerationLabel);
  t->computes(d_labels->gpVelocityStarLabel);

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
    auto interpolator = d_interpolator->clone(patch);
    std::vector<IntVector> nodeIndices(interpolator->size());
    std::vector<double> shapeFunction(interpolator->size());

    for(int body = 0; body < numBodies; body++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Get the particle subset needed for projection to the grid
      int numGhostNodes = 1;  // Linear interpolation only
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch, Ghost::AroundNodes,
                                                       numGhostNodes,
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
      new_dw->allocateAndPut(gVelocity_star, d_labels->gpVelocityStarLabel,  matlIndex, patch);
      gVelocity_star.initialize(Vector(0.0, 0.0, 0.0));

      NCVariable<Vector> gAcceleration;
      new_dw->allocateAndPut(gAcceleration,  d_labels->gpAccelerationLabel,  matlIndex, patch);
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
    //delete interpolator;

  }  //end  patch loop
}

/*------------------------------------------------------------------------------------------------
 *  Method:  scheduleContactMomentumExchangeAfterIntegration
 *  Purpose: Use the contact model to correct any contact loads (grid based contact)
 * ------------------------------------------------------------------------------------------------
 */
void 
Peridynamics::scheduleContactMomentumExchangeAfterIntegration(SchedulerP& sched,
                                                              const PatchSet* patches,
                                                              const MaterialSet* matls)
{
  cout_doing << "Doing schedule contact momentum exchange after integration: Peridynamics " 
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
        // cout_dbg << " node = " << c << " old_vel = " << gVelocity_old[c]
        //          << " mod_vel = " << gVelocity_star[c]
        //          << " acc = " << gAcceleration[c] << std::endl;
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

  t->requires(Task::NewDW, d_labels->pPositionStarLabel,     Ghost::None);
  t->requires(Task::NewDW, d_labels->pVelocityStarLabel,     Ghost::None);
  t->requires(Task::NewDW, d_labels->pDisplacementStarLabel, Ghost::None);

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
  cout_doing << "Doing update particle kinematics: Peridynamics " 
             << ":Processor : " << UintahParallelComponent::d_myworld->myrank() << ":"
             << __FILE__ << ":" << __LINE__ << std::endl;

  Matrix3 dummyMatrix(0.0);

  Uintah::delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for (int p=0; p<patches->size(); p++) {

    const Patch* patch = patches->get(p);

    // Create a copy of the LinearInterpolator
    auto interpolator = d_interpolator->clone(patch);
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

      constParticleVariable<Point> pPosition_star;
      new_dw->get(pPosition_star, d_labels->pPositionStarLabel, pset);

      constParticleVariable<Vector> pDisp_star;
      new_dw->get(pDisp_star, d_labels->pDisplacementStarLabel, pset);

      constParticleVariable<Vector> pVelocity_star;
      new_dw->get(pVelocity_star, d_labels->pVelocityStarLabel, pset);

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
	//pPosition_new[idx] = pPosition[idx] + vel*delT;
        //pDisp_new[idx] = pDisp[idx] + vel*delT;
        //pVelocity_new[idx] = pVelocity[idx] + acc*delT;
	pPosition_new[idx] = pPosition_star[idx];
        pDisp_new[idx] = pDisp_star[idx];
        pVelocity_new[idx] = pVelocity_star[idx];

        cout_dbg << " Particle " << idx << " position = " << pPosition_new[idx]
                 << " Displacement = " << pDisp_new[idx] << " Velocity = " << pVelocity_new[idx] << std::endl;

        dbg_extra << " Particle = " << idx << " x_err = " << pPosition_star[idx] - pPosition[idx] - vel*delT
                  << " u_err = " << pDisp_star[idx] - pDisp[idx] - vel*delT
                  << " v_err = " << pVelocity_star[idx] - pVelocity[idx] - acc*delT << std::endl;
      } // end particle loop

    }  // end of matl loop

    // remove the interpolator clone
    //delete interpolator;

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
  t->requires(Task::OldDW, d_labels->pSizeLabel,              Ghost::None);

  t->requires(Task::OldDW, d_labels->pHorizonLabel,           Ghost::None);
  t->requires(Task::OldDW, d_labels->pNeighborListLabel,      Ghost::None);
  t->requires(Task::OldDW, d_labels->pNeighborCountLabel,     Ghost::None);

  t->requires(Task::NewDW, d_labels->pPositionLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_labels->pVelocityLabel_preReloc, Ghost::None);

  t->computes(d_labels->pParticleIDLabel_preReloc);
  t->computes(d_labels->pMassLabel_preReloc);
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

  for (int p=0; p<patches->size(); p++) {

    const Patch* patch = patches->get(p);

    double kineticEnergy = 0.0;
    Vector centerOfMassPosition(0.0,0.0,0.0);
    Vector totalMomentum(0.0,0.0,0.0);

    int numBodies = d_sharedState->getNumPeridynamicsMatls();
    for (int body = 0; body < numBodies; body++) {

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(body);
      int matlIndex = peridynamic_matl->getDWIndex();

      // Not populating the delset, but we need this to satisfy Relocate
      //ParticleSubset* delset = scinew ParticleSubset(0, matlIndex, patch);
      //new_dw->deleteParticles(delset);      

      // Get the particle set in this patch
      ParticleSubset* pset = old_dw->getParticleSubset(matlIndex, patch);

      // Create a delete set that's needed during particle relocation
      ParticleSubset* delset = scinew ParticleSubset(0, matlIndex, patch);

      // Get the arrays of particle values needed for this task
      cout_dbg << "Before get particle id." << std::endl;
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);
      cout_dbg << "Got particle id." << std::endl;

      constParticleVariable<double> pMass;
      old_dw->get(pMass, d_labels->pMassLabel, pset);
      cout_dbg << "Got particle mass." << std::endl;

      constParticleVariable<Matrix3> pSize;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      cout_dbg << "Got particle size." << std::endl;

      constParticleVariable<double> pHorizon;
      old_dw->get(pHorizon, d_labels->pHorizonLabel, pset);
      cout_dbg << "Got particle horizon." << std::endl;

      constParticleVariable<int> pNeighborCount;
      old_dw->get(pNeighborCount, d_labels->pNeighborCountLabel, pset);
      cout_dbg << "Got particle neighbor count." << std::endl;

      constParticleVariable<NeighborList> pNeighborList;
      old_dw->get(pNeighborList, d_labels->pNeighborListLabel, pset);
      cout_dbg << "Got particle neighbor list." << std::endl;

      constParticleVariable<Point> pPosition_new;
      new_dw->get(pPosition_new, d_labels->pPositionLabel_preReloc, pset);
      cout_dbg << "Got particle position." << std::endl;

      constParticleVariable<Vector> pVelocity_new;
      new_dw->get(pVelocity_new, d_labels->pVelocityLabel_preReloc, pset);
      cout_dbg << "Got particle velocity." << std::endl;

      // Allocate the arrays of particle data that will be updated
      ParticleVariable<long64> pParticleID_new;
      new_dw->allocateAndPut(pParticleID_new, d_labels->pParticleIDLabel_preReloc, pset);

      ParticleVariable<double> pMass_new;
      new_dw->allocateAndPut(pMass_new, d_labels->pMassLabel_preReloc, pset);

      ParticleVariable<Matrix3> pSize_new;
      new_dw->allocateAndPut(pSize_new, d_labels->pSizeLabel_preReloc, pset);

      ParticleVariable<double> pHorizon_new;
      new_dw->allocateAndPut(pHorizon_new, d_labels->pHorizonLabel_preReloc, pset);

      ParticleVariable<int> pNeighborCount_new;
      new_dw->allocateAndPut(pNeighborCount_new, d_labels->pNeighborCountLabel_preReloc, pset);

      ParticleVariable<NeighborList> pNeighborList_new;
      new_dw->allocateAndPut(pNeighborList_new, d_labels->pNeighborListLabel_preReloc, pset);

      // Copy old data to new arrays
      pParticleID_new.copyData(pParticleID);
      pMass_new.copyData(pMass);
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

        // Delete particles whose mass is too small (**WARNING** hardcoded for now)
        if (pMass[idx] <= 1.0e-16) {
          std::cout << "Deleted particle with id " << idx << ":" << pParticleID[idx]
                    << " of material type " << matlIndex 
                    << " from patch " << patch << std::endl;
          delset->addParticle(idx);
        }


      } // end particle loop

      //  reduction variables
      new_dw->put(Uintah::sumvec_vartype(totalMomentum),        d_labels->TotalMomentumLabel);
      new_dw->put(Uintah::sum_vartype(kineticEnergy),           d_labels->KineticEnergyLabel);
      new_dw->put(Uintah::sumvec_vartype(centerOfMassPosition), d_labels->CenterOfMassPositionLabel);

      // particles to be deleted during relocation
      new_dw->deleteParticles(delset);

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

