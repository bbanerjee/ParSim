#include <CCA/Components/Peridynamics/PeridynamicsMaterialModel/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModel/DamageModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/Peridynamics/Peridynamics.h>
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
#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
#include <Core/Util/DebugStream.h>
#include <Core/Thread/Mutex.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

using namespace Vaango;

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern Mutex cerrLock;


/*! Construct */
Peridynamics::Peridynamics(const Uintah::ProcessorGroup* myworld) :
  Uintah::UintahParallelComponent(myworld)
{
  d_periLabels = SCIRun::scinew PeridynamicsLabel();
  d_periFlags = SCIRun::scinew PeridynamicsFlags(myworld);
  d_interpolator = SCIRun::scinew Uintah::LinearInterpolator();

  d_dataArchiver = 0;
  d_numGhostNodes = 1;
  d_numGhostParticles = 2;
  d_recompile = false;
}

/*! Delete */
Peridynamics::~Peridynamics()
{
  delete d_periLabels;
  delete d_periFlags;
  delete d_interpolator;
}

/*! Read input file and set up problem */
void 
Peridynamics::problemSetup(const Uintah::ProblemSpecP& prob_spec, 
                           const Uintah::ProblemSpecP& restart_prob_spec,
                           Uintah::GridP& grid,
                           Uintah::SimulationStateP& sharedState)
{
  std::cout << "Doing problemSetup: Peridynamics" << std::endl;

  d_sharedState = sharedState;
  dynamic_cast<Uintah::Scheduler*>(getPort("scheduler"))->setPositionVar(d_periLabels->pPositionLabel);
  
  d_dataArchiver = dynamic_cast<Uintah::Output*>(getPort("output"));
  if(!d_dataArchiver){
    throw Uintah::InternalError("Peridynamics:couldn't get output port", __FILE__, __LINE__);
  }

  Uintah::ProblemSpecP prob_spec_mat_ps = prob_spec->findBlockWithOutAttribute("MaterialProperties");
  Uintah::ProblemSpecP restart_mat_ps = prob_spec;
  if (restart_prob_spec) restart_mat_ps = restart_prob_spec;

  Uintah::ProblemSpecP peridynamic_soln_ps = restart_mat_ps->findBlock("Peridynamics");
  if (!peridynamic_soln_ps){
    std::ostringstream warn;
    warn<<"ERROR:Peridynamics:\n missing Peridynamics section in the input file\n";
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
 
  // Read all Peridynamics d_periFlags (look in PeridynamicsFlags.cc)
  d_periFlags->readPeridynamicsFlags(restart_mat_ps, d_dataArchiver);
  d_sharedState->setParticleGhostLayer(Uintah::Ghost::AroundNodes, d_numGhostParticles);

  // Creates Peridynamics material w/ constitutive models and damage models
  materialProblemSetup(restart_mat_ps);

  // Set up contact model
  d_contactModel = 
    ContactFactory::create(UintahParallelComponent::d_myworld, restart_mat_ps, sharedState, 
                           d_periLabels, d_periFlags);
}

void 
Peridynamics::outputProblemSpec(ProblemSpecP& root_ps)
{
  Uintah::ProblemSpecP root = root_ps->getRootNode();

  Uintah::ProblemSpecP d_periFlags_ps = root->appendChild("Peridynamics");
  d_periFlags->outputProblemSpec(d_periFlags_ps);

  Uintah::ProblemSpecP mat_ps = 0;
  mat_ps = root->findBlockWithOutAttribute("MaterialProperties");

  if (mat_ps == 0)
    mat_ps = root->appendChild("MaterialProperties");
    
  Uintah::ProblemSpecP peridynamic_ps = mat_ps->appendChild("Peridynamics");
  for (int i = 0; i < d_sharedState->getNumPeridynamicsMatls();i++) {
    PeridynamicsMaterial* mat = d_sharedState->getPeridynamicsMaterial(i);
    Uintah::ProblemSpecP cm_ps = mat->outputProblemSpec(peridynamic_ps);
  }

  d_contactModel->outputProblemSpec(peridynamic_ps);
}

void 
Peridynamics::scheduleInitialize(const Uintah::LevelP& level,
                                 Uintah::SchedulerP& sched)
{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::actuallyInitialize",
                                        this, &Peridynamics::actuallyInitialize);

  const Uintah::PatchSet* patches = level->eachPatch();
  Uintah::MaterialSubset* zeroth_matl = SCIRun::scinew Uintah::MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_periLabels->pPositionLabel);
  t->computes(d_periLabels->pMassLabel);
  t->computes(d_periLabels->pVolumeLabel);
  t->computes(d_periLabels->pDispLabel);
  t->computes(d_periLabels->pVelocityLabel);
  t->computes(d_periLabels->pExternalForceLabel);
  t->computes(d_periLabels->pParticleIDLabel);
  t->computes(d_periLabels->pVelGradLabel);
  t->computes(d_periLabels->pDefGradLabel);
  t->computes(d_periLabels->pDispGradLabel);
  t->computes(d_periLabels->pStressLabel);
  t->computes(d_sharedState->get_delt_label(),level.get_rep());
  t->computes(d_periLabels->pCellNAPIDLabel,zeroth_matl);

  int numPeridynamicsMats = d_sharedState->getNumPeridynamicsMatls();
  for(int m = 0; m < numPeridynamicsMats; m++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);

    // Add constitutive model computes
    PeridynamicsMaterialModel* cm = peridynamic_matl->getPeridynamicsMaterialModel();
    cm->addInitialComputesAndRequires(t, peridynamic_matl, patches);

    // Add damage model computes
    PeridynamicsFailureModel* fm = peridynamic_matl->getPeridynamicsFailureModel();
    fm->addInitialComputesAndRequires(t, peridynamic_matl, patches, d_periLabels);
  }

  sched->addTask(t, patches, d_sharedState->allPeridynamicsMaterials());

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) delete zeroth_matl; // shouln't happen, but...
}

/* _____________________________________________________________________
 Purpose:   Set variables that are normally set during the initialization
            phase, but get wiped clean when you restart
_____________________________________________________________________*/
void 
Peridynamics::restartInitialize()
{
  std::cout_doing<<"Doing restartInitialize\t\t\t\t\t Peridynamics"<<std::endl;
}

void 
Peridynamics::scheduleComputeStableTimestep(const Uintah::LevelP& level,
                                            Uintah::SchedulerP& sched)
{
  // Nothing to do here - delt is computed as a by-product of the
  // constitutive model
  // However, this task needs to do something in the case that Peridynamics
  // is being run on more than one level.
  Uintah::Task* t = 0;
  t = SCIRun::scinew   Uintah::Task("Peridynamics::actuallyComputeStableTimestep",
                            this, &Peridynamics::actuallyComputeStableTimestep);

  const Uintah::MaterialSet* peridynamic_matls = d_sharedState->allPeridynamicsMaterials();
  t->computes(d_sharedState->get_delt_label(),level.get_rep());
  sched->addTask(t,level->eachPatch(), peridynamic_matls);
}

void
Peridynamics::scheduleTimeAdvance(const Uintah::LevelP & level,
                                  Uintah::SchedulerP & sched)
{
  const Uintah::PatchSet* patches = level->eachPatch();
  const Uintah::MaterialSet* matls = d_sharedState->allPeridynamicsMaterials();
  const Uintah::MaterialSubset* peridynamic_matls_sub = matls->getUnion();

  scheduleInterpolateParticlesToGrid(sched, patches, matls);
  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleApplyContactLoads(sched, patches, matls);
  scheduleComputeInternalForce(sched, patches, matls);
  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  scheduleCorrectContactLoads(sched, patches, matls);
  scheduleSetGridBoundaryConditions(sched, patches, matls);
  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  scheduleComputeDamage(sched, patches, matls);
}

void 
Peridynamics::scheduleInterpolateParticlesToGrid(Uintah::SchedulerP& sched,
                                                 const Uintah::PatchSet* patches,
                                                 const Uintah::MaterialSet* matls)
{
  Uintah::Task* t = SCIRun::scinew   Uintah::Task("Peridynamics::interpolateParticlesToGrid",
                                          this, &Peridynamics::interpolateParticlesToGrid);

  t->requires(Uintah::Task::OldDW, d_periLabels->pPositionLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pMassLabel,     Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pVolumeLabel,   Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pVelocityLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pDefGradLabel,  Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::NewDW, d_periLabels->pExternalForceLabel_preReloc, 
              Uintah::Ghost::AroundNodes, d_numGhostParticles);

  t->computes(d_periLabels->gMassLabel);
  t->computes(d_periLabels->gVolumeLabel);
  t->computes(d_periLabels->gVelocityLabel);
  t->computes(d_periLabels->gExternalForceLabel);
  t->computes(d_periLabels->gMassLabel,     d_sharedState->getAllInOneMatl(), Uintah::Task::OutOfDomain);
  t->computes(d_periLabels->gVolumeLabel,   d_sharedState->getAllInOneMatl(), Uintah::Task::OutOfDomain);
  t->computes(d_periLabels->gVelocityLabel, d_sharedState->getAllInOneMatl(), Uintah::Task::OutOfDomain);

  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleApplyExternalLoads(Uintah::SchedulerP& sched,
                                         const Uintah::PatchSet* patches,
                                         const Uintah::MaterialSet* matls)
{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::applyExternalLoads",
                                        this, &Peridynamics::applyExternalLoads);
                  
  t->requires(Task::OldDW, d_periLabels->pPositionLabel,          Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pMassLabel,              Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pDispLabel,              Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pDefGradLabel,           Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pExternalForceLabel,     Ghost::None);

  t->computes(d_periLabels->pExternalForceLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleApplyContactLoads(Uintah::SchedulerP& sched,
                                        const Uintah::PatchSet* patches,
                                        const Uintah::MaterialSet* matls)
{
  d_contactModel->addComputesAndRequiresInterpolated(sched, patches, matls);
}

void 
Peridynamics::scheduleComputeInternalForce(Uintah::SchedulerP& sched,
                                           const Uintah::PatchSet* patches,
                                           const Uintah::MaterialSet* matls)
{
  Uintah::Task* t = SCIRun::scinew   Uintah::Task("Peridynamics::computeInternalForce",
                                          this, &Peridynamics::computeInternalForce);

  t->requires(Uintah::Task::OldDW, d_periLabels->pPositionLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pVolumeLabel,   Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pDefGradLabel,  Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pDispLabel,     Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::OldDW, d_periLabels->pVelocityLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);

  t->requires(Uintah::Task::NewDW, d_periLabels->gVolumeLabel,   Uintah::Ghost::None);
  t->requires(Uintah::Task::NewDW, d_periLabels->gVolumeLabel, d_sharedState->getAllInOneMatl(),
              Uintah::Task::OutOfDomain, Uintah::Ghost::None);

  int numMatls = d_sharedState->getNumPeridynamicsMatls();
  for(int m = 0; m < numMatls; m++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);
    const Uintah::MaterialSubset* matlset = peridynamic_matl->thisMaterial();

    t->computes(d_periLabels->pInternalForceLabel_preReloc, matlset);
    t->computes(d_periLabels->pStressLabel_preReloc, matlset);
    t->computes(d_periLabels->pDefGradLabel_preReloc, matlset);
    t->computes(d_periLabels->pVelGradLabel_preReloc, matlset);
  }

  t->computes(d_periLabels->StrainEnergyLabel);
  t->computes(d_sharedState->get_delt_label(),getLevel(patches));

  sched->addTask(t, patches, matls);

  scheduleComputeAccStrainEnergy(sched, patches, matls);
}

// Compute the accumulated strain energy
void 
Peridynamics::scheduleComputeAccStrainEnergy(Uintah::SchedulerP& sched,
                                             const Uintah::PatchSet* patches,
                                             const Uintah::MaterialSet* matls)
{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::computeAccStrainEnergy",
                                        this, &Peridynamics::computeAccStrainEnergy);
  t->requires(Uintah::Task::OldDW, d_periLabels->AccStrainEnergyLabel);
  t->requires(Uintah::Task::NewDW, d_periLabels->StrainEnergyLabel);
  t->computes(d_periLabels->AccStrainEnergyLabel);
  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleComputeAndIntegrateAcceleration(Uintah::SchedulerP& sched,
                                                      const Uintah::PatchSet* patches,
                                                      const Uintah::MaterialSet* matls)
{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::computeAndIntegrateAcceleration",
                                        this, &Peridynamics::computeAndIntegrateAcceleration);

  t->requires(Uintah::Task::OldDW, d_sharedState->get_delt_label() );

  t->requires(Uintah::Task::OldDW, d_periLabels->pMassLabel,          Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::NewDW, d_periLabels->pInternalForceLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::NewDW, d_periLabels->pExternalForceLabel, Uintah::Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Uintah::Task::NewDW, d_periLabels->pVelocityLabel,      Uintah::Ghost::AroundNodes, d_numGhostParticles);

  int numMatls = d_sharedState->getNumPeridynamicsMatls();
  for(int m = 0; m < numMatls; m++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);
    const Uintah::MaterialSubset* matlset = peridynamic_matl->thisMaterial();

    t->computes(d_periLabels->pVelocityStarLabel, matlset);
    t->computes(d_periLabels->pAccelerationLabel, matlset);
  }
  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleCorrectContactLoads(Uintah::SchedulerP& sched,
                                          const Uintah::PatchSet* patches,
                                          const Uintah::MaterialSet* matls)
{
  d_contactModel->addComputesAndRequiresIntegrated(sched, patches, matls);
}

void 
Peridynamics::scheduleSetGridBoundaryConditions(Uintah::SchedulerP& sched,
                                                const Uintah::PatchSet* patches,
                                                const Uintah::MaterialSet* matls)

{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::setGridBoundaryConditions",
                                        this, &Peridynamics::setGridBoundaryConditions);
                  
  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_sharedState->get_delt_label() );
  
  t->modifies(             d_periLabels->gAccelerationLabel,     mss);
  t->modifies(             d_periLabels->gVelocityStarLabel,     mss);
  t->requires(Task::NewDW, d_periLabels->gVelocityLabel,   Ghost::None);

  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleInterpolateToParticlesAndUpdate(Uintah::SchedulerP& sched,
                                                      const Uintah::PatchSet* patches,
                                                      const Uintah::MaterialSet* matls)

{
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::interpolateToParticlesAndUpdate",
                                        this, &Peridynamics::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );

  t->requires(Task::NewDW, d_periLabels->gAccelerationLabel, Uintah::Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::NewDW, d_periLabels->gVelocityStarLabel, Uintah::Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::OldDW, d_periLabels->pPositionLabel,     Uintah::Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pMassLabel,         Uintah::Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pParticleIDLabel,   Uintah::Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pVelocityLabel,     Uintah::Ghost::None);
  t->requires(Task::OldDW, d_periLabels->pDispLabel,         Uintah::Ghost::None);
  t->requires(Task::NewDW, d_periLabels->pDefGradLabel_preReloc, Uintah::Ghost::None);
  t->modifies(d_periLabels->pVolumeLabel_preReloc);

  t->computes(d_periLabels->pDispLabel_preReloc);
  t->computes(d_periLabels->pVelocityLabel_preReloc);
  t->computes(d_periLabels->pPositionLabel_preReloc);
  t->computes(d_periLabels->pParticleIDLabel_preReloc);
  t->computes(d_periLabels->pMassLabel_preReloc);

  //__________________________________
  //  reduction variables
  t->computes(d_periLabels->TotalMomentumLabel);
  t->computes(d_periLabels->KineticEnergyLabel);
  t->computes(d_periLabels->CenterOfMassPositionLabel);
  t->computes(d_periLabels->TotalMassLabel);

  sched->addTask(t, patches, matls);
}

void 
Peridynamics::scheduleComputecDamage(Uintah::SchedulerP& sched,
                                     const Uintah::PatchSet* patches,
                                     const Uintah::MaterialSet* matls)
{
  int numMatls = d_sharedState->getNumPeridynamicsMatls();
  Uintah::Task* t = SCIRun::scinew Uintah::Task("Peridynamics::computeDamage",
                                        this, &Peridynamics::computeDamage);
  for(int m = 0; m < numMatls; m++){
    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);

    PeridynamicsFailureModel* d_failureModel = peridynamic_matl->getPeridynamicsFailureModel();
    d_failureModel->addComputesAndRequires(t, peridynamic_matl, patches, d_periLabels);
  }

  sched->addTask(t, patches, matls);
}

void 
Peridynamics::actuallyInitialize(const Uintah::ProcessorGroup*,
                                 const Uintah::PatchSubset* patches,
                                 const Uintah::MaterialSubset* matls,
                                 Uintah::DataWarehouse*,
                                 Uintah::DataWarehouse* new_dw)
{
  particleIndex totalParticles=0;
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    Uintah::CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_periLabels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    for(int m=0;m<matls->size();m++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial( m );
      int indx = peridynamic_matl->getDWIndex();

      // Create particles
      particleIndex numParticles = peridynamic_matl->countParticles(patch);
      totalParticles+=numParticles;
      peridynamic_matl->createParticles(numParticles, cellNAPID, patch, new_dw);

      // Initialize constitutive model
      peridynamic_matl->getPeridynamicsMaterialModel()->initialize(patch, peridynamic_matl, new_dw);

      // Initialize damage model
      peridynamic_matl->getPeridynamicsFailureModel()->initialize(patch, peridynamic_matl, new_dw, d_periLabels);
    }
  }

  // Initialize the accumulated strain energy
  new_dw->put(max_vartype(0.0), d_periLabels->AccStrainEnergyLabel);
  new_dw->put(Uintah::sumlong_vartype(totalParticles), d_periLabels->partCountLabel);
}

void 
Peridynamics::actuallyComputeStableTimestep(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* ,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw)
{
  // Put something here to satisfy the need for a reduction operation in
  // the case that there are multiple levels present
  const Uintah::Level* level = getLevel(patches);
  new_dw->put(delt_vartype(999.0), d_periLabels->delTLabel, level);
}

void 
Peridynamics::interpolateParticlesToGrid(const Uintah::ProcessorGroup*,
                                         const Uintah::PatchSubset* patches,
                                         const Uintah::MaterialSubset* ,
                                         Uintah::DataWarehouse* old_dw,
                                         Uintah::DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Uintah::Patch* patch = patches->get(p);
    int numMatls = d_sharedState->getNumPeridynamicsMatls();

    // Use only linear interpolation for particle to grid 
    Uintah::ParticleInterpolator* interpolator = d_interpolator->clone(patch);
    std::vector<Uintah::IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    // Initialize global data
    Uintah::NCVariable<double>         gMassGlobal, gVolGlobal;
    Uintah::NCVariable<Uintah::Vector> gVelGlobal;
    new_dw->allocateAndPut(gMassGlobal, d_periLabels->gMassLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    new_dw->allocateAndPut(gVolGlobal, d_periLabels->gVolumeLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    new_dw->allocateAndPut(gVelGlobal, d_periLabels->gVelocityLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gMassGlobal.initialize(std::numeric_limits<double>::epsilon());
    gVolGlobal.initialize(std::numeric_limits<double>::epsilon());
    gVelGlobal.initialize(Uintah::Vector(0.0));

    for(int m = 0; m < numMatls; m++){

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);
      int dwi = peridynamic_matl->getDWIndex();

      // Create arrays for the particle data
      Uintah::constParticleVariable<Uintah::Point>  pPosition;
      Uintah::constParticleVariable<double> pMass, pVolume;
      Uintah::constParticleVariable<Uintah::Vector> pVelocity, pExtForce;
      Uintah::constParticleVariable<Uintah::Matrix3> pSize, pDefGrad_old;

      Uintah::ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch, Uintah::Ghost::AroundNodes, 
                                                               d_numGhostParticles, 
                                                               d_periLabels->pPositionLabel);

      old_dw->get(pPosition,     d_periLabels->pPositionLabel,      pset);
      old_dw->get(pMass,         d_periLabels->pMassLabel,          pset);
      old_dw->get(pVolume,       d_periLabels->pVolumeLabel,        pset);
      old_dw->get(pVelocity,     d_periLabels->pVelocityLabel,      pset);
      old_dw->get(pDefGrad_old,  d_periLabels->pDefGradLabel,       pset);
      old_dw->get(pSize,         d_periLabels->pSizeLabel,          pset);

      new_dw->get(pExtForce,     d_periLabels->pExtForceLabel_preReloc, pset);

      // Create arrays for the grid data
      NCVariable<double> gMass;
      NCVariable<double> gVol;
      NCVariable<Vector> gVel;
      NCVariable<Vector> gExtForce;

      new_dw->allocateAndPut(gMass,     d_periLabels->gMassLabel,          dwi, patch);
      new_dw->allocateAndPut(gVolume,   d_periLabels->gVolumeLabel,        dwi, patch);
      new_dw->allocateAndPut(gVelocity, d_periLabels->gVelocityLabel,      dwi, patch);
      new_dw->allocateAndPut(gExtForce, d_periLabels->gExternalForceLabel, dwi, patch);

      gMass.initialize(std::numeric_limits<double>::epsilon());
      gVolume.initialize(std::numeric_limits<double>::epsilon());
      gVelocity.initialize(Uintah::Vector(0,0,0));
      gExtForce.initialize(Uintah::Vector(0,0,0));

      // Interpolate particle data to Grid data.
      // This currently consists of the particle velocity and mass
      // Need to compute the lumped Global mass matrix and velocity
      // Uintah::Vector from the individual mass matrix and velocity vector
      // GridMass * GridVelocity =  S^T*M_D*ParticleVelocity
      Uintah::Vector totalMomentum(0.0,0.0,0.0);
      Uintah::Vector pMomentum;

      //loop over all particles in the patch:
      for (Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;
        interpolator->findCellAndWeights(pPosition[idx], ni, S, pSize[idx], pDefGrad_old[idx]);
        pMomentum = pVelocity[idx]*pMass[idx];
        totalMomentum += pMomentum;

        // Add each particles contribution to the local mass & velocity 
        // Must use the node indices
        Uintah::IntVector node;
        // Iterates through the nodes which receive information from the current particle
        for (unsigned int k = 0; k < ni.size(); k++) {
          node = ni[k];
          if(patch->containsNode(node)) {
            gMass[node] += pMass[idx]*S[k];
            gVelocity[node] += pMomentum*S[k];
            gVolume[node] += pVolume[idx]*S[k];
          }
        }

      } // End of particle loop

      for(Uintah::NodeIterator iter=patch->getExtraNodeIterator(); !iter.done();iter++){
        Uintah::IntVector c = *iter; 
        gMassGlobal[c] += gMass[c];
        gVolGlobal[c]  += gVolume[c];
        gVelGlobal[c]  += gVelocity[c];
        gVelocity[c] /= gMass[c];
      }

    }  // End loop over materials

    for(NodeIterator iter = patch->getNodeIterator(); !iter.done();iter++){
      Uintah::IntVector c = *iter;
      gVelGlobal[c] /= gMassGlobal[c];
    }

    delete interpolator;
  }  // End loop over patches
}

void 
Peridynamics::applyExternalLoads(const Uintah::ProcessorGroup* ,
                                 const Uintah::PatchSubset* patches,
                                 const Uintah::MaterialSubset*,
                                 Uintah::DataWarehouse* old_dw,
                                 Uintah::DataWarehouse* new_dw)
{
  // Get the current time
  double time = d_sharedState->getElapsedTime();

  // Loop thru patches to update external force vector
  for(int p=0;p<patches->size();p++){
    const Uintah::Patch* patch = patches->get(p);
    
    int numPeridynamicsMatls=d_sharedState->getNumPeridynamicsMatls();
    
    for(int m = 0; m < numPeridynamicsMatls; m++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial( m );
      int dwi = peridynamic_matl->getDWIndex();
      Uintah::ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Get the particle data
      Uintah::constParticleVariable<Uintah::Point>   pPosition;
      Uintah::constParticleVariable<Uintah::Matrix3> pSize;
      Uintah::constParticleVariable<Uintah::Matrix3> pDefGrad;
      Uintah::ParticleVariable<Uintah::Vector>       pExternalForce_new;

      old_dw->get(pPosition, d_periLabels->pPositionLabel, pset);
      old_dw->get(pSize,     d_periLabels->pSizeLabel,     pset);
      old_dw->get(pDefGrad,  d_periLabels->pDefGradLabel,  pset);

      new_dw->allocateAndPut(pExternalForce_new, d_periLabels->pExtForceLabel_preReloc,  pset);

      for(Uintah::ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        pExternalForce_new[*iter] = 0.;
      }
    } // matl loop
  }  // patch loop
}

void 
Peridynamics::computeInternalForce(const Uintah::ProcessorGroup*,
                                   const Uintah::PatchSubset* patches,
                                   const Uintah::MaterialSubset* ,
                                   Uintah::DataWarehouse* old_dw,
                                   Uintah::DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Uintah::Patch* patch = patches->get(p);

    Uintah::Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0/dx.x();
    oodx[1] = 1.0/dx.y();
    oodx[2] = 1.0/dx.z();
    Uintah::Matrix3 Id;
    Id.Identity();

    int numPeridynamicsMatls = d_sharedState->getNumPeridynamicsMatls();

    for(int m = 0; m < numPeridynamicsMatls; m++){

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial( m );
      int dwi = peridynamic_matl->getDWIndex();

      Uintah::constParticleVariable<Uintah::Point>   pPosition;
      Uintah::constParticleVariable<double>          pVolume;
      Uintah::constParticleVariable<Uintah::Matrix3> pStress;
      Uintah::constParticleVariable<Uintah::Matrix3> pSize;
      Uintah::constParticleVariable<Uintah::Matrix3> pDefGrad_old;
      Uintah::ParticleVariable<Uintah::Vector>       pInternalForce;
      Uintah::ParticleVariable<Uintah::Matrix3>      pStress_new;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch,
                                                       Ghost::AroundNodes, d_numGhostParticles,
                                                       d_periLabels->pPositionLabel);

      old_dw->get(pPosition,      d_periLabels->pPositionLabel,   pset);
      old_dw->get(pVolume,        d_periLabels->pVolumeLabel,     pset);
      old_dw->get(pStress,        d_periLabels->pStressLabel,     pset);
      old_dw->get(pSize,          d_periLabels->pSizeLabel,       pset);
      old_dw->get(pDefGrad_old,   d_periLabels->pDefGradLabel,    pset);

      new_dw->allocateAndPut(pStress_new,        d_periLabels->pStress_preReloc,             pset);
      new_dw->allocateAndPut(pInternalForce,     d_periLabels->pInternalForceLabel_preReloc, pset);
      pInternalForce.initialize(Uintah::Vector(0,0,0));

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;
        pStress_new[idx] = Uintah::Matrix3(0.0);
        pInternalForce[idx] = Uintah::Vector(0.0);
      } // end particle loop
    } // end matl loop
  } // end patch loop
  
}

void 
Peridynamics::computeAccStrainEnergy(const Uintah::ProcessorGroup*,
                                     const Uintah::PatchSubset*,
                                     const Uintah::MaterialSubset*,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw)
{
  // Get the totalStrainEnergy from the old datawarehouse
  Uintah::max_vartype accStrainEnergy;
  old_dw->get(accStrainEnergy, d_periLabels->AccStrainEnergyLabel);

  // Get the incremental strain energy from the new datawarehouse
  Uintah::sum_vartype incStrainEnergy;
  new_dw->get(incStrainEnergy, d_periLabels->StrainEnergyLabel);
  
  // Add the two a put into new dw
  double totalStrainEnergy = 
    (double) accStrainEnergy + (double) incStrainEnergy;
  new_dw->put(Uintah::max_vartype(totalStrainEnergy), d_periLabels->AccStrainEnergyLabel);
}

void 
Peridynamics::computeAndIntegrateAcceleration(const Uintah::ProcessorGroup*,
                                              const Uintah::PatchSubset* patches,
                                              const Uintah::MaterialSubset*,
                                              Uintah::DataWarehouse* old_dw,
                                              Uintah::DataWarehouse* new_dw)
{
  Uintah::delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );
 
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    Uintah::Vector gravity = d_periFlags->d_gravity;
    for(int m = 0; m < d_sharedState->getNumPeridynamicsMatls(); m++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial( m );
      int dwi = peridynamic_matl->getDWIndex();

      // Get required variables for this patch
      Uintah::constParticleVariable<Uintah::Vector> pInternalForce, pExternalForce, pVelocity;
      Uintah::constParticleVariable<double> pMass;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch,
                                                       Ghost::AroundNodes, d_numGhostParticles,
                                                       d_periLabels->pPositionLabel);

      new_dw->get(pInternalForce, d_periLabels->pInternalForceLabel, pset);
      new_dw->get(pExternalForce, d_periLabels->pExternalForceLabel, pset);
      new_dw->get(pMass,          d_periLabels->pMassLabel,          pset);
      new_dw->get(pVelocity,      d_periLabels->pVelocityLabel,      pset);

      // Create variables for the results
      Uintah::ParticleVariable<Uintah::Vector> pVelocity_star, pAcceleration;
      new_dw->allocateAndPut(pVelocity_star, d_periLabels->pVelocityStarLabel, pset);
      new_dw->allocateAndPut(pAcceleration,  d_periLabels->pAccelerationLabel, pset);
      pAcceleration.initialize(Uintah::Vector(0.,0.,0.));

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;
        pAcceleration[idx] = Uintah::Vector(0.0) + gravity;
        pVelocity_star[idx] = pVelocity[idx] + pAcceleration[idx] * delT;
      } // end particle loop
    } // end matls loop
  } // end patch loop
}

void 
Peridynamics::setGridBoundaryConditions(const Uintah::ProcessorGroup*,
                                        const Uintah::PatchSubset* patches,
                                        const Uintah::MaterialSubset* ,
                                        Uintah::DataWarehouse* old_dw,
                                        Uintah::DataWarehouse* new_dw)
{
  delt_vartype delT;            
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);

    int numPeridynamicsMatls = d_sharedState->getNumPeridynamicsMatls();

    for(int m = 0; m < numPeridynamicsMatls; m++){
      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial( m );
      int dwi = peridynamic_matl->getDWIndex();

      Uintah::NCVariable<Vector> gVelocity_star, gAcceleration;
      Uintah::constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(gAcceleration,  d_periLabels->gAccelerationLabel,  dwi,patch);
      new_dw->getModifiable(gVelocity_star, d_periLabels->gVelocityStarLabel,  dwi,patch);
      new_dw->get(gVelocity,                d_periLabels->gVelocityLabel,      dwi, patch, Ghost::None, 0);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      PeridynamicsBoundCond bc;
      bc.setBoundaryCondition(patch, dwi, "Velocity", gVelocity_star); 
      bc.setBoundaryCondition(patch, dwi, "Symmetric", gVelocity_star);

      // Now recompute acceleration as the difference between the velocity
      // interpolated to the grid (no bcs applied) and the new velocity_star
      for(Uintah::NodeIterator iter=patch->getExtraNodeIterator();!iter.done(); iter++){
        Uintah::IntVector c = *iter;
        gAcceleration[c] = (gVelocity_star[c] - gVelocity[c])/delT;
      } // node loop
    } // matl loop
  }  // patch loop
}

void 
Peridynamics::interpolateToParticlesAndUpdate(const Uintah::ProcessorGroup*,
                                              const Uintah::PatchSubset* patches,
                                              const Uintah::MaterialSubset* ,
                                              Uintah::DataWarehouse* old_dw,
                                              Uintah::DataWarehouse* new_dw)
{
  delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for(int p=0;p<patches->size();p++) {
    const Uintah::Patch* patch = patches->get(p);

    Uintah::ParticleInterpolator* interpolator = d_interpolator->clone(patch);
    std::vector<Uintah::IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively
    double totalMass = 0.0;
    double totalVolume = 0.0;
    double kineticEnergy = 0.0;
    Uintah::Vector centerOfMassPosition(0.0,0.0,0.0);
    Uintah::Vector totalMomentum(0.0,0.0,0.0);

    int numPeridynamicsMatls=d_sharedState->getNumPeridynamicsMatls();
    for(int m = 0; m < numPeridynamicsMatls; m++) {

      PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);
      int dwi = peridynamic_matl->getDWIndex();

      // Get the arrays of particle values to be changed
      Uintah::constParticleVariable<Uintah::Point> pPosition;
      Uintah::constParticleVariable<Uintah::Vector> pVelocity;
      Uintah::constParticleVariable<Uintah::Matrix3> pSize;
      Uintah::constParticleVariable<double> pMass;
      Uintah::constParticleVariable<double> pVolume;
      Uintah::constParticleVariable<long64> pids;
      Uintah::constParticleVariable<Uintah::Vector> pDisp;
      Uintah::constParticleVariable<Uintah::Matrix3> pDefGrad_new,pDefGrad_old;

      Uintah::ParticleVariable<Uintah::Point> pPosition_new;
      Uintah::ParticleVariable<Uintah::Vector> pVelocity_new;
      Uintah::ParticleVariable<Uintah::Matrix3> pSize_new;
      Uintah::ParticleVariable<double> pMass_new, pVolume_new;
      Uintah::ParticleVariable<long64> pids_new;
      Uintah::ParticleVariable<Uintah::Vector> pDisp_new;

      // Get the arrays of grid data on which the new part. values depend
      Uintah::constNCVariable<Uintah::Vector> gVelocity_star, gAcceleration;

      Uintah::ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(pPosition,    d_periLabels->pPositionLabel,                  pset);
      old_dw->get(pDisp,        d_periLabels->pDispLabel,                      pset);
      old_dw->get(pMass,        d_periLabels->pMassLabel,                      pset);
      old_dw->get(pVelocity,    d_periLabels->pVelocityLabel,                  pset);
      old_dw->get(pVolume,      d_periLabels->pVolumeLabel,                    pset);

      new_dw->get(pDefGrad_new,       d_periLabels->pDefGradLabel_preReloc,    pset);
      new_dw->getModifiable(pVolume,  d_periLabels->pVolumeLabel_preReloc,     pset);

      new_dw->allocateAndPut(pVelocity_new, d_periLabels->pVelocityLabel_preReloc,   pset);
      new_dw->allocateAndPut(pPosition_new, d_periLabels->pPositionLabel_preReloc,   pset);
      new_dw->allocateAndPut(pDisp_new,     d_periLabels->pDispLabel_preReloc,       pset);
      new_dw->allocateAndPut(pMass_new,     d_periLabels->pMassLabel_preReloc,       pset);

      //Carry forward ParticleID
      old_dw->get(pids,                d_periLabels->pParticleIDLabel,          pset);
      new_dw->allocateAndPut(pids_new, d_periLabels->pParticleIDLabel_preReloc, pset);
      pids_new.copyData(pids);

      //Carry forward ParticleSize
      old_dw->get(pSize,                d_periLabels->pSizeLabel,                pset);
      new_dw->allocateAndPut(pSize_new, d_periLabels->pSizeLabel_preReloc,       pset);
      pSize_new.copyData(pSize);

      new_dw->get(gVelocity_star,  d_periLabels->gVelocityStarLabel, dwi, patch,
                  Uintah::Ghost::AroundCells,d_numGhostParticles);
      new_dw->get(gAcceleration,   d_periLabels->gAccelerationLabel, dwi, patch,
                  Uintah::Ghost::AroundCells,d_numGhostParticles);

      // Loop over particles
      for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pPosition[idx], ni, S, pSize[idx], pDefGrad_new[idx]);

        // Accumulate the contribution from each surrounding vertex
        Uintah::Vector vel(0.0,0.0,0.0);
        Uintah::Vector acc(0.0,0.0,0.0);
        for (unsigned int kk = 0; kk < ni.size(); kk++) {
          Uintah::IntVector node = ni[kk];
          vel += gVelocity_star[node]  * S[kk];
          acc += gAcceleration[node]   * S[kk];
        }

        // Update the particle's position and velocity
	pPosition_new[idx]    = pPosition[idx] + vel*delT;
        pDisp_new[idx]        = pDisp[idx] + vel*delT;
        pVelocity_new[idx]    = pVelocity[idx] + acc*delT;

	pMass_new[idx] = pMass[idx];
	pVolume[idx]   = pVolume_new[idx];

        kineticEnergy += 0.5*pMass[idx]*pVelocity_new[idx].length2();
        centerOfMassPosition += (pPosition_new[idx]*pMass[idx]).asVector();
        totalMomentum  += pVelocity_new[idx]*pMass[idx];
        totalMass += pMass[idx];
        totalVolume += pVolume[idx];

      } // end particle loop

      //  reduction variables
      new_dw->put(Uintah::sum_vartype(totalMass),               d_periLabels->TotalMassLabel);
      new_dw->put(Uintah::sum_vartype(totalVolume),             d_periLabels->TotalVolumeDeformedLabel);
      new_dw->put(Uintah::sumvec_vartype(totalMomentum),        d_periLabels->TotalMomentumLabel);
      new_dw->put(Uintah::sum_vartype(kineticEnergy),           d_periLabels->KineticEnergyLabel);
      new_dw->put(Uintah::sumvec_vartype(centerOfMassPosition), d_periLabels->CenterOfMassPositionLabel);

    }  // end of matl loop

    delete interpolator;
  } // end of patch loop
  
}

void 
Peridynamics::computeDamage(const Uintah::ProcessorGroup*,
                            const Uintah::PatchSubset* patches,
                            const Uintah::MaterialSubset* ,
                            Uintah::DataWarehouse* old_dw,
                            Uintah::DataWarehouse* new_dw)
{
  for(int m = 0; m < d_sharedState->getNumPeridynamicsMatls(); m++){

    PeridynamicsMaterial* peridynamic_matl = d_sharedState->getPeridynamicsMaterial(m);

    PeridynamicsFailureModel* failureModel = peridynamic_matl->getPeridynamicsFailureModel();
    failureModel->computeDamage(patches, peridynamic_matl, old_dw, new_dw, d_periLabels);

  }
}

bool 
Peridynamics::needRecompile(double , double , const Uintah::GridP& )
{
  if(d_recompile){
    d_recompile = false;
    return true;
  }
  else{
    return false;
  }
}

void 
Peridynamics::materialProblemSetup(const Uintah::ProblemSpecP& prob_spec) 
{
  //Search for the MaterialProperties block and then get the Peridynamics section
  Uintah::ProblemSpecP mat_ps = prob_spec->findBlockWithOutAttribute("MaterialProperties");
  Uintah::ProblemSpecP peridynamics_mat_ps = mat_ps->findBlock("Peridynamics");

  for (Uintah::ProblemSpecP ps = peridynamics_mat_ps->findBlock("material"); ps != 0;
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
    PeridynamicsMaterial *mat = SCIRun::scinew PeridynamicsMaterial(ps, d_sharedState, d_periFlags);

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
