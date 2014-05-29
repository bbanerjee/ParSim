#include <CCA/Components/Peridynamics/unit_test/utBlank.h>

#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsDomainBoundCond.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
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

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern SCIRun::Mutex cerrLock;


/*! Construct */
utBlank::utBlank(const Uintah::ProcessorGroup* myworld) :
  Uintah::UintahParallelComponent(myworld)
{
  d_periLabels = scinew PeridynamicsLabel();
  d_periFlags = scinew PeridynamicsFlags(myworld);
  d_interpolator = scinew Uintah::LinearInterpolator();

  d_dataArchiver = 0;
  d_numGhostNodes = 1;
  d_numGhostParticles = 2;
  d_recompile = false;
}

/*! Delete */
utBlank::~utBlank()
{
  delete d_periLabels;
  delete d_periFlags;
  delete d_interpolator;
}

/*! Read input file and set up problem */
void 
utBlank::problemSetup(const Uintah::ProblemSpecP& prob_spec, 
                      const Uintah::ProblemSpecP& restart_prob_spec,
                      Uintah::GridP& grid,
                      Uintah::SimulationStateP& sharedState)
{
}

void 
utBlank::outputProblemSpec(Uintah::ProblemSpecP& root_ps)
{
}

void 
utBlank::scheduleInitialize(const Uintah::LevelP& level,
                            Uintah::SchedulerP& sched)
{
}

void 
utBlank::restartInitialize()
{
}

void 
utBlank::scheduleComputeStableTimestep(const Uintah::LevelP& level,
                                       Uintah::SchedulerP& sched)
{
}

void
utBlank::scheduleTimeAdvance(const Uintah::LevelP & level,
                             Uintah::SchedulerP & sched)
{
}

void 
utBlank::scheduleInterpolateParticlesToGrid(Uintah::SchedulerP& sched,
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleApplyExternalLoads(Uintah::SchedulerP& sched,
                                    const Uintah::PatchSet* patches,
                                    const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleApplyContactLoads(Uintah::SchedulerP& sched,
                                   const Uintah::PatchSet* patches,
                                   const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleComputeInternalForce(Uintah::SchedulerP& sched,
                                      const Uintah::PatchSet* patches,
                                      const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleComputeAccStrainEnergy(Uintah::SchedulerP& sched,
                                        const Uintah::PatchSet* patches,
                                        const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleComputeAndIntegrateAcceleration(Uintah::SchedulerP& sched,
                                                 const Uintah::PatchSet* patches,
                                                 const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleCorrectContactLoads(Uintah::SchedulerP& sched,
                                     const Uintah::PatchSet* patches,
                                     const Uintah::MaterialSet* matls)
{
}

void 
utBlank::scheduleSetGridBoundaryConditions(Uintah::SchedulerP& sched,
                                           const Uintah::PatchSet* patches,
                                           const Uintah::MaterialSet* matls)

{
}

void 
utBlank::scheduleUpdateParticleState(Uintah::SchedulerP& sched,
                                     const Uintah::PatchSet* patches,
                                     const Uintah::MaterialSet* matls)

{
}

void 
utBlank::scheduleComputeDamage(Uintah::SchedulerP& sched,
                               const Uintah::PatchSet* patches,
                               const Uintah::MaterialSet* matls)
{
}

void 
utBlank::actuallyInitialize(const Uintah::ProcessorGroup*,
                            const Uintah::PatchSubset* patches,
                            const Uintah::MaterialSubset* matls,
                            Uintah::DataWarehouse*,
                            Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::actuallyComputeStableTimestep(const Uintah::ProcessorGroup*,
                                       const Uintah::PatchSubset* patches,
                                       const Uintah::MaterialSubset* ,
                                       Uintah::DataWarehouse* old_dw,
                                       Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::interpolateParticlesToGrid(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset* ,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::applyExternalLoads(const Uintah::ProcessorGroup* ,
                            const Uintah::PatchSubset* patches,
                            const Uintah::MaterialSubset*,
                            Uintah::DataWarehouse* old_dw,
                            Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::computeInternalForce(const Uintah::ProcessorGroup*,
                              const Uintah::PatchSubset* patches,
                              const Uintah::MaterialSubset* ,
                              Uintah::DataWarehouse* old_dw,
                              Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::computeAccStrainEnergy(const Uintah::ProcessorGroup*,
                                const Uintah::PatchSubset*,
                                const Uintah::MaterialSubset*,
                                Uintah::DataWarehouse* old_dw,
                                Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::computeAndIntegrateAcceleration(const Uintah::ProcessorGroup*,
                                         const Uintah::PatchSubset* patches,
                                         const Uintah::MaterialSubset*,
                                         Uintah::DataWarehouse* old_dw,
                                         Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::setGridBoundaryConditions(const Uintah::ProcessorGroup*,
                                   const Uintah::PatchSubset* patches,
                                   const Uintah::MaterialSubset* ,
                                   Uintah::DataWarehouse* old_dw,
                                   Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::updateParticleState(const Uintah::ProcessorGroup*,
                             const Uintah::PatchSubset* patches,
                             const Uintah::MaterialSubset* ,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw)
{
}

void 
utBlank::computeDamage(const Uintah::ProcessorGroup*,
                       const Uintah::PatchSubset* patches,
                       const Uintah::MaterialSubset* ,
                       Uintah::DataWarehouse* old_dw,
                       Uintah::DataWarehouse* new_dw)
{
}

bool 
utBlank::needRecompile(double , double , const Uintah::GridP& )
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
utBlank::materialProblemSetup(const Uintah::ProblemSpecP& prob_spec) 
{
}
