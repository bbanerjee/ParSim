#ifndef Vaango_Peridynamics_H
#define Vaango_Peridynamics_H

#include <Core/Parallel/UintahParallelComponent.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Components/Peridynamics/PeridynamicsCommon.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Labels/PeridynamicsLabel.h>

#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>

#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   Peridynamics
    \brief   Parallel version of bond-based/state-based peridynamics
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class Peridynamics : public PeridynamicsCommon, 
                       public Uintah::SimulationInterface, 
                       public Uintah::UintahParallelComponent 
  {
  public:

    /*! Constructor/Destructor */
    Peridynamics(const Uintah::ProcessorGroup* myworld);
    virtual ~Peridynamics();

    
    /*! Problem spec reader and particle creation */
    virtual void problemSetup(const Uintah::ProblemSpecP& params, 
                              const Uintah::ProblemSpecP& restart_prob_spec, 
                              Uintah::GridP& grid,
                              Uintah::SimulationStateP& state);

    /*! Output writer for problem spec */
    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    /*! Schedule tasks */
    virtual void scheduleInitialize(const Uintah::LevelP& level, 
                                    Uintah::SchedulerP& sched);
    virtual void scheduleComputeStableTimestep(const Uintah::LevelP& level, 
                                               Uintah::SchedulerP& sched);
    virtual void scheduleTimeAdvance(const Uintah::LevelP& level, 
                                     Uintah::SchedulerP& sched);

    /*! For restarts */
    virtual void restartInitialize();

  protected:
 
    /*! Schedule tasks called by the public ones */
    virtual void scheduleInterpolateParticlesToGrid(Uintah::SchedulerP& sched, 
                                                    const Uintah::PatchSet* patches,
                                                    const Uintah::MaterialSet* matls);
    virtual void scheduleExMomInterpolated(Uintah::SchedulerP& sched, 
                                           const Uintah::PatchSet* patches,
                                           const Uintah::MaterialSet* matls);

    void scheduleComputeDeformationGradient(Uintah::SchedulerP& sched, 
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls);
  
    virtual void scheduleComputeStressTensor(Uintah::SchedulerP& sched, 
                                             const Uintah::PatchSet* patches,
                                             const Uintah::MaterialSet* matls);

    void scheduleComputeDamage(Uintah::SchedulerP& sched, 
                               const Uintah::PatchSet* patches,
                               const Uintah::MaterialSet* matls);

    virtual void scheduleComputeInternalForce(Uintah::SchedulerP& sched, 
                                              const Uintah::PatchSet* patches,
                                              const Uintah::MaterialSet* matls);

    virtual void scheduleComputeAndIntegrateAcceleration(Uintah::SchedulerP& sched,
                                                         const Uintah::PatchSet* patches,
                                                         const Uintah::MaterialSet* matls);

    virtual void scheduleExMomIntegrated(Uintah::SchedulerP& sched, 
                                         const Uintah::PatchSet* patches,
                                         const Uintah::MaterialSet* matls);

    void scheduleSetGridBoundaryConditions(Uintah::SchedulerP& sched, 
                                           const Uintah::PatchSet* patches,
                                           const Uintah::MaterialSet* matls);
                                                 
    void scheduleApplyExternalLoads(Uintah::SchedulerP& sched, 
                                    const Uintah::PatchSet* patches, 
                                    const Uintah::MaterialSet* matls);

    virtual void scheduleUpdateParticleState(Uintah::SchedulerP& sched, 
                                             const Uintah::PatchSet* patches,
                                             const Uintah::MaterialSet* matls);

    /*! Actual initialization */
    virtual void actuallyInitialize(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset* matls,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);

    /*! Actual stable timestep computation */
    void actuallyComputeStableTimestep(const Uintah::ProcessorGroup*,
                                       const Uintah::PatchSubset* patches,
                                       const Uintah::MaterialSubset* matls,
                                       Uintah::DataWarehouse* old_dw,
                                       Uintah::DataWarehouse* new_dw);

    /*! Interpolation of particles to grid for contact detection */
    virtual void interpolateParticlesToGrid(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* matls,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw);

    /*! Computation of deformation gradient */
    virtual void computeDeformationGradient(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* matls,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw);

    /*! Computation of stress tensor */
    virtual void computeStressTensor(const Uintah::ProcessorGroup*,
                                     const Uintah::PatchSubset* patches,
                                     const Uintah::MaterialSubset* matls,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw);

    /*! Computation of internal force */
    virtual void computeInternalForce(const Uintah::ProcessorGroup*,
                                      const Uintah::PatchSubset* patches,
                                      const Uintah::MaterialSubset* matls,
                                      Uintah::DataWarehouse* old_dw,
                                      Uintah::DataWarehouse* new_dw);

    /*! Integration of acceleration */
    virtual void computeAndIntegrateAcceleration(const Uintah::ProcessorGroup*,
                                                 const Uintah::PatchSubset* patches,
                                                 const Uintah::MaterialSubset* matls,
                                                 Uintah::DataWarehouse* old_dw,
                                                 Uintah::DataWarehouse* new_dw);

    /*! Grid boundary condition set up */
    void setGridBoundaryConditions(const Uintah::ProcessorGroup*,
                                   const Uintah::PatchSubset* patches,
                                   const Uintah::MaterialSubset* ,
                                   Uintah::DataWarehouse* old_dw,
                                   Uintah::DataWarehouse* new_dw);

    /*! Apply external loads to bodies inside the domain */
    void applyExternalLoads(const Uintah::ProcessorGroup*,
                            const Uintah::PatchSubset* patches,
                            const Uintah::MaterialSubset* ,
                            Uintah::DataWarehouse* old_dw,
                            Uintah::DataWarehouse* new_dw);

    /*! Update everything at the end of a timestep */
    virtual void updateParticleState(const Uintah::ProcessorGroup*,
                                     const Uintah::PatchSubset* patches,
                                     const Uintah::MaterialSubset* matls,
                                     Uintah::DataWarehouse* old_dw,
                                     Uintah::DataWarehouse* new_dw);

    /*! Need taskgraph recompile ? */  
    bool needRecompile(double time, double dt, const Uintah::GridP& grid);

  protected:
  
    Uintah::SimulationStateP d_sharedState;

    PeridynamicsLabel* d_periLabels;
    PeridynamicsFlags* d_periFlags;

    Uintah::Output* d_dataArchiver;

    int  d_numGhostNodes;  // Number of ghost nodes needed
    bool d_recompile;
  
  private:

    /*! Don't allow copying */
    Peridynamics(const Peridynamics&);
    Peridynamics& operator=(const Peridynamics&);
};
      
} // end namespace Vaango

#endif
