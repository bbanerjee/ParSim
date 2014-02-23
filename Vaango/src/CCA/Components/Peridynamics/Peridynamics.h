#ifndef Vaango_Peridynamics_H
#define Vaango_Peridynamics_H

#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <Core/Grid/SimulationStateP.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Components/MPM/Contact/Contact.h>

#include <Core/Parallel/UintahParallelComponent.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/ParticleInterpolator.h>

#include <Core/Geometry/Vector.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>


namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   Peridynamics
    \brief   Parallel version of bond-based/state-based peridynamics
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class Peridynamics : public Uintah::SimulationInterface, 
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

    virtual void scheduleApplyExternalLoads(Uintah::SchedulerP& sched, 
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls);

    virtual void scheduleApplyContactLoads(Uintah::SchedulerP& sched, 
                                           const Uintah::PatchSet* patches,
                                           const Uintah::MaterialSet* matls);

    virtual void scheduleComputeInternalForce(Uintah::SchedulerP& sched, 
                                              const Uintah::PatchSet* patches,
                                              const Uintah::MaterialSet* matls);

    virtual void scheduleComputeAccStrainEnergy(Uintah::SchedulerP& sched, 
                                                const Uintah::PatchSet* patches,
                                                const Uintah::MaterialSet* matls);

    virtual void scheduleComputeAndIntegrateAcceleration(Uintah::SchedulerP& sched,
                                                         const Uintah::PatchSet* patches,
                                                         const Uintah::MaterialSet* matls);

    virtual void scheduleCorrectContactLoads(Uintah::SchedulerP& sched, 
                                             const Uintah::PatchSet* patches,
                                             const Uintah::MaterialSet* matls);

    virtual void scheduleComputeDamage(Uintah::SchedulerP& sched, 
                                       const Uintah::PatchSet* patches,
                                       const Uintah::MaterialSet* matls);

    virtual void scheduleSetGridBoundaryConditions(Uintah::SchedulerP& sched, 
                                                   const Uintah::PatchSet* patches,
                                                   const Uintah::MaterialSet* matls);
                                                 
    virtual void scheduleUpdateParticleState(Uintah::SchedulerP& sched,
                                             const Uintah::PatchSet* patches,
                                             const Uintah::MaterialSet* matls);

  protected:

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
    //virtual void computeDeformationGradient(const Uintah::ProcessorGroup*,
    //                                        const Uintah::PatchSubset* patches,
    //                                        const Uintah::MaterialSubset* matls,
    //                                        Uintah::DataWarehouse* old_dw,
    //                                        Uintah::DataWarehouse* new_dw);

    /*! Computation of stress tensor */
    //virtual void computeStressTensor(const Uintah::ProcessorGroup*,
    //                                 const Uintah::PatchSubset* patches,
    //                                 const Uintah::MaterialSubset* matls,
    //                                 Uintah::DataWarehouse* old_dw,
    //                                 Uintah::DataWarehouse* new_dw);

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

    virtual void computeDamage(const Uintah::ProcessorGroup*,
                               const Uintah::PatchSubset* patches,
                               const Uintah::MaterialSubset* ,
                               Uintah::DataWarehouse* old_dw,
                               Uintah::DataWarehouse* new_dw);


    virtual void computeAccStrainEnergy(const Uintah::ProcessorGroup*,
                                        const Uintah::PatchSubset*,
                                        const Uintah::MaterialSubset*,
                                        Uintah::DataWarehouse* old_dw,
                                        Uintah::DataWarehouse* new_dw);

    /*! Need taskgraph recompile ? */  
    bool needRecompile(double time, double dt, const Uintah::GridP& grid);

    void materialProblemSetup(const Uintah::ProblemSpecP& prob_spec);

    template<typename T>
      void setParticleDefault(Uintah::ParticleVariable<T>& pvar,
                              const Uintah::VarLabel* label, 
                              Uintah::ParticleSubset* pset,
                              Uintah::DataWarehouse* new_dw,
                              const T& val)
    {
      new_dw->allocateAndPut(pvar, label, pset);
      Uintah::ParticleSubset::iterator iter = pset->begin();
      for (; iter != pset->end(); iter++) {
        pvar[*iter] = val;
      }
    }

  protected:
  
    Uintah::SimulationStateP d_sharedState;

    PeridynamicsLabel* d_periLabels;
    PeridynamicsFlags* d_periFlags;

    Uintah::ParticleInterpolator* d_interpolator;
    Uintah::Output* d_dataArchiver;
    Uintah::Contact* d_contactModel;

    int  d_numGhostNodes;      // Number of ghost nodes needed
    int  d_numGhostParticles;  // Number of ghost particles needed
    bool d_recompile;
  
  private:

    /*! Don't allow copying */
    Peridynamics(const Peridynamics&);
    Peridynamics& operator=(const Peridynamics&);
};
      
} // end namespace Vaango

#endif
