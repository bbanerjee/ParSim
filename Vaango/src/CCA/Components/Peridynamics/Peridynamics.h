/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef Vaango_Peridynamics_H
#define Vaango_Peridynamics_H

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsCommon.h>

#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/FamilyComputer/FamilyComputer.h>
#include <CCA/Components/Peridynamics/GradientComputer/PeridynamicsDefGradComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/BondInternalForceComputer.h>
#include <CCA/Components/Peridynamics/InternalForceComputer/ParticleInternalForceComputer.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>

#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MPMInterpolators/ParticleInterpolator.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleVariable.h>

#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class   Peridynamics
  \brief   Parallel version of bond-based/state-based peridynamics
  \author  Biswajit Banerjee
  \warning
*/
/////////////////////////////////////////////////////////////////////////////

class Peridynamics
  : public Uintah::SimulationCommon
  , public PeridynamicsCommon
{
public:
  /*! Constructor/Destructor */
  Peridynamics(const Uintah::ProcessorGroup* myworld,
               const Uintah::MaterialManagerP& mat_manager);

  virtual ~Peridynamics();

  /*! Don't allow copying/move */
  Peridynamics(const Peridynamics&) = delete;
  Peridynamics(Peridynamics&&)      = delete;
  Peridynamics&
  operator=(const Peridynamics&) = delete;
  Peridynamics&
  operator=(Peridynamics&&) = delete;

  /*! Problem spec reader and particle creation */
  virtual void
  problemSetup(const Uintah::ProblemSpecP& params,
               const Uintah::ProblemSpecP& restart_prob_spec,
               Uintah::GridP& grid,
               const std::string& input_ups_dir = "");

  /*! Output writer for problem spec */
  virtual void
  outputProblemSpec(Uintah::ProblemSpecP& ps);

  /*! Schedule tasks */
  virtual void
  scheduleInitialize(const Uintah::LevelP& level, Uintah::SchedulerP& sched);
  virtual void
  scheduleComputeStableTimestep(const Uintah::LevelP& level,
                                Uintah::SchedulerP& sched);
  virtual void
  scheduleTimeAdvance(const Uintah::LevelP& level, Uintah::SchedulerP& sched);

  /*! For restarts */

  virtual void
  scheduleRestartInitialize(const Uintah::LevelP& level,
                            Uintah::SchedulerP& sched)
  {
  }
  virtual void
  restartInitialize();

protected:
  /*! Schedule initialization of particle load BCs */
  void
  scheduleInitializeParticleLoadBCs(const Uintah::LevelP& level,
                                    Uintah::SchedulerP& sched);

  /*! Actual initialization */
  void
  actuallyInitialize(const Uintah::ProcessorGroup*,
                     const Uintah::PatchSubset* patches,
                     const Uintah::MaterialSubset* matls,
                     Uintah::DataWarehouse* old_dw,
                     Uintah::DataWarehouse* new_dw);

  /*!  Calculate the number of material points per load curve */
  void
  countSurfaceParticlesPerLoadCurve(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset*,
                                    Uintah::DataWarehouse*,
                                    Uintah::DataWarehouse* new_dw);

  /*!  Use the type of LoadBC to find the initial external force at each
   * particle */
  void
  initializeParticleLoadBC(const Uintah::ProcessorGroup*,
                           const Uintah::PatchSubset* patches,
                           const Uintah::MaterialSubset*,
                           Uintah::DataWarehouse*,
                           Uintah::DataWarehouse* new_dw);

  /*! Find list of neighbors */
  void
  findNeighborsInHorizon(const Uintah::ProcessorGroup*,
                         const Uintah::PatchSubset* patches,
                         const Uintah::MaterialSubset* matls,
                         Uintah::DataWarehouse*,
                         Uintah::DataWarehouse* new_dw);

  /*! Actual stable timestep computation */
  void
  actuallyComputeStableTimestep(const Uintah::ProcessorGroup*,
                                const Uintah::PatchSubset* patches,
                                const Uintah::MaterialSubset* matls,
                                Uintah::DataWarehouse* old_dw,
                                Uintah::DataWarehouse* new_dw);

  /*! Apply external loads to bodies inside the domain */
  void
  scheduleApplyExternalLoads(Uintah::SchedulerP& sched,
                             const Uintah::PatchSet* patches,
                             const Uintah::MaterialSet* matls);
  void
  applyExternalLoads(const Uintah::ProcessorGroup*,
                     const Uintah::PatchSubset* patches,
                     const Uintah::MaterialSubset*,
                     Uintah::DataWarehouse* old_dw,
                     Uintah::DataWarehouse* new_dw);

  /*! Interpolation of particles to grid for contact detection */
  void
  scheduleInterpolateParticlesToGrid(Uintah::SchedulerP& sched,
                                     const Uintah::PatchSet* patches,
                                     const Uintah::MaterialSet* matls);
  void
  interpolateParticlesToGrid(const Uintah::ProcessorGroup*,
                             const Uintah::PatchSubset* patches,
                             const Uintah::MaterialSubset* matls,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw);

  /*! Apply contact forces */
  void
  scheduleContactMomentumExchangeAfterInterpolate(
    Uintah::SchedulerP& sched,
    const Uintah::PatchSet* patches,
    const Uintah::MaterialSet* matls);

  /*! Computation of deformation gradient */
  void
  scheduleComputeDeformationGradient(Uintah::SchedulerP& sched,
                                     const Uintah::PatchSet* patches,
                                     const Uintah::MaterialSet* matls);
  void
  computeDeformationGradient(const Uintah::ProcessorGroup*,
                             const Uintah::PatchSubset* patches,
                             const Uintah::MaterialSubset* matls,
                             Uintah::DataWarehouse* old_dw,
                             Uintah::DataWarehouse* new_dw);

  /*! Computation of stress tensor */
  void
  scheduleComputeStressTensor(Uintah::SchedulerP& sched,
                              const Uintah::PatchSet* patches,
                              const Uintah::MaterialSet* matls);
  void
  computeStressTensor(const Uintah::ProcessorGroup*,
                      const Uintah::PatchSubset* patches,
                      const Uintah::MaterialSubset* matls,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::DataWarehouse* new_dw);

  /*! Computation of internal force */
  /*! Computation of grid internal force */
  /*! Computation of bond internal force */
  /*! Computation of particle internal force */
  void
  scheduleComputeInternalForce(Uintah::SchedulerP& sched,
                               const Uintah::PatchSet* patches,
                               const Uintah::MaterialSet* matls);
  void
  computeGridInternalForce(const Uintah::ProcessorGroup*,
                           const Uintah::PatchSubset* patches,
                           const Uintah::MaterialSubset* matls,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::DataWarehouse* new_dw);
  void
  computeBondInternalForce(const Uintah::ProcessorGroup*,
                           const Uintah::PatchSubset* patches,
                           const Uintah::MaterialSubset* matls,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::DataWarehouse* new_dw);
  void
  computeParticleInternalForce(const Uintah::ProcessorGroup*,
                               const Uintah::PatchSubset* patches,
                               const Uintah::MaterialSubset* matls,
                               Uintah::DataWarehouse* old_dw,
                               Uintah::DataWarehouse* new_dw);

  /*! Integration of grid acceleration */
  void
  scheduleComputeAndIntegrateGridAcceleration(Uintah::SchedulerP& sched,
                                              const Uintah::PatchSet* patches,
                                              const Uintah::MaterialSet* matls);
  void
  computeAndIntegrateGridAcceleration(const Uintah::ProcessorGroup*,
                                      const Uintah::PatchSubset* patches,
                                      const Uintah::MaterialSubset* matls,
                                      Uintah::DataWarehouse* old_dw,
                                      Uintah::DataWarehouse* new_dw);

  /*! Integration of particle acceleration */
  void
  scheduleComputeAndIntegrateParticleAcceleration(
    Uintah::SchedulerP& sched,
    const Uintah::PatchSet* patches,
    const Uintah::MaterialSet* matls);
  void
  computeAndIntegrateParticleAcceleration(const Uintah::ProcessorGroup*,
                                          const Uintah::PatchSubset* patches,
                                          const Uintah::MaterialSubset* matls,
                                          Uintah::DataWarehouse* old_dw,
                                          Uintah::DataWarehouse* new_dw);

  /*! Project particle acceleration and velocity to grid */
  void
  scheduleProjectParticleAccelerationToGrid(Uintah::SchedulerP& sched,
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls);
  void
  projectParticleAccelerationToGrid(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset* matls,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);

  /*! Correct contact forces */
  void
  scheduleContactMomentumExchangeAfterIntegration(
    Uintah::SchedulerP& sched,
    const Uintah::PatchSet* patches,
    const Uintah::MaterialSet* matls);

  /*! Grid boundary condition set up */
  void
  scheduleSetGridBoundaryConditions(Uintah::SchedulerP& sched,
                                    const Uintah::PatchSet* patches,
                                    const Uintah::MaterialSet* matls);
  void
  setGridBoundaryConditions(const Uintah::ProcessorGroup*,
                            const Uintah::PatchSubset* patches,
                            const Uintah::MaterialSubset*,
                            Uintah::DataWarehouse* old_dw,
                            Uintah::DataWarehouse* new_dw);

  /*! Update particle velocities and displacements from grid data */
  void
  scheduleUpdateParticleKinematics(Uintah::SchedulerP& sched,
                                   const Uintah::PatchSet* patches,
                                   const Uintah::MaterialSet* matls);
  void
  updateParticleKinematics(const Uintah::ProcessorGroup*,
                           const Uintah::PatchSubset* patches,
                           const Uintah::MaterialSubset* matls,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::DataWarehouse* new_dw);

  /*! Compute damage and remove bonds */
  void
  scheduleComputeDamage(Uintah::SchedulerP& sched,
                        const Uintah::PatchSet* patches,
                        const Uintah::MaterialSet* matls);
  void
  computeDamage(const Uintah::ProcessorGroup*,
                const Uintah::PatchSubset* patches,
                const Uintah::MaterialSubset*,
                Uintah::DataWarehouse* old_dw,
                Uintah::DataWarehouse* new_dw);

  /*! Finalize the particle state to get ready for the next time step */
  void
  scheduleFinalizeParticleState(Uintah::SchedulerP& sched,
                                const Uintah::PatchSet* patches,
                                const Uintah::MaterialSet* matls);
  void
  finalizeParticleState(const Uintah::ProcessorGroup*,
                        const Uintah::PatchSubset* patches,
                        const Uintah::MaterialSubset* matls,
                        Uintah::DataWarehouse* old_dw,
                        Uintah::DataWarehouse* new_dw);

  /*! Need taskgraph recompile ? */
  bool
  needRecompile(double time, double dt, const Uintah::GridP& grid);

  template<typename T>
  void
  setParticleDefault(Uintah::ParticleVariable<T>& pvar,
                     const Uintah::VarLabel* label,
                     Uintah::ParticleSubset* pset,
                     Uintah::DataWarehouse* new_dw,
                     const T& val)
  {
    new_dw->allocateAndPut(pvar, label, pset);
    for (auto idx : *pset) {
      pvar[idx] = val;
    }
  }

protected:
  Uintah::MaterialManagerP d_mat_manager{ nullptr };

  PeridynamicsLabel* d_pd_labels{ nullptr };
  PeridynamicsFlags* d_pd_flags{ nullptr };
  PeridynamicsDefGradComputer* d_defGradComputer{ nullptr };
  BondInternalForceComputer* d_bondIntForceComputer{ nullptr };
  ParticleInternalForceComputer* d_intForceComputer{ nullptr };
  FamilyComputer* d_familyComputer{ nullptr };
  ContactModelBase* d_contactModel{ nullptr };

  Uintah::ParticleInterpolator* d_interpolator{ nullptr };

  int d_numGhostNodes;     // Number of ghost nodes needed
  int d_numGhostParticles; // Number of ghost particles needed
  bool d_recompile;
  Uintah::MaterialSubset* d_loadCurveIndex{ nullptr };

private:
};

} // end namespace Vaango

#endif
