/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef UINTAH_HOMEBREW_IMP_MPM_H
#define UINTAH_HOMEBREW_IMP_MPM_H

#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <CCA/Components/MPM/Core/ImpMPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/GradientComputer/DeformationGradientComputer.h>
#include <CCA/Components/MPM/ImpMPMSolvers/Solver.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/SwitchingCriteria.h>

#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <sci_defs/petsc_defs.h>

#include <list>
#include <map>
#include <vector>

namespace Uintah {

class DataWarehouse;
class MPMLabel;
class ImpMPMLabel;
class ProcessorGroup;
class VarLabel;
class Task;
class ImplicitHeatConductionTasks;

class ParticleTempShape
{
public:
  double pTemperature;
  std::vector<IntVector> cellNodes;
  std::vector<double> shapeFnValues;
};

class ImpMPM
  : public SimulationCommon
  , public MPMCommon
{
public:
  ImpMPM(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);
  virtual ~ImpMPM();

  ImpMPM(const ImpMPM&) = delete;
  ImpMPM(ImpMPM&&)      = delete;
  ImpMPM&
  operator=(const ImpMPM&) = delete;
  ImpMPM&
  operator=(ImpMPM&&) = delete;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& mat_ps,
               GridP& grid);

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleSwitchInitialization(const LevelP& level, SchedulerP&);

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRefine(const PatchSet* patches, SchedulerP& scheduler);

  virtual void
  scheduleRefineInterface(const LevelP& fineLevel,
                          SchedulerP& scheduler,
                          bool needCoarse,
                          bool needFine);

  virtual void
  scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched);

  /// Schedule to mark flags for AMR regridding
  virtual void
  scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched);

  /// Schedule to mark initial flags for AMR regridding
  void
  scheduleInitialErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched);

  virtual double
  recomputeDelT(double dt) {
    return dt * d_mpm_flags->d_delTDecreaseFactor;
  }

  void
  scheduleSwitchTest(const LevelP& level, SchedulerP& sched);

  enum IntegratorType
  {
    Explicit,
    Implicit
  };

private:
  void
  scheduleInitializePressureBCs(const LevelP& level, SchedulerP&);

  void
  scheduleComputeDeformationGradient(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls,
                                     const bool recursion);

  void
  scheduleComputeStressTensor(SchedulerP&,
                              const PatchSet*,
                              const MaterialSet*,
                              const bool recursion);

  void
  scheduleComputeDeformationGradient(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls);

  void
  scheduleComputeStressTensor(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleFormStiffnessMatrix(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleComputeInternalForce(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  void
  scheduleFormQ(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleUpdateGridKinematics(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  void
  scheduleComputeParticleBodyForce(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls);

  void
  scheduleApplyExternalLoads(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleInterpolateParticlesToGrid(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSubset*,
                                     const MaterialSet*);

  void
  scheduleFindSurfaceParticles(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls);

  void
  scheduleDestroyMatrix(SchedulerP&,
                        const PatchSet*,
                        const MaterialSet*,
                        const bool recursion);

  void
  scheduleCreateMatrix(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleApplyBoundaryConditions(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  scheduleComputeContact(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleFindFixedDOF(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleSolveForDuCG(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleGetDisplacementIncrement(SchedulerP&,
                                   const PatchSet*,
                                   const MaterialSet*);

  void
  scheduleUpdateTotalDisplacement(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  scheduleComputeAcceleration(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleInterpolateToParticlesAndUpdate(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

  void
  scheduleInterpolateStressToGrid(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  scheduleIterate(SchedulerP&,
                  const LevelP&,
                  const PatchSet*,
                  const MaterialSet*);

  void
  scheduleCheckConvergence(SchedulerP&,
                           const LevelP&,
                           const PatchSet*,
                           const MaterialSet*);

private:
  friend class MPMICE;

  inline bool
  compare(double num1, double num2)
  {
    double EPSILON = 1.e-16;
    return (std::abs(num1 - num2) <= EPSILON);
  };

  void
  actuallyInitialize(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  printParticleCount(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  void
  initializePressureBC(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  initializeHeatFluxBC(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  actuallyComputeStableTimestep(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

  void
  computeParticleBodyForce(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  void
  applyExternalLoads(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  interpolateParticlesToGrid(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  void
  findSurfaceParticles(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  projectCCHeatSourceToNodes(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  void
  computeCCVolume(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  void
  rigidBody(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset* matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

  void
  destroyMatrix(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                DataWarehouse* old_dw,
                DataWarehouse* new_dw,
                bool recursion);

  void
  createMatrix(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  applyBoundaryConditions(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  void
  computeContact(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset* matls,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  findFixedDOF(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  computeDeformationGradient(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             bool recursion);

  void
  computeDeformationGradient(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  // This is for the computation with the 24 x 24 matrix
  void
  computeStressTensorImplicit(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw,
                              bool recursion);

  // No matrix calculations are performed.
  void
  computeStressTensorImplicit(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  void
  formStiffnessMatrix(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);
  //////////
  // Insert Documentation Here:
  void
  computeInternalForce(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  iterate(const ProcessorGroup* pg,
          const PatchSubset* patches,
          const MaterialSubset* matls,
          DataWarehouse* old_dw,
          DataWarehouse* new_dw,
          LevelP level,
          Scheduler* sched);

  void
  formQ(const ProcessorGroup*,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse* old_dw,
        DataWarehouse* new_dw);

  void
  solveForDuCG(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  solveForTemp(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  getDisplacementIncrement(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  void
  getTemperatureIncrement(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  void
  updateGridKinematics(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  checkConvergence(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  updateTotalDisplacement(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  computeAcceleration(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  void
  refine(const ProcessorGroup*,
         const PatchSubset* patches,
         const MaterialSubset* matls,
         DataWarehouse*,
         DataWarehouse* new_dw);

  void
  errorEstimate(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                DataWarehouse*,
                DataWarehouse* new_dw);

  void
  initialErrorEstimate(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse*,
                       DataWarehouse* new_dw);

  //////////
  // Insert Documentation Here:
  void
  interpolateStressToGrid(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

private:
  std::unique_ptr<MPMLabel> d_mpm_labels{ nullptr };
  std::unique_ptr<ImpMPMLabel> d_impmpmLabels{ nullptr };
  std::unique_ptr<ImpMPMFlags> d_mpm_flags{ nullptr };
  std::unique_ptr<DeformationGradientComputer> d_defGradComputer{ nullptr };
  std::unique_ptr<ImplicitHeatConductionTasks> d_heatConductionTasks{ nullptr };
  std::unique_ptr<Solver> d_solver{ nullptr };

  std::shared_ptr<MaterialSubset> d_oneMaterial{ nullptr };
  std::shared_ptr<MaterialSubset> d_loadCurveIndex{ nullptr };

  const PatchSet* d_perprocPatches{ nullptr };
  SwitchingCriteria* d_switchCriteria{ nullptr };

  // stuff for not having to recompile the iterative scheduler every timstep
  SchedulerP d_subsched;
  bool d_recompileSubsched{ true };

  int d_numGhostParticles{ 1 }; // Number of ghost particles needed.
  int d_numGhostNodes{ 1 };     // Number of ghost nodes     needed.
  int d_numIterations{ 0 };

  double d_nextOutputTime{ 0.0 };
  double d_SMALL_NUM_MPM{ 1.0e-200 };
  double d_initialDt{ 10000.0 };

  Vector d_contactDirections{ 1.0, 0.0, 0.0 }; // For rigid body contact
  std::string d_contactType{ "null" };
  bool d_rigidContact{ false };
  bool d_singleVelocityContact{ false };
  double d_contactStopTime{ 0.0 };                    // for rigid contact
  Vector d_velocityAfterContactStop{ 0.0, 0.0, 0.0 }; // for rigid contact

  // list of xminus, xplus, ...
  std::list<Patch::FaceType> d_boundaryTractionFaces;
};

} // end namespace Uintah

#endif
