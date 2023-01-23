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

#ifndef VAANGO_CCA_COMPONENTS_MPM_SERIALMPM_H
#define VAANGO_CCA_COMPONENTS_MPM_SERIALMPM_H

#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/SimulationCommon/SimulationCommon.h>
#include <Core/Grid/MaterialManagerP.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/SwitchingCriteria.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

// put here to avoid template problems
#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/GradientComputer/DeformationGradientComputer.h>
#include <CCA/Components/MPM/PhysicalBC/LoadCurve.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

namespace Uintah {

class ThermalContact;
class HeatConduction;
class AnalysisModule;
class SDInterfaceModel;
class SDInterfaceTasks;
class FluxBCModel;
class CZLabel;
class CohesiveZoneTasks;
class ScalarDiffusionTasks;

class SerialMPM
  : public SimulationCommon
  , public MPMCommon
{
public:
  std::unique_ptr<Contact> contactModel{ nullptr };
  std::unique_ptr<ThermalContact> thermalContactModel{ nullptr };
  std::unique_ptr<HeatConduction> heatConductionModel{ nullptr };
  std::unique_ptr<SDInterfaceTasks> sdInterfaceTasks{ nullptr };

public:
  SerialMPM(const ProcessorGroup* myworld, const MaterialManagerP matManager);

  virtual ~SerialMPM() noexcept(false);

  // No copy or move allowed
  SerialMPM(const SerialMPM&) = delete;
  SerialMPM(SerialMPM&&)      = delete;
  SerialMPM&
  operator=(const SerialMPM&) = delete;
  SerialMPM&
  operator=(SerialMPM&&) = delete;

  virtual double
  recomputeDelT(double delT);

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP&,
               MaterialManagerP&);

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

  virtual void
  scheduleDeleteGeometryObjects(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  restartInitialize();

  void
  schedulePrintParticleCount(const LevelP& level, SchedulerP& sched);

  void
  scheduleTotalParticleCount(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls);

  /*!
   * Schedule the initialization of the stress and deformation gradient
   * based on the body forces (which also have to be computed)
   */
  void
  scheduleInitializeStressAndDefGradFromBodyForce(const LevelP& level,
                                                  SchedulerP& sched);
  /*!
   * Actually initialize the body force acceleration
   */
  void
  initializeBodyForce(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset*,
                      DataWarehouse*,
                      DataWarehouse* new_dw);
  /*!
   * Actually initialize the stress and deformation gradient assuming linear
   * elastic behavior after computing the body force acceleration
   *
   * **WARNING** Assumes zero shear stresses and that body forces are aligned
   *             with coordinate directions
   */
  void
  initializeStressAndDefGradFromBodyForce(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset*,
                                          DataWarehouse*,
                                          DataWarehouse* new_dw);
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

  void
  setMPMLabel(MPMLabel* Mlb)
  {
    d_mpmLabels.reset(Mlb);
  };

  void
  setWithICE()
  {
    d_mpmFlags->d_withICE = true;
  };

public:
  enum IntegratorType
  {
    Explicit,
    Implicit,
    Fracture
  };

protected:
  friend class MPMICE;
  friend class SingleHydroMPM;

  virtual void
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
  totalParticleCount(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  //////////
  // Initialize particle data with a default values in the
  // new datawarehouse
  template<typename T>
  void
  setParticleDefault(ParticleVariable<T>& pvar,
                     const VarLabel* label,
                     ParticleSubset* pset,
                     DataWarehouse* new_dw,
                     const T& val);

  void
  printParticleLabels(vector<const VarLabel*> label,
                      DataWarehouse* dw,
                      int dwi,
                      const Patch* patch);

  void
  scheduleInitializePressureBCs(const LevelP& level, SchedulerP&);

  void
  scheduleInitializeMomentBCs(const LevelP& level, SchedulerP&);

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
  initializeMomentBC(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  scheduleRestartInitialize_alt(SchedulerP& sched, const LevelP& level);

  void
  restartInitialize_alt(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse*,
                        DataWarehouse* new_dw){};

  void
  restartInitializeTask(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset*,
                        DataWarehouse*,
                        DataWarehouse* new_dw);

  void
  actuallyComputeStableTimestep(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw);

  virtual void
  interpolateParticlesToGrid(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  virtual void
  computeNormals(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  findSurfaceParticles(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  computeLogisticRegression(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  virtual void
  computeSSPlusVp(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  virtual void
  computeSPlusSSPlusVp(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  virtual void
  addCohesiveZoneForces(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  void
  computeUnrotatedStressAndDeformationRate(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset*,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw);

  virtual void
  computeStressTensor(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  void
  computeRotatedStress(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  /*! Compute the velocity gradient ad deformation gradient */
  void
  computeDeformationGradient(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  /*! Compute the basic damage variables */
  void
  computeBasicDamage(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  /*! Update the erosion parameter if mass is to be removed */
  void
  updateErosionParameter(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  /*! Find particles that should be deleted */
  void
  findRogueParticles(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  //////////
  // Compute Accumulated Strain Energy
  void
  computeAccStrainEnergy(const ProcessorGroup*,
                         const PatchSubset*,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw);

  virtual void
  computeContactArea(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  virtual void
  computeInternalForce(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  virtual void
  computeAndIntegrateAcceleration(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  void
  setGridBoundaryConditions(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  /*!  Compute body forces for each particle.  Needed for rotating
    objects */
  void
  computeParticleBodyForce(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  //////////
  // This task is to be used for setting particle external force
  // and external heat rate.  I'm creating a separate task so that
  // user defined schemes for setting these can be implemented without
  // editing the core routines
  void
  applyExternalLoads(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  computeCurrentParticleSize(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  void
  addNewParticles(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  /*!  Convert the localized particles into particles of a new material
    with a different velocity field */
  void
  convertLocalizedParticles(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset* matls,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  virtual void
  interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  /*! Final particle update schedule and actual */
  virtual void
  scheduleFinalParticleUpdate(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  finalParticleUpdate(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  virtual void
  updateCohesiveZones(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  virtual void
  setPrescribedMotion(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  //////////
  // Allow blocks of particles to be moved according to a prescribed schedule:
  virtual void
  insertParticles(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);

  /*! Add new particles to the simulation based on criteria TBD: */
  virtual void
  scheduleAddParticles(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  addParticles(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  // Used to compute the particles initial physical size
  // for use in deformed particle visualization
  virtual void
  computeParticleScaleFactor(const ProcessorGroup*,
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

  virtual void
  scheduleInterpolateParticlesToGrid(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  virtual void
  scheduleComputeNormals(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls);

  virtual void
  scheduleFindSurfaceParticles(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls);

  virtual void
  scheduleComputeLogisticRegression(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls);

  virtual void
  scheduleAddCohesiveZoneForces(SchedulerP&,
                                const PatchSet*,
                                const MaterialSubset*,
                                const MaterialSubset*,
                                const MaterialSet*);

  virtual void
  scheduleComputeSSPlusVp(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeSPlusSSPlusVp(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  virtual void
  scheduleComputeHeatExchange(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleExMomInterpolated(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleUnrotateStressAndDeformationRate(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls);

  virtual void
  scheduleComputeStressTensor(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleRotateStress(SchedulerP& sched,
                       const PatchSet* patches,
                       const MaterialSet* matls);

  /*!----------------------------------------------------------------------
   * scheduleComputeXPICVelocities
   *-----------------------------------------------------------------------*/
  void
  scheduleComputeXPICVelocities(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls);
  /*!----------------------------------------------------------------------
   * computeParticleVelocityXPIC
   *-----------------------------------------------------------------------*/
  void
  computeParticleVelocityXPIC(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);
  /*!----------------------------------------------------------------------
   * computeGridVelocityXPIC
   *-----------------------------------------------------------------------*/
  void
  computeGridVelocityXPIC(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  void
  scheduleComputeDeformationGradient(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  scheduleComputeBasicDamage(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleUpdateErosionParameter(SchedulerP& sched,
                                 const PatchSet* patches,
                                 const MaterialSet* matls);

  void
  scheduleFindRogueParticles(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls);

  void
  scheduleComputeAccStrainEnergy(SchedulerP&,
                                 const PatchSet*,
                                 const MaterialSet*);

  virtual void
  scheduleComputeContactArea(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeInternalForce(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  virtual void
  scheduleComputeInternalHeatRate(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  virtual void
  scheduleComputeNodalHeatFlux(SchedulerP&,
                               const PatchSet*,
                               const MaterialSet*);

  virtual void
  scheduleSolveHeatEquations(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeAndIntegrateAcceleration(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

  virtual void
  scheduleIntegrateTemperatureRate(SchedulerP&,
                                   const PatchSet*,
                                   const MaterialSet*);

  virtual void
  scheduleExMomIntegrated(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleSetGridBoundaryConditions(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet* matls);

  void
  scheduleComputeParticleBodyForce(SchedulerP&,
                                   const PatchSet*,
                                   const MaterialSet*);

  void
  scheduleApplyExternalLoads(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleInterpolateToParticlesAndUpdate(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

  virtual void
  scheduleUpdateCohesiveZones(SchedulerP&,
                              const PatchSet*,
                              const MaterialSubset*,
                              const MaterialSubset*,
                              const MaterialSet*);

  virtual void
  scheduleSetPrescribedMotion(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleAddNewParticles(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleConvertLocalizedParticles(SchedulerP&,
                                    const PatchSet*,
                                    const MaterialSet*);

  void
  scheduleCheckNeedAddMPMMaterial(SchedulerP&,
                                  const PatchSet* patches,
                                  const MaterialSet*);

  virtual void
  scheduleInterpolateToParticlesAndUpdateMom1(SchedulerP&,
                                              const PatchSet*,
                                              const MaterialSet*);

  virtual void
  scheduleInterpolateParticleVelToGridMom(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*);

  virtual void
  scheduleInterpolateToParticlesAndUpdateMom2(SchedulerP&,
                                              const PatchSet*,
                                              const MaterialSet*);

  virtual void
  scheduleInsertParticles(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeParticleScaleFactor(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  virtual void
  interpolateToParticlesAndUpdateMom1(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* matls,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw);

  virtual void
  interpolateToParticlesAndUpdateMom2(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* matls,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw);

  virtual void
  interpolateParticleVelToGridMom(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw);

  void
  checkNeedAddMPMMaterial(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw);

  bool
  needRecompile(double time, double dt, const GridP& grid);

  void
  readPrescribedDeformations(string filename);

  void
  readInsertParticlesFile(string filename);

  virtual void
  scheduleSwitchTest(const LevelP& level, SchedulerP& sched);

protected:
  double d_nextOutputTime{ 0.0 };
  double d_SMALL_NUM_MPM{ 1.0e-200 };
  int d_numGhostParticles{ 1 }; // Number of ghost particles needed.
  int d_numGhostNodes{ 1 };     // Number of ghost nodes     needed.
  bool d_fracture{false};
  bool d_recompile{false};

  IntegratorType d_integrator;
  MaterialSubset* d_loadCurveIndex{nullptr};

  std::unique_ptr<MPMLabel> d_mpmLabels{nullptr};
  std::unique_ptr<MPMFlags> d_mpmFlags{nullptr};

  std::unique_ptr<DeformationGradientComputer> d_defGradComputer{nullptr};
  std::unique_ptr<ScalarDiffusionTasks> d_diffusionTasks{nullptr};
  std::unique_ptr<SwitchingCriteria> d_switchCriteria{nullptr};
  std::vector<std::unique_ptr<AnalysisModule>> d_analysisModules;

  MaterialManagerP d_materialManager{nullptr};
  Output* d_dataArchiver{nullptr};

  std::list<Patch::FaceType>
    d_boundaryTractionFaces; // list of xminus, xplus, yminus, ...
  std::vector<MPMPhysicalBC*> d_physicalBCs;

  std::vector<double> d_prescribedTimes; // These three are used only if
  std::vector<double> d_prescribedAngle; // d_prescribeDeformation
  std::vector<Vector>
    d_prescribedRotationAxis; // is "true".  It is "false" by default.
  std::vector<Matrix3> d_prescribedF;

  // The following are used iff the d_insertParticles flag is true.
  std::vector<double> d_IPTimes;
  std::vector<double> d_IPColor;
  std::vector<Vector> d_IPTranslate;
  std::vector<Vector> d_IPVelNew;

};

} // end namespace Uintah

#endif
