/*
 * The MIT License
 *
 * Copyright (c) 2018-2023 Biswajit Banerjee
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

#ifndef VAANGO_CCA_COMPONENTS_MPM_UofU_MPM_H
#define VAANGO_CCA_COMPONENTS_MPM_UofU_MPM_H

#include <CCA/Components/MPM/SerialMPM.h>

namespace Uintah {

class UofU_MPM
  : public SerialMPM
{
public:
  inline static const Matrix3 Identity{ 1.0, 0.0, 0.0, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0 };
  inline static const Matrix3 Zero{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  };

  UofU_MPM(const ProcessorGroup* myworld, const MaterialManagerP& matManager);

  UofU_MPM(const UofU_MPM&) = delete;
  UofU_MPM(UofU_MPM&&)      = delete;
  UofU_MPM&
  operator=(const UofU_MPM&) = delete;
  UofU_MPM&
  operator=(UofU_MPM&&) = delete;

  virtual ~UofU_MPM();

  virtual double
  recomputeDelT(double delT)
  {
    return delT * 0.1;
  }

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid);

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

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

  void
  setMPMLabel(MPMLabel* labels)
  {
    d_mpm_labels.reset(labels);
  };

  enum IntegratorType
  {
    Explicit,
    Implicit,
    Fracture
  };

protected:
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

  void
  checkGridVelocity(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  void
  computeVelocityAndDeformationGradient(const ProcessorGroup*,
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

  void
  computeAcceleration(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset*,
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

  virtual void
  setPrescribedMotion(const ProcessorGroup*,
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

  virtual void
  scheduleInterpolateParticlesToGrid(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  scheduleCheckGridVelocity(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleContactMomentumExchange(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls,
                                  const VarLabel* label);

  void
  scheduleComputeVelocityAndDeformationGradient(SchedulerP&,
                                                const PatchSet*,
                                                const MaterialSet*);

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

  void
  scheduleComputeAcceleration(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* matls);

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
  scheduleSetPrescribedMotion(SchedulerP&, const PatchSet*, const MaterialSet*);

  virtual void
  scheduleComputeParticleScaleFactor(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  readPrescribedDeformations(string filename);

  void
  readInsertParticlesFile(string filename);

};

} // end namespace Uintah

#endif
