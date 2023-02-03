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

#ifndef UINTAH_HOMEBREW_MPMICE_H
#define UINTAH_HOMEBREW_MPMICE_H

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/MPM/Core/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBC.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <CCA/Ports/SimulationInterface.h>
#include <CCA/Ports/SwitchingCriteria.h>

#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class ICE;
class ICELabel;
class MPMLabel;
class MPMICELabel;
class Output;

enum class MPMType
{
  STAND_MPMICE = 0,
  RIGID_MPMICE,
  SHELL_MPMICE,
  FRACTURE_MPMICE
};

class MPMICE final : public SimulationCommon
{
public:
  MPMICE(const ProcessorGroup* myworld,
         const MaterialManagerP& mat_manager,
         const MPMType& type,
         bool doAMR = false);

  virtual ~MPMICE();

  MPMICE(const MPMICE&) = delete;
  MPMICE(MPMICE&&)      = delete;
  MPMICE&
  operator=(const MPMICE&) = delete;
  MPMICE&
  operator=(MPMICE&&) = delete;

  double
  recomputeDelT(double delT) override;

  void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid) override;

  void
  outputProblemSpec(ProblemSpecP& ps) override;

  void
  scheduleInitialize(const LevelP& level, SchedulerP&) override;

  void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& sched) override;

  void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP&) override;

  // scheduleTimeAdvance version called by the AMR simulation controller.
  void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&) override;

  void
  scheduleFinalizeTimestep(const LevelP& level, SchedulerP&) override;

  void
  scheduleSwitchTest(const LevelP& level, SchedulerP& sched) override;

  //__________________________________
  //   AMR
  void
  scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                               SchedulerP& sched) override;

  void
  scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched) override;

  void
  scheduleRefineInterface(const LevelP& fineLevel,
                          SchedulerP& scheduler,
                          bool needOld,
                          bool needNew) override;

  void
  scheduleRefine(const PatchSet* patches, SchedulerP& sched) override;

  void
  scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched) override;

private:
  void
  scheduleInterpolateNCToCC_0(SchedulerP&,
                              const PatchSet*,
                              const MaterialSubset*,
                              const MaterialSet*);

  void
  scheduleCoarsenCC_0(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleCoarsenNCMass(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleComputeLagrangianValuesMPM(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSubset*,
                                     const MaterialSet*);

  void
  scheduleCoarsenLagrangianValuesMPM(SchedulerP&,
                                     const PatchSet*,
                                     const MaterialSet*);

  void
  scheduleInterpolateCCToNC(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleComputeCCVelAndTempRates(SchedulerP&,
                                   const PatchSet*,
                                   const MaterialSet*);

  void
  scheduleRefineCC(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleComputePressure(SchedulerP&,
                          const PatchSet*,
                          const MaterialSubset*,
                          const MaterialSubset*,
                          const MaterialSubset*,
                          const MaterialSet*);

  void
  scheduleInterpolatePressCCToPressNC(SchedulerP&,
                                      const PatchSet*,
                                      const MaterialSubset*,
                                      const MaterialSet*);

  void
  scheduleRefinePressCC(SchedulerP&,
                        const PatchSet*,
                        const MaterialSubset*,
                        const MaterialSet*);

  void
  scheduleInterpolatePAndGradP(SchedulerP&,
                               const PatchSet*,
                               const MaterialSubset*,
                               const MaterialSubset*,
                               const MaterialSubset*,
                               const MaterialSet*);

  // AMR
  void
  scheduleRefineVariableCC(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSet* matls,
                           const VarLabel* variable);

  template<typename T>
  void
  scheduleCoarsenVariableCC(SchedulerP& sched,
                            const PatchSet* patches,
                            const MaterialSet* matls,
                            const VarLabel* variable,
                            T defaultValue,
                            bool modifies,
                            const string& coarsenMethod);

  template<typename T>
  void
  scheduleCoarsenVariableNC(SchedulerP& sched,
                            const PatchSet* patches,
                            const MaterialSet* matls,
                            const VarLabel* variable,
                            T defaultValue,
                            bool modifies,
                            string coarsenMethod);

  void
  scheduleParticleRelocation(SchedulerP& sched,
                             const LevelP& level,
                             const MaterialSet* mpm_matls);

  //______________________________________________________________________
  //       A C T U A L   S T E P S :
  void
  actuallyInitialize(const ProcessorGroup*,
                     const PatchSubset* patch,
                     const MaterialSubset* matls,
                     DataWarehouse*,
                     DataWarehouse* new_dw);

  void
  interpolateNCToCC_0(const ProcessorGroup*,
                      const PatchSubset* patch,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  void
  computeLagrangianValuesMPM(const ProcessorGroup*,
                             const PatchSubset* patch,
                             const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  void
  computeEquilibrationPressure(const ProcessorGroup*,
                               const PatchSubset* patch,
                               const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse*,
                               const MaterialSubset* press_matl);

  void
  interpolateCCToNC(const ProcessorGroup*,
                    const PatchSubset* patch,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  void
  computeCCVelAndTempRates(const ProcessorGroup*,
                           const PatchSubset* patch,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  void
  interpolatePressCCToPressNC(const ProcessorGroup*,
                              const PatchSubset* patch,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  void
  interpolatePAndGradP(const ProcessorGroup*,
                       const PatchSubset* patch,
                       const MaterialSubset* matls,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw);

  void
  binaryPressureSearch(std::vector<constCCVariable<double>>& Temp,
                       std::vector<CCVariable<double>>& rho_micro,
                       std::vector<CCVariable<double>>& vol_frac,
                       std::vector<CCVariable<double>>& rho_CC_new,
                       std::vector<CCVariable<double>>& speedSound_new,
                       std::vector<double>& dp_drho,
                       std::vector<double>& dp_de,
                       std::vector<double>& press_eos,
                       constCCVariable<double>& press,
                       CCVariable<double>& press_new,
                       double press_ref,
                       std::vector<constCCVariable<double>>& cv,
                       std::vector<constCCVariable<double>>& gamma,
                       double convergence_crit,
                       int numALLMatls,
                       int& count,
                       double& sum,
                       IntVector c);

  // AMR
  void
  refine(const ProcessorGroup*,
         const PatchSubset* patches,
         const MaterialSubset* matls,
         DataWarehouse*,
         DataWarehouse* new_dw);

  template<typename T>
  void
  refineVariableCC(const ProcessorGroup*,
                   const PatchSubset* patch,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* variable);

  template<typename T>
  void
  coarsenDriver_stdNC(IntVector cl,
                      IntVector ch,
                      IntVector refinementRatio,
                      double ratio,
                      const Level* coarseLevel,
                      constNCVariable<T>& fine_q_NC,
                      NCVariable<T>& coarse_q_NC);

  template<typename T>
  void
  coarsenVariableCC(const ProcessorGroup*,
                    const PatchSubset* patch,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw,
                    const VarLabel* variable,
                    T defaultValue,
                    bool modifies,
                    string coarsenMethod);

  template<typename T>
  void
  coarsenVariableNC(const ProcessorGroup*,
                    const PatchSubset* patch,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw,
                    const VarLabel* variable,
                    T defaultValue,
                    bool modifies,
                    string coarsenMethod);

  void
  refineCoarseFineInterface(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* fine_old_dw,
                            DataWarehouse* fine_new_dw);

private:
  enum bctype
  {
    NONE = 0,
    FIXED,
    SYMMETRY,
    NEIGHBOR
  };

protected:
  std::unique_ptr<MPMICELabel> d_mpmice_labels{ nullptr };

  std::unique_ptr<SerialMPM> d_mpm{ nullptr };
  std::unique_ptr<ICE> d_ice{ nullptr };

  const MPMLabel* d_mpm_labels{ nullptr };
  const ICELabel* d_ice_labels{ nullptr };

  std::vector<std::unique_ptr<AnalysisModule>> d_analysisModules;

  SwitchingCriteria* d_switchCriteria{ nullptr };

  bool d_rigidMPM{ false };
  bool d_testForNegTemps_mpm{ true };
  bool do_mlmpmice{ false };
  bool d_useSimpleEquilibrationPressure{ false };

  int d_8or27{ 8 };
  int d_num_ghost_nodes{ 1 };

  double d_SMALL_NUM{ 1.0e-100 };
  double d_TINY_RHO{ 1.0e-12 };
  double d_convergence_tolerance{ 1.0e-16 };
};

} // End namespace Uintah

#endif
