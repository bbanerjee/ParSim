/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef VAANGO_CCA_COMPONENTS_ICE_ICE_H
#define VAANGO_CCA_COMPONENTS_ICE_ICE_H

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <CCA/Components/ICE/Advection/Advector.h>
#include <CCA/Components/ICE/Core/BoundaryCond.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/ICE/CustomBCs/LODI2.h>
#include <CCA/Components/ICE/TurbulenceModel/Turbulence.h>
#include <CCA/Components/ICE/customInitialize.h>

#include <CCA/Components/Models/MultiMatlExchange/ExchangeCoefficients.h>
#include <CCA/Components/Models/MultiMatlExchange/ExchangeModel.h>

#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>

#include <CCA/Ports/ModelInterface.h>

#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Grid/Variables/Utils.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/UintahMiscMath.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <sci_defs/hypre_defs.h>

#include <string>
#include <vector>

#define MAX_MATLS 16

using namespace Uintah::ExchangeModels;

namespace Uintah {

class ModelInterface;
class Turbulence;
class WallSheatStress;
class AnalysisModule;

// The following two structs are used by computeEquilibrationPressure to store
// debug information:
//
struct EqPress_dbgMatl
{
  int mat;
  double press_eos;
  double volFrac;
  double rhoMicro;
  double temp_CC;
  double rho_CC;
};

struct EqPress_dbg
{
  int count;
  double sumVolFrac;
  double press_new;
  double delPress;
  std::vector<EqPress_dbgMatl> matl;
};

class ICE : public SimulationCommon
{
public:
  ICE(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager);

  virtual ~ICE();

  ICE(const ICE&) = delete;
  ICE(ICE&&)      = delete;
  ICE&
  operator=(const ICE&) = delete;
  ICE&
  operator=(ICE&&) = delete;

  virtual double
  recomputeDelT(double delT);

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "");

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& sched);

  virtual void
  scheduleComputeStableTimestep(const LevelP&, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual void
  scheduleFinalizeTimestep(const LevelP& level, SchedulerP&);

  virtual void
  scheduleAnalysis(const LevelP& level, SchedulerP&);

  void
  scheduleComputePressure(SchedulerP&,
                          const PatchSet*,
                          const MaterialSubset*,
                          const MaterialSet*);

  void
  scheduleComputeTempFC(SchedulerP&,
                        const PatchSet*,
                        const MaterialSubset*,
                        const MaterialSubset*,
                        const MaterialSet*);

  void
  scheduleComputeVel_FC(SchedulerP&,
                        const PatchSet*,
                        const MaterialSubset*,
                        const MaterialSubset*,
                        const MaterialSubset*,
                        const MaterialSet*);

  void
  scheduleComputeDelPressAndUpdatePressCC(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSubset*,
                                          const MaterialSubset*,
                                          const MaterialSubset*,
                                          const MaterialSet*);

  void
  scheduleComputePressFC(SchedulerP&,
                         const PatchSet*,
                         const MaterialSubset*,
                         const MaterialSet*);

  void
  scheduleComputeThermoTransportProperties(SchedulerP&,
                                           const LevelP& level,
                                           const MaterialSet*);

  void
  scheduleVelTau_CC(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleViscousShearStress(SchedulerP&, const PatchSet*, const MaterialSet*);

  void
  scheduleAccumulateMomentumSourceSinks(SchedulerP&,
                                        const PatchSet*,
                                        const MaterialSubset*,
                                        const MaterialSubset*,
                                        const MaterialSubset*,
                                        const MaterialSet*);

  void
  scheduleAccumulateEnergySourceSinks(SchedulerP&,
                                      const PatchSet*,
                                      const MaterialSubset*,
                                      const MaterialSubset*,
                                      const MaterialSubset*,
                                      const MaterialSet*);

  void
  scheduleComputeLagrangianValues(SchedulerP&,
                                  const PatchSet*,
                                  const MaterialSet*);

  void
  scheduleComputeLagrangianSpecificVolume(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSubset*,
                                          const MaterialSubset*,
                                          const MaterialSubset*,
                                          const MaterialSet*);

  void
  scheduleMaxMach_on_Lodi_BC_Faces(SchedulerP&,
                                   const LevelP&,
                                   const MaterialSet*);

  void
  computesRequires_AMR_Refluxing(Task* t, const MaterialSet* ice_matls);

  void
  scheduleAdvectAndAdvanceInTime(SchedulerP&,
                                 const PatchSet*,
                                 const MaterialSubset*,
                                 const MaterialSet*);

  void
  scheduleConservedtoPrimitive_Vars(SchedulerP& sched,
                                    const PatchSet* patch_set,
                                    const MaterialSubset* ice_matlsub,
                                    const MaterialSet* ice_matls,
                                    const string& where);

  void
  scheduleTestConservation(SchedulerP&,
                           const PatchSet*,
                           const MaterialSubset*,
                           const MaterialSet*);

  void
  scheduleComputeTaskGraphIndex(SchedulerP& sched, const LevelP& level);

  //__________________________________
  //__________________________________
  //  I M P L I C I T   I C E

  void
  scheduleSetupMatrix(SchedulerP&,
                      const LevelP&,
                      const PatchSet*,
                      const MaterialSubset*,
                      const MaterialSet*);

  void
  scheduleSetupRHS(SchedulerP&,
                   const PatchSet*,
                   const MaterialSubset*,
                   const MaterialSet*,
                   bool insideOuterIterLoop,
                   const string& computes_or_modifies);

  void
  scheduleCompute_maxRHS(SchedulerP& sched,
                         const LevelP& level,
                         const MaterialSubset* one_matl,
                         const MaterialSet*);

  void
  scheduleUpdatePressure(SchedulerP&,
                         const LevelP&,
                         const PatchSet*,
                         const MaterialSubset*,
                         const MaterialSubset*,
                         const MaterialSubset*,
                         const MaterialSet*);

  void
  scheduleRecomputeVel_FC(SchedulerP& sched,
                          const PatchSet* patches,
                          const MaterialSubset*,
                          const MaterialSubset*,
                          const MaterialSubset*,
                          const MaterialSet*,
                          bool);

  void
  scheduleComputeDel_P(SchedulerP& sched,
                       const LevelP& level,
                       const PatchSet* patches,
                       const MaterialSubset* one_matl,
                       const MaterialSubset* press_matl,
                       const MaterialSet* all_matls);

  void
  scheduleImplicitPressureSolve(SchedulerP& sched,
                                const LevelP& level,
                                const PatchSet*,
                                const MaterialSubset* one_matl,
                                const MaterialSubset* press_matl,
                                const MaterialSubset* ice_matls,
                                const MaterialSubset* mpm_matls,
                                const MaterialSet* all_matls);

  //__________________________________
  //  I M P L I C I T   A M R I C E
  void
  scheduleCoarsen_delP(SchedulerP& sched,
                       const LevelP& level,
                       const MaterialSubset* press_matl,
                       const VarLabel* variable);

  void
  schedule_matrixBC_CFI_coarsePatch(SchedulerP& sched,
                                    const LevelP& coarseLevel,
                                    const MaterialSubset* one_matl,
                                    const MaterialSet* all_matls);

  void
  scheduleMultiLevelPressureSolve(SchedulerP& sched,
                                  const GridP grid,
                                  const PatchSet*,
                                  const MaterialSubset* one_matl,
                                  const MaterialSubset* press_matl,
                                  const MaterialSubset* ice_matls,
                                  const MaterialSubset* mpm_matls,
                                  const MaterialSet* all_matls);

  void
  scheduleZeroMatrix_UnderFinePatches(SchedulerP& sched,
                                      const LevelP& coarseLevel,
                                      const MaterialSubset* one_matl);

  void
  zeroMatrix_UnderFinePatches(const ProcessorGroup*,
                              const PatchSubset* coarsePatches,
                              const MaterialSubset*,
                              DataWarehouse*,
                              DataWarehouse* new_dw);

  void
  schedule_bogus_imp_delP(SchedulerP& sched,
                          const PatchSet* perProcPatches,
                          const MaterialSubset* press_matl,
                          const MaterialSet* all_matls);

  void
  scheduleAddReflux_RHS(SchedulerP& sched,
                        const LevelP& coarseLevel,
                        const MaterialSubset* one_matl,
                        const MaterialSet* all_matls,
                        const bool OnOff);

  //__________________________________
  //   M O D E L S
  void
  scheduleComputeModelSources(SchedulerP& sched,
                              const LevelP& level,
                              const MaterialSet* matls);

  void
  scheduleUpdateVolumeFraction(SchedulerP& sched,
                               const LevelP& level,
                               const MaterialSubset* press_matl,
                               const MaterialSet* matls);

  void
  scheduleComputeLagrangian_Transported_Vars(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet*);

  void
  setWithMPM()
  {
    d_with_mpm = true;
  };

  void
  setWithRigidMPM()
  {
    d_with_rigid_mpm = true;
  };

public:
  void
  actuallyInitialize(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* matls,
                     DataWarehouse*,
                     DataWarehouse* new_dw);

  void
  initializeSubTask_hydrostaticAdj(const ProcessorGroup*,
                                   const PatchSubset*,
                                   const MaterialSubset*,
                                   DataWarehouse*,
                                   DataWarehouse* new_dw);

  void
  actuallyComputeStableTimestep(const ProcessorGroup*,
                                const PatchSubset* patch,
                                const MaterialSubset* matls,
                                DataWarehouse*,
                                DataWarehouse*);

  void
  computeEquilibrationPressure(const ProcessorGroup*,
                               const PatchSubset* patch,
                               const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse*);

  void
  computeEquilPressure_1_matl(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

  void
  computeVel_FC(const ProcessorGroup*,
                const PatchSubset*,
                const MaterialSubset*,
                DataWarehouse*,
                DataWarehouse*);

  void
  updateVel_FC(const ProcessorGroup*,
               const PatchSubset*,
               const MaterialSubset*,
               DataWarehouse*,
               DataWarehouse*,
               bool);

  void
  computeTempFC(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset*,
                DataWarehouse*,
                DataWarehouse*);

  template<class T>
  void
  computeTempFace(CellIterator it,
                  IntVector adj_offset,
                  constCCVariable<double>& rho_CC,
                  constCCVariable<double>& Temp_CC,
                  T& Temp_FC);

  template<class T>
  void
  computeVelFace(int dir,
                 CellIterator it,
                 IntVector adj_offset,
                 double dx,
                 double delT,
                 double gravity,
                 constCCVariable<double>& rho_CC,
                 constCCVariable<double>& sp_vol_CC,
                 constCCVariable<Vector>& vel_CC,
                 constCCVariable<double>& press_CC,
                 T& vel_FC,
                 T& gradP_FC,
                 bool include_acc);

  template<class T>
  void
  updateVelFace(int dir,
                CellIterator it,
                IntVector adj_offset,
                double dx,
                double delT,
                constCCVariable<double>& sp_vol_CC,
                constCCVariable<double>& press_CC,
                T& vel_FC,
                T& grad_dp_FC);

  void
  computeDelPressAndUpdatePressCC(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse*,
                                  DataWarehouse*);

  void
  computePressFC(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset* matls,
                 DataWarehouse*,
                 DataWarehouse*);

  template<class T>
  void
  computePressFace(CellIterator it,
                   IntVector adj_offset,
                   constCCVariable<double>& sum_rho,
                   constCCVariable<double>& press_CC,
                   T& press_FC);

  void
  computeThermoTransportProperties(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ice_matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw);

  void
  VelTau_CC(const ProcessorGroup*,
            const PatchSubset* patches,
            const MaterialSubset* ice_matls,
            DataWarehouse* old_dw,
            DataWarehouse* new_dw);

  void
  computeVelTau_CCFace(const Patch* patch,
                       const Patch::FaceType face,
                       constCCVariable<Vector>& vel_CC,
                       CCVariable<Vector>& velTau_CC);

  void
  viscousShearStress(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* ice_matls,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw);

  void
  accumulateMomentumSourceSinks(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse*,
                                DataWarehouse*);

  void
  accumulateEnergySourceSinks(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse*,
                              DataWarehouse*);

  void
  computeLagrangianValues(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*,
                          DataWarehouse*);

  void
  computeLagrangianSpecificVolume(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse*,
                                  DataWarehouse*);

  void
  addExchangeToMomentumAndEnergy(const ProcessorGroup*,
                                 const PatchSubset*,
                                 const MaterialSubset*,
                                 DataWarehouse*,
                                 DataWarehouse*);

  void
  maxMach_on_Lodi_BC_Faces(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw);

  void
  advectAndAdvanceInTime(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* matls,
                         DataWarehouse*,
                         DataWarehouse*);

  void
  conservedtoPrimitive_Vars(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw);

  void
  computeTaskGraphIndex(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset*,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  //__________________________________
  //  I M P L I C I T   I C E
  void
  setupMatrix(const ProcessorGroup*,
              const PatchSubset* patches,
              const MaterialSubset*,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw);

  void
  setupRHS(const ProcessorGroup*,
           const PatchSubset* patches,
           const MaterialSubset*,
           DataWarehouse* old_dw,
           DataWarehouse* new_dw,
           bool insideOuterIterLoop,
           string computes_or_modifies);

  void
  compute_maxRHS(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse*,
                 DataWarehouse* new_dw);

  void
  updatePressure(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw);

  void
  computeDel_P(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset*,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw);

  void
  implicitPressureSolve(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset*,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw,
                        LevelP level,
                        const MaterialSubset*,
                        const MaterialSubset*);

  //__________________________________
  //  I M P L I C I T   A M R I C E
  void
  coarsen_delP(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse*,
               DataWarehouse* new_dw,
               const VarLabel* variable);

  void
  matrixCoarseLevelIterator(Patch::FaceType patchFace,
                            const Patch* coarsePatch,
                            const Patch* finePatch,
                            const Level* fineLevel,
                            CellIterator& iter,
                            bool& isRight_CP_FP_pair);

  void
  matrixBC_CFI_coarsePatch(const ProcessorGroup*,
                           const PatchSubset* coarsePatches,
                           const MaterialSubset* matls,
                           DataWarehouse*,
                           DataWarehouse* new_dw);

  void
  multiLevelPressureSolve(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw,
                          GridP grid,
                          const MaterialSubset*,
                          const MaterialSubset*);

  void
  bogus_imp_delP(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse*,
                 DataWarehouse* new_dw);

  void
  compute_refluxFluxes_RHS(const ProcessorGroup*,
                           const PatchSubset* coarsePatches,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw);

  void
  apply_refluxFluxes_RHS(const ProcessorGroup*,
                         const PatchSubset* coarsePatches,
                         const MaterialSubset*,
                         DataWarehouse*,
                         DataWarehouse* new_dw);

  //__________________________________
  //   M O D E L S

  void
  zeroModelSources(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse*,
                   DataWarehouse*);

  void
  updateVolumeFraction(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* matls,
                       DataWarehouse*,
                       DataWarehouse*);

  void
  computeLagrangian_Transported_Vars(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* matls,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);

  //__________________________________
  //   O T H E R

  void
  TestConservation(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse*,
                   DataWarehouse*);

  void
  hydrostaticPressureAdjustment(const Patch* patch,
                                const CCVariable<double>& rho_micro_CC,
                                CCVariable<double>& press_CC);

  IntVector
  upwindCell_X(const IntVector& c, const double& var, double is_logical_R_face);

  IntVector
  upwindCell_Y(const IntVector& c, const double& var, double is_logical_R_face);

  IntVector
  upwindCell_Z(const IntVector& c, const double& var, double is_logical_R_face);

  double
  getRefPress() const
  {
    return d_ref_press;
  }

  Vector
  getGravity() const
  {
    return d_gravity;
  }

  // debugging variables
  int d_dbgVar1{ 0 };
  int d_dbgVar2{ 0 };
  std::vector<IntVector> d_dbgIndices;

  // flags
  bool d_doRefluxing{ false };
  int d_surroundingMatl_indx{ -9 };
  bool d_impICE{ false };
  bool d_with_mpm{ false };
  bool d_with_rigid_mpm{ false };
  bool d_viscousFlow{ false };
  bool d_applyHydrostaticPress{ true };

  int d_max_iter_equilibration{ 100 };
  int d_max_iter_implicit{ 10 };
  int d_iters_before_timestep_recompute{ 100 };
  double d_outer_iter_tolerance{ 1.0e-30 };

  // ADD HEAT VARIABLES
  bool d_add_heat{ false };
  double d_add_heat_t_start{ 0.0 };
  double d_add_heat_t_final{ 1.0 };
  std::vector<int> d_add_heat_matls;
  std::vector<double> d_add_heat_coeff;

  double d_ref_press{ 0.0 };

public:
  // Particle state - communicated from MPM
  inline void
  setParticleGhostLayer(Ghost::GhostType type, int ngc)
  {
    particle_ghost_type  = type;
    particle_ghost_layer = ngc;
  }

  inline void
  getParticleGhostLayer(Ghost::GhostType& type, int& ngc)
  {
    type = particle_ghost_type;
    ngc  = particle_ghost_layer;
  }

private:
  //! so all components can know how many particle ghost cells to ask for
  Ghost::GhostType particle_ghost_type{ Ghost::None };
  int particle_ghost_layer{ 0 };

  // For AMR staff

protected:
  virtual void
  refineBoundaries(const Patch* patch,
                   CCVariable<double>& val,
                   DataWarehouse* new_dw,
                   const VarLabel* label,
                   int matl,
                   double factor);

  virtual void
  refineBoundaries(const Patch* patch,
                   CCVariable<Vector>& val,
                   DataWarehouse* new_dw,
                   const VarLabel* label,
                   int matl,
                   double factor);

  virtual void
  refineBoundaries(const Patch* patch,
                   SFCXVariable<double>& val,
                   DataWarehouse* new_dw,
                   const VarLabel* label,
                   int matl,
                   double factor);

  virtual void
  refineBoundaries(const Patch* patch,
                   SFCYVariable<double>& val,
                   DataWarehouse* new_dw,
                   const VarLabel* label,
                   int matl,
                   double factor);

  virtual void
  refineBoundaries(const Patch* patch,
                   SFCZVariable<double>& val,
                   DataWarehouse* new_dw,
                   const VarLabel* label,
                   int matl,
                   double factor);

  virtual void
  addRefineDependencies(Task* task, const VarLabel* var, int step, int nsteps);

  MaterialSubset* d_press_matl{ nullptr };
  MaterialSet* d_press_matlSet{ nullptr };

private:
#ifdef HAVE_HYPRE
  const VarLabel* hypre_solver_label;
#endif
  friend class MPMICE;
  friend class AMRICE;
  friend class impAMRICE;

  std::unique_ptr<ICELabel> d_ice_labels{ nullptr };
  SchedulerP d_subsched{ nullptr };

  std::string d_delT_scheme{ "aggressive" };
  bool d_recompileSubsched{ false };
  double d_EVIL_NUM{ -9.99e30 };
  double d_SMALL_NUM{ 1.0e-100 };
  double d_CFL{ d_EVIL_NUM };
  double d_delT_speedSoundKnob{ 1.0 };
  // used to modify the diffusion constribution to delT calc.
  double d_delT_diffusionKnob{ 1.0 };
  Vector d_gravity{ 0.0, 0.0, 0.0 };
  Vector d_fixedPressGrad = Vector(d_EVIL_NUM);
  Ghost::GhostType d_gn{ Ghost::None };
  Ghost::GhostType d_gac{ Ghost::AroundCells };

  //__________________________________
  // Misc
  std::unique_ptr<CustomBCDriver::customBC_globalVars> d_BC_globalVars;
  std::unique_ptr<customInitialize_basket> d_customInitialize_basket;

  std::unique_ptr<Advector> d_advector{ nullptr };
  int d_OrderOfAdvection{ 1 };
  bool d_useCompatibleFluxes{ true };
  bool d_clampSpecificVolume{ false };

  std::unique_ptr<Turbulence> d_turbulence{ nullptr };
  std::unique_ptr<WallShearStress> d_WallShearStressModel{ nullptr };

  std::vector<std::unique_ptr<AnalysisModule>> d_analysisModules;

  // exchange coefficients
  std::unique_ptr<ExchangeModel> d_exchModel{ nullptr };

  // flags for the conservation test
  struct conservationTest_flags
  {
    bool onOff;
    bool mass;
    bool momentum;
    bool energy;
    bool exchange;
  };
  std::unique_ptr<conservationTest_flags> d_conservationTest{ nullptr };

  //______________________________________________________________________
  //        models
  std::vector<std::unique_ptr<ModelInterface>> d_models;

  //______________________________________________________________________
  //      FUNCTIONS
  inline bool
  isEqual(const Vector& a, const Vector& b)
  {
    return (a.x() == b.x() && a.y() == b.y() && a.z() == b.z());
  };

  inline bool
  isEqual(const double a, const double b)
  {
    return a == b;
  };

  /*_____________________________________________________________________
    Purpose~  Returns if any CCVariable == value.  Useful for detecting
    uninitialized variables
    _____________________________________________________________________  */
  template<class T>
  bool
  isEqual(T value, CellIterator& iter, CCVariable<T>& q_CC, IntVector& cell)
  {
    for (; !iter.done(); iter++) {
      IntVector c = *iter;
      if (isEqual(q_CC[c], value)) {
        cell = c;
        return true;
      }
    }
    cell = IntVector(0, 0, 0);
    return false;
  }
};

} // End namespace Uintah

#endif
