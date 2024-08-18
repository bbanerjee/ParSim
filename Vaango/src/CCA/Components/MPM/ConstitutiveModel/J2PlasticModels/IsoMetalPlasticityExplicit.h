/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __SMALL_STRAIN_LARGE_ROTATION_PLASTIC_H__
#define __SMALL_STRAIN_LARGE_ROTATION_PLASTIC_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>

#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/DamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_MetalIso.h>
#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/FlowStressModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_Metal.h>
#include <CCA/Components/MPM/ConstitutiveModel/KinHardeningModels/KinematicHardeningModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MeltTempModels/MeltingTempModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/SpecHeatModels/SpecificHeatModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/StabilityModels/StabilityCheck.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/DeformationState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

#include <CCA/Ports/DataWarehouseP.h>

#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cmath>

namespace Uintah {

class MPMLabel;
class MPMFlags;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class IsoMetalPlasticityExplicit
  \brief Small strain/Large rotation Hypo-Elastic Plastic Constitutive Model
  \author Biswajit Banerjee \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n

  The rate of deformation and stress is rotated to material configuration
  before the updated values are calculated.  Failure is calculated and
  failed particles are tagged and their stresses are changed acoording to
  the failure model.

  Needs :
  1) Isotropic elastic moduli.
  2) Isotropic hardening flow rule in the form of a Plasticity Model.
  3) Kinematic hardening rule.
  4) Yield condition.
  5) Stability condition.
  6) Damage model.
  7) Shear modulus model.
  8) Melting temperature model.
  9) Specific heat model.

  \warning Only isotropic materials, von-Mises type yield conditions,
  associated flow rule
*/
/////////////////////////////////////////////////////////////////////////////

class IsoMetalPlasticityExplicit : public ConstitutiveModel
{

public:
  inline static const Matrix3 zero{
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  };
  inline static const Matrix3 one{
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0
  };
  static const double sqrtThreeTwo;
  static const double sqrtTwoThird;

  // Create datatype for storing model parameters
  struct CMData
  {
    double Bulk;       /*< Bulk modulus */
    double Shear;      /*< Shear Modulus */
    double CTE;        /*< Coeff. of thermal expansion */
    double Chi;        /*< Taylor-Quinney coefficient */
    double sigma_crit; /*< Critical stress */
  };

  // Create datatype for storing porosity parameters
  struct PorosityData
  {
    double f0;     /*< Initial mean porosity */
    double f0_std; /*< Initial standard deviation of porosity */
    double fc;     /*< Critical porosity */
    double fn;     /*< Volume fraction of void nucleating particles */
    double en;     /*< Mean strain for nucleation */
    double sn;     /*< Standard deviation of strain for nucleation */
    std::string porosityDist; /*< Initial porosity distribution*/
  };

  // Create datatype for storing damage parameters
  struct ScalarDamageData
  {
    double D0;     /*< Initial mean scalar damage */
    double D0_std; /*< Initial standard deviation of scalar damage */
    double Dc;     /*< Critical scalar damage */
    std::string scalarDamageDist; /*< Initial damage distrinution */
  };

  const VarLabel* pStrainRateLabel;
  const VarLabel* pPlasticStrainLabel;
  const VarLabel* pPlasticStrainRateLabel;
  const VarLabel* pEqStrainRateLabel;
  const VarLabel* pEqPlasticStrainLabel;
  const VarLabel* pEqPlasticStrainRateLabel;
  const VarLabel* pDamageLabel;
  const VarLabel* pPorosityLabel;
  const VarLabel* pLocalizedLabel;
  const VarLabel* pIntVarLabel;
  const VarLabel* pDStressDIntVarLabel;

  const VarLabel* pStrainRateLabel_preReloc;
  const VarLabel* pPlasticStrainLabel_preReloc;
  const VarLabel* pPlasticStrainRateLabel_preReloc;
  const VarLabel* pEqStrainRateLabel_preReloc;
  const VarLabel* pEqPlasticStrainLabel_preReloc;
  const VarLabel* pEqPlasticStrainRateLabel_preReloc;
  const VarLabel* pDamageLabel_preReloc;
  const VarLabel* pPorosityLabel_preReloc;
  const VarLabel* pLocalizedLabel_preReloc;
  const VarLabel* pIntVarLabel_preReloc;
  const VarLabel* pDStressDIntVarLabel_preReloc;

protected:
  CMData d_initialData;
  PorosityData d_porosity;
  ScalarDamageData d_scalarDam;

  double d_tol{ 1.0e-8 };
  double d_initialMaterialTemperature{ 294.0 };
  double d_isothermal{ 0.0 };
  bool d_doIsothermal{ false };
  bool d_useModifiedEOS{ false };
  bool d_evolvePorosity{ false };
  bool d_evolveDamage{ false };
  bool d_computeSpecificHeat{ false };
  bool d_checkTeplaFailureCriterion{ false };
  bool d_doMelting{ false };
  bool d_checkStressTriax{ false };

  // Erosion algorithms
  bool d_setStressToZero{ false };
  bool d_allowNoTension{ false };

  std::unique_ptr<Vaango::MPMEquationOfState> d_eos;
  std::unique_ptr<Vaango::ShearModulusModel> d_shear;
  std::unique_ptr<Vaango::KinematicHardeningModel> d_kinematic;
  std::unique_ptr<Vaango::ElasticModuli_MetalIso> d_elastic;
  std::unique_ptr<Vaango::IntVar_Metal> d_intvar;
  std::unique_ptr<Vaango::YieldCondition> d_yield;
  std::unique_ptr<MeltingTempModel> d_melt;
  std::unique_ptr<SpecificHeatModel> d_Cp;
  std::unique_ptr<FlowStressModel> d_flow;
  std::unique_ptr<DamageModel> d_damage;
  std::unique_ptr<StabilityCheck> d_stable;

public:
  ////////////////////////////////////////////////////////////////////////
  /*! \brief constructors */
  ////////////////////////////////////////////////////////////////////////
  IsoMetalPlasticityExplicit(ProblemSpecP& ps, MPMFlags* flag);
  IsoMetalPlasticityExplicit(const IsoMetalPlasticityExplicit* cm);
  IsoMetalPlasticityExplicit&
  operator=(const IsoMetalPlasticityExplicit& cm) = delete;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief destructor  */
  ////////////////////////////////////////////////////////////////////////
  ~IsoMetalPlasticityExplicit() override;

  ModelType
  modelType() const override
  {
    return ModelType::RATE_FORM;
  }

  void
  outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  // clone
  std::unique_ptr<ConstitutiveModel>
  clone() override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief initialize  each particle's constitutive model data */
  ////////////////////////////////////////////////////////////////////////
  virtual void
  initializeCMData(const Patch* patch,
                   const MPMMaterial* matl,
                   DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief compute stable timestep for this patch */
  ////////////////////////////////////////////////////////////////////////
  virtual void
  computeStableTimestep(const Patch* patch,
                        const MPMMaterial* matl,
                        DataWarehouse* new_dw);

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  virtual void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute stress at each particle in the patch

    The plastic work is converted into a rate of temperature increase
    using an equation of the form
    \f[
       \dot{T} = \frac{\chi}{\rho C_p}(\sigma:D^p)
    \f]
  */
  ////////////////////////////////////////////////////////////////////////
  virtual void
  computeStressTensor(const PatchSubset* patches,
                      const MPMMaterial* matl,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief carry forward CM data for RigidMPM */
  ////////////////////////////////////////////////////////////////////////
  void
  carryForward(const PatchSubset* patches,
               const MPMMaterial* matl,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  addRequiresDamageParameter(Task* task,
                             const MPMMaterial* matl,
                             const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  getDamageParameter(const Patch* patch,
                     ParticleVariable<int>& damage,
                     int dwi,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  allocateCMDataAddRequires(Task* task,
                            const MPMMaterial* matl,
                            const PatchSet* patch,
                            MPMLabel* lb) const override
  {
  }

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  allocateCMDataAdd(DataWarehouse* new_dw,
                    ParticleSubset* subset,
                    ParticleLabelVariableMap* newState,
                    ParticleSubset* delset,
                    DataWarehouse* old_dw) override
  {
  }

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  scheduleCheckNeedAddMPMMaterial(Task* task,
                                  const MPMMaterial* matl,
                                  const PatchSet* patches) const override
  {
  }

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  checkNeedAddMPMMaterial(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override
  {
  }

  ////////////////////////////////////////////////////////////////////////
  /*! \brief initialize  each particle's constitutive model data */
  ////////////////////////////////////////////////////////////////////////
  void
  addParticleState(std::vector<const VarLabel*>& from,
                   std::vector<const VarLabel*>& to) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Sockets for MPM-ICE */
  ////////////////////////////////////////////////////////////////////////
  double
  computeRhoMicroCM(double pressure,
                    const double p_ref,
                    const MPMMaterial* matl,
                    double temperature,
                    double rho_guess) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Sockets for MPM-ICE */
  ////////////////////////////////////////////////////////////////////////
  void
  computePressEOSCM(double rho_m,
                    double& press_eos,
                    double p_ref,
                    double& dp_drho,
                    double& ss_new,
                    const MPMMaterial* matl,
                    double temperature) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Sockets for MPM-ICE */
  ////////////////////////////////////////////////////////////////////////
  double
  getCompressibility() override;

protected:
  std::tuple<Vaango::Tensor::Matrix6Mandel,
             Vaango::Tensor::Vector6Mandel,
             Vaango::Tensor::Vector6Mandel,
             double>
  computeElasPlasTangentModulus(Vaango::Tensor::Matrix6Mandel& C_e,
                                std::vector<Matrix3>& dsigma_deta,
                                const ModelStateBase* state) const;

  ////////////////////////////////////////////////////////////////////////
  /*! compute stress at each particle in the patch */
  ////////////////////////////////////////////////////////////////////////
  void
  computeStressTensorExplicit(const PatchSubset* patches,
                              const MPMMaterial* matl,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw);

protected:
  void
  initializeLocalMPMLabels();

  void
  getInitialPorosityData(ProblemSpecP& ps);

  void
  getInitialDamageData(ProblemSpecP& ps);

  void
  setErosionAlgorithm();

protected:
  enum class Status
  {
    INVALID_VALUE,
    CONVERGED_IN_K,
    CONVERGED_IN_F,
    CONVERGED_IN_DELTA_GAMMA,
    CONVERGENCE_FAILURE
  };

  void
  computeSubstep();

  void
  updateAsElastic(particleIndex idx,
                  const MPMMaterial* matl,
                  const ModelStateBase& state_old,
                  const DeformationState& defState_new,
                  const Matrix3& sigma_trial,
                  const Matrix3& backStress_old,
                  Matrix3& pStress_new,
                  Matrix3& pBackStress_new);

  std::
    tuple<IsoMetalPlasticityExplicit::Status, std::string, int, double, double>
    updateAsPlastic(particleIndex idx,
                    const MPMMaterial* matl,
                    double delT,
                    const std::vector<Matrix3>& sigma_eta_old,
                    const ModelStateBase& state_old,
                    const Matrix3& backStress_old,
                    const DeformationState& defState_new,
                    double f_0,
                    const Matrix3& sigma_trial,
                    double rho_cur,
                    std::vector<Matrix3>& sigma_eta_new,
                    ModelStateBase& state_new,
                    Matrix3& pStress_new,
                    Matrix3& pBackStress_new,
                    double& pdTdt_new);

  void
  updateAsFluid(particleIndex idx,
                const MPMMaterial* matl,
                double delT,
                const DeformationState& defState_new,
                const ModelStateBase& state,
                const Matrix3& sigma_old,
                Matrix3& pStress_new,
                Matrix3& pBackStress_new);

  std::tuple<double,
             Matrix3,
             double,
             IsoMetalPlasticityExplicit::Status,
             std::string,
             int>
  doNewtonSolve(particleIndex idx,
                const MPMMaterial* matl,
                double delT,
                double f_0,
                const Matrix3& sigma_trial,
                const DeformationState& defState_new,
                const std::vector<Matrix3>& sigma_eta_old,
                const ModelStateBase& state_old,
                ModelStateBase& state) const;

  std::tuple<double, Matrix3, IsoMetalPlasticityExplicit::Status, std::string>
  computeDeltaGamma(particleIndex idx,
                    int iter,
                    double delT,
                    const std::vector<Matrix3>& sigma_eta_old,
                    const DeformationState& defState_new,
                    const ModelStateBase& state,
                    const Matrix3& sigma_k,
                    double f_k,
                    double Delta_gamma_old) const;
};

} // End namespace Uintah

#endif // __SMALL_STRAIN_LARGE_ROTATION_PLASTIC_H__
