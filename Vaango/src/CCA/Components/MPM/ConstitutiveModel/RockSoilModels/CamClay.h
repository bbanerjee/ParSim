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

#ifndef __CAM_CLAY_PLASTIC_H__
#define __CAM_CLAY_PLASTIC_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/BorjaEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_BorjaPressure.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>

#include <CCA/Ports/DataWarehouseP.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <cmath>

namespace Uintah {

class MPMLabel;
class MPMFlags;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class CamClay
  \brief Hyperelastic Cam Clay Plasticity Model
*/
/////////////////////////////////////////////////////////////////////////////

class CamClay : public ConstitutiveModel
{

public:
  const VarLabel* pStrainLabel;
  const VarLabel* pElasticStrainLabel;
  const VarLabel* pDeltaGammaLabel;

  const VarLabel* pStrainLabel_preReloc;
  const VarLabel* pElasticStrainLabel_preReloc;
  const VarLabel* pDeltaGammaLabel_preReloc;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief constructors */
  ////////////////////////////////////////////////////////////////////////
  explicit CamClay(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
  explicit CamClay(const CamClay* cm);
  CamClay&
  operator=(const CamClay& cm) = delete;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief destructor  */
  ////////////////////////////////////////////////////////////////////////
  ~CamClay() override;

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
  /*! \brief Computes and requires for the initialization of the model */
  ////////////////////////////////////////////////////////////////////////
  void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief initialize  each particle's constitutive model data */
  ////////////////////////////////////////////////////////////////////////
  void
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
  /*! \brief Computes ad requires for each time step */
  ////////////////////////////////////////////////////////////////////////
  void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute stress at each particle in the patch
  */
  ////////////////////////////////////////////////////////////////////////
  void
  computeStressTensor(const PatchSubset* patches,
                      const MPMMaterial* matl,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  addComputesAndRequires(Task*,
                         const MPMMaterial*,
                         const PatchSet*,
                         const bool,
                         const bool) const override{};

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
                             const PatchSet* patches) const override{};

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  getDamageParameter(const Patch* patch,
                     ParticleVariable<int>& damage,
                     int dwi,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw) override{};

  ////////////////////////////////////////////////////////////////////////
  /*! \brief allocate CM data requires */
  ////////////////////////////////////////////////////////////////////////
  void
  allocateCMDataAddRequires(Task* task,
                            const MPMMaterial* matl,
                            const PatchSet* patch,
                            Uintah::MPMLabel* lb) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief allocate cm data add */
  ////////////////////////////////////////////////////////////////////////
  void
  allocateCMDataAdd(DataWarehouse* new_dw,
                    ParticleSubset* subset,
                    ParticleLabelVariableMap* newState,
                    ParticleSubset* delset,
                    DataWarehouse* old_dw) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  scheduleCheckNeedAddMPMMaterial(Task* task,
                                  const MPMMaterial* matl,
                                  const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \brief Put documentation here. */
  ////////////////////////////////////////////////////////////////////////
  void
  checkNeedAddMPMMaterial(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override;

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

private:
  enum class SolveStatus
  {
    OK,
    INVALID_VALUE,
    CONVERGENCE_FAILURE
  };

  inline static const Matrix3 one{
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0
  };
  inline static const Matrix3 zero{ 0.0 };
  static const double sqrtThreeTwo;
  static const double sqrtTwoThird;

  Vaango::MPMEquationOfState* d_eos;
  Vaango::ShearModulusModel* d_shear;
  Vaango::YieldCondition* d_yield;
  std::shared_ptr<Vaango::IntVar_BorjaPressure> d_intvar;

  void
  initializeLocalMPMLabels();

  /* Do elastic update */
  std::tuple<Matrix3, Matrix3, Matrix3, double, double>
  elasticUpdate(const Matrix3& strain_elast_tr,
                const Matrix3& nn,
                double pDeltaGamma_old,
                double pPc_old,
                const Vaango::ModelState_CamClay& state);

  /* Do plastic update */
  std::
    tuple<Matrix3, Matrix3, Matrix3, double, double, SolveStatus, std::string>
    plasticUpdate(particleIndex idx,
                  const MPMMaterial* matl,
                  double ftrial,
                  const Matrix3& elasticStrain_old,
                  const Matrix3& strain_elast_tr,
                  double strain_elast_v_tr,
                  double strain_elast_s_tr,
                  const Matrix3& nn,
                  double pDeltaGamma_old,
                  double pPc_old,
                  Vaango::ModelState_CamClay& state);

  /* Compute the trial state */
  std::tuple<Matrix3, Matrix3, Vaango::ModelState_CamClay>
  computeTrialState(particleIndex idx,
                    const MPMMaterial* matl,
                    double delT,
                    const Matrix3& rateOfDef_new,
                    const Matrix3& elasticStrain_old,
                    double pVol_new,
                    double pMass,
                    double pPc,
                    double rho_0,
                    double rho_cur);

  /* Do the main Newton solve loop */
  std::tuple<CamClay::SolveStatus, int, double, double, std::string>
  doNewtonSolve(particleIndex idx,
                const MPMMaterial* matl,
                double ftrial,
                const Matrix3& strain_elast_tr,
                double strain_elast_v_tr,
                double strain_elast_s_tr,
                const Matrix3& elasticStrain_old,
                double& delgamma_old,
                double& strain_elast_v,
                double& strain_elast_s,
                Vaango::ModelState_CamClay& state);

  /* Compute delta gamma incremenet for each newton iteration*/
  std::tuple<double, std::vector<double>>
  computeDeltaGammaIncrement(Vaango::ModelState_CamClay& state,
                             double delgamma_k,
                             double dfdp,
                             double dfdq,
                             double rv,
                             double rs,
                             double rf);

  /* Check againt Borja-Tamagnini matrices*/
  Matrix3
  computeBorjaTamagniniAmatrix(double lambda_hat,
                               double kappa_hat,
                               double eps_v_e,
                               double eps_v0_e,
                               double p,
                               double p_c,
                               double q,
                               double M,
                               double mu,
                               double p0,
                               double alpha,
                               double eps_s_e,
                               double deltaPhi);
};

} // End namespace Uintah

#endif // __CAM_CLAY_PLASTIC_H__
