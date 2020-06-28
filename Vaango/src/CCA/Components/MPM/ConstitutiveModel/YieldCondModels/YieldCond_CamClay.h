/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __CAMCLAY_YIELD_MODEL_H__
#define __CAMCLAY_YIELD_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_BorjaPressure.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*!
 \class  YieldCond_CamClay
 \brief  The CAM-Clay yield condition

 The yield condition is given by
   Phi(p,q) = (q/M)^2 + p(p-p_c) <= 0
 where
   p = I_1/3 = 1/3 tr(sigma)
   q = sqrt(3 J2) = sqrt(3/2 s:s)
   sigma = Cauchy stress
   s = deviatoric stress = sigma - p I
   M = material constant = slope of critical state line
   p_c = consolidation pressure
         defined by the rate equation

         1/p_c dp_c/dt = 1/(lambdatilde - kappatilde) depsp_v/dt

         where lambdatilde = material constant
               kappatilde = material constant
               epsp_v = volumetric plastic strain

   p_c is computed using the Borja plasticity model and is
       treated as an internal variable.
*/

class YieldCond_CamClay : public YieldCondition
{

public:
  //! Constructor
  /*! Creates a YieldCond_CamClay function object */
  YieldCond_CamClay(Uintah::ProblemSpecP& ps, IntVar_BorjaPressure* intvar);
  YieldCond_CamClay(const YieldCond_CamClay* cm);
  YieldCond_CamClay& operator=(const YieldCond_CamClay&) = delete;

  //! Destructor
  ~YieldCond_CamClay() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["M"] = d_M;
    return params;
  }

  //--------------------------------------------------------------
  // Compute value of yield function
  //--------------------------------------------------------------
  std::pair<double, Util::YieldStatus>
  evalYieldCondition(const ModelStateBase* state) override;
  double
  evalYieldConditionMax(const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
  //--------------------------------------------------------------
  double
  df_dp(const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Compute df/dq  where q = sqrt(3 J_2), J_2 = 2nd invariant deviatoric stress
  //--------------------------------------------------------------
  double
  df_dq(const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dp)
  //--------------------------------------------------------------
  double
  d2f_dp_depsVol(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) override;

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dp)
  //--------------------------------------------------------------
  double
  d2f_dp_depsDev(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) override;

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dq)
  //--------------------------------------------------------------
  double
  d2f_dq_depsVol(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) override;

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dq)
  //--------------------------------------------------------------
  double
  d2f_dq_depsDev(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) override;

  //--------------------------------------------------------------
  // Compute df/depse_v
  //--------------------------------------------------------------
  double
  df_depsVol(const ModelStateBase* state,
             const MPMEquationOfState* eos,
             const ShearModulusModel* shear) override;

  //--------------------------------------------------------------
  // Compute df/depse_s
  //--------------------------------------------------------------
  double
  df_depsDev(const ModelStateBase* state,
             const MPMEquationOfState* eos,
             const ShearModulusModel* shear) override;

  double
  getInternalPoint(const ModelStateBase* state_old,
                   const ModelStateBase* state_new) override
  {
    return 0.0;
  }

  //================================================================================
  // Other options below.
  //================================================================================

  // Evaluate the yield function.
  double
  evalYieldCondition(const double p,
                     const double q,
                     const double p_c,
                     const double dummy1,
                     double& dummy2) override;

  // Evaluate yield condition (s = deviatoric stress = sigDev
  //                           p = state->pressure
  //                           p_c = state->yieldStress)
  double
  evalYieldCondition(const Uintah::Matrix3& sigDev,
                     const ModelStateBase* state) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
    with respect to \f$\sigma_{ij}\f$.
  */
  /////////////////////////////////////////////////////////////////////////
  Uintah::Matrix3
  df_dsigma(const Uintah::Matrix3& stress,
            const double p_c,
            const double dummy) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
    with respect to \f$s_{ij}\f$.

    This is for the associated flow rule with \f$s_{ij}\f$ being
    the deviatoric stress.
  */
  /////////////////////////////////////////////////////////////////////////
  Uintah::Matrix3
  df_dsigmaDev(const Uintah::Matrix3& stress,
               const double dummy1,
               const double dummy2) override;

  /*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
  Uintah::Matrix3
  df_dsigma(const Uintah::Matrix3& xi,
            const ModelStateBase* state) override;

  /*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
      where \f$s\f$ is deviatoric part of Cauchy stress and
      \f$\beta\f$ is the backstress */
  Uintah::Matrix3
  df_dxi(const Uintah::Matrix3& xi,
         const ModelStateBase* state) override;

  /* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
  std::pair<Uintah::Matrix3, Uintah::Matrix3>
  df_dsigmaDev_dbeta(const Uintah::Matrix3& xi,
                     const ModelStateBase* state) override;

  /*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$)*/
  double
  df_dplasticStrain(const Uintah::Matrix3& xi,
                    const double& d_sigy_dep,
                    const ModelStateBase* state) override;

  /*! Derivative with respect to the porosity (\f$\epsilon^p \f$)*/
  double
  df_dporosity(const Uintah::Matrix3& xi, const ModelStateBase* state) override;

  /*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
  double
  eval_h_alpha(const Uintah::Matrix3& xi, const ModelStateBase* state) override;

  /*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
  double
  eval_h_phi(const Uintah::Matrix3& xi,
             const double& factorA,
             const ModelStateBase* state) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the elastic-plastic tangent modulus.
  */
  /////////////////////////////////////////////////////////////////////////
  void
  computeElasPlasTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                const Uintah::Matrix3& sigma,
                                double sigY,
                                double dsigYdep,
                                double porosity,
                                double voidNuclFac,
                                Uintah::TangentModulusTensor& Cep) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the factor \f$h_1\f$ for plastic strain

    \f[
    h_1 = \frac{\sigma : f_{\sigma}}{\sigma_Y}
    \f]

    \return factor
  */
  /////////////////////////////////////////////////////////////////////////
  inline double
  computePlasticStrainFactor(double sigma_f_sigma, double sigma_Y);

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the continuum elasto-plastic tangent modulus
    assuming associated flow rule.

    \f[
    C_{ep} = C_{e} - \frac{(C_e:f_{\sigma})\otimes(f_{\sigma}:C_e)}
    {-f_q.h_q + f_{\sigma}:C_e:f_{\sigma}}
    \f]

    \return TangentModulusTensor \f$ C_{ep} \f$.
  */
  /////////////////////////////////////////////////////////////////////////
  void
  computeTangentModulus(const Uintah::TangentModulusTensor& Ce,
                        const Uintah::Matrix3& f_sigma,
                        double f_q1,
                        double h_q1,
                        Uintah::TangentModulusTensor& Cep);

private:
  double d_M; // slope of critical state line
  IntVar_BorjaPressure* d_intvar;

};

} // End namespace Uintah

#endif // __CAMCLAY_YIELD_MODEL_H__
