/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2018 Parresia Research Limited, New Zealand
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

#ifndef __ROUSSELIER_YIELD_MODEL_H__
#define __ROUSSELIER_YIELD_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

/*! \class  YieldCond_Rousselier
 *  \brief  Rousselier Yield Condition.
 *  \author Biswajit Banerjee
 *  \author C-SAFE and Department of Mechanical Engineering
 *  \author University of Utah
 *  \warning The stress tensor is the Cauchy stress and not the
 *           Kirchhoff stress.

 References:

 1) Bernauer, G. and Brocks, W., 2002, Fatigue Fract. Engg. Mater. Struct.,
 25, 363-384.

 The yield condition is given by
 \f[
 f(\sigma,k,T) =
 \frac{\sigma_{eq}}{1-phi} +
 D \sigma_1 phi \exp \left(\frac{(1/3)Tr(\sigma)}{(1-phi)\sigma_1}\right) -
 \sigma_f = 0
 \f]
 where \f$\f(\sigma,k,T)\f$ is the yield condition,
 \f$\sigma\f$ is the Cauchy stress,
 \f$k\f$ is a set of internal variable that evolve with time,
 \f$T\f$ is the temperature,
 \f$\sigma_{eq} = \sqrt{3J_2}\f$ is the von Mises equivalent stress given by
 \f$ \sigma_{eq} = \sqrt{\frac{3}{2}\sigma^{d}:\sigma^{d}}\f$ where
 \f$\sigma^{d}\f$ is the deviatoric part of the Cauchy stress,
 \f$\sigma_{f}\f$ is the flow stress,
 \f$D,\sigma_1\f$ are material constants, and
 \f$\phi\f$ is the porosity (void volume fraction).
*/

class YieldCond_Rousselier : public YieldCondition
{

public:
  /* Constants needed for Rousselier model */
  struct Params
  {
    double D;
    double sigma_1;
  };

  YieldCond_Rousselier(ProblemSpecP& ps);
  YieldCond_Rousselier(const YieldCond_Rousselier* cm);
  YieldCond_Rousselier&
  operator=(const YieldCond_Rousselier&) = delete;
  ~YieldCond_Rousselier() override = default;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["D"]       = d_params.D;
    params["sigma_1"] = d_params.sigma_1;
    return params;
  }

  //! Evaluate the yield function.
  double
  evalYieldCondition(const double equivStress,
                     const double flowStress,
                     const double traceOfCauchyStress,
                     const double porosity,
                     double& sig) override;

  double
  evalYieldCondition(const Uintah::Matrix3& xi,
                     const ModelStateBase* state) override;

  std::pair<double, Util::YieldStatus>
  evalYieldCondition(const ModelStateBase* state) override;

  double
  evalYieldConditionMax(const ModelStateBase* state) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the derivative of the yield function (f) 
    with respect to various quantities

    sigma = Cauchy stress
    sigmaDev = deviatoric stress
    xi = sigmaDev - beta
    beta =  backstress 
    p = volumetric stress = 1/3 Tr(sigma)
    q = sqrt(3 J_2), J_2 = 2nd invariant of sigmaDev
  */
  /////////////////////////////////////////////////////////////////////////
  void
  df_dsigma(const Uintah::Matrix3& stress,
            const double flowStress,
            const double porosity,
            Uintah::Matrix3& derivative) override;

  void
  df_dsigma(const Uintah::Matrix3& xi,
            const ModelStateBase* state,
            Uintah::Matrix3& df_dsigma) override;

  void
  df_dsigmaDev(const Uintah::Matrix3& stress,
               const double flowStress,
               const double porosity,
               Uintah::Matrix3& derivative) override;

  void
  df_dxi(const Uintah::Matrix3& xi,
         const ModelStateBase* state,
         Uintah::Matrix3& df_xi) override;

  void
  df_dsigmaDev_dbeta(const Uintah::Matrix3& xi,
                     const ModelStateBase* state,
                     Uintah::Matrix3& df_ds,
                     Uintah::Matrix3& df_dbeta) override;

  double
  df_dp(const ModelStateBase* state) override;

  double
  df_dq(const ModelStateBase* state) override;

  double
  df_dplasticStrain(const Uintah::Matrix3& xi,
                    const double& d_sigy_dep,
                    const ModelStateBase* state) override;

  double
  df_dporosity(const Uintah::Matrix3& xi, 
               const ModelStateBase* state) override;

  double
  d2f_dp_depsVol(const ModelStateBase* state,
                 const PressureModel* eos,
                 const ShearModulusModel* shear,
                 const InternalVariableModel* intvar) override;

  double
  d2f_dp_depsDev(const ModelStateBase* state,
                 const PressureModel* eos,
                 const ShearModulusModel* shear,
                 const InternalVariableModel* intvar) override;

  double
  d2f_dq_depsVol(const ModelStateBase* state,
                 const PressureModel* eos,
                 const ShearModulusModel* shear,
                 const InternalVariableModel* intvar) override;

  double
  d2f_dq_depsDev(const ModelStateBase* state,
                 const PressureModel* eos,
                 const ShearModulusModel* shear,
                 const InternalVariableModel* intvar) override;
  double
  df_depsVol(const ModelStateBase* state,
             const PressureModel* eos,
             const ShearModulusModel* shear,
             const InternalVariableModel* intvar) override;

  double
  df_depsDev(const ModelStateBase* state,
             const PressureModel* eos,
             const ShearModulusModel* shear,
             const InternalVariableModel* intvar) override;


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
  computeElasPlasTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                const Matrix3& sigma,
                                double sigY,
                                double dsigYdV,
                                double porosity,
                                double voidNuclFac,
                                Uintah::TangentModulusTensor& Cep) override;

  double
  getInternalPoint(const ModelStateBase* state_old,
                   const ModelStateBase* state_new) override
  {
    return 0.0;
  }

private:
  Params d_params;

  Uintah::Matrix3 
  df_dsigma(const Uintah::Matrix3& s_dev,
            double p,
            double phi) const;

  void
  computeTangentModulus(const Uintah::TangentModulusTensor& Ce,
                        const Uintah::Matrix3& f_sigma,
                        double f_q1,
                        double f_q2,
                        double h_q1,
                        double h_q2,
                        Uintah::TangentModulusTensor& Cep);

};

} // End namespace Vaango

#endif // __ROUSSELIER_YIELD_MODEL_H__
