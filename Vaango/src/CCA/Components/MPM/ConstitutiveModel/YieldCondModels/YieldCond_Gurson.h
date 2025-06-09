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

#ifndef __BB_GURSON_YIELD_MODEL_H__
#define __BB_GURSON_YIELD_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_Metal.h>
#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/FlowStressModel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

//////////////////////////////////////////////////////////////////////
/*!
  \class  YieldCond_Gurson
  \brief  Gurson-Tvergaard-Needleman Yield Condition.
  \author Biswajit Banerjee \n
  C-SAFE and Department of Mechanical Engineering
  University of Utah
  \warning The stress tensor is the Cauchy stress and not the
  Kirchhoff stress.

  References:

  1) Bernauer, G. and Brocks, W., 2002, Fatigue Fract. Engg. Mater. Struct.,
  25, 363-384.
  2) Ramaswamy, S. and Aravas, N., 1998, Comput. Methods Appl. Mech. Engrg.,
  163, 33-53.

  The yield condition is given by
  \f[
  f(\sigma,k,T) =
  \frac{\sigma_{eq}^2}{\sigma_f^2} +
  2 q_1 \phi_* \cosh \left(q_2 \frac{Tr(\sigma)}{2\sigma_f}\right) -
  (1+q_3 \phi_*^2) = 0
  \f]
  where \f$f(\sigma,k,T)\f$ is the yield condition,
  \f$\sigma\f$ is the Cauchy stress,
  \f$k\f$ is a set of internal variable that evolve with time,
  \f$T\f$ is the temperature,
  \f$\sigma_{eq}\f$ is the von Mises equivalent stress given by
  \f$ \sigma_{eq} = \sqrt{\frac{3}{2}\sigma^{d}:\sigma^{d}}\f$ where
  \f$\sigma^{d}\f$ is the deviatoric part of the Cauchy stress,
  \f$\sigma_{f}\f$ is the flow stress,
  \f$q_1,q_2,q_3\f$ are material constants, and
  \f$\phi_*\f$ is the porosity (damage) function.

  The damage function is given by
  \f$ \phi_* = \phi \f$ for \f$ \phi \le \phi_c \f$,
  \f$ \phi_* = \phi_c + k (\phi - \phi_c) \f$ for \f$ \phi > \phi_c \f$, where
  \f$ k \f$ is constant, and \f$ \phi \f$ is the porosity (void volume
  fraction).
*/
//////////////////////////////////////////////////////////////////////

class YieldCond_Gurson : public YieldCondition
{

public:
  /*! \struct CMData
    \brief Constants needed for GTN model */
  struct CMData
  {
    double q1;  /*< Constant q_1 */
    double q2;  /*< Constant q_2 */
    double q3;  /*< Constant q_3 */
    double k;   /*< Constant k */
    double f_c; /*< Critical void volume fraction */
  };

  /*! Constructor
    Creates a Gurson Yield Function object */
  YieldCond_Gurson(Uintah::ProblemSpecP& ps, 
                   IntVar_Metal* intvar,
                   const Uintah::FlowStressModel* flow);
  YieldCond_Gurson(const YieldCond_Gurson* cm);
  YieldCond_Gurson&
  operator=(const YieldCond_Gurson&) = delete;

  //! Destructor
  virtual ~YieldCond_Gurson() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["q1"]  = d_CM.q1;
    params["q2"]  = d_CM.q2;
    params["q3"]  = d_CM.q3;
    params["k"]   = d_CM.k;
    params["f_c"] = d_CM.f_c;
    return params;
  }

  //! Evaluate the yield function.
  std::pair<double, Util::YieldStatus>
  evalYieldCondition(const ModelStateBase* state) override;

  double
  computeYieldFunction(const ModelStateBase* state) const override;

  double
  evalYieldCondition(const Uintah::Matrix3& stress,
                     const ModelStateBase* state) override;

  double
  evalYieldConditionMax(const ModelStateBase* state) override;

  /*! Derivative with respect to the Cauchy stress (\f$\sigma \f$)*/
  Uintah::Matrix3
  df_dsigma(const ModelStateBase* state) override;

  Uintah::Matrix3
  df_dsigma(const Uintah::Matrix3& stress,
            const ModelStateBase* state) override;

  /*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
      where \f$s\f$ is deviatoric part of Cauchy stress and
      \f$\beta\f$ is the backstress */
  Uintah::Matrix3
  df_dxi(const Uintah::Matrix3& stress,
         const ModelStateBase* state) override;

  /* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
  std::pair<Uintah::Matrix3, Uintah::Matrix3>
  df_dsigmaDev_dbeta(const Uintah::Matrix3& stress,
                     const ModelStateBase* state) override;
  double
  df_dbeta_p(const Matrix3& stress, const ModelStateBase* state) const;

  //--------------------------------------------------------------
  // Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
  //--------------------------------------------------------------
  double
  df_dp([[maybe_unused]] const ModelStateBase* state) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute df/dq  where q = sqrt(3 J_2), J_2 = 2nd invariant deviatoric stress
  //--------------------------------------------------------------
  double
  df_dq([[maybe_unused]] const ModelStateBase* state) override
  {
    return 0.0;
  };

  /*! Derivative with respect to internal variables */
  void
  df_dintvar([[maybe_unused]] const ModelStateBase* state,
             [[maybe_unused]] MetalIntVar& df_dintvar) const override;

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dp)
  //--------------------------------------------------------------
  double
  d2f_dp_depsVol([[maybe_unused]] const ModelStateBase* state,
                 [[maybe_unused]] const MPMEquationOfState* eos,
                 [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dp)
  //--------------------------------------------------------------
  double
  d2f_dp_depsDev([[maybe_unused]] const ModelStateBase* state,
                 [[maybe_unused]] const MPMEquationOfState* eos,
                 [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dq)
  //--------------------------------------------------------------
  double
  d2f_dq_depsVol([[maybe_unused]] const ModelStateBase* state,
                 [[maybe_unused]] const MPMEquationOfState* eos,
                 [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dq)
  //--------------------------------------------------------------
  double
  d2f_dq_depsDev([[maybe_unused]] const ModelStateBase* state,
                 [[maybe_unused]] const MPMEquationOfState* eos,
                 [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute df/depse_v
  //--------------------------------------------------------------
  double
  df_depsVol([[maybe_unused]] const ModelStateBase* state,
             [[maybe_unused]] const MPMEquationOfState* eos,
             [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  //--------------------------------------------------------------
  // Compute df/depse_s
  //--------------------------------------------------------------
  double
  df_depsDev([[maybe_unused]] const ModelStateBase* state,
             [[maybe_unused]] const MPMEquationOfState* eos,
             [[maybe_unused]] const ShearModulusModel* shear) override
  {
    return 0.0;
  };

  double
  getInternalPoint([[maybe_unused]] const ModelStateBase* state_old,
                   [[maybe_unused]] const ModelStateBase* state_new) override
  {
    return 0.0;
  }

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
                        double f_q2,
                        double h_q1,
                        double h_q2,
                        Uintah::TangentModulusTensor& Cep);

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the elastic-plastic tangent modulus.
  */
  /////////////////////////////////////////////////////////////////////////
  void
  computeElasPlasTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                const Uintah::Matrix3& sigma,
                                double sigY,
                                double dsigYdV,
                                double porosity,
                                double voidNuclFac,
                                Uintah::TangentModulusTensor& Cep) override;

private:
  CMData d_CM;
  IntVar_Metal* d_intvar;
  const Uintah::FlowStressModel* d_flow;

  /* Compute the yield function */
  double
  computeYieldFunction(const Matrix3& stress,
                       const ModelStateBase* state) const;

  /*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$)*/
  double
  df_dplasticStrain(const ModelStateBase* state) const;

  /*! Derivative with respect to the porosity (\f$\phi\f$)*/
  double
  df_dporosity(const ModelStateBase* state) const;

};

} // End namespace Uintah

#endif // __BB_GURSON_YIELD_MODEL_H__
