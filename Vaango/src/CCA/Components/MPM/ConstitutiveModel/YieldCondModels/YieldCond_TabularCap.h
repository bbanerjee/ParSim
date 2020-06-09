/*
 * The MIT License
 *
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#ifndef __TABULAR_CAP_YIELD_CONDITION_MODEL_H__
#define __TABULAR_CAP_YIELD_CONDITION_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCondition.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularData.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/WeibParameters.h>

#include <Core/Grid/Variables/VarLabel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

namespace Vaango {

  using Polyline = std::vector<Uintah::Point>;

/*!
  \class  YieldCond_TabularCap
  \brief  The Tabular yield condition
*/

class YieldCond_TabularCap : public YieldCondition
{
friend std::ostream& operator<<(std::ostream& out, const YieldCond_TabularCap& yc);
public:

  YieldCond_TabularCap() = delete;
  YieldCond_TabularCap(const YieldCond_TabularCap&) = delete;
  ~YieldCond_TabularCap() = default;
  YieldCond_TabularCap& operator=(const YieldCond_TabularCap&) = delete;

  YieldCond_TabularCap(Uintah::ProblemSpecP& ps);
  YieldCond_TabularCap(const YieldCond_TabularCap* yc);

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["I1_min"] = -d_I1bar_max;
    params["I1_max"] = -d_I1bar_min;
    params["sqrtJ2_max"] = d_sqrtJ2_max;
    params["R"] = d_yield.capEllipticityRatio;
    return params;
  }

  //--------------------------------------------------------------
  // Compute value of yield function
  //--------------------------------------------------------------
  std::pair<double, Util::YieldStatus>
    evalYieldCondition(const ModelStateBase* state) override;
  double evalYieldConditionMax(const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Derivative with respect to the Cauchy stress (\f$\sigma \f$)
  //--------------------------------------------------------------
  void df_dsigma(const Uintah::Matrix3& xi, const ModelStateBase* state,
                      Uintah::Matrix3& df_dsigma) override;

  //--------------------------------------------------------------
  // Compute df/dp  where p = volumetric stress = 1/3 Tr(sigma)
  //--------------------------------------------------------------
  double df_dp(
    const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Compute df/dq  where q = sqrt(3 J_2), J_2 = 2nd invariant deviatoric stress
  //--------------------------------------------------------------
  double df_dq(
    const ModelStateBase* state) override;

  //--------------------------------------------------------------
  // Compute df/deps^p_v 
  //--------------------------------------------------------------
  double df_depsVol(
    const ModelStateBase* state, 
    const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;


  /**
   * Function: computeYieldSurfacePolylinePbarSqrtJ2
   *
   * Purpose: Compute a sequence of points representing the yield surface
   *          in pbar-sqrtJ2 space
   *
   * Inputs:
   *  state_old = old state
   *
   * Returns:
   *   std::vector<Point> 
   */
  Polyline
    computeYieldSurfacePolylinePbarSqrtJ2(const ModelStateBase* state) override;

  /**
   * Function: getUpdatedYieldConditionRange
   *
   * Purpose: Compute range of the yield surface in pbar-sqrtJ2 space
   *
   * Inputs:
   *  std::vector<Point> 
   *
   * Returns:
   *   std::array<double, 3>  = pbar_min, pbar_max, sqrtJ2_max
   */
  std::array<double, 3>
    getYieldConditionRange(const Polyline& yield_surface) override;

  /**
   * Function: getInternalPoint
   *
   * Purpose: Get a point that is inside the yield surface
   *
   * Inputs:
   *  state = state at the current time
   *
   * Returns:
   *   I1 = value of tr(stress) at a point inside the yield surface
   */
  double getInternalPoint(const ModelStateBase* state_old,
                          const ModelStateBase* state_trial) override
  {
    return 0.0;
  }

  /**
   * Functions: getClosestPoint and getClosestPointAndTangent
   *
   * Purpose: Get the point on the yield surface that is closest to a given
   * point (2D)
   *
   * Inputs:
   *  state = current state
   *  px = x-coordinate of point
   *  py = y-coordinate of point
   *
   * Outputs:
   *  cpx = x-coordinate of closest point on yield surface
   *  cpy = y-coordinate of closest point
   *  If Tangent:
   *  tx = x-component of tangent vector
   *  ty = y-component of tangent vector
   *
   * Returns:
   *   true - if the closest point can be found
   *   false - otherwise
   */
  bool getClosestPoint(const ModelStateBase* state, const double& px,
                       const double& py, double& cpx, double& cpy) override;
  bool getClosestPointAndTangent(const ModelStateBase* state, 
                                 const double& px, const double& py, 
                                 double& cpx, double& cpy,
                                 double& tx, double& ty) override;

  //================================================================================
  // Other options below.
  //================================================================================

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dp)
  //--------------------------------------------------------------
  double d2f_dp_depsVol(
    const ModelStateBase* state, const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dp)
  //--------------------------------------------------------------
  double d2f_dp_depsDev(
    const ModelStateBase* state, const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;

  //--------------------------------------------------------------
  // Compute d/depse_v(df/dq)
  //--------------------------------------------------------------
  double d2f_dq_depsVol(
    const ModelStateBase* state, const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;

  //--------------------------------------------------------------
  // Compute d/depse_s(df/dq)
  //--------------------------------------------------------------
  double d2f_dq_depsDev(
    const ModelStateBase* state, const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;

  //--------------------------------------------------------------
  // Compute df/depse_s
  //--------------------------------------------------------------
  double df_depsDev(
    const ModelStateBase* state, const PressureModel* eos,
    const ShearModulusModel* shear,
    const InternalVariableModel* intvar) override;
  // Evaluate the yield function.
  double evalYieldCondition(const double p, const double q, const double dummy0,
                            const double dummy1, double& dummy2) override;

  // Evaluate yield condition (s = deviatoric stress = sigDev
  //                           p = state->pressure
  //                           p_c = state->yieldStress)
  double evalYieldCondition(const Uintah::Matrix3& sigDev,
                            const ModelStateBase* state) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
    with respect to \f$\sigma_{ij}\f$.
  */
  /////////////////////////////////////////////////////////////////////////
  void df_dsigma(const Uintah::Matrix3& stress,
                                const double dummy1, const double dummy2,
                                Uintah::Matrix3& derivative) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the derivative of the yield function \f$(\Phi)\f$
    with respect to \f$s_{ij}\f$.

    This is for the associated flow rule with \f$s_{ij}\f$ being
    the deviatoric stress.
  */
  /////////////////////////////////////////////////////////////////////////
  void df_dsigmaDev(const Uintah::Matrix3& stress,
                                   const double dummy1, const double dummy2,
                                   Uintah::Matrix3& derivative) override;

  /*! Derivative with respect to the \f$xi\f$ where \f$\xi = s - \beta \f$
    where \f$s\f$ is deviatoric part of Cauchy stress and
    \f$\beta\f$ is the backstress */
  void df_dxi(const Uintah::Matrix3& xi, const ModelStateBase* state,
                   Uintah::Matrix3& df_xi) override;

  /* Derivative with respect to \f$ s \f$ and \f$ \beta \f$ */
  void df_dsigmaDev_dbeta(const Uintah::Matrix3& xi,
                           const ModelStateBase* state, Uintah::Matrix3& df_ds,
                           Uintah::Matrix3& df_dbeta) override;

  /*! Derivative with respect to the plastic strain (\f$\epsilon^p \f$)*/
  double df_dplasticStrain(const Uintah::Matrix3& xi, const double& d_sigy_dep,
                     const ModelStateBase* state) override;

  /*! Derivative with respect to the porosity (\f$\epsilon^p \f$)*/
  double df_dporosity(const Uintah::Matrix3& xi,
                      const ModelStateBase* state) override;

  /*! Compute h_alpha  where \f$d/dt(ep) = d/dt(gamma)~h_{\alpha}\f$ */
  double eval_h_alpha(const Uintah::Matrix3& xi,
                      const ModelStateBase* state) override;

  /*! Compute h_phi  where \f$d/dt(phi) = d/dt(gamma)~h_{\phi}\f$ */
  double eval_h_phi(const Uintah::Matrix3& xi, const double& factorA,
                    const ModelStateBase* state) override;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the elastic-plastic tangent modulus.
  */
  /////////////////////////////////////////////////////////////////////////
  void computeElasPlasTangentModulus(
    const Uintah::TangentModulusTensor& Ce, const Uintah::Matrix3& sigma,
    double sigY, double dsigYdep, double porosity, double voidNuclFac,
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
  inline double computePlasticStrainFactor(double sigma_f_sigma,
                                           double sigma_Y);

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
  void computeTangentModulus(const Uintah::TangentModulusTensor& Ce,
                             const Uintah::Matrix3& f_sigma, double f_q1,
                             double h_q1, Uintah::TangentModulusTensor& Cep);


  /* Compute points on the cap */
  void computeCapPoints(double X_bar, Polyline& p_q_all);

private:
  /**
   *  These are the parameters that are read from the input file
   */
  struct YieldFunctionParameters
  {
    TabularData table;
    double capEllipticityRatio;
    YieldFunctionParameters() = default;
    YieldFunctionParameters(Uintah::ProblemSpecP& ps)
      : table(ps)
    {
      table.setup();
      ps->require("cap_ellipticity_ratio", capEllipticityRatio);
    }
    YieldFunctionParameters(const YieldFunctionParameters& yf)
    {
      table = yf.table;
      capEllipticityRatio = yf.capEllipticityRatio;
    }
    YieldFunctionParameters& operator=(const YieldFunctionParameters& yf)
    {
      if (this != &yf) {
        table = yf.table;
        capEllipticityRatio = yf.capEllipticityRatio;
      }
      return *this;
    }
  };

  YieldFunctionParameters d_yield;
  double d_I1bar_min;
  double d_I1bar_max;
  double d_sqrtJ2_max;
  Polyline  d_polyline;
  std::vector<Uintah::Vector> d_normals;

  /* Some helpers and checks */
  void checkInputParameters();
  void setYieldConditionRange();
  void saveAsPolyline();
  void computeNormals();

  /* Find the closest point */
  Uintah::Point getClosestPoint(const Polyline& polyline,
                                const double& p_bar, const double& sqrtJ2);
  Uintah::Point getClosestPointTable(const ModelState_TabularCap* state,
                                     const Uintah::Point& z_r_pt);
  Uintah::Point getClosestPointSpline(const ModelState_TabularCap* state,
                                      const Uintah::Point& z_r_pt);
  std::tuple<Uintah::Point, Uintah::Vector>
    getClosestPointSplineNewton(const ModelState_TabularCap* state, 
                                const Uintah::Point& z_r_pt);

  /* Compute the height of the elliptical cap */
  double computeEllipseHeight(const Polyline& p_q_points, double p_cap);

  /* Convert yield function data to z_rprime coordinates */
  void convertToZRprime(const double& sqrtKG, const Polyline& p_q_points, 
                        Polyline& z_r_points) const;

  /* Check the location of the trial point by computing distances */
  Util::YieldStatus checkClosestPointDistance(const ModelState_TabularCap* state);
};

} // End namespace Uintah

#endif // __TABULAR_CAP_YIELD_CONDITION_MODEL_H__
