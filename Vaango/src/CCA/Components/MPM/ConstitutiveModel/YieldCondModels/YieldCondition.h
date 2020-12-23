/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __BB_YIELD_CONDITION_H__
#define __BB_YIELD_CONDITION_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulusModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>
#include <Core/Grid/Task.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <Core/ProblemSpec/ProblemSpec.h>

namespace Vaango {

using ParameterDict = std::map<std::string, double>;
using Polyline      = std::vector<Uintah::Point>;

/*! \class YieldCondition
 *  \brief A generic wrapper for various yield conditions
 *  \author Biswajit Banerjee,
 *  \author C-SAFE and Department of Mechanical Engineering,
 *  \author University of Utah.
 *  \warning Mixing and matching yield conditions with damage and plasticity
 *           models should be done with care.  No checks are provided to stop
 *           the user from using the wrong combination of models.
 *
 * Provides an abstract base class for various yield conditions used
 * in the plasticity and damage models
*/
class YieldCondition
{
public:
  //! Construct a yield condition.
  /*! This is an abstract base class. */
  YieldCondition();

  //! Destructor of yield condition.
  /*! Virtual to ensure correct behavior */
  virtual ~YieldCondition();

  virtual void
  outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;

  virtual void
  addParticleState(std::vector<const VarLabel*>& from,
                   std::vector<const VarLabel*>& to){};

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Get/compute the parameters of the yield condition model
   */
  /////////////////////////////////////////////////////////////////////////
  virtual std::map<std::string, double>
  getParameters() const = 0;
  virtual void
  computeModelParameters(double factor){};

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate the yield function \f$(f)\f$.
    If \f$f \le 0\f$ the state is elastic.
    If \f$f > 0\f$ the state is plastic.
  */
  /////////////////////////////////////////////////////////////////////////
  virtual std::pair<double, Util::YieldStatus>
  evalYieldCondition(const ModelStateBase* state) = 0;

  virtual double
  computeYieldFunction(const ModelStateBase* state) const = 0;

  virtual double
  evalYieldCondition(const Uintah::Matrix3& stress,
                     const ModelStateBase* state) = 0;

  virtual double
  evalYieldConditionMax(const ModelStateBase* state) = 0;

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
  virtual Uintah::Matrix3
  df_dsigma(const ModelStateBase* state) = 0;

  virtual Uintah::Matrix3
  df_dsigma(const Uintah::Matrix3& stress, const ModelStateBase* state) = 0;

  virtual Uintah::Matrix3
  df_dxi(const Uintah::Matrix3& stress, const ModelStateBase* state) = 0;

  virtual std::pair<Uintah::Matrix3, Uintah::Matrix3>
  df_dsigmaDev_dbeta(const Uintah::Matrix3& stress,
                     const ModelStateBase* state) = 0;

  virtual double
  df_dp(const ModelStateBase* state) = 0;

  virtual double
  df_dq(const ModelStateBase* state) = 0;

  /* Derivatives of the yield function wrt internal variables */
  virtual void
  df_dintvar(const ModelStateBase* state,
             MetalIntVar& df_dintvar) const 
  {
    std::ostringstream out;
    out << "**ERROR** df_dintvar MetalIntVar argument should not be "
        << "called by the yield condition model currently in use.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  virtual void
  df_dintvar(const ModelStateBase* state,
             ArenaIntVar& df_dintvar) const
  {
    std::ostringstream out;
    out << "**ERROR** df_dintvar ArenaIntVar argument should not be "
        << "called by the yield condition model currently in use.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  virtual void
  df_dintvar(const ModelStateBase* state,
             BorjaIntVar& df_dintvar) const
  {
    std::ostringstream out;
    out << "**ERROR** df_dintvar BorjaIntVar argument should not be "
        << "called by the yield condition model currently in use.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  virtual void
  df_dintvar(const ModelStateBase* state,
             SoilBrannonIntVar& df_dintvar) const
  {
    std::ostringstream out;
    out << "**ERROR** df_dintvar SoilBrannonIntVar argument should not be "
        << "called by the yield condition model currently in use.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  virtual void
  df_dintvar(const ModelStateBase* state,
             TabularCapIntVar& df_dintvar) const
  {
    std::ostringstream out;
    out << "**ERROR** df_dintvar TabularCapIntVar argument should not be "
        << "called by the yield condition model currently in use.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  /* Compute d/depse_v(df/dp) */
  virtual double
  d2f_dp_depsVol(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) = 0;

  /* Compute d/depse_s(df/dp) */
  virtual double
  d2f_dp_depsDev(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) = 0;

  /* Compute d/depse_v(df/dq) */
  virtual double
  d2f_dq_depsVol(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) = 0;

  /* Compute d/depse_s(df/dq) */
  virtual double
  d2f_dq_depsDev(const ModelStateBase* state,
                 const MPMEquationOfState* eos,
                 const ShearModulusModel* shear) = 0;

  /* Compute df/depse_v */
  virtual double
  df_depsVol(const ModelStateBase* state,
             const MPMEquationOfState* eos,
             const ShearModulusModel* shear) = 0;

  /* Compute df/depse_s */
  virtual double
  df_depsDev(const ModelStateBase* state,
             const MPMEquationOfState* eos,
             const ShearModulusModel* shear) = 0;

  /////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the elastic-plastic tangent modulus.
  */
  /////////////////////////////////////////////////////////////////////////
  virtual void
  computeElasPlasTangentModulus(const Uintah::TangentModulusTensor& Ce,
                                const Uintah::Matrix3& sigma,
                                double sigY,
                                double dsigYdV,
                                double porosity,
                                double voidNuclFac,
                                Uintah::TangentModulusTensor& Cep) = 0;

  /*! Compute continuum elastic-plastic tangent modulus.
     df_dsigma = r */
  void
  computeElasPlasTangentModulus(const Uintah::Matrix3& r,
                                const Uintah::Matrix3& df_ds,
                                const Uintah::Matrix3& h_beta,
                                const Uintah::Matrix3& df_dbeta,
                                const double& h_alpha,
                                const double& df_dep,
                                const double& h_phi,
                                const double& df_phi,
                                const double& J,
                                const double& dp_dJ,
                                const ModelStateBase* state,
                                Uintah::TangentModulusTensor& Cep);

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
  virtual Polyline
  computeYieldSurfacePolylinePbarSqrtJ2(const ModelStateBase* state_old)
  {
    Polyline dummy;
    return dummy;
  }

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
  virtual std::array<double, 3>
  getYieldConditionRange(const Polyline& yield_surface)
  {
    std::array<double, 3> dummy;
    return dummy;
  }

  /**
   * Function: getInternalPoint
   *
   * Purpose: Get a point that is inside the yield surface
   *
   * Inputs:
   *  state_old = old state
   *  state_new = new state
   *
   * Returns:
   *   I1 = value of tr(stress) at a point inside the yield surface
   */
  virtual double
  getInternalPoint(const ModelStateBase* state_old,
                   const ModelStateBase* state_new) = 0;

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
  virtual bool
  getClosestPoint(const ModelStateBase* state,
                  const double& px,
                  const double& py,
                  double& cpx,
                  double& cpy)
  {
    return false;
  }
  virtual bool
  getClosestPointAndTangent(const ModelStateBase* state,
                            const double& px,
                            const double& py,
                            double& cpx,
                            double& cpy,
                            double& tx,
                            double& ty)
  {
    return false;
  }
  virtual bool
  getClosestPointAndTangent(const ModelStateBase* state,
                            const Polyline& z_r_table,
                            const Util::PolylineKDTree& z_r_index, 
                            const double& px,
                            const double& py,
                            double& cpx,
                            double& cpy,
                            double& tx,
                            double& ty)
  {
    return false;
  }

  /**
   * These are needed for keeping track of point-to-point material variability
   */
  virtual void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* matl,
                                const PatchSet* patch) const {};

  virtual void
  initializeLocalVariables(const Patch* patch,
                           ParticleSubset* pset,
                           DataWarehouse* new_dw,
                           constParticleVariable<double>& pVolume){};

  virtual void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches) const {};

  virtual void
  copyLocalVariables(ParticleSubset* pset,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw){};

  virtual std::vector<std::string>
  getLocalVariableLabels() const
  {
    std::vector<std::string> pYieldParamLabels;
    pYieldParamLabels.emplace_back("None");
    return pYieldParamLabels;
  }

  virtual std::vector<Uintah::constParticleVariable<double>>
  getLocalVariables(Uintah::ParticleSubset* pset, Uintah::DataWarehouse* old_dw)
  {
    constParticleVariable<double> pNull;
    std::vector<constParticleVariable<double>> pYieldParams;
    pYieldParams.emplace_back(pNull);
    return pYieldParams;
  }

  /**
   *  This is used to scale the yield parameters
   */
  virtual void
  updateLocalVariables(ParticleSubset* pset,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw,
                       constParticleVariable<double>& pCoherence_old,
                       const ParticleVariable<double>& pCoherence_new){};

};

} // End namespace Vaango

#endif // __BB_YIELD_CONDITION_H__
