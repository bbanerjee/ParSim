/*
 * The MIT License
 *
 * Copyright (c) 2015-2106 Parresia Research Limited, New Zealand
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

#ifndef ___ELASTIC_MODULI_ARENA_MODEL_H__
#define ___ELASTIC_MODULI_ARENA_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/AirEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/GraniteEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/WaterEOS.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <limits>

namespace Vaango {

/*! \class ElasticModuli_Arena
 *  \brief The updated elasticity for Arenisca3
 *  \author Biswajit Banerjee,
 *          Caveat below by Michael Homel.
 *
 *  Purpose:
 *  Compute the nonlinear elastic tangent stiffness as a function of the
 * pressure
 *  plastic strain, and fluid parameters.
 *
 * Caveat:
 *  To be thermodynamically consistent, the shear modulus in an isotropic model
 *  must be constant, but the bulk modulus can depend on pressure.  However,
 * this
 *  leads to a Poisson's ratio that approaches 0.5 at high pressures, which is
 *  inconsistent with experimental data for the Poisson's ratio, inferred from
 * the
 *  Young's modulus.  Induced anisotropy is likely the cause of the discrepency,
 *  but it may be better to allow the shear modulus to vary so the Poisson's
 * ratio
 *  remains reasonable.
 *
 *  If the user has specified a nonzero value of nu1 and nu2, the shear modulus
 * will
 *  vary with pressure so the drained Poisson's ratio transitions from nu1 to
 * nu1+nu2 as
 *  the bulk modulus varies from b0 to b0+b1.  The fluid model further affects
 * the
 *  bulk modulus, but does not alter the shear modulus, so the pore fluid does
 *  increase the Poisson's ratio.
 *
 */
class ElasticModuli_Arena : public ElasticModuliModel
{

private:
  /* Tangent bulk modulus parameters */
  struct BulkModulusParameters
  {
    double b0;
    double b1;
    double b2;
    double b3;
    double b4;
  };

  /* Tangent shear modulus parameters */
  struct ShearModulusParameters
  {
    double G0;
    double nu1;
    double nu2;
  };

  BulkModulusParameters d_bulk;
  ShearModulusParameters d_shear;

  /* Tangent bulk modulus models for air, water, granite */
  AirEOS d_air;
  WaterEOS d_water;
  GraniteEOS d_granite;

  void checkInputParameters();

  void computeDrainedModuli(
    const double& I1_eff_bar, const double& ev_p_bar, double& KK,
    double& GG) const override; 

  void computePartialSaturatedModuli(
    const double& I1_eff_bar, const double& pw_bar, const double& ev_p_bar,
    const double& phi, const double& S_w, double& KK,
    double& GG) const override; 


public:
  /*! Construct a constant elasticity model. */
  ElasticModuli_Arena(Uintah::ProblemSpecP& ps);

  /*! Construct a copy of constant elasticity model. */
  ElasticModuli_Arena(const ElasticModuli_Arena& smm) = default;

  /*! Destructor of constant elasticity model.   */
  ~ElasticModuli_Arena() override = default;

  ElasticModuli_Arena& operator=(const ElasticModuli_Arena& smm) = delete;

  // Polymorphic cloning method (preferred over operator=(const T*))
  [[nodiscard]] std::unique_ptr<ElasticModuliModel>
  clone() const
  {
    return std::make_unique<ElasticModuli_Arena>(*this);
  }

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["b0"] = d_bulk.b0;
    params["b1"] = d_bulk.b1;
    params["b2"] = d_bulk.b2;
    params["b3"] = d_bulk.b3;
    params["b4"] = d_bulk.b4;
    params["G0"] = d_shear.G0;
    params["nu1"] = d_shear.nu1;
    params["nu2"] = d_shear.nu2;
    return params;
  }

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const override;
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const override;

  ElasticModuli getElasticModuliLowerBound() const override
  {
    return getInitialElasticModuli();
  }
  ElasticModuli getElasticModuliUpperBound() const override
  {
    return ElasticModuli(std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
  }

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> computeDModuliDIntVar(const ModelStateBase* state) const override;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const override;
};
} // End namespace Vaango

#endif // __ELASTIC_MODULI_ARENA_MODEL_H__
