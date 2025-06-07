/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/IsoNonlinHypoelastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <vector>

using namespace Vaango;
using Uintah::Matrix3;
using Uintah::DeformationState;

IsoNonlinHypoelastic::IsoNonlinHypoelastic() 
{
  d_elastic = nullptr;
  d_intvar = nullptr;
}

IsoNonlinHypoelastic::IsoNonlinHypoelastic(const ElasticModuliModel* elastic,
                                           const InternalVariableModel* intvar)
{
  d_elastic = elastic;
  d_intvar = intvar;
}

Matrix3
IsoNonlinHypoelastic::computeStress(double delT,
                                    const Matrix3& stress_old,
                                    const DeformationState* deformState,
                                    const ModelStateBase* modelState) const
{
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(modelState);
  double kappa = moduli.bulkModulus;
  double mu = moduli.shearModulus;
  double lambda = kappa - (2.0 / 3.0) * mu;

  Matrix3 D = deformState->D;
  Matrix3 stress_rate = 
    Vaango::Util::Identity * (lambda * D.Trace()) + D * (2.0 * mu);

  //std::cout << "delT = " << delT << "\n";
  //std::cout << "sig_old = " << stress_old << "\n";
  //std::cout << "stress_rate = " << stress_rate << "\n";
  //std::cout << "K = " << kappa << " G = " << mu << "\n";
  //std::cout << "stress_inc = " << stress_rate * delT << "\n";

  Matrix3 stress_new = stress_old + stress_rate * delT;
  //std::cout << "sig_new = " << stress_new << "\n";
  return stress_new;
}

std::vector<Matrix3>
IsoNonlinHypoelastic::computeDStressDIntVar(double delT,
                                            const std::vector<Matrix3>& derivStress_old,
                                            const DeformationState* deformState,
                                            const ModelStateBase* modelState) const
{
  std::vector<ElasticModuli> derivModuli = d_elastic->computeDModuliDIntVar(modelState);
  std::vector<Matrix3> dStress;
  int intVarID = 0;
  for (auto deriv : derivModuli) {
    double derivKappa = deriv.bulkModulus;
    double derivMu = deriv.shearModulus;
    double derivLambda = derivKappa - (2.0 / 3.0) * derivMu;

    Matrix3 D = deformState->D;
    Matrix3 derivStress_rate = 
      Vaango::Util::Identity * (derivLambda * D.Trace()) + D * (2.0 * derivMu);
    
    Matrix3 derivStress_new = derivStress_old[intVarID] + derivStress_rate * delT;
    dStress.push_back(derivStress_new);
    ++intVarID;
  }
  return dStress;
}

std::pair<Matrix3, std::vector<Matrix3>>
IsoNonlinHypoelastic::computeStressAndDStressDIntVar(double delT,
                                                     const Matrix3& stress_old,
                                                     const std::vector<Matrix3>& derivStress_old,
                                                     const DeformationState* deformState,
                                                     const ModelStateBase* modelState) const
{
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(modelState);
  double kappa = moduli.bulkModulus;
  double mu = moduli.shearModulus;
  double lambda = kappa - (2.0 / 3.0) * mu;

  Matrix3 D = deformState->D;
  Matrix3 stress_rate = 
    Vaango::Util::Identity * (lambda * D.Trace()) + D * (2.0 * mu);

  Matrix3 stress_new = stress_old + stress_rate * delT;

  std::vector<ElasticModuli> derivModuli = d_elastic->computeDModuliDIntVar(modelState);
  std::vector<Matrix3> derivStress_new;
  int intVarID = 0;
  for (auto deriv : derivModuli) {
    double derivKappa = deriv.bulkModulus;
    double derivMu = deriv.shearModulus;
    double derivLambda = derivKappa - (2.0 / 3.0) * derivMu;

    Matrix3 dStress_rate = 
      Vaango::Util::Identity * (derivLambda * D.Trace()) + D * (2.0 * derivMu);
    
    Matrix3 dStress_new = derivStress_old[intVarID] + dStress_rate * delT;
    derivStress_new.push_back(dStress_new);
    ++intVarID;
  }
  return std::make_pair(stress_new, derivStress_new);
}

Tensor::Matrix6Mandel
IsoNonlinHypoelastic::computeElasticTangentModulus(const ModelStateBase* state) const
{
  return d_elastic->computeElasticTangentModulus(state);
}
