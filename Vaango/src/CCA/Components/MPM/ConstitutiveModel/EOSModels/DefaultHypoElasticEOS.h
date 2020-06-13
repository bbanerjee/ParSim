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

#ifndef __DEFAULT_HYPOELASTIC_EOS_MODEL_H__
#define __DEFAULT_HYPOELASTIC_EOS_MODEL_H__

#include "MPMEquationOfState.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

////////////////////////////////////////////////////////////////////////////
/*!
  \class DefaultHypoElasticEOS
  \brief Not really an equation of state but just an isotropic
  hypoelastic pressure calculator based on bulk modulus.
  \author Biswajit Banerjee, \n
  C-SAFE and Department of Mechanical Engineering, \n
  University of Utah \n

  The equation of state is given by
  \f[
  p = Tr(D) K \Delta T
  \f]
  where \n
  \f$p\f$ = pressure\n
  \f$D\f$ = rate of deformation tensor\n
  \f$K\f$ = bulk modulus\n
  \f$\Delta T\f$ = time increment
*/
////////////////////////////////////////////////////////////////////////////

class DefaultHypoElasticEOS : public MPMEquationOfState
{
public:
  // constructors
  DefaultHypoElasticEOS(); // This constructor is used when there is
                           // no equation_of_state tag in the input
                           // file  ** WARNING **
  DefaultHypoElasticEOS(ProblemSpecP& ps);
  DefaultHypoElasticEOS(const DefaultHypoElasticEOS* cm);
  DefaultHypoElasticEOS& operator=(const DefaultHypoElasticEOS& cm) = delete;

  // destructor
  ~DefaultHypoElasticEOS() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  std::map<std::string, double> getParameters() const override;

  //////////
  // Calculate the pressure using a equation of state
  double computePressure(const MPMMaterial* matl, const ModelStateBase* state,
                         const Matrix3& deformGrad,
                         const Matrix3& rateOfDeformation,
                         const double& delT) override;

  double eval_dp_dJ(const MPMMaterial* matl, const double& detF,
                    const ModelStateBase* state) override;

  // Compute pressure (option 1)
  double computePressure(const double& rho_orig,
                         const double& rho_cur) override;

  // Compute pressure (option 2)
  void computePressure(const double& rho_orig, const double& rho_cur,
                       double& pressure, double& dp_drho,
                       double& csquared) override;

  // Compute bulk modulus
  double computeInitialBulkModulus() override;
  double computeBulkModulus(const double& rho_orig,
                            const double& rho_cur) override;
  double computeBulkModulus(const ModelStateBase* state) override;

  // Compute strain energy
  double computeStrainEnergy(const ModelStateBase* state) override;
  double computeStrainEnergy(const double& rho_orig,
                             const double& rho_cur) override;

  // Compute density given pressure
  double computeDensity(const double& rho_orig,
                        const double& pressure) override;

  double computeDpDepse_v(const Vaango::ModelStateBase*) const override;
  double computeDpDepse_s(const Vaango::ModelStateBase*) const override;
  double computeElasticVolumetricStrain(const double& pp,
                                        const double& p0) override;

  double computeExpElasticVolumetricStrain(const double& pp,
                                           const double& p0) override;
  double computeDerivExpElasticVolumetricStrain(const double& pp,
                                                const double& p0,
                                                double& exp_eps_e_v) override;
};

} // End namespace Uintah

#endif // __DEFAULT_HYPOELASTIC_EOS_MODEL_H__
