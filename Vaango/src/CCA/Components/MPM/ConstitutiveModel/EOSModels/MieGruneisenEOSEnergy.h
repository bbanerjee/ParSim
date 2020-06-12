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

#ifndef __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__
#define __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__

#include "MPMEquationOfState.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

////////////////////////////////////////////////////////////////////////////
/*!
  \class MieGruneisenEOSEnergy

  \brief A Mie-Gruneisen type equation of state model

  \author Biswajit Banerjee \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n

  Reference:

  Zocher, Maudlin, Chen, Flower-Maudlin, 2000,
  European Congress on Computational Methods in Applied Science
  and Engineering, ECOMAS 2000, Barcelona)


  The equation of state is given by
  \f[
  p = \frac{\rho_0 C_0^2 \zeta
            \left[1 + \left(1-\frac{\Gamma_0}{2}\right)\zeta\right]}
           {\left[1 - (S_{\alpha} - 1) \zeta\right]^2 + \Gamma_0 C_p T}
  \f]
  where
  \f$ p\f$ = pressure \n
  \f$ C_0 \f$= bulk speed of sound \n
  \f$ \zeta = (\rho/\rho_0 - 1)\f$ \n
  where \f$\rho\f$ = current density \n
  \f$\rho_0\f$ = initial density \n
  \f$ E\f$ = internal energy = \f$C_p T\f$ \n
  where \f$C_p\f$ = specfic heat at constant pressure \n
  \f$T\f$ = temperature \n
  \f$\Gamma_0\f$ = Gruneisen's gamma at reference state \n
  \f$S_{\alpha}\f$ = linear Hugoniot slope coefficient

  \ Modified by Jim Guilkey to be energy based.

*/
////////////////////////////////////////////////////////////////////////////

class MieGruneisenEOSEnergy : public MPMEquationOfState
{
public:
  struct CMData
  {
    double C_0;
    double Gamma_0;
    double S_1;
    double S_2;
    double S_3;
  };

  // constructors
  MieGruneisenEOSEnergy(ProblemSpecP& ps);
  MieGruneisenEOSEnergy(const MieGruneisenEOSEnergy* cm);
  MieGruneisenEOSEnergy& operator=(const MieGruneisenEOSEnergy& cm) = delete;

  // Special operator for computing internal energy
  double operator()(double eta) const;

  // destructor
  ~MieGruneisenEOSEnergy() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  std::map<std::string, double> getParameters() const override;

  /////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure using a equation of state */
  /////////////////////////////////////////////////////////////////////////
  double computePressure(const MPMMaterial* matl, const ModelStateBase* state,
                         const Matrix3& deformGrad,
                         const Matrix3& rateOfDeformation,
                         const double& delT) override;

  // Calculate rate of temperature change due to compression/expansion
  double computeIsentropicTemperatureRate(const double T, const double rho_0,
                                          const double rho_cur,
                                          const double Dtrace) override;

  double eval_dp_dJ(const MPMMaterial* matl, const double& delF,
                    const ModelStateBase* state) override;

  // Compute pressure (option 1)
  double computePressure(const double& rho_orig,
                         const double& rho_cur) override;

  // Compute pressure (option 2)
  void computePressure(const double& rho_orig, const double& rho_cur,
                       double& pressure, double& dp_drho,
                       double& csquared) override;

  // Compute bulk modulus
  double computeBulkModulus(const double& rho_orig,
                            const double& rho_cur) override;
  double computeBulkModulus(const ModelStateBase* state) override;

  // Compute strain energy
  double computeStrainEnergy(const double& rho_orig,
                             const double& rho_cur) override;

  // Compute density given pressure
  double computeDensity(const double& rho_orig,
                        const double& pressure) override;

private:
  typedef double (MieGruneisenEOSEnergy::*pFuncPtr)(const double&,
                                                    const double&);
  typedef double (MieGruneisenEOSEnergy::*dpdJFuncPtr)(const double&,
                                                       const double&);

  // Find root of p(eta) - p0 = 0 using Ridder's method
  double findEtaRidder(pFuncPtr pFunc, const double& rho_orig, const double& p0,
                       double& etamin, double& etamax, const double& tolerance,
                       const int& maxIter);

  // Find root of p(eta) - p0 = 0 using Newton's method
  double findEtaNewton(pFuncPtr pFunc, dpdJFuncPtr dpdJFunc,
                       const double& rho_orig, const double& p0,
                       const double& J0, const double& tolerance,
                       const int& maxIter);

  // Compute p for compressive volumetric deformations
  double pCompression(const double& rho_orig, const double& eta);

  // Compute dp/dJ for compressive volumetric deformations
  double dpdJCompression(const double& rho_orig, const double& eta);

  // Compute p for tensile volumetric deformations
  double pTension(const double& rho_orig, const double& eta);

  // Compute dp/dJ for tensile volumetric deformations
  double dpdJTension(const double& rho_orig, const double& eta);

private:
  CMData d_const;

};

} // End namespace Uintah

#endif // __MIE_GRUNEISEN_EOS_ENERGY_MODEL_H__
