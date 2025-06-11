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

#ifndef __nullptr_PRESSURE_TEMPLATED_MODEL_H__
#define __nullptr_PRESSURE_TEMPLATED_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EquationOfStateT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_MetalT.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class EOS_NullT

  \brief A dummy EOS

*/
////////////////////////////////////////////////////////////////////////////

class EOS_NullT : public EquationOfStateT<EOS_NullT, ModelState_MetalT>
{

public:
  EOS_NullT([[maybe_unused]] Uintah::ProblemSpecP& ps) {}
  EOS_NullT([[maybe_unused]] const EOS_NullT* cm) {}
  EOS_NullT&
  operator=(const EOS_NullT& cm) = delete;
  ~EOS_NullT()                   = default;

  void
  l_outputProblemSpec([[maybe_unused]] Uintah::ProblemSpecP& ps)
  {
  }

  /*! Get parameters */
  std::map<std::string, double>
  l_getParameters() const
  {
    std::map<std::string, double> params;
    return params;
  }

  /////////////////////////////////////////////////////////////////////////
  /*! Calculate the pressure using a equation of state */
  /////////////////////////////////////////////////////////////////////////
  double
  l_computePressure([[maybe_unused]] const Uintah::MPMMaterial* matl,
                    [[maybe_unused]] const ModelState_MetalT* state,
                    [[maybe_unused]] const Uintah::Matrix3& deformGrad,
                    [[maybe_unused]] const Uintah::Matrix3& rateOfDeformation,
                    [[maybe_unused]] const double& delT)
  {
    return 0.0;
  }

  double
  l_eval_dp_dJ([[maybe_unused]] const Uintah::MPMMaterial* matl,
               [[maybe_unused]] const double& delF,
               [[maybe_unused]] const ModelState_MetalT* state)
  {
    return 0.0;
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_v
      where epse_v = tr(epse)
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  l_computeDpDepse_v([[maybe_unused]] const ModelState_MetalT* state) const
  {
    return 0.0;
  }

  ////////////////////////////////////////////////////////////////////////
  /*! Calculate the derivative of p with respect to epse_s
      where epse_s = sqrt{2}{3} ||ee||
            ee = epse - 1/3 tr(epse) I
            epse = total elastic strain */
  ////////////////////////////////////////////////////////////////////////
  double
  l_computeDpDepse_s([[maybe_unused]] const ModelState_MetalT* state) const
  {
    return 0.0;
  }

  // Calculate rate of temperature change due to compression/expansion
  double
  l_computeIsentropicTemperatureRate([[maybe_unused]] const double T,
                                     [[maybe_unused]] const double rho_0,
                                     [[maybe_unused]] const double rho_cur,
                                     [[maybe_unused]] const double Dtrace)
  {
    return 0.0;
  }

  // Compute pressure (option 1)
  double
  l_computePressure([[maybe_unused]] const double& rho_orig, [[maybe_unused]] const double& rho_cur)
  {
    return 0.0;
  }

  // Compute pressure (option 2)
  void
  l_computePressure([[maybe_unused]] const double& rho_orig,
                    [[maybe_unused]] const double& rho_cur,
                    [[maybe_unused]] double& pressure,
                    [[maybe_unused]] double& dp_drho,
                    [[maybe_unused]] double& csquared)
  {
  }

  // Compute bulk modulus
  void
  setInitialBulkModulus()
  {
  }
  double
  l_computeInitialBulkModulus()
  {
    return 0.0;
  }
  double
  l_computeBulkModulus([[maybe_unused]] const ModelState_MetalT* state)
  {
    return 0.0;
  }
  double
  l_computeBulkModulus([[maybe_unused]] const double& rho_orig, [[maybe_unused]] const double& rho_cur)
  {
    return 0.0;
  }

  // Compute strain energy
  double
  l_computeStrainEnergy([[maybe_unused]] const ModelState_MetalT* state)
  {
    return 0.0;
  }
  double
  l_computeStrainEnergy([[maybe_unused]] const double& rho_orig, [[maybe_unused]] const double& rho_cur)
  {
    return 0.0;
  }

  // Compute density given pressure
  double
  l_computeDensity([[maybe_unused]] const double& rho_orig, [[maybe_unused]] const double& pressure)
  {
    return 0.0;
  }

  double
  l_computeElasticVolumetricStrain([[maybe_unused]] const double& pp, [[maybe_unused]] const double& p0)
  {
    return 0.0;
  }
  double
  l_computeExpElasticVolumetricStrain([[maybe_unused]] const double& pp, [[maybe_unused]] const double& p0)
  {
    return 0.0;
  }
  double
  l_computeDerivExpElasticVolumetricStrain([[maybe_unused]] const double& pp,
                                           [[maybe_unused]] const double& p0,
                                           [[maybe_unused]] double& exp_eps_e_v)
  {
    return 0.0;
  }
};

} // End namespace Vaango

#endif // __nullptr_PRESSURE_TEMPLATED_MODEL_H__
