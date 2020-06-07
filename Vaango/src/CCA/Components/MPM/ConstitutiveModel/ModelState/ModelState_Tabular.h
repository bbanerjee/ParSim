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

#ifndef __MODEL_STATE_TABULAR_H__
#define __MODEL_STATE_TABULAR_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Default.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_Tabular
  \brief A structure that stores the state data that is specialized for
         the Tabular plasticity model.
         ** Derived from PlasticityState:ModelState
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_Tabular : public ModelState_Default
{

public:
  static const Uintah::Matrix3 Identity;
  static const double sqrtTwo;
  static const double sqrtThree;

  Uintah::long64 particleID;

  double I1;      // I1= Tr(sigma)
  double J2;
  double sqrt_J2; // sqrt(J2)
  double zz;      // Lode coordinate 'z'
  double rr;      // Lode coordinate 'r'
  double ep_v;      // ep_v = Tr(ep) : Volumetric part of the plastic strain
  double ep_eq;     // The equivalent plastic strain computed from the current
                    // plastic strain (This quantity can decrease)
  double ep_cum_eq; // The cumulative equivalent plastic strain
                    // (This quantity increases monotonically)

  Uintah::Matrix3 stressTensor; 
  Uintah::Matrix3 deviatoricStressTensor; 
  Uintah::Matrix3 elasticStrainTensor; 
  Uintah::Matrix3 plasticStrainTensor;

  ModelState_Tabular();

  ModelState_Tabular(const ModelState_Tabular& state) = default;
  ModelState_Tabular(const ModelState_Tabular* state);

  ~ModelState_Tabular() = default;

  ModelState_Tabular& operator=(const ModelState_Tabular& state) = default;
  ModelState_Tabular* operator=(const ModelState_Tabular* state);

  virtual 
  size_t numStateVar() const override
  {
    auto numBase = ModelState_Default::numStateVar();
    auto numThis = 16u;
    return numBase + numThis;
  }

  void updateStressInvariants();
  void updatePlasticStrainInvariants();

  friend std::ostream& operator<<(std::ostream& os,
                                  const ModelState_Tabular& state)
  {
    os << "ParticleID = " << state.particleID << "\n"
       << "\t sigma = " << state.stressTensor << "\n"
       << "\t dev(sigma) = " << state.deviatoricStressTensor << "\n"
       << "\t I1 = " << state.I1 << ", J2 = " << state.J2
       << ", sqrt_J2 = " << state.sqrt_J2 << ", r = " << state.rr
       << ", z = " << state.zz << "\n"
       << "\t eps_e = " << state.elasticStrainTensor << "\n"
       << "\t eps_p = " << state.plasticStrainTensor << "\n"
       << "\t evp = " << state.ep_v << " ep_eq = " << state.ep_eq 
       << "\t ep_cum_eq = " << state.ep_cum_eq << "\n"
       << "\t K = " << state.bulkModulus << ", G = " << state.shearModulus
       << "\n";
    return os;
  }
};

} // End namespace Uintah

#endif // __MODEL_STATE_TABULAR_H__
