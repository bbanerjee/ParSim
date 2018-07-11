/*
 * The MIT License
 *
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

#ifndef __MODEL_STATE_TABULAR_CAP_H__
#define __MODEL_STATE_TABULAR_CAP_H__

#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>
#include <vector>
#include <iterator>

namespace Vaango {

  using Polyline = std::vector<Uintah::Point>;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_TabularCap
  \brief A structure that stores the state data that is specialized for
         the Tabular plasticity model with cap.
         ** Derived from PlasticityState:ModelState
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_TabularCap : public ModelState_Tabular
{

public:
  double capX;             // Hydrostatic strength in I1 space
  Polyline yield_f_pts;    // Polyline representing yield function with cap
                           // in pbar-q space

  ModelState_TabularCap();

  ModelState_TabularCap(const ModelState_TabularCap& state) = default;
  ModelState_TabularCap(const ModelState_TabularCap* state);

  ~ModelState_TabularCap() = default;

  ModelState_TabularCap& operator=(const ModelState_TabularCap& state) = default;
  ModelState_TabularCap* operator=(const ModelState_TabularCap* state);

  friend std::ostream& operator<<(std::ostream& os,
                                  const ModelState_TabularCap& state)
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
       << ", X = " << state.capX << "\n"
       << "\t Points = " ;
    std::copy(state.yield_f_pts.begin(), state.yield_f_pts.end(),
              std::ostream_iterator<Uintah::Point>(os, " "));
    os << "\n";
    return os;
  }
};

} // End namespace Uintah

#endif // __MODEL_STATE_TABULAR_CAP_H__
