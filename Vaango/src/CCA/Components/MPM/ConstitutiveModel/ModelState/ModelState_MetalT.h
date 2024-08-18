/*
 * The MIT License
 *
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __DERIVED_MODEL_STATE_METAL_TEMPLATED_H__
#define __DERIVED_MODEL_STATE_METAL_TEMPLATED_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateT.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_MetalT
  \brief A structure that stores the state data that is specialized for
         the Borja model.
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_MetalT : public ModelStateT<ModelState_MetalT>
{

public:

  ModelState_MetalT();

  ModelState_MetalT(const ModelState_MetalT& state);
  ModelState_MetalT(const ModelState_MetalT* state);

  ~ModelState_MetalT();

  ModelState_MetalT& operator=(const ModelState_MetalT& state);
  ModelState_MetalT* operator=(const ModelState_MetalT* state);
  void copyLocalState(const ModelState_MetalT* state);

  size_t numLocalStateVar() const 
  {
    auto numThis = 0u;
    return numThis;
  }

  void updateLocalStressInvariants(const Uintah::Matrix3& stress) {}
  void updateLocalStressInvariants() {}

  std::pair<Uintah::Matrix3, Uintah::Matrix3> 
  updateLocalStrainScalars(const Uintah::Matrix3& strain,
                           const Uintah::Matrix3& strain_trial) 
  {
    double epse_v = strain.Trace(); 
    double epse_v_tr = strain_trial.Trace(); 
    auto eps_dev = strain - Vaango::Util::Identity * (epse_v / 3.0);
    auto eps_tr_dev = strain_trial - Vaango::Util::Identity * (epse_v_tr / 3.0);
    return std::make_pair(eps_dev, eps_tr_dev);
  }
};

} // End namespace Uintah

#endif // __DERIVED_MODEL_STATE_METAL_TEMPLATED_H__
