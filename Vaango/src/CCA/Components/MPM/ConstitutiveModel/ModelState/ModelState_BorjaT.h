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

#ifndef __DERIVED_MODEL_STATE_CAMCLAY_BORJA_H__
#define __DERIVED_MODEL_STATE_CAMCLAY_BORJA_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateT.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_BorjaT
  \brief A structure that stores the state data that is specialized for
         the Borja model.
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_BorjaT : public ModelStateT<ModelState_BorjaT>
{

public:
  double p_c;  // consolidation pressure
  double p_c0; // consolidation pressure at the beginning of time step

  double epse_v;    // volumetric elastic strain = tr(epse)
  double epse_s;    // deviatoric elastic strain = sqrt(2/3) ||ee||
                    //  ee = epse - 1/3 epse_v I
  double epse_v_tr; // trial volumetric elastic strain
  double epse_s_tr; // trial deviatoric elastic strain

  Uintah::Matrix3 elasticStrainTensor;
  Uintah::Matrix3 elasticStrainTensorTrial;

  ModelState_BorjaT();

  ModelState_BorjaT(const ModelState_BorjaT& state);
  ModelState_BorjaT(const ModelState_BorjaT* state);

  ~ModelState_BorjaT();

  ModelState_BorjaT& operator=(const ModelState_BorjaT& state);
  ModelState_BorjaT* operator=(const ModelState_BorjaT* state);
  void copyLocalState(const ModelState_BorjaT* state);

  size_t numLocalStateVar() const 
  {
    auto numThis = 8u;
    return numThis;
  }

  void updateLocalStressInvariants(const Uintah::Matrix3& stress) {}
  void updateLocalStressInvariants() {}

  std::pair<Uintah::Matrix3, Uintah::Matrix3> 
  updateLocalStrainScalars(const Uintah::Matrix3& strain,
                           const Uintah::Matrix3& strain_trial) 
  {
    elasticStrainTensor = strain;
    elasticStrainTensorTrial = strain_trial;
    epse_v = strain.Trace(); 
    epse_v_tr = strain_trial.Trace(); 
    auto eps_dev = strain - Vaango::Util::Identity * (epse_v / 3.0);
    auto eps_tr_dev = strain_trial - Vaango::Util::Identity * (epse_v_tr / 3.0);
    auto ee = eps_dev.Contract(eps_dev);
    auto ee_tr = eps_tr_dev.Contract(eps_tr_dev);
    epse_s = std::sqrt(2.0 * ee / 3.0);
    epse_s_tr = std::sqrt(2.0 * ee_tr / 3.0);
    
    return std::make_pair(eps_dev, eps_tr_dev);
  }
};

} // End namespace Uintah

#endif // __DERIVED_MODEL_STATE_CAMCLAY_BORJA_H__
