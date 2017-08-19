/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef __DERIVED_MODEL_STATE_CAMCLAY_DATA_H__
#define __DERIVED_MODEL_STATE_CAMCLAY_DATA_H__

#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Default.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_CamClay
  \brief A structure that stores the state data that is specialized for
         the CamClay model.
         ** Derived from PlasticityState:ModelState
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_CamClay : public ModelState_Default
{

public:
  double p_c;  // consolidation pressure
  double p_c0; // consolidation pressure at the beginning of time step

  double p; // pressure = tr(sigma)/3
  double q; // shear = sqrt(3J2); J2 = 1/2 s:s; s = sigma - p I

  double epse_v;    // volumetric elastic strain = tr(epse)
  double epse_s;    // deviatoric elastic strain = sqrt(2/3) ||ee||
                    //  ee = epse - 1/3 epse_v I
  double epse_v_tr; // trial volumetric elastic strain
  double epse_s_tr; // trial deviatoric elastic strain

  Uintah::Matrix3 elasticStrainTensor;
  Uintah::Matrix3 elasticStrainTensorTrial;

  ModelState_CamClay();

  ModelState_CamClay(const ModelState_CamClay& state);
  ModelState_CamClay(const ModelState_CamClay* state);

  ~ModelState_CamClay() override;

  ModelState_CamClay& operator=(const ModelState_CamClay& state);
  ModelState_CamClay* operator=(const ModelState_CamClay* state);
};

} // End namespace Uintah

#endif // __DERIVED_MODEL_STATE_CAMCLAY_DATA_H__
