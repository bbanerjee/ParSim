/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __MODEL_STATE_ARENISCA3_H__
#define __MODEL_STATE_ARENISCA3_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_Arenisca3
  \brief A structure that stores the state data that is specialized for
         the Arenisca3 model.
         ** Derived from PlasticityState:ModelState
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_Arenisca3 : public ModelStateBase
{

public:
  double I1;                           // I1 = Tr(sigma)
  double sqrt_J2;                      // sqrt(J2)
  Uintah::Matrix3 plasticStrainTensor; // The tensor form of plastic strain
  double kappa;                        // The cap kappa parameter
  double capX; // The cap hydrostatic compressive strength X
  double zeta; // The back stress parameter

  ModelState_Arenisca3();

  ModelState_Arenisca3(const ModelState_Arenisca3& state);
  ModelState_Arenisca3(const ModelState_Arenisca3* state);

  ~ModelState_Arenisca3() override;

  ModelState_Arenisca3& operator=(const ModelState_Arenisca3& state);
  ModelState_Arenisca3* operator=(const ModelState_Arenisca3* state);

  virtual 
  size_t numStateVar() const override
  {
    auto numBase = ModelStateBase::numStateVar();
    auto numThis = 6u;
    return numBase + numThis;
  }
};

} // End namespace Uintah

#endif // __MODEL_STATE_ARENISCA3_H__
