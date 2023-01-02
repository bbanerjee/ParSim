/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __DERIVED_MODEL_STATE_SOIL_MODEL_BRANNON_DATA_H__
#define __DERIVED_MODEL_STATE_SOIL_MODEL_BRANNON_DATA_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState_SoilModelBrannon
  \brief A structure that stores the state data that is specialized for
         the SoilModelBrannon model.
         ** Derived from PlasticityState:ModelState
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelState_SoilModelBrannon : public ModelStateBase
{

public:
  Uintah::Matrix3 plasticStrainTensor; // The tensor form of plastic strain
  double kappa;                        // The cap kappa parameter
  double CR;                           // Initial cap radius
  double maxX;        // The max. cap hydrostatic compressive strength X
  double eps_v;       // Volumetrc strain
  double delta_eps_v; // Change in Volumetrc strain
  double scale_eps_v; // Scale factor for Volumetrc strain

  ModelState_SoilModelBrannon();

  ModelState_SoilModelBrannon(const ModelState_SoilModelBrannon& state);
  ModelState_SoilModelBrannon(const ModelState_SoilModelBrannon* state);

  ~ModelState_SoilModelBrannon() override;

  ModelState_SoilModelBrannon& operator=(
    const ModelState_SoilModelBrannon& state);
  ModelState_SoilModelBrannon* operator=(
    const ModelState_SoilModelBrannon* state);

  virtual 
  size_t numStateVar() const override
  {
    auto numBase = ModelStateBase::numStateVar();
    auto numThis = 7u;
    return numBase + numThis;
  }
};

} // End namespace Uintah

#endif // __DERIVED_MODEL_STATE_SOIL_MODEL_BRANNON_DATA_H__
