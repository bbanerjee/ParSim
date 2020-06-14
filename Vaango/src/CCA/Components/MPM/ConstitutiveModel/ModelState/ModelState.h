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

#ifndef __TEMPLATE_MODEL_STATE_H__
#define __TEMPLATE_MODEL_STATE_H__

#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Math/Matrix3.h>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelState
  \brief Base class for structure that stores the model state
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

template <typename DerivedT>
class ModelState
{
public:

  double yieldStress;
  double strainRate;
  double plasticStrainRate;
  double plasticStrain;
  double pressure;
  double temperature;
  double initialTemperature;
  double density;
  double initialDensity;
  double volume;
  double initialVolume;
  double bulkModulus;
  double initialBulkModulus;
  double shearModulus;
  double initialShearModulus;
  double meltingTemp;
  double initialMeltTemp;
  double specificHeat;
  double porosity;
  double energy;
  double I1;  // tr(stress)
  double J2;  // 1/2 dev(stress):dev(stress)
  double p;   // I1/3
  double q;   // sqrt(3 J2)
  const Uintah::Matrix3* backStress;

  void copyState(const DerivedT* state)
  {
    yieldStress = state->yieldStress;
    strainRate = state->strainRate;
    plasticStrainRate = state->plasticStrainRate;
    plasticStrain = state->plasticStrain;
    pressure = state->pressure;
    temperature = state->temperature;
    initialTemperature = state->initialTemperature;
    density = state->density;
    initialDensity = state->initialDensity;
    volume = state->volume;
    initialVolume = state->initialVolume;
    bulkModulus = state->bulkModulus;
    initialBulkModulus = state->initialBulkModulus;
    shearModulus = state->shearModulus;
    initialShearModulus = state->initialShearModulus;
    meltingTemp = state->meltingTemp;
    initialMeltTemp = state->initialMeltTemp;
    specificHeat = state->specificHeat;
    porosity = state->porosity;
    energy = state->energy;
    I1 = state->I1;
    J2 = state->J2;
    p = state->p;
    q = state->q;
    backStress = state->backStress;

    derived()->copyLocalState(state);
  }

  size_t numStateVar() const
  {
    return derived()->numLocalStateVar() + 25u;
  }

  /* Returns deviatoric stress */
  Uintah::Matrix3 updateStressInvariants(const Uintah::Matrix3& stress)
  {
    I1 = stress.Trace(); 
    p = I1 / 3.0;
    pressure = p;
    Uintah::Matrix3 devStress = stress - Vaango::Util::Identity * p; 
    J2 = 0.5 * devStress.Contract(devStress);
    q = std::sqrt(3.0 * J2);

    derived()->updateLocalStressInvariants(stress);

    return devStress;
  }

  void updateStressInvariants() 
  {
    derived()->updateLocalStressInvariants();
  }

  /* Returns deviatoric strain */
  std::pair<Uintah::Matrix3, Uintah::Matrix3> 
  updateStrainScalars(const Uintah::Matrix3& strain,
                      const Uintah::Matrix3& strain_trial)
  {
    return derived()->updateLocalStrainScalars(strain, strain_trial);
  }

private:

  ModelState()
  {
    yieldStress = 0.0;
    strainRate = 0.0;
    plasticStrainRate = 0.0;
    plasticStrain = 0.0;
    pressure = 0.0;
    temperature = 0.0;
    initialTemperature = 0.0;
    density = 0.0;
    initialDensity = 0.0;
    volume = 0.0;
    initialVolume = 0.0;
    bulkModulus = 0.0;
    initialBulkModulus = 0.0;
    shearModulus = 0.0;
    initialShearModulus = 0.0;
    meltingTemp = 0.0;
    initialMeltTemp = 0.0;
    specificHeat = 0.0;
    porosity = 0.0;
    energy = 0.0;
    I1 = 0.0;
    J2 = 0.0;
    p = 0.0;
    q = 0.0;
    backStress = nullptr;
  }

  DerivedT* derived()
  {
    return static_cast<DerivedT*>(this);
  }

  const DerivedT* derived() const
  {
    return static_cast<const DerivedT*>(this);
  }

  friend DerivedT;

};

} // End namespace Vaango

#endif // __TEMPLATE_MODEL_STATE_H__
