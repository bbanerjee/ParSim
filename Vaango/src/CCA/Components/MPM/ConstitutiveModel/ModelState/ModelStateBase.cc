/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>

using namespace Vaango;

ModelStateBase::ModelStateBase()
{
  eqStrainRate = 0.0;
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
  yieldStress = 0.0;
  lambdaIncPlastic = 0.0;
  eqPlasticStrainRate = 0.0;
  eqPlasticStrain = 0.0;
  porosity = 0.0;
  energy = 0.0;
  I1 = 0.0;
  J2 = 0.0;
  p = 0.0;
  q = 0.0;
  devStress = Vaango::Util::Zero;
  backStress = Vaango::Util::Zero;
  strainRate = Vaango::Util::Zero;
  plasticFlowDirection = Vaango::Util::Zero;
  d_plasticStrain = Vaango::Util::Zero;
  d_stress = Vaango::Util::Zero;
}

ModelStateBase::ModelStateBase(const ModelStateBase& state)
{
  *this = &state;
}

ModelStateBase::ModelStateBase(const ModelStateBase* state)
{
  *this = state;
}

ModelStateBase::~ModelStateBase() = default;

ModelStateBase&
ModelStateBase::operator=(const ModelStateBase& state_in)
{
  const ModelStateBase* state = &state_in;
  if (this == state)
    return *this;

  *this = state;
  return *this;
}

ModelStateBase*
ModelStateBase::operator=(const ModelStateBase* state)
{
  if (this == state)
    return this;
  eqStrainRate = state->eqStrainRate;
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
  yieldStress = state->yieldStress;
  lambdaIncPlastic = state->lambdaIncPlastic;
  eqPlasticStrainRate = state->eqPlasticStrainRate;
  eqPlasticStrain = state->eqPlasticStrain;
  porosity = state->porosity;
  energy = state->energy;
  I1 = state->I1;
  J2 = state->J2;
  p = state->p;
  q = state->q;
  devStress = state->devStress;
  backStress = state->backStress;
  strainRate = state->strainRate;
  plasticFlowDirection = state->plasticFlowDirection;
  d_plasticStrain = state->d_plasticStrain;
  d_stress = state->d_stress;
  return this;
}

/* Set the stress and get deviatoric stress */
void
ModelStateBase::setStress(const Uintah::Matrix3& stress) 
{
  d_stress = stress;
  devStress = updateStressInvariants(stress); 
}

/* Get the stress */
Uintah::Matrix3 
ModelStateBase::getStress() const 
{
  return d_stress;
}

/* Set the plastic strain */
void 
ModelStateBase::setPlasticStrain(const Uintah::Matrix3& ep) 
{
  d_plasticStrain = ep;
  //eqPlasticStrain = ep.Norm();
}

/* Get the plastic strain */
Uintah::Matrix3 
ModelStateBase::getPlasticStrain() const 
{
  return d_plasticStrain;
}

Uintah::Matrix3
ModelStateBase::updateStressInvariants(const Uintah::Matrix3& stress)
{
  I1 = stress.Trace(); 
  p = I1 / 3.0;
  pressure = p;

  Uintah::Matrix3 s = stress - Vaango::Util::Identity * p; 
  J2 = 0.5 * s.Contract(s);
  q = std::sqrt(3.0 * J2);
  return s;
}
