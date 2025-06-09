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
#include "ModelStateBase.h"

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

ModelStateBase&
ModelStateBase::operator=(const ModelStateBase& state)
{
  if (this == &state)
    return *this;

  this->eqStrainRate = state.eqStrainRate;
  this->pressure = state.pressure;
  this->temperature = state.temperature;
  this->initialTemperature = state.initialTemperature;
  this->density = state.density;
  this->initialDensity = state.initialDensity;
  this->volume = state.volume;
  this->initialVolume = state.initialVolume;
  this->bulkModulus = state.bulkModulus;
  this->initialBulkModulus = state.initialBulkModulus;
  this->shearModulus = state.shearModulus;
  this->initialShearModulus = state.initialShearModulus;
  this->meltingTemp = state.meltingTemp;
  this->initialMeltTemp = state.initialMeltTemp;
  this->specificHeat = state.specificHeat;
  this->yieldStress = state.yieldStress;
  this->lambdaIncPlastic = state.lambdaIncPlastic;
  this->eqPlasticStrainRate = state.eqPlasticStrainRate;
  this->eqPlasticStrain = state.eqPlasticStrain;
  this->porosity = state.porosity;
  this->energy = state.energy;
  this->I1 = state.I1;
  this->J2 = state.J2;
  this->p = state.p;
  this->q = state.q;
  this->devStress = state.devStress;
  this->backStress = state.backStress;
  this->strainRate = state.strainRate;
  this->plasticFlowDirection = state.plasticFlowDirection;
  this->d_plasticStrain = state.d_plasticStrain;
  this->d_stress = state.d_stress;
  return *this;
}

ModelStateBase&
Vaango::ModelStateBase::operator=(ModelStateBase&& state) noexcept
{
  if (this == &state) { // Correct self-assignment check for reference
    return *this;
  }
  this->eqStrainRate = std::move(state.eqStrainRate);
  this->pressure = std::move(state.pressure);
  this->temperature = std::move(state.temperature);
  this->initialTemperature = std::move(state.initialTemperature);
  this->density = std::move(state.density);
  this->initialDensity = std::move(state.initialDensity);
  this->volume = std::move(state.volume);
  this->initialVolume = std::move(state.initialVolume);
  this->bulkModulus = std::move(state.bulkModulus);
  this->initialBulkModulus = std::move(state.initialBulkModulus);
  this->shearModulus = std::move(state.shearModulus);
  this->initialShearModulus = std::move(state.initialShearModulus);
  this->meltingTemp = std::move(state.meltingTemp);
  this->initialMeltTemp = std::move(state.initialMeltTemp);
  this->specificHeat = std::move(state.specificHeat);
  this->yieldStress = std::move(state.yieldStress);
  this->lambdaIncPlastic = std::move(state.lambdaIncPlastic);
  this->eqPlasticStrainRate = std::move(state.eqPlasticStrainRate);
  this->eqPlasticStrain = std::move(state.eqPlasticStrain);
  this->porosity = std::move(state.porosity);
  this->energy = std::move(state.energy);
  this->I1 = std::move(state.I1);
  this->J2 = std::move(state.J2);
  this->p = std::move(state.p);
  this->q = std::move(state.q);
  this->devStress = std::move(state.devStress);
  this->backStress = std::move(state.backStress);
  this->strainRate = std::move(state.strainRate);
  this->plasticFlowDirection = std::move(state.plasticFlowDirection);
  this->d_plasticStrain = std::move(state.d_plasticStrain);
  this->d_stress = std::move(state.d_stress);
  return *this;
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
