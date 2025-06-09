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

#ifndef __MODEL_STATE_BASE_H__
#define __MODEL_STATE_BASE_H__

#include <Core/Math/Matrix3.h>

#include <memory>
#include <utility>

namespace Vaango {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ModelStateBase
  \brief Base class for structure that stores the model state
  \author Biswajit Banerjee \n
*/
/////////////////////////////////////////////////////////////////////////////

class ModelStateBase
{
public:
  double eqStrainRate;
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
  double yieldStress;
  double lambdaIncPlastic; // magnitude of plastic strain rate * delT
  double eqPlasticStrainRate;
  double eqPlasticStrain;
  double porosity;
  double energy;
  double I1;  // tr(stress)
  double J2;  // 1/2 dev(stress):dev(stress)
  double p;   // I1/3
  double q;   // sqrt(3 J2)

  Uintah::Matrix3 devStress;
  Uintah::Matrix3 backStress;
  Uintah::Matrix3 strainRate;
  Uintah::Matrix3 plasticFlowDirection;

private:

  Uintah::Matrix3 d_plasticStrain;
  Uintah::Matrix3 d_stress;

public:

  ModelStateBase();

  // Virtual copy constructor
  // We'll rely on the clone method for polymorphic copying.
  ModelStateBase(const ModelStateBase& state) = default;

  // Virtual move constructor
  ModelStateBase(ModelStateBase&& other) noexcept = default;

  virtual ~ModelStateBase() = default;

  // Virtual copy assignment operator (standard for polymorphic assignment)
  virtual ModelStateBase&
  operator=(const ModelStateBase& state);

  // Virtual move assignment operator
  virtual ModelStateBase&
  operator=(ModelStateBase&& state) noexcept;

  // Polymorphic cloning method (preferred over operator=(const T*))
  [[nodiscard]] virtual std::unique_ptr<ModelStateBase>
  clone() const
  {
    return std::make_unique<ModelStateBase>(*this);
  }

  virtual
  size_t numStateVar() const
  {
    auto numThis = 25u;
    return numThis;
  }

  /* Set the stress and get deviatoric stress */
  void setStress(const Uintah::Matrix3& stress);

  /* Get the stress */
  Uintah::Matrix3 getStress() const;

  /* Set the plastic strain */
  void setPlasticStrain(const Uintah::Matrix3& ep);

  /* Get the plastic strain */
  Uintah::Matrix3 getPlasticStrain() const;

  virtual 
  void updateStressInvariants() {}

private:

  /* Returns deviatoric stress */
  Uintah::Matrix3 updateStressInvariants(const Uintah::Matrix3& stress);

};

} // End namespace Uintah

#endif // __MODEL_STATE_BASE_H__
