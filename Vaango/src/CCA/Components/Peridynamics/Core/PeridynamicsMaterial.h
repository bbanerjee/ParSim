/*
 * The MIT License
 *
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

#ifndef __VAANGO_PERIDYNAMICS_MATERIAL_H__
#define __VAANGO_PERIDYNAMICS_MATERIAL_H__

#include <Core/Grid/Material.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <vector>

namespace Uintah {
class VarLabel;
class Patch;
class DataWarehouse;
class GeometryObject;
} // namespace Uintah

namespace Vaango {

class PeridynamicsFlags;
class PeridynamicsLabel;
class PeridynamicsDamageModel;
class PeridynamicsMaterialModel;
class ParticleCreator;

class PeridynamicsMaterial : public Uintah::Material
{

public:
  // Default Constructor
  PeridynamicsMaterial();

  // Standard Peridynamics Material Constructor
  PeridynamicsMaterial(Uintah::ProblemSpecP& ps,
                       Uintah::MaterialManagerP& mat_manager,
                       PeridynamicsFlags* flags,
                       bool isRestart);

  ~PeridynamicsMaterial() override;

  // Prevent copying/move of this class
  PeridynamicsMaterial(const PeridynamicsMaterial& mpmm) = delete;
  PeridynamicsMaterial(PeridynamicsMaterial&& mpmm)      = delete;
  PeridynamicsMaterial&
  operator=(const PeridynamicsMaterial& mpmm) = delete;
  PeridynamicsMaterial&
  operator=(PeridynamicsMaterial&& mpmm) = delete;

  virtual Uintah::ProblemSpecP
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  using VarLabelVector = std::vector<const Uintah::VarLabel*>;
  virtual void
  registerParticleState(std::vector<VarLabelVector>& state,
                        std::vector<VarLabelVector>& state_preReloc);

  // Return correct constitutive model pointer for this material
  PeridynamicsMaterialModel*
  getMaterialModel() const;

  // Return correct basic damage model pointer for this material
  PeridynamicsDamageModel*
  getDamageModel() const;

  double
  getInitialDensity() const
  {
    return d_density;
  }

  Uintah::particleIndex
  countParticles(const Uintah::Patch* patch);
  void
  createParticles(Uintah::particleIndex numParticles,
                  Uintah::CCVariable<short int>& cellNAPID,
                  const Uintah::Patch* patch,
                  Uintah::DataWarehouse* new_dw);

  ParticleCreator*
  getParticleCreator();

private:
  PeridynamicsLabel* d_varLabel{nullptr};
  PeridynamicsMaterialModel* d_materialModel{nullptr};
  PeridynamicsDamageModel* d_damageModel{nullptr};
  ParticleCreator* d_particle_creator{nullptr};

  double d_density{0.0};

  std::vector<Uintah::GeometryObject*> d_geom_objs;

  ///////////////////////////////////////////////////////////////////////////
  //
  // The standard set of initialization actions except particlecreator
  //
  void
  standardInitialization(Uintah::ProblemSpecP& ps,
                         const Uintah::MaterialManagerP& mat_manager,
                         PeridynamicsFlags* flags,
                         bool isRestart);
};

} // End namespace Vaango

#endif // __VAANGO_PERIDYNAMICS_MATERIAL_H__
