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

#ifndef __CCA_COMPONENTS_MPM_MPMMATERIAL_H__
#define __CCA_COMPONENTS_MPM_MPMMATERIAL_H__

// Do not EVER put a #include for anything in another CCA/Components in here.
// (#includes of other MPM files is ok.  However, if you #include'd ARCHES
// or something, then a circular dependency would be created.)

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

// This is to avoid circular dependencies between MPMMaterial and
// BasicDamageModel
// Better design needed? --BB
namespace Vaango {
class BasicDamageModel;
}

namespace Uintah {

class Patch;
class Grid;
class DataWarehouse;
class VarLabel;
class GeometryObject;
class ConstitutiveModel;
class MPMLabel;
class ParticleCreator;
class ScalarDiffusionModel;

class MPMMaterial final : public Material
{
public:
  // Default Constructor
  MPMMaterial();

  // Standard MPM Material Constructor
  explicit MPMMaterial(ProblemSpecP& ps,
                       MaterialManagerP& mats,
                       MPMFlags* flags,
                       bool isRestart);

  ~MPMMaterial() override;

  // Prevent copying of this class
  MPMMaterial(const MPMMaterial& mpmm) = delete;
  auto
  operator=(const MPMMaterial& mpmm) -> MPMMaterial& = delete;

  // void registerParticleState(MaterialManager* ss) override;
  using VarLabelVector = std::vector<const VarLabel*>;
  void
  registerParticleState(std::vector<VarLabelVector>& state,
                        std::vector<VarLabelVector>& state_preReloc);

  auto
  outputProblemSpec(ProblemSpecP& ps) -> ProblemSpecP override;

  /*!  Create a copy of the material without the associated geometry */
  void
  copyWithoutGeom(ProblemSpecP& ps, const MPMMaterial* mat, MPMFlags* flags);

  void
  deleteGeomObjects()
  {
    d_geom_objs.clear();
  }

  [[nodiscard]] auto
  nullGeomObject() const -> int;

  [[nodiscard]] auto
  getConstitutiveModel() const -> ConstitutiveModel*;

  [[nodiscard]] auto
  getBasicDamageModel() const -> Vaango::BasicDamageModel*;

  [[nodiscard]] auto
  getScalarDiffusionModel() const -> ScalarDiffusionModel*;

  auto
  createParticles(CCVariable<short int>& cellNAPID,
                  const Patch*,
                  DataWarehouse* new_dw) -> particleIndex;

  auto
  getParticleCreator() -> ParticleCreator*;

  [[nodiscard]] auto
  getInitialDensity() const -> double;

  // Get the specific heats at room temperature
  [[nodiscard]] auto
  getInitialCp() const -> double;
  [[nodiscard]] auto
  getInitialCv() const -> double;

  // for temperature dependent plasticity models
  [[nodiscard]] auto
  getRoomTemperature() const -> double;
  [[nodiscard]] auto
  getMeltTemperature() const -> double;

  [[nodiscard]] auto
  getIsRigid() const -> bool
  {
    return d_is_rigid;
  }
  void
  setIsRigid(bool is_rigid)
  {
    d_is_rigid = is_rigid;
  }

  [[nodiscard]] auto
  doBasicDamage() const -> bool
  {
    return d_doBasicDamage;
  }

  [[nodiscard]] auto
  getIncludeFlowWork() const -> bool
  {
    return d_includeFlowWork;
  }

  [[nodiscard]] auto
  getSpecificHeat() const -> double
  {
    return d_specificHeat;
  }
  [[nodiscard]] auto
  getThermalConductivity() const -> double
  {
    return d_thermalConductivity;
  }

  // For scalar diffusion
  auto
  doConcReduction() -> bool
  {
    return d_doConcReduction;
  };

  // For activating particle insertion
  [[nodiscard]] auto
  getActivationTime() const -> double
  {
    return d_activation_time;
  }

  [[nodiscard]] auto
  isActive() const -> bool
  {
    return d_is_active;
  }

  void
  setActive(bool flag)
  {
    d_is_active = flag;
  }

  // Rigid material that transmitss force
  [[nodiscard]] auto
  isRigidForceTransmittingMaterial() const -> bool
  {
    return d_is_force_transmitting_material;
  }

  void
  setAsRigidForceTransmittingMaterial(bool flag)
  {
    d_is_force_transmitting_material = flag;
  }

  // For hydromechanical coupling
  [[nodiscard]] auto
  getWaterDensity() const -> double
  {
    return d_waterdensity;
  }

  [[nodiscard]] auto
  getPorosity() const -> double
  {
    return d_porosity;
  }

  [[nodiscard]] auto
  getPermeability() const -> double
  {
    return d_permeability;
  }

  [[nodiscard]] auto
  getInitialPorePressure() const -> double
  {
    return d_initial_porepressure;
  }

  // For MPMICE
  [[nodiscard]] auto
  getGamma() const -> double;

  void
  initializeCCVariables(CCVariable<double>& rhom,
                        CCVariable<double>& rhC,
                        CCVariable<double>& temp,
                        CCVariable<Vector>& vCC,
                        CCVariable<double>& vfCC,
                        const Patch* patch);

  void
  initializeDummyCCVariables(CCVariable<double>& rhom,
                             CCVariable<double>& rhC,
                             CCVariable<double>& temp,
                             CCVariable<Vector>& vCC,
                             CCVariable<double>& vfCC,
                             const Patch* patch);

private:
  // The standard set of initialization actions except particlecreator
  void
  standardInitialization(ProblemSpecP& ps,
                         MaterialManagerP& mats,
                         MPMFlags* flags,
                         const bool isRestart);

private:
  const MPMFlags* d_flags{ nullptr };
  std::unique_ptr<MPMLabel> d_lb{ nullptr };
  std::unique_ptr<ConstitutiveModel> d_cm{ nullptr };
  std::unique_ptr<ParticleCreator> d_particle_creator{ nullptr };
  std::unique_ptr<ScalarDiffusionModel> d_sdm{ nullptr };
  std::unique_ptr<Vaango::BasicDamageModel> d_basicDamageModel{ nullptr };
  std::vector<std::shared_ptr<GeometryObject>> d_geom_objs;

  bool d_doBasicDamage{ false };
  bool d_is_rigid{ false }; // for implicit rigid body contact
  bool d_includeFlowWork{ false };

  double d_density{ 0.0 };
  double d_specificHeat{ 0.0 };
  double d_thermalConductivity{ 0.0 };

  // Specific heats at constant pressure and constant volume
  // (values at room temperature - [273.15 + 20] K)
  double d_Cp{ 0.0 }, d_Cv{ 0.0 };

  // for temperature dependent plasticity models
  double d_troom{ 0.0 };
  double d_tmelt{ 0.0 };

  // For scalar diffusion
  bool d_doConcReduction{ false };

  // For rigid-body force transmission
  bool d_is_force_transmitting_material{ false };

  // For activating inserted particles
  bool d_is_active{ true };
  double d_activation_time{ 0.0 };

  // For hydromechanical coupling
  double d_waterdensity{ 0.0 }, d_porosity{ 0.0 }, d_permeability{ 0.0 },
    d_initial_porepressure{ 0.0 };

  // For triangulated surfaces
  bool d_all_triangle_geometry{ true };
};

} // End namespace Uintah

#endif // __CCA_COMPONENTS_MPM_MPMMATERIAL_H__
