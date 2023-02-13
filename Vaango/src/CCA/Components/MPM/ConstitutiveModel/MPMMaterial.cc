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

//  MPMMaterial.cc

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreatorFactory.h>
#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ScalarDiffusionModel.h>
#include <CCA/Components/MPM/ReactionDiffusion/ScalarDiffusionModelFactory.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/IntVector.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/NullGeometryPiece.h>
#include <Core/GeometryPiece/TriGeometryPiece.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <iostream>
#include <list>
#include <string>

#define d_TINY_RHO 1.0e-12 // also defined  ICE.cc and ICEMaterial.cc

namespace Uintah {

// Default constructor
MPMMaterial::MPMMaterial()
{
  d_lb = std::make_unique<MPMLabel>();
}

// Standard Constructor
MPMMaterial::MPMMaterial(ProblemSpecP& ps,
                         MaterialManagerP& matManager,
                         MPMFlags* flags,
                         bool isRestart)
  : Material(ps)
  , d_flags(flags)
{
  d_lb = std::make_unique<MPMLabel>();
  standardInitialization(ps, matManager, flags, isRestart);
  d_particle_creator = ParticleCreatorFactory::create(ps, this, flags);
}

void
MPMMaterial::standardInitialization(ProblemSpecP& ps,
                                    MaterialManagerP& matManager,
                                    MPMFlags* flags,
                                    bool isRestart)

{
  // Follow the layout of the input file
  // Steps:
  // 1.  Determine the type of constitutive model and create it.
  // 1.1  Determine if scalar diffusion is used and the type of
  //      scalar diffusion model and create it.
  //      Added for reactive flow component.
  // 2.  Get the general properties of the material such as
  //     density, thermal_conductivity, specific_heat.
  // 3.  Loop through all of the geometry pieces that make up a single
  //     geometry object.
  // 4.  Within the geometry object, assign the boundary conditions
  //     to the object.
  // 5.  Assign the velocity field.

  // Step 0 -- create the basic damage model
  d_doBasicDamage = false;
  ps->get("do_basic_damage", d_doBasicDamage);
  if (d_doBasicDamage) {
    d_basicDamageModel = std::make_unique<Vaango::BasicDamageModel>(ps, flags);
    d_basicDamageModel->setMaterialManager(matManager);
  }

  // Step 1 -- create the constitutive gmodel.
  d_cm = ConstitutiveModelFactory::create(ps, flags);
  d_cm->setMaterialManager(matManager);
  if (!d_cm) {
    std::ostringstream desc;
    desc << "An error occured in the ConstitutiveModelFactory that has \n"
         << "slipped through the existing bullet proofing.\n";
    throw ParameterNotFound(desc.str(), __FILE__, __LINE__);
  }

  // Step 1.1 -- check if scalar diffusion is used and
  // create the scalar diffusion model.
  if (flags->d_doScalarDiffusion) {
    d_sdm = ScalarDiffusionModelFactory::create(ps, matManager, flags);
  }

  // Step 2 -- get the general material properties
  ps->require("density", d_density);
  ps->require("thermal_conductivity", d_thermalConductivity);
  ps->require("specific_heat", d_specificHeat);

  // Assume the the centered specific heat is C_v
  d_Cv = d_specificHeat;

  // Set C_p = C_v if not C_p data are entered
  d_Cp = d_Cv;
  ps->get("C_p", d_Cp);

  d_troom = 294.0;
  d_tmelt = 1.0e6;
  ps->get("room_temp", d_troom);
  ps->get("melt_temp", d_tmelt);

  // This is currently only used in the implicit code, but should
  // be put to use in the explicit code as well.
  d_is_rigid = false;
  ps->get("is_rigid", d_is_rigid);

  d_includeFlowWork = false;
  ps->get("includeFlowWork", d_includeFlowWork);

  // For scalar diffusion
  // This is used for the autocycleflux boundary conditions
  d_doConcReduction = false;
  ps->get("do_conc_reduction", d_doConcReduction);

  // For activating materials if needed
  d_is_active = true;
  ps->get("is_active", d_is_active);
  d_activation_time = 0.0;
  ps->get("activation_time", d_activation_time);

  // Material is force transmitting (moves according to sum of forces)
  d_is_force_transmitting_material = false;
  ps->get("is_force_transmitting_material", d_is_force_transmitting_material);
  flags->d_reductionVars->sumTransmittedForce = true;

  // Also use for Darcy momentum exchange model
  ps->get("permeability", d_permeability);

  // For MPM hydro-mechanical coupling
  if (flags->d_coupledFlow) {
    // Rigid material does not require porosity and permeability
    if (!ps->findBlockWithAttributeValue(
          "constitutive_model", "type", "rigid")) {
      ps->require("water_density", d_waterdensity);
      ps->require("porosity", d_porosity);
      d_initial_porepressure = 0.0;
      ps->get("initial_pore_pressure", d_initial_porepressure);
    }
  }

  // Step 3 -- Loop through all of the pieces in this geometry object
  // int piece_num = 0;
  std::list<GeometryObject::DataItem> geom_obj_data;
  geom_obj_data.emplace_back("res", GeometryObject::IntVector);
  geom_obj_data.emplace_back("temperature", GeometryObject::Double);
  geom_obj_data.emplace_back("velocity", GeometryObject::Vector);
  geom_obj_data.emplace_back("affineTransformation_A0", GeometryObject::Vector);
  geom_obj_data.emplace_back("affineTransformation_A1", GeometryObject::Vector);
  geom_obj_data.emplace_back("affineTransformation_A2", GeometryObject::Vector);
  geom_obj_data.emplace_back("affineTransformation_b", GeometryObject::Vector);
  geom_obj_data.emplace_back("volumeFraction", GeometryObject::Double);

  if (flags->d_withColor) {
    geom_obj_data.emplace_back("color", GeometryObject::Double);
  }

  // ReactiveFlow Diffusion Component
  if (flags->d_doScalarDiffusion) {
    geom_obj_data.emplace_back("concentration", GeometryObject::Double);
  }

  // If this ia a restart, don't create geometry
  if (isRestart) {
    return;
  }

  // Create the geometry pieces
  for (ProblemSpecP geom_obj_ps = ps->findBlock("geom_object");
       geom_obj_ps != 0;
       geom_obj_ps = geom_obj_ps->findNextBlock("geom_object")) {
    std::vector<GeometryPieceP> pieces;
    GeometryPieceFactory::create(geom_obj_ps, pieces);

    // Tag if a triangle geometry piece exists in this set of objects
    for (const auto& geom_piece : pieces) {
      if (!static_cast<TriGeometryPiece*>(geom_piece.get())) {
        d_all_triangle_geometry = false;
      }
    }

    GeometryPieceP mainpiece;
    if (pieces.size() == 0) {
      throw ParameterNotFound(
        "No piece specified in geom_object", __FILE__, __LINE__);
    } else if (pieces.size() > 1) {
      mainpiece = std::make_shared<UnionGeometryPiece>(pieces);
    } else {
      mainpiece = pieces[0];
    }

    d_geom_objs.push_back(
      std::make_shared<GeometryObject>(mainpiece, geom_obj_ps, geom_obj_data));
  }
}

MPMMaterial::~MPMMaterial() = default;

void
MPMMaterial::registerParticleState(std::vector<VarLabelVector>& state,
                                   std::vector<VarLabelVector>& state_preReloc)
{
  state.push_back(d_particle_creator->returnParticleState());
  state_preReloc.push_back(d_particle_creator->returnParticleStatePreReloc());
}

auto
MPMMaterial::outputProblemSpec(ProblemSpecP& ps) -> ProblemSpecP
{
  ProblemSpecP mpm_ps = Material::outputProblemSpec(ps);
  mpm_ps->appendElement("density", d_density);
  mpm_ps->appendElement("thermal_conductivity", d_thermalConductivity);
  mpm_ps->appendElement("specific_heat", d_specificHeat);
  mpm_ps->appendElement("C_p", d_Cp);
  mpm_ps->appendElement("room_temp", d_troom);
  mpm_ps->appendElement("melt_temp", d_tmelt);
  mpm_ps->appendElement("is_rigid", d_is_rigid);

  mpm_ps->appendElement("includeFlowWork", d_includeFlowWork);

  mpm_ps->appendElement("do_conc_reduction", d_doConcReduction);

  mpm_ps->appendElement("do_basic_damage", d_doBasicDamage);
  if (d_doBasicDamage) {
    d_basicDamageModel->outputProblemSpecDamage(mpm_ps);
  }
  mpm_ps->appendElement("is_force_transmitting_material",
                        d_is_force_transmitting_material);
  mpm_ps->appendElement("is_active", d_is_active);
  mpm_ps->appendElement("activation_time", d_activation_time);

  mpm_ps->appendElement("permeability", d_permeability);

  // For MPM hydro-mechanical coupling
  if (d_flags->d_coupledFlow) {
    mpm_ps->appendElement("water_density", d_waterdensity);
    mpm_ps->appendElement("porosity", d_porosity);
    mpm_ps->appendElement("initial_pore_pressure", d_initial_porepressure);
  }

  d_cm->outputProblemSpec(mpm_ps);

  if (getScalarDiffusionModel()) {
    d_sdm->outputProblemSpec(mpm_ps);
  }

  for (const auto& obj : d_geom_objs) {
    obj->outputProblemSpec(mpm_ps);
  }

  return mpm_ps;
}

void
MPMMaterial::copyWithoutGeom(ProblemSpecP& ps,
                             const MPMMaterial* mat,
                             MPMFlags* flags)
{
  d_doBasicDamage = mat->d_doBasicDamage;
  if (d_doBasicDamage) {
    d_basicDamageModel = mat->d_basicDamageModel->clone();
  }

  d_cm                  = mat->d_cm->clone();
  d_density             = mat->d_density;
  d_thermalConductivity = mat->d_thermalConductivity;
  d_specificHeat        = mat->d_specificHeat;
  d_Cv                  = mat->d_Cv;
  d_Cp                  = mat->d_Cp;
  d_troom               = mat->d_troom;
  d_tmelt               = mat->d_tmelt;
  d_is_rigid            = mat->d_is_rigid;

  d_includeFlowWork = mat->d_includeFlowWork;
  d_doConcReduction = mat->d_doConcReduction;

  d_is_force_transmitting_material = mat->d_is_force_transmitting_material;
  d_is_active                      = mat->d_is_active;
  d_activation_time                = mat->d_activation_time;

  d_permeability = mat->d_permeability;
  if (d_flags->d_coupledFlow) {
    d_waterdensity         = mat->d_waterdensity;
    d_porosity             = mat->d_porosity;
    d_initial_porepressure = mat->d_initial_porepressure;
  }

  // Check to see which ParticleCreator object we need
  d_particle_creator = ParticleCreatorFactory::create(ps, this, flags);
}

auto
MPMMaterial::getConstitutiveModel() const -> ConstitutiveModel*
{
  // Return the pointer to the constitutive model associated
  // with this material

  return d_cm.get();
}

auto
MPMMaterial::getBasicDamageModel() const -> Vaango::BasicDamageModel*
{
  // Return the pointer to the basic damage model associated
  // with this material

  return d_basicDamageModel.get();
}

auto
MPMMaterial::getScalarDiffusionModel() const -> ScalarDiffusionModel*
{
  return d_sdm.get();
}

auto
MPMMaterial::createParticles(CCVariable<short int>& cellNAPID,
                             const Patch* patch,
                             DataWarehouse* new_dw) -> particleIndex
{
  return d_particle_creator->createParticles(
    this, cellNAPID, patch, new_dw, d_geom_objs);
}

auto
MPMMaterial::getParticleCreator() -> ParticleCreator*
{
  return d_particle_creator.get();
}

auto
MPMMaterial::getInitialDensity() const -> double
{
  return d_density;
}

auto
MPMMaterial::getInitialCp() const -> double
{
  return d_Cp;
}

auto
MPMMaterial::getInitialCv() const -> double
{
  return d_Cv;
}

auto
MPMMaterial::getRoomTemperature() const -> double
{
  return d_troom;
}

auto
MPMMaterial::getMeltTemperature() const -> double
{
  return d_tmelt;
}

auto
MPMMaterial::nullGeomObject() const -> int
{
  int count = 0;
  for (const auto& object : d_geom_objs) {
    GeometryPieceP piece = object->getPiece();
    auto* null_piece     = dynamic_cast<NullGeometryPiece*>(piece.get());
    if (null_piece) {
      return count;
    }
    count++;
  }
  return -1;
}

/* ---------------------------------------------------------------------
   Function~  MPMMaterial::initializeCells--
   Notes:  This function initializeCCVariables.  Reasonable values for
   CC Variables need to be present in all the cells and evolve, even though
   there is no mass.  This is essentially the same routine that is in
   ICEMaterial.cc
   _____________________________________________________________________*/
void
MPMMaterial::initializeCCVariables(CCVariable<double>& rho_micro,
                                   CCVariable<double>& rho_CC,
                                   CCVariable<double>& temp,
                                   CCVariable<Vector>& vel_CC,
                                   CCVariable<double>& vol_frac_CC,
                                   const Patch* patch)
{
  // initialize to -9 so bullet proofing will catch it any cell that
  // isn't initialized
  vel_CC.initialize(Vector(0., 0., 0.));
  rho_micro.initialize(-9.0);
  rho_CC.initialize(-9.0);
  temp.initialize(-9.0);
  vol_frac_CC.initialize(0.0);
  Vector dx = patch->dCell();

  for (auto& d_geom_obj : d_geom_objs) {
    GeometryPieceP piece = d_geom_obj->getPiece();
    IntVector ppc        = d_geom_obj->getInitialData_IntVector("res");
    Vector dxpp          = patch->dCell() / ppc;
    Vector dcorner       = dxpp * 0.5;
    double totalppc      = ppc.x() * ppc.y() * ppc.z();

    // Find the bounds of a region a little bigger than the piece's BBox.
    Box bb = piece->getBoundingBox();

    Point bb_low(bb.lower().x() - 3.0 * dx.x(),
                 bb.lower().y() - 3.0 * dx.y(),
                 bb.lower().z() - 3.0 * dx.z());

    Point bb_up(bb.upper().x() + 3.0 * dx.x(),
                bb.upper().y() + 3.0 * dx.y(),
                bb.upper().z() + 3.0 * dx.z());

    for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
         iter++) {
      IntVector c = *iter;

      Point lower = patch->nodePosition(*iter) + dcorner;
      int count   = 0;

      for (int ix = 0; ix < ppc.x(); ix++) {
        for (int iy = 0; iy < ppc.y(); iy++) {
          for (int iz = 0; iz < ppc.z(); iz++) {
            IntVector idx(ix, iy, iz);
            Point p = lower + dxpp * idx;
            if (piece->inside(p)) {
              count++;
            }
          }
        }
      }

      double ups_volFrac = d_geom_obj->getInitialData_double("volumeFraction");
      if (ups_volFrac == -1.0) {
        vol_frac_CC[c] +=
          count / totalppc; // there can be contributions from multiple objects
      } else {
        vol_frac_CC[c] = ups_volFrac * count / (totalppc);
      }

      rho_micro[c] = getInitialDensity();
      rho_CC[c]    = rho_micro[c] * vol_frac_CC[c] + d_TINY_RHO;

      // these values of temp_CC and vel_CC are only used away from the mpm
      // objects
      // on the first timestep in interpolateNC_CC_0.  We just need reasonable
      // values
      temp[c] = 300.0;

      Point pd = patch->cellPosition(c);

      bool inside_bb_lo =
        (pd.x() > bb_low.x() && pd.y() > bb_low.y() && pd.z() > bb_low.z());
      bool inside_bb_up =
        (pd.x() < bb_up.x() && pd.y() < bb_up.y() && pd.z() < bb_up.z());

      // warning:  If two objects share the same cell then the last object sets
      // the values.
      //  This isn't a big deal since interpolateNC_CC_0 will only use these
      //  values on the first
      //  timmestep
      if (inside_bb_lo && inside_bb_up) {
        vel_CC[c] = d_geom_obj->getInitialData_Vector("velocity");
        temp[c]   = d_geom_obj->getInitialData_double("temperature");
      }

    } // Loop over cells
  }   // Loop over geom_objects
}
//______________________________________________________________________
//
void
MPMMaterial::initializeDummyCCVariables(CCVariable<double>& rho_micro,
                                        CCVariable<double>& rho_CC,
                                        CCVariable<double>& temp,
                                        CCVariable<Vector>& vel_CC,
                                        CCVariable<double>& vol_frac_CC,
                                        const Patch*)
{
  vel_CC.initialize(Vector(0., 0., 0.));
  rho_micro.initialize(d_density);
  rho_CC.initialize(d_TINY_RHO);
  temp.initialize(d_troom);
  vol_frac_CC.initialize(1.0);
}

} // end namespace Uintah