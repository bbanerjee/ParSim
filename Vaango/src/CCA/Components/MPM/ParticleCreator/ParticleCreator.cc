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

#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ScalarDiffusionModel.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>

#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Components/MPM/MMS/MMS.h>

#include <CCA/Components/MPM/PhysicalBC/CrackBC.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/DynamicSDFGeometry.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/CylinderGeometryPiece.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <fstream>
#include <iostream>
#include <memory>

namespace Uintah {

ParticleCreator::ParticleCreator(MPMMaterial* matl, MPMFlags* flags)
{
  d_matl                 = matl;
  d_mpm_labels          = std::make_unique<MPMLabel>();
  d_amrmpm_labels       = std::make_unique<AMRMPMLabel>();
  d_hydrompm_labels     = std::make_unique<HydroMPMLabel>();
  d_useLoadCurves       = flags->d_useLoadCurves;
  d_useLoadCurvesVector = flags->d_useLoadCurvesVector;
  d_withColor           = flags->d_withColor;
  d_artificialViscosity = flags->d_artificialViscosity;
  d_computeScaleFactor  = flags->d_computeScaleFactor;
  d_doScalarDiffusion   = flags->d_doScalarDiffusion;
  d_useCPTI             = flags->d_useCPTI;
  d_flags               = flags;

  registerPermanentParticleState(matl);
}

auto
ParticleCreator::createParticles(MPMMaterial* matl,
                                 CCVariable<short int>& cellNAPID,
                                 const Patch* patch,
                                 DataWarehouse* new_dw,
                                 const VecGeometryObjectSP& geom_objs)
  -> particleIndex
{
  ObjectVars obj_vars;
  particleIndex numParticles = 0;
  for (const auto& geom : geom_objs) {
    numParticles += countAndCreateParticles(patch, geom.get(), obj_vars);
  }

  ParticleVars pvars;
  int dwi = matl->getDWIndex();
  [[maybe_unused]] ParticleSubset* pset =
    allocateVariables(numParticles, dwi, patch, new_dw, pvars);

  particleIndex current_particle_index = 0;

  for (const auto& obj_sp : geom_objs) {
    GeometryObject* obj_ptr = obj_sp.get();
    GeometryPieceP piece = obj_sp->getPiece();

    Box b = piece->getBoundingBox().intersect(patch->getExtraBox());
    if (b.degenerate()) {
      continue;
    }

    Vector dxpp = patch->dCell() / obj_sp->getInitialData_IntVector("res");

    // Special case exception for SpecialGeomPieces
    // (includes FileGP, SmoothGP, AbaqusGP, CorrugatedGP etc.)
    auto* sgp = dynamic_cast<SpecialGeomPiece*>(piece.get());

    // Set up the local SpecialGeometryPiece particle data
    SGPData<double> sgp_vol_data;
    SGPData<double> sgp_temp_data;
    SGPData<Vector> sgp_force_data;
    SGPData<Vector> sgp_fiber_data;
    SGPData<Vector> sgp_velocity_data;
    SGPData<Matrix3> sgp_size_data;
    SGPData<double> sgp_color_data;
    SGPData<double> sgp_concentration_data;
    SGPData<Vector> sgp_area_data;
    SGPData<double> sgp_poscharge_data;
    SGPData<double> sgp_negcharge_data;
    SGPData<double> sgp_permittivity_data;

    if (sgp) {

      std::cout << std::format(
        "Created a special geometry of type {} with #particles = {}\n",
        piece->getType(), numParticles);

      sgp_vol_data.initialize(sgp, obj_vars, "p.volume", obj_ptr);
      sgp_temp_data.initialize(sgp, obj_vars, "p.temperature", obj_ptr);
      sgp_force_data.initialize(sgp, obj_vars, "p.externalforce", obj_ptr);
      sgp_fiber_data.initialize(sgp, obj_vars, "p.fiberdirs", obj_ptr);
      sgp_velocity_data.initialize(sgp, obj_vars, "p.velocity", obj_ptr);
      sgp_size_data.initialize(sgp, obj_vars, "p.size", obj_ptr);

      if (d_withColor) {
        sgp_color_data.initialize(sgp, obj_vars, "p.color", obj_ptr);
      }

      if (d_doScalarDiffusion) {
        sgp_concentration_data.initialize(
          sgp, obj_vars, "p.concentration", obj_ptr);
        sgp_area_data.initialize(sgp, obj_vars, "p.area", obj_ptr);
      }

      if (d_withGaussSolver) {
        sgp_poscharge_data.initialize(sgp, obj_vars, "p.poscharge", obj_ptr);
        sgp_negcharge_data.initialize(sgp, obj_vars, "p.negcharge", obj_ptr);
        sgp_permittivity_data.initialize(
          sgp, obj_vars, "p.permittivity", obj_ptr);
      }

    } else {
      std::cout << "Created a geometry with #particles = " << numParticles
                << std::endl;
    }

    Vector eps(1e-10, 1e-10, 1e-10);
    for (const auto& point : obj_vars.points.at(obj_ptr)) {
      IntVector cell_idx;
      if (!patch->findCell(point, cell_idx)) {
        if (!patch->findCell(point - eps, cell_idx)) {
          continue;
        }
      }
      if (!patch->containsPoint(point) && !patch->containsPoint(point - eps)) {
        continue;
      }

      particleIndex pidx = current_particle_index++;

      //std::cout << std::format("Point[{}] = {} Cell = {}\n", pidx, point, cell_idx);
      initializeParticle(
        patch, obj_ptr, matl, point, cell_idx, pidx, cellNAPID, pvars);

      // Again, everything below exists for SpecialGeometryPiece only
      if (sgp) {

        if (auto vol = sgp_vol_data.get_and_advance()) {
          pvars.pVolume[pidx] = *vol;
          pvars.pMass[pidx]   = matl->getInitialDensity() * pvars.pVolume[pidx];
        }

        if (auto temp = sgp_temp_data.get_and_advance()) {
          pvars.pTemperature[pidx] = *temp;
        }

        if (auto force = sgp_force_data.get_and_advance()) {
          pvars.pExternalForce[pidx] = *force;
        }

        if (auto vel = sgp_velocity_data.get_and_advance()) {
          pvars.pVelocity[pidx] = *vel;
        }

        if (auto fiber_dir = sgp_fiber_data.get_and_advance()) {
          pvars.pFiberDir[pidx] = *fiber_dir;
        }

        if (auto size_val =
              sgp_size_data.get_and_advance()) { // Renamed 'size' to 'size_val'
                                                 // to avoid conflict
          Vector dxcc       = patch->dCell();
          pvars.pSize[pidx] = *size_val;

          // Calculate volume based on type
          if (!d_useCPTI) { // CPDI and others
            pvars.pVolume[pidx] = std::abs(pvars.pSize[pidx].Determinant());
          } else { // CPTI
            pvars.pVolume[pidx] =
              std::abs(pvars.pSize[pidx].Determinant() / 6.0);
          }
          pvars.pMass[pidx] = matl->getInitialDensity() * pvars.pVolume[pidx];

          // Modify pSize (R-vectors) normalized by cell spacing
          Matrix3 cell_size_norm(1.0 / dxcc.x(),
                                 0.0,
                                 0.0,
                                 0.0,
                                 1.0 / dxcc.y(),
                                 0.0,
                                 0.0,
                                 0.0,
                                 1.0 / dxcc.z());
          pvars.pSize[pidx] = pvars.pSize[pidx] * cell_size_norm;
        }

        if (d_withColor) {
          if (auto color = sgp_color_data.get_and_advance()) {
            pvars.pColor[pidx] = *color;
          }
        }

        if (d_doScalarDiffusion) {
          if (auto area = sgp_area_data.get_and_advance()) {
            pvars.pArea[pidx] = *area;
          }
          if (auto conc = sgp_concentration_data.get_and_advance()) {
            pvars.pConcentration[pidx] = *conc;
          }
        }

        if (d_withGaussSolver) {
          if (auto pos_charge = sgp_poscharge_data.get_and_advance()) {
            pvars.pPosCharge[pidx] = *pos_charge;
          }
          if (auto neg_charge = sgp_negcharge_data.get_and_advance()) {
            pvars.pNegCharge[pidx] = *neg_charge;
          }
          if (auto permittivity = sgp_permittivity_data.get_and_advance()) {
            pvars.pPermittivity[pidx] = *permittivity;
          }
        }
      }

      // If the particle is on the surface and if there is
      // a physical BC attached to it then mark with the
      // physical BC pointer
      if (d_useLoadCurves) {
        pvars.pLoadCurveID[pidx] =
          checkForSurface(piece, point, dxpp) ? getLoadCurveID(point, dxpp) : 0;
      }
      //for (const auto pidx : *pset) {
      //  std::cout << std::format(
      //    "pidx = {}, volume = {}\n", pidx, pvars.pVolume[pidx]);
      //}
    }
  }

  //for (const auto pidx: *pset) {
  //  std::cout << std::format("pidx = {}, volume = {}\n", pidx, pvars.pVolume[pidx]);
  //}
  return numParticles;
}

// Get the LoadCurveID applicable for this material point
// WARNING : Should be called only once per particle during a simulation
// because it updates the number of particles to which a BC is applied.
auto
ParticleCreator::getLoadCurveID(const Point& pp, const Vector& dxpp) -> int
{
  int ret = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();

    // std::cerr << " BC Type = " << bcType << std::endl;
    if (bcType == "Pressure") {
      auto* pbc = dynamic_cast<PressureBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
        // std::cout << "\t surface particle; flagged material pt" << std::endl;
        ret = pbc->loadCurveID();
      }
    } else if (bcType == "Velocity") {
      auto* vbc = dynamic_cast<VelocityBC*>(bc.get());
      if (vbc->flagMaterialPoint(pp, dxpp)) {
        // std::cout << "\t surface particle; flagged material pt" << std::endl;
        ret = vbc->loadCurveID();
      }
    } else if (bcType == "Moment") {
      auto* pbc = dynamic_cast<MomentBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
        ret = pbc->loadCurveID();
      }
    } else if (bcType == "HeatFlux") {
      auto* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
      if (hfbc->flagMaterialPoint(pp, dxpp)) {
        ret = hfbc->loadCurveID();
      }
    }
  }
  return ret;
}

// Print MPM physical boundary condition information
void
ParticleCreator::printPhysicalBCs()
{
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();
    if (bcType == "Pressure") {
      auto* pbc = dynamic_cast<PressureBC*>(bc.get());
      std::cerr << *pbc << std::endl;
    }
    if (bcType == "Velocity") {
      auto* vbc = dynamic_cast<VelocityBC*>(bc.get());
      std::cerr << *vbc << std::endl;
    }
    if (bcType == "Moment") {
      auto* pbc = dynamic_cast<MomentBC*>(bc.get());
      std::cerr << *pbc << std::endl;
    }
    if (bcType == "HeatFlux") {
      auto* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
      std::cerr << *hfbc << std::endl;
    }
  }
}

void
ParticleCreator::applyForceBC(const Vector& dxpp,
                              const Point& pp,
                              const double& pMass,
                              Vector& pExtForce)
{
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();

    // std::cerr << " BC Type = " << bcType << std::endl;
    if (bcType == "Force") {
      auto* fbc = dynamic_cast<ForceBC*>(bc.get());

      Box fbcBox(fbc->getLowerRange() - dxpp, fbc->getUpperRange() + dxpp);

      // std::cerr << "BC Box = " << bcBox << " Point = " << pp << std::endl;
      if (fbcBox.contains(pp)) {
        pExtForce = fbc->getForceDensity() * pMass;
        // std::cerr << "External Force on Particle = " << pExtForce
        //      << " Force Density = " << bc->getForceDensity()
        //      << " Particle Mass = " << pMass << std::endl;
      }
    }
  }
}

auto
ParticleCreator::allocateVariables(particleIndex numParticles,
                                   int dwi,
                                   const Patch* patch,
                                   DataWarehouse* new_dw,
                                   ParticleVars& pvars) -> ParticleSubset*
{

  ParticleSubset* subset =
    new_dw->createParticleSubset(numParticles, dwi, patch);
  new_dw->allocateAndPut(pvars.position, d_mpm_labels->pXLabel, subset);
  new_dw->allocateAndPut(pvars.pDisp, d_mpm_labels->pDispLabel, subset);
  new_dw->allocateAndPut(pvars.pVelocity, d_mpm_labels->pVelocityLabel, subset);
  new_dw->allocateAndPut(pvars.pAcc, d_mpm_labels->pAccelerationLabel, subset);
  new_dw->allocateAndPut(
    pvars.pExternalForce, d_mpm_labels->pExternalForceLabel, subset);
  new_dw->allocateAndPut(pvars.pMass, d_mpm_labels->pMassLabel, subset);
  new_dw->allocateAndPut(pvars.pVolume, d_mpm_labels->pVolumeLabel, subset);
  new_dw->allocateAndPut(
    pvars.pTemperature, d_mpm_labels->pTemperatureLabel, subset);
  new_dw->allocateAndPut(
    pvars.pParticleID, d_mpm_labels->pParticleIDLabel, subset);

  if (d_flags->d_enableDEM) {
    new_dw->allocateAndPut(pvars.pX0, d_mpm_labels->pX0Label, subset);
    new_dw->allocateAndPut(
      pvars.pRigidBodyID, d_mpm_labels->pRigidBodyIDLabel, subset);
    new_dw->allocateAndPut(
      pvars.pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, subset);
    new_dw->allocateAndPut(pvars.pTorque, d_mpm_labels->pTorqueLabel, subset);
    new_dw->allocateAndPut(
      pvars.pOrientation, d_mpm_labels->pOrientationLabel, subset);
    new_dw->allocateAndPut(pvars.pRadius, d_mpm_labels->pRadiusLabel, subset);
    new_dw->allocateAndPut(
      pvars.pInertiaTensor, d_mpm_labels->pInertiaTensorLabel, subset);
  }

  new_dw->allocateAndPut(pvars.pSize, d_mpm_labels->pSizeLabel, subset);
  new_dw->allocateAndPut(pvars.pFiberDir, d_mpm_labels->pFiberDirLabel, subset);
  // for thermal stress
  new_dw->allocateAndPut(
    pvars.pTempPrevious, d_mpm_labels->pTempPreviousLabel, subset);

  if (d_useLoadCurves) {
    new_dw->allocateAndPut(
      pvars.pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, subset);
  }
  if (d_withColor) {
    new_dw->allocateAndPut(pvars.pColor, d_mpm_labels->pColorLabel, subset);
  }
  if (d_artificialViscosity) {
    new_dw->allocateAndPut(pvars.p_q, d_mpm_labels->p_qLabel, subset);
  }

  // For AMR
  new_dw->allocateAndPut(pvars.pRefined, d_mpm_labels->pRefinedLabel, subset);
  if (d_flags->d_AMR) {
    new_dw->allocateAndPut(
      pvars.pLastLevel, d_mpm_labels->pLastLevelLabel, subset);
  }

  // For body force calculation
  new_dw->allocateAndPut(
    pvars.pBodyForceAcc, d_mpm_labels->pBodyForceAccLabel, subset);
  new_dw->allocateAndPut(
    pvars.pCoriolisImportance, d_mpm_labels->pCoriolisImportanceLabel, subset);

  // For switching between implicit and explicit
  new_dw->allocateAndPut(
    pvars.pExternalHeatFlux, d_mpm_labels->pExternalHeatFluxLabel, subset);

  // For friction contact
  new_dw->allocateAndPut(pvars.pSurface, d_mpm_labels->pSurfLabel, subset);

  return subset;
}

void
ParticleCreator::createPoints(const Patch* patch,
                              GeometryObject* obj,
                              ObjectVars& obj_vars)
{
  GeometryPieceP piece = obj->getPiece();
  Box b2               = patch->getExtraBox();
  IntVector ppc        = obj->getInitialData_IntVector("res");
  Vector dxpp          = patch->dCell() / ppc;
  Vector dcorner       = dxpp * 0.5;

  // Affine transformation for making conforming particle distributions
  // to be used in the conforming CPDI simulations. The input vectors are
  // optional and if you do not wish to use the affine transformation, just do
  // not define them in the input file.
  Vector affineTrans_A0 = obj->getInitialData_Vector("affineTransformation_A0");
  Vector affineTrans_A1 = obj->getInitialData_Vector("affineTransformation_A1");
  Vector affineTrans_A2 = obj->getInitialData_Vector("affineTransformation_A2");
  Vector affineTrans_b  = obj->getInitialData_Vector("affineTransformation_b");
  Matrix3 affineTrans_A(affineTrans_A0[0],
                        affineTrans_A0[1],
                        affineTrans_A0[2],
                        affineTrans_A1[0],
                        affineTrans_A1[1],
                        affineTrans_A1[2],
                        affineTrans_A2[0],
                        affineTrans_A2[1],
                        affineTrans_A2[2]);

  // AMR stuff
  const Level* curLevel = patch->getLevel();
  bool hasFiner         = curLevel->hasFinerLevel();
  Level* fineLevel      = nullptr;
  if (hasFiner) {
    fineLevel = (Level*)curLevel->getFinerLevel().get_rep();
  }

  // Actual computation
  for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
    Point lower = patch->nodePosition(*iter) + dcorner;
    IntVector c = *iter;

    if (hasFiner) { // Don't create particles if a finer level exists here
      const Point CC         = patch->cellPosition(c);
      bool includeExtraCells = false;
      const Patch* patchExists =
        fineLevel->getPatchFromPoint(CC, includeExtraCells);
      if (patchExists != nullptr) {
        continue;
      }
    }

    for (int ix = 0; ix < ppc.x(); ix++) {
      for (int iy = 0; iy < ppc.y(); iy++) {
        for (int iz = 0; iz < ppc.z(); iz++) {

          IntVector idx(ix, iy, iz);
          Point p = lower + dxpp * idx;
          if (!b2.contains(p)) {
            throw InternalError(
              "Particle created outside of patch?", __FILE__, __LINE__);
          }
          if (piece->inside(p)) {
            Vector p1(p(0), p(1), p(2));
            p1   = affineTrans_A * p1 + affineTrans_b;
            p(0) = p1[0];
            p(1) = p1[1];
            p(2) = p1[2];

            // Adds a key to obj_vars if it doesn't exist
            obj_vars.points[obj].push_back(p);
          }
        } // z
      }   // y
    }     // x
  }       // iterator

  /*
  //  *** This part is associated with CBDI_CompressiveCylinder.ups input file.
  //      It creates conforming particle distribution to be used in the
  simulation.
  //      To use that you need to uncomment the following commands to create the
  //      conforming particle distribution and comment above commands that are
  used
  //      to create non-conforming particle distributions.

    geompoints::key_type key(patch,obj);
    int resolutionPart=1;
    int nPar1=180*resolutionPart;
    int nPar2=16*resolutionPart;

    for(int ix=1;ix < nPar1+1; ix++){
      double ttemp,rtemp;
      ttemp=(ix-0.5)*2.0*3.14159265358979/nPar1;
      for(int iy=1;iy < nPar2+1; iy++){
          rtemp=0.75+(iy-0.5)*0.5/nPar2;
          Point p(rtemp*cos(ttemp),rtemp*sin(ttemp),0.5);
          if(patch->containsPoint(p)){
            object_points[key].push_back(p);
          }
      }
    }
  */
}

void
ParticleCreator::initializeParticle(const Patch* patch,
                                    GeometryObject* obj,
                                    MPMMaterial* matl,
                                    Point p,
                                    IntVector cell_idx,
                                    particleIndex i,
                                    CCVariable<short int>& cellNAPID,
                                    ParticleVars& pvars)
{
  IntVector ppc = obj->getInitialData_IntVector("res");
  Vector dxpp   = patch->dCell() / obj->getInitialData_IntVector("res");
  Vector dxcc   = patch->dCell();

  // Affine transformation for making conforming particle distributions
  //  to be used in the conforming CPDI simulations. The input vectors are
  //  optional and if you do not want to use affine transformations, just do
  //  not define them in the input file.

  Vector affineTrans_A0 = obj->getInitialData_Vector("affineTransformation_A0");
  Vector affineTrans_A1 = obj->getInitialData_Vector("affineTransformation_A1");
  Vector affineTrans_A2 = obj->getInitialData_Vector("affineTransformation_A2");
  // Vector affineTrans_b= obj->getInitialData_Vector("affineTransformation_b");
  Matrix3 affineTrans_A(affineTrans_A0[0],
                        affineTrans_A0[1],
                        affineTrans_A0[2],
                        affineTrans_A1[0],
                        affineTrans_A1[1],
                        affineTrans_A1[2],
                        affineTrans_A2[0],
                        affineTrans_A2[1],
                        affineTrans_A2[2]);
  // The size matrix is used for storing particle domain sizes
  // (Rvectors for CPDI and CPTI) normalized by the grid spacing
  Matrix3 size(1. / ((double)ppc.x()),
               0.,
               0.,
               0.,
               1. / ((double)ppc.y()),
               0.,
               0.,
               0.,
               1. / ((double)ppc.z()));
  size = affineTrans_A * size;

  /*
  //  *** This part is associated with CBDI_CompressiveCylinder.ups input file.
  //      It determines particle domain sizes for the conforming particle
  distribution,
  //      which is used in the simulation.
  //      To activate that you need to uncomment the following commands to
  determine
  //      particle domain sizes for the conforming particle distribution, and
  //      comment above commands that are used to determine particle domain
  sizes for
  //      non-conforming particle distributions.

    int resolutionPart=1;
    int nPar1=180*resolutionPart;
    int nPar2=16*resolutionPart;
    double pi=3.14159265358979;
    double rtemp=sqrt(p.x()*p.x()+p.y()*p.y());
    Matrix3 size(2.*pi/nPar1*p.y()/dxcc[0],2.*0.25/nPar2/rtemp*p.x()/dxcc[1],0.,
                -2.*pi/nPar1*p.x()/dxcc[0],2.*0.25/nPar2/rtemp*p.y()/dxcc[1],0.,
                                        0., 0.,1.);
  */

  pvars.pTemperature[i] = obj->getInitialData_double("temperature");

  // For AMR
  const Level* curLevel = patch->getLevel();
  pvars.pRefined[i]     = curLevel->getIndex();

  // MMS
  // std::cout << "mms_type = " << d_flags->d_mmsType << "\n";
  string mms_type = d_flags->d_mmsType;
  if (mms_type != "none") {
    MMS MMSObject;
    MMSObject.initializeParticleForMMS(pvars.position,
                                       pvars.pVelocity,
                                       pvars.pSize,
                                       pvars.pDisp,
                                       pvars.pMass,
                                       pvars.pVolume,
                                       p,
                                       dxcc,
                                       size,
                                       patch,
                                       d_flags,
                                       i);
  } else {
    if (std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z())) {
      std::ostringstream err;
      err << "**ERROR** Particle created with NaN position: " << p << "\n"
          << "  Geometry Object index: " << matl->getGeometryObjectIndex(obj) << "\n";
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
    pvars.position[i] = p;
    if (d_flags->d_axisymmetric) {
      // assume unit radian extent in the circumferential direction
      pvars.pVolume[i] = p.x() *
                         (size(0, 0) * size(1, 1) - size(0, 1) * size(1, 0)) *
                         dxcc.x() * dxcc.y();
      // std::cout << "dx_cc = " << dxcc << " size = " << size
      //           << " x = " << p.x() << " vol = " << pvars.pVolume[i] << "\n";
    } else {
      // standard voxel volume
      pvars.pVolume[i] = size.Determinant() * dxcc.x() * dxcc.y() * dxcc.z();
    }

    pvars.pSize[i]     = size;
    pvars.pDisp[i]     = Vector(0., 0., 0.);
    pvars.pVelocity[i] = obj->getInitialData_Vector("velocity");
    pvars.pAcc[i]      = Vector(0., 0., 0.);

    // std::cout << "i = " << i << " px = " << pvars.position[i]
    //           << " size = " << pvars.pSize[i] << "\n";
    double vol_frac_CC = 1.0;
    try {
      if (obj->getInitialData_double("volumeFraction") == -1.0) {
        vol_frac_CC    = 1.0;
        pvars.pMass[i] = matl->getInitialDensity() * pvars.pVolume[i];
      } else {
        vol_frac_CC = obj->getInitialData_double("volumeFraction");
        pvars.pMass[i] =
          matl->getInitialDensity() * pvars.pVolume[i] * vol_frac_CC;
      }
    } catch (...) {
      vol_frac_CC    = 1.0;
      pvars.pMass[i] = matl->getInitialDensity() * pvars.pVolume[i];
    }
  }

  if (d_withColor) {
    pvars.pColor[i] = obj->getInitialData_double("color");
  }
  if (d_artificialViscosity) {
    pvars.p_q[i] = 0.;
  }
  if (d_flags->d_AMR) {
    pvars.pLastLevel[i] = curLevel->getID();
  }

  pvars.pTempPrevious[i] = pvars.pTemperature[i];

  Vector pExtForce(0, 0, 0);
  applyForceBC(dxpp, p, pvars.pMass[i], pExtForce);

  pvars.pExternalForce[i] = pExtForce;
  pvars.pFiberDir[i]      = matl->getConstitutiveModel()->getInitialFiberDir();

  // AMR
  if (d_flags->d_AMR) {
    pvars.pLastLevel[i] = curLevel->getID();
  }

  // For body forces
  pvars.pBodyForceAcc[i]       = Vector(0.0, 0.0, 0.0); // Init to zero
  pvars.pCoriolisImportance[i] = 0.0;

  // For switching between explicit and implicit MPM
  pvars.pExternalHeatFlux[i] = 0.0;

  // For friction contact
  GeometryPieceP piece = obj->getPiece();
  pvars.pSurface[i]    = checkForSurface2(piece, p, dxpp);

  if (d_flags->d_enableDEM) {
    pvars.pX0[i] = p;
    // Rigid Body ID
    int objIdx = matl->getGeometryObjectIndex(obj);
    long64 rbID =
      ((long64)matl->getDWIndex() << 32) | (long64)(objIdx >= 0 ? objIdx : 0);
    pvars.pRigidBodyID[i] = rbID;

    // Initialize other DEM variables
    pvars.pAngularVelocity[i] = Vector(0, 0, 0);
    pvars.pTorque[i]          = Vector(0, 0, 0);
    pvars.pOrientation[i]     = Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    pvars.pRadius[i] = matl->getParticleRadius();

    // For discrete material, check if Master or Slave
    if (d_matl->isDiscrete()) {
      Box box      = piece->getBoundingBox();
      Point center = box.lower() + (box.upper() - box.lower()) * 0.5;
      double distSq = (p - center).length2();
      if (distSq > 1.0e-9) {
        // Slave/Dummy particle
        pvars.pMass[i]   = 0.0;
        pvars.pVolume[i] = 0.0;
        pvars.pRadius[i] = 0.0;
      } else {
        // Master particle: use total volume of geometry piece if available
        double vol = 0;
        if (auto* sphere = dynamic_cast<SphereGeometryPiece*>(piece.get())) {
          vol              = sphere->volume();
          pvars.pRadius[i] = sphere->radius();
        } else if (auto* box_piece =
                     dynamic_cast<BoxGeometryPiece*>(piece.get())) {
          vol              = box_piece->volume();
          pvars.pRadius[i] = 0.5 * box_piece->smallestSide();
        } else if (auto* cyl =
                     dynamic_cast<CylinderGeometryPiece*>(piece.get())) {
          vol              = cyl->volume();
          pvars.pRadius[i] = cyl->radius();
        } else if (auto* sdf =
                     dynamic_cast<DynamicSDFGeometry*>(piece.get())) {
          Box b = sdf->getBoundingBox();
          Vector diag = b.upper() - b.lower();
          vol = diag.x() * diag.y() * diag.z();
          pvars.pRadius[i] = 0.5 * std::min({diag.x(), diag.y(), diag.z()});
        } else {
          // Fallback for other geometry pieces
          Box b = piece->getBoundingBox();
          Vector diag = b.upper() - b.lower();
          vol = diag.x() * diag.y() * diag.z();
          if (pvars.pRadius[i] == 0) {
            pvars.pRadius[i] = matl->getParticleRadius();
          }
        }
        if (vol > 0) {
          pvars.pVolume[i] = vol;
          pvars.pMass[i]   = matl->getInitialDensity() * vol;
        }
      }
    }

    double I = 0.4 * pvars.pMass[i] * pvars.pRadius[i] * pvars.pRadius[i];
    pvars.pInertiaTensor[i] = Matrix3(I, 0, 0, 0, I, 0, 0, 0, I);
  }

  // Cell ids
  ASSERT(cell_idx.x() <= 0xffff && cell_idx.y() <= 0xffff &&
         cell_idx.z() <= 0xffff);

  long64 cellID = ((long64)cell_idx.x() << 16) | ((long64)cell_idx.y() << 32) |
                  ((long64)cell_idx.z() << 48);

  short int& myCellNAPID = cellNAPID[cell_idx];
  pvars.pParticleID[i]   = (cellID | (long64)myCellNAPID);
  ASSERT(myCellNAPID < 0x7fff);
  myCellNAPID++;
}

auto
ParticleCreator::countAndCreateParticles(const Patch* patch,
                                         GeometryObject* obj,
                                         ObjectVars& obj_vars) -> particleIndex
{
  GeometryPieceP piece = obj->getPiece();
  Box b1               = piece->getBoundingBox();
  Box b2               = patch->getExtraBox();
  Box b                = b1.intersect(b2);
  if (b.degenerate()) {
    return 0;
  }

  // Create pairs for the particle variables created by
  // the SpecialGeometryPieces
  auto vol_pair   = std::make_pair("p.volume", obj);
  auto temp_pair  = std::make_pair("p.temperature", obj);
  auto color_pair = std::make_pair("p.color", obj);
  auto conc_pair  = std::make_pair("p.concentration", obj);
  auto pos_pair   = std::make_pair("p.poscharge", obj);
  auto neg_pair   = std::make_pair("p.negcharge", obj);
  auto perm_pair  = std::make_pair("p.permittivity", obj);
  auto force_pair = std::make_pair("p.externalforce", obj);
  auto fiber_pair = std::make_pair("p.fiber", obj);
  auto vel_pair   = std::make_pair("p.velocity", obj);
  auto area_pair  = std::make_pair("p.area", obj);
  auto size_pair  = std::make_pair("p.size", obj);

  // If the material is discrete, we want exactly one master particle
  // at the center of the object and several slave particles as geometric proxies.
  if (d_matl->isDiscrete()) {
    Box box      = piece->getBoundingBox();
    Point center = box.lower() + (box.upper() - box.lower()) * 0.5;

    Vector eps(1e-10, 1e-10, 1e-10);

    // Master particle at center
    Point m_center = center;
    if (!patch->containsPoint(m_center)) {
      if (patch->containsPoint(m_center - eps)) {
        m_center = m_center - eps * 0.1;
      } else if (patch->containsPoint(m_center + eps)) {
        m_center = m_center + eps * 0.1;
      }
    }

    if (patch->containsPoint(m_center)) {
      obj_vars.points[obj].push_back(m_center);
    }

    // Add surface dummy particles for visualization.
    // These particles have zero mass and do not participate in contact,
    // but follow the rigid body's motion.
    if (auto* cyl = dynamic_cast<CylinderGeometryPiece*>(piece.get())) {
      Point bottom = cyl->bottom();
      Point top    = cyl->top();
      double rad   = cyl->radius();
      Vector axis  = top - bottom;
      double h     = axis.length();
      axis.normalize();

      // Find two vectors orthogonal to the axis
      Vector v1, v2;
      axis.findOrthogonal(v1, v2);

      int n_theta = 16;
      int n_z     = 8;
      for (int i_z = 0; i_z <= n_z; ++i_z) {
        double z = (double)i_z / (double)n_z * h;
        Point p_z = bottom + axis * z;
        for (int i_t = 0; i_t < n_theta; ++i_t) {
          double theta = (double)i_t / (double)n_theta * 2.0 * M_PI;
          Point p = p_z + (v1 * cos(theta) + v2 * sin(theta)) * rad;
          if (patch->containsPoint(p)) {
            obj_vars.points[obj].push_back(p);
          }
        }
      }
      // Caps
      int n_r = 3;
      for (int i_r = 1; i_r < n_r; ++i_r) {
        double r = (double)i_r / (double)n_r * rad;
        for (int i_t = 0; i_t < n_theta; ++i_t) {
          double theta = (double)i_t / (double)n_theta * 2.0 * M_PI;
          Vector radial = (v1 * cos(theta) + v2 * sin(theta)) * r;
          Point p_bot = bottom + radial;
          Point p_top = top + radial;
          if (patch->containsPoint(p_bot)) obj_vars.points[obj].push_back(p_bot);
          if (patch->containsPoint(p_top)) obj_vars.points[obj].push_back(p_top);
        }
      }
    } else if (auto* sphere = dynamic_cast<SphereGeometryPiece*>(piece.get())) {
      Point center = sphere->origin();
      double rad   = sphere->radius();
      int n_theta  = 16;
      int n_phi    = 8;
      for (int i_p = 1; i_p < n_phi; ++i_p) {
        double phi = (double)i_p / (double)n_phi * M_PI;
        double sin_phi = sin(phi);
        double cos_phi = cos(phi);
        for (int i_t = 0; i_t < n_theta; ++i_t) {
          double theta = (double)i_t / (double)n_theta * 2.0 * M_PI;
          Point p = center + Vector(rad * sin_phi * cos(theta),
                                    rad * sin_phi * sin(theta),
                                    rad * cos_phi);
          if (patch->containsPoint(p)) {
            obj_vars.points[obj].push_back(p);
          }
        }
      }
      // Poles
      Point p_north = center + Vector(0, 0, rad);
      Point p_south = center + Vector(0, 0, -rad);
      if (patch->containsPoint(p_north)) obj_vars.points[obj].push_back(p_north);
      if (patch->containsPoint(p_south)) obj_vars.points[obj].push_back(p_south);
    } else {
      // Fallback: Sparse grid points + Corners
      // (Original slave particle logic already handled corners below)
    }

    // Slave particles at corners (Geometric Proxies)
    // These ensure the object is visible to neighbor patches.
    Point corners[8];
    corners[0] = box.lower();
    corners[1] = Point(box.upper().x(), box.lower().y(), box.lower().z());
    corners[2] = Point(box.upper().x(), box.upper().y(), box.lower().z());
    corners[3] = Point(box.lower().x(), box.upper().y(), box.lower().z());
    corners[4] = Point(box.lower().x(), box.lower().y(), box.upper().z());
    corners[5] = Point(box.upper().x(), box.lower().y(), box.upper().z());
    corners[6] = box.upper();
    corners[7] = Point(box.lower().x(), box.upper().y(), box.upper().z());

    for (int k = 0; k < 8; ++k) {
      Point p = corners[k];
      Vector to_center = center - p;
      if (to_center.length() > 1e-9) {
        to_center.normalize();
        p = p + to_center * 1e-9;
      }

      if (patch->containsPoint(p) || 
          patch->containsPoint(p - eps)) {
        // Ensure uniqueness (important for 2D where corners overlap)
        bool exists = false;
        for (const auto& existing_p : obj_vars.points[obj]) {
          if ((existing_p - p).length2() < 1.0e-16) {
            exists = true;
            break;
          }
        }
        if (!exists) {
          obj_vars.points[obj].push_back(p);
        }
      }
    }

    return static_cast<particleIndex>(obj_vars.points[obj].size());
  }

  // If the object is a SpecialGeomPiece (e.g. FileGeometryPiece or
  // SmoothCylGeomPiece) then use the particle creators in that
  // class to do the counting
  auto* sgp = dynamic_cast<SpecialGeomPiece*>(piece.get());
  if (sgp) {
    int numPts = 0;
    auto* fgp  = dynamic_cast<FileGeometryPiece*>(piece.get());
    if (fgp) {
      if (d_useCPTI) {
        proc0cout << "*** Reading CPTI file ***" << std::endl;
      }
      fgp->readPoints(patch->getID());
      numPts = fgp->returnPointCount();
      proc0cout << "Number of points read from file = " << numPts << std::endl;
    } else {
      Vector dxpp = patch->dCell() / obj->getInitialData_IntVector("res");
      double dx   = Min(Min(dxpp.x(), dxpp.y()), dxpp.z());
      sgp->setParticleSpacing(dx);
      sgp->setCellSize(patch->dCell());
      numPts = sgp->createPoints();
      proc0cout << "Special Geom Piece: Number of points created = " << numPts
                << std::endl;
    }

    auto points         = sgp->getPoints();
    auto pVols          = sgp->getScalar("p.volume");
    auto pTemps         = sgp->getScalar("p.temperature");
    auto pColors        = sgp->getScalar("p.color");
    auto pConcentration = sgp->getScalar("p.concentration");
    auto pPosCharge     = sgp->getScalar("p.poscharge");
    auto pNegCharge     = sgp->getScalar("p.negcharge");
    auto pPermittivity  = sgp->getScalar("p.permittivity");
    auto pForces        = sgp->getVector("p.externalforce");
    auto pFiberDirs     = sgp->getVector("p.fiberdir");
    auto pVelocities    = sgp->getVector("p.velocity");
    auto pAreas         = sgp->getVector("p.area");
    auto pSizes         = sgp->getTensor("p.size");

    Point p;
    IntVector cell_idx;

    for (int ii = 0; ii < numPts; ++ii) {
      p = points->at(ii);
      //IntVector l(patch->getCellLowIndex());
      //IntVector h(patch->getCellHighIndex());
      //IntVector c = patch->getLevel()->getCellIndex(p);
      //std::cout << std::format("Point = {}, cell index = {}, low = {}, hi = {}\n", p, c, l, h);
      if (patch->findCell(p, cell_idx)) {
        if (patch->containsPoint(p)) {
          obj_vars.points[obj].push_back(p);

          if (pVols) {
            obj_vars.scalars[vol_pair].push_back(pVols->at(ii));
          }
          if (pTemps) {
            obj_vars.scalars[temp_pair].push_back(pTemps->at(ii));
          }
          if (pColors) {
            obj_vars.scalars[color_pair].push_back(pColors->at(ii));
          }
          if (pConcentration) {
            obj_vars.scalars[conc_pair].push_back(pConcentration->at(ii));
          }
          if (pPosCharge) {
            obj_vars.scalars[pos_pair].push_back(pPosCharge->at(ii));
          }
          if (pNegCharge) {
            obj_vars.scalars[neg_pair].push_back(pNegCharge->at(ii));
          }
          if (pPermittivity) {
            obj_vars.scalars[perm_pair].push_back(pPermittivity->at(ii));
          }
          if (pForces) {
            obj_vars.vectors[force_pair].push_back(pForces->at(ii));
          }
          if (pFiberDirs) {
            obj_vars.vectors[fiber_pair].push_back(pFiberDirs->at(ii));
          }
          if (pVelocities) {
            obj_vars.vectors[vel_pair].push_back(pVelocities->at(ii));
          }
          if (pAreas) {
            obj_vars.vectors[area_pair].push_back(pAreas->at(ii));
          }
          if (pSizes) {
            obj_vars.tensors[size_pair].push_back(pSizes->at(ii));
          }
        }
      } // patch contains cell
    }
    // sgp->deletePoints();
    // sgp->deleteVolume();
  } else {
    createPoints(patch, obj, obj_vars);
  }

  return static_cast<particleIndex>(obj_vars.points[obj].size());
}

auto
ParticleCreator::returnParticleState() -> vector<const VarLabel*>
{
  return particle_state;
}

auto
ParticleCreator::returnParticleStatePreReloc() -> vector<const VarLabel*>
{
  return particle_state_preReloc;
}

void
ParticleCreator::registerPermanentParticleState(MPMMaterial* matl)
{
  particle_state.push_back(d_mpm_labels->pDispLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pDispLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pVelocityLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pVelocityLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pAccelerationLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pAccelerationLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pExternalForceLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pExtForceLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pMassLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pMassLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pVolumeLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pVolumeLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pTemperatureLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pTemperatureLabel_preReloc);

  // for thermal stress
  particle_state.push_back(d_mpm_labels->pTempPreviousLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pTempPreviousLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pParticleIDLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pParticleIDLabel_preReloc);

  if (d_withColor) {
    particle_state.push_back(d_mpm_labels->pColorLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pColorLabel_preReloc);
  }

  particle_state.push_back(d_mpm_labels->pSizeLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pSizeLabel_preReloc);

  if (d_useLoadCurves) {
    particle_state.push_back(d_mpm_labels->pLoadCurveIDLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pLoadCurveIDLabel_preReloc);
  }

  particle_state.push_back(d_mpm_labels->pDispGradLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pDispGradLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pDefGradLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pDefGradLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pStressLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pStressLabel_preReloc);

  if (d_artificialViscosity) {
    particle_state.push_back(d_mpm_labels->p_qLabel);
    particle_state_preReloc.push_back(d_mpm_labels->p_qLabel_preReloc);
  }

  if (d_computeScaleFactor) {
    particle_state.push_back(d_mpm_labels->pScaleFactorLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pScaleFactorLabel_preReloc);
  }

  matl->getConstitutiveModel()->addParticleState(particle_state,
                                                 particle_state_preReloc);

  if (matl->doBasicDamage()) {
    matl->getBasicDamageModel()->addParticleState(particle_state,
                                                  particle_state_preReloc);
  }

  // For AMR
  if (d_flags->d_refineParticles) {
    particle_state.push_back(d_mpm_labels->pRefinedLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pRefinedLabel_preReloc);
  }

  if (d_flags->d_AMR) {
    particle_state.push_back(d_mpm_labels->pLastLevelLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pLastLevelLabel_preReloc);
  }

  // For body forces
  particle_state.push_back(d_mpm_labels->pBodyForceAccLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pBodyForceAccLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pCoriolisImportanceLabel);
  particle_state_preReloc.push_back(
    d_mpm_labels->pCoriolisImportanceLabel_preReloc);

  // For switching between implicit and explicit MPM
  particle_state.push_back(d_mpm_labels->pExternalHeatFluxLabel);
  particle_state.push_back(d_mpm_labels->pVelGradLabel);

  particle_state_preReloc.push_back(
    d_mpm_labels->pExternalHeatFluxLabel_preReloc);
  particle_state_preReloc.push_back(d_mpm_labels->pVelGradLabel_preReloc);

  particle_state.push_back(d_mpm_labels->pSurfLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pSurfLabel_preReloc);

  if (d_flags->d_enableDEM) {
    particle_state.push_back(d_mpm_labels->pX0Label);
    particle_state_preReloc.push_back(d_mpm_labels->pX0Label_preReloc);

    particle_state.push_back(d_mpm_labels->pRigidBodyIDLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pRigidBodyIDLabel_preReloc);

    particle_state.push_back(d_mpm_labels->pAngularVelocityLabel);
    particle_state_preReloc.push_back(
      d_mpm_labels->pAngularVelocityLabel_preReloc);

    particle_state.push_back(d_mpm_labels->pTorqueLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pTorqueLabel_preReloc);

    particle_state.push_back(d_mpm_labels->pOrientationLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pOrientationLabel_preReloc);

    particle_state.push_back(d_mpm_labels->pRadiusLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pRadiusLabel_preReloc);

    particle_state.push_back(d_mpm_labels->pInertiaTensorLabel);
    particle_state_preReloc.push_back(d_mpm_labels->pInertiaTensorLabel_preReloc);
  }

  // For scalar diffusion
  if (d_doScalarDiffusion) {
    particle_state.push_back(d_mpm_labels->diffusion->pConcentration);
    particle_state_preReloc.push_back(
      d_mpm_labels->diffusion->pConcentration_preReloc);

    particle_state.push_back(d_mpm_labels->diffusion->pConcPrevious);
    particle_state_preReloc.push_back(
      d_mpm_labels->diffusion->pConcPrevious_preReloc);

    particle_state.push_back(d_mpm_labels->diffusion->pGradConcentration);
    particle_state_preReloc.push_back(
      d_mpm_labels->diffusion->pGradConcentration_preReloc);

    particle_state.push_back(d_mpm_labels->diffusion->pExternalScalarFlux);
    particle_state_preReloc.push_back(
      d_mpm_labels->diffusion->pExternalScalarFlux_preReloc);

    particle_state.push_back(d_mpm_labels->diffusion->pArea);
    particle_state_preReloc.push_back(d_mpm_labels->diffusion->pArea_preReloc);

    matl->getScalarDiffusionModel()->addParticleState(particle_state,
                                                      particle_state_preReloc);
  }

  if (d_coupledFlow && !matl->getIsRigid()) {
    // if (d_coupledflow ) {
    particle_state.push_back(d_hydrompm_labels->pFluidMassLabel);
    particle_state.push_back(d_hydrompm_labels->pSolidMassLabel);
    particle_state.push_back(d_hydrompm_labels->pPorePressureLabel);
    particle_state.push_back(d_hydrompm_labels->pPorosityLabel);

    // Error Cannot find in relocateParticles ???

    particle_state_preReloc.push_back(
      d_hydrompm_labels->pFluidMassLabel_preReloc);
    particle_state_preReloc.push_back(
      d_hydrompm_labels->pSolidMassLabel_preReloc);
    particle_state_preReloc.push_back(
      d_hydrompm_labels->pPorePressureLabel_preReloc);
    particle_state_preReloc.push_back(
      d_hydrompm_labels->pPorosityLabel_preReloc);

    if (d_flags->d_integratorType == "explicit") {
      particle_state.push_back(d_hydrompm_labels->pFluidVelocityLabel);
      particle_state_preReloc.push_back(
        d_hydrompm_labels->pFluidVelocityLabel_preReloc);
    }
  }
}

auto
ParticleCreator::checkForSurface(const GeometryPieceP piece,
                                 const Point p,
                                 const Vector dxpp) -> int
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface
  // std::cout << "GeometryPiece = " << piece->getType() << std::endl;
  // std::cout << " Point = " << p << " box = " << dxpp << std::endl;

  int ss = 0;
  // Check to the left (-x)
  if (!piece->inside(p - Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check to the right (+x)
  if (!piece->inside(p + Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check behind (-y)
  if (!piece->inside(p - Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
  // Check in front (+y)
  if (!piece->inside(p + Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
#if 1
  // Check below (-z)
  if (!piece->inside(p - Vector(0., 0., dxpp.z()))) {
    ss++;
  }
  // Check above (+z)
  if (!piece->inside(p + Vector(0., 0., dxpp.z()))) {
    ss++;
  }
#endif

  if (ss > 0) {
    return 1;
  } else {
    return 0;
  }
}

auto
ParticleCreator::checkForSurface2(const GeometryPieceP piece,
                                  const Point p,
                                  const Vector dxpp) -> double
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface

  int ss = 0;
  // Check to the left (-x)
  if (!piece->inside(p - Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check to the right (+x)
  if (!piece->inside(p + Vector(dxpp.x(), 0., 0.))) {
    ss++;
  }
  // Check behind (-y)
  if (!piece->inside(p - Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
  // Check in front (+y)
  if (!piece->inside(p + Vector(0., dxpp.y(), 0.))) {
    ss++;
  }
  // if (d_flags->d_ndim==3) {
  //  Check below (-z)
  if (!piece->inside(p - Vector(0., 0., dxpp.z()))) {
    ss++;
  }
  // Check above (+z)
  if (!piece->inside(p + Vector(0., 0., dxpp.z()))) {
    ss++;
  }
  //}

  if (ss > 0) {
    return 1.0;
  } else {
    return 0.0;
  }
}

} // end namespace Uintah