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

#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/Core/MPMMaterial.h>

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ScalarDiffusionModel.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>

#include <CCA/Components/MPM/Core/AMRMPMLabel.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMDiffusionLabel.h>

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
#include <Core/GeometryPiece/SpecialGeomPiece.h>

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

particleIndex
ParticleCreator::createParticles(MPMMaterial* matl,
                                 CCVariable<short int>& cellNAPID,
                                 const Patch* patch,
                                 DataWarehouse* new_dw,
                                 const VecGeometryObjectSP& geom_objs)
{
  ObjectVars obj_vars;
  particleIndex numParticles = 0;
  for (const auto& geom : geom_objs) {
    numParticles += countAndCreateParticles(patch, geom.get(), obj_vars);
  }

  ParticleVars pvars;
  int dwi = matl->getDWIndex();
  allocateVariables(numParticles, dwi, patch, new_dw, pvars);

  particleIndex start = 0;

  for (auto& obj : geom_objs) {
    particleIndex count  = 0;
    GeometryPieceP piece = obj->getPiece();
    Box b1               = piece->getBoundingBox();
    Box b2               = patch->getExtraBox();
    Box b                = b1.intersect(b2);
    if (b.degenerate()) {
      count = 0;
      continue;
    }

    Vector dxpp = patch->dCell() / obj->getInitialData_IntVector("res");

    // Special case exception for SpecialGeomPieces
    // (includes FileGP, SmoothGP, AbaqusGP, CorrugatedGP etc.)
    SpecialGeomPiece* sgp = dynamic_cast<SpecialGeomPiece*>(piece.get());

    // Set up pointers to SpecialGeometryPiece particle data
    const std::vector<double>* pVolumes        = nullptr;
    const std::vector<double>* pTemperatures   = nullptr;
    const std::vector<double>* pColors         = nullptr;
    const std::vector<double>* pConcentrations = nullptr;
    const std::vector<double>* pPosCharges     = nullptr;
    const std::vector<double>* pNegCharges     = nullptr;
    const std::vector<double>* pPermittivities = nullptr;
    const std::vector<Vector>* pForces         = nullptr;
    const std::vector<Vector>* pFiberDirs      = nullptr;
    const std::vector<Vector>* pVelocities     = nullptr;
    const std::vector<Vector>* pAreas          = nullptr;
    const std::vector<Matrix3>* pSizes         = nullptr;

    // Create pairs for the particle variables created by
    // the SpecialGeometryPieces
    auto vol_pair   = std::make_pair("p.volume", obj.get());
    auto temp_pair  = std::make_pair("p.temperature", obj.get());
    auto color_pair = std::make_pair("p.color", obj.get());
    auto conc_pair  = std::make_pair("p.concentration", obj.get());
    auto pos_pair   = std::make_pair("p.poscharge", obj.get());
    auto neg_pair   = std::make_pair("p.negcharge", obj.get());
    auto perm_pair  = std::make_pair("p.permittivity", obj.get());
    auto force_pair = std::make_pair("p.externalforce", obj.get());
    auto fiber_pair = std::make_pair("p.fiber", obj.get());
    auto vel_pair   = std::make_pair("p.velocity", obj.get());
    auto area_pair  = std::make_pair("p.area", obj.get());
    auto size_pair  = std::make_pair("p.size", obj.get());

    // Set up iterators for SpecialGeometryObject data
    std::vector<double>::const_iterator sgp_vol_iter;
    std::vector<double>::const_iterator sgp_temp_iter;
    std::vector<double>::const_iterator sgp_color_iter;
    std::vector<double>::const_iterator sgp_concentration_iter;
    std::vector<double>::const_iterator sgp_poscharge_iter;
    std::vector<double>::const_iterator sgp_negcharge_iter;
    std::vector<double>::const_iterator sgp_permittivity_iter;
    std::vector<Vector>::const_iterator sgp_force_iter;
    std::vector<Vector>::const_iterator sgp_fiber_iter;
    std::vector<Vector>::const_iterator sgp_velocity_iter;
    std::vector<Vector>::const_iterator sgp_area_iter;
    std::vector<Matrix3>::const_iterator sgp_size_iter;

    if (sgp) {

      std::cout << "Created a special geometry with #particles = "
                << numParticles << std::endl;

      if ((pVolumes = sgp->getScalar("p.volume"))) {
        sgp_vol_iter = obj_vars.scalars.at(vol_pair).begin();
      }

      if ((pTemperatures = sgp->getScalar("p.temperature"))) {
        sgp_temp_iter = obj_vars.scalars.at(temp_pair).begin();
      }

      if ((pForces = sgp->getVector("p.externalforce"))) {
        sgp_force_iter = obj_vars.vectors.at(force_pair).begin();
      }

      if ((pFiberDirs = sgp->getVector("p.fiberdirs"))) {
        sgp_fiber_iter = obj_vars.vectors.at(fiber_pair).begin();
      }

      if ((pVelocities = sgp->getVector("p.velocity"))) {
        sgp_velocity_iter = obj_vars.vectors.at(vel_pair).begin();
      }

      if ((pSizes = sgp->getTensor("p.size"))) {
        sgp_size_iter = obj_vars.tensors.at(size_pair).begin();
      }

      if (d_withColor) {
        if ((pColors = sgp->getScalar("p.color"))) {
          sgp_color_iter = obj_vars.scalars.at(color_pair).begin();
        }
      }

      if (d_doScalarDiffusion) {
        if ((pConcentrations = sgp->getScalar("p.concentration"))) {
          sgp_concentration_iter = obj_vars.scalars.at(conc_pair).begin();
        }

        if ((pAreas = sgp->getVector("p.area"))) {
          sgp_area_iter = obj_vars.vectors.at(area_pair).begin();
        }
      }

      if (d_withGaussSolver) {
        if ((pPosCharges = sgp->getScalar("p.poscharge"))) {
          sgp_poscharge_iter = obj_vars.scalars.at(pos_pair).begin();
        }
        if ((pNegCharges = sgp->getScalar("p.negcharge"))) {
          sgp_negcharge_iter = obj_vars.scalars.at(neg_pair).begin();
        }
        if ((pPermittivities = sgp->getScalar("p.permittivity"))) {
          sgp_permittivity_iter = obj_vars.scalars.at(perm_pair).begin();
        }
      }

    } else {
      std::cout << "Created a geometry with #particles = " << numParticles
                << std::endl;
    }

    for (auto point : obj_vars.points.at(obj.get())) {
      IntVector cell_idx;
      if (!patch->findCell(point, cell_idx)) {
        continue;
      }

      if (!patch->containsPoint(point)) {
        continue;
      }

      particleIndex pidx = start + count;

      // std::cout << "Point["<<pidx<<"]="<<point<<" Cell = "<<cell_idx<<endl;
      initializeParticle(patch,
                         obj.get(),
                         matl,
                         point,
                         cell_idx,
                         pidx,
                         cellNAPID,
                         pvars);

      // Again, everything below exists for SpecialGeometryPiece only
      if (sgp) {

        if (sgp_vol_iter != obj_vars.scalars.at(vol_pair).end()) {
          pvars.pVolume[pidx] = *sgp_vol_iter;
          pvars.pMass[pidx]   = matl->getInitialDensity() * pvars.pVolume[pidx];
          ++sgp_vol_iter;
        }

        if (sgp_temp_iter != obj_vars.scalars.at(temp_pair).end()) {
          pvars.pTemperature[pidx] = *sgp_temp_iter;
          ++sgp_temp_iter;
        }

        if (sgp_force_iter != obj_vars.vectors.at(force_pair).end()) {
          pvars.pExternalForce[pidx] = *sgp_force_iter;
          ++sgp_force_iter;
        }

        if (sgp_velocity_iter != obj_vars.vectors.at(vel_pair).end()) {
          pvars.pVelocity[pidx] = *sgp_velocity_iter;
          ++sgp_velocity_iter;
        }

        if (sgp_fiber_iter != obj_vars.vectors.at(fiber_pair).end()) {
          pvars.pFiberDir[pidx] = *sgp_fiber_iter;
          ++sgp_fiber_iter;
        }

        if (!d_useCPTI) {
          // CPDI and others
          // Read pSize from file
          if (sgp_size_iter != obj_vars.tensors.at(size_pair).end()) {
            Vector dxcc       = patch->dCell();
            pvars.pSize[pidx] = *sgp_size_iter;
            // Calculate CPDI hexahedron volume from pSize
            // (if volume not passed from FileGeometryPiece)
            pvars.pVolume[pidx] = std::abs(pvars.pSize[pidx].Determinant());
            pvars.pMass[pidx] = matl->getInitialDensity() * pvars.pVolume[pidx];

            // Modify pSize (CPDI R-vectors) to be normalized by cell spacing
            Matrix3 size(1. / ((double)dxcc.x()),
                         0.,
                         0.,
                         0.,
                         1. / ((double)dxcc.y()),
                         0.,
                         0.,
                         0.,
                         1. / ((double)dxcc.z()));
            pvars.pSize[pidx] = pvars.pSize[pidx] * size;
            ++sgp_size_iter;
          }
        } else {
          // CPTI
          // Read pSize from file
          if (sgp_size_iter != obj_vars.tensors.at(size_pair).end()) {
            Vector dxcc       = patch->dCell();
            pvars.pSize[pidx] = *sgp_size_iter;
            // Calculate CPTI tetrahedron volume from pSize
            // (if volume not passed from FileGeometryPiece)
            pvars.pVolume[pidx] =
              std::abs(pvars.pSize[pidx].Determinant() / 6.0);
            pvars.pMass[pidx] = matl->getInitialDensity() * pvars.pVolume[pidx];

            // Modify pSize (CPTI R-vectors) to be normalized by cell spacing
            Matrix3 size(1. / ((double)dxcc.x()),
                         0.,
                         0.,
                         0.,
                         1. / ((double)dxcc.y()),
                         0.,
                         0.,
                         0.,
                         1. / ((double)dxcc.z()));
            pvars.pSize[pidx] = pvars.pSize[pidx] * size;
            ++sgp_size_iter;
          }
        }

        if (sgp_color_iter != obj_vars.scalars.at(color_pair).end()) {
          pvars.pColor[pidx] = *sgp_color_iter;
          ++sgp_color_iter;
        }

        if (sgp_area_iter != obj_vars.vectors.at(area_pair).end()) {
          pvars.pArea[pidx] = *sgp_area_iter;
          ++sgp_area_iter;
        }

        if (sgp_concentration_iter != obj_vars.scalars.at(conc_pair).end()) {
          pvars.pConcentration[pidx] = *sgp_concentration_iter;
          ++sgp_concentration_iter;
        }
      }

      // If the particle is on the surface and if there is
      // a physical BC attached to it then mark with the
      // physical BC pointer
      if (d_useLoadCurves) {
        if (checkForSurface(piece, point, dxpp)) {
          pvars.pLoadCurveID[pidx] = getLoadCurveID(point, dxpp);
          // std::cout << " Particle: " << pidx << " use_load_curves = " <<
          // d_useLoadCurves << std::endl; std::cout << "\t surface particle;
          // Load curve id = " << pvars.pLoadCurveID[pidx] << std::endl;
        } else {
          pvars.pLoadCurveID[pidx] = 0;
          // std::cout << "\t not surface particle; Load curve id = " <<
          // pLoadCurveID[pidx] << std::endl;
        }
      }
      count++;
    }
    start += count;
  }
  return numParticles;
}

// Get the LoadCurveID applicable for this material point
// WARNING : Should be called only once per particle during a simulation
// because it updates the number of particles to which a BC is applied.
int
ParticleCreator::getLoadCurveID(const Point& pp, const Vector& dxpp)
{
  int ret = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();

    // std::cerr << " BC Type = " << bcType << std::endl;
    if (bcType == "Pressure") {
      PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
        // std::cout << "\t surface particle; flagged material pt" << std::endl;
        ret = pbc->loadCurveID();
      }
    } else if (bcType == "Velocity") {
      VelocityBC* vbc = dynamic_cast<VelocityBC*>(bc.get());
      if (vbc->flagMaterialPoint(pp, dxpp)) {
        // std::cout << "\t surface particle; flagged material pt" << std::endl;
        ret = vbc->loadCurveID();
      }
    } else if (bcType == "Moment") {
      MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
        ret = pbc->loadCurveID();
      }
    } else if (bcType == "HeatFlux") {
      HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
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
      PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
      std::cerr << *pbc << std::endl;
    }
    if (bcType == "Velocity") {
      VelocityBC* vbc = dynamic_cast<VelocityBC*>(bc.get());
      std::cerr << *vbc << std::endl;
    }
    if (bcType == "Moment") {
      MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
      std::cerr << *pbc << std::endl;
    }
    if (bcType == "HeatFlux") {
      HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
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
      ForceBC* fbc = dynamic_cast<ForceBC*>(bc.get());

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

ParticleSubset*
ParticleCreator::allocateVariables(particleIndex numParticles,
                                   int dwi,
                                   const Patch* patch,
                                   DataWarehouse* new_dw,
                                   ParticleVars& pvars)
{

  ParticleSubset* subset =
    new_dw->createParticleSubset(numParticles, dwi, patch);
  new_dw->allocateAndPut(pvars.position, d_mpm_labels->pXLabel, subset);
  new_dw->allocateAndPut(pvars.pDisp, d_mpm_labels->pDispLabel, subset);
  new_dw->allocateAndPut(pvars.pVelocity, d_mpm_labels->pVelocityLabel, subset);
  new_dw->allocateAndPut(pvars.pAcc, d_mpm_labels->pAccelerationLabel, subset);
  new_dw->allocateAndPut(pvars.pExternalForce,
                         d_mpm_labels->pExternalForceLabel,
                         subset);
  new_dw->allocateAndPut(pvars.pMass, d_mpm_labels->pMassLabel, subset);
  new_dw->allocateAndPut(pvars.pVolume, d_mpm_labels->pVolumeLabel, subset);
  new_dw->allocateAndPut(pvars.pTemperature,
                         d_mpm_labels->pTemperatureLabel,
                         subset);
  new_dw->allocateAndPut(pvars.pParticleID,
                         d_mpm_labels->pParticleIDLabel,
                         subset);
  new_dw->allocateAndPut(pvars.pSize, d_mpm_labels->pSizeLabel, subset);
  new_dw->allocateAndPut(pvars.pFiberDir, d_mpm_labels->pFiberDirLabel, subset);
  // for thermal stress
  new_dw->allocateAndPut(pvars.pTempPrevious,
                         d_mpm_labels->pTempPreviousLabel,
                         subset);

  if (d_useLoadCurves) {
    new_dw->allocateAndPut(pvars.pLoadCurveID,
                           d_mpm_labels->pLoadCurveIDLabel,
                           subset);
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
    new_dw->allocateAndPut(pvars.pLastLevel,
                           d_mpm_labels->pLastLevelLabel,
                           subset);
  }

  // For body force calculation
  new_dw->allocateAndPut(pvars.pBodyForceAcc,
                         d_mpm_labels->pBodyForceAccLabel,
                         subset);
  new_dw->allocateAndPut(pvars.pCoriolisImportance,
                         d_mpm_labels->pCoriolisImportanceLabel,
                         subset);

  // For switching between implicit and explicit
  new_dw->allocateAndPut(pvars.pExternalHeatFlux,
                         d_mpm_labels->pExternalHeatFluxLabel,
                         subset);

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
  Level* fineLevel      = 0;
  if (hasFiner) {
    fineLevel = (Level*)curLevel->getFinerLevel().get_rep();
  }
  for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
    Point lower = patch->nodePosition(*iter) + dcorner;
    IntVector c = *iter;

    if (hasFiner) { // Don't create particles if a finer level exists here
      const Point CC         = patch->cellPosition(c);
      bool includeExtraCells = false;
      const Patch* patchExists =
        fineLevel->getPatchFromPoint(CC, includeExtraCells);
      if (patchExists != 0) {
        continue;
      }
    }

    for (int ix = 0; ix < ppc.x(); ix++) {
      for (int iy = 0; iy < ppc.y(); iy++) {
        for (int iz = 0; iz < ppc.z(); iz++) {

          IntVector idx(ix, iy, iz);
          Point p = lower + dxpp * idx;
          if (!b2.contains(p)) {
            throw InternalError("Particle created outside of patch?",
                                __FILE__,
                                __LINE__);
          }
          if (piece->inside(p)) {
            Vector p1(p(0), p(1), p(2));
            p1   = affineTrans_A * p1 + affineTrans_b;
            p(0) = p1[0];
            p(1) = p1[1];
            p(2) = p1[2];
            obj_vars.points.at(obj).push_back(p);
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

particleIndex
ParticleCreator::countAndCreateParticles(const Patch* patch,
                                         GeometryObject* obj,
                                         ObjectVars& obj_vars)
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

  // If the object is a SpecialGeomPiece (e.g. FileGeometryPiece or
  // SmoothCylGeomPiece) then use the particle creators in that
  // class to do the counting
  SpecialGeomPiece* sgp = dynamic_cast<SpecialGeomPiece*>(piece.get());
  if (sgp) {
    int numPts             = 0;
    FileGeometryPiece* fgp = dynamic_cast<FileGeometryPiece*>(piece.get());
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

vector<const VarLabel*>
ParticleCreator::returnParticleState()
{
  return particle_state;
}

vector<const VarLabel*>
ParticleCreator::returnParticleStatePreReloc()
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

  // For friction contact
  particle_state.push_back(d_mpm_labels->pSurfLabel);
  particle_state_preReloc.push_back(d_mpm_labels->pSurfLabel_preReloc);

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

  if (d_withGaussSolver) {
    particle_state.push_back(d_amrmpm_labels->pPosChargeLabel);
    particle_state_preReloc.push_back(
      d_amrmpm_labels->pPosChargeLabel_preReloc);

    particle_state.push_back(d_amrmpm_labels->pNegChargeLabel);
    particle_state_preReloc.push_back(
      d_amrmpm_labels->pNegChargeLabel_preReloc);

    particle_state.push_back(d_amrmpm_labels->pPosChargeGradLabel);
    particle_state_preReloc.push_back(
      d_amrmpm_labels->pPosChargeGradLabel_preReloc);

    particle_state.push_back(d_amrmpm_labels->pNegChargeGradLabel);
    particle_state_preReloc.push_back(
      d_amrmpm_labels->pNegChargeGradLabel_preReloc);

    particle_state.push_back(d_amrmpm_labels->pPermittivityLabel);
    particle_state_preReloc.push_back(
      d_amrmpm_labels->pPermittivityLabel_preReloc);
  }
}

int
ParticleCreator::checkForSurface(const GeometryPieceP piece,
                                 const Point p,
                                 const Vector dxpp)
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

double
ParticleCreator::checkForSurface2(const GeometryPieceP piece,
                                  const Point p,
                                  const Vector dxpp)
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