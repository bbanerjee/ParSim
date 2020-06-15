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

#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/SmoothGeomPiece.h>
#include <Core/Labels/MPMLabel.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/CrackBC.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <CCA/Components/MPM/MMS/MMS.h>
#include <fstream>
#include <iostream>

using namespace Uintah;
using std::vector;
using std::cerr;
using std::ofstream;

ParticleCreator::ParticleCreator(MPMMaterial* matl, 
                                 MPMFlags* flags)
:d_lock("Particle Creator lock")
{
  d_lb = scinew MPMLabel();
  d_flags = flags;

  d_useLoadCurves = flags->d_useLoadCurves;
  d_withColor = flags->d_withColor;
  d_artificialViscosity = flags->d_artificialViscosity;
  d_computeScaleFactor = flags->d_computeScaleFactor;
  d_doScalarDiffusion = flags->d_doScalarDiffusion;
  d_useCPTI = flags->d_useCPTI;

  registerPermanentParticleState(matl);
}

ParticleCreator::~ParticleCreator()
{
  if (d_lb) {
    delete d_lb;
  }
}

particleIndex
ParticleCreator::createParticles(MPMMaterial* matl,
                                 CCVariable<short int>& cellNAPID,
                                 const Patch* patch,DataWarehouse* new_dw,
                                 vector<GeometryObject*>& d_geom_objs)
{
  ObjectVars vars;
  particleIndex numParticles = 0;
  for (auto geom : d_geom_objs) {
    numParticles += countAndCreateParticles(patch, geom, vars);
  }
  
  ParticleVars pvars;
  int dwi = matl->getDWIndex();
  allocateVariables(numParticles,dwi,patch,new_dw, pvars);

  particleIndex start = 0;
  
  for (auto obj : d_geom_objs) {
    particleIndex count = 0;
    GeometryPieceP piece = obj->getPiece();
    Box b1 = piece->getBoundingBox();
    Box b2 = patch->getExtraBox();
    Box b = b1.intersect(b2);
    if(b.degenerate()) {
      count = 0;
      continue;
    }

    Vector dxpp = patch->dCell()/obj->getInitialData_IntVector("res");    

    // Special case exception for SmoothGeomPieces and FileGeometryPieces
    SmoothGeomPiece *sgp = dynamic_cast<SmoothGeomPiece*>(piece.get_rep());
    vector<double>* volumes       = 0;
    vector<double>* temperatures  = 0;
    vector<double>* colors        = 0;
    vector<Vector>* pforces       = 0;
    vector<Vector>* pFiberDirs    = 0;
    vector<Vector>* pvelocities   = 0;    // gcd adds and new change name
    vector<Matrix3>* pSizes       = 0;
    if (sgp){

      std::cout << "Created a special geometry with #particles = " << numParticles << std::endl;
      volumes      = sgp->getVolume();
      temperatures = sgp->getTemperature();
      pforces      = sgp->getForces();
      pFiberDirs   = sgp->getFiberDirs();
      pvelocities  = sgp->getVelocity();  // gcd adds and new change name
      pSizes       = sgp->getSize();

      if(d_withColor){
        colors      = sgp->getColors();
      }
    } else {
      std::cout << "Created a geometry with #particles = " << numParticles << std::endl;
    }

    // For getting particle volumes (if they exist)
    vector<double>::const_iterator voliter;
    if (volumes) {
      if (!volumes->empty()) voliter = vars.d_object_vols[obj].begin();
    }

    // For getting particle temps (if they exist)
    vector<double>::const_iterator tempiter;
    if (temperatures) {
      if (!temperatures->empty()) tempiter = vars.d_object_temps[obj].begin();
    }

    // For getting particle external forces (if they exist)
    vector<Vector>::const_iterator forceiter;
    if (pforces) {
      if (!pforces->empty()) forceiter = vars.d_object_forces[obj].begin();
    }

    // For getting particle fiber directions (if they exist)
    vector<Vector>::const_iterator fiberiter;
    if (pFiberDirs) {
      if (!pFiberDirs->empty()) fiberiter = vars.d_object_fibers[obj].begin();
    }
    
    // For getting particle velocities (if they exist)   // gcd adds
    vector<Vector>::const_iterator velocityiter;
    if (pvelocities) {                             // new change name
      if (!pvelocities->empty()) velocityiter =
              vars.d_object_velocity[obj].begin();  // new change name
    }                                                    // end gcd adds
    
    // For getting particle sizes (if they exist)
    vector<Matrix3>::const_iterator sizeiter;
    if (pSizes) {
      if (!pSizes->empty()) sizeiter = vars.d_object_size[obj].begin();
      if (d_flags->d_AMR) {
        cerr << "WARNING:  The particle size when using smooth or file\n"; 
        cerr << "geom pieces needs some work when used with AMR" << endl;
      }
    }

    // For getting particles colors (if they exist)
    vector<double>::const_iterator coloriter;
    if (colors) {
      if (!colors->empty()) coloriter = vars.d_object_colors[obj].begin();
    }

    for(auto point : vars.d_object_points[obj]) {
      IntVector cell_idx;
      if (!patch->findCell(point,cell_idx)) continue;

      if (!patch->containsPoint(point)) continue;
      
      particleIndex pidx = start+count;      

      //std::cout << "Point["<<pidx<<"]="<<point<<" Cell = "<<cell_idx<<endl;
 
      initializeParticle(patch, obj,matl,point,cell_idx,pidx,cellNAPID, pvars);

      // Again, everything below exists for FileGeometryPiece only, where
      // a user can describe the geometry as a series of points in a file.
      // One can also describe any of the fields below in the file as well.
      // See FileGeometryPiece for usage.

      if (temperatures) {
        if (!temperatures->empty()) {
          pvars.pTemperature[pidx] = *tempiter;
          ++tempiter;
        }
      }

      if (pforces) {                           
        if (!pforces->empty()) {
          pvars.pExternalForce[pidx] = *forceiter;
          ++forceiter;
        }
      }

      if (pvelocities) {                           // gcd adds and change name 
        if (!pvelocities->empty()) {               // and change name
          pvars.pVelocity[pidx] = *velocityiter;
          ++velocityiter;
        }
      }                                         // end gcd adds

      if (pFiberDirs) {
        if (!pFiberDirs->empty()) {
          pvars.pFiberDir[pidx] = *fiberiter;
          ++fiberiter;
        }
      }
      
      if (volumes) {
        if (!volumes->empty()) {
          pvars.pVolume[pidx] = *voliter;
          pvars.pMass[pidx] = matl->getInitialDensity()*pvars.pVolume[pidx];
          ++voliter;
        }
      }

      // CPDI
      if (pSizes && (d_flags->d_interpolatorType=="cpdi")) {

        // Read pSize from file
        if (!pSizes->empty()) {
          Vector dxcc = patch->dCell(); 
          pvars.pSize[pidx] = *sizeiter;
          if (volumes->empty()) {
            // Calculate CPDI hexahedron volume from pSize 
            // (if volume not passed from FileGeometryPiece)
            pvars.pVolume[pidx] = std::abs(pvars.pSize[pidx].Determinant());
            pvars.pMass[pidx] = matl->getInitialDensity()*pvars.pVolume[pidx];
          }

          // Modify pSize (CPDI R-vectors) to be normalized by cell spacing
          Matrix3 size(1./((double) dxcc.x()),0.,0.,
                       0.,1./((double) dxcc.y()),0.,
                       0.,0.,1./((double) dxcc.z()));
          pvars.pSize[pidx]= pvars.pSize[pidx]*size;
          ++sizeiter;
        }
      }

      // CPTI
      if (pSizes && d_useCPTI) {

        // Read pSize from file
        if (!pSizes->empty()) {
          Vector dxcc = patch->dCell(); 
          pvars.pSize[pidx] = *sizeiter;
          if (volumes->empty()) {
            // Calculate CPTI tetrahedron volume from pSize 
            // (if volume not passed from FileGeometryPiece)
            pvars.pVolume[pidx] = std::abs(pvars.pSize[pidx].Determinant()/6.0);
            pvars.pMass[pidx] = matl->getInitialDensity()*pvars.pVolume[pidx];
          }

          // Modify pSize (CPTI R-vectors) to be normalized by cell spacing
          Matrix3 size(1./((double) dxcc.x()),0.,0.,
                       0.,1./((double) dxcc.y()),0.,
                       0.,0.,1./((double) dxcc.z()));
          pvars.pSize[pidx]= pvars.pSize[pidx]*size;
          ++sizeiter;
        }
      }

      if (colors) {
        if (!colors->empty()) {
          pvars.pColor[pidx] = *coloriter;
          ++coloriter;
        }
      }

      // If the particle is on the surface and if there is
      // a physical BC attached to it then mark with the 
      // physical BC pointer
      if (d_useLoadCurves) {
        if (checkForSurface(piece,point,dxpp)) {
          pvars.pLoadCurveID[pidx] = getLoadCurveID(point, dxpp);
          //std::cout << " Particle: " << pidx << " use_load_curves = " << d_useLoadCurves << std::endl;
          //std::cout << "\t surface particle; Load curve id = " << pvars.pLoadCurveID[pidx] << std::endl;
        } else {
          pvars.pLoadCurveID[pidx] = 0;
          //std::cout << "\t not surface particle; Load curve id = " << pLoadCurveID[pidx] << std::endl;
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
int ParticleCreator::getLoadCurveID(const Point& pp, const Vector& dxpp)
{
  int ret=0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();
        
    //cerr << " BC Type = " << bcType << endl;
    if (bcType == "Pressure") {
      PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
         //std::cout << "\t surface particle; flagged material pt" << std::endl;
         ret = pbc->loadCurveID(); 
      } 
    } else if (bcType == "Velocity") {
      VelocityBC* vbc = dynamic_cast<VelocityBC*>(bc.get());
      if (vbc->flagMaterialPoint(pp, dxpp)) {
         //std::cout << "\t surface particle; flagged material pt" << std::endl;
         ret = vbc->loadCurveID(); 
      } 
    } else if (bcType == "Moment") {
      MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
      if (pbc->flagMaterialPoint(pp, dxpp)) {
         ret = pbc->loadCurveID(); 
      }
    }
    else if (bcType == "HeatFlux") {      
      HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
      if (hfbc->flagMaterialPoint(pp, dxpp)) {
        ret = hfbc->loadCurveID(); 
      }
    }
  }
  return ret;
}

// Print MPM physical boundary condition information
void ParticleCreator::printPhysicalBCs()
{
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    string bcType = bc->getType();
    if (bcType == "Pressure") {
      PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
      cerr << *pbc << endl;
    }
    if (bcType == "Velocity") {
      VelocityBC* vbc = dynamic_cast<VelocityBC*>(bc.get());
      cerr << *vbc << endl;
    }
    if (bcType == "Moment") {
      MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
      cerr << *pbc << endl;
    }
    if (bcType == "HeatFlux") {
      HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(bc.get());
      cerr << *hfbc << endl;
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
        
    //cerr << " BC Type = " << bcType << endl;
    if (bcType == "Force") {
      ForceBC* fbc = dynamic_cast<ForceBC*>(bc.get());

      Box fbcBox(fbc->getLowerRange()-dxpp,fbc->getUpperRange()+dxpp);

      //cerr << "BC Box = " << bcBox << " Point = " << pp << endl;
      if(fbcBox.contains(pp)) {
        pExtForce = fbc->getForceDensity() * pMass;
        //cerr << "External Force on Particle = " << pExtForce 
        //     << " Force Density = " << bc->getForceDensity() 
        //     << " Particle Mass = " << pMass << endl;
      }
    } 
  }
}

ParticleSubset* 
ParticleCreator::allocateVariables(particleIndex numParticles, 
                                   int dwi, const Patch* patch,
                                   DataWarehouse* new_dw,
                                   ParticleVars& pvars)
{

  ParticleSubset* subset = new_dw->createParticleSubset(numParticles,dwi,
                                                        patch);
  new_dw->allocateAndPut(pvars.position,       d_lb->pXLabel,             subset);
  new_dw->allocateAndPut(pvars.pDisp,          d_lb->pDispLabel,          subset);
  new_dw->allocateAndPut(pvars.pVelocity,      d_lb->pVelocityLabel,      subset); 
  new_dw->allocateAndPut(pvars.pAcc,           d_lb->pAccelerationLabel,  subset);
  new_dw->allocateAndPut(pvars.pExternalForce, d_lb->pExternalForceLabel, subset);
  new_dw->allocateAndPut(pvars.pMass,          d_lb->pMassLabel,          subset);
  new_dw->allocateAndPut(pvars.pVolume,        d_lb->pVolumeLabel,        subset);
  new_dw->allocateAndPut(pvars.pTemperature,   d_lb->pTemperatureLabel,   subset);
  new_dw->allocateAndPut(pvars.pParticleID,    d_lb->pParticleIDLabel,    subset);
  new_dw->allocateAndPut(pvars.pSize,          d_lb->pSizeLabel,          subset);
  new_dw->allocateAndPut(pvars.pFiberDir,      d_lb->pFiberDirLabel,      subset); 
  // for thermal stress
  new_dw->allocateAndPut(pvars.pTempPrevious,  d_lb->pTempPreviousLabel,  subset); 
  
  if (d_useLoadCurves) {
    new_dw->allocateAndPut(pvars.pLoadCurveID, d_lb->pLoadCurveIDLabel,   subset); 
  }
  if(d_withColor){
     new_dw->allocateAndPut(pvars.pColor,      d_lb->pColorLabel,         subset);
  }
  if(d_artificialViscosity){
     new_dw->allocateAndPut(pvars.p_q,         d_lb->p_qLabel,            subset);
  }

  // For AMR
  new_dw->allocateAndPut(pvars.pRefined,       d_lb->pRefinedLabel,       subset);
  if (d_flags->d_AMR) {
     new_dw->allocateAndPut(pvars.pLastLevel,  d_lb->pLastLevelLabel,     subset);
  }

  // For body force calculation
  new_dw->allocateAndPut(pvars.pBodyForceAcc,       d_lb->pBodyForceAccLabel,       subset);
  new_dw->allocateAndPut(pvars.pCoriolisImportance, d_lb->pCoriolisImportanceLabel, subset);

  // For switching between implicit and explicit
  new_dw->allocateAndPut(pvars.pExternalHeatFlux, d_lb->pExternalHeatFluxLabel, subset);

  // For friction contact
  new_dw->allocateAndPut(pvars.pSurface,       d_lb->pSurfLabel,         subset);

  return subset;
}


void ParticleCreator::createPoints(const Patch* patch, GeometryObject* obj,
                                   ObjectVars& vars)
{

  GeometryPieceP piece = obj->getPiece();
  Box b2 = patch->getExtraBox();
  IntVector ppc = obj->getInitialData_IntVector("res");
  Vector dxpp = patch->dCell()/ppc;
  Vector dcorner = dxpp*0.5;

  // Affine transformation for making conforming particle distributions
  // to be used in the conforming CPDI simulations. The input vectors are
  // optional and if you do not wish to use the affine transformation, just do
  // not define them in the input file.
  Vector affineTrans_A0=obj->getInitialData_Vector("affineTransformation_A0");
  Vector affineTrans_A1=obj->getInitialData_Vector("affineTransformation_A1");
  Vector affineTrans_A2=obj->getInitialData_Vector("affineTransformation_A2");
  Vector affineTrans_b= obj->getInitialData_Vector("affineTransformation_b");
  Matrix3 affineTrans_A(
          affineTrans_A0[0],affineTrans_A0[1],affineTrans_A0[2],
          affineTrans_A1[0],affineTrans_A1[1],affineTrans_A1[2],
          affineTrans_A2[0],affineTrans_A2[1],affineTrans_A2[2]);

  // AMR stuff
  const Level* curLevel = patch->getLevel();
  bool hasFiner = curLevel->hasFinerLevel();
  Level* fineLevel=0;
  if(hasFiner){
    fineLevel = (Level*) curLevel->getFinerLevel().get_rep();
  }
  for(CellIterator iter = patch->getCellIterator(); !iter.done(); iter++){
    Point lower = patch->nodePosition(*iter) + dcorner;
    IntVector c = *iter;
    
    if(hasFiner){ // Don't create particles if a finer level exists here
      const Point CC = patch->cellPosition(c);
      bool includeExtraCells= false;
      const Patch* patchExists = fineLevel->getPatchFromPoint(CC,includeExtraCells);
      if(patchExists != 0){
       continue;
      }
    }

    for(int ix=0;ix < ppc.x(); ix++){
      for(int iy=0;iy < ppc.y(); iy++){
        for(int iz=0;iz < ppc.z(); iz++){
        
          IntVector idx(ix, iy, iz);
          Point p = lower + dxpp*idx;
          if (!b2.contains(p)){
            throw InternalError("Particle created outside of patch?", __FILE__, __LINE__);
          }
          if (piece->inside(p)){ 
            Vector p1(p(0),p(1),p(2));
            p1=affineTrans_A*p1+affineTrans_b;
            p(0)=p1[0];
            p(1)=p1[1];
            p(2)=p1[2];
            vars.d_object_points[obj].push_back(p);
          }
        }  // z
      }  // y
    }  // x
  }  // iterator

/*
//  *** This part is associated with CBDI_CompressiveCylinder.ups input file.
//      It creates conforming particle distribution to be used in the simulation.
//      To use that you need to uncomment the following commands to create the
//      conforming particle distribution and comment above commands that are used
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
          d_object_points[key].push_back(p);
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
  Vector dxpp = patch->dCell()/obj->getInitialData_IntVector("res");
  Vector dxcc = patch->dCell();

  // Affine transformation for making conforming particle distributions
  //  to be used in the conforming CPDI simulations. The input vectors are
  //  optional and if you do not want to use affine transformations, just do
  //  not define them in the input file.

  Vector affineTrans_A0=obj->getInitialData_Vector("affineTransformation_A0");
  Vector affineTrans_A1=obj->getInitialData_Vector("affineTransformation_A1");
  Vector affineTrans_A2=obj->getInitialData_Vector("affineTransformation_A2");
  //Vector affineTrans_b= obj->getInitialData_Vector("affineTransformation_b");
  Matrix3 affineTrans_A(
          affineTrans_A0[0],affineTrans_A0[1],affineTrans_A0[2],
          affineTrans_A1[0],affineTrans_A1[1],affineTrans_A1[2],
          affineTrans_A2[0],affineTrans_A2[1],affineTrans_A2[2]);
  // The size matrix is used for storing particle domain sizes 
  // (Rvectors for CPDI and CPTI) normalized by the grid spacing
  Matrix3 size(1./((double) ppc.x()),0.,0.,
               0.,1./((double) ppc.y()),0.,
               0.,0.,1./((double) ppc.z()));
  size=affineTrans_A*size;

/*
//  *** This part is associated with CBDI_CompressiveCylinder.ups input file.
//      It determines particle domain sizes for the conforming particle distribution,
//      which is used in the simulation.
//      To activate that you need to uncomment the following commands to determine
//      particle domain sizes for the conforming particle distribution, and
//      comment above commands that are used to determine particle domain sizes for
//      non-conforming particle distributions.

  int resolutionPart=1;
  int nPar1=180*resolutionPart;
  int nPar2=16*resolutionPart;
  double pi=3.14159265358979;
  double rtemp=sqrt(p.x()*p.x()+p.y()*p.y());
  Matrix3 size(2.*pi/nPar1*p.y()/dxcc[0],2.*0.25/nPar2/rtemp*p.x()/dxcc[1],0.,
              -2.*pi/nPar1*p.x()/dxcc[0],2.*0.25/nPar2/rtemp*p.y()/dxcc[1],0.,
                                      0.,                               0.,1.);
*/

  pvars.pTemperature[i] = obj->getInitialData_double("temperature");


  // For AMR
  const Level* curLevel = patch->getLevel();
  pvars.pRefined[i]     = curLevel->getIndex();

  //MMS
  //std::cout << "mms_type = " << d_flags->d_mmsType << "\n";
  string mms_type = d_flags->d_mmsType;
  if(mms_type != "none") {
    MMS MMSObject;
    MMSObject.initializeParticleForMMS(pvars.position, pvars.pVelocity, pvars.pSize,
                                       pvars.pDisp, pvars.pMass, pvars.pVolume,
                                       p, dxcc, size, patch, d_flags, i);
  }  else {
    pvars.position[i] = p;
     if(d_flags->d_axisymmetric){
      // assume unit radian extent in the circumferential direction
      pvars.pVolume[i] = p.x()*
                   (size(0,0)*size(1,1)-size(0,1)*size(1,0))*dxcc.x()*dxcc.y();
      //std::cout << "dx_cc = " << dxcc << " size = " << size 
      //          << " x = " << p.x() << " vol = " << pvars.pVolume[i] << "\n";
     } else {
     // standard voxel volume
     pvars.pVolume[i]  = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
    }

    pvars.pSize[i]      = size;
    pvars.pDisp[i]      = Vector(0.,0.,0.);
    pvars.pVelocity[i]  = obj->getInitialData_Vector("velocity");
    pvars.pAcc[i]       = Vector(0.,0.,0.);

    //std::cout << "i = " << i << " px = " << pvars.position[i]
    //          << " size = " << pvars.pSize[i] << "\n";
    double vol_frac_CC = 1.0;
    try {
     if(obj->getInitialData_double("volumeFraction") == -1.0) {    
      vol_frac_CC = 1.0;
      pvars.pMass[i]      = matl->getInitialDensity()*pvars.pVolume[i];
     } else {
      vol_frac_CC = obj->getInitialData_double("volumeFraction");
      pvars.pMass[i]      = matl->getInitialDensity()*pvars.pVolume[i]*vol_frac_CC;
     }
    } catch (...) {
      vol_frac_CC = 1.0;       
      pvars.pMass[i]      = matl->getInitialDensity()*pvars.pVolume[i];
    }
  }
  
  if(d_withColor){
    pvars.pColor[i] = obj->getInitialData_double("color");
  }
  if(d_artificialViscosity){
    pvars.p_q[i] = 0.;
  }
  if(d_flags->d_AMR){
    pvars.pLastLevel[i] = curLevel->getID();
  }
  
  pvars.pTempPrevious[i]  = pvars.pTemperature[i];

  Vector pExtForce(0,0,0);
  applyForceBC(dxpp, p, pvars.pMass[i], pExtForce);
  
  pvars.pExternalForce[i] = pExtForce;
  pvars.pFiberDir[i]      = matl->getConstitutiveModel()->getInitialFiberDir();

  // AMR
  if (d_flags->d_AMR) {
    pvars.pLastLevel[i] = curLevel->getID();
  }

  // For body forces
  pvars.pBodyForceAcc[i] = Vector(0.0, 0.0, 0.0); // Init to zero
  pvars.pCoriolisImportance[i] = 0.0;

  // For switching between explicit and implicit MPM
  pvars.pExternalHeatFlux[i] = 0.0;

  // For friction contact
  GeometryPieceP piece = obj->getPiece();
  pvars.pSurface[i] = checkForSurface2(piece, p, dxpp);

  // Cell ids
  ASSERT(cell_idx.x() <= 0xffff && 
         cell_idx.y() <= 0xffff && 
         cell_idx.z() <= 0xffff);
         
  long64 cellID = ((long64)cell_idx.x() << 16) | 
                  ((long64)cell_idx.y() << 32) | 
                  ((long64)cell_idx.z() << 48);
                  
  short int& myCellNAPID = cellNAPID[cell_idx];
  pvars.pParticleID[i] = (cellID | (long64) myCellNAPID);
  ASSERT(myCellNAPID < 0x7fff);
  myCellNAPID++;
}

particleIndex 
ParticleCreator::countAndCreateParticles(const Patch* patch, 
                                         GeometryObject* obj,
                                         ObjectVars& vars)
{
  GeometryPieceP piece = obj->getPiece();
  Box b1 = piece->getBoundingBox();
  Box b2 = patch->getExtraBox();
  Box b = b1.intersect(b2);
  if(b.degenerate()) return 0;
  
  // If the object is a SmoothGeomPiece (e.g. FileGeometryPiece or
  // SmoothCylGeomPiece) then use the particle creators in that 
  // class to do the counting d
  SmoothGeomPiece   *sgp = dynamic_cast<SmoothGeomPiece*>(piece.get_rep());
  if (sgp) {
    int numPts = 0;
    FileGeometryPiece *fgp = dynamic_cast<FileGeometryPiece*>(piece.get_rep());
    if(fgp){
      if (d_useCPTI) {
        proc0cout << "*** Reading CPTI file ***" << endl;
      }
      fgp->readPoints(patch->getID());
      numPts = fgp->returnPointCount();
      proc0cout << "Number of points read from file = " << numPts << std::endl;
    } else {
      Vector dxpp = patch->dCell()/obj->getInitialData_IntVector("res");    
      double dx   = Min(Min(dxpp.x(),dxpp.y()), dxpp.z());
      sgp->setParticleSpacing(dx);
      sgp->setCellSize(patch->dCell());
      numPts = sgp->createPoints();
      proc0cout << "Special Geom Piece: Number of points created = " << numPts << std::endl;
    }
    vector<Point>* points      = sgp->getPoints();
    vector<double>* vols       = sgp->getVolume();
    vector<double>* temps      = sgp->getTemperature();
    vector<double>* colors     = sgp->getColors();
    vector<Vector>* pforces    = sgp->getForces();
    vector<Vector>* pFiberDirs = sgp->getFiberDirs();
    vector<Vector>* pvelocities= sgp->getVelocity();
    vector<Matrix3>* pSizes    = sgp->getSize();
    Point p;
    IntVector cell_idx;
    
    for (int ii = 0; ii < numPts; ++ii) {
      p = points->at(ii);
      if (patch->findCell(p,cell_idx)) {
        if (patch->containsPoint(p)) {
          vars.d_object_points[obj].push_back(p);
          
          if (!vols->empty()) {
            double vol = vols->at(ii); 
            vars.d_object_vols[obj].push_back(vol);
          }
          if (!temps->empty()) {
            double temp = temps->at(ii); 
            vars.d_object_temps[obj].push_back(temp);
          }
          if (!pforces->empty()) {
            Vector pforce = pforces->at(ii); 
            vars.d_object_forces[obj].push_back(pforce);
          }
          if (!pFiberDirs->empty()) {
            Vector pfiber = pFiberDirs->at(ii); 
            vars.d_object_fibers[obj].push_back(pfiber);
          }
          if (!pvelocities->empty()) {
            Vector pvel = pvelocities->at(ii); 
            vars.d_object_velocity[obj].push_back(pvel);
          }
          if (!pSizes->empty()) {
            Matrix3 psz = pSizes->at(ii); 
            vars.d_object_size[obj].push_back(psz);
          }
          if (!colors->empty()) {
            double color = colors->at(ii); 
            vars.d_object_colors[obj].push_back(color);
          }
        } 
      }  // patch contains cell
    }
    //sgp->deletePoints();
    //sgp->deleteVolume();
  } else {
    createPoints(patch, obj, vars);
  }
  
  return (particleIndex) vars.d_object_points[obj].size();
}

vector<const VarLabel* > ParticleCreator::returnParticleState()
{
  return particle_state;
}


vector<const VarLabel* > ParticleCreator::returnParticleStatePreReloc()
{
  return particle_state_preReloc;
}

void ParticleCreator::registerPermanentParticleState(MPMMaterial* matl)
{
  particle_state.push_back(d_lb->pDispLabel);
  particle_state_preReloc.push_back(d_lb->pDispLabel_preReloc);

  particle_state.push_back(d_lb->pVelocityLabel);
  particle_state_preReloc.push_back(d_lb->pVelocityLabel_preReloc);

  particle_state.push_back(d_lb->pAccelerationLabel);
  particle_state_preReloc.push_back(d_lb->pAccelerationLabel_preReloc);

  particle_state.push_back(d_lb->pExternalForceLabel);
  particle_state_preReloc.push_back(d_lb->pExtForceLabel_preReloc);

  particle_state.push_back(d_lb->pMassLabel);
  particle_state_preReloc.push_back(d_lb->pMassLabel_preReloc);

  particle_state.push_back(d_lb->pVolumeLabel);
  particle_state_preReloc.push_back(d_lb->pVolumeLabel_preReloc);

  particle_state.push_back(d_lb->pTemperatureLabel);
  particle_state_preReloc.push_back(d_lb->pTemperatureLabel_preReloc);

  // for thermal stress
  particle_state.push_back(d_lb->pTempPreviousLabel);
  particle_state_preReloc.push_back(d_lb->pTempPreviousLabel_preReloc);
  
  particle_state.push_back(d_lb->pParticleIDLabel);
  particle_state_preReloc.push_back(d_lb->pParticleIDLabel_preReloc);
  

  if (d_withColor){
    particle_state.push_back(d_lb->pColorLabel);
    particle_state_preReloc.push_back(d_lb->pColorLabel_preReloc);
  }

  particle_state.push_back(d_lb->pSizeLabel);
  particle_state_preReloc.push_back(d_lb->pSizeLabel_preReloc);

  if (d_useLoadCurves) {
    particle_state.push_back(d_lb->pLoadCurveIDLabel);
    particle_state_preReloc.push_back(d_lb->pLoadCurveIDLabel_preReloc);
  }

  particle_state.push_back(d_lb->pDispGradLabel);
  particle_state_preReloc.push_back(d_lb->pDispGradLabel_preReloc);

  particle_state.push_back(d_lb->pDefGradLabel);
  particle_state_preReloc.push_back(d_lb->pDefGradLabel_preReloc);

  particle_state.push_back(d_lb->pStressLabel);
  particle_state_preReloc.push_back(d_lb->pStressLabel_preReloc);

  if (d_artificialViscosity) {
    particle_state.push_back(d_lb->p_qLabel);
    particle_state_preReloc.push_back(d_lb->p_qLabel_preReloc);
  }

  if (d_computeScaleFactor) {
    particle_state.push_back(d_lb->pScaleFactorLabel);
    particle_state_preReloc.push_back(d_lb->pScaleFactorLabel_preReloc);
  }

  matl->getConstitutiveModel()->addParticleState(particle_state,
                                                 particle_state_preReloc);

  if (matl->d_doBasicDamage) {
    matl->getBasicDamageModel()->addParticleState(particle_state,
                                                  particle_state_preReloc);
  }

  // For AMR
  if (d_flags->d_refineParticles) {
    particle_state.push_back(d_lb->pRefinedLabel);
    particle_state_preReloc.push_back(d_lb->pRefinedLabel_preReloc);
  }

  if (d_flags->d_AMR) {
    particle_state.push_back(d_lb->pLastLevelLabel);
    particle_state_preReloc.push_back(d_lb->pLastLevelLabel_preReloc);
  }

  // For body forces
  particle_state.push_back(d_lb->pBodyForceAccLabel);
  particle_state_preReloc.push_back(d_lb->pBodyForceAccLabel_preReloc);

  particle_state.push_back(d_lb->pCoriolisImportanceLabel);
  particle_state_preReloc.push_back(d_lb->pCoriolisImportanceLabel_preReloc);

  // For switching between implicit and explicit MPM
  particle_state.push_back(d_lb->pExternalHeatFluxLabel);
  particle_state.push_back(d_lb->pVelGradLabel);

  particle_state_preReloc.push_back(d_lb->pExternalHeatFluxLabel_preReloc);
  particle_state_preReloc.push_back(d_lb->pVelGradLabel_preReloc);

  // For friction contact
  particle_state.push_back(d_lb->pSurfLabel);
  particle_state_preReloc.push_back(d_lb->pSurfLabel_preReloc);
}

int
ParticleCreator::checkForSurface( const GeometryPieceP piece, const Point p,
                                  const Vector dxpp )
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface
  //std::cout << "GeometryPiece = " << piece->getType() << std::endl;
  //std::cout << " Point = " << p << " box = " << dxpp << std::endl;
  
  int ss = 0;
  // Check to the left (-x)
  if(!piece->inside(p-Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check to the right (+x)
  if(!piece->inside(p+Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check behind (-y)
  if(!piece->inside(p-Vector(0.,dxpp.y(),0.)))
    ss++;
  // Check in front (+y)
  if(!piece->inside(p+Vector(0.,dxpp.y(),0.)))
    ss++;
#if 1
  // Check below (-z)
  if(!piece->inside(p-Vector(0.,0.,dxpp.z())))
    ss++;
  // Check above (+z)
  if(!piece->inside(p+Vector(0.,0.,dxpp.z())))
    ss++;
#endif

  if(ss>0){
    return 1;
  }
  else {
    return 0;
  }
}

double
ParticleCreator::checkForSurface2(const GeometryPieceP piece, const Point p,
                                  const Vector dxpp )
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface

  int ss = 0;
  // Check to the left (-x)
  if(!piece->inside(p-Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check to the right (+x)
  if(!piece->inside(p+Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check behind (-y)
  if(!piece->inside(p-Vector(0.,dxpp.y(),0.)))
    ss++;
  // Check in front (+y)
  if(!piece->inside(p+Vector(0.,dxpp.y(),0.)))
    ss++;
  //if (d_flags->d_ndim==3) {
    // Check below (-z)
    if(!piece->inside(p-Vector(0.,0.,dxpp.z())))
      ss++;
    // Check above (+z)
    if(!piece->inside(p+Vector(0.,0.,dxpp.z())))
      ss++;
  //}

  if(ss>0){
    return 1.0;
  } else {
    return 0.0;
  }
}

/*! For particle material change operation */
void ParticleCreator::allocateVariablesAddRequires(Task* task, 
                                                   const MPMMaterial* ,
                                                   const PatchSet* ) const
{
  d_lock.writeLock();
  Ghost::GhostType  gn = Ghost::None;
  task->requires(Task::OldDW,d_lb->pParticleIDLabel,        gn);
  task->requires(Task::OldDW,d_lb->pXLabel,                 gn);
  task->requires(Task::OldDW,d_lb->pMassLabel,              gn);
  //task->requires(Task::OldDW,d_lb->pVolumeLabel,          gn);
  task->requires(Task::OldDW,d_lb->pTemperatureLabel,       gn);
  task->requires(Task::OldDW,d_lb->pDispLabel,              gn);
  task->requires(Task::OldDW,d_lb->pVelocityLabel,          gn);
  task->requires(Task::OldDW,d_lb->pAccelerationLabel,      gn);
  task->requires(Task::NewDW,d_lb->pExtForceLabel_preReloc, gn);
  //task->requires(Task::OldDW,d_lb->pExternalForceLabel,   gn);
  task->requires(Task::NewDW,d_lb->pVolumeLabel_preReloc,   gn);
  task->requires(Task::OldDW,d_lb->pSizeLabel,              gn);
  // for thermal stress
  task->requires(Task::OldDW,d_lb->pTempPreviousLabel, gn); 

  if (d_useLoadCurves){
    task->requires(Task::OldDW,d_lb->pLoadCurveIDLabel, gn);
  }
  if (d_withColor){
    task->requires(Task::OldDW,d_lb->pColorLabel,       gn);
  }
  if(d_artificialViscosity){
    task->requires(Task::OldDW,d_lb->p_qLabel,          gn);
  }

  // For body forces
  task->requires(Task::OldDW, d_lb->pBodyForceAccLabel, gn);
  task->requires(Task::OldDW, d_lb->pCoriolisImportanceLabel, gn);

  // For switching between explicit and implicit MPM
  task->requires(Task::OldDW, d_lb->pExternalHeatFluxLabel, gn);

  d_lock.writeUnlock();
}


/*! For particle material change operation */
void ParticleCreator::allocateVariablesAdd(DataWarehouse* new_dw,
                                           ParticleSubset* addset,
                                           ParticleLabelVariableMap* newState,
                                           ParticleSubset* delset,
                                           DataWarehouse* old_dw)
{
  d_lock.writeLock();

  ParticleVars pvars;
  ParticleSubset::iterator n,o;

  constParticleVariable<Point>  o_position;
  constParticleVariable<Vector> o_disp, o_velocity, o_acc;
  constParticleVariable<Vector> o_external_force;
  constParticleVariable<double> o_mass;
  constParticleVariable<double> o_volume;
  constParticleVariable<double> o_temperature;
  constParticleVariable<double> o_sp_vol;
  constParticleVariable<long64> o_particleID;
  constParticleVariable<Matrix3> o_size;
  constParticleVariable<int>    o_loadcurve;
  constParticleVariable<double> o_tempPrevious; // for thermal stress
  constParticleVariable<double> o_color;
  constParticleVariable<double> o_q;

  constParticleVariable<Vector> o_BodyForceAcc;
  constParticleVariable<double> o_CoriolisImportance;

  constParticleVariable<double>  o_ExternalHeatFlux;
  constParticleVariable<Matrix3> o_VelGrad;

  constParticleVariable<Matrix3> o_DispGrad;
  constParticleVariable<double>  o_dTdt;
  
  new_dw->allocateTemporary(pvars.position,       addset);
  new_dw->allocateTemporary(pvars.pDisp,          addset);
  new_dw->allocateTemporary(pvars.pVelocity,      addset); 
  new_dw->allocateTemporary(pvars.pAcc,           addset);
  new_dw->allocateTemporary(pvars.pExternalForce, addset);
  new_dw->allocateTemporary(pvars.pMass,          addset);
  new_dw->allocateTemporary(pvars.pVolume,        addset);
  new_dw->allocateTemporary(pvars.pTemperature,   addset);
  new_dw->allocateTemporary(pvars.pParticleID,    addset);
  new_dw->allocateTemporary(pvars.pSize,          addset);
  new_dw->allocateTemporary(pvars.pLoadCurveID,   addset); 
  new_dw->allocateTemporary(pvars.pTempPrevious,  addset);

  new_dw->allocateTemporary(pvars.pBodyForceAcc,        addset);
  new_dw->allocateTemporary(pvars.pCoriolisImportance,  addset);

  new_dw->allocateTemporary(pvars.pExternalHeatFlux,    addset);

  old_dw->get(o_position,       d_lb->pXLabel,                delset);
  old_dw->get(o_mass,           d_lb->pMassLabel,             delset);
  old_dw->get(o_particleID,     d_lb->pParticleIDLabel,       delset);
  old_dw->get(o_temperature,    d_lb->pTemperatureLabel,      delset);
  old_dw->get(o_disp,           d_lb->pDispLabel,             delset);
  old_dw->get(o_velocity,       d_lb->pVelocityLabel,         delset);
  old_dw->get(o_acc,            d_lb->pAccelerationLabel,     delset);
  new_dw->get(o_external_force, d_lb->pExtForceLabel_preReloc,delset);
  //old_dw->get(o_external_force,d_lb->pExternalForceLabel,   delset);
  new_dw->get(o_volume,         d_lb->pVolumeLabel_preReloc,  delset);
  //old_dw->get(o_volume,       d_lb->pVolumeLabel,           delset);
  old_dw->get(o_size,           d_lb->pSizeLabel,             delset);
  old_dw->get(o_tempPrevious,   d_lb->pTempPreviousLabel,     delset);
  
  if (d_useLoadCurves){ 
    old_dw->get(o_loadcurve,    d_lb->pLoadCurveIDLabel,      delset);
  }
  if(d_withColor){
    new_dw->allocateTemporary(pvars.pColor,         addset); 
    old_dw->get(o_color,        d_lb->pColorLabel,            delset);
  }
  if(d_artificialViscosity){
    new_dw->allocateTemporary(pvars.p_q,         addset); 
    old_dw->get(o_q,        d_lb->p_qLabel,            delset);
  }
   
  // For body force 
  old_dw->get(o_BodyForceAcc, d_lb->pBodyForceAccLabel, delset);
  old_dw->get(o_CoriolisImportance, d_lb->pCoriolisImportanceLabel, delset);
  new_dw->allocateTemporary(pvars.pBodyForceAcc, addset);
  new_dw->allocateTemporary(pvars.pCoriolisImportance, addset);

  // For switching between explicit and implicit MPM
  old_dw->get(o_ExternalHeatFlux, d_lb->pExternalHeatFluxLabel, delset);
  new_dw->allocateTemporary(pvars.pExternalHeatFlux, addset);

  n = addset->begin();
  for (o=delset->begin(); o != delset->end(); o++, n++) {
    pvars.position[*n]      = o_position[*o];
    pvars.pDisp[*n]         = o_disp[*o];
    pvars.pVelocity[*n]     = o_velocity[*o];
    pvars.pAcc[*n]          = o_acc[*o];
    pvars.pExternalForce[*n]= o_external_force[*o];
    pvars.pMass[*n]         = o_mass[*o];
    pvars.pVolume[*n]       = o_volume[*o];
    pvars.pTemperature[*n]  = o_temperature[*o];
    pvars.pParticleID[*n]   = o_particleID[*o];
    pvars.pSize[*n]         = o_size[*o];
    pvars.pTempPrevious[*n] = o_tempPrevious[*o];  // for thermal stress
    if (d_useLoadCurves){ 
      pvars.pLoadCurveID[*n]= o_loadcurve[*o];
    }
    if (d_withColor){
      pvars.pColor[*n]      = o_color[*o];
    }
    if(d_artificialViscosity){
      pvars.p_q[*n]      = o_q[*o];
    }

    // For body force 
    pvars.pBodyForceAcc[*n] = o_BodyForceAcc[*o];
    pvars.pCoriolisImportance[*n] = o_CoriolisImportance[*o];

    // For switching between explicit and implicit MPM
    pvars.pExternalHeatFlux[*n] = o_ExternalHeatFlux[*o];
  }

  (*newState)[d_lb->pXLabel]              = pvars.position.clone();
  (*newState)[d_lb->pDispLabel]           = pvars.pDisp.clone();
  (*newState)[d_lb->pVelocityLabel]       = pvars.pVelocity.clone();
  (*newState)[d_lb->pAccelerationLabel]   = pvars.pAcc.clone();
  (*newState)[d_lb->pExternalForceLabel]  = pvars.pExternalForce.clone();
  (*newState)[d_lb->pMassLabel]           = pvars.pMass.clone();
  (*newState)[d_lb->pVolumeLabel]         = pvars.pVolume.clone();
  (*newState)[d_lb->pTemperatureLabel]    = pvars.pTemperature.clone();
  (*newState)[d_lb->pParticleIDLabel]     = pvars.pParticleID.clone();
  (*newState)[d_lb->pSizeLabel]           = pvars.pSize.clone();
  (*newState)[d_lb->pTempPreviousLabel]   = pvars.pTempPrevious.clone(); // for thermal stress
  
  if (d_useLoadCurves){ 
    (*newState)[d_lb->pLoadCurveIDLabel]= pvars.pLoadCurveID.clone();
  }
  if(d_withColor){
    (*newState)[d_lb->pColorLabel]      = pvars.pColor.clone();
  }
  if(d_artificialViscosity){
    (*newState)[d_lb->p_qLabel]         = pvars.p_q.clone();
  }

  // For body force
  (*newState)[d_lb->pBodyForceAccLabel] = pvars.pBodyForceAcc.clone();
  (*newState)[d_lb->pCoriolisImportanceLabel] = pvars.pCoriolisImportance.clone();

  // For switching between explicit and implicit MPM
  (*newState)[d_lb->pExternalHeatFluxLabel] = pvars.pExternalHeatFlux.clone();

  d_lock.writeUnlock();
}

