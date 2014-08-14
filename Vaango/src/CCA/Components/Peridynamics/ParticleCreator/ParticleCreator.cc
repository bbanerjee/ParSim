#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleLoadBCFactory.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticlePressureBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleForceBC.h>
#include <CCA/Components/Peridynamics/ParticleBC/ParticleNormalForceBC.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/SmoothGeomPiece.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

using namespace Vaango;

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "PDPartDoing:+,PDPartDebug:+".....
//  bash     : export SCI_DEBUG="PDPartDoing:+,PDPartDebug:+" )
//  default is OFF
using Uintah::DebugStream;
static DebugStream cout_doing("PDPartDoing", false);
static DebugStream cout_dbg("PDPartDebug", false);


ParticleCreator::ParticleCreator(PeridynamicsMaterial* matl, 
                                 PeridynamicsFlags* flags)
  : d_lock("Particle Creator lock")
{
  d_varLabel = scinew PeridynamicsLabel();
  d_flags = flags;
  registerPermanentParticleState(matl);
}

ParticleCreator::~ParticleCreator()
{
  delete d_varLabel;
}


particleIndex 
ParticleCreator::countParticles(const Uintah::Patch* patch,
                                std::vector<Uintah::GeometryObject*>& d_geom_objs)
{
  cout_doing << "\t Counting particles in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  d_lock.writeLock();
  particleIndex sum = 0;
  for (auto geom = d_geom_objs.begin(); geom != d_geom_objs.end(); ++geom){ 
    cout_dbg << "\t\t Calling countAndCreatePoints for patch " << patch << std::endl;
    sum += countAndCreateParticles(patch,*geom);
  }
  
  d_lock.writeUnlock();
  return sum;
}


particleIndex 
ParticleCreator::countAndCreateParticles(const Uintah::Patch* patch, 
                                         Uintah::GeometryObject* obj)
{
  cout_doing << "\t\t Counting and creating particles in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  GeometryPoints::key_type pointKey(patch, obj);
  GeometryScalars::key_type scalarKey(patch, obj);
  GeometryVectors::key_type vectorKey(patch, obj);

  Uintah::GeometryPieceP piece = obj->getPiece();
  Uintah::Box b1 = piece->getBoundingBox();
  Uintah::Box b2 = patch->getExtraBox();
  Uintah::Box b = b1.intersect(b2);
  if (b.degenerate()) return 0;
  
  // If the object is a FileGeometryPiece (e.g. FileGeometryPiece or
  // SmoothGeomPiece) then use the particle creators in that 
  // class to do the counting 
  Uintah::SmoothGeomPiece *sgp = dynamic_cast<Uintah::SmoothGeomPiece*>(piece.get_rep());
  if (sgp) {
    int numPts = 0;
    
    // Check whether it is actually a FileGeometryPiece
    Uintah::FileGeometryPiece* fgp = dynamic_cast<Uintah::FileGeometryPiece*>(sgp);
    if (fgp) {
      fgp->readPoints(patch->getID());
      numPts = fgp->returnPointCount();
    } else {
      SCIRun::Vector dxpp = patch->dCell()/obj->getInitialData_IntVector("res");
      double dx   = std::min(std::min(dxpp.x(),dxpp.y()), dxpp.z());
      sgp->setParticleSpacing(dx);
      numPts = sgp->createPoints();
      std::cout << "Smooth Geom Piece: Number of points created = " << numPts << std::endl;
    }

    std::vector<SCIRun::Point>* points = sgp->getPoints();
    std::vector<double>* vols = sgp->getVolume();
    std::vector<SCIRun::Vector>* pvelocities= sgp->getVelocity();
    std::vector<SCIRun::Vector>* pforces= sgp->getForces();

    SCIRun::Point pp;
    SCIRun::IntVector cell_idx;
    
    for (int ii = 0; ii < numPts; ++ii) {
      pp = points->at(ii);
      if (patch->findCell(pp, cell_idx)) {
        if (patch->containsPoint(pp)) {
          d_object_points[pointKey].push_back(pp);
          
          if (!vols->empty()) {
            double vol = vols->at(ii); 
            d_object_vols[scalarKey].push_back(vol);
          }
          if (!pvelocities->empty()) {
	    SCIRun::Vector pvel = pvelocities->at(ii); 
            d_object_velocity[vectorKey].push_back(pvel);
          }
          if (!pforces->empty()) {
	    SCIRun::Vector pforce = pforces->at(ii); 
            d_object_forces[vectorKey].push_back(pforce);
          }
        }  // patch contains cell
      }
    }
  } else { // Not a SmoothGeomPiece or FileGeometryPiece
    cout_dbg << "\t\t Calling createPoints for patch " << patch << std::endl;
    createPoints(patch,obj);
  }
  
  return (particleIndex) d_object_points[pointKey].size();
}

void 
ParticleCreator::createPoints(const Uintah::Patch* patch, 
                              Uintah::GeometryObject* obj)
{
  cout_doing << "\t\t Creating points in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  GeometryPoints::key_type key(patch, obj);

  Uintah::GeometryPieceP geom_piece = obj->getPiece();
  SCIRun::Box box_with_extra_cells = patch->getExtraBox();

  SCIRun::IntVector particles_per_cell = obj->getInitialData_IntVector("res");
  SCIRun::Vector dxpp = patch->dCell()/particles_per_cell;
  SCIRun::Vector dcorner = dxpp*0.5;

  // Iterate through cells in patch
  for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
  
    SCIRun::IntVector cell = *iter;
    SCIRun::Point lower = patch->nodePosition(cell) + dcorner;

    for (int ix = 0; ix < particles_per_cell.x(); ix++) {
      for (int iy = 0; iy < particles_per_cell.y(); iy++) {
        for (int iz = 0; iz < particles_per_cell.z(); iz++) {

          SCIRun::IntVector idx(ix, iy, iz);
          SCIRun::Point point = lower + dxpp*idx;
  
          if (!box_with_extra_cells.contains(point)) {
             throw Uintah::InternalError("Particle created outside of patch ?",
                                 __FILE__, __LINE__);
          }

          if (geom_piece->inside(point)) {
            d_object_points[key].push_back(point);
          }
          
        } 
      } 
    } 
  } // end cell iterator

  cout_dbg << "\t\t Number of points created in patch " << patch << " = " << d_object_points[key].size()
           << std::endl;

}

Uintah::ParticleSubset* 
ParticleCreator::createParticles(PeridynamicsMaterial* matl,
                                 particleIndex numParticles,
                                 Uintah::CCVariable<short int>& cellNAPID,
                                 const Uintah::Patch* patch,
                                 Uintah::DataWarehouse* new_dw,
                                 std::vector<Uintah::GeometryObject*>& d_geom_objs)
{
  cout_doing << "\t\t\t Creating particles: Peridynamics: " << __FILE__ << ":" << __LINE__ << std::endl;

  d_lock.writeLock();

  // Allocate the space for particle variables associated with this material
  int matlIndex = matl->getDWIndex();
  Uintah::ParticleSubset* subset = allocateVariables(numParticles, matlIndex, patch, new_dw);

  // Loop through geometry objects
  particleIndex start = 0;
  for (auto obj = d_geom_objs.begin(); obj != d_geom_objs.end(); ++obj) {
    particleIndex count = 0;

    // Get the geometry piece and its bounding box
    Uintah::GeometryPieceP piece = (*obj)->getPiece();
    Uintah::Box b1 = piece->getBoundingBox();
    Uintah::Box b2 = patch->getExtraBox();
    Uintah::Box b = b1.intersect(b2);
    if(b.degenerate()) {
      count = 0;
      continue;
    }

    // Set up arrays for particle volumes etc.
    std::vector<double>*         pVolArray = 0;   // particle volumes
    std::vector<SCIRun::Vector>* pVelArray = 0;   // particle velocities 
    std::vector<SCIRun::Vector>* pForceArray = 0; // particle forces 

    // Special case exception for SmoothGeomPieces and FileGeometryPieces
    // FileGeometryPieces are derived from SmoothGeomPiece and contain the particle data in a file
    // while smooth geometry pieces generate these same data.  
    // **WARNING** Not sure what the effect is on Abaqus type input files.
    Uintah::SmoothGeomPiece *gp = dynamic_cast<Uintah::SmoothGeomPiece*>(piece.get_rep());
    if (gp){
      pVolArray  = gp->getVolume();
      pVelArray  = gp->getVelocity();  
      pForceArray = gp->getForces();
    }

    // Get the distance from the actual boundary of the object to the
    // particle nearest the boundary.  This is needed later for locating
    // boundary particles
    SCIRun::Vector dxpp = patch->dCell()/(*obj)->getInitialData_IntVector("res");    

    // Set up iterator for getting particle volumes, velocities, etc. (if they exist)
    GeometryScalars::key_type scalarKey(patch, *obj);
    GeometryVectors::key_type vectorKey(patch, *obj);
    auto voliter = d_object_vols[scalarKey].begin();
    auto velocityiter = d_object_velocity[vectorKey].begin();  
    auto forceiter = d_object_forces[vectorKey].begin();  
    
    // Loop through object points
    GeometryPoints::key_type pointKey(patch, *obj);
    for ( auto itr = d_object_points[pointKey].begin(); itr!=d_object_points[pointKey].end(); ++itr) {

      SCIRun::Point pPosition = *itr;
      SCIRun::IntVector cell_idx;
      if (!patch->findCell(pPosition, cell_idx)) continue;
      if (!patch->containsPoint(pPosition)) continue;
      
      particleIndex pidx = start+count;      
      cout_dbg << "\t\t\t CreateParticles: Point["<<pidx<<"]="<<pPosition<<" Cell = "<<cell_idx<<std::endl;
 
      initializeParticle(patch, obj, matl, pPosition, cell_idx, pidx, cellNAPID);
      
      if (pVolArray) {
        if (!pVolArray->empty()) {
          d_pvolume[pidx] = *voliter;
          d_pmass[pidx] = matl->getInitialDensity()*d_pvolume[pidx];
          ++voliter;
        }
      }

      if (pVelArray) {                           
        if (!pVelArray->empty()) {               
          d_pvelocity[pidx] = *velocityiter;
          ++velocityiter;
        }
      }                                         

      if (pForceArray) {                           
        if (!pForceArray->empty()) {               
          d_pexternalforce[pidx] = *forceiter;
          ++forceiter;
        }
      }                                         

      // If the particle is on the surface and if there is
      // a physical BC attached to it then mark with the 
      // physical BC pointer
      if (d_flags->d_useLoadCurves) {
        if (checkForSurface(piece, pPosition, dxpp)) {
          d_pLoadCurveID[pidx] = getLoadCurveID(pPosition, dxpp);
        } else {
          d_pLoadCurveID[pidx] = 0;
        }
      }

      count++;
    }
    start += count;
  }
  d_lock.writeUnlock();
  return subset;
}

Uintah::ParticleSubset* 
ParticleCreator::allocateVariables(particleIndex numParticles, 
                                   int matlIndex, 
                                   const Uintah::Patch* patch,
                                   Uintah::DataWarehouse* new_dw)
{
  cout_doing << "\t\t\t Allocating variables in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  Uintah::ParticleSubset* subset = new_dw->createParticleSubset(numParticles,matlIndex,
                                                        patch);
  new_dw->allocateAndPut(d_pparticleID,    d_varLabel->pParticleIDLabel,    subset);
  new_dw->allocateAndPut(d_position,       d_varLabel->pPositionLabel,      subset);
  new_dw->allocateAndPut(d_pmass,          d_varLabel->pMassLabel,          subset);
  new_dw->allocateAndPut(d_psize,          d_varLabel->pSizeLabel,          subset);
  new_dw->allocateAndPut(d_pvolume,        d_varLabel->pVolumeLabel,        subset);
  new_dw->allocateAndPut(d_pdisp,          d_varLabel->pDisplacementLabel,  subset);
  new_dw->allocateAndPut(d_pvelocity,      d_varLabel->pVelocityLabel,      subset); 
  new_dw->allocateAndPut(d_pexternalforce, d_varLabel->pExternalForceLabel, subset); 
  new_dw->allocateAndPut(d_pHorizon,       d_varLabel->pHorizonLabel,       subset);

  if (d_flags->d_useLoadCurves) {
    new_dw->allocateAndPut(d_pLoadCurveID, d_varLabel->pLoadCurveIDLabel,   subset);
  }

  return subset;
}

void ParticleCreator::allocateVariablesAddRequires(Uintah::Task* task, 
                                                   const PeridynamicsMaterial* ,
                                                   const Uintah::PatchSet* ) const
{
  cout_doing << "\t Scheduling task variables in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;

  d_lock.writeLock();
  Uintah::Ghost::GhostType  gn = Uintah::Ghost::None;
  task->requires(Uintah::Task::OldDW, d_varLabel->pParticleIDLabel,  gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pPositionLabel,    gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pMassLabel,        gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pSizeLabel,        gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pDisplacementLabel, gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pVelocityLabel,    gn);

  task->requires(Uintah::Task::NewDW, d_varLabel->pVolumeLabel_preReloc, gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pHorizonLabel,  gn);

  if (d_flags->d_useLoadCurves){
    task->requires(Uintah::Task::OldDW, d_varLabel->pLoadCurveIDLabel, gn);
  }

  d_lock.writeUnlock();
}

void 
ParticleCreator::initializeParticle(const Uintah::Patch* patch,
                                    std::vector<Uintah::GeometryObject*>::const_iterator obj,
                                    PeridynamicsMaterial* matl,
                                    SCIRun::Point pPosition,
                                    SCIRun::IntVector cell_idx,
                                    particleIndex pidx,
                                    Uintah::CCVariable<short int>& cellNAPID)
{
  cout_doing << "\t\t\t\t Initializing particles in particle creator: Peridynamics: " 
             << __FILE__ << ":" << __LINE__ << std::endl;


  SCIRun::IntVector ppc = (*obj)->getInitialData_IntVector("res");
  cout_dbg << "\t\t\t\t Particles per cell = " << ppc << std::endl; 

  SCIRun::Vector dxcc = patch->dCell();
  cout_dbg << "\t\t\t\t Cell dimensions = " << dxcc << std::endl; 

  Uintah::Matrix3 size(1./((double) ppc.x()),0.,0.,
                       0.,1./((double) ppc.y()),0.,
                       0.,0.,1./((double) ppc.z()));
  cout_dbg << "\t\t\t\t Particle size = " << size << std::endl; 
  cout_dbg << "\t\t\t\t Particle location = " << pPosition << std::endl; 
  cout_dbg << "\t\t\t\t Particle index = " << pidx << std::endl; 

  d_position[pidx] = pPosition;
  d_pvolume[pidx]  = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  d_psize[pidx]      = size;
  d_pvelocity[pidx]  = (*obj)->getInitialData_Vector("velocity");
  d_pmass[pidx]      = matl->getInitialDensity()*d_pvolume[pidx];
  d_pdisp[pidx]      = SCIRun::Vector(0.,0.,0.);
  
  // Compute the max length of the side of a cell and set the horizon accordingly
  double maxCellEdge = std::max(std::max(dxcc.x(), dxcc.y()), dxcc.z());
  d_pHorizon[pidx] = maxCellEdge*d_flags->d_numCellsInHorizon;
  
  // Initialize external force (updated later using ParticleLoadBCs)
  SCIRun::Vector pExtForce(0.0, 0.0, 0.0);
  d_pexternalforce[pidx] = pExtForce;

  ASSERT(cell_idx.x() <= 0xffff && 
         cell_idx.y() <= 0xffff && 
         cell_idx.z() <= 0xffff);
         
  Uintah::long64 cellID = ((Uintah::long64)cell_idx.x() << 16) | 
                          ((Uintah::long64)cell_idx.y() << 32) | 
                          ((Uintah::long64)cell_idx.z() << 48);
                  
  short int& myCellNAPID = cellNAPID[cell_idx];
  d_pparticleID[pidx] = (cellID | (Uintah::long64) myCellNAPID);
  ASSERT(myCellNAPID < 0x7fff);
  myCellNAPID++;
}

int
ParticleCreator::checkForSurface( const Uintah::GeometryPieceP piece, const SCIRun::Point p,
                                  const SCIRun::Vector dxpp )
{

  //  Check the candidate points which surround the point just passed
  //   in.  If any of those points are not also inside the object
  //  the current point is on the surface

  int ss = 0;
  // Check to the left (-x)
  if(!piece->inside(p-SCIRun::Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check to the right (+x)
  if(!piece->inside(p+SCIRun::Vector(dxpp.x(),0.,0.)))
    ss++;
  // Check behind (-y)
  if(!piece->inside(p-SCIRun::Vector(0.,dxpp.y(),0.)))
    ss++;
  // Check in front (+y)
  if(!piece->inside(p+SCIRun::Vector(0.,dxpp.y(),0.)))
    ss++;
  // Check below (-z)
  if(!piece->inside(p-SCIRun::Vector(0.,0.,dxpp.z())))
    ss++;
  // Check above (+z)
  if(!piece->inside(p+SCIRun::Vector(0.,0.,dxpp.z())))
    ss++;

  if(ss>0){
    return 1;
  }
  else {
    return 0;
  }
}

// Get the LoadCurveID applicable for this material point
// WARNING : Should be called only once per particle during a simulation 
// because it updates the number of particles to which a BC is applied.
int 
ParticleCreator::getLoadCurveID(const SCIRun::Point& pp, const SCIRun::Vector& dxpp)
{
  int ret=0;
  for (unsigned int ii = 0; ii < ParticleLoadBCFactory::particleLoadBCs.size(); ii++){

    std::string bcs_type = ParticleLoadBCFactory::particleLoadBCs[ii]->getType();

    //cerr << " BC Type = " << bcs_type << endl;
    if (bcs_type == "Pressure") {
      ParticlePressureBC* bc =
        dynamic_cast<ParticlePressureBC*>(ParticleLoadBCFactory::particleLoadBCs[ii]);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
         ret = bc->loadCurveID();
      }
    }
    else if (bcs_type == "Force") {
      ParticleForceBC* bc =
        dynamic_cast<ParticleForceBC*>(ParticleLoadBCFactory::particleLoadBCs[ii]);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
        ret = bc->loadCurveID();
      }
    }
    else if (bcs_type == "NormalForce") {
      ParticleNormalForceBC* bc =
        dynamic_cast<ParticleNormalForceBC*>(ParticleLoadBCFactory::particleLoadBCs[ii]);
      if (bc->flagSurfaceParticle(pp, dxpp)) {
        ret = bc->loadCurveID();
      }
    }
  }
  return ret;
}


std::vector<const Uintah::VarLabel* > ParticleCreator::returnParticleState()
{
  return particle_state;
}


std::vector<const Uintah::VarLabel* > ParticleCreator::returnParticleStatePreReloc()
{
  return particle_state_preReloc;
}

void ParticleCreator::registerPermanentParticleState(PeridynamicsMaterial* matl)
{
  d_lock.writeLock();
  particle_state.push_back(d_varLabel->pDisplacementLabel);
  particle_state_preReloc.push_back(d_varLabel->pDisplacementLabel_preReloc);

  particle_state.push_back(d_varLabel->pVelocityLabel);
  particle_state_preReloc.push_back(d_varLabel->pVelocityLabel_preReloc);

  particle_state.push_back(d_varLabel->pExternalForceLabel);
  particle_state_preReloc.push_back(d_varLabel->pExternalForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pMassLabel);
  particle_state_preReloc.push_back(d_varLabel->pMassLabel_preReloc);

  particle_state.push_back(d_varLabel->pVolumeLabel);
  particle_state_preReloc.push_back(d_varLabel->pVolumeLabel_preReloc);

  particle_state.push_back(d_varLabel->pParticleIDLabel);
  particle_state_preReloc.push_back(d_varLabel->pParticleIDLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pSizeLabel);
  particle_state_preReloc.push_back(d_varLabel->pSizeLabel_preReloc);

  particle_state.push_back(d_varLabel->pDefGradLabel);
  particle_state_preReloc.push_back(d_varLabel->pDefGradLabel_preReloc);

  particle_state.push_back(d_varLabel->pStressLabel);
  particle_state_preReloc.push_back(d_varLabel->pStressLabel_preReloc);

  particle_state.push_back(d_varLabel->pPK1StressLabel);
  particle_state_preReloc.push_back(d_varLabel->pPK1StressLabel_preReloc);

  particle_state.push_back(d_varLabel->pHorizonLabel);
  particle_state_preReloc.push_back(d_varLabel->pHorizonLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pNeighborListLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborListLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pNeighborConnLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborConnLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pNeighborCountLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborCountLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pNeighborBondForceLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborBondForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pNeighborBondEnergyLabel);
  particle_state_preReloc.push_back(d_varLabel->pNeighborBondEnergyLabel_preReloc);

  particle_state.push_back(d_varLabel->pInternalForceLabel);
  particle_state_preReloc.push_back(d_varLabel->pInternalForceLabel_preReloc);

  particle_state.push_back(d_varLabel->pDamageLabel);
  particle_state_preReloc.push_back(d_varLabel->pDamageLabel_preReloc);

  if (d_flags->d_useLoadCurves) {
    particle_state.push_back(d_varLabel->pLoadCurveIDLabel);
    particle_state_preReloc.push_back(d_varLabel->pLoadCurveIDLabel_preReloc);
  }

  matl->getMaterialModel()->addParticleState(particle_state,
                                             particle_state_preReloc);

  matl->getDamageModel()->addParticleState(particle_state,
                                           particle_state_preReloc);

  d_lock.writeUnlock();
}

