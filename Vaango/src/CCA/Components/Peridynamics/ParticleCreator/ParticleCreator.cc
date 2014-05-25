#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarLabel.h>

#include <fstream>
#include <iostream>

using namespace Vaango;

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

Uintah::ParticleSubset* 
ParticleCreator::createParticles(PeridynamicsMaterial* matl,
                                 particleIndex numParticles,
                                 Uintah::CCVariable<short int>& cellNAPID,
                                 const Uintah::Patch* patch,
                                 Uintah::DataWarehouse* new_dw,
                                 std::vector<Uintah::GeometryObject*>& d_geom_objs)
{
  d_lock.writeLock();

  int dwi = matl->getDWIndex();
  Uintah::ParticleSubset* subset = allocateVariables(numParticles,dwi,patch,new_dw);

  particleIndex start = 0;
  
  std::vector<Uintah::GeometryObject*>::const_iterator obj;
  for (obj = d_geom_objs.begin(); obj != d_geom_objs.end(); ++obj) {
    particleIndex count = 0;
    Uintah::GeometryPieceP piece = (*obj)->getPiece();
    Uintah::Box b1 = piece->getBoundingBox();
    Uintah::Box b2 = patch->getExtraBox();
    Uintah::Box b = b1.intersect(b2);
    if(b.degenerate()) {
      count = 0;
      continue;
    }

    //SCIRun::Vector dxpp = patch->dCell()/(*obj)->getInitialData_IntVector("res");    

    // Special case exception for FileGeometryPieces and FileGeometryPieces
    Uintah::FileGeometryPiece *fgp = dynamic_cast<Uintah::FileGeometryPiece*>(piece.get_rep());
    std::vector<double>* volumes       = 0;
    std::vector<SCIRun::Vector>* pvelocities   = 0;    // gcd adds and new change name
    if (fgp){
      volumes      = fgp->getVolume();
      pvelocities  = fgp->getVelocity();  // gcd adds and new change name
    }

    // For getting particle volumes (if they exist)
    std::vector<double>::const_iterator voliter;
    geomvols::key_type volkey(patch,*obj);
    if (volumes) {
      if (!volumes->empty()) voliter = d_object_vols[volkey].begin();
    }

    // For getting particle velocities (if they exist)   // gcd adds
    std::vector<SCIRun::Vector>::const_iterator velocityiter;
    geomvecs::key_type pvelocitykey(patch,*obj);
    if (pvelocities) {                             // new change name
      if (!pvelocities->empty()) velocityiter =
              d_object_velocity[pvelocitykey].begin();  // new change name
    }                                                    // end gcd adds
    
    std::vector<SCIRun::Point>::const_iterator itr;
    geompoints::key_type key(patch,*obj);
    for(itr=d_object_points[key].begin();itr!=d_object_points[key].end();++itr){
      SCIRun::IntVector cell_idx;
      if (!patch->findCell(*itr,cell_idx)) continue;

      if (!patch->containsPoint(*itr)) continue;
      
      particleIndex pidx = start+count;      
      //cerr << "Point["<<pidx<<"]="<<*itr<<" Cell = "<<cell_idx<<endl;
 
      initializeParticle(patch,obj,matl,*itr,cell_idx,pidx,cellNAPID);
      
      if (volumes) {
        if (!volumes->empty()) {
          pvolume[pidx] = *voliter;
          pmass[pidx] = matl->getInitialDensity()*pvolume[pidx];
          ++voliter;
        }
      }

      if (pvelocities) {                           // gcd adds and change name 
        if (!pvelocities->empty()) {               // and change name
          pvelocity[pidx] = *velocityiter;
          ++velocityiter;
        }
      }                                         // end gcd adds

      count++;
    }
    start += count;
  }
  d_lock.writeUnlock();
  return subset;
}

Uintah::ParticleSubset* 
ParticleCreator::allocateVariables(particleIndex numParticles, 
                                   int dwi, const Uintah::Patch* patch,
                                   Uintah::DataWarehouse* new_dw)
{

  Uintah::ParticleSubset* subset = new_dw->createParticleSubset(numParticles,dwi,
                                                        patch);
  new_dw->allocateAndPut(pparticleID,    d_varLabel->pParticleIDLabel,    subset);
  new_dw->allocateAndPut(px,             d_varLabel->pXLabel,             subset);
  new_dw->allocateAndPut(position,       d_varLabel->pPositionLabel,      subset);
  new_dw->allocateAndPut(pmass,          d_varLabel->pMassLabel,          subset);
  new_dw->allocateAndPut(psize,          d_varLabel->pSizeLabel,          subset);
  new_dw->allocateAndPut(pvolume,        d_varLabel->pVolumeLabel,        subset);
  new_dw->allocateAndPut(pdisp,          d_varLabel->pDisplacementLabel,  subset);
  new_dw->allocateAndPut(pvelocity,      d_varLabel->pVelocityLabel,      subset); 

  new_dw->allocateAndPut(pHorizon,       d_varLabel->pHorizonLabel,       subset);
  
  return subset;
}

void ParticleCreator::allocateVariablesAddRequires(Uintah::Task* task, 
                                                   const PeridynamicsMaterial* ,
                                                   const Uintah::PatchSet* ) const
{
  d_lock.writeLock();
  Uintah::Ghost::GhostType  gn = Uintah::Ghost::None;
  task->requires(Uintah::Task::OldDW, d_varLabel->pParticleIDLabel,  gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pPositionLabel,    gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pXLabel,           gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pMassLabel,        gn);
  //task->requires(Task::OldDW,d_varLabel->pVolumeLabel,    gn);
  task->requires(Uintah::Task::NewDW, d_varLabel->pVolumeLabel_preReloc,   gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pDisplacementLabel, gn);
  task->requires(Uintah::Task::OldDW, d_varLabel->pVelocityLabel,    gn);

  task->requires(Uintah::Task::OldDW, d_varLabel->pHorizonLabel,  gn);
  d_lock.writeUnlock();
}


void ParticleCreator::createPoints(const Uintah::Patch* patch, Uintah::GeometryObject* obj)
{

}


void 
ParticleCreator::initializeParticle(const Uintah::Patch* patch,
                                    std::vector<Uintah::GeometryObject*>::const_iterator obj,
                                    PeridynamicsMaterial* matl,
                                    SCIRun::Point p,
                                    SCIRun::IntVector cell_idx,
                                    particleIndex i,
                                    Uintah::CCVariable<short int>& cellNAPID)
{
  SCIRun::IntVector ppc = (*obj)->getInitialData_IntVector("res");
  //SCIRun::Vector dxpp = patch->dCell()/(*obj)->getInitialData_IntVector("res");
  SCIRun::Vector dxcc = patch->dCell();
  Uintah::Matrix3 size(1./((double) ppc.x()),0.,0.,
                       0.,1./((double) ppc.y()),0.,
                       0.,0.,1./((double) ppc.z()));

  px[i] = p;
  position[i] = p;
  pvolume[i]  = size.Determinant()*dxcc.x()*dxcc.y()*dxcc.z();
  psize[i]      = size;
  pvelocity[i]  = (*obj)->getInitialData_Vector("velocity");
  pmass[i]      = matl->getInitialDensity()*pvolume[i];
  pdisp[i]      = SCIRun::Vector(0.,0.,0.);
  
  // Compute the max length of the side of a cell
  double maxCellEdge = std::max(std::max(dxcc.x(), dxcc.y()), dxcc.z());
  pHorizon[i] = maxCellEdge*d_flags->d_numCellsInHorizon;
  
  ASSERT(cell_idx.x() <= 0xffff && 
         cell_idx.y() <= 0xffff && 
         cell_idx.z() <= 0xffff);
         
  Uintah::long64 cellID = ((Uintah::long64)cell_idx.x() << 16) | 
                  ((Uintah::long64)cell_idx.y() << 32) | 
                  ((Uintah::long64)cell_idx.z() << 48);
                  
  short int& myCellNAPID = cellNAPID[cell_idx];
  pparticleID[i] = (cellID | (Uintah::long64) myCellNAPID);
  ASSERT(myCellNAPID < 0x7fff);
  myCellNAPID++;
}

particleIndex 
ParticleCreator::countParticles(const Uintah::Patch* patch,
                                std::vector<Uintah::GeometryObject*>& d_geom_objs)
{
  d_lock.writeLock();
  particleIndex sum = 0;
  std::vector<Uintah::GeometryObject*>::const_iterator geom;
  for (geom=d_geom_objs.begin(); geom != d_geom_objs.end(); ++geom){ 
    sum += countAndCreateParticles(patch,*geom);
  }
  
  d_lock.writeUnlock();
  return sum;
}


particleIndex 
ParticleCreator::countAndCreateParticles(const Uintah::Patch* patch, 
                                         Uintah::GeometryObject* obj)
{
  geompoints::key_type key(patch,obj);
  geomvols::key_type   volkey(patch,obj);
  geomvecs::key_type   pvelocitykey(patch,obj);
  Uintah::GeometryPieceP piece = obj->getPiece();
  Uintah::Box b1 = piece->getBoundingBox();
  Uintah::Box b2 = patch->getExtraBox();
  Uintah::Box b = b1.intersect(b2);
  if(b.degenerate()) return 0;
  
  // If the object is a FileGeometryPiece (e.g. FileGeometryPiece or
  // SmoothCylGeomPiece) then use the particle creators in that 
  // class to do the counting d
  Uintah::FileGeometryPiece *fgp = dynamic_cast<Uintah::FileGeometryPiece*>(piece.get_rep());
  if (fgp) {
    fgp->readPoints(patch->getID());
    int numPts = fgp->returnPointCount();

    std::vector<SCIRun::Point>* points = fgp->getPoints();
    std::vector<double>* vols = fgp->getVolume();
    std::vector<SCIRun::Vector>* pvelocities= fgp->getVelocity();
    SCIRun::Point p;
    SCIRun::IntVector cell_idx;
    
    for (int ii = 0; ii < numPts; ++ii) {
      p = points->at(ii);
      if (patch->findCell(p,cell_idx)) {
        if (patch->containsPoint(p)) {
          d_object_points[key].push_back(p);
          
          if (!vols->empty()) {
            double vol = vols->at(ii); 
            d_object_vols[volkey].push_back(vol);
          }
          if (!pvelocities->empty()) {
	    SCIRun::Vector pvel = pvelocities->at(ii); 
            d_object_velocity[pvelocitykey].push_back(pvel);
          }
        }  // patch contains cell
      }
    }
  } else {
    createPoints(patch,obj);
  }
  
  return (particleIndex) d_object_points[key].size();
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

  particle_state.push_back(d_varLabel->pMassLabel);
  particle_state_preReloc.push_back(d_varLabel->pMassLabel_preReloc);

  particle_state.push_back(d_varLabel->pVolumeLabel);
  particle_state_preReloc.push_back(d_varLabel->pVolumeLabel_preReloc);

  particle_state.push_back(d_varLabel->pParticleIDLabel);
  particle_state_preReloc.push_back(d_varLabel->pParticleIDLabel_preReloc);
  
  particle_state.push_back(d_varLabel->pDispGradLabel);
  particle_state_preReloc.push_back(d_varLabel->pDispGradLabel_preReloc);

  particle_state.push_back(d_varLabel->pVelGradLabel);
  particle_state_preReloc.push_back(d_varLabel->pVelGradLabel_preReloc);

  particle_state.push_back(d_varLabel->pDefGradLabel);
  particle_state_preReloc.push_back(d_varLabel->pDefGradLabel_preReloc);

  particle_state.push_back(d_varLabel->pStressLabel);
  particle_state_preReloc.push_back(d_varLabel->pStressLabel_preReloc);

  matl->getMaterialModel()->addParticleState(particle_state,
                                             particle_state_preReloc);

  matl->getDamageModel()->addParticleState(particle_state,
                                           particle_state_preReloc);

  d_lock.writeUnlock();
}

