#include <CCA/Components/Peridynamics/PeridynamicsDomainBoundCond.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <vector>
#include <iostream>

using namespace Vaango;

PeridynamicsDomainBoundCond::PeridynamicsDomainBoundCond()
{
}

PeridynamicsDomainBoundCond::~PeridynamicsDomainBoundCond()
{
}

void 
PeridynamicsDomainBoundCond::setBoundaryCondition(const Uintah::Patch* patch,int dwi,
                                                  const std::string& type, 
                                                  Uintah::NCVariable<SCIRun::Vector>& variable,
                                                  std::string interp_type)
{
  for(Uintah::Patch::FaceType face = Uintah::Patch::startFace;
      face <= Uintah::Patch::endFace; face=Uintah::Patch::nextFace(face)){
    SCIRun::IntVector oneCell = patch->faceDirection(face);

    if (patch->getBCType(face) == Uintah::Patch::None) {
      int numChildren = patch->getBCDataArray(face)->getNumberChildren(dwi);
      SCIRun::IntVector l(0,0,0),h(0,0,0),off(0,0,0);
      if(interp_type=="gimp" || interp_type=="3rdorderBS" || interp_type=="cpdi"){
        patch->getFaceExtraNodes(face,0,l,h);
      }
      for (int child = 0; child < numChildren; child++) {
        Uintah::Iterator nbound_ptr;
        Uintah::Iterator nu;        // not used;

        if (type == "Velocity"){
          const  Uintah::BoundCondBase* bcb = 
             patch->getArrayBCValues(face,dwi,"Velocity",nu,nbound_ptr,child);

          const Uintah::BoundCond<SCIRun::Vector>* bc = 
             dynamic_cast<const Uintah::BoundCond<SCIRun::Vector>*>(bcb); 
          if (bc != 0) {
            if (bc->getBCType__NEW() == "Dirichlet") {
              SCIRun::Vector bcv = bc->getValue();
              for (nbound_ptr.reset();!nbound_ptr.done();nbound_ptr++){ 
                SCIRun::IntVector nd = *nbound_ptr;
                variable[nd] = bcv;
              }
              if(interp_type=="gimp" || interp_type=="3rdorderBS" 
                                     || interp_type=="cpdi"){
                for(Uintah::NodeIterator it(l,h); !it.done(); it++) {
                  SCIRun::IntVector nd = *it;
                  variable[nd] = bcv;
                }
              }
            }
            delete bc;
          } else
          delete bcb;

        } else if (type == "Symmetric"){
          const Uintah::BoundCondBase* bcb =
            patch->getArrayBCValues(face,dwi,"Symmetric",nu,nbound_ptr,child);

          if (bcb->getBCType__NEW() == "symmetry") {
            if (face == Uintah::Patch::xplus || face == Uintah::Patch::xminus){
              if(interp_type=="linear"){
               for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++) {
                SCIRun::IntVector nd = *nbound_ptr;
                variable[nd] = SCIRun::Vector(0.,variable[nd].y(), variable[nd].z());
               }
              } // linear
              if(interp_type=="gimp" || interp_type=="cpdi" 
                                     || interp_type=="3rdorderBS"){
                SCIRun::IntVector off = SCIRun::IntVector(1,0,0);
                SCIRun::IntVector L(0,0,0),H(0,0,0);
                SCIRun::IntVector inner;
                if(face==Uintah::Patch::xminus){
                  L = l+off; H = h+off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(0.,variable[nd].y(), variable[nd].z());
                  }
                } else if(face==Uintah::Patch::xplus){
                  L = l-off; H = h-off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(0.,variable[nd].y(), variable[nd].z());
                  }
                }
                if(face==Uintah::Patch::xminus){
                  inner = SCIRun::IntVector(2,0,0);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { //extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(-variable[nd+inner].x(),
                                           variable[nd+inner].y(), 
                                           variable[nd+inner].z());
                  }
                } else if(face==Uintah::Patch::xplus){
                  inner = SCIRun::IntVector(-2,0,0);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { //extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(-variable[nd+inner].x(),
                                           variable[nd+inner].y(), 
                                           variable[nd+inner].z());
                  }
                }
              }  // cpdi, gimp or 3rdorderBS
            } // xplus/xminus faces

            if (face == Uintah::Patch::yplus || face == Uintah::Patch::yminus){
              if(interp_type=="linear"){
               for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++){
                SCIRun::IntVector nd = *nbound_ptr;
                variable[nd] = SCIRun::Vector(variable[nd].x(),0.,variable[nd].z());
               }
              } // linear
              if(interp_type=="gimp" || interp_type=="cpdi" 
                                     || interp_type=="3rdorderBS"){
                SCIRun::IntVector off = SCIRun::IntVector(0,1,0);
                SCIRun::IntVector L(0,0,0),H(0,0,0);
                SCIRun::IntVector inner;
                if(face==Uintah::Patch::yminus){
                  L = l+off; H = h+off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(variable[nd].x(),0.,variable[nd].z());
                  }
                } else if(face==Uintah::Patch::yplus){
                  L = l-off; H = h-off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(variable[nd].x(),0.,variable[nd].z());
                  }
                }
                if(face==Uintah::Patch::yminus){
                  inner = SCIRun::IntVector(0,2,0);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { // extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(variable[nd+inner].x(),
                                         -variable[nd+inner].y(),
                                          variable[nd+inner].z());
                  }
                } else if(face==Uintah::Patch::yplus){
                  inner = SCIRun::IntVector(0,-2,0);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { // extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(variable[nd+inner].x(),
                                         -variable[nd+inner].y(),
                                          variable[nd+inner].z());
                  }
                }
              } // cpdi or gimp
            }  // yplus/yminus faces
            if (face == Uintah::Patch::zplus || face == Uintah::Patch::zminus){
              if(interp_type=="linear"){
               for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++){
                SCIRun::IntVector nd = *nbound_ptr;
                variable[nd] = SCIRun::Vector(variable[nd].x(), variable[nd].y(),0.);
               }
              } // linear
              if(interp_type=="gimp" || interp_type=="cpdi" 
                                     || interp_type=="3rdorderBS"){
                SCIRun::IntVector off = SCIRun::IntVector(0,0,1);
                SCIRun::IntVector L(0,0,0),H(0,0,0);
                SCIRun::IntVector inner;
                if(face==Uintah::Patch::zminus){
                  L = l+off; H = h+off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(variable[nd].x(), variable[nd].y(),0.);
                  }
                } else if(face==Uintah::Patch::zplus){
                  L = l-off; H = h-off;
                  for(Uintah::NodeIterator it(L,H); !it.done(); it++){//bndy face nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd]=SCIRun::Vector(variable[nd].x(), variable[nd].y(),0.);
                  }
                }
                if(l.z()==-1 || h.z()==3){
                 if(face==Uintah::Patch::zminus){
                  inner = SCIRun::IntVector(0,0,2);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { // extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(variable[nd+inner].x(),
                                          variable[nd+inner].y(),
                                         -variable[nd+inner].z());
                  }
                 } else if(face==Uintah::Patch::zplus){
                  inner = SCIRun::IntVector(0,0,-2);
                  for(Uintah::NodeIterator it(l,h); !it.done(); it++) { // extra nodes
                    SCIRun::IntVector nd = *it;
                    variable[nd] = SCIRun::Vector(variable[nd+inner].x(),
                                          variable[nd+inner].y(),
                                         -variable[nd+inner].z());
                  }
                 }
                }
              } // cpdi or gimp
            } // zplus/zminus
            delete bcb;
          } else{
            delete bcb;
          }
        }
      }
    } else
      continue;
  }
}

void 
PeridynamicsDomainBoundCond::setBoundaryCondition(const Uintah::Patch* patch, int dwi,
                                                  const std::string& type, 
                                                  Uintah::NCVariable<double>& variable,
                                                  std::string interp_type)
{
  for(Uintah::Patch::FaceType face = Uintah::Patch::startFace;
      face <= Uintah::Patch::endFace; face=Uintah::Patch::nextFace(face)){
    SCIRun::IntVector oneCell = patch->faceDirection(face);
    if (patch->getBCType(face) == Uintah::Patch::None) {
      int numChildren = patch->getBCDataArray(face)->getNumberChildren(dwi);
      SCIRun::IntVector l(0,0,0),h(0,0,0);
      if(interp_type=="gimp" || interp_type=="3rdorderBS" || interp_type=="cpdi"){
        patch->getFaceExtraNodes(face,0,l,h);
      }
      for (int child = 0; child < numChildren; child++) {
        Uintah::Iterator nbound_ptr;
        Uintah::Iterator nu;  // not used

        if(type=="Pressure" || type=="Temperature"){
          const Uintah::BoundCondBase *bcb = 
            patch->getArrayBCValues(face,dwi,type,nu,nbound_ptr, child);

          const Uintah::BoundCond<double>* bc = 
            dynamic_cast<const Uintah::BoundCond<double>*>(bcb);
          
          if (bc != 0) {
            if (bc->getBCType__NEW() == "Dirichlet") {
              double bcv = bc->getValue();
              for (nbound_ptr.reset(); !nbound_ptr.done();nbound_ptr++){
                SCIRun::IntVector nd = *nbound_ptr;
                variable[nd] = bcv;
              }
              if(interp_type=="gimp" || interp_type=="3rdorderBS" || interp_type=="cpdi"){
                for(Uintah::NodeIterator it(l,h); !it.done(); it++) {
                  SCIRun::IntVector nd = *it;
                  variable[nd] = bcv;
                }
              }
            }
            
            if (bc->getBCType__NEW() == "Neumann"){
              SCIRun::Vector deltax = patch->dCell();
              double dx = -9;
              SCIRun::IntVector off(-9,-9,-9);
              if (face == Uintah::Patch::xplus){
                dx = deltax.x();
                off=SCIRun::IntVector(1,0,0);
              }
              else if (face == Uintah::Patch::xminus){
                dx = deltax.x();
                off=SCIRun::IntVector(-1,0,0);
              }
              else if (face == Uintah::Patch::yplus){
                dx = deltax.y();
                off=SCIRun::IntVector(0,1,0);
              }
              else if (face == Uintah::Patch::yminus){
                dx = deltax.y();
                off=SCIRun::IntVector(0,-1,0);
              }
              else if (face == Uintah::Patch::zplus){
                dx = deltax.z();
                off=SCIRun::IntVector(0,0,1);
              }
              else if (face == Uintah::Patch::zminus){
                dx = deltax.z();
                off=SCIRun::IntVector(0,0,-1);
              }
              
              double gradv = bc->getValue();

              for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
		SCIRun::IntVector nd = *nbound_ptr;
		variable[nd] = variable[nd-off] - gradv*dx;
	      }

              for(Uintah::NodeIterator it(l,h); !it.done(); it++) {
                SCIRun::IntVector nd = *it;
                variable[nd] = variable[nd-off] - gradv*dx;
              }
            }
            
            delete bc;
          } else
          delete bcb;
        }
      }  // child
    } else
      continue;
  }
}
