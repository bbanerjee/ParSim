#include <Peridynamics/PeriDEMBond.h>
#include <Core/Util/Utility.h>
#include <iostream>

namespace pd {

void
PeriDEMBond::applyBondForce()
{
  if (!isAlive) { // not alive
    return;
  }

  // (1) calculate current bond vector
  // get current global coordinates of the initProjectorLocal
  Vec currProjectorGlobal =
    demParticle->localToGlobal(initProjectorLocal) + demParticle->currentPosition();
  currBondVec = periPoint->currentPosition() - currProjectorGlobal;

  // (2) check if bond is alive
  // at present, use the same criterioin as the peri-bond used in pd

  // do not allow the peri-dem-bond broken at inclusion example
  REAL stretch = (vnormL2(currBondVec) - vnormL2(initBondVec)) / vnormL2(initBondVec);
  if (stretch > util::getParam<REAL>("bondStretchLimit")) {
    isAlive = false;
    return; // do not need to calculate forces
  }

  // (3) calculate bond force and apply bond force to peri-particle and
  // dem-particle
  // normal vector of bond w.r.t. initBondVec  (nn).b = (b.n)n
  Vec bondn = dot(currBondVec, normalize(initBondVec)) * normalize(initBondVec);        
  // tangent vector of bond w.r.t. initBondVec (I- nn).b = b - (b.n)n
  Vec bondt = currBondVec - bondn; 

  // only one k, kn=kt, which makes the reaction force is
  // more likely in velocity direction
  REAL k_periBndry = dem::plan_gravity * 40; 
  //	REAL kt_periBndry =
  // plan_gravity/(1+util::getParam<REAL>("bondStretchLimit"));
  // //
  // just for test, July 15, 2014

  REAL stretch_normal = vnormL2(bondn) - vnormL2(initBondVec);
  if (stretch_normal > 0)
    k_periBndry =
      dem::plan_gravity /
      (1 + util::getParam<REAL>("bondStretchLimit")) * 10;

  Vec fn =
    (bondn - initBondVec) *
    k_periBndry; // force is pointing from the projector to the peri-point
  Vec ft =
    bondt *
    k_periBndry; // force is pointing from the projector to the peri-point
  // //std::cout << "fn: " << fn.x() << ", " << fn.y() << ", " << fn.z()
  //<< std::endl;
  // //std::cout << "ft: " << ft.x() << ", " << ft.y() << ", " << ft.z()
  //<< std::endl << std::endl;
  // apply forces to peri-point
  periPoint->addAccelerationByForce(-fn - ft);
  // apply forces to dem-particle
  demParticle->addForce(fn + ft);
  //	demParticle->addMoment(
  //(currProjectorGlobal-demParticle->currentPosition())*(fn+ft) );

} // end applyBondForce

/*
    // this is only for spherical DEM particles, approximation method, but save
a lot of computations,
    // no, just save one minute for 100 time steps
    void PeriDEMBond::applyBondForce(){
        if(!isAlive){ // not alive
            return;
        }

        // (1) calculate current bond vector
        // get current global coordinates of the initProjectorLocal
        REAL currBondL =
vnormL2(periPoint->currentPosition()-demParticle->currentPosition())
- demParticle->radiusA();
        REAL initBondL = vnormL2(initBondVec);

        // (2) check if bond is alive
        // at present, use the same criterioin as the peri-bond used in
pd

// do not allow the peri-dem-bond broken at inclusion example
        REAL stretch = (currBondL - initBondL)/initBondL;
        if( stretch >
util::getParam<REAL>("bondStretchLimit")){
            isAlive = false;
            return;	// do not need to calculate forces
        }

        // (3) calculate bond force and apply bond force to peri-particle and
dem-particle
        REAL kn_periBndry =
util::getParam<REAL>("periYoung")*1e-1;	// just in value
        Vec fn = (currBondL-initBondL)*kn_periBndry*normalize(initBondVec);
// force is
pointing from the projector to the peri-point

        // apply forces to peri-point
        periPoint->addAccelerationByForce(-fn);
        // apply forces to dem-particle
        demParticle->addForce(fn);
//	demParticle->addMoment(
(currProjectorGlobal-demParticle->currentPosition())*(fn+ft) );

    } // end applyBondForce
*/

/*
    void PeriDEMBond::applyBondForce(){

        Vec peri_posi = periPoint->currentPosition();
        Vec dem_posi  = demParticle->currentPosition();
        Vec r_vec = dem_posi-peri_posi;

        REAL dist = vnormL2(r_vec);
        REAL radi = demParticle->radiusA();
        if(dist > radi)
            return;
        REAL maxAllowableRelativeOverlap =
util::getParam<REAL>("maxAllowableRelativeOverlap");
        if(dist < (1-maxAllowableRelativeOverlap)*radi)
            dist = (1-maxAllowableRelativeOverlap)*radi;

        Vec f_dem =
(radi-dist)*normalize(r_vec)*util::getParam<REAL>("periYoung")*1e-2;
//
forces by peri-points on DEM particle
////std::cout << "f_dem: " << f_dem.x() << ", " << f_dem.y() << ", " <<
f_dem.z() << std::endl;
        // apply forces to peri-point
        periPoint->addAccelerationByForce(f_dem);
        // apply forces to dem-particle
        demParticle->addForce(f_dem);

    } // end applyBondForce
*/

// this is used to test the coupled force model, October 10, 2014
// in this test model, the sand-peri-points will move along the dem-particle
void
PeriDEMBond::applyBondBoundary()
{
  if (!isAlive) { // not alive
    return;
  }

  // now initProjectorLocal is the local coordinate of the peri-point as the
  // bond created
  dem::Vec curr_posi_global =
    demParticle->localToGlobal(initProjectorLocal) + demParticle->currentPosition();
  periPoint->setCurrentPositionition(curr_posi_global);

} // end applyBondBoundary

} // end pd
