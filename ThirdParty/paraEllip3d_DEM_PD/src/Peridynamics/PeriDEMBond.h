#ifndef PERIDEMBOND_H
#define PERIDEMBOND_H

#include <Core/Const/const.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/DEMParticle.h>
#include <InputOutput/InputParameter.h>
#include <Peridynamics/PeriParticle.h>
#include <boost/mpi.hpp>

// the bond between peri-points and the sand particles
// July 14, 2014
using dem::Vec;
using dem::DEMParticle;

namespace pd {

class PeriDEMBond
{

public:
  PeriDEMBond()
  {
    initProjectorLocal = 0;
    initBondVec = 0;
    currBondVec = initBondVec;
    isAlive = true;
    periPoint = NULL;
    demParticle = NULL;
  }

  PeriDEMBond(const Vec& alocal, DEMParticle* dem_pt,
              PeriParticle* peri_pt)
  {
    initProjectorLocal = alocal;

    Vec aglobal = dem_pt->currentPosition() + dem_pt->localToGlobal(alocal);

    initBondVec = peri_pt->currentPosition() - aglobal;
    currBondVec = initBondVec;
    isAlive = true;
    periPoint = peri_pt;
    demParticle = dem_pt;

    periPoint->addDEMId(demParticle->getId());
  }

  ~PeriDEMBond()
  {
    periPoint->eraseDEMId(demParticle->getId());
    periPoint = NULL;
    demParticle = NULL;
  }

  //    void eraseDEMIDFromPeri(){
  //	periPoint->eraseDEMID(demParticle->getID());
  //    }

  /*
  // this is used to test the coupled force model, October 10, 2014
  // in this test model, the sand-peri-points will move along the dem-particle
      PeriDEMBond(const Vec& alocal, particle* dem_pt,
  PeriParticle* peri_pt){
        initProjectorLocal = alocal;	// now alocal the peri-point itself
  instead of
  the projector

        isAlive = true;
        periPoint = peri_pt;
        demParticle = dem_pt;
      }
  */

  void applyBondForce();
  void applyBondBoundary(); // this is used to test the coupled force model,
                            // October 10, 2014
  // in this test model, the sand-peri-points will move along the dem-particle
  bool getIsAlive() { return isAlive; }

private:
  // bond has two ends, one end is the peri-point,
  // the other is the projector of this peri-point on the surface of the sand
  // particle,
  // this initProjectorLocal is the position of this projector point in the
  // local coordinates
  // of this sand particle
  Vec initProjectorLocal;
  Vec initBondVec; // initial bond vector pointing from the projector to the
                   // peri-point
  Vec currBondVec; // current bond vector pointing from the current surface
                   // point to the peri-point
  bool isAlive;    // the state if the bond is alive
  PeriParticle* periPoint;
  DEMParticle* demParticle;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& initProjectorLocal;
    ar& initBondVec;
    ar& currBondVec;
    ar& isAlive;
    ar& periPoint;
    ar& demParticle;
  }
};

} // end pd

#endif
