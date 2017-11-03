#ifndef BOND_H
#define BOND_H

#include <Core/Math/Matrix.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/PeriContainers.h>
#include <InputOutput/InputParameter.h>
#include <boost/mpi.hpp>

#include <iostream>
#include <string>

namespace pd {

class PeriBond
{

public:
  //-------------------------------------------------------------------------
  // Default Constructor
  PeriBond();

  // Overload Constructor
  PeriBond(REAL, PeriParticleP, PeriParticleP);

  // Destructor
  ~PeriBond();

  //-------------------------------------------------------------------------
  // Accessor Functions
  bool getIsAlive() const;
  // getisAlive - returns state of the PeriBond
  // @return - state of the PeriBond

  bool getIsRecv() const { return isRecv; }
  void setIsRecv() { isRecv = true; }

  REAL getWeight() const;
  // getweight - returns weight of the PeriBond
  // @return - weight of the PeriBond

  REAL getInitLength() const;
  // getinitLength - returns initial length of the PeriBond
  // @return - initial length of the PeriBond

  REAL volume(bool) const;

  dem::Vec getXi(bool) const;

  dem::Vec getEta(bool) const;

  dem::Vec getEtaHalf(bool, const REAL) const;

  PeriParticleP getPt1() const { return pt1; }

  PeriParticleP getPt2() const { return pt2; }

  dem::Matrix getMicroK(
    const bool) const; // get the contribution of K from one single bond

  dem::Matrix getMicroN(const bool, const bool) const;

  dem::Matrix getMicroNHalf(const bool, const bool, const REAL) const;

  dem::Matrix getMicroNDeltaU(const bool, const bool, const REAL) const;

  //-------------------------------------------------------------------------
  // Mutator Functions
  void setIsAlive(bool);
  // setisAlive - sets state of the PeriBond
  // @param - state of the PeriBond

  void setWeight(REAL);
  // setweight - sets weight of the PeriBond
  // @param - weight of the PeriBond

  void setInitLength(REAL);
  // setinitLength - sets initial length of the PeriBond
  // @param - initial length of the PeriBond

  void setAliveFalse() { isAlive = false; }

  //-------------------------------------------------------------------------
  // Utility Functions

  REAL calcCurrentLength();

  void checkIfAlive();

private:
  // Member Variables
  bool isAlive;      // if the PeriBond is alive or not
  bool isRecv;       // is this peri-bonds between recvPeriParticle
  REAL weight;       // influence function
  REAL initLength;   // initial PeriBond length
  PeriParticleP pt1; // what for? store address of the particles it belongs to.
  PeriParticleP pt2;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& isAlive;
    ar& isRecv;
    ar& weight;
    ar& initLength;
    //      ar & pt1;
    //      ar & pt2;
  }
}; // end PeriBond

} // end pd
#endif
