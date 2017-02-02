#ifndef PERIBOUNDARYBOND_H
#define PERIBOUNDARYBOND_H

#include <Core/Types/realtypes.h>
#include <Core/Math/Vec.h>
#include <Peridynamics/PeriParticle.h>
#include <InputOutput/Parameter.h>
#include <boost/mpi.hpp>

// the bond between peri-points and the boundaries
// July 14, 2014
namespace dem {

class PeriBoundaryBond {

public:
  PeriBoundaryBond() {
    initBoundaryProjector = 0;
    initBondVec = 0;
    currBondVec = initBondVec;
    isAlive = true;
    periPoint = 0;
  }

  PeriBoundaryBond(const Vec &a, periDynamics::PeriParticle *pt) {
    initBoundaryProjector = a;
    initBondVec = pt->getCurrPosition() - a;
    currBondVec = initBondVec;
    isAlive = true;
    periPoint = pt;
  }

  void applyBondForce(REAL bndry_coord, int bndry_type);
  // bndry_coord is the coordinate of the boundary, int is boundary type

private:
  // bond has two ends, one end is the peri-point,
  // the other is the projector of this peri-point in the boundary,
  // this initBoundaryProjector is the position of this projector point in the
  // initial
  Vec initBoundaryProjector;
  Vec initBondVec; // the initial bond vector pointing from the projector to the
                   // peri-point
  Vec currBondVec; // the current bond vector pointing from the current boundary
                   // point to the peri-point
  bool isAlive; // the state if the bond is alive
  periDynamics::PeriParticle *periPoint;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &initBoundaryProjector;
    ar &initBondVec;
    ar &currBondVec;
    ar &isAlive;
    ar &periPoint;
  }
};

} // end dem

#endif
