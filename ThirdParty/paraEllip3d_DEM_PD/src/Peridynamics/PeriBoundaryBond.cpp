#include <Peridynamics/PeriBoundaryBond.h>
#include <Core/Util/Utility.h>
#include <iostream>
#include <stdlib.h>

namespace dem {

void
PeriBoundaryBond::applyBondForce(REAL bndry_coord, int bndry_type)
{
  if (!isAlive) { // not alive
    return;
  }

  // (1) calculate the current bond vector
  // bndry_coord is the coordinate of the boundary,
  // bndry_type refers the type of the boundary, following the boundary type in
  // the assembly.h
  // 3 ---- back boundary, x_min
  // 1 ---- front boundary, x_max
  // 4 ---- left boundary, y_min
  // 2 ---- right boundary, y_max
  // 6 ---- bottom boundary, z_min
  // 5 ---- top boundary, z_max

  Vec currBoundaryPoint = initBoundaryProjector;
  switch (bndry_type) {
    case 3:
      currBoundaryPoint.setX(bndry_coord);
      break;
    case 1:
      currBoundaryPoint.setX(bndry_coord);
      break;
    case 4:
      currBoundaryPoint.setY(bndry_coord);
      break;
    case 2:
      currBoundaryPoint.setY(bndry_coord);
      break;
    case 6:
      currBoundaryPoint.setZ(bndry_coord);
      break;
    case 5:
      currBoundaryPoint.setZ(bndry_coord);
      break;
    default:
      std::cout << "boundary type should be between 1 to 6..." << std::endl;
      exit(-1);
      break;
  } // end switch

  currBondVec = periPoint->currentPosition() - currBoundaryPoint;

  // (2) check bond if alive
  // at present, use the same criterioin as the peri-bond used in periDynamics
  REAL stretch = (vfabs(currBondVec) - vfabs(initBondVec)) / vfabs(initBondVec);
  if (stretch > util::getParam<REAL>("bondStretchLimit") ||
      stretch <
        -2.0 * util::getParam<REAL>("bondStretchLimit")) {
    isAlive = false;
    return; // do not need to calculate forces
  }

  // (3) calculate bond force and apply bond force to peri-particle and boundary
  REAL kn_periBndry =
    util::getParam<REAL>("periYoung"); // just in value
  REAL kt_periBndry = util::getParam<REAL>("periYoung"); // just for test, July 15, 2014

  Vec bondn =
    currBondVec % initBondVec / vfabs(initBondVec) *
    normalize(initBondVec);        // normal vector of bond w.r.t. initBondVec
  Vec bondt = currBondVec - bondn; // tangent vector of bond w.r.t. initBondVec

  Vec fn =
    (bondn - initBondVec) *
    kn_periBndry; // force is pointing from the projector to the peri-point
  Vec ft =
    bondt *
    kt_periBndry; // force is pointing from the projector to the peri-point

  // apply forces to peri-point
  periPoint->addAccelerationByForce(-fn - ft);
  // apply forces to boundary
  // do not apply at present, since boundary force is only used to adjust the
  // motion of boundaries. July 15, 2014

} // end applyBondForce

} // end dem
