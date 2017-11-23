#include <Core/Const/Constants.h>
#include <Core/Math/root6.h>
#include <Core/Math/root6_old.h>
#include <Core/Util/Utility.h>
#include <DiscreteElements/DEMContact.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

//#define MINDLIN_ASSUMED
//#define MINDLIN_KNOWN

//#define USE_ROOT6_OLD

namespace dem {

DEMContact::DEMContact()
{
  d_p1 = nullptr;
  d_p2 = nullptr;
  d_penetration = 0;
  d_contactRadius = 0;
  d_point1 = d_point2 = 0;
  d_radius1 = d_radius2 = 0;
  d_normalDirection = d_tangentDirection = 0;

  d_isInContact = false;
  d_tangentLoadingActive = d_prevTangentLoadingActive = true;
  d_normalForce = d_prevNormalForce = 0;
  d_tangentForce = d_prevTangentForce = 0;
  d_tangentDisplacement = d_prevTangentDisplacement = 0;
  d_tangentDisplacementStart = 0;
  d_tangentSlidingActive = d_prevTangentSlidingActive = false;
  d_tangentForcePeak = 0;

  d_cohesionForce = 0;
  d_spinResist = 0;

  d_E0 = d_G0 = d_R0 = 0;
}

DEMContact::DEMContact(DEMParticle* t1, DEMParticle* t2)
{
  d_p1 = t1;
  d_p2 = t2;
  d_penetration = 0;
  d_contactRadius = 0;
  d_point1 = d_point2 = 0;
  d_radius1 = d_radius2 = 0;
  d_normalDirection = d_tangentDirection = 0;

  d_isInContact = false;
  d_tangentLoadingActive = d_prevTangentLoadingActive = true;
  d_normalForce = d_prevNormalForce = 0;
  d_tangentForce = d_prevTangentForce = 0;
  d_tangentDisplacement = d_prevTangentDisplacement = 0;
  d_tangentDisplacementStart = 0;
  d_tangentSlidingActive = d_prevTangentSlidingActive = false;
  d_tangentForcePeak = 0;

  d_cohesionForce = 0;
  d_spinResist = 0;

  d_E0 = d_G0 = d_R0 = 0;
}

bool
DEMContact::isRedundant(const DEMContact& other) const
{
  std::size_t id1 = getP1()->getId();
  std::size_t id2 = getP2()->getId();
  std::size_t oId1 = (other.getP1())->getId();
  std::size_t oId2 = (other.getP2())->getId();

  return ((id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2));
}

bool
DEMContact::operator==(const DEMContact& other) const
{
  std::size_t id1 = getP1()->getId();
  std::size_t id2 = getP2()->getId();
  std::size_t oId1 = (other.getP1())->getId();
  std::size_t oId2 = (other.getP2())->getId();

  return ((id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2));
}

DEMParticle*
DEMContact::getP1() const
{
  return d_p1;
}

DEMParticle*
DEMContact::getP2() const
{
  return d_p2;
}

bool
DEMContact::isOverlapped(REAL minRelativeOverlap, REAL measurableOverlap)
{
  // v[0] is the point on p2, v[1] is the point on p1
  REAL coef1[10], coef2[10];
  d_p1->getGlobalCoef(coef1); 
  d_p2->getGlobalCoef(coef2);

  d_radius1 = d_p1->computeRadius(d_point1);
  d_radius2 = d_p2->computeRadius(d_point2);
  auto maxRadius = std::max(d_radius1, d_radius2);


  Vec v[2];
#ifdef USE_ROOT6_OLD
  bool b1 = root6_old(coef1, coef2, v[0], d_radius1, d_p1->getId(), d_p2->getId());
  bool b2 = root6_old(coef2, coef1, v[1], d_radius2, d_p2->getId(), d_p1->getId());
#else
  bool b1 = root6(coef1, coef2, v[0], d_radius1, d_p1->getId(), d_p2->getId());
  bool b2 = root6(coef2, coef1, v[1], d_radius2, d_p2->getId(), d_p1->getId());
#endif
  d_point1 = v[1];
  d_point2 = v[0];
  d_penetration = vnormL2(d_point1 - d_point2);
  /*
  if ((d_p1->getId() == 2 && d_p2->getId() == 94) ||
      (d_p1->getId() == 94 && d_p2->getId() == 2)) {
    auto printCoef = [](const auto& coef) {
      std::ostringstream os;
      for (int ii = 0; ii < 10; ii++) {
        os << coef[ii] << " ";
      }
      return os.str();
    };
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout << " MPI_Rank = " << world_rank << " Particles "
    //          << " p1 = " << d_p1->getId() << " p2 = " << d_p2->getId()
    //          << " d_penetration = " << d_penetration << "\n"
    //          << " \t Pos1 = " << d_p1->currentPosition()
    //          << " \t Pos2 = " << d_p2->currentPosition() << "\n"
    //          << " \t b1 = " << b1 << " b2 = " << b2 << "\n"
    //          << " \t coef1 = [" << printCoef(coef1) << "] \n"
    //          << " \t coef2 = [" << printCoef(coef2) << "] \n"
    //          << " \t point1 = " << d_point1 << " point2 = " << d_point2 << "\n"
    //          << " \t point1_old = " << v2_old
    //          << " point2_old = " << v1_old
    //          << " radius1 = " << d_radius1 << " radius2 = " << d_radius2
    //          << "\n";
  }
  */

  if (b1 && b2 &&
      d_penetration / (2 * maxRadius) > minRelativeOverlap &&
      nearbyint(d_penetration / measurableOverlap) >= 1) { 
    // a strict detection method
    d_isInContact = true;
    return true;
  } else {
    d_isInContact = false;
    return false;
  }
}

void
DEMContact::checkinPreviousContactTangents(ContactTangentArray& contactTangentVec)
{
  for (auto& tangent : contactTangentVec) {
    if (tangent.ct_particle1 == d_p1->getId() && 
        tangent.ct_particle2 == d_p2->getId()) {
      d_prevTangentForce = tangent.ct_tangentForce;
      d_prevTangentDisplacement = tangent.ct_tangentDisplacement;
      d_prevTangentLoadingActive = tangent.ct_tangentLoadingActive;
      d_tangentDisplacementStart = tangent.ct_tangentDispStart;
      d_tangentForcePeak = tangent.ct_tangentForcePeak;
      d_prevTangentSlidingActive = tangent.ct_tangentSlidingActive;
      /*
      if ((it.d_particle1 == 2 && it.d_particle2 == 94) ||
          (it.d_particle1 == 94 && it.d_particle2 == 2)) {
        //std::cout << "DEMContact tangents: " << std::setprecision(16)
        //          << " TangentForce =  " << d_prevTangentForce
        //          << " TangentDisp =  " << d_prevTangentDisplacement << "\n";
      }
      */
      break;
    }
  }
}

void
DEMContact::checkoutContactTangents(ContactTangentArray& contactTangentVec)
{
  contactTangentVec.push_back(DEMContactTangent(d_p1->getId(), d_p2->getId(), d_tangentForce,
                                     d_tangentDisplacement, d_tangentLoadingActive, d_tangentDisplacementStart,
                                     d_tangentForcePeak, d_tangentSlidingActive));
}

// isOverlapped() has been called in findContact() in dem.cpp and
// information recorded,
// now this function is called by internalForce() in dem.cpp.
void
DEMContact::computeContactForces(std::size_t iteration)
{
  if (!d_isInContact) {
    d_isInContact = false;
    d_tangentLoadingActive = false;
    d_tangentForcePeak = 0;
    d_normalForce = 0;
    d_tangentForce = 0;
    d_tangentDisplacement = 0; // total value
    d_normalDirection = 0;
    d_tangentDirection = 0;

    d_penetration = 0;
    d_contactRadius = 0;
    d_radius1 = d_radius2 = 0;
    d_spinResist = 0;
    return;
  }

  /*
  if ((d_p1->getId() == 2 && d_p2->getId() == 94) ||
      (d_p1->getId() == 94 && d_p2->getId() == 2))  {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout << "Before DEMContact: MPI_rank = " << world_rank 
    //          << " iter = " << iteration << std::setprecision(16)
    //          << " id = " << d_p1->getId() << " " << d_p2->getId() << "\n\t"
    //          << " f1 = " << d_p1->force() << "\n\t"
    //          << " f2 = " << d_p2->force() << "\n\t"
    //          << " d_penetration=" << d_penetration << "\n\t" << "\n";
  }
  */

  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");
  REAL maxAllowableRelativeOverlap = util::getParam<REAL>("maxAllowableRelativeOverlap");
  REAL minMeasurableOverlap = util::getParam<REAL>("minMeasurableOverlap");
  REAL contactCohesion = util::getParam<REAL>("contactCohesion");
  REAL contactDamp = util::getParam<REAL>("contactDamp");
  REAL contactFric = util::getParam<REAL>("contactFric");

  // obtain normal force, using absolute equation instead of stiffness method
  d_p1->setContactNum(d_p1->getNumBoundaryContacts() + 1);
  d_p2->setContactNum(d_p2->getNumBoundaryContacts() + 1);
  d_p1->setInContact(true);
  d_p2->setInContact(true);

  d_R0 = d_radius1 * d_radius2 / (d_radius1 + d_radius2);
  d_E0 = 0.5 * young / (1 - poisson * poisson);
  REAL allowedOverlap = 2.0 * fmin(d_radius1, d_radius2) * maxAllowableRelativeOverlap;
  if (d_penetration > allowedOverlap) {
    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);
    inf << " DEMContact.cpp: iter=" << std::setw(8) << iteration
        << " ptcl1=" << std::setw(8) << getP1()->getId()
        << " ptcl2=" << std::setw(8) << getP2()->getId()
        << " d_penetration=" << std::setw(OWID) << d_penetration
        << " allow=" << std::setw(OWID) << allowedOverlap << std::endl;
    MPI_Status status;
    //int length = OWID * 2 + 8 * 3 + 19 + 7 * 3 + 8 + 1;
    MPI_File_write_shared(overlapInf, const_cast<char*>(inf.str().c_str()),
                          inf.str().length(), MPI_CHAR, &status);

    d_penetration = allowedOverlap;
  }

  d_penetration = nearbyint(d_penetration / minMeasurableOverlap) * minMeasurableOverlap;
  d_contactRadius = sqrt(d_penetration * d_R0);
  // d_normalDirection points from particle 1 to particle 2
  d_normalDirection = normalize(d_point1 - d_point2);
  // normalForce pointing to particle 1 // pow(d_penetration, 1.5)
  d_normalForce = -sqrt(d_penetration * d_penetration * d_penetration) * sqrt(d_R0) * 4 *
                  d_E0 / 3 * d_normalDirection;

  // apply cohesion force
  d_cohesionForce = Pi * (d_penetration * d_R0) * contactCohesion * d_normalDirection;

  // obtain normal damping force
  Vec cp = (d_point1 + d_point2) / 2;
  Vec veloc1 =
    d_p1->currentVelocity() + cross(d_p1->currentAngularVelocity(), (cp - d_p1->currentPosition()));
  Vec veloc2 =
    d_p2->currentVelocity() + cross(d_p2->currentAngularVelocity(), (cp - d_p2->currentPosition()));
  REAL m1 = getP1()->mass();
  REAL m2 = getP2()->mass();
  REAL kn = pow(6 * vnormL2(d_normalForce) * d_R0 * pow(d_E0, 2), 1.0 / 3.0);
  REAL dampCritical = 2 * sqrt(m1 * m2 / (m1 + m2) * kn); // critical damping
  Vec contactDampingForce = contactDamp * dampCritical *
                        dot(veloc1 - veloc2, d_normalDirection) * d_normalDirection;

  d_vibraTimeStep = 2.0 * sqrt(m1 * m2 / (m1 + m2) / kn);
  Vec relativeVel = veloc1 - veloc2;
  d_impactTimeStep =
    (relativeVel.lengthSq() < std::numeric_limits<double>::min())
      ? std::numeric_limits<double>::max()
      : allowedOverlap / fabs(dot(relativeVel, d_normalDirection));

  // obtain tangential force
  if (contactFric != 0) {
    d_G0 = young / 2 / (1 + poisson);
    // RelaDispInc points along point1's displacement relative to point2
    Vec RelaDispInc = (veloc1 - veloc2) * timeStep;
    Vec tangentDispInc = RelaDispInc - dot(RelaDispInc, d_normalDirection) * d_normalDirection;
    // prevTangentDisp read by checkinPreviousContactTangents()
    d_tangentDisplacement = d_prevTangentDisplacement + tangentDispInc;
    if (vnormL2(d_tangentDisplacement) == 0) {
      d_tangentDirection = 0;
    } else {
      // tangentDirc points along Tangentential forces exerted on particle 1
      d_tangentDirection = normalize(-d_tangentDisplacement);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // linear friction model
    REAL fP = contactFric * vnormL2(d_normalForce);
    REAL ks = 4 * d_G0 * d_contactRadius / (2 - poisson);
    // d_prevTangentForce read by CheckinPreTangent()
    d_tangentForce = d_prevTangentForce + ks * (-tangentDispInc);
    if (vnormL2(d_tangentForce) > fP) {
      d_tangentForce = fP * d_tangentDirection;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MINDLIN_ASSUMED
    computeTangentForceMindlinAssumed(contactFric, poisson, tangentDispInc);
#endif
#ifdef MINDLIN_KNOWN
    computeTangentForceMindlinKnown(contactFric, poisson, tangentDispInc);
#endif
  }

  // apply forces
  Vec totalForce = d_normalForce +  d_tangentForce + d_cohesionForce
                   - contactDampingForce;
  Vec momentArm = cp - d_p1->currentPosition();
  Vec totalMoment = cross(momentArm, (d_normalForce + d_tangentForce - contactDampingForce));

  // Update the forces and moments
  //d_p1->addForce(totalForce);
  //d_p2->addForce(-totalForce);
  //d_p1->addMoment(totalMoment);
  //d_p2->addMoment(-totalMoment);
  d_p1->addForceIDMap(totalForce, d_p2->getId());
  d_p2->addForceIDMap(-totalForce, d_p1->getId());
  d_p1->addMomentIDMap(totalMoment, d_p2->getId());
  d_p2->addMomentIDMap(-totalMoment, d_p1->getId());

  /*
  if ((d_p1->getId() == 2 && d_p2->getId() == 94) ||
      (d_p1->getId() == 94 && d_p2->getId() == 2))  {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout << "After DEMContact: MPI_rank = " << world_rank 
    //          << " iter = " << iteration << std::setprecision(16)
    //          << " id = " << d_p1->getId() << " " << d_p2->getId() << "\n\t"
    //          << " f1 = " << d_p1->force() << "\n\t"
    //          << " f2 = " << d_p2->force() << "\n\t"
    //          << " d_penetration=" << d_penetration << "\n\t" << "\n";
  }
  */
  /*
  {
    //std::cout << " cohesionForce=" << d_cohesionForce << "\n\t"
              << " normalForce=" << d_normalForce << "\n\t"
              << " tangentForce =" << d_tangentForce << "\n\t"
              << " dampingForce = " << contactDampingForce << "\n\t"
              << " totalForce = " << totalForce << "\n\t"
              << " totalMoment = " << totalMoment << "\n\t"
              << " accumulated time=" << iteration * timeStep << "\n"
              << "\t cp = " << cp << "\n\t"
              << " d_point1 = " << d_point1 << " d_point2 = " << d_point2
              << "\n\t"
              << " V1 = " << d_p1->currentVelocity()
              << " V2 = " << d_p2->currentVelocity() << "\n\t"
              << " W1 = " << d_p1->currentAngularVelocity()
              << " W2 = " << d_p2->currentAngularVelocity() << "\n\t"
              << " P1 = " << d_p1->currentPosition()
              << " P2 = " << d_p2->currentPosition() << "\n\t"
              << " veloc1 = " << veloc1 << " veloc2 = " << veloc2 << "\n";
  }
  */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mindlin's model (loading/unloading condition assumed)
// This model is not recommended as it is impossible to strictly determine
// loading/unloading condition
// unless load is known (the case of pure moment rotation).
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
DEMContact::computeTangentForceMindlinAssumed(const REAL& contactFric,
                                           const REAL& poisson,
                                           const Vec& tangentDispInc)
{
  REAL val = 0, ks = 0;
  REAL fP = contactFric * vnormL2(d_normalForce);
  d_tangentLoadingActive = (dot(d_prevTangentDisplacement, tangentDispInc) >= 0);

  if (d_tangentLoadingActive) {        // loading
    if (!d_prevTangentLoadingActive) { // pre-step is unloading
      val = 8 * d_G0 * d_contactRadius * vnormL2(tangentDispInc) /
            (3 * (2 - poisson) * fP);
      d_tangentDisplacementStart = d_prevTangentDisplacement;
    } else // pre-step is loading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tangentDisplacement - d_tangentDisplacementStart) /
            (3 * (2 - poisson) * fP);

    if (val > 1.0)
      d_tangentForce = fP * d_tangentDirection;
    else {
      ks = 4 * d_G0 * d_contactRadius / (2 - poisson) * sqrt(1 - val);
      // incremental method
      d_tangentForce =
        d_prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
      // total value method: tangentForce = fP*(1-pow(1-val, 1.5))*d_tangentDirection;
    }
  } else {                  // unloading
    if (d_prevTangentLoadingActive) { // pre-step is loading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tangentDisplacement - d_tangentDisplacementStart) /
            (3 * (2 - poisson) * fP);
      d_tangentForcePeak = vnormL2(d_prevTangentForce);
    } else // pre-step is unloading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tangentDisplacement - d_tangentDisplacementStart) /
            (3 * (2 - poisson) * fP);

    if (val > 1.0 || d_tangentForcePeak > fP)
      d_tangentForce = fP * d_tangentDirection;
    else {
      ks = 2 * sqrt(2) * d_G0 * d_contactRadius / (2 - poisson) *
           sqrt(1 + pow(1 - d_tangentForcePeak / fP, 2.0 / 3.0) + val);
      // incremental method
      d_tangentForce =
        d_prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
      // total value method: d_tangentForce = (d_tangentForcePeak-2*fP*(1-sqrt(2)/4*pow(1+
      // pow(1-d_tangentForcePeak/fP,2.0/3.0) + val,1.5)))*d_tangentDirection;
    }
  }

  if (vnormL2(d_tangentForce) > fP)
    d_tangentForce = fP * d_tangentDirection;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mindlin's model (loading/unloading condition known for pure moment rotation
// case)
// As loading/unloading condition is known, both incremental and total value
// method work well.
// Herein sliding history is incorporated.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
DEMContact::computeTangentForceMindlinKnown(const REAL& contactFric,
                                         const REAL& poisson,
                                         const Vec& tangentDispInc)
{
  REAL val = 0, ks = 0;
  REAL fP = contactFric * vnormL2(d_normalForce);
  if (d_prevTangentSlidingActive)
    val =
      8 * d_G0 * d_contactRadius * vnormL2(tangentDispInc) / (3 * (2 - poisson) * fP);
  else
    val = 8 * d_G0 * d_contactRadius * vnormL2(d_tangentDisplacement - d_tangentDisplacementStart) /
          (3 * (2 - poisson) * fP);

  if (iteration > 10000 &&
      iteration < 11000) { // loading (and possible sliding)
    if (val > 1.0) {
      d_tangentForce = fP * d_tangentDirection;
      d_tangentSlidingActive = true;
    } else {
      if (!d_prevTangentSlidingActive) {
        ks = 4 * d_G0 * d_contactRadius / (2 - poisson) * sqrt(1 - val);
        d_tangentForce =
          d_prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        d_tangentSlidingActive = false;
      } else {
        if (vnormL2(d_tangentForce) > vnormL2(d_prevTangentForce))
          d_tangentSlidingActive = true;
        else
          d_tangentSlidingActive = false;
      }
    }
    d_tangentForcePeak = vnormL2(d_tangentForce);
  } else { // (possible sliding and) unloading
    if (val > 1.0 || d_tangentForcePeak > fP) {
      d_tangentForce = fP * d_tangentDirection;
      d_tangentSlidingActive = true;
    } else {
      if (!d_prevTangentSlidingActive) {
        ks = 2 * sqrt(2) * d_G0 * d_contactRadius / (2 - poisson) *
             sqrt(1 + pow(1 - d_tangentForcePeak / fP, 2.0 / 3.0) + val);
        d_tangentForce =
          d_prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        d_tangentSlidingActive = false;
      } else {
        if (vnormL2(d_tangentForce) > vnormL2(d_prevTangentForce))
          d_tangentSlidingActive = true;
        else {
          d_tangentSlidingActive = false;
          d_tangentDisplacementStart = d_tangentDisplacement;
        }
      }
    }
  }

  /*
        //std::cout<< "DEMContact.h: iter="<iteration
                 << " prevTangentSlide=" << d_prevTangentSlidingActive
                 << " tangentSlide=" << d_tangentSlidingActive
                 << " val=" << val
                 << " ks=" << ks
                 << " tangentDispInc.x=" << tangentDispInc.x()
                 << " d_prevTangentForce=" << vnormL2(d_prevTangentForce)
                 << " d_tangentForce" << vnormL2(d_tangentForce)
                 << std::endl;
        */

  if (vnormL2(d_tangentForce) > fP)
    d_tangentForce = fP * d_tangentDirection;
}

} // namespace dem ends
