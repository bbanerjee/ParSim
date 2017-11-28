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
  d_vibrationTimeStep = std::numeric_limits<double>::max();
  d_impactTimeStep = d_vibrationTimeStep;
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

  d_vibrationTimeStep = std::numeric_limits<double>::max();
  d_impactTimeStep = d_vibrationTimeStep;
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
DEMContact::isOverlapped(REAL minRelativeOverlap, REAL measurableOverlap,
                         std::size_t iteration)
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

  d_isInContact = false;

  if (b1 && b2 &&
      d_penetration / (2 * maxRadius) > minRelativeOverlap &&
      nearbyint(d_penetration / measurableOverlap) >= 1) { 
    // a strict detection method
    //if (d_p1->getId() == 284 || d_p2->getId() == 284) {
    //  std::cout << "Iteration = " << iteration
    //            << " p1 = (" << d_p1->getId() << "," 
    //            << static_cast<int>(d_p1->getType()) << ")"
    //            << " p2 = (" << d_p2->getId() << "," 
    //            << static_cast<int>(d_p2->getType()) << ")"
    //            << " b1 = " << std::boolalpha << b1
    //            << " b2 = " << std::boolalpha << b2
    //            << " penetration = " << d_penetration
    //            << " maxRadius = " << maxRadius << "\n";
    //}
    d_isInContact = true;
  } 

  return d_isInContact;
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
int
DEMContact::computeContactForces(REAL timeStep, std::size_t iteration,
                                 REAL stiffness, REAL shearModulus,
                                 REAL cohesion, REAL damping,
                                 REAL friction, REAL maxOverlapFactor,
                                 REAL minMeasurableOverlap)
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
    return 0;
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

  // obtain normal force, using absolute equation instead of stiffness method
  d_p1->setContactNum(d_p1->getNumBoundaryContacts() + 1);
  d_p2->setContactNum(d_p2->getNumBoundaryContacts() + 1);
  d_p1->setInContact(true);
  d_p2->setInContact(true);

  d_R0 = d_radius1 * d_radius2 / (d_radius1 + d_radius2);
  d_E0 = stiffness;
  auto poissonRatio = 1 - shearModulus/stiffness;

  REAL allowedOverlap = 2.0 * std::min(d_radius1, d_radius2) * maxOverlapFactor;
  int excessiveOverlap = 0;
  if (d_penetration > allowedOverlap) {
    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);
    inf << " DEMContact.cpp: " << std::setw(4) << __LINE__ << ":"
        << " iter=" << std::setw(8) << iteration
        << " ptcl1=(" << std::setw(8) << getP1()->getId()
        << ", " << std::setw(2) << static_cast<int>(getP1()->getType()) << ")"
        << " ptcl2=(" << std::setw(8) << getP2()->getId()
        << ", " << std::setw(2) << static_cast<int>(getP2()->getType()) << ")"
        << " d_penetration=" << std::setw(OWID) << d_penetration
        << " allow=" << std::setw(OWID) << allowedOverlap << std::endl;
    MPI_Status status;
    //int length = OWID * 2 + 8 * 3 + 19 + 7 * 3 + 8 + 1;
    MPI_File_write_shared(overlapInf, const_cast<char*>(inf.str().c_str()),
                          inf.str().length(), MPI_CHAR, &status);

    d_penetration = allowedOverlap;
    excessiveOverlap = 1;
  }

  d_penetration = 
    nearbyint(d_penetration / minMeasurableOverlap) * minMeasurableOverlap;

  d_contactRadius = std::sqrt(d_penetration * d_R0);
  // d_normalDirection points from particle 1 to particle 2
  d_normalDirection = normalize(d_point1 - d_point2);
  // normalForce pointing to particle 1 // pow(d_penetration, 1.5)
  d_normalForce = -std::sqrt(d_penetration * d_penetration * d_penetration) * 
                   std::sqrt(d_R0) * 4 * d_E0 / 3 * d_normalDirection;
  auto kn = computeNormalStiffness(d_normalForce, d_R0, d_E0);

  d_cohesionForce = Pi * (d_penetration * d_R0) * cohesion * d_normalDirection;

  REAL m1 = getP1()->mass();
  REAL m2 = getP2()->mass();

  Vec midPoint = (d_point1 + d_point2) / 2;
  Vec pVelocity1 = computeVelocity(d_p1, midPoint);
  Vec pVelocity2 = computeVelocity(d_p2, midPoint);
  Vec relativeVel = pVelocity1 - pVelocity2;

  auto contactDampingForce = computeDampingForce(m1, m2, relativeVel,
                                                 damping, kn, d_normalDirection);

  updateTimestep(m1, m2, relativeVel, kn, d_normalDirection, allowedOverlap);

  if (friction != 0) {
    computeTangentForce(timeStep, iteration, shearModulus, poissonRatio,
                        friction, d_contactRadius, relativeVel,
                        d_normalDirection, d_normalForce);
  }

  updateForceAndMoment(d_normalForce, d_tangentForce, d_cohesionForce, 
                       contactDampingForce, midPoint, d_p1, d_p2);

  return excessiveOverlap;
}

// Compute the damping force
Vec
DEMContact::computeDampingForce(REAL pMass1, REAL pMass2, 
                                const Vec& relativeVel,
                                REAL dampingCoeff, REAL normalStiffness,
                                const Vec& contactNormal)
{
  // critical damping
  REAL dampCritical = 
    2 * std::sqrt(pMass1 * pMass2 / (pMass1 + pMass2) * normalStiffness); 
  Vec contactDampingForce = dampingCoeff * dampCritical *
    dot(relativeVel, contactNormal) * contactNormal;
  return contactDampingForce;
}

// Update time step sizes
void
DEMContact::updateTimestep(REAL pMass1, REAL pMass2, 
                           const Vec& relativeVel, 
                           REAL normalStiffness, const Vec& contactNormal,
                           REAL allowedOverlap)
{
  d_vibrationTimeStep = 
    2 * std::sqrt(pMass1 * pMass2 / ((pMass1 + pMass2) * normalStiffness));
  d_impactTimeStep =
    (relativeVel.lengthSq() < std::numeric_limits<double>::min())
      ? std::numeric_limits<double>::max()
      : allowedOverlap / std::abs(dot(relativeVel, contactNormal));
  //std::cout << "normalStiffness = " << normalStiffness
  //          << " t_v = " << d_vibrationTimeStep
  //          << " t_i = " << d_impactTimeStep << "\n";
}

// Compute the tangential force
void
DEMContact::computeTangentForce(REAL timeStep, std::size_t iteration,
                                REAL shearModulus, 
                                REAL poissonRatio,
                                REAL friction, 
                                REAL contactRadius,
                                const Vec& relativeVel, 
                                const Vec& contactNormal, 
                                const Vec& normalForce)
{
  // disp points along point1's displacement relative to point2
  Vec dispInc = relativeVel * timeStep;
  Vec tangentDispInc = dispInc - dot(dispInc, contactNormal) * contactNormal;

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
  REAL fp = friction * vnormL2(normalForce);
  REAL ks = 4 * shearModulus * contactRadius / (2 - poissonRatio);

  // d_prevTangentForce read by CheckinPreTangent()
  d_tangentForce = d_prevTangentForce + ks * (-tangentDispInc);
  if (vnormL2(d_tangentForce) > fp) {
    d_tangentForce = fp * d_tangentDirection;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  #ifdef MINDLIN_ASSUMED
    computeTangentForceMindlinAssumed(friction, poissonRatio, tangentDispInc);
  #endif
  #ifdef MINDLIN_KNOWN
    computeTangentForceMindlinKnown(friction, poissonRatio, tangentDispInc, iteration);
  #endif
}

// Update the forces and moments
void
DEMContact::updateForceAndMoment(const Vec& normalForce, const Vec& tangentForce,
                                 const Vec& cohesionForce, const Vec& dampingForce,
                                 const Vec& momentCenter,
                                 DEMParticle* p1, DEMParticle* p2) 
{
  Vec totalForce = normalForce + tangentForce + cohesionForce - dampingForce;
  Vec totalMoment = cross((momentCenter - p1->currentPosition()), 
                          (normalForce + tangentForce - dampingForce));

  auto particle1 = p1->getId();
  auto particle2 = p2->getId();
  p1->addForceIDMap(totalForce, particle2);
  p2->addForceIDMap(-totalForce, particle1);
  p1->addMomentIDMap(totalMoment, particle2);
  p2->addMomentIDMap(-totalMoment, particle1);

  //std::cout << " cohesionForce=" << d_cohesionForce << "\n\t"
  //          << " normalForce=" << d_normalForce << "\n\t"
  //          << " tangentForce =" << d_tangentForce << "\n\t"
  //          << " dampingForce = " << contactDampingForce << "\n\t"
  //          << " totalForce = " << totalForce << "\n\t"
  //          << " totalMoment = " << totalMoment << "\n\t"
  //          << "\t cp = " << cp << "\n\t";
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
                                            const Vec& tangentDispInc,
                                            std::size_t iteration)
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
