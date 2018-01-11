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
}

DEMContact::DEMContact(DEMParticle* t1, DEMParticle* t2)
{
  d_p1 = t1;
  d_p2 = t2;
}

DEMContact::DEMContact(DEMParticle* t1, DEMParticle* t2,
                       const DEMContactData& data)
{
  d_p1 = t1;
  d_p2 = t2;
  d_data = data;
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

  d_data.radius1 = d_p1->computeRadius(d_data.point1);
  d_data.radius2 = d_p2->computeRadius(d_data.point2);
  auto maxRadius = std::max(d_data.radius1, d_data.radius2);

  Vec v[2];
#ifdef USE_ROOT6_OLD
  bool b1 = root6_old(coef1, coef2, v[0], d_data.radius1, d_p1->getId(), d_p2->getId());
  bool b2 = root6_old(coef2, coef1, v[1], d_data.radius2, d_p2->getId(), d_p1->getId());
#else
  bool b1 = root6(coef1, coef2, v[0], d_data.radius1, d_p1->getId(), d_p2->getId());
  bool b2 = root6(coef2, coef1, v[1], d_data.radius2, d_p2->getId(), d_p1->getId());
#endif
  d_data.point1 = v[1];
  d_data.point2 = v[0];
  d_data.penetration = vnormL2(d_data.point1 - d_data.point2);
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
    //          << " d_penetration = " << d_data.penetration << "\n"
    //          << " \t Pos1 = " << d_p1->currentPosition()
    //          << " \t Pos2 = " << d_p2->currentPosition() << "\n"
    //          << " \t b1 = " << b1 << " b2 = " << b2 << "\n"
    //          << " \t coef1 = [" << printCoef(coef1) << "] \n"
    //          << " \t coef2 = [" << printCoef(coef2) << "] \n"
    //          << " \t point1 = " << d_data.point1 << " point2 = " << d_data.point2 << "\n"
    //          << " \t point1_old = " << v2_old
    //          << " point2_old = " << v1_old
    //          << " radius1 = " << d_data.radius1 << " radius2 = " << d_data.radius2
    //          << "\n";
  }
  */

  d_data.isInContact = false;

  if (b1 && b2 &&
      d_data.penetration / (2 * maxRadius) > minRelativeOverlap &&
      nearbyint(d_data.penetration / measurableOverlap) >= 1) { 
    // a strict detection method
    //if (d_p1->getId() == 284 || d_p2->getId() == 284) {
    //  std::cout << "Iteration = " << iteration
    //            << " p1 = (" << d_p1->getId() << "," 
    //            << static_cast<int>(d_p1->getType()) << ")"
    //            << " p2 = (" << d_p2->getId() << "," 
    //            << static_cast<int>(d_p2->getType()) << ")"
    //            << " b1 = " << std::boolalpha << b1
    //            << " b2 = " << std::boolalpha << b2
    //            << " penetration = " << d_data.penetration
    //            << " maxRadius = " << maxRadius << "\n";
    //}
    d_data.isInContact = true;
  } 

  return d_data.isInContact;
}

void
DEMContact::checkinPreviousContactTangents(DEMContactTangentArray& contactTangentVec)
{
  for (auto& tangent : contactTangentVec) {
    if (tangent.ct_particle1 == d_p1->getId() && 
        tangent.ct_particle2 == d_p2->getId()) {
      d_data.prevTangentForce = tangent.ct_tangentForce;
      d_data.prevTangentDisplacement = tangent.ct_tangentDisplacement;
      d_data.prevTangentLoadingActive = tangent.ct_tangentLoadingActive;
      d_data.tangentDisplacementStart = tangent.ct_tangentDispStart;
      d_data.tangentForcePeak = tangent.ct_tangentForcePeak;
      d_data.prevTangentSlidingActive = tangent.ct_tangentSlidingActive;
      /*
      if ((it.d_particle1 == 2 && it.d_particle2 == 94) ||
          (it.d_particle1 == 94 && it.d_particle2 == 2)) {
        //std::cout << "DEMContact tangents: " << std::setprecision(16)
        //          << " TangentForce =  " << d_data.prevTangentForce
        //          << " TangentDisp =  " << d_data.prevTangentDisplacement << "\n";
      }
      */
      break;
    }
  }
}

void
DEMContact::checkoutContactTangents(DEMContactTangentArray& contactTangentVec)
{
  contactTangentVec.push_back(DEMContactTangent(d_p1->getId(), d_p2->getId(), 
    d_data.tangentForce, d_data.tangentDisplacement, d_data.tangentLoadingActive, 
    d_data.tangentDisplacementStart, d_data.tangentForcePeak, 
    d_data.tangentSlidingActive));
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
  if (!d_data.isInContact) {
    d_data.isInContact = false;
    d_data.tangentLoadingActive = false;
    d_data.tangentForcePeak = 0;
    d_data.normalForce = 0;
    d_data.tangentForce = 0;
    d_data.tangentDisplacement = 0; // total value
    d_data.normalDirection = 0;
    d_data.tangentDirection = 0;

    d_data.penetration = 0;
    d_data.contactRadius = 0;
    d_data.radius1 = d_data.radius2 = 0;
    d_data.spinResist = 0;
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
    //          << " d_penetration=" << d_data.penetration << "\n\t" << "\n";
  }
  */

  // obtain normal force, using absolute equation instead of stiffness method
  d_p1->setContactNum(d_p1->getNumBoundaryContacts() + 1);
  d_p2->setContactNum(d_p2->getNumBoundaryContacts() + 1);
  d_p1->setInContact(true);
  d_p2->setInContact(true);

  d_data.R0 = d_data.radius1 * d_data.radius2 / (d_data.radius1 + d_data.radius2);
  d_data.E0 = stiffness;
  auto poissonRatio = 1 - shearModulus/stiffness;

  REAL allowedOverlap = 2.0 * std::min(d_data.radius1, d_data.radius2) * maxOverlapFactor;
  int excessiveOverlap = 0;
  if (d_data.penetration > allowedOverlap) {
    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);
    inf << " DEMContact.cpp: " << std::setw(4) << __LINE__ << ":"
        << " iter=" << std::setw(8) << iteration
        << " ptcl1=(" << std::setw(8) << getP1()->getId()
        << ", " << std::setw(2) << static_cast<int>(getP1()->getType()) << ")"
        << " ptcl2=(" << std::setw(8) << getP2()->getId()
        << ", " << std::setw(2) << static_cast<int>(getP2()->getType()) << ")"
        << " d_penetration=" << std::setw(OWID) << d_data.penetration
        << " allow=" << std::setw(OWID) << allowedOverlap << std::endl;
    MPI_Status status;
    //int length = OWID * 2 + 8 * 3 + 19 + 7 * 3 + 8 + 1;
    MPI_File_write_shared(overlapInf, const_cast<char*>(inf.str().c_str()),
                          inf.str().length(), MPI_CHAR, &status);

    //d_data.penetration = allowedOverlap;
    excessiveOverlap = 1;
  }

  //d_data.penetration = 
  //  nearbyint(d_data.penetration / minMeasurableOverlap) * minMeasurableOverlap;

  d_data.contactRadius = std::sqrt(d_data.penetration * d_data.R0);
  // d_normalDirection points from particle 1 to particle 2
  d_data.normalDirection = normalize(d_data.point1 - d_data.point2);
  // normalForce pointing to particle 1 // pow(d_data.penetration, 1.5)
  d_data.normalForce = -std::sqrt(d_data.penetration * d_data.penetration * d_data.penetration) * 
                   std::sqrt(d_data.R0) * 4 * d_data.E0 / 3 * d_data.normalDirection;
  auto kn = computeNormalStiffness(d_data.normalForce, d_data.R0, d_data.E0);

  d_data.cohesionForce = Pi * (d_data.penetration * d_data.R0) * cohesion * d_data.normalDirection;

  REAL m1 = getP1()->mass();
  REAL m2 = getP2()->mass();

  Vec midPoint = (d_data.point1 + d_data.point2) / 2;
  Vec pVelocity1 = computeVelocity(d_p1, midPoint);
  Vec pVelocity2 = computeVelocity(d_p2, midPoint);
  Vec relativeVel = pVelocity1 - pVelocity2;

  auto contactDampingForce = computeDampingForce(m1, m2, relativeVel,
                                                 damping, kn, d_data.normalDirection);

  updateTimestep(m1, m2, relativeVel, kn, d_data.normalDirection, allowedOverlap);

  if (friction != 0) {
    computeTangentForce(timeStep, iteration, shearModulus, poissonRatio,
                        friction, d_data.contactRadius, relativeVel,
                        d_data.normalDirection, d_data.normalForce);
  }

  updateForceAndMoment(d_data.normalForce, d_data.tangentForce, d_data.cohesionForce, 
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
  d_data.vibrationTimeStep = 
    2 * std::sqrt(pMass1 * pMass2 / ((pMass1 + pMass2) * normalStiffness));
  d_data.impactTimeStep =
    (relativeVel.lengthSq() < std::numeric_limits<double>::min())
      ? std::numeric_limits<double>::max()
      : allowedOverlap / std::abs(dot(relativeVel, contactNormal));
  //std::cout << "normalStiffness = " << normalStiffness
  //          << " t_v = " << d_data.vibrationTimeStep
  //          << " t_i = " << d_data.impactTimeStep << "\n";
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
  d_data.tangentDisplacement = d_data.prevTangentDisplacement + tangentDispInc;
  if (vnormL2(d_data.tangentDisplacement) == 0) {
    d_data.tangentDirection = 0;
  } else {
    // tangentDirc points along Tangentential forces exerted on particle 1
    d_data.tangentDirection = normalize(-d_data.tangentDisplacement);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // linear friction model
  REAL fp = friction * vnormL2(normalForce);
  REAL ks = 4 * shearModulus * contactRadius / (2 - poissonRatio);

  // d_prevTangentForce read by CheckinPreTangent()
  d_data.tangentForce = d_data.prevTangentForce + ks * (-tangentDispInc);
  if (vnormL2(d_data.tangentForce) > fp) {
    d_data.tangentForce = fp * d_data.tangentDirection;
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

  auto particle2 = p2->getId();
  p1->addForceIDMap(totalForce, particle2);
  p1->addMomentIDMap(totalMoment, particle2);

  auto particle1 = p1->getId();
  p2->addForceIDMap(-totalForce, particle1);
  p2->addMomentIDMap(-totalMoment, particle1);

  //std::cout << " cohesionForce=" << d_data.cohesionForce << "\n\t"
  //          << " normalForce=" << d_data.normalForce << "\n\t"
  //          << " tangentForce =" << d_data.tangentForce << "\n\t"
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
  REAL fP = contactFric * vnormL2(d_data.normalForce);
  d_data.tangentLoadingActive = (dot(d_data.prevTangentDisplacement, tangentDispInc) >= 0);

  if (d_data.tangentLoadingActive) {        // loading
    if (!d_data.prevTangentLoadingActive) { // pre-step is unloading
      val = 8 * d_data.G0 * d_data.contactRadius * vnormL2(tangentDispInc) /
        (3 * (2 - poisson) * fP);
      d_data.tangentDisplacementStart = d_data.prevTangentDisplacement;
    } else { // pre-step is loading
      val = 8 * d_data.G0 * d_data.contactRadius * 
        vnormL2(d_data.tangentDisplacement - d_data.tangentDisplacementStart) /
        (3 * (2 - poisson) * fP);
    }

    if (val > 1.0)
      d_data.tangentForce = fP * d_data.tangentDirection;
    else {
      ks = 4 * d_data.G0 * d_data.contactRadius / (2 - poisson) * sqrt(1 - val);
      // incremental method
      // tangentDispInc determines signs
      d_data.tangentForce = d_data.prevTangentForce + ks * (-tangentDispInc); 
      // total value method: tangentForce = fP*(1-pow(1-val, 1.5))*d_data.tangentDirection;
    }
  } else {                  // unloading
    if (d_data.prevTangentLoadingActive) { // pre-step is loading
      val = 8 * d_data.G0 * d_data.contactRadius * 
        vnormL2(d_data.tangentDisplacement - d_data.tangentDisplacementStart) /
        (3 * (2 - poisson) * fP);
      d_data.tangentForcePeak = vnormL2(d_data.prevTangentForce);
    } else // pre-step is unloading
      val = 8 * d_data.G0 * d_data.contactRadius * 
        vnormL2(d_data.tangentDisplacement - d_data.tangentDisplacementStart) /
        (3 * (2 - poisson) * fP);

    if (val > 1.0 || d_data.tangentForcePeak > fP)
      d_data.tangentForce = fP * d_data.tangentDirection;
    else {
      ks = 2 * sqrt(2) * d_data.G0 * d_data.contactRadius / (2 - poisson) *
           sqrt(1 + pow(1 - d_data.tangentForcePeak / fP, 2.0 / 3.0) + val);
      // incremental method
      d_data.tangentForce =
        d_data.prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
      // total value method: d_data.tangentForce = (d_data.tangentForcePeak-2*fP*(1-sqrt(2)/4*pow(1+
      // pow(1-d_data.tangentForcePeak/fP,2.0/3.0) + val,1.5)))*d_data.tangentDirection;
    }
  }

  if (vnormL2(d_data.tangentForce) > fP)
    d_data.tangentForce = fP * d_data.tangentDirection;
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
  REAL fP = contactFric * vnormL2(d_data.normalForce);
  if (d_data.prevTangentSlidingActive) {
    val = 8 * d_data.G0 * d_data.contactRadius * 
      vnormL2(tangentDispInc) / (3 * (2 - poisson) * fP);
  } else {
    val = 8 * d_data.G0 * d_data.contactRadius * 
      vnormL2(d_data.tangentDisplacement - d_data.tangentDisplacementStart) /
      (3 * (2 - poisson) * fP);
  }

  if (iteration > 10000 &&
      iteration < 11000) { // loading (and possible sliding)
    if (val > 1.0) {
      d_data.tangentForce = fP * d_data.tangentDirection;
      d_data.tangentSlidingActive = true;
    } else {
      if (!d_data.prevTangentSlidingActive) {
        ks = 4 * d_data.G0 * d_data.contactRadius / (2 - poisson) * sqrt(1 - val);
        d_data.tangentForce =
          d_data.prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        d_data.tangentSlidingActive = false;
      } else {
        if (vnormL2(d_data.tangentForce) > vnormL2(d_data.prevTangentForce))
          d_data.tangentSlidingActive = true;
        else
          d_data.tangentSlidingActive = false;
      }
    }
    d_data.tangentForcePeak = vnormL2(d_data.tangentForce);
  } else { // (possible sliding and) unloading
    if (val > 1.0 || d_data.tangentForcePeak > fP) {
      d_data.tangentForce = fP * d_data.tangentDirection;
      d_data.tangentSlidingActive = true;
    } else {
      if (!d_data.prevTangentSlidingActive) {
        ks = 2 * sqrt(2) * d_data.G0 * d_data.contactRadius / (2 - poisson) *
             sqrt(1 + pow(1 - d_data.tangentForcePeak / fP, 2.0 / 3.0) + val);
        d_data.tangentForce =
          d_data.prevTangentForce + ks * (-tangentDispInc); // tangentDispInc determines signs
        d_data.tangentSlidingActive = false;
      } else {
        if (vnormL2(d_data.tangentForce) > vnormL2(d_data.prevTangentForce))
          d_data.tangentSlidingActive = true;
        else {
          d_data.tangentSlidingActive = false;
          d_data.tangentDisplacementStart = d_data.tangentDisplacement;
        }
      }
    }
  }

  /*
        //std::cout<< "DEMContact.h: iter="<iteration
                 << " prevTangentSlide=" << d_data.prevTangentSlidingActive
                 << " tangentSlide=" << d_data.tangentSlidingActive
                 << " val=" << val
                 << " ks=" << ks
                 << " tangentDispInc.x=" << tangentDispInc.x()
                 << " d_data.prevTangentForce=" << vnormL2(d_data.prevTangentForce)
                 << " d_data.tangentForce" << vnormL2(d_data.tangentForce)
                 << std::endl;
        */

  if (vnormL2(d_data.tangentForce) > fP)
    d_data.tangentForce = fP * d_data.tangentDirection;
}

void 
DEMContact::write(std::stringstream& dataStream) const
{
  auto point1 = getPoint1();
  auto point2 = getPoint2();
  auto center = (point1 + point2)/2;
  dataStream << std::setw(OWID) << getP1()->getId() 
             << std::setw(OWID) << getP2()->getId() 
             << std::setw(OWID) << point1.x() 
             << std::setw(OWID) << point1.y() 
             << std::setw(OWID) << point1.z() 
             << std::setw(OWID) << point2.x() 
             << std::setw(OWID) << point2.y() 
             << std::setw(OWID) << point2.z() 
             << std::setw(OWID) << radius1()
             << std::setw(OWID) << radius2() 
             << std::setw(OWID) << getPenetration() 
             << std::setw(OWID) << getTangentDisplacement()
             << std::setw(OWID) << getContactRadius() 
             << std::setw(OWID) << getR0() 
             << std::setw(OWID) << getE0() 
             << std::setw(OWID) << getNormalForceMagnitude() 
             << std::setw(OWID) << getTangentForceMagnitude()
             << std::setw(OWID) << center.x() 
             << std::setw(OWID) << center.y() 
             << std::setw(OWID) << center.z() 
             << std::setw(OWID) << getNormalForce().x() 
             << std::setw(OWID) << getNormalForce().y() 
             << std::setw(OWID) << getNormalForce().z()
             << std::setw(OWID) << getTangentForce().x() 
             << std::setw(OWID) << getTangentForce().y() 
             << std::setw(OWID) << getTangentForce().z() 
             << std::setw(OWID) << getVibrationTimeStep() 
             << std::setw(OWID) << getImpactTimeStep() << std::endl;
}

void 
DEMContact::write(zen::XmlOut& xml) const
{
  zen::XmlElement& root = xml.ref();
  zen::XmlElement& child = root.addChild("contact");
  zen::XmlOut xml_child(child);
  xml_child.attribute("particle1", getP1()->getId());
  xml_child.attribute("particle2", getP2()->getId());

  d_data.write(xml_child);
}

DEMContactData::DEMContactData()
{
  penetration = 0;
  contactRadius = 0;
  point1 = point2 = 0;
  radius1 = radius2 = 0;
  normalDirection = tangentDirection = 0;

  isInContact = false;
  tangentLoadingActive = prevTangentLoadingActive = true;
  normalForce = prevNormalForce = 0;
  tangentForce = prevTangentForce = 0;
  tangentDisplacement = prevTangentDisplacement = 0;
  tangentDisplacementStart = 0;
  tangentSlidingActive = prevTangentSlidingActive = false;
  tangentForcePeak = 0;

  cohesionForce = 0;
  spinResist = 0;

  E0 = G0 = R0 = 0;
  vibrationTimeStep = std::numeric_limits<double>::max();
  impactTimeStep = vibrationTimeStep;
}

void
DEMContactData::write(zen::XmlOut& xml) const
{
  writeVector("point1", point1, xml);
  writeVector("point2", point2, xml);
  writeVector("normal_force", normalForce, xml);
  writeVector("normal_force_prev", prevNormalForce, xml);
  writeVector("tangent_force", tangentForce, xml);
  writeVector("tangent_force_prev", prevTangentForce, xml);
  writeVector("tangent_displacement", tangentDisplacement, xml);
  writeVector("tangent_displacement_start", tangentDisplacementStart, xml);
  writeVector("tangent_displacement_prev", prevTangentDisplacement, xml);
  writeVector("normal_direction", normalDirection, xml);
  writeVector("tangent_direction", tangentDirection, xml);
  writeVector("cohesion_force", cohesionForce, xml);
  writeVector("spin_resistance", spinResist, xml);

  writeScalar("radius1", radius1, xml);
  writeScalar("radius2", radius2, xml);
  writeScalar("penetration", penetration, xml);
  writeScalar("contact_radius", contactRadius , xml);
  writeScalar("R0", R0 , xml);
  writeScalar("E0", E0 , xml);
  writeScalar("G0", G0 , xml);
  writeScalar("tangent_force_peak", tangentForcePeak, xml);
  writeScalar("vibration_timestep", vibrationTimeStep, xml); 
  writeScalar("impact_timestep", impactTimeStep, xml);

  writeBoolean("is_in_contact", isInContact, xml);
  writeBoolean("tangent_loading_active", tangentLoadingActive, xml);
  writeBoolean("tangent_sliding_active", tangentSlidingActive, xml);
  writeBoolean("tangent_loading_active_prev", prevTangentLoadingActive, xml);
  writeBoolean("tangent_sliding_active_prev", prevTangentSlidingActive, xml);
}

void 
DEMContactData::writeBoolean(const std::string& label,
                             bool variable,
                             zen::XmlOut& xml) const
{
  xml[label](variable);
}

void 
DEMContactData::writeScalar(const std::string& label,
                            double variable,
                            zen::XmlOut& xml) const
{
  xml[label](variable);
}

void 
DEMContactData::writeVector(const std::string& label,
                            const Vec& variable,
                            zen::XmlOut& xml) const
{
  std::ostringstream stream;
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.precision(dem::OPREC);
  stream << "["
        << variable.x() << ", " << variable.y() << ", " << variable.z()
        << "]";
  xml[label](stream.str());
}

void
DEMContactData::read(const zen::XmlIn& xml) 
{
  readValue<Vec>(xml, "point1", point1);
  readValue<Vec>(xml, "point2", point2);
  readValue<Vec>(xml, "normal_force", normalForce);
  readValue<Vec>(xml, "normal_force_prev", prevNormalForce);
  readValue<Vec>(xml, "tangent_force", tangentForce);
  readValue<Vec>(xml, "tangent_force_prev", prevTangentForce);
  readValue<Vec>(xml, "tangent_displacement", tangentDisplacement);
  readValue<Vec>(xml, "tangent_displacement_start", tangentDisplacementStart);
  readValue<Vec>(xml, "tangent_displacement_prev", prevTangentDisplacement);
  readValue<Vec>(xml, "normal_direction", normalDirection);
  readValue<Vec>(xml, "tangent_direction", tangentDirection);
  readValue<Vec>(xml, "cohesion_force", cohesionForce);
  readValue<Vec>(xml, "spin_resistance", spinResist);

  readValue<REAL>(xml, "radius1", radius1);
  readValue<REAL>(xml, "radius2", radius2);
  readValue<REAL>(xml, "penetration", penetration);
  readValue<REAL>(xml, "contact_radius", contactRadius);
  readValue<REAL>(xml, "R0", R0);
  readValue<REAL>(xml, "E0", E0);
  readValue<REAL>(xml, "G0", G0);
  readValue<REAL>(xml, "tangent_force_peak", tangentForcePeak);
  readValue<REAL>(xml, "vibration_timestep", vibrationTimeStep);
  readValue<REAL>(xml, "impact_timestep", impactTimeStep);

  readValue<bool>(xml, "is_in_contact", isInContact);
  readValue<bool>(xml, "tangent_loading_active", tangentLoadingActive);
  readValue<bool>(xml, "tangent_sliding_active", tangentSlidingActive);
  readValue<bool>(xml, "tangent_loading_active_prev", prevTangentLoadingActive);
  readValue<bool>(xml, "tangent_sliding_active_prev", prevTangentSlidingActive);
}

template <typename T>
bool
DEMContactData::readValue(const zen::XmlIn& ps, 
                          const std::string& label,
                          T& value)
{
  // Check the label exists
  auto prop_ps = ps[label];
  if (!prop_ps) {
    std::cerr << "**ERROR** DEMContactData reader: " 
              << label << " not found \n";
    return false;
  }

  // Get the data into a string and remove "[" and "]"
  std::string data;
  prop_ps(data);
  data.erase(std::remove(data.begin(), data.end(), '['), data.end());
  data.erase(std::remove(data.begin(), data.end(), ']'), data.end());

  // Convert string to correct type
  //std::cout << "data str = " << data << std::endl;
  value = Ellip3D::IOUtil::convert<T>(data);
  //std::cout << "data value = " << value << std::endl;
  return true;
}

} // namespace dem ends
