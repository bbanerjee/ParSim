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
  d_penetr = 0;
  d_contactRadius = 0;
  d_point1 = d_point2 = 0;
  d_radius1 = d_radius2 = 0;
  d_normalDirc = d_tgtDirc = 0;

  d_isInContact = false;
  d_tgtLoading = d_prevTgtLoading = true;
  d_normalForce = d_prevNormalForce = 0;
  d_tgtForce = d_prevTgtForce = 0;
  d_tgtDisp = d_prevTgtDisp = 0;
  d_tgtDispStart = 0;
  d_tgtSlide = d_prevTgtSlide = false;
  d_tgtPeak = 0;

  d_cohesionForce = 0;
  d_spinResist = 0;

  d_E0 = d_G0 = d_R0 = 0;
}

DEMContact::DEMContact(DEMParticle* t1, DEMParticle* t2)
{
  d_p1 = t1;
  d_p2 = t2;
  d_penetr = 0;
  d_contactRadius = 0;
  d_point1 = d_point2 = 0;
  d_radius1 = d_radius2 = 0;
  d_normalDirc = d_tgtDirc = 0;

  d_isInContact = false;
  d_tgtLoading = d_prevTgtLoading = true;
  d_normalForce = d_prevNormalForce = 0;
  d_tgtForce = d_prevTgtForce = 0;
  d_tgtDisp = d_prevTgtDisp = 0;
  d_tgtDispStart = 0;
  d_tgtSlide = d_prevTgtSlide = false;
  d_tgtPeak = 0;

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
DEMContact::isOverlapped()
{
  // v[0] is the point on p2, v[1] is the point on p1
  REAL coef1[10], coef2[10];
  d_p1->getGlobalCoef(coef1); 
  d_p2->getGlobalCoef(coef2);

  d_radius1 = d_p1->getRadius(d_point1);
  d_radius2 = d_p2->getRadius(d_point2);


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
  d_penetr = vnormL2(d_point1 - d_point2);
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
    //          << " penetr = " << d_penetr << "\n"
    //          << " \t Pos1 = " << d_p1->currentPosition()
    //          << " \t Pos2 = " << d_p2->currentPosition() << "\n"
    //          << " \t b1 = " << b1 << " b2 = " << b2 << "\n"
    //          << " \t coef1 = [" << printCoef(coef1) << "] \n"
    //          << " \t coef2 = [" << printCoef(coef2) << "] \n"
    //          << " \t point1 = " << d_point1 << " point2 = " << d_point2 << "\n"
              //<< " \t point1_old = " << v2_old
              //<< " point2_old = " << v1_old
    //          << " radius1 = " << d_radius1 << " radius2 = " << d_radius2
    //          << "\n";
  }
  */

  if (b1 && b2 &&
      d_penetr / (2.0 * fmax(d_radius1, d_radius2)) >
        util::getParam<REAL>("minRelaOverlap") &&
      nearbyint(d_penetr / util::getParam<REAL>("measureOverlap")) >=
        1) { // a strict detection method
    d_isInContact = true;
    return true;
  } else {
    d_isInContact = false;
    return false;
  }
}

void
DEMContact::checkinPrevTgt(ContactTangentArray& contactTgtVec)
{
  for (auto& it : contactTgtVec) {
    if (it.d_ptcl1 == d_p1->getId() && it.d_ptcl2 == d_p2->getId()) {
      d_prevTgtForce = it.d_tgtForce;
      d_prevTgtDisp = it.d_tgtDisp;
      d_prevTgtLoading = it.d_tgtLoading;
      d_tgtDispStart = it.d_tgtDispStart;
      d_tgtPeak = it.d_tgtPeak;
      d_prevTgtSlide = it.d_tgtSlide;
      /*
      if ((it.d_ptcl1 == 2 && it.d_ptcl2 == 94) ||
          (it.d_ptcl1 == 94 && it.d_ptcl2 == 2)) {
        //std::cout << "DEMContact tangents: " << std::setprecision(16)
        //          << " TgtForce =  " << d_prevTgtForce
        //          << " TgtDisp =  " << d_prevTgtDisp << "\n";
      }
      */
      break;
    }
  }
}

void
DEMContact::checkoutTgt(ContactTangentArray& contactTgtVec)
{
  contactTgtVec.push_back(ContactTgt(d_p1->getId(), d_p2->getId(), d_tgtForce,
                                     d_tgtDisp, d_tgtLoading, d_tgtDispStart,
                                     d_tgtPeak, d_tgtSlide));
}

// isOverlapped() has been called in findContact() in dem.cpp and
// information recorded,
// now this function is called by internalForce() in dem.cpp.
void
DEMContact::contactForce()
{
  if (!d_isInContact) {
    d_isInContact = false;
    d_tgtLoading = false;
    d_tgtPeak = 0;
    d_normalForce = 0;
    d_tgtForce = 0;
    d_tgtDisp = 0; // total value
    d_normalDirc = 0;
    d_tgtDirc = 0;

    d_penetr = 0;
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
    //          << " f1 = " << d_p1->getForce() << "\n\t"
    //          << " f2 = " << d_p2->getForce() << "\n\t"
    //          << " penetr=" << d_penetr << "\n\t" << "\n";
  }
  */

  REAL young = util::getParam<REAL>("young");
  REAL poisson = util::getParam<REAL>("poisson");
  REAL maxRelaOverlap = util::getParam<REAL>("maxRelaOverlap");
  REAL measureOverlap = util::getParam<REAL>("measureOverlap");
  REAL contactCohesion = util::getParam<REAL>("contactCohesion");
  REAL contactDamp = util::getParam<REAL>("contactDamp");
  REAL contactFric = util::getParam<REAL>("contactFric");

  // obtain normal force, using absolute equation instead of stiffness method
  d_p1->setContactNum(d_p1->getContactNum() + 1);
  d_p2->setContactNum(d_p2->getContactNum() + 1);
  d_p1->setInContact(true);
  d_p2->setInContact(true);

  d_R0 = d_radius1 * d_radius2 / (d_radius1 + d_radius2);
  d_E0 = 0.5 * young / (1 - poisson * poisson);
  REAL allowedOverlap = 2.0 * fmin(d_radius1, d_radius2) * maxRelaOverlap;
  if (d_penetr > allowedOverlap) {
    std::stringstream inf;
    inf.setf(std::ios::scientific, std::ios::floatfield);
    inf << " DEMContact.cpp: iter=" << std::setw(8) << iteration
        << " ptcl1=" << std::setw(8) << getP1()->getId()
        << " ptcl2=" << std::setw(8) << getP2()->getId()
        << " penetr=" << std::setw(OWID) << d_penetr
        << " allow=" << std::setw(OWID) << allowedOverlap << std::endl;
    MPI_Status status;
    int length = OWID * 2 + 8 * 3 + 19 + 7 * 3 + 8 + 1;
    MPI_File_write_shared(overlapInf, const_cast<char*>(inf.str().c_str()),
                          length, MPI_CHAR, &status);

    d_penetr = allowedOverlap;
  }

  d_penetr = nearbyint(d_penetr / measureOverlap) * measureOverlap;
  d_contactRadius = sqrt(d_penetr * d_R0);
  // d_normalDirc points from particle 1 to particle 2
  d_normalDirc = normalize(d_point1 - d_point2);
  // normalForce pointing to particle 1 // pow(d_penetr, 1.5)
  d_normalForce = -sqrt(d_penetr * d_penetr * d_penetr) * sqrt(d_R0) * 4 *
                  d_E0 / 3 * d_normalDirc;

  // apply cohesion force
  d_cohesionForce = Pi * (d_penetr * d_R0) * contactCohesion * d_normalDirc;

  // obtain normal damping force
  Vec cp = (d_point1 + d_point2) / 2;
  Vec veloc1 =
    d_p1->currentVelocity() + cross(d_p1->currentOmega(), (cp - d_p1->currentPosition()));
  Vec veloc2 =
    d_p2->currentVelocity() + cross(d_p2->currentOmega(), (cp - d_p2->currentPosition()));
  REAL m1 = getP1()->getMass();
  REAL m2 = getP2()->getMass();
  REAL kn = pow(6 * vnormL2(d_normalForce) * d_R0 * pow(d_E0, 2), 1.0 / 3.0);
  REAL dampCritical = 2 * sqrt(m1 * m2 / (m1 + m2) * kn); // critical damping
  Vec cntDampingForce = contactDamp * dampCritical *
                        dot(veloc1 - veloc2, d_normalDirc) * d_normalDirc;

  d_vibraTimeStep = 2.0 * sqrt(m1 * m2 / (m1 + m2) / kn);
  Vec relativeVel = veloc1 - veloc2;
  d_impactTimeStep =
    (relativeVel.lengthSq() < std::numeric_limits<double>::min())
      ? std::numeric_limits<double>::max()
      : allowedOverlap / fabs(dot(relativeVel, d_normalDirc));

  // obtain tangential force
  if (contactFric != 0) {
    d_G0 = young / 2 / (1 + poisson);
    // RelaDispInc points along point1's displacement relative to point2
    Vec RelaDispInc = (veloc1 - veloc2) * timeStep;
    Vec tgtDispInc = RelaDispInc - dot(RelaDispInc, d_normalDirc) * d_normalDirc;
    // prevTgtDisp read by checkinPrevTgt()
    d_tgtDisp = d_prevTgtDisp + tgtDispInc;
    if (vnormL2(d_tgtDisp) == 0) {
      d_tgtDirc = 0;
    } else {
      // tgtDirc points along Tgtential forces exerted on particle 1
      d_tgtDirc = normalize(-d_tgtDisp);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // linear friction model
    REAL fP = contactFric * vnormL2(d_normalForce);
    REAL ks = 4 * d_G0 * d_contactRadius / (2 - poisson);
    // d_prevTgtForce read by CheckinPreTgt()
    d_tgtForce = d_prevTgtForce + ks * (-tgtDispInc);
    if (vnormL2(d_tgtForce) > fP) {
      d_tgtForce = fP * d_tgtDirc;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MINDLIN_ASSUMED
    computeTangentForceMindlinAssumed(contactFric, poisson, tgtDispInc);
#endif
#ifdef MINDLIN_KNOWN
    computeTangentForceMindlinKnown(contactFric, poisson, tgtDispInc);
#endif
  }

  // apply forces
  Vec totalForce = d_normalForce +  d_tgtForce + d_cohesionForce
                   - cntDampingForce;
  Vec momentArm = cp - d_p1->currentPosition();
  Vec totalMoment = cross(momentArm, (d_normalForce + d_tgtForce - cntDampingForce));

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
    //          << " f1 = " << d_p1->getForce() << "\n\t"
    //          << " f2 = " << d_p2->getForce() << "\n\t"
    //          << " penetr=" << d_penetr << "\n\t" << "\n";
  }
  */
  /*
    //std::cout << " cohesionForce=" << d_cohesionForce << "\n\t"
              << " normalForce=" << d_normalForce << "\n\t"
              << " tgtForce =" << d_tgtForce << "\n\t"
              << " dampingForce = " << cntDampingForce << "\n\t"
              << " totalForce = " << totalForce << "\n\t"
              << " totalMoment = " << totalMoment << "\n\t"
              << " accumulated time=" << iteration * timeStep << "\n"
              << "\t cp = " << cp << "\n\t"
              << " d_point1 = " << d_point1 << " d_point2 = " << d_point2
              << "\n\t"
              << " V1 = " << d_p1->currentVelocity()
              << " V2 = " << d_p2->currentVelocity() << "\n\t"
              << " W1 = " << d_p1->currentOmega()
              << " W2 = " << d_p2->currentOmega() << "\n\t"
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
                                           const Vec& tgtDispInc)
{
  REAL val = 0, ks = 0;
  REAL fP = contactFric * vnormL2(d_normalForce);
  d_tgtLoading = (dot(d_prevTgtDisp, tgtDispInc) >= 0);

  if (d_tgtLoading) {        // loading
    if (!d_prevTgtLoading) { // pre-step is unloading
      val = 8 * d_G0 * d_contactRadius * vnormL2(tgtDispInc) /
            (3 * (2 - poisson) * fP);
      d_tgtDispStart = d_prevTgtDisp;
    } else // pre-step is loading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tgtDisp - d_tgtDispStart) /
            (3 * (2 - poisson) * fP);

    if (val > 1.0)
      d_tgtForce = fP * d_tgtDirc;
    else {
      ks = 4 * d_G0 * d_contactRadius / (2 - poisson) * sqrt(1 - val);
      // incremental method
      d_tgtForce =
        d_prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
      // total value method: tgtForce = fP*(1-pow(1-val, 1.5))*d_tgtDirc;
    }
  } else {                  // unloading
    if (d_prevTgtLoading) { // pre-step is loading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tgtDisp - d_tgtDispStart) /
            (3 * (2 - poisson) * fP);
      d_tgtPeak = vnormL2(d_prevTgtForce);
    } else // pre-step is unloading
      val = 8 * d_G0 * d_contactRadius * vnormL2(d_tgtDisp - d_tgtDispStart) /
            (3 * (2 - poisson) * fP);

    if (val > 1.0 || d_tgtPeak > fP)
      d_tgtForce = fP * d_tgtDirc;
    else {
      ks = 2 * sqrt(2) * d_G0 * d_contactRadius / (2 - poisson) *
           sqrt(1 + pow(1 - d_tgtPeak / fP, 2.0 / 3.0) + val);
      // incremental method
      d_tgtForce =
        d_prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
      // total value method: d_tgtForce = (d_tgtPeak-2*fP*(1-sqrt(2)/4*pow(1+
      // pow(1-d_tgtPeak/fP,2.0/3.0) + val,1.5)))*d_tgtDirc;
    }
  }

  if (vnormL2(d_tgtForce) > fP)
    d_tgtForce = fP * d_tgtDirc;
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
                                         const Vec& tgtDispInc)
{
  REAL val = 0, ks = 0;
  REAL fP = contactFric * vnormL2(d_normalForce);
  if (d_prevTgtSlide)
    val =
      8 * d_G0 * d_contactRadius * vnormL2(tgtDispInc) / (3 * (2 - poisson) * fP);
  else
    val = 8 * d_G0 * d_contactRadius * vnormL2(d_tgtDisp - d_tgtDispStart) /
          (3 * (2 - poisson) * fP);

  if (iteration > 10000 &&
      iteration < 11000) { // loading (and possible sliding)
    if (val > 1.0) {
      d_tgtForce = fP * d_tgtDirc;
      d_tgtSlide = true;
    } else {
      if (!d_prevTgtSlide) {
        ks = 4 * d_G0 * d_contactRadius / (2 - poisson) * sqrt(1 - val);
        d_tgtForce =
          d_prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
        d_tgtSlide = false;
      } else {
        if (vnormL2(d_tgtForce) > vnormL2(d_prevTgtForce))
          d_tgtSlide = true;
        else
          d_tgtSlide = false;
      }
    }
    d_tgtPeak = vnormL2(d_tgtForce);
  } else { // (possible sliding and) unloading
    if (val > 1.0 || d_tgtPeak > fP) {
      d_tgtForce = fP * d_tgtDirc;
      d_tgtSlide = true;
    } else {
      if (!d_prevTgtSlide) {
        ks = 2 * sqrt(2) * d_G0 * d_contactRadius / (2 - poisson) *
             sqrt(1 + pow(1 - d_tgtPeak / fP, 2.0 / 3.0) + val);
        d_tgtForce =
          d_prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
        d_tgtSlide = false;
      } else {
        if (vnormL2(d_tgtForce) > vnormL2(d_prevTgtForce))
          d_tgtSlide = true;
        else {
          d_tgtSlide = false;
          d_tgtDispStart = d_tgtDisp;
        }
      }
    }
  }

  /*
        //std::cout<< "DEMContact.h: iter="<iteration
                 << " prevTgtSlide=" << d_prevTgtSlide
                 << " tgtSlide=" << d_tgtSlide
                 << " val=" << val
                 << " ks=" << ks
                 << " tgtDispInc.x=" << tgtDispInc.x()
                 << " d_prevTgtForce=" << vnormL2(d_prevTgtForce)
                 << " d_tgtForce" << vnormL2(d_tgtForce)
                 << std::endl;
        */

  if (vnormL2(d_tgtForce) > fP)
    d_tgtForce = fP * d_tgtDirc;
}

} // namespace dem ends
