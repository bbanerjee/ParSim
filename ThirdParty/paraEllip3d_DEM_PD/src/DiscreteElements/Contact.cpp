#include <DiscreteElements/Contact.h>
#include <InputOutput/Parameter.h>
#include <Core/Const/const.h>
#include <Core/Math/root6.h>
#include <iostream>
#include <fstream>
#include <iomanip>

//#define MINDLIN_ASSUMED
//#define MINDLIN_KNOWN

namespace dem {

Contact::Contact() {
  p1 = NULL;
  p2 = NULL;
  penetr = 0;
  contactRadius = 0;
  point1 = point2 = 0;
  radius1 = radius2 = 0;
  normalDirc = tgtDirc = 0;

  isInContact = false;
  tgtLoading = prevTgtLoading = true;
  normalForce = prevNormalForce = 0;
  tgtForce = prevTgtForce = 0;
  tgtDisp = prevTgtDisp = 0;
  tgtDispStart = 0;
  tgtSlide = prevTgtSlide = false;
  tgtPeak = 0;

  cohesionForce = 0;
  spinResist = 0;

  E0 = G0 = R0 = 0;
}

Contact::Contact(Particle* t1, Particle* t2) {
  p1 = t1;
  p2 = t2;
  penetr = 0;
  contactRadius = 0;
  point1 = point2 = 0;
  radius1 = radius2 = 0;
  normalDirc = tgtDirc = 0;

  isInContact = false;
  tgtLoading = prevTgtLoading = true;
  normalForce = prevNormalForce = 0;
  tgtForce = prevTgtForce = 0;
  tgtDisp = prevTgtDisp = 0;
  tgtDispStart = 0;
  tgtSlide = prevTgtSlide = false;
  tgtPeak = 0;

  cohesionForce = 0;
  spinResist = 0;

  E0 = G0 = R0 = 0;
}

bool Contact::isRedundant(const Contact &other) const {
  std::size_t id1 = getP1()->getId();
  std::size_t id2 = getP2()->getId();
  std::size_t oId1 = (other.getP1())->getId();
  std::size_t oId2 = (other.getP2())->getId();

  return ((id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2));
}

bool Contact::operator==(const Contact &other) const {
  std::size_t id1 = getP1()->getId();
  std::size_t id2 = getP2()->getId();
  std::size_t oId1 = (other.getP1())->getId();
  std::size_t oId2 = (other.getP2())->getId();

  return ((id2 == oId1 && id1 == oId2) || (id1 == oId1 && id2 == oId2));
}

Particle* Contact::getP1() const { return p1; }

Particle* Contact::getP2() const { return p2; }

bool Contact::isOverlapped() {
  REAL coef1[10], coef2[10];
  p1->getGlobalCoef(coef1); // v[0] is the point on p2, v[1] is the point on p1
  p2->getGlobalCoef(coef2);
  Vec v[2];
  bool b1 = root6(coef1, coef2, v[0]);
  bool b2 = root6(coef2, coef1, v[1]);
  point1 = v[1];
  point2 = v[0];
  radius1 = p1->getRadius(point1);
  radius2 = p2->getRadius(point2);
  penetr = vfabs(point1 - point2);

  if (b1 && b2 &&
      penetr / (2.0 * fmax(radius1, radius2)) >
          dem::Parameter::getSingleton().parameter["minRelaOverlap"] &&
      nearbyint(penetr /
                dem::Parameter::getSingleton().parameter["measureOverlap"]) >=
          1) { // a strict detection method
    isInContact = true;
    return true;
  } else {
    isInContact = false;
    return false;
  }
}

void Contact::checkinPrevTgt(ContactTangentArray &contactTgtVec) {
  for (ContactTangentArray::iterator it = contactTgtVec.begin();
       it != contactTgtVec.end(); ++it) {
    if (it->ptcl1 == p1->getId() && it->ptcl2 == p2->getId()) {
      prevTgtForce = it->tgtForce;
      prevTgtDisp = it->tgtDisp;
      prevTgtLoading = it->tgtLoading;
      tgtDispStart = it->tgtDispStart;
      tgtPeak = it->tgtPeak;
      prevTgtSlide = it->tgtSlide;
      break;
    }
  }
}

void Contact::checkoutTgt(ContactTangentArray &contactTgtVec) {
  contactTgtVec.push_back(
      ContactTgt(p1->getId(), p2->getId(), tgtForce, tgtDisp, tgtLoading,
                 tgtDispStart, tgtPeak, tgtSlide));
}

void Contact::contactForce() {
  // isOverlapped() has been called in findContact() in assembly.cpp and
  // information recorded,
  // now this function is called by internalForce() in assembly.cpp.

  if (isInContact) {

    REAL young = dem::Parameter::getSingleton().parameter["young"];
    REAL poisson = dem::Parameter::getSingleton().parameter["poisson"];
    REAL maxRelaOverlap =
        dem::Parameter::getSingleton().parameter["maxRelaOverlap"];
    REAL measureOverlap =
        dem::Parameter::getSingleton().parameter["measureOverlap"];
    REAL contactCohesion =
        dem::Parameter::getSingleton().parameter["contactCohesion"];
    REAL contactDamp = dem::Parameter::getSingleton().parameter["contactDamp"];
    REAL contactFric = dem::Parameter::getSingleton().parameter["contactFric"];

    // obtain normal force, using absolute equation instead of stiffness method
    p1->setContactNum(p1->getContactNum() + 1);
    p2->setContactNum(p2->getContactNum() + 1);
    p1->setInContact(true);
    p2->setInContact(true);

    R0 = radius1 * radius2 / (radius1 + radius2);
    E0 = 0.5 * young / (1 - poisson * poisson);
    REAL allowedOverlap = 2.0 * fmin(radius1, radius2) * maxRelaOverlap;
    if (penetr > allowedOverlap) {
      std::stringstream inf;
      inf.setf(std::ios::scientific, std::ios::floatfield);
      inf << " Contact.cpp: iter=" << std::setw(8) << iteration
          << " ptcl1=" << std::setw(8) << getP1()->getId()
          << " ptcl2=" << std::setw(8) << getP2()->getId()
          << " penetr=" << std::setw(OWID) << penetr
          << " allow=" << std::setw(OWID) << allowedOverlap << std::endl;
      MPI_Status status;
      int length = OWID * 2 + 8 * 3 + 19 + 7 * 3 + 8 + 1;
      MPI_File_write_shared(overlapInf, const_cast<char *>(inf.str().c_str()),
                            length, MPI_CHAR, &status);

      penetr = allowedOverlap;
    }

    penetr = nearbyint(penetr / measureOverlap) * measureOverlap;
    contactRadius = sqrt(penetr * R0);
    normalDirc = normalize(
        point1 - point2);     // normalDirc points from particle 1 to particle 2
    normalForce = -sqrt(penetr * penetr * penetr) * sqrt(R0) * 4 * E0 / 3 *
                  normalDirc; // normalForce pointing to particle 1
                              // pow(penetr, 1.5)

    // apply cohesion force
    cohesionForce = Pi * (penetr * R0) * contactCohesion * normalDirc;
    p1->addForce(cohesionForce);
    p2->addForce(-cohesionForce);

    // apply normal force
    p1->addForce(normalForce);
    p2->addForce(-normalForce);
    p1->addMoment(((point1 + point2) / 2 - p1->getCurrPos()) % normalForce);
    p2->addMoment(((point1 + point2) / 2 - p2->getCurrPos()) % (-normalForce));

    /*
      std::cout<< "Contact.h: iter=" << iteration
	       << " penetr=" << penetr
	       << " cohesionForce=" << vfabs(cohesionForce)
	       << " normalForce=" << vfabs(normalForce)
	       << " accumulated time=" << iteration * timeStep
	       << std::endl;
      */

    // obtain normal damping force
    Vec cp = (point1 + point2) / 2;
    Vec veloc1 =
        p1->getCurrVeloc() + p1->getCurrOmga() % (cp - p1->getCurrPos());
    Vec veloc2 =
        p2->getCurrVeloc() + p2->getCurrOmga() % (cp - p2->getCurrPos());
    REAL m1 = getP1()->getMass();
    REAL m2 = getP2()->getMass();
    REAL kn = pow(6 * vfabs(normalForce) * R0 * pow(E0, 2), 1.0 / 3.0);
    REAL dampCritical = 2 * sqrt(m1 * m2 / (m1 + m2) * kn); // critical damping
    Vec cntDampingForce = contactDamp * dampCritical *
                          ((veloc1 - veloc2) * normalDirc) * normalDirc;
    vibraTimeStep = 2.0 * sqrt(m1 * m2 / (m1 + m2) / kn);
    impactTimeStep = allowedOverlap / fabs((veloc1 - veloc2) * normalDirc);

    // apply normal damping force
    p1->addForce(-cntDampingForce);
    p2->addForce(cntDampingForce);
    p1->addMoment(((point1 + point2) / 2 - p1->getCurrPos()) %
                  (-cntDampingForce));
    p2->addMoment(((point1 + point2) / 2 - p2->getCurrPos()) % cntDampingForce);

    if (contactFric != 0) {
      // obtain tangential force
      G0 = young / 2 / (1 + poisson); // RelaDispInc points along point1's
                                      // displacement relative to point2
      Vec RelaDispInc = (veloc1 - veloc2) * timeStep;
      Vec tgtDispInc = RelaDispInc - (RelaDispInc * normalDirc) * normalDirc;
      tgtDisp =
          prevTgtDisp + tgtDispInc; // prevTgtDisp read by checkinPrevTgt()
      if (vfabs(tgtDisp) == 0)
        tgtDirc = 0;
      else
        tgtDirc = normalize(-tgtDisp); // tgtDirc points along Tgtential forces
                                       // exerted on particle 1

      REAL fP = 0;
      REAL ks = 0;

      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // linear friction model
      fP = contactFric * vfabs(normalForce);
      ks = 4 * G0 * contactRadius / (2 - poisson);
      tgtForce = prevTgtForce +
                 ks * (-tgtDispInc); // prevTgtForce read by CheckinPreTgt()
      if (vfabs(tgtForce) > fP)
        tgtForce = fP * tgtDirc;
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mindlin's model (loading/unloading condition assumed)
// This model is not recommended as it is impossible to strictly determine
// loading/unloading condition
// unless load is known (the case of pure moment rotation).
#ifdef MINDLIN_ASSUMED
      REAL val = 0;
      fP = contactFric * vfabs(normalForce);
      tgtLoading = (prevTgtDisp * tgtDispInc >= 0);

      if (tgtLoading) {        // loading
        if (!prevTgtLoading) { // pre-step is unloading
          val = 8 * G0 * contactRadius * vfabs(tgtDispInc) /
                (3 * (2 - poisson) * fP);
          tgtDispStart = prevTgtDisp;
        } else // pre-step is loading
          val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart) /
                (3 * (2 - poisson) * fP);

        if (val > 1.0)
          tgtForce = fP * tgtDirc;
        else {
          ks = 4 * G0 * contactRadius / (2 - poisson) * sqrt(1 - val);
          //incremental method
          tgtForce =
              prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
          //total value method: tgtForce = fP*(1-pow(1-val, 1.5))*tgtDirc;
        }
      } else {                // unloading
        if (prevTgtLoading) { // pre-step is loading
          val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart) /
                (3 * (2 - poisson) * fP);
          tgtPeak = vfabs(prevTgtForce);
        } else // pre-step is unloading
          val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart) /
                (3 * (2 - poisson) * fP);

        if (val > 1.0 || tgtPeak > fP)
          tgtForce = fP * tgtDirc;
        else {
          ks = 2 * sqrt(2) * G0 * contactRadius / (2 - poisson) *
               sqrt(1 + pow(1 - tgtPeak / fP, 2.0 / 3.0) + val);
          //incremental method
          tgtForce =
              prevTgtForce + ks * (-tgtDispInc); // tgtDispInc determines signs
          //total value method: tgtForce = (tgtPeak-2*fP*(1-sqrt(2)/4*pow(1+
          //pow(1-tgtPeak/fP,2.0/3.0) + val,1.5)))*tgtDirc;
        }
      }

      if (vfabs(tgtForce) > fP)
        tgtForce = fP * tgtDirc;
#endif
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mindlin's model (loading/unloading condition known for pure moment rotation
// case)
// As loading/unloading condition is known, both incremental and total value
// method work well.
// Herein sliding history is incorporated.
#ifdef MINDLIN_KNOWN
      REAL val = 0;
      fP = contactFric * vfabs(normalForce);
      if (prevTgtSlide)
        val = 8 * G0 * contactRadius * vfabs(tgtDispInc) /
              (3 * (2 - poisson) * fP);
      else
        val = 8 * G0 * contactRadius * vfabs(tgtDisp - tgtDispStart) /
              (3 * (2 - poisson) * fP);

      if (iteration > 10000 &&
          iteration < 11000) { // loading (and possible sliding)
        if (val > 1.0) {
          tgtForce = fP * tgtDirc;
          tgtSlide = true;
        } else {
          if (!prevTgtSlide) {
            ks = 4 * G0 * contactRadius / (2 - poisson) * sqrt(1 - val);
            tgtForce = prevTgtForce +
                       ks * (-tgtDispInc); // tgtDispInc determines signs
            tgtSlide = false;
          } else {
            if (vfabs(tgtForce) > vfabs(prevTgtForce))
              tgtSlide = true;
            else
              tgtSlide = false;
          }
        }
        tgtPeak = vfabs(tgtForce);
      } else { // (possible sliding and) unloading
        if (val > 1.0 || tgtPeak > fP) {
          tgtForce = fP * tgtDirc;
          tgtSlide = true;
        } else {
          if (!prevTgtSlide) {
            ks = 2 * sqrt(2) * G0 * contactRadius / (2 - poisson) *
                 sqrt(1 + pow(1 - tgtPeak / fP, 2.0 / 3.0) + val);
            tgtForce = prevTgtForce +
                       ks * (-tgtDispInc); // tgtDispInc determines signs
            tgtSlide = false;
          } else {
            if (vfabs(tgtForce) > vfabs(prevTgtForce))
              tgtSlide = true;
            else {
              tgtSlide = false;
              tgtDispStart = tgtDisp;
            }
          }
        }
      }

      /*
     	std::cout<< "Contact.h: iter="<iteration
     		 << " prevTgtSlide=" << prevTgtSlide
     		 << " tgtSlide=" << tgtSlide
     		 << " val=" << val
     		 << " ks=" << ks
     		 << " tgtDispInc.x=" << tgtDispInc.getX()
     		 << " prevTgtForce=" << vfabs(prevTgtForce)
     		 << " tgtForce" << vfabs(tgtForce)
     		 << std::endl;
     	*/

      if (vfabs(tgtForce) > fP)
        tgtForce = fP * tgtDirc;
#endif
      /////////////////////////////////////////////////////////////////////////////////////////////////////////

      // apply tangential force
      p1->addForce(tgtForce);
      p2->addForce(-tgtForce);
      p1->addMoment(((point1 + point2) / 2 - p1->getCurrPos()) % tgtForce);
      p2->addMoment(-((point1 + point2) / 2 - p2->getCurrPos()) % tgtForce);
    }

  } else {
    isInContact = false;
    tgtLoading = false;
    tgtPeak = 0;
    normalForce = 0;
    tgtForce = 0;
    tgtDisp = 0; //total value
    normalDirc = 0;
    tgtDirc = 0;

    penetr = 0;
    contactRadius = 0;
    radius1 = radius2 = 0;
    spinResist = 0;
  }

}

} // namespace dem ends
