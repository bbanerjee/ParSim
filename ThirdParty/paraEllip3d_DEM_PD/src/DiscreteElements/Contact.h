#ifndef CONTACT_H
#define CONTACT_H

#include <Core/Types/realtypes.h>
#include <DiscreteElements/Particle.h>
#include <cstddef>
#include <vector>
#include <boost/functional/hash.hpp>

namespace dem {

class ContactTgt {

public:
  std::size_t ptcl1;
  std::size_t ptcl2;
  Vec tgtForce;
  Vec tgtDisp;
  bool tgtLoading;
  Vec tgtDispStart;
  REAL tgtPeak;
  bool tgtSlide;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &ptcl1;
    ar &ptcl2;
    ar &tgtForce;
    ar &tgtDisp;
    ar &tgtLoading;
    ar &tgtDispStart;
    ar &tgtPeak;
    ar &tgtSlide;
  }

public:
  ContactTgt()
      : ptcl1(0), ptcl2(0), tgtForce(0), tgtDisp(0), tgtLoading(0),
        tgtDispStart(0), tgtPeak(0), tgtSlide(false) {}

  ContactTgt(std::size_t _ptcl1, std::size_t _ptcl2, Vec _tf, Vec _td, bool _tl,
             Vec _tds, REAL _tp, bool _ts)
      : ptcl1(_ptcl1), ptcl2(_ptcl2), tgtForce(_tf), tgtDisp(_td),
        tgtLoading(_tl), tgtDispStart(_tds), tgtPeak(_tp), tgtSlide(_ts) {}

};

class Contact {
public:
  Contact();
  Contact(Particle *t1, Particle *t2);

  Particle *getP1() const;
  Particle *getP2() const;
  Vec getPoint1() const { return point1; }
  Vec getPoint2() const { return point2; }
  REAL getRadius1() const { return radius1; }
  REAL getRadius2() const { return radius2; }
  REAL getR0() const { return R0; }
  REAL getE0() const { return E0; }
  REAL getVibraTimeStep() const { return vibraTimeStep; }
  REAL getImpactTimeStep() const { return impactTimeStep; }

  bool isOverlapped();
  void contactForce(); // calculate normal and tangential force of contact
  REAL getNormalForce() const { return vfabs(normalForce); }
  REAL getTgtForce() const { return vfabs(tgtForce); }
  REAL getPenetration() const { return penetr; }
  REAL getContactRadius() const { return contactRadius; }
  REAL getTgtDisp() const {
    return vfabs(tgtDisp);
  } // total value during a process of contact
  void checkoutTgt(std::vector<ContactTgt> &contactTgtVec);
  void checkinPrevTgt(std::vector<ContactTgt> &contactTgtVec);
  Vec normalForceVec() const { return normalForce; }
  Vec tgtForceVec() const { return tgtForce; }
  bool isRedundant(const Contact &other) const;
  bool operator==(const Contact &other) const;

private:
  Particle *p1;       // particle 1
  Particle *p2;       // particle 2
  REAL penetr;        // penetr
  REAL contactRadius; // radius of contact surface
  Vec point1;         // point1 on particle 1, innermost to particle 2
  Vec point2;         // point2 on particle 2, innermost to particle 1
  REAL radius1;       // radius of osculating circles at point1
  REAL radius2;       // radius of osculating circles at point2
  Vec normalDirc; // normal direction, pointing from particle 1 to particle 2
  Vec tgtDirc;    // tangential direction

  bool isInContact; // are p1 and p1 in contact
  bool tgtLoading;  // tangential loading or unloading
  Vec normalForce;  // pointing from particle 2 to paticle 1
  Vec tgtForce; // TgtrDirc points along tangential forces exerted on particle 1
  Vec tgtDisp;  // tangential relative displacment total vector
  Vec tgtDispStart; // displacement start value for each loading-unloading loop
  bool tgtSlide;    // tangential silde or not

  bool prevTgtLoading; // previous loading-unloading status
  Vec prevNormalForce;
  Vec prevTgtForce;
  Vec prevTgtDisp; // previous tangential relative displacment total vector
  bool prevTgtSlide;
  REAL tgtPeak;

  Vec cohesionForce; // cohesion force between particles
  Vec spinResist;

  REAL E0;
  REAL G0;
  REAL R0;
  REAL vibraTimeStep;
  REAL impactTimeStep;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &p1;
    ar &p2;
    ar &penetr;
    ar &contactRadius;
    ar &point1;
    ar &point2;
    ar &radius1;
    ar &radius2;
    ar &normalDirc;
    ar &tgtDirc;
    ar &isInContact;
    ar &tgtLoading;
    ar &normalForce;
    ar &tgtForce;
    ar &tgtDisp;
    ar &tgtDispStart;
    ar &tgtSlide;
    ar &prevTgtLoading;
    ar &prevNormalForce;
    ar &prevTgtForce;
    ar &prevTgtDisp;
    ar &prevTgtSlide;
    ar &tgtPeak;
    ar &cohesionForce;
    ar &spinResist;
    ar &E0;
    ar &G0;
    ar &R0;
    ar &vibraTimeStep;
    ar &impactTimeStep;
  }

public:
  friend std::size_t hash_value(const Contact &c) {
    boost::hash<std::size_t> hasher;
    return hasher(c.getP1()->getId() * c.getP2()->getId());
  }

};

}

#endif
