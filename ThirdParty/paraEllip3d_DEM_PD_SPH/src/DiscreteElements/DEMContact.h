#ifndef CONTACT_H
#define CONTACT_H

#include <Core/Types/RealTypes.h>
#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/DEMParticle.h>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <vector>

namespace dem {

class ContactTgt
{

public:
  std::size_t d_ptcl1;
  std::size_t d_ptcl2;
  Vec d_tgtForce;
  Vec d_tgtDisp;
  bool d_tgtLoading;
  Vec d_tgtDispStart;
  REAL d_tgtPeak;
  bool d_tgtSlide;

  ContactTgt()
    : d_ptcl1(0)
    , d_ptcl2(0)
    , d_tgtForce(0)
    , d_tgtDisp(0)
    , d_tgtLoading(0)
    , d_tgtDispStart(0)
    , d_tgtPeak(0)
    , d_tgtSlide(false)
  {
  }

  ContactTgt(std::size_t _ptcl1, std::size_t _ptcl2, Vec _tf, Vec _td, bool _tl,
             Vec _tds, REAL _tp, bool _ts)
    : d_ptcl1(_ptcl1)
    , d_ptcl2(_ptcl2)
    , d_tgtForce(_tf)
    , d_tgtDisp(_td)
    , d_tgtLoading(_tl)
    , d_tgtDispStart(_tds)
    , d_tgtPeak(_tp)
    , d_tgtSlide(_ts)
  {
  }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_ptcl1;
    ar& d_ptcl2;
    ar& d_tgtForce;
    ar& d_tgtDisp;
    ar& d_tgtLoading;
    ar& d_tgtDispStart;
    ar& d_tgtPeak;
    ar& d_tgtSlide;
  }
};

class DEMContact
{
public:
  DEMContact();
  DEMContact(DEMParticle* t1, DEMParticle* t2);

  DEMParticle* getP1() const;
  DEMParticle* getP2() const;
  Vec getPoint1() const { return d_point1; }
  Vec getPoint2() const { return d_point2; }
  REAL getRadius1() const { return d_radius1; }
  REAL getRadius2() const { return d_radius2; }
  REAL getR0() const { return d_R0; }
  REAL getE0() const { return d_E0; }
  REAL getVibraTimeStep() const { return d_vibraTimeStep; }
  REAL getImpactTimeStep() const { return d_impactTimeStep; }

  bool isOverlapped();
  void contactForce(); // calculate normal and tangential force of contact
  REAL getNormalForce() const { return vnormL2(d_normalForce); }
  REAL getTgtForce() const { return vnormL2(d_tgtForce); }
  REAL getPenetration() const { return d_penetr; }
  REAL getContactRadius() const { return d_contactRadius; }
  REAL getTgtDisp() const
  {
    return vnormL2(d_tgtDisp);
  } // total value during a process of contact
  void checkoutTgt(std::vector<ContactTgt>& contactTgtVec);
  void checkinPrevTgt(std::vector<ContactTgt>& contactTgtVec);
  Vec normalForceVec() const { return d_normalForce; }
  Vec tgtForceVec() const { return d_tgtForce; }
  bool isRedundant(const DEMContact& other) const;
  bool operator==(const DEMContact& other) const;

private:
  DEMParticle* d_p1;       // particle 1
  DEMParticle* d_p2;       // particle 2
  REAL d_penetr;        // d_penetration
  REAL d_contactRadius; // radius of contact surface
  Vec d_point1;         // point1 on particle 1, innermost to particle 2
  Vec d_point2;         // point2 on particle 2, innermost to particle 1
  REAL d_radius1;       // radius of osculating circles at point1
  REAL d_radius2;       // radius of osculating circles at point2
  Vec d_normalDirc; // normal direction, pointing from particle 1 to particle 2
  Vec d_tgtDirc;    // tangential direction

  bool d_isInContact; // are p1 and p1 in contact
  bool d_tgtLoading;  // tangential loading or unloading
  Vec d_normalForce;  // pointing from particle 2 to paticle 1
  Vec
    d_tgtForce; // TgtrDirc points along tangential forces exerted on particle 1
  Vec d_tgtDisp; // tangential relative displacment total vector
  Vec
    d_tgtDispStart; // displacement start value for each loading-unloading loop
  bool d_tgtSlide;  // tangential silde or not

  bool d_prevTgtLoading; // previous loading-unloading status
  Vec d_prevNormalForce;
  Vec d_prevTgtForce;
  Vec d_prevTgtDisp; // previous tangential relative displacment total vector
  bool d_prevTgtSlide;
  REAL d_tgtPeak;

  Vec d_cohesionForce; // cohesion force between particles
  Vec d_spinResist;

  REAL d_E0;
  REAL d_G0;
  REAL d_R0;
  REAL d_vibraTimeStep;
  REAL d_impactTimeStep;

  void computeTangentForceMindlinAssumed(const REAL& contactFric,
                                         const REAL& poisson,
                                         const Vec& tgtDispInc);
  void computeTangentForceMindlinKnown(const REAL& contactFric,
                                       const REAL& poisson,
                                       const Vec& tgtDispInc);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_p1;
    ar& d_p2;
    ar& d_penetr;
    ar& d_contactRadius;
    ar& d_point1;
    ar& d_point2;
    ar& d_radius1;
    ar& d_radius2;
    ar& d_normalDirc;
    ar& d_tgtDirc;
    ar& d_isInContact;
    ar& d_tgtLoading;
    ar& d_normalForce;
    ar& d_tgtForce;
    ar& d_tgtDisp;
    ar& d_tgtDispStart;
    ar& d_tgtSlide;
    ar& d_prevTgtLoading;
    ar& d_prevNormalForce;
    ar& d_prevTgtForce;
    ar& d_prevTgtDisp;
    ar& d_prevTgtSlide;
    ar& d_tgtPeak;
    ar& d_cohesionForce;
    ar& d_spinResist;
    ar& d_E0;
    ar& d_G0;
    ar& d_R0;
    ar& d_vibraTimeStep;
    ar& d_impactTimeStep;
  }

public:
  friend std::size_t hash_value(const DEMContact& c)
  {
    boost::hash<std::size_t> hasher;
    return hasher(c.getP1()->getId() * c.getP2()->getId());
  }
};
}

#endif
