#ifndef ELLIP3D_DEM_CONTACT_H
#define ELLIP3D_DEM_CONTACT_H

#include <Core/Types/RealTypes.h>
#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/DEMContactTangent.h>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <vector>

namespace dem {

class DEMContact
{
public:
  DEMContact();
  DEMContact(DEMParticle* t1, DEMParticle* t2);

  DEMParticle* getP1() const;
  DEMParticle* getP2() const;
  Vec getPoint1() const { return d_point1; }
  Vec getPoint2() const { return d_point2; }
  REAL radius1() const { return d_radius1; }
  REAL radius2() const { return d_radius2; }
  REAL getR0() const { return d_R0; }
  REAL getE0() const { return d_E0; }
  REAL getVibraTimeStep() const { return d_vibraTimeStep; }
  REAL getImpactTimeStep() const { return d_impactTimeStep; }

  bool isOverlapped();
  void computeContactForces(); // calculate normal and tangential force of contact
  REAL getNormalForceMagnitude() const { return vnormL2(d_normalForce); }
  REAL getTangentForceMagnitude() const { return vnormL2(d_tangentForce); }
  REAL getPenetration() const { return d_penetration; }
  REAL getContactRadius() const { return d_contactRadius; }
  REAL getTangentDisplacement() const
  {
    return vnormL2(d_tangentDisplacement);
  } // total value during a process of contact
  void checkoutContactTangents(std::vector<DEMContactTangent>& contactTangentVec);
  void checkinPreviousContactTangents(std::vector<DEMContactTangent>& contactTangentVec);
  Vec getNormalForce() const { return d_normalForce; }
  Vec getTangentForce() const { return d_tangentForce; }
  bool isRedundant(const DEMContact& other) const;
  bool operator==(const DEMContact& other) const;

private:
  DEMParticle* d_p1;       // particle 1
  DEMParticle* d_p2;       // particle 2
  REAL d_penetration;        // d_penetration
  REAL d_contactRadius; // radius of contact surface
  Vec d_point1;         // point1 on particle 1, innermost to particle 2
  Vec d_point2;         // point2 on particle 2, innermost to particle 1
  REAL d_radius1;       // radius of osculating circles at point1
  REAL d_radius2;       // radius of osculating circles at point2
  Vec d_normalDirection; // normal direction, pointing from particle 1 to particle 2
  Vec d_tangentDirection;    // tangential direction

  bool d_isInContact; // are p1 and p1 in contact
  bool d_tangentLoadingActive;  // tangential loading or unloading
  Vec d_normalForce;  // pointing from particle 2 to paticle 1
  Vec
    d_tangentForce; // TangentrDirc points along tangential forces exerted on particle 1
  Vec d_tangentDisplacement; // tangential relative displacment total vector
  Vec
    d_tangentDisplacementStart; // displacement start value for each loading-unloading loop
  bool d_tangentSlidingActive;  // tangential silde or not

  bool d_prevTangentLoadingActive; // previous loading-unloading status
  Vec d_prevNormalForce;
  Vec d_prevTangentForce;
  Vec d_prevTangentDisplacement; // previous tangential relative displacment total vector
  bool d_prevTangentSlidingActive;
  REAL d_tangentForcePeak;

  Vec d_cohesionForce; // cohesion force between particles
  Vec d_spinResist;

  REAL d_E0;
  REAL d_G0;
  REAL d_R0;
  REAL d_vibraTimeStep;
  REAL d_impactTimeStep;

  void computeTangentForceMindlinAssumed(const REAL& contactFric,
                                         const REAL& poisson,
                                         const Vec& tangentDispInc);
  void computeTangentForceMindlinKnown(const REAL& contactFric,
                                       const REAL& poisson,
                                       const Vec& tangentDispInc);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_p1;
    ar& d_p2;
    ar& d_penetration;
    ar& d_contactRadius;
    ar& d_point1;
    ar& d_point2;
    ar& d_radius1;
    ar& d_radius2;
    ar& d_normalDirection;
    ar& d_tangentDirection;
    ar& d_isInContact;
    ar& d_tangentLoadingActive;
    ar& d_normalForce;
    ar& d_tangentForce;
    ar& d_tangentDisplacement;
    ar& d_tangentDisplacementStart;
    ar& d_tangentSlidingActive;
    ar& d_prevTangentLoadingActive;
    ar& d_prevNormalForce;
    ar& d_prevTangentForce;
    ar& d_prevTangentDisplacement;
    ar& d_prevTangentSlidingActive;
    ar& d_tangentForcePeak;
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
