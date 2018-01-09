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

struct DEMContactData {

public:

  DEMContactData();
  ~DEMContactData() = default;

  void write(zen::XmlOut& xml) const;
  void read(const zen::XmlIn& xml); 

  Vec point1;                   // point1 on particle 1, innermost to particle 2
  Vec point2;                   // point2 on particle 2, innermost to particle 1
  Vec normalForce;              // pointing from particle 2 to paticle 1
  Vec prevNormalForce;
  Vec tangentForce;             // TangentrDirc points along tangential forces 
                                // exerted on particle 1
  Vec prevTangentForce;
  Vec tangentDisplacement;      // tangential relative displacment total vector

  Vec tangentDisplacementStart; // displacement start value for each 
                                // loading-unloading loop
  Vec prevTangentDisplacement;  // previous tangential relative displacment 
                                // total vector
  Vec normalDirection;          // normal direction, pointing from particle 1 
                                // to particle 2
  Vec tangentDirection;         // tangential direction
  Vec cohesionForce;            // cohesion force between particles
  Vec spinResist;

  REAL radius1;            // radius of osculating circles at point1
  REAL radius2;            // radius of osculating circles at point2
  REAL penetration;        // penetration
  REAL contactRadius;      // radius of contact surface
  REAL R0;
  REAL E0;
  REAL G0;
  REAL tangentForcePeak;
  REAL vibrationTimeStep;
  REAL impactTimeStep;

  bool isInContact;              // are p1 and p1 in contact
  bool tangentLoadingActive;     // tangential loading or unloading
  bool tangentSlidingActive;     // tangential silde or not
  bool prevTangentLoadingActive; // previous loading-unloading status
  bool prevTangentSlidingActive;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& penetration;
    ar& contactRadius;
    ar& point1;
    ar& point2;
    ar& radius1;
    ar& radius2;
    ar& normalDirection;
    ar& tangentDirection;
    ar& isInContact;
    ar& tangentLoadingActive;
    ar& normalForce;
    ar& tangentForce;
    ar& tangentDisplacement;
    ar& tangentDisplacementStart;
    ar& tangentSlidingActive;
    ar& prevTangentLoadingActive;
    ar& prevNormalForce;
    ar& prevTangentForce;
    ar& prevTangentDisplacement;
    ar& prevTangentSlidingActive;
    ar& tangentForcePeak;
    ar& cohesionForce;
    ar& spinResist;
    ar& E0;
    ar& G0;
    ar& R0;
    ar& vibrationTimeStep;
    ar& impactTimeStep;
  }

private:

  void writeBoolean(const std::string& label,
                    bool variable,
                    zen::XmlOut& xml) const;

  void writeScalar(const std::string& label,
                   double variable,
                   zen::XmlOut& xml) const;

  void writeVector(const std::string& label,
                   const Vec& variable,
                   zen::XmlOut& xml) const;

  template <typename T>
  bool readValue(const zen::XmlIn& ps, 
                 const std::string& label,
                 T& value);
};

class DEMContact
{
public:
  DEMContact();
  DEMContact(DEMParticle* t1, DEMParticle* t2);
  DEMContact(DEMParticle* t1, DEMParticle* t2, const DEMContactData& data);

  DEMParticle* getP1() const;
  DEMParticle* getP2() const;
  Vec getPoint1() const { return d_data.point1; }
  Vec getPoint2() const { return d_data.point2; }
  REAL radius1() const { return d_data.radius1; }
  REAL radius2() const { return d_data.radius2; }
  REAL getR0() const { return d_data.R0; }
  REAL getE0() const { return d_data.E0; }
  REAL getVibrationTimeStep() const { return d_data.vibrationTimeStep; }
  REAL getImpactTimeStep() const { return d_data.impactTimeStep; }

  bool isOverlapped(REAL minRelativeOverlap, REAL measurableOverlap,
                    std::size_t iteration);

  // calculate normal and tangential force of contact
  int computeContactForces(REAL timestep, std::size_t iteration,
                           REAL stiffness, REAL shearModulus,
                           REAL cohesion, REAL damping,
                           REAL friction, REAL maxOverlapFactor,
                           REAL minMeasurableOverlap);

  inline
  Vec computeVelocity(DEMParticle* particle, const Vec& midPoint) {
    Vec velocity = particle->currentVelocity() + 
      cross(particle->currentAngularVelocity(), 
            midPoint - particle->currentPosition());
    return velocity;
  }

  inline 
  REAL computeNormalStiffness(const Vec& normalForce, REAL R, REAL E) {
    REAL kn = std::pow(6 * vnormL2(normalForce) * R * E * E, 1.0 / 3.0);
    return kn;
  }

  Vec computeDampingForce(REAL pMass1, REAL pMass2, const Vec& relativeVel,
                          REAL dampingCoeff, 
                          REAL normalStiffness, const Vec& contactNormal);

  void updateTimestep(REAL pMass1, REAL pMass2, const Vec& relativeVel,
                      REAL normalStiffness, const Vec& contactNormal,
                      REAL allowedOverlap);

  void computeTangentForce(REAL timestep, std::size_t iteration,
                           REAL shearModulus, REAL poissonRatio, REAL friction, 
                           REAL contactRadius, const Vec& relativeVel, 
                           const Vec& contactNormal, const Vec& normalForce);

  // Update the forces and moments
  void updateForceAndMoment(const Vec& normalForce, const Vec& tangentForce,
                            const Vec& cohesionForce, const Vec& dampingForce,
                            const Vec& momentCenter,
                            DEMParticle* p1, DEMParticle* p2);

  REAL getNormalForceMagnitude() const { return vnormL2(d_data.normalForce); }
  REAL getTangentForceMagnitude() const { return vnormL2(d_data.tangentForce); }
  REAL getPenetration() const { return d_data.penetration; }
  REAL getContactRadius() const { return d_data.contactRadius; }
  REAL getTangentDisplacement() const
  {
    return vnormL2(d_data.tangentDisplacement);
  } // total value during a process of contact
  Vec getTangentDisplacementVec() const
  {
    return d_data.tangentDisplacement;
  } 

  void checkoutContactTangents(std::vector<DEMContactTangent>& contactTangentVec);
  void checkinPreviousContactTangents(std::vector<DEMContactTangent>& contactTangentVec);
  Vec getNormalForce() const { return d_data.normalForce; }
  Vec getTangentForce() const { return d_data.tangentForce; }
  bool isRedundant(const DEMContact& other) const;
  bool operator==(const DEMContact& other) const;

  void write(std::stringstream& str) const;
  void write(zen::XmlOut& xml) const;

private:
  DEMParticle* d_p1;            // particle 1
  DEMParticle* d_p2;            // particle 2
  DEMContactData d_data;

  void computeTangentForceMindlinAssumed(const REAL& contactFric,
                                         const REAL& poisson,
                                         const Vec& tangentDispInc);
  void computeTangentForceMindlinKnown(const REAL& contactFric,
                                       const REAL& poisson,
                                       const Vec& tangentDispInc,
                                       std::size_t iteration);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_p1;
    ar& d_p2;
    ar& d_data;
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
