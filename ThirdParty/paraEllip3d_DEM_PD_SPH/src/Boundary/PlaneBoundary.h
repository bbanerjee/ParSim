#ifndef PLANE_BOUNDARY_H
#define PLANE_BOUNDARY_H

#include <Boundary/Boundary.h>
#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

///////////////////////////////////////
class PlaneBoundary : public Boundary
{
public:

  PlaneBoundary(BoundaryID id = BoundaryID::NONE, 
                Boundary::BoundaryType tp = Boundary::BoundaryType::NONE, 
                EdgeCount en = 0)
    : Boundary(id, tp, en)
    , d_direction(0)
    , d_position(0)
    , d_previousPosition(0)
    , d_velocity(0)
    , d_previousVelocity(0)
  {
  }

  PlaneBoundary(Boundary::BoundaryType type, std::ifstream& ifs);
  PlaneBoundary(Boundary::BoundaryType type, BoundaryID id, const XMLProblemSpec& ps);
  PlaneBoundary(Boundary::BoundaryType type, BoundaryID id, const JsonProblemSpec& ps);

  Vec getDirection() const { return d_direction; }
  Vec getPosition() const override { return d_position; }
  Vec getVelocity() const override { return d_velocity; }
  Vec getPreviousPosition() const override { return d_previousPosition; }
  Vec getPreviousVelocity() const override { return d_previousVelocity; }

  void setDirection(Vec dir) { d_direction = dir; }
  void setPosition(Vec pnt) override { d_position = pnt; }
  void setVelocity(Vec vel) override { d_velocity = vel; }

  REAL distanceToBdry(Vec pos) const
  {
    return dot((pos - d_position) , normalize(d_direction));
  }
  REAL distanceToBdry(Vec pos, Plane pn) const
  {
    return dot(pos - pn.getPosition() , normalize(pn.getDirection()));
  }

  void print(std::ostream& os) override;
  void printContactInfo(std::ostream& os) override;

  void updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ) override;
  void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ) override;
  void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ) override;
  void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY,
                         REAL areaZ) override;
  void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ,
                          REAL sigmaX, REAL sigmaY) override;
  void findBdryContact(DEMParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap) override;

private:
  Vec d_direction;
  Vec d_position;
  Vec d_previousPosition;
  Vec d_velocity;
  Vec d_previousVelocity;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& d_direction;
    ar& d_position;
    ar& d_previousPosition;
    ar& d_velocity;
    ar& d_previousVelocity;
  }

};

} // namespace dem ends

#endif
