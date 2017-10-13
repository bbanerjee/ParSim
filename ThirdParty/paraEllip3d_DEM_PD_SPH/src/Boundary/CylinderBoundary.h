#ifndef CYLINDER_BOUNDARY_H
#define CYLINDER_BOUNDARY_H

#include <Boundary/Boundary.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Cylinder.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

///////////////////////////////////////
class CylinderBoundary : public Boundary
{

public:
  CylinderBoundary()
    : Boundary()
    , direction(0)
    , point(0)
    , prevPoint(0)
    , velocity(0)
    , prevVelocity(0)
    , radius(0)
  {
  }

  CylinderBoundary(Boundary::BoundaryType type, std::ifstream& ifs);

  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const XMLProblemSpec& ps);
  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const JsonProblemSpec& ps);

  Vec getDirection() const { return direction; }
  Vec getPosition() const override { return point; }
  Vec getVelocity() const override { return velocity; }
  Vec getPreviousPosition() const override { return prevPoint; }
  Vec getPreviousVelocity() const override { return prevVelocity; }
  REAL getRadius() const { return radius; }

  void setDirection(Vec dir) { direction = dir; }
  void setPosition(Vec pnt) override { point = pnt; }
  void setVelocity(Vec vel) override { velocity = vel; }

  void print(std::ostream& os) override
  {
    Boundary::print(os);
    os << std::setw(OWID) << direction.x() << std::setw(OWID) << direction.y()
       << std::setw(OWID) << direction.z() << std::setw(OWID) << point.x()
       << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
       << std::setw(OWID) << radius << std::endl
       << std::endl;
  }

  void findBdryContact(DEMParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTangentMap) override;

private:
  Vec direction;
  Vec point;
  Vec prevPoint;
  Vec velocity;
  Vec prevVelocity;
  REAL radius;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& direction;
    ar& point;
    ar& radius;
  }

};

} // namespace dem ends

#endif
