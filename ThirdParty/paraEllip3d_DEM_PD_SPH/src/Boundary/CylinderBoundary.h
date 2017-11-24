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
    , d_direction(0)
    , d_point(0)
    , d_prevPoint(0)
    , d_velocity(0)
    , d_prevVelocity(0)
    , d_radius(0)
  {
  }

  CylinderBoundary(Boundary::BoundaryType type, std::ifstream& ifs);

  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const XMLProblemSpec& ps);
  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const JsonProblemSpec& ps);

  Vec getDirection() const { return d_direction; }
  Vec getPosition() const override { return d_point; }
  Vec getVelocity() const override { return d_velocity; }
  Vec previousPosition() const override { return d_prevPoint; }
  Vec previousVelocity() const override { return d_prevVelocity; }
  REAL radius() const { return d_radius; }

  void setDirection(Vec dir) { d_direction = dir; }
  void setPosition(Vec pnt) override { d_point = pnt; }
  void setVelocity(Vec vel) override { d_velocity = vel; }

  void print(std::ostream& os) override
  {
    Boundary::print(os);
    os << std::setw(OWID) << d_direction.x() << std::setw(OWID) << d_direction.y()
       << std::setw(OWID) << d_direction.z() << std::setw(OWID) << d_point.x()
       << std::setw(OWID) << d_point.y() << std::setw(OWID) << d_point.z()
       << std::setw(OWID) << d_radius << std::endl
       << std::endl;
  }

  void findBoundaryContacts(DEMParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTangentMap,
                     std::size_t iteration) override;

private:
  Vec d_direction;
  Vec d_point;
  Vec d_prevPoint;
  Vec d_velocity;
  Vec d_prevVelocity;
  REAL d_radius;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& d_direction;
    ar& d_point;
    ar& d_radius;
  }

};

} // namespace dem ends

#endif
