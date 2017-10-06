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
    , direc(0)
    , point(0)
    , prevPoint(0)
    , veloc(0)
    , prevVeloc(0)
    , radius(0)
  {
  }

  CylinderBoundary(Boundary::BoundaryType type, std::ifstream& ifs);

  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const XMLProblemSpec& ps);
  CylinderBoundary(Boundary::BoundaryType tp, BoundaryID id, const JsonProblemSpec& ps);

  Vec getDirec() const { return direc; }
  Vec getPoint() const override { return point; }
  Vec getVeloc() const override { return veloc; }
  Vec getPrevPoint() const override { return prevPoint; }
  Vec getPrevVeloc() const override { return prevVeloc; }
  REAL getRadius() const { return radius; }

  void setDirec(Vec dir) { direc = dir; }
  void setPoint(Vec pnt) override { point = pnt; }
  void setVeloc(Vec vel) override { veloc = vel; }

  void print(std::ostream& os) override
  {
    Boundary::print(os);
    os << std::setw(OWID) << direc.x() << std::setw(OWID) << direc.y()
       << std::setw(OWID) << direc.z() << std::setw(OWID) << point.x()
       << std::setw(OWID) << point.y() << std::setw(OWID) << point.z()
       << std::setw(OWID) << radius << std::endl
       << std::endl;
  }

  void findBdryContact(DEMParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap) override;

private:
  Vec direc;
  Vec point;
  Vec prevPoint;
  Vec veloc;
  Vec prevVeloc;
  REAL radius;

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& direc;
    ar& point;
    ar& radius;
  }

};

} // namespace dem ends

#endif
