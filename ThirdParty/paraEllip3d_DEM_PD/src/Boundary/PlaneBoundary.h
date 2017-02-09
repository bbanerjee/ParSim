#ifndef PLANE_BOUNDARY_H
#define PLANE_BOUNDARY_H

#include <Boundary/Boundary.h>
#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/Containers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>

namespace dem {

class Particle; // forward declaration, only use pointer to class Particle

///////////////////////////////////////
class PlaneBoundary : public Boundary
{
private:
  Vec direc;
  Vec point;
  Vec prevPoint;
  Vec veloc;
  Vec prevVeloc;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& boost::serialization::base_object<Boundary>(*this);
    ar& direc;
    ar& point;
    ar& prevPoint;
    ar& veloc;
    ar& prevVeloc;
  }

public:
  PlaneBoundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0)
    : Boundary(i, tp, en)
    , direc(0)
    , point(0)
    , prevPoint(0)
    , veloc(0)
    , prevVeloc(0)
  {
  }

  PlaneBoundary(std::size_t type, std::ifstream& ifs);

  Vec getDirec() const { return direc; }
  Vec getPoint() const { return point; }
  Vec getVeloc() const { return veloc; }
  Vec getPrevPoint() const { return prevPoint; }
  Vec getPrevVeloc() const { return prevVeloc; }

  void setDirec(Vec dir) { direc = dir; }
  void setPoint(Vec pnt) { point = pnt; }
  void setVeloc(Vec vel) { veloc = vel; }

  REAL distanceToBdry(Vec pos) const
  {
    return (pos - point) * normalize(direc);
  }
  REAL distanceToBdry(Vec pos, Plane pn) const
  {
    return (pos - pn.getPoint()) * normalize(pn.getDirec());
  }

  void print(std::ostream& os);
  void printContactInfo(std::ostream& os);

  void updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ);
  void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ,
                          REAL sigmaX, REAL sigmaY);
  void findBdryContact(ParticlePArray& ptcls);
  void boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap);
};

} // namespace dem ends

#endif
