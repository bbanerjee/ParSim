#ifndef PLANE_BOUNDARY_H
#define PLANE_BOUNDARY_H

#include <Boundary/Boundary.h>
#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>

namespace dem {

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

///////////////////////////////////////
class PlaneBoundary : public Boundary
{
private:
  Vec direc;
  Vec point;
  Vec prevPoint;
  Vec veloc;
  Vec prevVeloc;

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
  PlaneBoundary(BoundaryId id = 0, BoundaryType tp = 0, EdgeCount en = 0)
    : Boundary(id, tp, en)
    , direc(0)
    , point(0)
    , prevPoint(0)
    , veloc(0)
    , prevVeloc(0)
  {
  }

  PlaneBoundary(BoundaryType type, std::ifstream& ifs);
  PlaneBoundary(BoundaryId id, BoundaryType type, const XMLProblemSpec& ps);
  PlaneBoundary(BoundaryId id, BoundaryType type, const JsonProblemSpec& ps);

  Vec getDirec() const { return direc; }
  Vec getPoint() const override { return point; }
  Vec getVeloc() const override { return veloc; }
  Vec getPrevPoint() const override { return prevPoint; }
  Vec getPrevVeloc() const override { return prevVeloc; }

  void setDirec(Vec dir) { direc = dir; }
  void setPoint(Vec pnt) override { point = pnt; }
  void setVeloc(Vec vel) override { veloc = vel; }

  REAL distanceToBdry(Vec pos) const
  {
    return dot((pos - point) , normalize(direc));
  }
  REAL distanceToBdry(Vec pos, Plane pn) const
  {
    return dot(pos - pn.getPoint() , normalize(pn.getDirec()));
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
  void findBdryContact(ParticlePArray& ptcls) override;
  void boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap) override;
};

} // namespace dem ends

#endif
