#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/Containers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/Containers.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>
#include <vector>

namespace dem {

class Particle; // forward declaration, only use pointer to class Particle

///////////////////////////////////////
class Boundary
{ // abstract base class
protected:
  std::size_t id;
  std::size_t type;

  // extra edges that are necessary to define a finite plane
  // e.g., side wall of a top-open container
  std::size_t extraNum;
  std::vector<Plane> extraEdge;

  ParticlePArray possParticle;
  BoundaryContactArray contactInfo;
  std::size_t contactNum;
  Vec normal;
  Vec tangt;
  REAL penetr;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& id;
    ar& type;
    ar& extraNum;
    ar& extraEdge;
    ar& possParticle;
    ar& contactInfo;
    ar& contactNum;
    ar& normal;
    ar& tangt;
    ar& penetr;
  }

public:
  Boundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0)
    : id(i)
    , type(tp)
    , extraNum(en)
    , contactNum(0)
    , normal(0)
    , tangt(0)
    , penetr(0)
  {
  }

  Boundary(std::size_t type, std::ifstream& ifs);
  virtual ~Boundary() = default; // polymorphic base class requires a virtual destructor

  std::size_t getId() { return id; }
  std::size_t getType() { return type; }
  ParticlePArray& getPossParticle() { return possParticle; }
  BoundaryContactArray& getContactInfo() { return contactInfo; }
  std::size_t getContactNum() const { return contactNum; }
  Vec getNormalForce() const { return normal; }
  Vec getTangtForce() const { return tangt; }
  REAL getAvgPenetr() const { return penetr; }

  virtual void print(std::ostream& os);
  virtual void printContactInfo(std::ostream& os);
  virtual void findBdryContact(ParticlePArray& ptcls) = 0;
  virtual void boundaryForce(BoundaryTangentArrayMap& boundaryTgtMap) = 0;
  virtual void updateStatForce();
  void clearStatForce();
  void clearContactInfo();

  virtual void updateIsotropic(REAL simga, REAL areaX, REAL areaY, REAL areaZ)
  {
  }
  virtual void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
  virtual void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
  virtual void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY, REAL areaZ)
  {
  }
  virtual void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY,
                                  REAL areaZ, REAL sigmaX, REAL sigmaY)
  {
  }
  virtual Vec getPoint() const = 0;
  virtual Vec getVeloc() const = 0;
  virtual Vec getPrevPoint() const = 0;
  virtual Vec getPrevVeloc() const = 0;
  virtual void setPoint(Vec pnt) = 0;
  virtual void setVeloc(Vec vel) = 0;
};

} // namespace dem ends

#endif
