#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/DEMContainers.h>
#include <InputOutput/json/json.hpp>
#include <InputOutput/zenxml/xml.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>
#include <vector>

namespace dem {

using BoundaryId = std::size_t;
using BoundaryType = std::size_t;
using EdgeCount = std::size_t;
using EdgeArray = std::vector<Plane>;
using ContactCount = std::size_t;
using XMLProblemSpec = zen::XmlIn;
using JsonProblemSpec = nlohmann::json;

class DEMParticle; // forward declaration, only use pointer to class DEMParticle

///////////////////////////////////////
class Boundary
{ // abstract base class
protected:
  BoundaryId b_id;
  BoundaryType b_type;

  // extra edges that are necessary to define a finite plane
  // e.g., side wall of a top-open container
  EdgeCount b_extraNum;
  EdgeArray b_extraEdge;

  ParticlePArray possParticle;
  BoundaryContactArray contactInfo;
  ContactCount contactNum;
  Vec normal;
  Vec tangt;
  REAL penetr;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& b_id;
    ar& b_type;
    ar& b_extraNum;
    ar& b_extraEdge;
    ar& possParticle;
    ar& contactInfo;
    ar& contactNum;
    ar& normal;
    ar& tangt;
    ar& penetr;
  }

public:
  Boundary(BoundaryId id, BoundaryType tp = 0, EdgeCount en = 0)
    : b_id(id)
    , b_type(tp)
    , b_extraNum(en)
    , contactNum(0)
    , normal(0)
    , tangt(0)
    , penetr(0)
  {
  }

  Boundary();
  virtual ~Boundary() =
    default; // polymorphic base class requires a virtual destructor

  BoundaryId getId() { return b_id; }
  BoundaryType getType() { return b_type; }
  ParticlePArray& getPossParticle() { return possParticle; }
  BoundaryContactArray& getContactInfo() { return contactInfo; }
  ContactCount getContactNum() const { return contactNum; }
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
