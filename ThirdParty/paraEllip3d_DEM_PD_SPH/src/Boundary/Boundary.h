#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <Boundary/BoundaryContact.h>
#include <Boundary/BoundaryTangent.h>
#include <Boundary/BoundaryContainers.h>
#include <Core/Geometry/Plane.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <DiscreteElements/DEMContainers.h>
#include <InputOutput/json/json.hpp>
#include <InputOutput/zenxml/xml.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <iostream>
#include <utility>
#include <vector>

namespace dem {

using EdgeCount = std::size_t;
using EdgeArray = std::vector<Plane>;
using ContactCount = std::size_t;
using XMLProblemSpec = zen::XmlIn;
using JsonProblemSpec = nlohmann::json;

class DEMParticle; 

class Boundary
{ 
public:

  enum class BoundaryType
  {
    NONE = 0,
    PLANE = 1,
    CYLINDER = 2
  };

  enum class BoundaryID
  {
    NONE = 0,
    XMINUS = 1,
    XPLUS  = 2,
    YMINUS = 3,
    YPLUS  = 4,
    ZMINUS = 5,
    ZPLUS  = 6
  };

  enum BoundaryFlag
  {
    ONLY_BOTTOM_BOUNDARY,
    NO_TOP_BOUNDARY,
    ALL_BOUNDARIES
  };

  static BoundaryType getBoundaryType(const std::string& str)
  {
    if (str == "plane")
      return BoundaryType::PLANE;
    else if (str == "cylinder")
      return BoundaryType::CYLINDER;
    else
      return BoundaryType::NONE;
  }

  static BoundaryType getBoundaryType(int type)
  {
    return static_cast<BoundaryType>(type);
  }

  static BoundaryID getBoundaryID(const std::string& str)
  {
    if (str == "x-")
      return BoundaryID::XMINUS;
    else if (str == "x+")
      return BoundaryID::XPLUS;
    else if (str == "y-")
      return BoundaryID::YMINUS;
    else if (str == "y+")
      return BoundaryID::YPLUS;
    else if (str == "z-")
      return BoundaryID::ZMINUS;
    else if (str == "z+")
      return BoundaryID::ZPLUS;
    else
      return BoundaryID::NONE;
  }

  static BoundaryID getBoundaryID(int id)
  {
    return static_cast<BoundaryID>(id);
  }

  static std::string getBoundaryIDStr(BoundaryID id)
  {
    if (id == BoundaryID::XMINUS)
     return "x-";
    else if (id == BoundaryID::XPLUS)
     return "x+";
    else if (id == BoundaryID::YMINUS)
     return "y-";
    else if (id == BoundaryID::YPLUS)
     return "y+";
    else if (id == BoundaryID::ZMINUS)
     return "z-";
    else if (id == BoundaryID::ZPLUS)
     return "z+";
    else
     return "none";
  }

  static BoundaryFlag getBoundaryFlag(std::size_t boundaryNum)
  {
    if (boundaryNum == 1) 
      return BoundaryFlag::ONLY_BOTTOM_BOUNDARY;
    else if (boundaryNum == 5)
      return BoundaryFlag::NO_TOP_BOUNDARY;
    else
      return BoundaryFlag::ALL_BOUNDARIES;
  }

public:

  Boundary(BoundaryID id, BoundaryType tp = BoundaryType::NONE, EdgeCount en = 0)
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
  virtual ~Boundary() = default;

  BoundaryID getId() { return b_id; }
  BoundaryType getType() { return b_type; }
  DEMParticlePArray& getPossParticle() { return possParticle; }
  BoundaryContactArray& getContactInfo() { return contactInfo; }
  ContactCount getContactNum() const { return contactNum; }
  Vec getNormalForce() const { return normal; }
  Vec getTangtForce() const { return tangt; }
  REAL getAvgPenetr() const { return penetr; }

  virtual void print(std::ostream& os);
  virtual void printContactInfo(std::ostream& os);
  virtual void findBdryContact(DEMParticlePArray& ptcls) = 0;
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

protected:

  BoundaryID b_id;
  BoundaryType b_type;

  // extra edges that are necessary to define a finite plane
  // e.g., side wall of a top-open container
  EdgeCount b_extraNum;
  EdgeArray b_extraEdge;

  DEMParticlePArray possParticle;
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

};

} // namespace dem ends

#endif
