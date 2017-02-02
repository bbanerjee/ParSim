#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <Core/Types/realtypes.h>
#include <Core/Math/Vec.h>
#include <Core/Geometry/Cylinder.h>
#include <DiscreteElements/Containers.h>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstddef>
#include <cstdlib>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace dem {

class Particle; // forward declaration, only use pointer to class Particle

/////////////////////////////////////
class BoundaryTgt {
public:
  std::size_t particleId;
  Vec tgtForce;
  Vec tgtDisp;
  bool tgtLoading;
  Vec tgtDispStart;
  REAL tgtPeak;

public:
  BoundaryTgt()
      : particleId(0), tgtForce(0), tgtDisp(0), tgtLoading(false),
        tgtDispStart(0), tgtPeak(0) {}

  BoundaryTgt(std::size_t _particleId, Vec _v1, Vec _v2, bool _b, Vec _v3,
              REAL _tp)
      : particleId(_particleId), tgtForce(_v1), tgtDisp(_v2), tgtLoading(_b),
        tgtDispStart(_v3), tgtPeak(_tp) {}

};

/////////////////////////////////////
class BdryContact {
public:
  Particle *ptcl;
  Vec point;
  Vec normal;
  Vec tangt;
  REAL penetr;

public:
  BdryContact() : ptcl(NULL), point(0), normal(0), tangt(0), penetr(0) {}

  BdryContact(Particle *p, Vec pt, Vec nm, Vec tg, REAL pntr)
      : ptcl(p), point(pt), normal(nm), tangt(tg), penetr(pntr) {}

  void print(std::ostream &os) {
    os << std::setw(OWID) << point.getX() << std::setw(OWID) << point.getY()
       << std::setw(OWID) << point.getZ() << std::setw(OWID) << normal.getX()
       << std::setw(OWID) << normal.getY() << std::setw(OWID) << normal.getZ()
       << std::setw(OWID) << tangt.getX() << std::setw(OWID) << tangt.getY()
       << std::setw(OWID) << tangt.getZ() << std::setw(OWID) << penetr
       << std::endl;
  }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &ptcl;
    ar &point;
    ar &normal;
    ar &tangt;
    ar &penetr;
  }
};

/////////////////////////////////////
class Plane {
public:
  Vec direc;
  Vec point;

public:
  Plane() : direc(0), point(0) {}

  Plane(Vec dir, Vec pt) : direc(dir), point(pt) {}

  void print(std::ostream &os) {
    os << std::setw(OWID) << direc.getX() << std::setw(OWID) << direc.getY()
       << std::setw(OWID) << direc.getZ() << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY() << std::setw(OWID) << point.getZ()
       << std::endl;
  }

  Vec getDirec() const { return direc; }
  Vec getPoint() const { return point; }

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &direc;
    ar &point;
  }
};

///////////////////////////////////////
class Boundary { // abstract base class
protected:
  std::size_t id;
  std::size_t type;

  // extra edges that are necessary to define a finite plane
  // e.g., side wall of a top-open container
  std::size_t extraNum;
  std::vector<Plane> extraEdge;

  ParticlePArray possParticle;
  std::vector<BdryContact> contactInfo;
  std::size_t contactNum;
  Vec normal;
  Vec tangt;
  REAL penetr;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &id;
    ar &type;
    ar &extraNum;
    ar &extraEdge;
    ar &possParticle;
    ar &contactInfo;
    ar &contactNum;
    ar &normal;
    ar &tangt;
    ar &penetr;
  }

public:
  Boundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0)
      : id(i), type(tp), extraNum(en), contactNum(0), normal(0), tangt(0),
        penetr(0) {}

  Boundary(std::size_t type, std::ifstream &ifs);
  virtual ~Boundary() {
  } // polymorphic base class requires a virtual destructor

  std::size_t getId() { return id; }
  std::size_t getType() { return type; }
  ParticlePArray& getPossParticle() { return possParticle; }
  std::vector<BdryContact> &getContactInfo() { return contactInfo; }
  std::size_t getContactNum() const { return contactNum; }
  Vec getNormalForce() const { return normal; }
  Vec getTangtForce() const { return tangt; }
  REAL getAvgPenetr() const { return penetr; }

  virtual void print(std::ostream &os);
  virtual void printContactInfo(std::ostream &os);
  virtual void findBdryContact(ParticlePArray &ptcls) = 0;
  virtual void boundaryForce(
      std::map<std::size_t, std::vector<BoundaryTgt> > &boundaryTgtMap) = 0;
  virtual void updateStatForce();
  void clearStatForce();
  void clearContactInfo();

  virtual void updateIsotropic(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {
  }
  virtual void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
  virtual void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ) {}
  virtual void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY,
                                 REAL areaZ) {}
  virtual void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY,
                                  REAL areaZ, REAL sigmaX, REAL sigmaY) {}
  virtual Vec getPoint() const = 0;
  virtual Vec getVeloc() const = 0;
  virtual Vec getPrevPoint() const = 0;
  virtual Vec getPrevVeloc() const = 0;
  virtual void setPoint(Vec pnt) = 0;
  virtual void setVeloc(Vec vel) = 0;
};

///////////////////////////////////////
class planeBoundary : public Boundary {
private:
  Vec direc;
  Vec point;
  Vec prevPoint;
  Vec veloc;
  Vec prevVeloc;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<Boundary>(*this);
    ar &direc;
    ar &point;
    ar &prevPoint;
    ar &veloc;
    ar &prevVeloc;
  }

public:
  planeBoundary(std::size_t i = 0, std::size_t tp = 0, std::size_t en = 0)
      : Boundary(i, tp, en), direc(0), point(0), prevPoint(0), veloc(0),
        prevVeloc(0) {}

  planeBoundary(std::size_t type, std::ifstream &ifs);

  Vec getDirec() const { return direc; }
  Vec getPoint() const { return point; }
  Vec getVeloc() const { return veloc; }
  Vec getPrevPoint() const { return prevPoint; }
  Vec getPrevVeloc() const { return prevVeloc; }

  void setDirec(Vec dir) { direc = dir; }
  void setPoint(Vec pnt) { point = pnt; }
  void setVeloc(Vec vel) { veloc = vel; }

  REAL distanceToBdry(Vec pos) const {
    return (pos - point) * normalize(direc);
  }
  REAL distanceToBdry(Vec pos, Plane pn) const {
    return (pos - pn.getPoint()) * normalize(pn.getDirec());
  }

  void print(std::ostream &os);
  void printContactInfo(std::ostream &os);

  void updateIsotropic(REAL sigma, REAL areaX, REAL areaY, REAL areaZ);
  void updateOdometer(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updateTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updatePlaneStrain(REAL simga, REAL areaX, REAL areaY, REAL areaZ);
  void updateTrueTriaxial(REAL simga, REAL areaX, REAL areaY, REAL areaZ,
                          REAL sigmaX, REAL sigmaY);
  void findBdryContact(ParticlePArray &ptcls);
  void boundaryForce(
      std::map<std::size_t, std::vector<BoundaryTgt> > &boundaryTgtMap);

};

///////////////////////////////////////
class cylinderBoundary : public Boundary {
private:
  Vec direc;
  Vec point;
  Vec prevPoint;
  Vec veloc;
  Vec prevVeloc;
  REAL radius;

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &boost::serialization::base_object<Boundary>(*this);
    ar &direc;
    ar &point;
    ar &radius;
  }

public:
  cylinderBoundary()
      : Boundary(), direc(0), point(0), prevPoint(0), veloc(0), prevVeloc(0),
        radius(0) {}

  cylinderBoundary(std::size_t type, std::ifstream &ifs);

  Vec getDirec() const { return direc; }
  Vec getPoint() const { return point; }
  Vec getVeloc() const { return veloc; }
  Vec getPrevPoint() const { return prevPoint; }
  Vec getPrevVeloc() const { return prevVeloc; }
  REAL getRadius() const { return radius; }

  void setDirec(Vec dir) { direc = dir; }
  void setPoint(Vec pnt) { point = pnt; }
  void setVeloc(Vec vel) { veloc = vel; }

  void print(std::ostream &os) {
    Boundary::print(os);
    os << std::setw(OWID) << direc.getX() << std::setw(OWID) << direc.getY()
       << std::setw(OWID) << direc.getZ() << std::setw(OWID) << point.getX()
       << std::setw(OWID) << point.getY() << std::setw(OWID) << point.getZ()
       << std::setw(OWID) << radius << std::endl << std::endl;
  }

  void findBdryContact(ParticlePArray &ptcls);
  void boundaryForce(
      std::map<std::size_t, std::vector<BoundaryTgt> > &boundaryTgtMap);

};

} // namespace dem ends

#endif
