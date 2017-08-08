#ifndef PARTICLE_H
#define PARTICLE_H

#include <Boundary/Boundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/Cylinder.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/Gradation.h>
#include <InputOutput/InputParameter.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <cstddef>
#include <map>
#include <vector>

namespace dem {

class DEMParticle
{

private:
  // types of individual particle:
  //   0 - free particle
  //   1 - fixed particle
  //   2 - special case 2 (pure moment): translate first, then rotate only,
  // MNT_START needs to be defined
  //   3 - special case 3 (displacemental ellipsoidal pile): translate in
  // vertical direction only
  //   4 - special case 4 (impacting ellipsoidal penetrator): impact with inital
  // velocity in vertical direction only
  //   5 - free boundary particle
  //   6 - translate only, no rotation
  //  10 - ghost particle
  std::size_t d_id;
  std::size_t d_type;
  REAL d_a, d_b, d_c; // three semi-axle length, must satisfy a >= b >= c
  REAL d_young; // note: a(currDirecA), b(currDirecB), c(currDirecC) corresponds
                // to x, y, z in local frame, respectively
  REAL d_poisson;
  Vec d_currPos; // particle center
  Vec d_prevPos;
  Vec d_currDirecA, d_currDirecB,
    d_currDirecC; // direction of the three axles, in radian
  Vec d_prevDirecA, d_prevDirecB, d_prevDirecC;
  Vec d_currVeloc; // the velocity of the mass center
  Vec d_prevVeloc;
  Vec d_currOmga; // angular velocity in global frame!
  Vec d_prevOmga;
  Vec d_force;
  std::map<size_t, Vec> d_forceIDMap;
  Vec d_prevForce;
  Vec d_moment;
  std::map<size_t, Vec> d_momentIDMap;
  Vec d_prevMoment;
  Vec d_constForce;
  Vec d_constMoment;
  REAL d_density; // specific gravity
  REAL d_mass;
  REAL d_volume;
  Vec d_momentJ;      // moment of inertia in local body-fixed frame
  REAL d_coef[10];    // particle's coefficients in global coordinates
  REAL d_kinetEnergy; // kinetic energy
  std::size_t d_contactNum;
  bool d_inContact; // in contact with other particle or boundary
  std::vector<std::vector<REAL>> d_fluidGrid;

public:
  DEMParticle();
  DEMParticle(std::size_t n, std::size_t type, Vec center, REAL r, REAL young,
           REAL poisson);
  DEMParticle(std::size_t n, std::size_t type, Vec center, REAL a, REAL b, REAL c,
           REAL young, REAL poisson);
  DEMParticle(std::size_t n, std::size_t type, Vec center, Gradation& grad,
           REAL young, REAL poisson);
  DEMParticle(std::size_t n, std::size_t type, Vec dim, Vec position, Vec dirca,
           Vec dircb, Vec dircc, REAL young, REAL poisson);

  std::size_t getId() const { return d_id; }
  std::size_t getType() const { return d_type; }
  REAL getA() const { return d_a; }
  REAL getB() const { return d_b; }
  REAL getC() const { return d_c; }
  REAL getYoung() const { return d_young; }
  REAL getPoisson() const { return d_poisson; };
  REAL getVolume() const { return d_volume; }
  REAL getMass() const { return d_mass; }
  REAL getDensity() const { return d_density; }
  Vec currentPosition() const { return d_currPos; }
  Vec getPrevPosition() const { return d_prevPos; }
  Vec getCurrDirecA() const { return d_currDirecA; }
  Vec getCurrDirecB() const { return d_currDirecB; }
  Vec getCurrDirecC() const { return d_currDirecC; }
  Vec getPrevDirecA() const { return d_prevDirecA; }
  Vec getPrevDirecB() const { return d_prevDirecB; }
  Vec getPrevDirecC() const { return d_prevDirecC; }
  Vec currentVel() const { return d_currVeloc; }
  Vec getPrevVeloc() const { return d_prevVeloc; }
  Vec currentOmega() const { return d_currOmga; }
  Vec getPrevOmga() const { return d_prevOmga; }
  Vec getForce() const { return d_force; }
  std::map<size_t, Vec> getForceIDMap() const { return d_forceIDMap; }
  Vec getMoment() const { return d_moment; }
  std::map<size_t, Vec> getMomentIDMap() const { return d_momentIDMap; }
  Vec getAccel() const { return d_force / d_mass; }
  Vec getConstForce() const { return d_constForce; }
  Vec getConstMoment() const { return d_constMoment; }
  Vec getmomentJ() const { return d_momentJ; }
  bool isInContact() const { return d_inContact; }
  std::size_t getContactNum() const { return d_contactNum; }

  REAL getRadius(Vec v) const;
  REAL getTransEnergy() const;
  REAL getRotatEnergy() const;
  REAL getKinetEnergy() const;
  REAL getPotenEnergy(REAL ref) const;

  void setId(std::size_t n) { d_id = n; }
  void setType(std::size_t n) { d_type = n; }
  void setA(REAL dd) { d_a = dd; }
  void setB(REAL dd) { d_b = dd; }
  void setC(REAL dd) { d_c = dd; }
  void expand(REAL percent)
  {
    d_a *= (1 + percent);
    d_b *= (1 + percent);
    d_c *= (1 + percent);
  }
  void setCurrPos(Vec vv) { d_currPos = vv; }
  void setPrevPos(Vec vv) { d_prevPos = vv; }
  void setCurrDirecA(Vec vv) { d_currDirecA = vv; }
  void setCurrDirecB(Vec vv) { d_currDirecB = vv; }
  void setCurrDirecC(Vec vv) { d_currDirecC = vv; }
  void setPrevDirecA(Vec vv) { d_prevDirecA = vv; }
  void setPrevDirecB(Vec vv) { d_prevDirecB = vv; }
  void setPrevDirecC(Vec vv) { d_prevDirecC = vv; }
  void setCurrVeloc(Vec vv) { d_currVeloc = vv; }
  void setPrevVeloc(Vec vv) { d_prevVeloc = vv; }
  void setCurrOmega(Vec vv) { d_currOmga = vv; }
  void setPrevOmega(Vec vv) { d_prevOmga = vv; }
  void setForce(Vec vv)
  {
    d_force = vv;
  }
  void setMoment(Vec vv) { d_moment = vv; }
  void setConstForce(Vec vv) { d_constForce = vv; }
  void setConstMoment(Vec vv) { d_constMoment = vv; }
  void setmomentJ(Vec v) { d_momentJ = v; }
  void setMass(REAL d) { d_mass = d; }
  void setDensity(REAL dn) { d_density = dn; }
  void setInContact(bool value) { d_inContact = value; }
  void setContactNum(std::size_t num) { d_contactNum = num; }

  void clearContactForce();
  void addForce(Vec vv)
  {
    d_force += vv;
  }
  void addForceIDMap(Vec vv, size_t id) {
    d_forceIDMap[id] = vv;
  }
  void addMoment(Vec vv) { d_moment += vv; }
  void addMomentIDMap(Vec vv, size_t id) {
    d_momentIDMap[id] = vv;
  }
  void update();

  Vec globalToLocal(Vec input) const;
  Vec localToGlobal(Vec input) const;

  Vec globalToLocalPrev(Vec input) const; // based on previous step
  Vec localToGlobalPrev(Vec input) const;

  // update global coefficients in the following form based on
  // position/dimensions/orientations
  // a0 x^2 + a1 y^2 + a2 z^2 + a3 xy + a4 yz + a5 zx + a6 x + a7 y + a8 z + a9
  // = 0
  void globalCoef();
  void getGlobalCoef(REAL coef[]) const; // retrieve global coeffs into coef[]
  REAL surfaceError(Vec pt) const;

  // v is the point the line passing through, dirc is the unit vector parallel
  // to the line
  bool intersectWithLine(Vec v, Vec dirc, Vec rt[]) const;

  // find the point on plane which is deepest into a particles, px + qy + rz + s
  // = 0 is the equation
  // of the plane, true means intersection; false means no intersection.
  bool nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, Vec& ptnp) const;

  // calculate the normal force between particle and a plane rigid boundary
  void planeRBForce(PlaneBoundary* plane,
                    BoundaryTangentArrayMap& BoundarytgtMap,
                    BoundaryTangentArray& vtmp);

  // calculate the normal force between particle and a cylinder wall
  Vec cylinderRBForce(std::size_t boundaryId, const Cylinder& S, int side);
  void clearFluidGrid();
  void recordFluidGrid(std::size_t i, std::size_t j, std::size_t k,
                       REAL volFrac);
  std::vector<std::vector<REAL>>& getFluidGrid() { return d_fluidGrid; }

  // Check if the DEM particle contains a point + buffer
  // and returns a local coordinate of the point
  bool containsPoint(const dem::Vec& point,
                     const dem::Vec& dem_point,
                     const REAL& bufferLength,
                     dem::Vec& localCoord,
                     bool& insideGhostLayer) const;

private:
  void init();

private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_id;
    ar& d_type;
    ar& d_a;
    ar& d_b;
    ar& d_c;
    ar& d_young;
    ar& d_poisson;
    ar& d_currPos;
    ar& d_prevPos;
    ar& d_currDirecA;
    ar& d_currDirecB;
    ar& d_currDirecC;
    ar& d_prevDirecA;
    ar& d_prevDirecB;
    ar& d_prevDirecC;
    ar& d_currVeloc;
    ar& d_prevVeloc;
    ar& d_currOmga;
    ar& d_prevOmga;
    ar& d_force;
    ar& d_forceIDMap;
    ar& d_prevForce;
    ar& d_moment;
    ar& d_momentIDMap;
    ar& d_prevMoment;
    ar& d_constForce;
    ar& d_constMoment;
    ar& d_density;
    ar& d_mass;
    ar& d_volume;
    ar& d_momentJ;
    ar& d_coef;
    ar& d_kinetEnergy;
    ar& d_contactNum;
    ar& d_inContact;
    ar& d_fluidGrid;
  }
};

} // namespace dem ends

#endif
