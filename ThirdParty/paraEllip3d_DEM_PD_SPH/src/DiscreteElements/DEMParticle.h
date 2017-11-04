#ifndef PARTICLE_H
#define PARTICLE_H

#include <Boundary/Boundary.h>
#include <Boundary/PlaneBoundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/Cylinder.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
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
public:
  enum class DEMParticleShape
  {
    NONE = 0,
    ELLIPSOID = 1,
    SPHERE = 2,
    CUBE = 3,
    POLYELLIPSOID = 4
  };

  static std::string getDEMParticleShape(DEMParticleShape shape);
  static DEMParticleShape getDEMParticleShape(const std::string& shape);

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
  enum class DEMParticleType
  {
    FREE = 0,
    FIXED = 1,
    ROTATE_ONLY = 2,
    TRANSLATE_Z_ONLY = 3,
    IMPACT_Z_ONLY = 4,
    BOUNDARY_FREE = 5,
    TRANSLATE_ONLY = 6,
    GHOST = 10
  };

  static std::string getDEMParticleType(DEMParticleType type);
  static DEMParticleType getDEMParticleType(const std::string& type);

  DEMParticle();
  DEMParticle(std::size_t n, DEMParticleShape shape, DEMParticleType type, 
              Vec center, REAL r, REAL young, REAL poisson);
  DEMParticle(std::size_t n, DEMParticleShape shape, DEMParticleType type, 
              Vec center, REAL a, REAL b, REAL c, REAL young, REAL poisson);
  DEMParticle(std::size_t n, DEMParticleShape shape, DEMParticleType type, 
              Vec center, Gradation& grad, REAL young, REAL poisson);
  DEMParticle(std::size_t n, DEMParticleShape shape, DEMParticleType type, 
              Vec dim, Vec position, Vec dirca, Vec dircb, Vec dircc, 
              REAL young, REAL poisson);

  std::size_t getId() const { return d_id; }
  const DEMParticleShape& getShape() const { return d_shape; }
  const DEMParticleType&  getType() const { return d_type; }
  REAL radiusA() const { return d_a; }
  REAL radiusB() const { return d_b; }
  REAL radiusC() const { return d_c; }
  REAL youngsModulus() const { return d_young; }
  REAL poissonsRatio() const { return d_poisson; };
  REAL volume() const { return d_volume; }
  REAL mass() const { return d_mass; }
  REAL getTotalMass() const { return d_totalMassAllParticles; }
  REAL density() const { return d_density; }
  Vec currentPosition() const { return d_currPos; }
  Vec previousPosition() const { return d_prevPos; }
  Vec currentAnglesAxisA() const { return d_currDirecA; }
  Vec currentAnglesAxisB() const { return d_currDirecB; }
  Vec currentAnglesAxisC() const { return d_currDirecC; }
  Vec previousAnglesAxisA() const { return d_prevDirecA; }
  Vec previousAnglesAxisB() const { return d_prevDirecB; }
  Vec previousAnglesAxisC() const { return d_prevDirecC; }
  Vec currentVelocity() const { return d_currentVelocity; }
  Vec previousVelocity() const { return d_previousVelocity; }
  Vec currentAngularVelocity() const { return d_currOmga; }
  Vec previousAngularVelocity() const { return d_prevOmga; }
  Vec force() const { return d_force; }
  std::map<size_t, Vec> forceIDMap() const { return d_forceIDMap; }
  Vec moment() const { return d_moment; }
  std::map<size_t, Vec> momentIDMap() const { return d_momentIDMap; }
  Vec acceleration() const { return d_force / d_mass; }
  Vec getConstForce() const { return d_constForce; }
  Vec getConstMoment() const { return d_constMoment; }
  Vec getmomentJ() const { return d_momentJ; }
  bool isInContact() const { return d_inContact; }
  std::size_t getNumBoundaryContacts() const { return d_contactNum; }

  REAL computeRadius(Vec v) const;
  REAL computeTranslationalEnergy() const;
  REAL computeRotationalEnergy() const;
  REAL computeKineticEnergy() const;
  REAL computePotentialEnergy(REAL ref) const;

  void setId(std::size_t n) { d_id = n; }
  void setShape(std::size_t n) { d_shape = static_cast<DEMParticleShape>(n); }
  void setType(std::size_t n) { d_type = static_cast<DEMParticleType>(n); }
  void setRadiusA(REAL dd) { d_a = dd; }
  void setRadiusB(REAL dd) { d_b = dd; }
  void setRadiusC(REAL dd) { d_c = dd; }
  void expand(REAL percent)
  {
    d_a *= (1 + percent);
    d_b *= (1 + percent);
    d_c *= (1 + percent);
  }
  void setCurrentPosition(Vec vv) { d_currPos = vv; }
  void setPreviousPosition(Vec vv) { d_prevPos = vv; }
  void setCurrentAnglesAxisA(Vec vv) { d_currDirecA = vv; }
  void setCurrentAnglesAxisB(Vec vv) { d_currDirecB = vv; }
  void setCurrentAnglesAxisC(Vec vv) { d_currDirecC = vv; }
  void setPreviousAnglesAxisA(Vec vv) { d_prevDirecA = vv; }
  void setPreviousAnglesAxisB(Vec vv) { d_prevDirecB = vv; }
  void setPreviousAnglesAxisC(Vec vv) { d_prevDirecC = vv; }
  void setCurrentVelocity(Vec vv) { d_currentVelocity = vv; }
  void setPreviousVelocity(Vec vv) { d_previousVelocity = vv; }
  void setCurrentAngularVelocity(Vec vv) { d_currOmga = vv; }
  void setPreviousAngularVelocity(Vec vv) { d_prevOmga = vv; }
  void setForce(Vec vv)
  {
    d_force = vv;
  }
  void setMoment(Vec vv) { d_moment = vv; }
  void setConstForce(Vec vv) { d_constForce = vv; }
  void setConstMoment(Vec vv) { d_constMoment = vv; }
  void setmomentJ(Vec v) { d_momentJ = v; }
  void setMass(REAL d) { d_mass = d; }
  void setTotalMass(REAL mass) { d_totalMassAllParticles = mass; }
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
  void computeGlobalCoef();
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
                    BoundaryTangentArrayMap& BoundarytangentMap,
                    BoundaryTangentArray& vtmp);

  // calculate the normal force between particle and a cylinder wall
  Vec cylinderRBForce(std::size_t BoundaryID, const Cylinder& S, int side);
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

  REAL shortestDistToBoundary(const Vec& point) const;

  void dragForce();

  friend std::ostream& operator<<(std::ostream& os, const DEMParticle& pp)
  {
    os << "ID = " << pp.d_id << " Pos: " << pp.d_currPos 
       << " Rad: (" << pp.d_a << ", " << pp.d_b << ", " << pp.d_c << ")"
       << " Ax_a: " << pp.d_currDirecA
       << " Ax_b: " << pp.d_currDirecB
       << " Ax_c: " << pp.d_currDirecC << "\n";
    return os;
  }

private:
  
  std::size_t d_id;
  DEMParticleShape d_shape;
  DEMParticleType d_type;

  // three semi-axle length, must satisfy a >= b >= c
  // note: a(currDirecA), b(currDirecB), c(currDirecC) corresponds
  // to x, y, z in local frame, respectively
  REAL d_a, d_b, d_c; 
  REAL d_young; 
  REAL d_poisson;

  // particle center
  Vec d_currPos; 
  Vec d_prevPos;

  // direction angles of the three axles, in radian
  Vec d_currDirecA, d_currDirecB, d_currDirecC; 
  Vec d_prevDirecA, d_prevDirecB, d_prevDirecC;

  // the velocity of the mass center
  Vec d_currentVelocity; 
  Vec d_previousVelocity;

  // angular velocity in global frame!
  Vec d_currOmga; 
  Vec d_prevOmga;
  Vec d_force;
  std::map<size_t, Vec> d_forceIDMap;

  Vec d_prevForce;
  Vec d_moment;
  std::map<size_t, Vec> d_momentIDMap;

  Vec d_prevMoment;
  Vec d_constForce;
  Vec d_constMoment;

  // specific gravity
  REAL d_density; 
  REAL d_mass;
  REAL d_volume;

  // total mass of all particles (needs to be kept at ecah particle
  // for traction BCs where a = f/m is solved)
  REAL d_totalMassAllParticles; 

  // moment of inertia in local body-fixed frame
  Vec d_momentJ;      

  // particle's coefficients in global coordinates
  REAL d_coef[10];    

  REAL d_kineticEnergy;
  std::size_t d_contactNum;

  // in contact with other particle or boundary
  bool d_inContact; 
  std::vector<std::vector<REAL>> d_fluidGrid;

  void init();

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_id;
    ar& d_shape;
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
    ar& d_currentVelocity;
    ar& d_previousVelocity;
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
    ar& d_totalMassAllParticles;
    ar& d_momentJ;
    ar& d_coef;
    ar& d_kineticEnergy;
    ar& d_contactNum;
    ar& d_inContact;
    ar& d_fluidGrid;
  }
};

} // namespace dem ends

#endif
