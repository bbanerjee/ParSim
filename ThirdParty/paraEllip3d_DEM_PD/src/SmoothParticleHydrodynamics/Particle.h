#ifndef PARTICLE_H
#define PARTICLE_H

#include "Parameter.h"
#include "realtypes.h"
#include "Vec.h"
#include "Gradation.h"
#include "Rectangle.h"
#include "Cylinder.h"
#include "Boundary.h"
#include "SPHParticle.h"
#include <cstddef>
#include <map>
#include <vector>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace dem {
  
  class Particle {

  private:
    // types of individual particle:
    //   0 - free particle
    //   1 - fixed particle
    //   2 - special case 2 (pure moment): translate first, then rotate only, MNT_START needs to be defined
    //   3 - special case 3 (displacemental ellipsoidal pile): translate in vertical direction only
    //   4 - special case 4 (impacting ellipsoidal penetrator): impact with inital velocity in vertical direction only
    //   5 - free boundary particle
    //   6 - translate only, no rotation
    //  10 - ghost particle
    std::size_t  id;
    std::size_t  type;            
    REAL a, b, c;    // three semi-axle length, must satisfy a >= b >= c
    REAL young;      // note: a(currDirecA), b(currDirecB), c(currDirecC) corresponds to x, y, z in local frame, respectively
    REAL poisson;
    Vec  currPos;    // particle center
    Vec  prevPos;
    Vec  currDirecA, currDirecB, currDirecC; // direction of the three axles, in radian
    Vec  prevDirecA, prevDirecB, prevDirecC;
    Vec  currVeloc;  // the velocity of the mass center
    Vec  prevVeloc;
    Vec  currOmga;   // angular velocity in global frame!
    Vec  prevOmga;
    Vec  force;
    Vec  prevForce;
    Vec  moment;
    Vec  prevMoment;
    Vec  constForce;
    Vec  constMoment;
    REAL density;    // specific gravity
    REAL mass;
    REAL volume;
    Vec  momentJ;    // moment of inertia in local body-fixed frame
    REAL coef[10];   // particle's coefficients in global coordinates
    REAL kinetEnergy;// kinetic energy
    std::size_t  contactNum;
    bool inContact;  // in contact with other particle or boundary
    std::vector< std::vector<REAL> > fluidGrid;

  public:
    std::vector<sph::SPHParticle*> SPHGhostParticleVec;

  public:
    Particle();
    Particle(std::size_t n, std::size_t type, Vec center, REAL r, REAL young, REAL poisson);
    Particle(std::size_t n, std::size_t type, Vec center, REAL a, REAL b, REAL c, REAL young, REAL poisson);
    Particle(std::size_t n, std::size_t type, Vec center, Gradation& grad, REAL young, REAL poisson);
    Particle(std::size_t n, std::size_t type, Vec dim, Vec position, Vec dirca, Vec dircb, Vec dircc, REAL young, REAL poisson);
    
    std::size_t  getId() const {return id;}
    std::size_t  getType() const {return type;}
    REAL getA() const {return a;}
    REAL getB() const {return b;}
    REAL getC() const {return c;}
    REAL getYoung() const {return young;}
    REAL getPoisson() const {return poisson;};
    REAL getVolume() const {return volume;}
    REAL getMass() const {return mass;}
    REAL getDensity() const {return density;}
    Vec  getCurrPos() const {return currPos;}
    Vec  getPrevPos() const {return prevPos;}
    Vec  getCurrDirecA() const {return currDirecA;}
    Vec  getCurrDirecB() const {return currDirecB;}
    Vec  getCurrDirecC() const {return currDirecC;}
    Vec  getPrevDirecA() const {return prevDirecA;}
    Vec  getPrevDirecB() const {return prevDirecB;}
    Vec  getPrevDirecC() const {return prevDirecC;}
    Vec  getCurrVeloc() const {return currVeloc;}
    Vec  getPrevVeloc() const {return prevVeloc;}
    Vec  getCurrOmga() const {return currOmga;}
    Vec  getPrevOmga() const {return prevOmga;}
    Vec  getForce() const {return force;}
    Vec  getMoment() const {return moment;}
    Vec  getAccel() const {return force/mass;}
    Vec  getConstForce() const {return constForce;}
    Vec  getConstMoment() const {return constMoment;}
    Vec  getmomentJ() const {return momentJ;}
    bool isInContact() const {return inContact;}
    std::size_t  getContactNum() const {return contactNum;}

    REAL getRadius(Vec v) const;
    REAL getTransEnergy() const;
    REAL getRotatEnergy() const;
    REAL getKinetEnergy() const;
    REAL getPotenEnergy(REAL ref) const;
    
    void setId(std::size_t n) {id = n;}
    void setType(std::size_t n) {type = n;}
    void setA(REAL dd) {a = dd;}
    void setB(REAL dd) {b = dd;}
    void setC(REAL dd) {c = dd;}
    void expand(REAL percent) {a *=  (1+percent); b *=  (1+percent); c *=  (1+percent);}
    void setCurrPos(Vec vv) {currPos = vv;}
    void setPrevPos(Vec vv) {prevPos = vv;}
    void setCurrDirecA(Vec vv) {currDirecA = vv;}
    void setCurrDirecB(Vec vv) {currDirecB = vv;}
    void setCurrDirecC(Vec vv) {currDirecC = vv;}
    void setPrevDirecA(Vec vv) {prevDirecA = vv;}
    void setPrevDirecB(Vec vv) {prevDirecB = vv;}
    void setPrevDirecC(Vec vv) {prevDirecC = vv;}
    void setCurrVeloc(Vec vv) {currVeloc = vv;}
    void setPrevVeloc(Vec vv) {prevVeloc = vv;}
    void setCurrOmga(Vec vv) {currOmga = vv;}
    void setPrevOmga(Vec vv) {prevOmga = vv;}
    void setForce(Vec vv) {force = vv;}
    void setMoment(Vec vv) {moment = vv;}
    void setConstForce(Vec vv) {constForce = vv;}
    void setConstMoment(Vec vv) {constMoment = vv;}
    void setmomentJ(Vec v) {momentJ = v;}
    void setMass(REAL d) {mass = d;}
    void setDensity(REAL dn) {density = dn;}
    void setInContact(bool value) {inContact = value;}
    void setContactNum(std::size_t num) {contactNum = num;}

    void clearContactForce();
    void addForce(Vec vv) {force += vv;}
    void addMoment(Vec vv) {moment += vv;}
    void update();
    void dragForce();

    Vec globalToLocal(Vec input) const;
    Vec localToGlobal(Vec input) const;
    
    // update global coefficients in the following form based on position/dimensions/orientations
    // a0 x^2 + a1 y^2 + a2 z^2 + a3 xy + a4 yz + a5 zx + a6 x + a7 y + a8 z + a9 = 0
    void globalCoef();  
    void getGlobalCoef(REAL coef[]) const; // retrieve global coeffs into coef[]
    REAL surfaceError(Vec pt) const;
    
    // v is the point the line passing through, dirc is the unit vector parallel to the line
    bool intersectWithLine(Vec v, Vec dirc, Vec rt[]) const;
    
    // find the point on plane which is deepest into a particles, px + qy + rz + s = 0 is the equation 
    // of the plane, true means intersection; false means no intersection.
    bool nearestPTOnPlane(REAL p, REAL q, REAL r, REAL s, Vec &ptnp) const;
    
    // calculate the normal force between particle and a plane rigid boundary
    void planeRBForce(planeBoundary *plane,
		      std::map<std::size_t,std::vector<BoundaryTgt> > &BoundarytgtMap,
		      std::vector<BoundaryTgt> &vtmp);
    
    // calculate the normal force between particle and a cylinder wall
    Vec cylinderRBForce(std::size_t boundaryId, const Cylinder &S, int side);
    void clearFluidGrid();
    void recordFluidGrid(std::size_t i, std::size_t j, std::size_t k);
    std::vector< std::vector<REAL> > & getFluidGrid() { return fluidGrid; }

    // sph 
    void setDemParticleInSPHParticle();
    void setNULLDemParticleInSPHParticle();
    void updateSPHGhostParticle();
    
  private:
    void init();    

  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & id;
      ar & type;
      ar & a & b & c;
      ar & young;
      ar & poisson;
      ar & currPos;
      ar & prevPos;
      ar & currDirecA & currDirecB & currDirecC;
      ar & prevDirecA & prevDirecB & prevDirecC;
      ar & currVeloc;
      ar & prevVeloc;
      ar & currOmga;
      ar & prevOmga;
      ar & force;
      ar & prevForce;
      ar & moment;
      ar & prevMoment;
      ar & constForce;
      ar & constMoment;
      ar & density;
      ar & mass;
      ar & volume;
      ar & momentJ;
      ar & coef;
      ar & kinetEnergy;
      ar & contactNum;
      ar & inContact;
      ar & fluidGrid;
      ar & SPHGhostParticleVec;
    }  
    
  };
  
} // namespace dem ends

#endif
