#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H


#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
#include <Core/Types/integertypes.h>
#include <InputOutput/InputParameter.h>

#include <iostream>
#include <vector>
#include <math.h>
#include <boost/mpi.hpp>

namespace dem {
  class DEMParticle;
}

namespace sph {

enum class SPHParticleType : char {
  FREE,
  GHOST,
  BOUNDARY,
  NONE
};

class SPHParticle {

public:

  // Default Constructor
  SPHParticle();
  SPHParticle(ParticleID id, REAL mass, REAL density, REAL x, REAL y, REAL z, 
              const dem::Vec& local, SPHParticleType type);

  ~SPHParticle() {};

  // initialize volume, pressure, mu, displacement, velocity and acceleration
  void initialize();	
  REAL calculateVolume() { return d_mass/d_density; }
  REAL calculatePressure();
  REAL calculateViscosity();

  void setDensityDot(REAL a) { d_densityDot = a; }
  void setVelocityDot(const dem::Vec& a) { d_acceleration = a; }
  void setDensityDotVelocityDotZero() {
    d_densityDot = 0; d_acceleration = 0; d_velocityCorrection = 0;
  }
  void setCurrPosition(const dem::Vec& a) { d_currPos = a; }
  void setCurrPositionX(REAL a) { d_currPos.setX(a); }
  void setCurrPositionY(REAL a) { d_currPos.setY(a); }
  void setCurrPositionInitial() { d_currPos = d_initialPos; }
  void setCurrVelocity(const dem::Vec& a) { d_velocity = a; }
  void setType(SPHParticleType a) { d_type = a; }
  void setDemParticle(dem::DEMParticle* p) { d_demParticle = p; }
  void setNULLDemParticle() { d_demParticle = nullptr; }
  void addVelocityDot(const dem::Vec& a) { d_acceleration += a; }
  void addVelocityCorrection(const dem::Vec& a) { d_velocityCorrection += a; }
  void addDensityDot(REAL a) { d_densityDot += a; }
  void addCurrPositionX(REAL a) {d_currPos.setX(d_currPos.x()+a);}

  ParticleID getId() const { return d_id; }
  REAL getMass() const {return d_mass;}
  REAL getDensity() const {return d_density;}
  REAL getDensityDot() const {return d_densityDot;}

  // everytime when using getVolume(), 
  // make sure that volume, pressure, viscosity  have been calculated!!!!
  REAL getVolume() const {return d_volume;}	
  REAL getPressure() const {return d_pressure;}	
  REAL getViscosity() const {return d_mu;}	

  REAL getKineticEnergy() const {return 0.5*d_mass*dot(d_velocity,d_velocity);}
  dem::Vec getInitPosition() const {return d_initialPos;}
  dem::Vec currentPosition() const {return d_currPos;}
  dem::Vec getLocalPosition() const{return d_localCoords;}
  dem::Vec getTrialPosition() const;
  dem::Vec getDisplacement() const {return d_currPos-d_initialPos;}
  dem::Vec getVelocity() const {return d_velocity;}
  dem::Vec getVelocityDot() const {return d_acceleration;}
  SPHParticleType getType() const {return d_type;}
  dem::DEMParticle* getDEMParticle() {return d_demParticle;}

  void fixXYZ() {
    d_acceleration = 0;
    d_velocityCorrection = 0;
  }
  void fixYandZ() {
    fixY();
    fixZ();
  }
  void fixZ() {
    d_acceleration.setZ(0);
    d_velocityCorrection.setZ(0);
  }
  void fixY() {
    d_acceleration.setY(0);
    d_velocityCorrection.setY(0);
  }
  void fixZinApplyBoundary(REAL fix_z) {
    fixZ();
    d_velocity.setZ(0); 
    d_currPos.setZ(fix_z);
  }

  void update();	// update density, velocity and positions
  void updateDensity(const REAL& delT);	// only update density
  void updatePosition(const REAL& delT);	// only update position
  void updateVelocity(const REAL& delT);

  void updatePositionDensityLeapFrog(const REAL& delT);
  void initialVelocityLeapFrog(const REAL& delT);


  private:

  ParticleID d_id;
  REAL d_mass;		
  REAL d_density;	
  REAL d_volume;	
  REAL d_pressure;

  // dynamic viscosity, mu=density*dem::nu; 
  // need to be calculated before being used!!!
  REAL d_mu;	
  REAL d_densityDot;

  dem::Vec d_initialPos;
  dem::Vec d_currPos;	

  dem::Vec d_velocity;
  dem::Vec d_acceleration;

  // velocity correction, the delta_a term, by Monanghan's paper(1994)
  dem::Vec d_velocityCorrection;	

  // for the ghost point only, the local coordinates in the dem particle
  dem::Vec d_localCoords;	

  // variable for linked list searching
  // particle type: 1, free particle; 2, ghost particle; 3, boundary particle
  // if is ghost particle, then pointer to its dem particle, we don't have 
  // this in paraDEM-SPH 
  // there is two tricks in the implementation of MPI: 
  // (1) before send sph particles, demParticle = NULL; 
  //     after receive, assign demParticle. 
  // (2) When delete dem particles, should also free the memory of 
  //     SPHGhostParticleVec
  SPHParticleType d_type;	
  dem::DEMParticle* d_demParticle;	

  friend class boost::serialization::access;
  template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & d_id;
      ar & d_mass;
      ar & d_density;
      ar & d_volume;
      ar & d_pressure;
      ar & d_mu;
      ar & d_densityDot;
      ar & d_currPos;
      ar & d_initialPos;
      ar & d_velocity;
      ar & d_acceleration;
      ar & d_velocityCorrection;
      ar & d_localCoords;
      ar & d_type;
      ar & d_demParticle;
  }
   
}; // end particle


} // end namespace sph

#endif
