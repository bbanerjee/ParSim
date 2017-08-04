#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H


#include <iostream>
#include <vector>
#include <math.h>
#include <boost/mpi.hpp>

#include "Vec.h"
#include "realtypes.h"
#include "Parameter.h"

namespace dem{
class Particle;
}


namespace sph{

class SPHParticle {

public:
    // Default Constructor
    SPHParticle();
    SPHParticle(REAL mass, REAL density, REAL x, REAL y, REAL z, int t);
    SPHParticle(REAL mass, REAL density, REAL x, REAL y, REAL z, dem::Vec local, int t);

    ~SPHParticle() {};

    void calculateParticleVolume() {volume=mass/density;}
//    void calculateParticlePressure() {pressure= dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);}
    void calculateParticlePressure() {pressure= dem::Parameter::getSingleton().parameter["P0"]*(std::pow(density/dem::Parameter::getSingleton().parameter["SPHInitialDensity"], dem::Parameter::getSingleton().parameter["gamma"])-1);}
    void calculateParticleViscosity() {mu=density*dem::Parameter::getSingleton().parameter["nu"];}	// dynamic viscosity

    void setDensityDot(REAL a) {densityDot = a;}
    void setVelocityDot(const dem::Vec& a) {velocityDot = a;}
    void setDensityDotVelocityDotZero() {densityDot = 0; velocityDot = 0; velocityCorrection = 0;}
    void setCurrPosition(dem::Vec a) {curr_x = a;}
    void setCurrPositionX(REAL a) {curr_x.setX(a);}
    void setCurrPositionY(REAL a) {curr_x.setY(a);}
    void setCurrPositionInitial() {curr_x = initial_X;}
    void setCurrVelocity(dem::Vec a) {velocity = a;}
    void setType(int a) {type=a;}
    void setDemParticle(dem::Particle* p) {demParticle = p;}
    void setNULLDemParticle() {demParticle = NULL;}
    void addVelocityDot(const dem::Vec& a) {velocityDot = velocityDot+a;}
    void addVelocityCorrection(const dem::Vec& a) {velocityCorrection = velocityCorrection+a;}
    void addDensityDot(REAL a) {densityDot = densityDot+a;}
    void addCurrPositionX(REAL a) {curr_x.setX(curr_x.getX()+a);}

    REAL getParticleMass() const {return mass;}
    REAL getParticleDensity() const {return density;}
    REAL getDensityDot() const {return densityDot;}
    REAL getParticleVolume() const {return volume;}	// everytime when use getParticleVolume(), make sure that volume has been calculated!!!!
    REAL getParticlePressure() const {return pressure;}	// everytime when use getParticlePressure(), make sure that pressure has been calculated!!!!
    REAL getParticleViscosity() const {return mu;}	// everytime when use getParticleViscosity(), make sure that viscosity has been calculated!!!!
    REAL getKineticEnergy() const {return 0.5*mass*velocity*velocity;}
    dem::Vec getInitPosition() const {return initial_X;}
    dem::Vec getCurrPosition() const {return curr_x;}
    dem::Vec getLocalPosition() const{return local_X;}
    dem::Vec getTrialPosition() const {return curr_x+velocity*(dem::Parameter::getSingleton().parameter["timeStep"]);}	// return the trial position, used for the boundary handling model based on momentum conservation
    dem::Vec getDisplacement() const {return curr_x-initial_X;}
    dem::Vec getVelocity() const {return velocity;}
    dem::Vec getVelocityDot() const {return velocityDot;}
    int getType() const {return type;}
    dem::Particle* getDemParticle() {return demParticle;}
	
    void fixYandZ() {velocityDot.setY(0); velocityDot.setZ(0); velocityCorrection.setY(0); velocityCorrection.setZ(0);}
    void fixZ() {velocityDot.setZ(0);velocityCorrection.setZ(0);}
    void fixY() {velocityDot.setY(0);velocityCorrection.setY(0);}
    void fixXYZ() {velocityDot=0;velocityCorrection=0;}
    void fixZinApplyBoundary(REAL fix_z) {velocityDot.setZ(0); velocityCorrection.setZ(0); velocity.setZ(0); curr_x.setZ(fix_z);}
    void initial();	// initial displacement, velocity and acceleration
    void updateParticle();	// update density, velocity and positions
    void updateParticleDensity();	// only update density

    void updateParticlePositionDensityLeapFrog();
    void updateParticleVelocityLeapFrog();
    void initialParticleVelocityLeapFrog();

	
private:

    REAL mass;		// mass
    REAL density;	// mass density
    REAL volume;	// volume, v=m/density, need to be calculated before being used!!!
    REAL pressure;	// pressure, need to be calculated before being used!!!
    REAL mu;	// dynamic viscosity, mu=density*dem::nu; need to be calculated before being used!!!
    REAL densityDot;	// density dot

    dem::Vec curr_x;	// current position
    dem::Vec initial_X;	// initial position

    dem::Vec velocity;	// velocity of SPH particle
    dem::Vec velocityDot;	// velocity dot, acceleration
    dem::Vec velocityCorrection;	// velocity correction, the delta_a term, by Monanghan's paper(1994)

    dem::Vec local_X;	// for the ghost point only, the local coordinates in the dem particle

    // variable for linked list searching
    int type;	// particle type: 1, free particle; 2, ghost particle; 3, boundary particle
    dem::Particle* demParticle;	// if is ghost particle, then pointer to its dem particle, we don't have this in paraDEM-SPH
				// there is two tricks in the implementation of MPI: (1) before send sph particles, demParticle = NULL; 
				//   after receive, assign demParticle. (2) When delete dem particles, should also free the memory of SPHGhostParticleVec

    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & mass;
      ar & density;
      ar & volume;
      ar & pressure;
      ar & mu;
      ar & densityDot;
      ar & curr_x;
      ar & initial_X;
      ar & velocity;
      ar & velocityDot;
      ar & velocityCorrection;
      ar & local_X;
      ar & type;
      ar & demParticle;
    }
   
}; // end particle


} // end namespace sph

#endif
