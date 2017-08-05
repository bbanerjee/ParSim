// Function Definitions
#include "SPHParticle.h"
#include "Particle.h"

namespace sph {

    SPHParticle::SPHParticle()
	: mass(0), density(0), volume(0), pressure(0), mu(0), densityDot(0),
	  initial_X(0), curr_x(0), velocity(0), velocityDot(0), local_X(0), type(0) {}

    SPHParticle::SPHParticle(REAL m, REAL rho, REAL x, REAL y, REAL z, int t){
	mass = m;
    	density = rho;
    	volume = mass/density;
//    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);
   	pressure = dem::Parameter::getSingleton().parameter["P0"]*(std::pow(density/dem::Parameter::getSingleton().parameter["SPHInitialDensity"], dem::Parameter::getSingleton().parameter["gamma"])-1);;
    	mu = density*dem::Parameter::getSingleton().parameter["nu"];
    	densityDot = 0;

    	initial_X = dem::Vec(x,y,z);
	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

	local_X = 0;	// free SPH point
	
	if(t!=1 && t!=3){
	    std::cout << "Error in creating SPH free/boundary particle!" << std::endl;
	    exit(-1);
	}
	type = t;
	demParticle = NULL;
	
    } // end SPHParticle()

    SPHParticle::SPHParticle(REAL m, REAL rho, REAL x, REAL y, REAL z, dem::Vec local, int t){
	mass = m;
    	density = rho;
    	volume = mass/density;
//    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);
    	pressure = dem::Parameter::getSingleton().parameter["P0"]*(std::pow(density/dem::Parameter::getSingleton().parameter["SPHInitialDensity"], dem::Parameter::getSingleton().parameter["gamma"])-1);
    	mu = density*dem::Parameter::getSingleton().parameter["nu"];
    	densityDot = 0;

    	initial_X = dem::Vec(x,y,z);
	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

	local_X = local;	// ghost SPH point

	if(t!=2){
	    std::cout << "Error in creating ghost SPH particle!" << std::endl;
	    exit(-1);
	}
	type = 2;
	demParticle = NULL;
	
    } // end SPHParticle()
	
    void SPHParticle::initial(){

    	volume = mass/density;
//    	pressure = dem::P0*(std::pow(density/dem::SPHInitialDensity, dem::gamma)-1);
    	pressure = dem::Parameter::getSingleton().parameter["P0"]*(std::pow(density/dem::Parameter::getSingleton().parameter["SPHInitialDensity"], dem::Parameter::getSingleton().parameter["gamma"])-1);
    	mu = density*dem::Parameter::getSingleton().parameter["nu"];
    	densityDot = 0;

	curr_x = initial_X;
    	velocity = 0;
    	velocityDot = 0;

    } // end initial()

    void SPHParticle::updateParticle(){	// here the forward Euler time integration is used

//if(fabs(densityDot)>100){
//    densityDot = 0.1*densityDot;
//}
	density = density+densityDot*(dem::Parameter::getSingleton().parameter["timeStep"]);
	velocity = velocity+velocityDot*(dem::Parameter::getSingleton().parameter["timeStep"])-(dem::Parameter::getSingleton().parameter["sphDamping"])*velocity;
	curr_x = curr_x+(velocity+velocityCorrection)*(dem::Parameter::getSingleton().parameter["timeStep"]);
//	curr_x = curr_x+velocity*(dem::Parameter::getSingleton().parameter["timeStep"]);

    } // end updateParticle()

    void SPHParticle::updateParticleDensity(){	// here the forward Euler time integration is used

	density = density+densityDot*(dem::Parameter::getSingleton().parameter["timeStep"]);

    } // end updateParticle()


    void SPHParticle::updateParticlePositionDensityLeapFrog(){	// update position and density based on equation (4.1)

	curr_x = curr_x+(velocity+velocityCorrection)*(dem::Parameter::getSingleton().parameter["timeStep"]);
	density = density+densityDot*(dem::Parameter::getSingleton().parameter["timeStep"]);

    } // end updateParticlePositionDensityLeapFrog()

    void SPHParticle::updateParticleVelocityLeapFrog(){	// update velocity based on equation (4.2)

	velocity = velocity+velocityDot*(dem::Parameter::getSingleton().parameter["timeStep"]);

    } // end updateParticleVelocityLeapFrog()

    void SPHParticle::initialParticleVelocityLeapFrog(){	// update velocity based on equation (4.3)

	velocity = velocity+velocityDot*(dem::Parameter::getSingleton().parameter["timeStep"])*0.5;

    } // end initialParticleVelocityLeapFrog()



} // end namespace sph

