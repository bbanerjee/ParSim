#ifndef PERIPARTICLE_H
#define PERIPARTICLE_H

#include <algorithm>
#include <iostream>
#include <vector>

#include <Core/Math/Matrix.h>
#include <Core/Math/Vec.h>
#include <Core/Types/RealTypes.h>
#include <Core/Util/Utility.h>
#include <Peridynamics/PeriContainers.h>
#include <InputOutput/InputParameter.h>
#include <Peridynamics/PeriBond.h>
#include <Peridynamics/globfuncs.h>
#include <boost/mpi.hpp>

using ParticleID = std::uint64_t;

namespace pd {

class PeriDomain;

class PeriParticle
{

public:
  // Default Constructor
  PeriParticle();
  PeriParticle(ParticleID id, REAL x, REAL y, REAL z);
  PeriParticle(const PeriParticle&);

  ~PeriParticle();

  // void calcKinv();	// calculate dem::Matrix the inverse of K dem::Matrix
  // void calcAcceleration();	// calculate the acceleration of the particle
  // void calcStress();		// calculate stress

  ParticleID getId() {return d_id;}

  // setVolume - sets the volume of the particle
  // @param newParticleVolume - volume of the particle
  void setVolume(REAL newParticleVolume) { d_volume = newParticleVolume; } 

  // The mass should be set at the beginning of the simulation when the
  // volume is the initial volume
  void setMass(REAL initialVolume) {d_mass = util::getParam<REAL>("periDensity") * initialVolume;}
  void setInitPosition(const dem::Vec& pos) { d_initPosition = pos; }

  REAL mass() const { return d_mass; }
  REAL volume() const { return d_volume; }
  REAL getHorizonSize() const { return d_horizonSize; }
  REAL getSigma11() const { return d_sigma11; }
  REAL getSigma12() const { return d_sigma12; }
  REAL getSigma13() const { return d_sigma13; }
  REAL getSigma21() const { return d_sigma21; }
  REAL getSigma22() const { return d_sigma22; }
  REAL getSigma23() const { return d_sigma23; }
  REAL getSigma31() const { return d_sigma31; }
  REAL getSigma32() const { return d_sigma32; }
  REAL getSigma33() const { return d_sigma33; }
  int getBondsNumber() const { return d_bondVec.size(); }
  dem::Vec getInitPosition() const { return d_initPosition; }
  dem::Vec currentPosition() const { return d_initPosition + d_displacement; }
  dem::Vec previousPosition() const { return d_initPosition + d_prevDisp; }
  dem::Vec getDisplacement() const { return d_displacement; }
  dem::Vec getVelocity() const { return d_velocity; }
  bool getIsAlive() const { return d_isAlive; }
  dem::Matrix getSigma() const { return d_sigma; }
  dem::Matrix getDeformationGradient() const { return d_deformationGradient; }
  dem::Matrix getParticleKinv() const { return d_Kinv; }
  //    dem::Matrix getIsv() const {return isv;}
  REAL getIsv() const { return d_isv11; }
  dem::Vec getVelocityHalf() const { return d_velocityHalf; }

  dem::Vec accelerationeration() const { return d_acceleration; }

  void checkParticleAlive();
  void replaceHorizonSizeIfLarger(
    REAL tmp); // replace this->horizonSize if tmp is larger
  void prescribeBottomDisplacement(REAL disp)
  {
    d_displacement.setZ(disp);
    d_velocity.setZ(0.0);
  }
  void prescribeDisplacementX(REAL disp)
  {
    d_displacement.setX(disp);
    d_velocity.setX(0.0);
  }
  void prescribeDisplacementY(REAL disp)
  {
    d_displacement.setY(disp);
    d_velocity.setY(0.0);
  }
  void prescribeDisplacement(dem::Vec disp)
  {
    d_displacement = disp;
    d_velocity = dem::Vec(0.0);
  }
  void prescribeTopDisplacement(REAL disp)
  {
    d_displacement.setZ(disp);
    d_velocity.setZ(0.0);
  }
  void setInitVelocity(const dem::Vec& tmp) { d_velocity = tmp; }
  void setInitIsv(REAL isv_tmp) { d_isv11 = isv_tmp; }
  void constructMatrixMember(); // construct these Matrix members

  void setAcceleration(dem::Vec newAcceleration)
  {
    d_acceleration = newAcceleration;
  }
  void setCurrentPositionition(dem::Vec curr_posi)
  {
    d_displacement = curr_posi - d_initPosition;
  }
  void addAcceleration(dem::Vec acce_add) { d_acceleration += acce_add; }
  void addAccelerationByForce(dem::Vec force)
  {
    d_acceleration +=
      force / (d_volume *
               util::getParam<REAL>("periDensity") *
               util::getParam<REAL>("massScale"));
  }
  // add acceleration based on force, July 15, 2014
  //    void pushBackNeighborVec(PeriParticle* pt) {neighborVec.push_back(pt);}
  void pushBackBondVec(PeriBondP bt) { d_bondVec.push_back(bt); }
  void setAliveFalse() { d_isAlive = false; }
  void setStressZero() { d_sigma = dem::zeros(3, 3); }
  void setAccelerationZero() { d_acceleration = 0; }
  void setParticleKinv(dem::Matrix& tmp) { d_Kinv = tmp; }
  void calcParticleKinv();
  void calcParticleStress();
  void calcParticleAcceleration();

  void updateDisplacement(); // update displacement and velocityHalf
  void updateVelocity();     // update velocity
  void initial();            // initial displacement, velocity and acceleration
  void eraseRecvPeriBonds();
  void assignSigma();
  void assignKinv();
  void clearPeriBonds() { d_bondVec.clear(); }
  void releaseBondVec();

  void addDEMId(int id) { d_BondedDEMParticleID.push_back(id); }
  void eraseDEMId(int id)
  {
    d_BondedDEMParticleID.erase(
      std::remove(d_BondedDEMParticleID.begin(), d_BondedDEMParticleID.end(), id),
      d_BondedDEMParticleID.end());
  }
  bool isBonded(int id)
  {
    return (std::find(d_BondedDEMParticleID.begin(), d_BondedDEMParticleID.end(),
                      id) != d_BondedDEMParticleID.end());
  }

private:
  bool d_isAlive; // if the peri-particle is alive

  // REAL particleDensity;	// particle density ==> defined globally, same
  // material, usually
  // dem::dem::Vec currPosition;	// current position vector of particle,
  // no
  // need,
  // save spaces
  ParticleID d_id;       // Unique particle ID
  dem::Vec d_initPosition; // initial position vector of particle

  REAL d_mass;             // particle mass
  REAL d_volume;   // particle volume
  dem::Vec d_displacement; // particle displacement
  dem::Vec d_prevDisp;     // displacement at previous step
  dem::Vec d_velocity;     // particle velocity
  dem::Vec d_velocityHalf; // velocity at half time-step, for Velocity-Verlet
                         // integration
  dem::Vec d_acceleration; // particle acceleration

  dem::Matrix d_sigma;               // Cauchy stress
  dem::Matrix d_deformationGradient; // deformation gradient tensor
  dem::Matrix d_deformationGradientHalf;
  dem::Matrix d_Kinv; // the inverse K Matrix
                    //    dem::Matrix isv;		// a 1x5 matrix
  REAL d_isv11; // actually, only isv(1,1) is used, and it may be changed from
              // each time step
  dem::Matrix d_tangentModulus;

  REAL d_horizonSize; // used to get neighbor list

  //    std::vector<PeriParticle*> neighborVec;	// neighbor list of this
  // particle
  PeriBondPArray d_bondVec; // Bonds connected to this particle,
  // peri-bonds will be constructed after scattering in each cpu

  // in order to keep stress values after gathering
  REAL d_sigma11;
  REAL d_sigma12;
  REAL d_sigma13;

  REAL d_sigma21;
  REAL d_sigma22;
  REAL d_sigma23;

  REAL d_sigma31;
  REAL d_sigma32;
  REAL d_sigma33;

  // in order to keep Kinv values after commuPeriParticle
  REAL d_Kinv11;
  REAL d_Kinv12;
  REAL d_Kinv13;

  REAL d_Kinv21;
  REAL d_Kinv22;
  REAL d_Kinv23;

  REAL d_Kinv31;
  REAL d_Kinv32;
  REAL d_Kinv33;

  std::vector<int> d_BondedDEMParticleID; // the ID of DEM particles to which this
                                        // peri-point is bonded

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_id;
    ar& d_isAlive;
    ar& d_initPosition;
    ar& d_mass;
    ar& d_volume;
    ar& d_displacement;
    ar& d_prevDisp;
    ar& d_velocity;
    ar& d_velocityHalf;
    ar& d_acceleration;
    //      ar & sigma;
    //      ar & deformationGradient;
    //      ar & deformationGradientHalf;
    //      ar & Kinv;
    ar& d_isv11;
    //      ar & tangentModulus;
    ar& d_horizonSize;
    //      ar & bondVec;
    ar& d_sigma11;
    ar& d_sigma12;
    ar& d_sigma13;
    ar& d_sigma21;
    ar& d_sigma22;
    ar& d_sigma23;
    ar& d_sigma31;
    ar& d_sigma32;
    ar& d_sigma33;
    ar& d_Kinv11;
    ar& d_Kinv12;
    ar& d_Kinv13;
    ar& d_Kinv21;
    ar& d_Kinv22;
    ar& d_Kinv23;
    ar& d_Kinv31;
    ar& d_Kinv32;
    ar& d_Kinv33;
    ar& d_BondedDEMParticleID;
  }

}; // end particle

} // end pd

#endif
