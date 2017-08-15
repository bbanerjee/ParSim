#ifndef PERIPARTICLE_H
#define PERIPARTICLE_H

#include <algorithm>
#include <iostream>
#include <vector>

#include <Core/Math/Matrix.h>
#include <Core/Math/Vec.h>
#include <Core/Types/realtypes.h>
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
  void setVolume(REAL newParticleVolume) { volume = newParticleVolume; } 

  // The mass should be set at the beginning of the simulation when the
  // volume is the initial volume
  void setMass(REAL initialVolume) {mass = util::getParam<REAL>("periDensity") * initialVolume;}
  void setInitPosition(const dem::Vec& pos) { initPosition = pos; }

  REAL getMass() const { return mass; }
  REAL getVolume() const { return volume; }
  REAL getHorizonSize() const { return horizonSize; }
  REAL getSigma11() const { return sigma11; }
  REAL getSigma12() const { return sigma12; }
  REAL getSigma13() const { return sigma13; }
  REAL getSigma21() const { return sigma21; }
  REAL getSigma22() const { return sigma22; }
  REAL getSigma23() const { return sigma23; }
  REAL getSigma31() const { return sigma31; }
  REAL getSigma32() const { return sigma32; }
  REAL getSigma33() const { return sigma33; }
  int getBondsNumber() const { return bondVec.size(); }
  dem::Vec getInitPosition() const { return initPosition; }
  dem::Vec currentPosition() const { return initPosition + displacement; }
  dem::Vec getPrevPosition() const { return initPosition + prevDisp; }
  dem::Vec getDisplacement() const { return displacement; }
  dem::Vec getVelocity() const { return velocity; }
  bool getIsAlive() const { return isAlive; }
  dem::Matrix getSigma() const { return sigma; }
  dem::Matrix getDeformationGradient() const { return deformationGradient; }
  dem::Matrix getParticleKinv() const { return Kinv; }
  //    dem::Matrix getIsv() const {return isv;}
  REAL getIsv() const { return isv11; }
  dem::Vec getVelocityHalf() const { return velocityHalf; }

  dem::Vec getAcceleration() const { return acceleration; }

  void checkParticleAlive();
  void replaceHorizonSizeIfLarger(
    REAL tmp); // replace this->horizonSize if tmp is larger
  void prescribeBottomDisplacement(REAL disp)
  {
    displacement.setZ(disp);
    velocity.setZ(0.0);
  }
  void prescribeDisplacementX(REAL disp)
  {
    displacement.setX(disp);
    velocity.setX(0.0);
  }
  void prescribeDisplacementY(REAL disp)
  {
    displacement.setY(disp);
    velocity.setY(0.0);
  }
  void prescribeDisplacement(dem::Vec disp)
  {
    displacement = disp;
    velocity = dem::Vec(0.0);
  }
  void prescribeTopDisplacement(REAL disp)
  {
    displacement.setZ(disp);
    velocity.setZ(0.0);
  }
  void setInitVelocity(const dem::Vec& tmp) { velocity = tmp; }
  void setInitIsv(REAL isv_tmp) { isv11 = isv_tmp; }
  void constructMatrixMember(); // construct these Matrix members

  void setAcceleration(dem::Vec newAcceleration)
  {
    acceleration = newAcceleration;
  }
  void setCurrPosition(dem::Vec curr_posi)
  {
    displacement = curr_posi - initPosition;
  }
  void addAcceleration(dem::Vec acce_add) { acceleration += acce_add; }
  void addAccelerationByForce(dem::Vec force)
  {
    acceleration +=
      force / (volume *
               util::getParam<REAL>("periDensity") *
               util::getParam<REAL>("massScale"));
  }
  // add acceleration based on force, July 15, 2014
  //    void pushBackNeighborVec(PeriParticle* pt) {neighborVec.push_back(pt);}
  void pushBackBondVec(PeriBondP bt) { bondVec.push_back(bt); }
  void setAliveFalse() { isAlive = false; }
  void setStressZero() { sigma = dem::zeros(3, 3); }
  void setAccelerationZero() { acceleration = 0; }
  void setParticleKinv(dem::Matrix& tmp) { Kinv = tmp; }
  void calcParticleKinv();
  void calcParticleStress();
  void calcParticleAcceleration();

  void updateDisplacement(); // update displacement and velocityHalf
  void updateVelocity();     // update velocity
  void initial();            // initial displacement, velocity and acceleration
  void eraseRecvPeriBonds();
  void assignSigma();
  void assignKinv();
  void clearPeriBonds() { bondVec.clear(); }
  void releaseBondVec();

  void addDEMId(int id) { BondedDEMParticleID.push_back(id); }
  void eraseDEMId(int id)
  {
    BondedDEMParticleID.erase(
      std::remove(BondedDEMParticleID.begin(), BondedDEMParticleID.end(), id),
      BondedDEMParticleID.end());
  }
  bool isBonded(int id)
  {
    return (std::find(BondedDEMParticleID.begin(), BondedDEMParticleID.end(),
                      id) != BondedDEMParticleID.end());
  }

private:
  bool isAlive; // if the peri-particle is alive

  // REAL particleDensity;	// particle density ==> defined globally, same
  // material, usually
  // dem::dem::Vec currPosition;	// current position vector of particle,
  // no
  // need,
  // save spaces
  ParticleID d_id;       // Unique particle ID
  dem::Vec initPosition; // initial position vector of particle

  REAL mass;             // particle mass
  REAL volume;   // particle volume
  dem::Vec displacement; // particle displacement
  dem::Vec prevDisp;     // displacement at previous step
  dem::Vec velocity;     // particle velocity
  dem::Vec velocityHalf; // velocity at half time-step, for Velocity-Verlet
                         // integration
  dem::Vec acceleration; // particle acceleration

  dem::Matrix sigma;               // Cauchy stress
  dem::Matrix deformationGradient; // deformation gradient tensor
  dem::Matrix deformationGradientHalf;
  dem::Matrix Kinv; // the inverse K Matrix
                    //    dem::Matrix isv;		// a 1x5 matrix
  REAL isv11; // actually, only isv(1,1) is used, and it may be changed from
              // each time step
  dem::Matrix tangentModulus;

  REAL horizonSize; // used to get neighbor list

  //    std::vector<PeriParticle*> neighborVec;	// neighbor list of this
  // particle
  PeriBondPArray bondVec; // Bonds connected to this particle,
  // peri-bonds will be constructed after scattering in each cpu

  // in order to keep stress values after gathering
  REAL sigma11;
  REAL sigma12;
  REAL sigma13;

  REAL sigma21;
  REAL sigma22;
  REAL sigma23;

  REAL sigma31;
  REAL sigma32;
  REAL sigma33;

  // in order to keep Kinv values after commuPeriParticle
  REAL Kinv11;
  REAL Kinv12;
  REAL Kinv13;

  REAL Kinv21;
  REAL Kinv22;
  REAL Kinv23;

  REAL Kinv31;
  REAL Kinv32;
  REAL Kinv33;

  std::vector<int> BondedDEMParticleID; // the ID of DEM particles to which this
                                        // peri-point is bonded

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& d_id;
    ar& isAlive;
    ar& initPosition;
    ar& mass;
    ar& volume;
    ar& displacement;
    ar& prevDisp;
    ar& velocity;
    ar& velocityHalf;
    ar& acceleration;
    //      ar & sigma;
    //      ar & deformationGradient;
    //      ar & deformationGradientHalf;
    //      ar & Kinv;
    ar& isv11;
    //      ar & tangentModulus;
    ar& horizonSize;
    //      ar & bondVec;
    ar& sigma11;
    ar& sigma12;
    ar& sigma13;
    ar& sigma21;
    ar& sigma22;
    ar& sigma23;
    ar& sigma31;
    ar& sigma32;
    ar& sigma33;
    ar& Kinv11;
    ar& Kinv12;
    ar& Kinv13;
    ar& Kinv21;
    ar& Kinv22;
    ar& Kinv23;
    ar& Kinv31;
    ar& Kinv32;
    ar& Kinv33;
    ar& BondedDEMParticleID;
  }

}; // end particle

} // end pd

#endif
