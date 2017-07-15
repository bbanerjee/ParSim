#ifndef PERIDYNAMICS_H
#define PERIDYNAMICS_H

#include <Peridynamics/PeriBoundaryBond.h>
#include <Peridynamics/PeriContainers.h>
#include <Peridynamics/PeriDEMBond.h>
#include <Peridynamics/PeriElement.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/globfuncs.h>

#define isProc0_macro ( dem::Peridynamics::getMPIRank() == 0 )
#define proc0cout if( isProc0_macro ) std::cout
#define proc0cerr if( isProc0_macro ) std::cerr

namespace pd {

class Peridynamics
{
public:

  Peridynamics();
  ~Peridynamics();

  // MPI methods
  int getMPIRank() const { return mpiRank; }
  boost::mpi::communicator getMPIWorld() const { return boostWorld; }

  // Accessor methods
  const PeriParticlePArray& getPeriParticleVec() const
  {
    return periParticleVec;
  }
  const PeriParticlePArray& getRecvPeriParticleVec() const
  {
    return recvPeriParticleVec;
  }

  // Pure peridynamics only
  // initializes the velocity, displacement and acceleration
  // @param const std::string& - name of the input file
  void initializePeridynamics(const std::string&); 

  // readData - reads controlling parameters, particle
  // positions and mesh connectivities
  // @param const std::string& - reference of the input file name
  void readPeriDynamicsData(const std::string&); 

  // calcParticleVolume - calculates the particle
  // volume associated with each particle
  void calcParticleVolume();   
  void calcParticleVolumeForHex();   
  void calcParticleVolumeForTet();   

  // calcHorizonSize - calculates the horizon size for
  // each peri-particle
  void calcHorizonSize(); 
  void calcHorizonSizeForHex(); 
  void calcHorizonSizeForTet(); 

  // neibhor - searches and constructs the particle
  // neighborlists
  void constructNeighbor(); 

  void findPeriParticleInBox(const Box& container, 
                             const PeriParticlePArray& allParticles,
                             PeriParticlePArray& foundParticles);

  void removePeriParticleOutBox(const Box& container,
                                PeriParticlePArray& periParticleVec);

  void commuPeriParticle();
  void releaseRecvPeriParticle();
  void releasePeriBondVec();
  void releaseGatheredPeriParticle();
  void migratePeriParticle();
  void gatherPeriParticle();
  void findFixedPeriParticles();    // for all cpus
  void findBoundaryPeriParticles(); // for all cpus
  void applyPeriBoundaryCondition();

  void updatePeriGrid();

  void openPeriProgress(std::ofstream& ofs, const std::string& str);
  void printPeriProgress(std::ofstream& ofs, const int iframe) const;
  void printPeriProgressHalf(std::ofstream& ofs, const int iframe) const;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  int getnPeriParticle() const
  {
    return nPeriParticle;
  } // getnPeriParticle - returns number of peri-particles

  void setInitIsv();
  void prescribeEssentialBoundaryCondition(
    const int); // apply the fixed boundary condition on the bottom particles
  void solve(const std::string&); // solve - solve this problem
  // @param const std::string& - output file name for tecplot visualization
  void runFirstHalfStep();  // run step1 and step2 in Page 5 of Houfu's notes,
                            // equations (15) and (16)
  void runSecondHalfStep(); // run step3 in page 5 of Houfu's notes, equation
                            // (17)
  void findRecvPeriBonds(); // find peri-bonds between periParticleVec and
                            // recvPeriParticleVec
  void findPeriDEMBonds();  // find sand-peri bonds in each cpu, i.e.
                            // periParticleVec and ParticleVec
  void clearPeriDEMBonds();
  void eraseBrokenPeriDEMBonds();
  void writeMesh(const std::string&); // writeMesh - outputs the mesh, used for
                                      // problem checking
  // @param char * - reference of the output file name
  void writeMeshCheckVolume(const std::string&); // writeMeshCheckVolume -
                                                 // outputs the mesh and volume,
                                                 // will
  // @param (std::ofstream, int) - (reference of the output file name, frame
  // index)
  void printPeriDomain(const std::string&) const;
  void printRecvPeriDomain(const std::string&) const;
  void printPeriParticle(const std::string& str) const;

  void printPeriDomainSphere(
    const std::string&) const; // print stress in spherical coordinates

  void constructPeriMatrix();     // construct Matrix members in periParticleVec
  void constructRecvPeriMatrix(); // construct Matrix members in
                                  // recvPeriParticleVec, construction here
                                  // since currently
  // the pointer array in Matrix cannot be transfered well between cpus
  void calcDeformationGradient(); // calcDeformationGradient - calculates the
  // deformation gradient for each peri-particle
  void calcParticleKinv();     // calculate the inverse of K for all particles
  void calcRecvParticleKinv(); // calculate the inverse of K for all
                               // recvPeriParticles
  void calcParticleStress();   // calculate the Cauchy Stress for all particles
  void calcParticleAcceleration(); // calculate the acceleration for all
                                   // particles
  void checkBondParticleAlive(); // for each particle, check the state(alive or
                                 // not) of each surrounding bond;
  // if there's no alive bond for a particle, then this particle is disabled.
  void ApplyExternalForce(int istep);
  void writeDisplacementData(const std::string&, const std::string&,
                             const std::string&);

  void removeInsidePeriParticles(); // delete those peri-points that are inside
                                    // sand particles
  void constructBoundarySandPeriBonds(); // construct boundary bonds, sand
                                         // bonds, July 14, 2014
  void applyCoupledForces(); // check if peri-boundary bonds and peri-dem bonds
                             // are still alive, July 15, 2014
  // if so, apply coupled forces to peri-points and sand particles/boundaries
  void applyCoupledBoundary(); // this is used to test the coupled force model,
                               // October 10, 2014
  // in this test model, the sand-peri-points will move along the dem-particle
  void rigidInclusion();
  void pullOutPeri();


private:

  // MPI data
  boost::mpi::communicator boostWorld;
  MPI_Comm mpiWorld, cartComm;
  int mpiProcX, mpiProcY, mpiProcZ;
  int mpiRank, mpiSize, mpiTag, mpiCoords[3];
  int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
  int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2;
  int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2;
  int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2;
  int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2;
  int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;

  // stream
  std::ofstream periProgInf;
  std::ofstream periProgInfHalf;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// pd part
  /////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int nPeriParticle; // number of PeriParticles in the domain, only for master
                     // cpu
  int nele;          // number of elements in the mesh, only for master cpu

  int ndim; // dimension of the problem, only for master cpu
            //  int nsteps;	// number of total steps
            //  int printInterval;	// print interval

  PeriElementArray connectivity; // mesh connectivity, only for master cpu

  REAL point_interval; // for all cpus, broadcast in scatterDEMPeriParticle()
  REAL maxHorizonSize; // the maximum horizon size of all peri-points, for all
                       // cpus

  PeriParticlePArray
    allPeriParticleVecInitial; // only for master cpu
  // in simulations, we have peri-points and dem sands, for convience of domain
  // generation, we input a cuboid domain
  // filling with peri-points, this particle vector stores all these
  // peri-points. Then since we have sand particles also
  // in the same assembly, we need to remove those peri-points that are within
  // the sand particles to generate the periParticleVec
  // that are used for the calculation. July 14, 2014
  PeriParticlePArray allPeriParticleVec; // only for master cpu
  PeriParticlePArray
    periParticleVec; // coordinates of all the particles in this cpu, local
                     // peri-points in current cpu
  PeriParticlePArray interfacePeriParticleVec;
  PeriParticlePArray
    fixedPeriParticleVec; // this is the peri-points that are
                          // inside the rigid
                          // inclusion, fix these peri-points
  PeriParticlePArray
    outerfacePeriParticleVec; // this is only for the hollow
                              // spherical example
  //  PeriParticlePArray bottomBoundaryVec;	// particles
  // that are in the bottom boundary
  //  PeriParticlePArray topBoundaryVec;	// particles that
  // are in the bottom boundary
  //  PeriParticlePArray cubicTopBoundaryVec;	// particles
  // that are in the top boundary of the cubic
  //  PeriBondPArray totalBondVec;	// all the Bonds in the
  // domain

  // for all cpus, these bonds, including boundary bonds, peri-bonds, and
  // peri-DEM-bonds will
  // constructed after scattering within each cpu
  PeriBoundaryBondPArray bottomBoundaryBondVec;
  PeriBoundaryBondPArray topBoundaryBondVec;
  PeriBoundaryBondPArray leftBoundaryBondVec;
  PeriBoundaryBondPArray rightBoundaryBondVec;
  PeriBoundaryBondPArray frontBoundaryBondVec;
  PeriBoundaryBondPArray backBoundaryBondVec;

  PeriDEMBondPArray
    periDEMBondVec; // the bonds that between the sand particle and
  // the peri-points that are near to this particle, in each cpu

  PeriParticlePArray rperiParticleX1,
    rperiParticleX2; // r stands for received
  PeriParticlePArray rperiParticleY1, rperiParticleY2;
  PeriParticlePArray rperiParticleZ1, rperiParticleZ2;
  PeriParticlePArray rperiParticleX1Y1, rperiParticleX1Y2,
    rperiParticleX1Z1, rperiParticleX1Z2;
  PeriParticlePArray rperiParticleX2Y1, rperiParticleX2Y2,
    rperiParticleX2Z1, rperiParticleX2Z2;
  PeriParticlePArray rperiParticleY1Z1, rperiParticleY1Z2,
    rperiParticleY2Z1, rperiParticleY2Z2;
  PeriParticlePArray rperiParticleX1Y1Z1, rperiParticleX1Y1Z2,
    rperiParticleX1Y2Z1, rperiParticleX1Y2Z2;
  PeriParticlePArray rperiParticleX2Y1Z1, rperiParticleX2Y1Z2,
    rperiParticleX2Y2Z1, rperiParticleX2Y2Z2;
  PeriParticlePArray
    recvPeriParticleVec; // received particles per process
  PeriParticlePArray mergePeriParticleVec;

  PeriParticlePArray
    bottomBoundaryVec; // particles that are in the bottom boundary
  PeriParticlePArray
    frontBoundaryVec; // particles that are in the front y boundary
  PeriParticlePArray
    leftBoundaryVec; // particles that are in the left x boundary

  PeriParticlePArray
    topBoundaryInnerVec; // particles that are in the bottom boundary
  PeriParticlePArray
    topBoundaryEdgeVec; // particles that are in the bottom boundary
  PeriParticlePArray
    topBoundaryCornerVec; // particles that are in the bottom boundary

  PeriBondPArray
    recvPeriBondVec; // peri-bonds between recvPeriParticle in each cpu
  PeriBondPArray
    periBondVec; // peri-bonds between periParticleVec in each cpu

  //  // for elasticity verification purpose.
  //  int Uzindex[39]; // particle index along the z direction.[x = 0, y = 0]
  //  int Uyindex[5];  // particle index along the y direction.[x = 0, z = 0]
  //  int Uxindex[5];  // particle index along the x direction.[y = 0, z = 0]

};

} // namespace pd

#endif
