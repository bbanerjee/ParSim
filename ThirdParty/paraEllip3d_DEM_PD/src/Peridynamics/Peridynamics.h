#ifndef PERIDYNAMICS_H
#define PERIDYNAMICS_H

#include <Peridynamics/PeriBoundaryBond.h>
#include <Peridynamics/PeriContainers.h>
#include <Peridynamics/PeriDEMBond.h>
#include <Peridynamics/PeriElement.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/globfuncs.h>
#include <DiscreteElements/Patch.h>

/*
#define isProc0_macro ( dem::Peridynamics::getMPIRank() == 0 )
#define proc0cout if( isProc0_macro ) std::cout
#define proc0cerr if( isProc0_macro ) std::cerr
*/

using PeriPatch = dem::Patch<pd::PeriParticlePArray>;

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

  void setGrid(dem::Box cont) { d_periGrid = cont; }

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

  // calculate the inverse of K for all particles
  void calcParticleKinv();     

  // apply the fixed boundary condition on the bottom particles
  void prescribeEssentialBoundaryCondition(const int); 

  // calculate the Cauchy Stress for all particles
  void calcParticleStress();   

  // calculate the acceleration for all particles
  void calcParticleAcceleration(); 

  // Apply traction bounary conditions
  void applyTractionBoundary(int);

  // Scatter the peridynamics particles
  void scatterPeriParticle(const dem::Box& allContainer);

  void findPeriParticleInBox(const dem::Box& container, 
                             const PeriParticlePArray& allParticles,
                             PeriParticlePArray& foundParticles);

  // construct Matrix members in periParticleVec
  void constructPeriMatrix();

  void commuPeriParticle(int iteration,
                         const double& maxDEMParticleRadius);

  void removePeriParticleOutBox(const dem::Box& container,
                                PeriParticlePArray& periParticleVec);

  // find peri-bonds between periParticleVec and recvPeriParticleVec
  void updatePeridynamicBonds(PeriParticlePArray& recvParticles,
                              PeriBondPArray& recvPeriBondVec,
                              PeriParticlePArray& periParticleVec);
  void findRecvPeriBonds(); 

  // find sand-peri bonds in each cpu, i.e. periParticleVec and ParticleVec
  void findPeriDEMBonds(dem::ParticlePArray mergeParticles);  

  void clearPeriDEMBonds();
  void eraseBrokenPeriDEMBonds();

  bool isBdryProcess();

  void releaseRecvPeriParticle();
  void releasePeriBondVec();
  void releaseGatheredPeriParticle();

  void migratePeriParticle(int iteration);
  void gatherPeriParticle();

  void findFixedPeriParticles();    // for all cpus
  void findBoundaryPeriParticles(); // for all cpus
  void applyPeriBoundaryCondition();

  void updatePeriGrid(const PeriParticlePArray& particles);

  void openPeriProgress(std::ofstream& ofs, const std::string& str);
  void printPeriProgress(std::ofstream& ofs, const int iframe) const;
  void printPeriProgressHalf(std::ofstream& ofs, const int iframe) const;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  int getnPeriParticle() const
  {
    return nPeriParticle;
  } // getnPeriParticle - returns number of peri-particles

  void setInitIsv();

  // solve - solve this problem
  // @param const std::string& - output file name for tecplot visualization
  void solvePurePeridynamics(const std::string&, int printInterval); 

  // run step1 and step2 in Page 5 of Houfu's notes,
  // equations (15) and (16)
  void runFirstHalfStep();  

  // run step3 in page 5 of Houfu's notes, equation (17)
  void runSecondHalfStep(); 

  // writeMesh - outputs the mesh, used for
  // problem checking
  // @param char * - reference of the output file name
  void writeMesh(const std::string&); 

  // writeMeshCheckVolume -
  // outputs the mesh and volume,
  // will
  // @param (std::ofstream, int) - (reference of the output file name, frame
  // index)
  void writeMeshCheckVolume(const std::string&); 

  void writeParticleTecplot(std::ofstream& ofs, const int iframe) const;

  void printPeriDomain(const std::string&) const;
  void printRecvPeriDomain(const std::string&) const;
  void printPeriParticle(const std::string& str) const;

  // print stress in spherical coordinates
  void printPeriDomainSphere(const std::string&) const; 

  // construct Matrix members in
  // recvPeriParticleVec, construction here
  // since currently
  // the pointer array in Matrix cannot be transfered well between cpus
  void constructRecvPeriMatrix(); 

  // calcDeformationGradient - calculates the
  // deformation gradient for each peri-particle
  void calcDeformationGradient(); 

  // calculate the inverse of K for all
  // recvPeriParticles
  void calcRecvParticleKinv(); 

  // for each particle, check the state(alive or
  // not) of each surrounding bond;
  // if there's no alive bond for a particle, then this particle is disabled.
  void checkBondParticleAlive(); 

  void ApplyExternalForce(int istep);
  void writeDisplacementData(const std::string&, const std::string&,
                             const std::string&);

  // delete those peri-points that are inside sand particles
  void removeInsidePeriParticles(const dem::ParticlePArray& particles);     

  // construct boundary bonds, sand bonds, July 14, 2014
  void constructBoundarySandPeriBonds(); 

  // check if peri-boundary bonds and peri-dem bonds
  // are still alive, July 15, 2014
  // if so, apply coupled forces to peri-points and sand particles/boundaries
  void applyCoupledForces(); 

  // this is used to test the coupled force model,
  // October 10, 2014
  // in this test model, the sand-peri-points will move along the dem-particle
  void applyCoupledBoundary(); 

  void rigidInclusion();
  void pullOutPeri();

private:

  // MPI data
  boost::mpi::communicator boostWorld;
  MPI_Comm mpiWorld, cartComm;
  int mpiProcX, mpiProcY, mpiProcZ;
  int mpiRank, mpiSize, mpiTag;

  dem::IntVec d_mpiProcs;
  dem::IntVec d_mpiCoords;

  std::unique_ptr<PeriPatch> d_patchP;
  void createPatch(int iteration, const REAL& ghostWidth);
  void updatePatch(int iteration, const REAL& ghostWidth);

  // the maximum horizon size of all peri-points, for all cpus
  REAL d_maxHorizonSize; 
  REAL d_maxDistBetweenParticles;
  dem::Box d_periGrid;

  int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
  int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2;
  int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2;
  int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2;
  int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2;
  int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;

  // stream
  std::ofstream periProgInf;
  std::ofstream periProgInfHalf;

  // number of PeriParticles in the domain, only for master cpu
  int nPeriParticle; 

  // number of elements in the mesh, only for master cpu
  int nele;          

  // dimension of the problem, only for master cpu
  int ndim; 

  // mesh connectivity, only for master cpu
  PeriElementArray d_connectivity; 

  // only for master cpu
  // in simulations, we have peri-points and dem sands, for convience of domain
  // generation, we input a cuboid domain
  // filling with peri-points, this particle vector stores all these
  // peri-points. Then since we have sand particles also
  // in the same dem, we need to remove those peri-points that are within
  // the sand particles to generate the periParticleVec
  // that are used for the calculation. July 14, 2014
  PeriParticlePArray d_allPeriParticlesInitial; 

  // only for master cpu
  PeriParticlePArray allPeriParticleVec; 

  // coordinates of all the particles in this cpu, local 
  // peri-points in current cpu
  PeriParticlePArray periParticleVec; 
  PeriParticlePArray periParticleVecInitial; 

  PeriParticlePArray interfacePeriParticleVec;

  // this is the peri-points that are
  // inside the rigid
  // inclusion, fix these peri-points
  PeriParticlePArray fixedPeriParticleVec; 

  // this is only for the hollow
  // spherical example
  PeriParticlePArray outerfacePeriParticleVec; 

  // for all cpus, these bonds, including boundary bonds, peri-bonds, and
  // peri-DEM-bonds will
  // constructed after scattering within each cpu
  PeriBoundaryBondPArray bottomBoundaryBondVec;
  PeriBoundaryBondPArray topBoundaryBondVec;
  PeriBoundaryBondPArray leftBoundaryBondVec;
  PeriBoundaryBondPArray rightBoundaryBondVec;
  PeriBoundaryBondPArray frontBoundaryBondVec;
  PeriBoundaryBondPArray backBoundaryBondVec;

  // the bonds that between the sand particle and
  // the peri-points that are near to this particle, in each cpu
  PeriDEMBondPArray periDEMBondVec; 

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

  // particles that are in the top boundary
  PeriParticlePArray topBoundaryVec; 
  // particles that are in the bottom boundary
  PeriParticlePArray bottomBoundaryVec; 
  // particles that are in the front y boundary
  PeriParticlePArray frontBoundaryVec; 
  // particles that are in the left x boundary
  PeriParticlePArray leftBoundaryVec; 

  // particles that are in the bottom boundary
  PeriParticlePArray topBoundaryInnerVec; 
  // particles that are in the bottom boundary
  PeriParticlePArray topBoundaryEdgeVec; 
  // particles that are in the bottom boundary
  PeriParticlePArray topBoundaryCornerVec; 

  // particles that are in the top boundary of the cubic
  PeriParticlePArray cubicTopBoundaryVec;	

  // peri-bonds between recvPeriParticle in each cpu
  PeriBondPArray recvPeriBondVec; 
  // peri-bonds between periParticleVec in each cpu
  PeriBondPArray periBondVec; 

  //  // for elasticity verification purpose.
  int Uzindex[39]; // particle index along the z direction.[x = 0, y = 0]
  int Uyindex[5];  // particle index along the y direction.[x = 0, z = 0]
  int Uxindex[5];  // particle index along the x direction.[y = 0, z = 0]

};

} // namespace pd

#endif
