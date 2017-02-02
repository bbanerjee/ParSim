#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <Core/Types/realtypes.h>
#include <InputOutput/Parameter.h>
#include <Core/Math/Vec.h>
#include <DiscreteElements/Gradation.h>
#include <DiscreteElements/Particle.h>
#include <DiscreteElements/Contact.h>
#include <DiscreteElements/Boundary.h>
#include <DiscreteElements/Containers.h>
#include <Core/Geometry/Rectangle.h>
#include <Core/Geometry/Cylinder.h>
#include <DiscreteElements/Spring.h>
#include <FluidDynamics/Fluid.h>
#include <cstddef>
#include <map>
#include <vector>
#include <fstream>
#include <memory>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unique_ptr.hpp>
//include peridynamics files
#include <Peridynamics/globfuncs.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/PeriBoundaryBond.h>
#include <Peridynamics/PeriDEMBond.h>

namespace dem {

class Assembly {
private:

  // particles property
  Gradation gradation;            // particles gradation, broadcast among processes for once
  ParticlePArray allParticleVec;  // all particles, only meaningful to root process
  ParticlePArray particleVec;     // particles per process
  std::size_t trimHistoryNum;     // historical maximum numbering before trimming,
                                  // only meaningful to root process

  ContactArray contactVec;           // contacts per process
  ContactTangentArray contactTgtVec; // tangential contact force and displacement per process
  std::size_t allContactNum;         // estimated total contact number, only meaningful
                                     // to root process

  MembraneParticlePArray memBoundary;  // membrane particle boundaries
  SpringUPArray springVec;             // springs connecting membrane particles

  // container property
  Rectangle allContainer; // whole container, broadcast among processes for once
  Rectangle container;    // container per process
  Rectangle cavity;       // cavity inside container
  Rectangle grid;         // adaptive compute grid, broadcast among processes for once,
                          // updated per process

  // boundaries property
  BoundaryPArray boundaryVec;       // rigid boundaries, broadcast among processes upon changed.
  BoundaryPArray mergeBoundaryVec;  // rigid boundaries with stats from all processes
  BoundaryPArray cavityBoundaryVec; // rigid cavity boundaries
  std::map<std::size_t, BoundaryTangentArray >
      boundaryTgtMap; // particle-boundary contact tangential info

  // fluid property
  Fluid fluid;

  // average data
  REAL avgNormal; // only meaningful to root process
  REAL avgShear;  // only meaningful to root process
  REAL avgPenetr; // only meaningful to root process

  // energy data
  REAL transEnergy; // only meaningful to root process
  REAL rotatEnergy; // only meaningful to root process
  REAL kinetEnergy; // only meaningful to root process
  REAL graviEnergy; // only meaningful to root process
  REAL mechaEnergy; // only meaningful to root process

  // time step
  REAL vibraTimeStep;  // meaningful to all processes
  REAL impactTimeStep; // meaningful to all processes

  // MPI data
  boost::mpi::communicator boostWorld;
  MPI_Comm mpiWorld, cartComm;
  std::vector<std::size_t> bdryProcess;
  int mpiProcX, mpiProcY, mpiProcZ;
  int mpiRank, mpiSize, mpiTag, mpiCoords[3];
  int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
  int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2;
  int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2;
  int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2;
  int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2;
  int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;
  ParticlePArray rParticleX1, rParticleX2; // r stands for received
  ParticlePArray rParticleY1, rParticleY2;
  ParticlePArray rParticleZ1, rParticleZ2;
  ParticlePArray rParticleX1Y1, rParticleX1Y2, rParticleX1Z1,
      rParticleX1Z2;
  ParticlePArray rParticleX2Y1, rParticleX2Y2, rParticleX2Z1,
      rParticleX2Z2;
  ParticlePArray rParticleY1Z1, rParticleY1Z2, rParticleY2Z1,
      rParticleY2Z2;
  ParticlePArray rParticleX1Y1Z1, rParticleX1Y1Z2, rParticleX1Y2Z1,
      rParticleX1Y2Z2;
  ParticlePArray rParticleX2Y1Z1, rParticleX2Y1Z2, rParticleX2Y2Z1,
      rParticleX2Y2Z2;
  ParticlePArray recvParticleVec;  // received particles per process
  ParticlePArray mergeParticleVec; // merged particles per process

  // stream
  std::ofstream progressInf;
  std::ofstream balancedInf;
  std::ofstream periProgInf;
  std::ofstream periProgInfHalf;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// periDynamics part
  /////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int nPeriParticle; // number of PeriParticles in the domain, only for master
                     // cpu
  int nele; // number of elements in the mesh, only for master cpu

  int ndim;           // dimension of the problem, only for master cpu
                      //  int nsteps;	// number of total steps
                      //  int printInterval;	// print interval
  int **connectivity; // mesh connectivity, only for master cpu

  REAL point_interval; // for all cpus, broadcast in scatterDEMPeriParticle()
  REAL maxHorizonSize; // the maximum horizon size of all peri-points, for all
                       // cpus

  std::vector<periDynamics::PeriParticle *>
      allPeriParticleVecInitial; // only for master cpu
  // in simulations, we have peri-points and dem sands, for convience of domain
  // generation, we input a cuboid domain
  // filling with peri-points, this particle vector stores all these
  // peri-points. Then since we have sand particles also
  // in the same assembly, we need to remove those peri-points that are within
  // the sand particles to generate the periParticleVec
                          // that are used for the calculation. July 14, 2014
  std::vector<periDynamics::PeriParticle *>
      allPeriParticleVec; // only for master cpu
  std::vector<periDynamics::PeriParticle *>
      periParticleVec; // coordinates of all the particles in this cpu, local
                       // peri-points in current cpu
  std::vector<periDynamics::PeriParticle *> interfacePeriParticleVec;
  std::vector<periDynamics::PeriParticle *>
      fixedPeriParticleVec; // this is the peri-points that are inside the rigid
                            // inclusion, fix these peri-points
  std::vector<periDynamics::PeriParticle *>
      outerfacePeriParticleVec; // this is only for the hollow spherical example
  //  std::vector<periDynamics::PeriParticle*> bottomBoundaryVec;	// particles
  // that are in the bottom boundary
  //  std::vector<periDynamics::PeriParticle*> topBoundaryVec;	// particles that
  // are in the bottom boundary
  //  std::vector<periDynamics::PeriParticle*> cubicTopBoundaryVec;	// particles
  // that are in the top boundary of the cubic
  //  std::vector<periDynamics::PeriBond*> totalBondVec;	// all the Bonds in the
  // domain

  // for all cpus, these bonds, including boundary bonds, peri-bonds, and
  // peri-DEM-bonds will
  // constructed after scattering within each cpu
  std::vector<PeriBoundaryBond *> bottomBoundaryBondVec;
  std::vector<PeriBoundaryBond *> topBoundaryBondVec;
  std::vector<PeriBoundaryBond *> leftBoundaryBondVec;
  std::vector<PeriBoundaryBond *> rightBoundaryBondVec;
  std::vector<PeriBoundaryBond *> frontBoundaryBondVec;
  std::vector<PeriBoundaryBond *> backBoundaryBondVec;

  std::vector<PeriDEMBond *>
      periDEMBondVec; // the bonds that between the sand particle and
  // the peri-points that are near to this particle, in each cpu

  std::vector<periDynamics::PeriParticle *> rperiParticleX1,
      rperiParticleX2; // r stands for received
  std::vector<periDynamics::PeriParticle *> rperiParticleY1, rperiParticleY2;
  std::vector<periDynamics::PeriParticle *> rperiParticleZ1, rperiParticleZ2;
  std::vector<periDynamics::PeriParticle *> rperiParticleX1Y1,
      rperiParticleX1Y2, rperiParticleX1Z1, rperiParticleX1Z2;
  std::vector<periDynamics::PeriParticle *> rperiParticleX2Y1,
      rperiParticleX2Y2, rperiParticleX2Z1, rperiParticleX2Z2;
  std::vector<periDynamics::PeriParticle *> rperiParticleY1Z1,
      rperiParticleY1Z2, rperiParticleY2Z1, rperiParticleY2Z2;
  std::vector<periDynamics::PeriParticle *> rperiParticleX1Y1Z1,
      rperiParticleX1Y1Z2, rperiParticleX1Y2Z1, rperiParticleX1Y2Z2;
  std::vector<periDynamics::PeriParticle *> rperiParticleX2Y1Z1,
      rperiParticleX2Y1Z2, rperiParticleX2Y2Z1, rperiParticleX2Y2Z2;
  std::vector<periDynamics::PeriParticle *>
      recvPeriParticleVec; // received particles per process
  std::vector<periDynamics::PeriParticle *> mergePeriParticleVec;

  std::vector<periDynamics::PeriParticle *>
      bottomBoundaryVec; // particles that are in the bottom boundary
  std::vector<periDynamics::PeriParticle *>
      frontBoundaryVec;  // particles that are in the front y boundary
  std::vector<periDynamics::PeriParticle *>
      leftBoundaryVec;   // particles that are in the left x boundary

  std::vector<periDynamics::PeriParticle *>
      topBoundaryInnerVec;  // particles that are in the bottom boundary
  std::vector<periDynamics::PeriParticle *>
      topBoundaryEdgeVec;   // particles that are in the bottom boundary
  std::vector<periDynamics::PeriParticle *>
      topBoundaryCornerVec; // particles that are in the bottom boundary

  std::vector<periDynamics::PeriBond *>
      recvPeriBondVec; // peri-bonds between recvPeriParticle in each cpu
  std::vector<periDynamics::PeriBond *>
      periBondVec;     // peri-bonds between periParticleVec in each cpu

  //  // for elasticity verification purpose.
  //  int Uzindex[39]; // particle index along the z direction.[x = 0, y = 0]
  //  int Uyindex[5];  // particle index along the y direction.[x = 0, z = 0]
  //  int Uxindex[5];  // particle index along the x direction.[y = 0, z = 0]

public:
  Assembly()
      : trimHistoryNum(0), allContactNum(0), avgNormal(0), avgShear(0),
        avgPenetr(0), transEnergy(0), rotatEnergy(0), kinetEnergy(0),
        graviEnergy(0), mechaEnergy(0), vibraTimeStep(0), impactTimeStep(0),
        nPeriParticle(0) {}

  ~Assembly() {
    // release memory pointed to by pointers in the container
    //for (auto it = allParticleVec.begin(); it != allParticleVec.end(); ++it)
    //  delete (*it);

    //for (auto it = particleVec.begin(); it != particleVec.end(); ++it)
    //  delete (*it);

    //for (auto it = boundaryVec.begin(); it != boundaryVec.end(); ++it)
    //  delete (*it);

    //for (auto it = cavityBoundaryVec.begin(); it != cavityBoundaryVec.end(); ++it)
    //  delete (*it);

    //for (auto it = springVec.begin(); it != springVec.end(); ++it)
    //  delete (*it);

    // in case of consecutive simulations
    allParticleVec.clear();
    particleVec.clear();
    boundaryVec.clear();
    cavityBoundaryVec.clear();
    springVec.clear();

    // periDynamics part
    // free the spaces of these pointer vector
    for (std::vector<periDynamics::PeriParticle *>::iterator pt =
             allPeriParticleVecInitial.begin();
         pt != allPeriParticleVecInitial.end(); pt++) {
      delete (*pt);
    }
    delete connectivity;
    allPeriParticleVecInitial.clear();
    allPeriParticleVec.clear(); // this is part of allPeriParticleVecInitial
    interfacePeriParticleVec.clear();
    outerfacePeriParticleVec.clear();
    //    bottomBoundaryVec.clear();
    //    cubicTopBoundaryVec.clear();

    // delete memories in current cpus
    for (std::vector<periDynamics::PeriParticle *>::iterator pt =
             periParticleVec.begin();
         pt != periParticleVec.end(); pt++) {
      delete (*pt);
    }
    periParticleVec.clear();
    fixedPeriParticleVec.clear(); // this is part of periParticleVec
    bottomBoundaryVec.clear();
    frontBoundaryVec.clear();
    leftBoundaryVec.clear();
    topBoundaryInnerVec.clear();
    topBoundaryEdgeVec.clear();
    topBoundaryCornerVec.clear();

    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             topBoundaryBondVec.begin();
         pt != topBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }
    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             bottomBoundaryBondVec.begin();
         pt != bottomBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }
    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             leftBoundaryBondVec.begin();
         pt != leftBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }
    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             rightBoundaryBondVec.begin();
         pt != rightBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }
    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             frontBoundaryBondVec.begin();
         pt != frontBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }
    for (std::vector<dem::PeriBoundaryBond *>::iterator pt =
             backBoundaryBondVec.begin();
         pt != backBoundaryBondVec.end(); pt++) {
      delete (*pt);
    }

    topBoundaryBondVec.clear();
    bottomBoundaryBondVec.clear();
    leftBoundaryBondVec.clear();
    rightBoundaryBondVec.clear();
    frontBoundaryBondVec.clear();
    backBoundaryBondVec.clear();

    for (std::vector<dem::PeriDEMBond *>::iterator pt = periDEMBondVec.begin();
         pt != periDEMBondVec.end(); pt++) {
      delete (*pt);
    }
    periDEMBondVec.clear();

    for (std::vector<periDynamics::PeriBond *>::iterator pt =
             recvPeriBondVec.begin();
         pt != recvPeriBondVec.end(); pt++) {
      delete (*pt);
    }
    recvPeriBondVec.clear();

    for (std::vector<periDynamics::PeriBond *>::iterator pt =
             periBondVec.begin();
         pt != periBondVec.end(); pt++) {
      delete (*pt);
    }
    periBondVec.clear();
  } // ~Assembly()

  void setCommunicator(boost::mpi::communicator &comm);
  void setContainer(Rectangle cont) { allContainer = cont; }
  void setGrid(Rectangle cont) { grid = cont; }
  void setGradation(Gradation grad) { gradation = grad; }

  void tuneMassPercentage();
  void depositIntoContainer();
  void resumeDepositIntoContainer();
  void expandCavityParticle();
  void resumeExpandCavityParticle();
  void generateParticle(std::size_t particleLayers, const char *genParticle);
  void buildBoundary(std::size_t boundaryNum, const char *boundaryFile);
  void trimOnly();
  void trim(bool toRebuild, const char *inputParticle, const char *trmParticle);
  void deposit(const char *inputBoundary, const char *inputParticle);
  void proceedFromPreset();
  void coupleWithSonicFluid();

  void isotropic();
  void odometer();
  void triaxial();
  void planeStrain();
  void trueTriaxial();
  bool tractionErrorTol(REAL sigma, std::string type, REAL sigmaX = 0,
                        REAL sigmaY = 0);
  void getStartDimension(REAL &distX, REAL &distY, REAL &distZ);

  void setCavity(Rectangle cav) { cavity = cav; }

  void readParticle(const char *str);
  void readBoundary(const char *str);
  void scatterParticle();
  void scatterDEMPeriParticle();
  void commuParticle();
  void commuPeriParticle();
  bool isBdryProcess();
  void releaseRecvParticle();
  void releaseRecvPeriParticle();
  void releasePeriBondVec();
  void releaseGatheredParticle();
  void releaseGatheredPeriParticle();
  void releaseGatheredContact();
  void migrateParticle();
  void migratePeriParticle();
  void removeParticleOutRectangle();
  void removePeriParticleOutRectangle();
  void gatherParticle();
  void gatherPeriParticle();
  void gatherBdryContact();
  void findFixedPeriParticles();    // for all cpus
  void findBoundaryPeriParticles(); // for all cpus
  void applyPeriBoundaryCondition();
  void applyTractionBoundary(int);

  void updateGrid();
  void updateGridMinX();
  void updateGridMaxX();
  void updateGridMinY();
  void updateGridMaxY();
  void updateGridMinZ();
  void updateGridMaxZ();
  void updatePeriGrid();

  void openDepositProg(std::ofstream &ofs, const char *str);
  void printDepositProg(std::ofstream &ofs);
  void openCompressProg(std::ofstream &ofs, const char *str);
  void openPeriProgress(std::ofstream &ofs, const char *str);
  void printCompressProg(std::ofstream &ofs, REAL distX, REAL distY,
                         REAL distZ);
  void printPeriProgress(std::ofstream &ofs, const int iframe) const;
  void printPeriProgressHalf(std::ofstream &ofs, const int iframe) const;
  void openParticleProg(std::ofstream &ofs, const char *str);
  void closeProg(std::ofstream &ofs);

  void trimCavity(bool toRebuild, const char *Particlefile,
                  const char *cavParticle);
  void readCavityBoundary(const char *boundaryfile);
  void buildCavityBoundary(std::size_t existMaxId, const char *boundaryfile);
  void findContact();          // detect and resolve contact between particles
  void findBdryContact();      // find particles on boundaries
  void findParticleOnCavity(); // find particle on cavity boundaries

  void clearContactForce(); // clear forces and moments for all particles
  void internalForce();     // calculate inter-particle forces
  void springForce();
  void boundaryForce(); // calcualte forces between rigid boundaries and
                        // particles
  void cavityBoundaryForce();
  void updateParticle(); // update motion of particles

  REAL ellipPileForce();  // for force pile only
  void ellipPileUpdate(); // for force pile only

  Vec ellipPileDimn();
  REAL ellipPileTipZ();
  REAL ellipPilePeneVol();

  void updateBoundary(REAL simga, std::string type, REAL sigmaX = 0,
                      REAL sigmaY = 0);

  REAL getMass() const;
  REAL getAvgPenetr() const;
  REAL getParticleVolume() const;

  void calcTimeStep();
  void calcVibraTimeStep();
  void calcImpactTimeStep();
  void calcContactNum();

  REAL getAvgTransVelocity() const;
  REAL getAvgRotatVelocity() const;
  REAL getAvgForce() const;
  REAL getAvgMoment() const;

  void calcTransEnergy();
  void calcRotatEnergy();
  void calcKinetEnergy();
  void calcGraviEnergy(REAL ref);
  void calcMechaEnergy();
  void gatherEnergy();

  void setTrimHistoryNum(std::size_t n) { trimHistoryNum = n; }
  void printParticle(const char *str) const; // print all particles
  void printBdryContact(
      const char *str) const;                // print all boundary contact info
  void printParticle(
      const char *str,
      ParticlePArray &particleVec) const; // print particles info
  void printMemParticle(const char *str) const;    // print membrane particles
  void plotSpring(const char *str) const; // print springs in Tecplot format
  void plotBoundary(const char *str) const;
  void plotGrid(const char *str) const;
  void plotCavity(const char *str) const;
  void checkMembrane(std::vector<REAL> &vx) const;
  void printContact(char *str) const;        // print contacts information
  void printBoundary(const char *str) const; // print rigid boundaries info
  void printCavityBoundary(const char *str) const; // print cavity boundaries
  void printCavityParticle(std::size_t total, const char *str) const;

  // continue to deposit after a cavity is created inside the particle
  // assemblage
  void depositAfterCavity(std::size_t total_steps, std::size_t snapNum,
                          std::size_t interval, const char *iniptclfile,
                          const char *inibdryfile, const char *inicavefile,
                          const char *Particlefile, const char *contactfile,
                          const char *progressfile, const char *debugfile);

  // create a specimen by depositing particles into particle boundaries
  void deposit_PtclBdry(Gradation &grad, std::size_t freetype, REAL rsize,
                        std::size_t total_steps, std::size_t snapNum,
                        std::size_t interval, const char *iniptclfile,
                        const char *Particlefile, const char *contactfile,
                        const char *progressfile, const char *debugfile);

  // scale the assembly with particle boundaries from deposited state until it
  // reaches steady state
  void scale_PtclBdry(
      std::size_t total_steps = 50000, // total_steps
      std::size_t snapNum = 100,       // number of snapNum
      std::size_t interval = 10,       // print interval
      REAL dimn = 0.05,       // dimension of particle-composed-boundary
      REAL rsize = 1.0,       // relative container size
      const char *iniptclfile =
          "dep_particle_end", // input file, initial particles
      const char *Particlefile =
          "scl_particle", // output file, resulted particles, including snapNum
      const char *contactfile =
          "scl_contact",  // output file, resulted contacts, including snapNum
      const char *progressfile =
          "scl_progress", // output file, statistical info
      const char *debugfile = "scl_debug"); // output file, debug info

  // generate particles in space for particle boundaries
  void generate_p(Gradation &grad, const char *str, std::size_t freetype,
                  REAL rsize, REAL ht);

  void deGravitation(std::size_t total_steps, std::size_t snapNum,
                     std::size_t interval, bool toRebuild,
                     const char *iniptclfile, const char *Particlefile,
                     const char *contactfile, const char *progressfile,
                     const char *debugfile);

  // actual deposit function for particle boundaries
  void deposit_p(
      std::size_t total_steps = 50000, // total_steps
      std::size_t snapNum = 100,       // number of snapNum
      std::size_t interval = 10,       // print interval
      REAL dimn = 0.05,       // dimension of particle-composed-boundary
      REAL rsize = 1.0,       // relative container size
      const char *iniptclfile =
          "flo_particle_end", // input file, initial particles
      const char *Particlefile =
          "dep_particle", // output file, resulted particles, including snapNum
      const char *contactfile =
          "dep_contact",  // output file, resulted contacts, including snapNum
      const char *progressfile =
          "dep_progress", // output file, statistical info
      const char *debugfile = "dep_debug"); // output file, debug info

  //squeeze paticles inside a container by moving the boundaries
  void squeeze(
      std::size_t total_steps = 20000, // total_steps
      std::size_t init_steps = 5000,   // initial_steps to reach equilibrium
      std::size_t snapNum = 100,       // number of snapNum
      std::size_t interval = 10,       // print interval
      int flag = -1,                   // -1 squeeze; +1 loosen
      const char *iniptclfile =
          "flo_particle_end",          // input file, initial particles
      const char *inibdryfile =
          "dep_boundary_ini",          // input file, initial boundaries
      const char *Particlefile =
          "dep_particle", // output file, resulted particles, including snapNum
      const char *boundaryfile =
          "dep_boundary", // output file, resulted boundaries
      const char *contactfile =
          "dep_contact",  // output file, resulted contacts, including snapNum
      const char *progressfile =
          "dep_progress", // output file, statistical info
      const char *debugfile = "dep_debug"); // output file, debug info

  void deposit_repose(std::size_t interval, const char *inibdryfile,
                      const char *Particlefile, const char *contactfile,
                      const char *progressfile, const char *debugfile);

  void angleOfRepose(std::size_t interval, const char *inibdryfile,
                     const char *Particlefile, const char *contactfile,
                     const char *progressfile, const char *debugfile);

  REAL getPtclMinX(const ParticlePArray &particleVec) const;
  REAL getPtclMaxX(const ParticlePArray &particleVec) const;
  REAL getPtclMinY(const ParticlePArray &particleVec) const;
  REAL getPtclMaxY(const ParticlePArray &particleVec) const;
  REAL getPtclMinZ(const ParticlePArray &particleVec) const;
  REAL getPtclMaxZ(const ParticlePArray &particleVec) const;

  void collapse(std::size_t total_steps, std::size_t snapNum,
                std::size_t interval, const char *iniptclfile,
                const char *initboundary, const char *Particlefile,
                const char *contactfile, const char *progressfile,
                const char *debugfile);

  void createMemParticle(REAL rRadius, bool toRebuild, const char *Particlefile,
                         const char *allParticle);

  void iso_MemBdry(std::size_t total_steps, std::size_t snapNum,
                   std::size_t interval, REAL sigma3, REAL rRadius,
                   bool toRebuild, const char *iniptclfile,
                   const char *Particlefile, const char *contactfile,
                   const char *progressfile, const char *debugfile);

  void TrimPtclBdryByHeight(REAL height, const char *iniptclfile,
                            const char *Particlefile);

  void applyParticleBoundary(std::size_t total_steps = 100000,
                             std::size_t snapNum = 100,
                             std::size_t nterval = 10, REAL sigma = 1.0e+4,
                             const char *iniptclfile = "cre_particle",
                             const char *inibdryfile = "cre_bounary",
                             const char *Particlefile = "iso_particle",
                             const char *boundaryfile = "iso_boundary",
                             const char *contactfile = "iso_contact",
                             const char *progressfile = "iso_progress",
                             const char *balancedfile = "iso_balanced",
                             const char *debugfile = "iso_debug");

  // The confining pressure is 500kPa. This function initializes triaxial
  // compression test.
  void triaxialPtclBdryIni(std::size_t total_steps = 10000,
                           std::size_t snapNum = 100, std::size_t interval = 10,
                           REAL sigma = 5.0e+5,
                           const char *iniptclfile = "ini_particle_ini",
                           const char *inibdryfile = "ini_boundary_ini",
                           const char *Particlefile = "ini_particle",
                           const char *boundaryfile = "ini_boundary",
                           const char *contactfile = "ini_contact",
                           const char *progressfile = "ini_progress",
                           const char *debugfile = "ini_debug");

  // The confining pressure is 500kPa. This function performs triaxial
  // compression test.
  // Displacement boundaries are used in axial direction.
  void triaxialPtclBdry(std::size_t total_steps = 100000,
                        std::size_t snapNum = 100, std::size_t interval = 10,
                        const char *iniptclfile = "iso_particle_100k",
                        const char *inibdryfile = "iso_boundary_100k",
                        const char *Particlefile = "tri_particle",
                        const char *boundaryfile = "tri_boundary",
                        const char *contactfile = "tri_contact",
                        const char *progressfile = "tri_progress",
                        const char *balancedfile = "tri_balanced",
                        const char *debugfile = "tri_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // A rectangular pile is then drived into the particles using displacement
  // control.
  void rectPile_Disp(std::size_t total_steps = 50000, std::size_t snapNum = 100,
                     std::size_t interval = 10,
                     const char *iniptclfile = "pile_particle_ini",
                     const char *inibdryfile = "pile_boundary_ini",
                     const char *Particlefile = "pile_particle",
                     const char *boundaryfile = "pile_boundary",
                     const char *contactfile = "pile_contact",
                     const char *progressfile = "pile_progress",
                     const char *debugfile = "pile_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // An ellipsoidal pile is then drived into the particles using displacement
  // control.
  void ellipPile_Disp(std::size_t total_steps = 50000,
                      std::size_t snapNum = 100, std::size_t interval = 10,
                      REAL dimn = 0.05, REAL rsize = 1.0,
                      const char *iniptclfile = "pile_particle_ini",
                      const char *Particlefile = "pile_particle",
                      const char *contactfile = "pile_contact",
                      const char *progressfile = "pile_progress",
                      const char *debugfile = "pile_debug");

  // The specimen has been deposited with gravitation within rigid boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial
  // velocity.
  void ellipPile_Impact(std::size_t total_steps = 50000,
                        std::size_t snapNum = 100, std::size_t interval = 10,
                        REAL dimn = 0.05,
                        const char *iniptclfile = "ipt_particle_ini",
                        const char *inibdryfile = "dep_boundary_ini",
                        const char *Particlefile = "ipt_particle",
                        const char *contactfile = "ipt_contact",
                        const char *progressfile = "ipt_progress",
                        const char *debugfile = "ipt_debug");

  // The specimen has been deposited with gravitation within particle
  // boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial
  // velocity.
  void ellipPile_Impact_p(std::size_t total_steps = 50000,
                          std::size_t snapNum = 100, std::size_t interval = 10,
                          REAL dimn = 0.05,
                          const char *iniptclfile = "ipt_particle_ini",
                          const char *Particlefile = "ipt_particle",
                          const char *contactfile = "ipt_contact",
                          const char *progressfile = "ipt_progress",
                          const char *debugfile = "ipt_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // An ellipsoidal pile is then drived into the particles using force control.
  void ellipPile_Force(std::size_t total_steps = 50000,
                       std::size_t snapNum = 100, std::size_t interval = 10,
                       REAL dimn = 0.05, REAL force = 1.0e+4,
                       std::size_t division = 100,
                       const char *iniptclfile = "pile_particle_ini",
                       const char *Particlefile = "pile_particle",
                       const char *contactfile = "pile_contact",
                       const char *progressfile = "pile_progress",
                       const char *balancedfile = "pile_balanced",
                       const char *debugfile = "pile_debug");

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// periDynamics part
  /////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int getnPeriParticle() const {
    return nPeriParticle;
  } // getnPeriParticle - returns number of peri-particles

  void setInitIsv();
  void initialPeriDynamics(const char *); // initial - initializes the velocity,
                                          // displacement and acceleration
  // @param const char* - name of the input file
  void prescribeEssentialBoundaryCondition(
      const int); // apply the fixed boundary condition on the bottom particles
  void solve(const char *); // solve - solve this problem
  // @param const char* - output file name for tecplot visualization
  void runFirstHalfStep(); // run step1 and step2 in Page 5 of Houfu's notes,
                           // equations (15) and (16)
  void runSecondHalfStep(); // run step3 in page 5 of Houfu's notes, equation
                            // (17)
  void constructNeighbor(); // neibhor - searches and constructs the particle
                            // neighborlists
  void findRecvPeriBonds(); // find peri-bonds between periParticleVec and
                            // recvPeriParticleVec
  void findPeriDEMBonds(); // find sand-peri bonds in each cpu, i.e.
                           // periParticleVec and ParticleVec
  void clearPeriDEMBonds();
  void eraseBrokenPeriDEMBonds();
  void readPeriDynamicsData(
      const char *); // readData - reads controlling parameters, particle
                     // positions and mesh connectivities
  // @param char * - reference of the input file name
  void writeMesh(
      const char *); // writeMesh - outputs the mesh, used for problem checking
                     // @param char * - reference of the output file name
  void writeMeshCheckVolume(
      const char *); // writeMeshCheckVolume - outputs the mesh and volume, will
                     // be used for volume comuptation checking
                     // @param char * - reference of the output file name
  void writeParticleTecplot(
      std::ofstream &,
      const int) const; // writeParticleTecplot - outputs the information of all
                        // particles, for tecplot visualization
  // @param (std::ofstream, int) - (reference of the output file name, frame
  // index)
  void printPeriDomain(const char *) const;
  void printRecvPeriDomain(const char *) const;
  void printPeriParticle(const char *str) const;

  void printPeriDomainSphere(
      const char *) const; // print stress in spherical coordinates

  void constructPeriMatrix(); // construct Matrix members in periParticleVec
  void constructRecvPeriMatrix(); // construct Matrix members in
                                  // recvPeriParticleVec, construction here
                                  // since currently
  // the pointer array in Matrix cannot be transfered well between cpus
  void calcDeformationGradient(); // calcDeformationGradient - calculates the
                                  // deformation gradient for each peri-particle
  void calcHorizonSize(); // calcHorizonSize - calculates the horizon size for
                          // each peri-particle
  void calcParticleVolume(); // calcParticleVolume - calculates the particle
                             // volume associated with each particle
  void calcParticleKinv();   // calculate the inverse of K for all particles
  void calcRecvParticleKinv(); // calculate the inverse of K for all
                               // recvPeriParticles
  void calcParticleStress(); // calculate the Cauchy Stress for all particles
  void calcParticleAcceleration(); // calculate the acceleration for all
                                   // particles
  void checkBondParticleAlive(); // for each particle, check the state(alive or
                                 // not) of each surrounding bond;
  // if there's no alive bond for a particle, then this particle is disabled.
  void ApplyExternalForce(int istep);
  void writeDisplacementData(const char *, const char *, const char *);

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

public:
  void findParticleInRectangle(const Rectangle &container,
                               const ParticlePArray &allParticle,
                               ParticlePArray &foundParticle);

  void findPeriParticleInRectangle(
      const Rectangle &container,
      const std::vector<periDynamics::PeriParticle *> &allParticle,
      std::vector<periDynamics::PeriParticle *> &foundParticle);

};

} // namespace dem

#endif
