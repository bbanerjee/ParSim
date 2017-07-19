#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <Boundary/Boundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/Cylinder.h>
#include <Core/Math/Vec.h>
#include <Core/Math/IntVec.h>
#include <Core/Types/realtypes.h>
#include <DiscreteElements/Contact.h>
#include <DiscreteElements/Containers.h>
#include <DiscreteElements/Gradation.h>
#include <DiscreteElements/Particle.h>
#include <DiscreteElements/Spring.h>
#include <DiscreteElements/Patch.h>
#include <FluidDynamics/Fluid.h>
#include <InputOutput/Parameter.h>
#include <InputOutput/Output.h>
#include <InputOutput/OutputTecplot.h>
#include <InputOutput/OutputVTK.h>
#include <boost/mpi.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <cstddef>
#include <fstream>
#include <map>
#include <memory>
#include <vector>
// include peridynamics files
#include <Peridynamics/PeriBoundaryBond.h>
#include <Peridynamics/PeriContainers.h>
#include <Peridynamics/PeriDEMBond.h>
#include <Peridynamics/PeriElement.h>
#include <Peridynamics/PeriParticle.h>
#include <Peridynamics/globfuncs.h>

#define isProc0_macro ( dem::Assembly::getMPIRank() == 0 )
#define proc0cout if( isProc0_macro ) std::cout
#define proc0cerr if( isProc0_macro ) std::cerr

using pd::PeriElementArray;
using pd::PeriBoundaryBondPArray;
using pd::PeriDEMBondPArray;

namespace dem {

class Assembly
{
public:
  Assembly()
    : trimHistoryNum(0)
    , allContactNum(0)
    , avgNormal(0)
    , avgShear(0)
    , avgPenetr(0)
    , transEnergy(0)
    , rotatEnergy(0)
    , kinetEnergy(0)
    , graviEnergy(0)
    , mechaEnergy(0)
    , vibraTimeStep(0)
    , impactTimeStep(0)
    , nPeriParticle(0)
  {
  }

  ~Assembly()
  {
    // in case of consecutive simulations
    allParticleVec.clear();
    particleVec.clear();
    boundaryVec.clear();
    cavityBoundaryVec.clear();
    springVec.clear();

    // pd part
    allPeriParticleVecInitial.clear();
    allPeriParticleVec.clear(); // this is part of allPeriParticleVecInitial
    interfacePeriParticleVec.clear();
    outerfacePeriParticleVec.clear();
    //    bottomBoundaryVec.clear();
    //    cubicTopBoundaryVec.clear();

    // delete memories in current cpus
    periParticleVec.clear();
    fixedPeriParticleVec.clear(); // this is part of periParticleVec
    bottomBoundaryVec.clear();
    frontBoundaryVec.clear();
    leftBoundaryVec.clear();
    topBoundaryInnerVec.clear();
    topBoundaryEdgeVec.clear();
    topBoundaryCornerVec.clear();

    topBoundaryBondVec.clear();
    bottomBoundaryBondVec.clear();
    leftBoundaryBondVec.clear();
    rightBoundaryBondVec.clear();
    frontBoundaryBondVec.clear();
    backBoundaryBondVec.clear();

    periDEMBondVec.clear();
    recvPeriBondVec.clear();
    periBondVec.clear();
  } // ~Assembly()

  // Accessor methods
  int getMPIRank() const { return mpiRank; }
  boost::mpi::communicator getMPIWorld() const { return boostWorld; }
  const ParticlePArray& getAllParticleVec() const { return allParticleVec; }
  const ParticlePArray& getParticleVec() const { return particleVec; }
  const pd::PeriParticlePArray& getPeriParticleVec() const
  {
    return periParticleVec;
  }
  const pd::PeriParticlePArray& getRecvPeriParticleVec() const
  {
    return recvPeriParticleVec;
  }
  const Gradation& getGradation() const { return gradation; }
  const Box& getAllContainer() const { return allContainer; }
  Fluid& getFluid() { return fluid; }

  void setCommunicator(boost::mpi::communicator& comm);
  void setContainer(Box cont) { allContainer = cont; }
  void setGrid(Box cont) { grid = cont; }
  void setGradation(Gradation grad) { gradation = grad; }

  void generateParticle(std::size_t particleLayers,
                        const std::string& genParticle);
  void buildBoundary(std::size_t boundaryNum, const std::string& boundaryFile);
  void trim(bool toRebuild, const std::string& inputParticle,
            const std::string& trmParticle);
  void deposit(const std::string& boundaryFile,
               const std::string& particleFile);

  bool tractionErrorTol(REAL sigma, std::string type, REAL sigmaX = 0,
                        REAL sigmaY = 0);
  void getStartDimension(REAL& distX, REAL& distY, REAL& distZ);

  void setCavity(Box cav) { cavity = cav; }

  void readBoundary(const std::string& fileName);
  void readParticles(const std::string& fileName);

  void scatterParticle();
  void scatterDEMPeriParticle();
  void commuParticle();
  void commuParticle(const int& iteration);
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
  void removeParticleOutBox();
  void removePeriParticleOutBox();
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

  void createOutputWriter(const std::string& outputFolder, const int& iter) {
    bool writeVTK = true;
    if (writeVTK) {
      d_writer = std::make_unique<OutputVTK>(outputFolder, iter);
    } else {
      d_writer = std::make_unique<OutputTecplot>(outputFolder, iter);
    }
  }

  void updateFileNames(const int& iter, const std::string& extension) {
    d_writer->updateFileNames(iter, extension);
  }
  void updateFileNames(const int& iter) {
    d_writer->updateFileNames(iter);
  }
  std::string getParticleFileName() const {
    return d_writer->getParticleFileName();
  }

  void openDepositProg(std::ofstream& ofs, const std::string& str);
  void printDepositProg(std::ofstream& ofs);
  void openCompressProg(std::ofstream& ofs, const std::string& str);
  void openPeriProgress(std::ofstream& ofs, const std::string& str);
  void printCompressProg(std::ofstream& ofs, REAL distX, REAL distY,
                         REAL distZ);
  void printPeriProgress(std::ofstream& ofs, const int iframe) const;
  void printPeriProgressHalf(std::ofstream& ofs, const int iframe) const;
  void openParticleProg(std::ofstream& ofs, const std::string& str);
  void closeProg(std::ofstream& ofs);

  void trimCavity(bool toRebuild, const std::string& Particlefile,
                  const std::string& cavParticle);
  void readCavityBoundary(const std::string& boundaryfile);
  void buildCavityBoundary(std::size_t existMaxId,
                           const std::string& boundaryfile);
  void findContact(); // detect and resolve contact between particles
  void findContactSingleThread();
  void findContactMultiThread(int numThreads);
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
  void plotParticle() const; // print all particles
  void plotParticle(ParticlePArray& particleVec) const; // print particles info
  void printParticle(const std::string& fileName) const; // print all particles
  void printParticle(const std::string& fileName, ParticlePArray& particleVec) const; // print particles info
  void printBdryContact() const; // print all boundary contact info
  void printMemParticle(
    const std::string& str) const; // print membrane particles
  void plotSpring(
    const std::string& str) const; // print springs in Tecplot format
  void plotBoundary() const;
  void printBoundary() const; // print rigid boundaries info
  void plotGrid() const;
  void plotCavity(const std::string& str) const;
  void checkMembrane(std::vector<REAL>& vx) const;
  void printContact(const std::string& str) const; // print contacts information
  void printCavityBoundary(
    const std::string& str) const; // print cavity boundaries
  void printCavityParticle(std::size_t total, const std::string& str) const;

  // continue to deposit after a cavity is created inside the particle
  // assemblage
  void depositAfterCavity(std::size_t total_steps, std::size_t snapNum,
                          std::size_t interval, const std::string& iniptclfile,
                          const std::string& inibdryfile,
                          const std::string& inicavefile,
                          const std::string& Particlefile,
                          const std::string& contactfile,
                          const std::string& progressfile,
                          const std::string& debugfile);

  // create a specimen by depositing particles into particle boundaries
  void deposit_PtclBdry(Gradation& grad, std::size_t freetype, REAL rsize,
                        std::size_t total_steps, std::size_t snapNum,
                        std::size_t interval, const std::string& iniptclfile,
                        const std::string& Particlefile,
                        const std::string& contactfile,
                        const std::string& progressfile,
                        const std::string& debugfile);

  // scale the assembly with particle boundaries from deposited state until it
  // reaches steady state
  void scale_PtclBdry(
    std::size_t total_steps = 50000, // total_steps
    std::size_t snapNum = 100,       // number of snapNum
    std::size_t interval = 10,       // print interval
    REAL dimn = 0.05,                // dimension of particle-composed-boundary
    REAL rsize = 1.0,                // relative container size
    const std::string& iniptclfile =
      "dep_particle_end", // input file, initial particles
    const std::string& Particlefile =
      "scl_particle", // output file, resulted particles, including snapNum
    const std::string& contactfile =
      "scl_contact", // output file, resulted contacts, including snapNum
    const std::string& progressfile =
      "scl_progress", // output file, statistical info
    const std::string& debugfile = "scl_debug"); // output file, debug info

  // generate particles in space for particle boundaries
  void generate_p(Gradation& grad, const std::string& str, std::size_t freetype,
                  REAL rsize, REAL ht);

  void deGravitation(std::size_t total_steps, std::size_t snapNum,
                     std::size_t interval, bool toRebuild,
                     const std::string& iniptclfile,
                     const std::string& Particlefile,
                     const std::string& contactfile,
                     const std::string& progressfile,
                     const std::string& debugfile);

  // actual deposit function for particle boundaries
  void deposit_p(
    std::size_t total_steps = 50000, // total_steps
    std::size_t snapNum = 100,       // number of snapNum
    std::size_t interval = 10,       // print interval
    REAL dimn = 0.05,                // dimension of particle-composed-boundary
    REAL rsize = 1.0,                // relative container size
    const std::string& iniptclfile =
      "flo_particle_end", // input file, initial particles
    const std::string& Particlefile =
      "dep_particle", // output file, resulted particles, including snapNum
    const std::string& contactfile =
      "dep_contact", // output file, resulted contacts, including snapNum
    const std::string& progressfile =
      "dep_progress", // output file, statistical info
    const std::string& debugfile = "dep_debug"); // output file, debug info

  // squeeze paticles inside a container by moving the boundaries
  void squeeze(
    std::size_t total_steps = 20000, // total_steps
    std::size_t init_steps = 5000,   // initial_steps to reach equilibrium
    std::size_t snapNum = 100,       // number of snapNum
    std::size_t interval = 10,       // print interval
    int flag = -1,                   // -1 squeeze; +1 loosen
    const std::string& iniptclfile =
      "flo_particle_end", // input file, initial particles
    const std::string& inibdryfile =
      "dep_boundary_ini", // input file, initial boundaries
    const std::string& Particlefile =
      "dep_particle", // output file, resulted particles, including snapNum
    const std::string& boundaryfile =
      "dep_boundary", // output file, resulted boundaries
    const std::string& contactfile =
      "dep_contact", // output file, resulted contacts, including snapNum
    const std::string& progressfile =
      "dep_progress", // output file, statistical info
    const std::string& debugfile = "dep_debug"); // output file, debug info

  void deposit_repose(std::size_t interval, const std::string& inibdryfile,
                      const std::string& Particlefile,
                      const std::string& contactfile,
                      const std::string& progressfile,
                      const std::string& debugfile);

  void angleOfRepose(std::size_t interval, const std::string& inibdryfile,
                     const std::string& Particlefile,
                     const std::string& contactfile,
                     const std::string& progressfile,
                     const std::string& debugfile);

  REAL getPtclMinX(const ParticlePArray& particleVec) const;
  REAL getPtclMaxX(const ParticlePArray& particleVec) const;
  REAL getPtclMinY(const ParticlePArray& particleVec) const;
  REAL getPtclMaxY(const ParticlePArray& particleVec) const;
  REAL getPtclMinZ(const ParticlePArray& particleVec) const;
  REAL getPtclMaxZ(const ParticlePArray& particleVec) const;

  void collapse(std::size_t total_steps, std::size_t snapNum,
                std::size_t interval, const std::string& iniptclfile,
                const std::string& initboundary,
                const std::string& Particlefile, const std::string& contactfile,
                const std::string& progressfile, const std::string& debugfile);

  void createMemParticle(REAL rRadius, bool toRebuild,
                         const std::string& Particlefile,
                         const std::string& allParticle);

  void iso_MemBdry(std::size_t total_steps, std::size_t snapNum,
                   std::size_t interval, REAL sigma3, REAL rRadius,
                   bool toRebuild, const std::string& iniptclfile,
                   const std::string& Particlefile,
                   const std::string& contactfile,
                   const std::string& progressfile,
                   const std::string& debugfile);

  void TrimPtclBdryByHeight(REAL height, const std::string& iniptclfile,
                            const std::string& Particlefile);

  void applyParticleBoundary(std::size_t total_steps = 100000,
                             std::size_t snapNum = 100,
                             std::size_t nterval = 10, REAL sigma = 1.0e+4,
                             const std::string& iniptclfile = "cre_particle",
                             const std::string& inibdryfile = "cre_bounary",
                             const std::string& Particlefile = "iso_particle",
                             const std::string& boundaryfile = "iso_boundary",
                             const std::string& contactfile = "iso_contact",
                             const std::string& progressfile = "iso_progress",
                             const std::string& balancedfile = "iso_balanced",
                             const std::string& debugfile = "iso_debug");

  // The confining pressure is 500kPa. This function initializes triaxial
  // compression test.
  void triaxialPtclBdryIni(std::size_t total_steps = 10000,
                           std::size_t snapNum = 100, std::size_t interval = 10,
                           REAL sigma = 5.0e+5,
                           const std::string& iniptclfile = "ini_particle_ini",
                           const std::string& inibdryfile = "ini_boundary_ini",
                           const std::string& Particlefile = "ini_particle",
                           const std::string& boundaryfile = "ini_boundary",
                           const std::string& contactfile = "ini_contact",
                           const std::string& progressfile = "ini_progress",
                           const std::string& debugfile = "ini_debug");

  // The confining pressure is 500kPa. This function performs triaxial
  // compression test.
  // Displacement boundaries are used in axial direction.
  void triaxialPtclBdry(std::size_t total_steps = 100000,
                        std::size_t snapNum = 100, std::size_t interval = 10,
                        const std::string& iniptclfile = "iso_particle_100k",
                        const std::string& inibdryfile = "iso_boundary_100k",
                        const std::string& Particlefile = "tri_particle",
                        const std::string& boundaryfile = "tri_boundary",
                        const std::string& contactfile = "tri_contact",
                        const std::string& progressfile = "tri_progress",
                        const std::string& balancedfile = "tri_balanced",
                        const std::string& debugfile = "tri_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // A rectangular pile is then drived into the particles using displacement
  // control.
  void rectPile_Disp(std::size_t total_steps = 50000, std::size_t snapNum = 100,
                     std::size_t interval = 10,
                     const std::string& iniptclfile = "pile_particle_ini",
                     const std::string& inibdryfile = "pile_boundary_ini",
                     const std::string& Particlefile = "pile_particle",
                     const std::string& boundaryfile = "pile_boundary",
                     const std::string& contactfile = "pile_contact",
                     const std::string& progressfile = "pile_progress",
                     const std::string& debugfile = "pile_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // An ellipsoidal pile is then drived into the particles using displacement
  // control.
  void ellipPile_Disp(std::size_t total_steps = 50000,
                      std::size_t snapNum = 100, std::size_t interval = 10,
                      REAL dimn = 0.05, REAL rsize = 1.0,
                      const std::string& iniptclfile = "pile_particle_ini",
                      const std::string& Particlefile = "pile_particle",
                      const std::string& contactfile = "pile_contact",
                      const std::string& progressfile = "pile_progress",
                      const std::string& debugfile = "pile_debug");

  // The specimen has been deposited with gravitation within rigid boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial
  // velocity.
  void ellipPile_Impact(std::size_t total_steps = 50000,
                        std::size_t snapNum = 100, std::size_t interval = 10,
                        REAL dimn = 0.05,
                        const std::string& iniptclfile = "ipt_particle_ini",
                        const std::string& inibdryfile = "dep_boundary_ini",
                        const std::string& Particlefile = "ipt_particle",
                        const std::string& contactfile = "ipt_contact",
                        const std::string& progressfile = "ipt_progress",
                        const std::string& debugfile = "ipt_debug");

  // The specimen has been deposited with gravitation within particle
  // boundaries.
  // An ellipsoidal penetrator is then impacted into the particles with initial
  // velocity.
  void ellipPile_Impact_p(std::size_t total_steps = 50000,
                          std::size_t snapNum = 100, std::size_t interval = 10,
                          REAL dimn = 0.05,
                          const std::string& iniptclfile = "ipt_particle_ini",
                          const std::string& Particlefile = "ipt_particle",
                          const std::string& contactfile = "ipt_contact",
                          const std::string& progressfile = "ipt_progress",
                          const std::string& debugfile = "ipt_debug");

  // The specimen has been deposited with gravitation within boundaries composed
  // of particles.
  // An ellipsoidal pile is then drived into the particles using force control.
  void ellipPile_Force(std::size_t total_steps = 50000,
                       std::size_t snapNum = 100, std::size_t interval = 10,
                       REAL dimn = 0.05, REAL force = 1.0e+4,
                       std::size_t division = 100,
                       const std::string& iniptclfile = "pile_particle_ini",
                       const std::string& Particlefile = "pile_particle",
                       const std::string& contactfile = "pile_contact",
                       const std::string& progressfile = "pile_progress",
                       const std::string& balancedfile = "pile_balanced",
                       const std::string& debugfile = "pile_debug");

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// pd part
  /////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int getnPeriParticle() const
  {
    return nPeriParticle;
  } // getnPeriParticle - returns number of peri-particles

  void setInitIsv();
  void initialPeriDynamics(
    const std::string&); // initial - initializes the velocity,
                         // displacement and acceleration
  // @param const std::string& - name of the input file
  void prescribeEssentialBoundaryCondition(
    const int); // apply the fixed boundary condition on the bottom particles
  void solve(const std::string&); // solve - solve this problem
  // @param const std::string& - output file name for tecplot visualization
  void runFirstHalfStep();  // run step1 and step2 in Page 5 of Houfu's notes,
                            // equations (15) and (16)
  void runSecondHalfStep(); // run step3 in page 5 of Houfu's notes, equation
                            // (17)
  void constructNeighbor(); // neibhor - searches and constructs the particle
                            // neighborlists
  void findRecvPeriBonds(); // find peri-bonds between periParticleVec and
                            // recvPeriParticleVec
  void findPeriDEMBonds();  // find sand-peri bonds in each cpu, i.e.
                            // periParticleVec and ParticleVec
  void clearPeriDEMBonds();
  void eraseBrokenPeriDEMBonds();
  void readPeriDynamicsData(
    const std::string&); // readData - reads controlling parameters, particle
                         // positions and mesh connectivities
  // @param char * - reference of the input file name
  void writeMesh(const std::string&); // writeMesh - outputs the mesh, used for
                                      // problem checking
  // @param char * - reference of the output file name
  void writeMeshCheckVolume(const std::string&); // writeMeshCheckVolume -
                                                 // outputs the mesh and volume,
                                                 // will
  // be used for volume comuptation checking
  // @param char * - reference of the output file name
  void writeParticleTecplot(
    std::ofstream&,
    const int) const; // writeParticleTecplot - outputs the information of all
                      // particles, for tecplot visualization
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
  void calcHorizonSize(); // calcHorizonSize - calculates the horizon size for
                          // each peri-particle
  void calcParticleVolume();   // calcParticleVolume - calculates the particle
                               // volume associated with each particle
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

  void findParticleInBox(const Box& container,
                         const ParticlePArray& allParticle,
                         ParticlePArray& foundParticle);

  void findPeriParticleInBox(
    const Box& container, const pd::PeriParticlePArray& allParticle,
    pd::PeriParticlePArray& foundParticle);

private:

  // The output writer pointer
  std::unique_ptr<Output> d_writer;

  // particles property
  Gradation
    gradation; // particles gradation, broadcast among processes for once
  ParticlePArray
    allParticleVec;           // all particles, only meaningful to root process
  ParticlePArray particleVec; // particles per process
  std::size_t trimHistoryNum; // historical maximum numbering before trimming,
                              // only meaningful to root process

  ContactArray contactVec; // contacts per process
  ContactTangentArray
    contactTgtVec; // tangential contact force and displacement per process
  std::size_t allContactNum; // estimated total contact number, only meaningful
                             // to root process

  MembraneParticlePArray memBoundary; // membrane particle boundaries
  SpringUPArray springVec;            // springs connecting membrane particles

  // container property
  Box allContainer; // whole container, broadcast among processes for once
  Box container;    // container per process
  Box cavity;       // cavity inside container
  Box grid; // adaptive compute grid, broadcast among processes for once,
            // updated per process

  // boundaries property
  BoundaryPArray
    boundaryVec; // rigid boundaries, broadcast among processes upon changed.
  BoundaryPArray
    mergeBoundaryVec; // rigid boundaries with stats from all processes
  BoundaryPArray cavityBoundaryVec; // rigid cavity boundaries
  std::map<std::size_t, BoundaryTangentArray>
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
  int mpiRank, mpiSize, mpiTag;

  IntVec d_mpiProcs;
  IntVec d_mpiCoords;

  std::unique_ptr<Patch> d_patchP;
  void createPatch(int iteration, const REAL& ghostWidth);
  void updatePatch(int iteration, const REAL& ghostWidth);

  int rankX1, rankX2, rankY1, rankY2, rankZ1, rankZ2;
  int rankX1Y1, rankX1Y2, rankX1Z1, rankX1Z2;
  int rankX2Y1, rankX2Y2, rankX2Z1, rankX2Z2;
  int rankY1Z1, rankY1Z2, rankY2Z1, rankY2Z2;
  int rankX1Y1Z1, rankX1Y1Z2, rankX1Y2Z1, rankX1Y2Z2;
  int rankX2Y1Z1, rankX2Y1Z2, rankX2Y2Z1, rankX2Y2Z2;
  ParticlePArray rParticleX1, rParticleX2; // r stands for received
  ParticlePArray rParticleY1, rParticleY2;
  ParticlePArray rParticleZ1, rParticleZ2;
  ParticlePArray rParticleX1Y1, rParticleX1Y2, rParticleX1Z1, rParticleX1Z2;
  ParticlePArray rParticleX2Y1, rParticleX2Y2, rParticleX2Z1, rParticleX2Z2;
  ParticlePArray rParticleY1Z1, rParticleY1Z2, rParticleY2Z1, rParticleY2Z2;
  ParticlePArray rParticleX1Y1Z1, rParticleX1Y1Z2, rParticleX1Y2Z1,
    rParticleX1Y2Z2;
  ParticlePArray rParticleX2Y1Z1, rParticleX2Y1Z2, rParticleX2Y2Z1,
    rParticleX2Y2Z2;
  ParticlePHashMap sentParticleVec;  // received particles per process
  ParticlePArray recvParticleVec;  // received particles per process
  ParticlePArray mergeParticleVec; // merged particles per process

  // stream
  std::ofstream progressInf;
  std::ofstream balancedInf;
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

  pd::PeriParticlePArray
    allPeriParticleVecInitial; // only for master cpu
  // in simulations, we have peri-points and dem sands, for convience of domain
  // generation, we input a cuboid domain
  // filling with peri-points, this particle vector stores all these
  // peri-points. Then since we have sand particles also
  // in the same assembly, we need to remove those peri-points that are within
  // the sand particles to generate the periParticleVec
  // that are used for the calculation. July 14, 2014
  pd::PeriParticlePArray allPeriParticleVec; // only for master cpu
  pd::PeriParticlePArray
    periParticleVec; // coordinates of all the particles in this cpu, local
                     // peri-points in current cpu
  pd::PeriParticlePArray interfacePeriParticleVec;
  pd::PeriParticlePArray
    fixedPeriParticleVec; // this is the peri-points that are
                          // inside the rigid
                          // inclusion, fix these peri-points
  pd::PeriParticlePArray
    outerfacePeriParticleVec; // this is only for the hollow
                              // spherical example
  //  pd::PeriParticlePArray bottomBoundaryVec;	// particles
  // that are in the bottom boundary
  //  pd::PeriParticlePArray topBoundaryVec;	// particles that
  // are in the bottom boundary
  //  pd::PeriParticlePArray cubicTopBoundaryVec;	// particles
  // that are in the top boundary of the cubic
  //  pd::PeriBondPArray totalBondVec;	// all the Bonds in the
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

  pd::PeriParticlePArray rperiParticleX1,
    rperiParticleX2; // r stands for received
  pd::PeriParticlePArray rperiParticleY1, rperiParticleY2;
  pd::PeriParticlePArray rperiParticleZ1, rperiParticleZ2;
  pd::PeriParticlePArray rperiParticleX1Y1, rperiParticleX1Y2,
    rperiParticleX1Z1, rperiParticleX1Z2;
  pd::PeriParticlePArray rperiParticleX2Y1, rperiParticleX2Y2,
    rperiParticleX2Z1, rperiParticleX2Z2;
  pd::PeriParticlePArray rperiParticleY1Z1, rperiParticleY1Z2,
    rperiParticleY2Z1, rperiParticleY2Z2;
  pd::PeriParticlePArray rperiParticleX1Y1Z1, rperiParticleX1Y1Z2,
    rperiParticleX1Y2Z1, rperiParticleX1Y2Z2;
  pd::PeriParticlePArray rperiParticleX2Y1Z1, rperiParticleX2Y1Z2,
    rperiParticleX2Y2Z1, rperiParticleX2Y2Z2;
  pd::PeriParticlePArray
    recvPeriParticleVec; // received particles per process
  pd::PeriParticlePArray mergePeriParticleVec;

  pd::PeriParticlePArray
    bottomBoundaryVec; // particles that are in the bottom boundary
  pd::PeriParticlePArray
    frontBoundaryVec; // particles that are in the front y boundary
  pd::PeriParticlePArray
    leftBoundaryVec; // particles that are in the left x boundary

  pd::PeriParticlePArray
    topBoundaryInnerVec; // particles that are in the bottom boundary
  pd::PeriParticlePArray
    topBoundaryEdgeVec; // particles that are in the bottom boundary
  pd::PeriParticlePArray
    topBoundaryCornerVec; // particles that are in the bottom boundary

  pd::PeriBondPArray
    recvPeriBondVec; // peri-bonds between recvPeriParticle in each cpu
  pd::PeriBondPArray
    periBondVec; // peri-bonds between periParticleVec in each cpu

  //  // for elasticity verification purpose.
  //  int Uzindex[39]; // particle index along the z direction.[x = 0, y = 0]
  //  int Uyindex[5];  // particle index along the y direction.[x = 0, z = 0]
  //  int Uxindex[5];  // particle index along the x direction.[y = 0, z = 0]

};

} // namespace dem

#endif
