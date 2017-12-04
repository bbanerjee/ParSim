#ifndef DEM_H
#define DEM_H

#include <Boundary/Boundary.h>
#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Geometry/Cylinder.h>
#include <Core/Math/Vec.h>
#include <Core/Math/IntVec.h>
#include <Core/Types/RealTypes.h>
#include <DiscreteElements/DEMBoundaryConditions.h>
#include <DiscreteElements/DEMContact.h>
#include <DiscreteElements/DEMContainers.h>
#include <DiscreteElements/Gradation.h>
#include <DiscreteElements/DEMParticle.h>
#include <DiscreteElements/Spring.h>
#include <Core/Parallel/Patch.h>
#include <FluidDynamics/Fluid.h>
#include <InputOutput/InputParameter.h>
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

#define isProc0_macro ( dem::DiscreteElements::getMPIRank() == 0 )
#define proc0cout if( isProc0_macro ) std::cout
#define proc0cerr if( isProc0_macro ) std::cerr

namespace dem {

class DiscreteElements
{
public:

  DiscreteElements()
    : d_minParticleRadius(1.0e20)
    , d_maxParticleRadius(-1.0e20)
    , d_trimHistoryNum(0)
    , d_allContactNum(0)
    , d_avgNormalForce(0)
    , d_avgShearForce(0)
    , d_avgPenetration(0)
    , d_overlapCount(0)
    , d_translationalEnergy(0)
    , d_rotationalEnergy(0)
    , d_kineticEnergy(0)
    , d_gravitationalEnergy(0)
    , d_mechanicalEnergy(0)
    , d_vibrationTimeStep(0)
    , d_impactTimeStep(0)
  {
  }

  ~DiscreteElements()
  {
    // in case of consecutive simulations
    d_allDEMParticles.clear();
    d_patchParticles.clear();
    d_boundaries.clear();
    d_cavityBoundaries.clear();
    d_springs.clear();
  } // ~DiscreteElements()

  // Accessor methods
  boost::mpi::communicator getMPIWorld() const { return boostWorld; }
  static MPI_Comm getMPIWorldComm() { return s_mpiWorld; }
  static MPI_Comm getMPICartComm() { return s_cartComm; }
  static int getMPIRank() { return s_mpiRank; }
  static int getMPISize() { return s_mpiSize; }
  static IntVec getMPIProcs() { return s_mpiProcs; }
  static IntVec getMPICoords() { return s_mpiCoords; }

  const DEMParticlePArray& getAllDEMParticleVec() const { return d_allDEMParticles; }
  DEMParticlePArray& getModifiableAllParticleVec() { return d_allDEMParticles; }

  const DEMParticlePArray& getDEMParticleVec() const { return d_patchParticles; }
  DEMParticlePArray& getModifiableParticleVec() { return d_patchParticles; }
  DEMParticlePArray& getMergedParticleVec() { return d_mergedParticles; }

  std::size_t numOverlappingParticles() const
  {
    return d_overlapCount;
  }

  const Gradation& getGradation() const { return d_gradation; }
  const Box& getSpatialDomain() const { return d_spatialDomain; }
  Box& getModifiableSpatialDomain() { return d_spatialDomain; }
  Fluid& getFluid() { return fluid; }

  void setCommunicator(boost::mpi::communicator& comm);
  void setSpatialDomain(Box cont) { d_spatialDomain = cont; }
  void setPatchBox(Box cont) { d_demPatchBox = cont; }
  void setGradation(Gradation grad) { d_gradation = grad; }
  void setAllDEMParticleVec(const DEMParticlePArray& particles) {
    d_allDEMParticles = particles;
  }

  /**
   * Set min and max particle radius
   */
  void setMinMaxParticleRadius()
  {
    if (d_gradation.getPtclMaxRadius() == -1) {
      auto minmaxIter = std::minmax_element(d_allDEMParticles.begin(), 
                                            d_allDEMParticles.end(),
      [](const DEMParticleP& p1, const DEMParticleP& p2){
        auto max_p1 = std::max({p1->radiusA(), p1->radiusB(), p1->radiusC()});
        auto max_p2 = std::max({p2->radiusA(), p2->radiusB(), p2->radiusC()});
        return max_p1 < max_p2;
      });
      d_minParticleRadius = std::min({(*minmaxIter.first)->radiusA(),
                                      (*minmaxIter.first)->radiusB(),
                                      (*minmaxIter.first)->radiusC()});
      d_maxParticleRadius = std::max({(*minmaxIter.second)->radiusA(),
                                      (*minmaxIter.second)->radiusB(),
                                      (*minmaxIter.second)->radiusC()});
    } else {
      d_minParticleRadius = d_gradation.getPtclMinRadius();
      d_maxParticleRadius = d_gradation.getPtclMaxRadius();
    }
  }

  /**
   * 1) Reads the bounding box from an input boundary file
   * 2) Creates internal boundaries form DEM contact calculations
   * 3) Creates particles
   * 4) Writes the particle data and internal boundary data to the
   *    output files
   */
  void createAndSaveParticlesAndBoundaries(const std::string& boundaryFilename,
                                           const std::string& particleFilename);

  void trim(bool toRebuild, const std::string& inputParticle,
            const std::string& trmParticle);
  void deposit(const std::string& boundaryFilename,
               const std::string& particleFilename);

  bool areBoundaryTractionsEquilibrated(REAL time);
  bool areBoundaryTractionsEquilibrated(REAL sigma, std::string type, 
                                        REAL sigmaX = 0, REAL sigmaY = 0);

  void getStartDimension(REAL& distX, REAL& distY, REAL& distZ);

  void setCavity(Box cav) { cavity = cav; }

  void readBoundary(const std::string& fileName);
  void updateBoundaryAreas(BoundaryPArray& boundaries);

  void readParticles(const std::string& fileName);
  void readBoundaryConditions(const std::string& fileName);

  void scatterParticles();
  void scatterDEMPeriParticle();
  void communicateGhostParticles();
  void communicateGhostParticles(std::size_t iteration);
  bool isBoundaryProcess();
  void releaseReceivedParticles();
  void releaseGatheredParticle();
  void releaseGatheredContact();
  void migrateParticles(std::size_t iteration);
  void removeParticleOutBox();
  void gatherParticles();
  void gatherBoundaryContacts();
  std::size_t findBoundaryPeriParticles(); // for all cpus
  void applyTractionBoundary(int);

  void updatePatchBox();
  void updatePatchBoxMinX();
  void updatePatchBoxMaxX();
  void updatePatchBoxMinY();
  void updatePatchBoxMaxY();
  void updatePatchBoxMinZ();
  void updatePatchBoxMaxZ();
  void updatePeriPatchGrid();

  void createOutputWriter(const std::string& outputFolder, const int& iter) {
    bool writeVTK = true;
    if (writeVTK) {
      d_writer = std::make_unique<OutputVTK<DEMParticlePArray>>(outputFolder, iter);
    } else {
      d_writer = std::make_unique<OutputTecplot<DEMParticlePArray>>(outputFolder, iter);
    }
  }

  void updateFileNames(const int& iter, const std::string& extension) {
    d_writer->updateFilenames(iter, extension);
  }
  void updateFileNames(const int& iter) {
    d_writer->updateFilenames(iter);
  }
  std::string getParticleFileName() const {
    return d_writer->getParticleFilename();
  }

  void openProgressOutputFile(std::ofstream& ofs, const std::string& str);
  void openPeriProgress(std::ofstream& ofs, const std::string& str);
  void appendToProgressOutputFile(std::ofstream& ofs, 
                                  std::size_t iteration, REAL timeStep,
                                  REAL distX = 1, REAL distY = 1,
                                  REAL distZ = 1);
  void openParticleProg(std::ofstream& ofs, const std::string& str);
  void closeProgressOutputFile(std::ofstream& ofs);

  void trimCavity(bool toRebuild, const std::string& Particlefile,
                  const std::string& cavParticle);
  void readCavityBoundary(const std::string& boundaryfile);
  void buildCavityBoundary(std::size_t existMaxId,
                           const std::string& boundaryfile);

  // detect and resolve contact between particles
  void findContact(std::size_t iteration); 
  void findContactSingleThread(REAL minOverlap, REAL measOverlap,
                               std::size_t iteration);
  void findContactMultiThread(int numThreads, 
                              REAL minOverlap, REAL measOverlap,
                              std::size_t iteration);
  std::size_t findBoundaryContacts(std::size_t iteration);      // find particles on boundaries
  void findParticleOnCavity(); // find particle on cavity boundaries

  void initializeForces();

  void clearContactForce(); // clear forces and moments for all particles
  void applyBodyForce(); // Apply body forces

  // calculate inter-particle forces
  void internalForce(REAL timeStep, std::size_t iteration);     

  void springForce();

  // calculate forces between rigid boundaries and // particles
  void boundaryForce(REAL timeStep, std::size_t iteration); 
  void cavityBoundaryForce();

  // update motion of particles
  void updateParticles(REAL timeStep, std::size_t iteration); 

  REAL ellipPileForce();  // for force pile only
  void ellipPileUpdate(); // for force pile only

  Vec ellipPileDimn();
  REAL ellipPileTipZ();
  REAL ellipPilePeneVol();

  void updateBoundary(REAL time, REAL delT, REAL mass);
  void updateBoundary(REAL simga, std::string type, REAL sigmaX = 0,
                      REAL sigmaY = 0);

  REAL getTotalMassFromPatchParticleData() const;
  REAL mass() const;
  REAL getAvgPenetration() const;
  REAL volume() const;

  REAL calcTimeStep(REAL curTimeStep);
  void calcVibrationTimeStep();
  void calcImpactTimeStep();
  void calcContactNum();

  REAL getAvgTranslationalVelocity() const;
  REAL getAvgRotationalVelocity() const;
  REAL getAvgForce() const;
  REAL getAvgMoment() const;

  void calcTranslationalEnergy();
  void calcRotationalEnergy();
  void calcKineticEnergy();
  void calcGravitationalEnergy(REAL ref);
  void calcMechanicalEnergy();
  void gatherEnergy();

  void setTrimHistoryNum(std::size_t n) { d_trimHistoryNum = n; }
  void writeParticlesToFile(int frame) const; // print all particles
  void writeParticlesToFile(DEMParticlePArray& d_patchParticles, int frame) const; // print particles info
  void printParticlesCSV(const std::string& folderName,
                         const std::string& fileName, 
                         int frame) const;
  void printParticlesCSV(const std::string& folderName,
                         const std::string& fileName, 
                         const DEMParticlePArray& particles, 
                         int frame) const;
  void printParticlesXML(const std::string& folderName,
                         const std::string& fileName, 
                         int frame) const;
  void printParticlesXML(const std::string& folderName,
                         const std::string& fileName, 
                         const DEMParticlePArray& particles, 
                         int frame) const;
  void printBoundaryContacts() const; // print all boundary contact info
  void printMemParticle(
    const std::string& str) const; // print membrane particles
  void plotSpring(
    const std::string& str) const; // print springs in Tecplot format
  void writeBoundaryToFile() const;
  void writeBoundaryToFile(const OrientedBox& domain) const;
  void printBoundary() const; // print rigid boundaries info
  void writePatchGridToFile() const;
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

  // scale the dem with particle boundaries from deposited state until it
  // reaches steady state
  void scale_PtclBdry(
    std::size_t total_steps = 50000, // total_steps
    std::size_t snapNum = 100,       // number of snapNum
    std::size_t interval = 10,       // print interval
    REAL dimn = 0.05,                // dimension of particle-composed-boundary
    REAL rsize = 1.0,                // relative spatial domain size
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
    REAL rsize = 1.0,                // relative spatial domain size
    const std::string& iniptclfile =
      "flo_particle_end", // input file, initial particles
    const std::string& Particlefile =
      "dep_particle", // output file, resulted particles, including snapNum
    const std::string& contactfile =
      "dep_contact", // output file, resulted contacts, including snapNum
    const std::string& progressfile =
      "dep_progress", // output file, statistical info
    const std::string& debugfile = "dep_debug"); // output file, debug info

  // squeeze paticles inside a spatial domain by moving the boundaries
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

  REAL getPtclMinX(const DEMParticlePArray& d_patchParticles) const;
  REAL getPtclMaxX(const DEMParticlePArray& d_patchParticles) const;
  REAL getPtclMinY(const DEMParticlePArray& d_patchParticles) const;
  REAL getPtclMaxY(const DEMParticlePArray& d_patchParticles) const;
  REAL getPtclMinZ(const DEMParticlePArray& d_patchParticles) const;
  REAL getPtclMaxZ(const DEMParticlePArray& d_patchParticles) const;

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

  void constructBoundarySandPeriBonds(); // construct boundary bonds, sand
                                         // bonds, July 14, 2014

  void applyCoupledForces(); // check if peri-boundary bonds and peri-dem bonds
                             // are still alive, July 15, 2014

  // if so, apply coupled forces to peri-points and sand particles/boundaries
  void applyCoupledBoundary(); // this is used to test the coupled force model,
                               // October 10, 2014

  void findParticleInBox(const Box& spatialDomain,
                         const DEMParticlePArray& allParticle,
                         DEMParticlePArray& foundParticle);

  void dragForce();

  void applyParticleBC(double time, const Box& spatialDomain,
                       OrientedBox& modifiableDomain,
                       DEMParticlePArray& particles) {
    d_bc.applyParticleBC(time, spatialDomain, modifiableDomain, particles);
  }

  void applyPatchParticleBC(double time, OrientedBox& modifiableDomain) {
    applyParticleBC(time, d_spatialDomain, modifiableDomain, d_patchParticles);
  }

  void allowPatchDomainResize(Boundary::BoundaryID face) {
    d_allowPatchDomainResize[face] = true;
  }
  bool patchDomainResizeAllowed(Boundary::BoundaryID face) const {
    if (d_allowPatchDomainResize.find(face) != d_allowPatchDomainResize.end()) {
      return true;
    }
    return false;
  }

  void resizePatchBox();

  void switchParticleType(const DEMParticle::DEMParticleType& from,
                          const DEMParticle::DEMParticleType& to,
                          DEMParticlePArray& particles) {
    for (auto& particle : particles) {
      if (particle->getType() == from) {
        particle->setType(to);
      }
    }
  }

private:

  std::map<Boundary::BoundaryID, bool> d_allowPatchDomainResize;

  // The output writer pointer
  std::unique_ptr<Output> d_writer;

  // The boundary condition data
  DEMBoundaryConditions d_bc;

  // Particle gradation; broadcast among processes for once
  Gradation d_gradation; 
  REAL d_minParticleRadius;
  REAL d_maxParticleRadius;

  // all particles, only meaningful to root process
  DEMParticlePArray d_allDEMParticles;           
  DEMParticlePArray d_patchParticles; // particles per process

  // historical maximum numbering before trimming,
  // only meaningful to root process
  std::size_t d_trimHistoryNum; 

  ContactArray d_contacts; // contacts per process

  // tangential contact force and displacement per process
  ContactTangentArray d_contactTangents; 

  // estimated total contact number, only meaningful to root process
  std::size_t d_allContactNum; 

  MembraneParticlePArray d_membraneBoundaries; // membrane particle boundaries
  SpringUPArray d_springs;            // springs connecting membrane particles

  // whole spatial domain, broadcast among processes for once
  Box d_spatialDomain; 
  Box d_patchDomain;    // spatial domain per process
  Box cavity;       // cavity inside spatial domain

  // adaptive compute patchBox, broadcast among processes for once,
  // updated per process
  Box d_demPatchBox; 

  // rigid boundaries, broadcast among processes upon changed.
  BoundaryPArray d_boundaries; 

  // rigid boundaries with stats from all processes
  BoundaryPArray d_mergedBoundaries; 

  BoundaryPArray d_cavityBoundaries; // rigid cavity boundaries

  // particle-boundary contact tangential info
  std::map<std::size_t, BoundaryTangentArray> d_boundaryTangentMap; 

  // fluid property
  Fluid fluid;

  // average data
  REAL d_avgNormalForce; // only meaningful to root process
  REAL d_avgShearForce;  // only meaningful to root process
  REAL d_avgPenetration; // only meaningful to root process
  std::size_t d_overlapCount;

  // energy data
  REAL d_translationalEnergy; // only meaningful to root process
  REAL d_rotationalEnergy; // only meaningful to root process
  REAL d_kineticEnergy; // only meaningful to root process
  REAL d_gravitationalEnergy; // only meaningful to root process
  REAL d_mechanicalEnergy; // only meaningful to root process

  // time step
  REAL d_vibrationTimeStep;  // meaningful to all processes
  REAL d_impactTimeStep; // meaningful to all processes

  // MPI data
  boost::mpi::communicator boostWorld;
  static MPI_Comm s_mpiWorld;
  static MPI_Comm s_cartComm;
  static int s_mpiRank;
  static int s_mpiSize;
  static IntVec s_mpiProcs;
  static IntVec s_mpiCoords;
  int mpiTag;
  std::vector<std::size_t> bdryProcess;

  std::unique_ptr<Patch<DEMParticlePArray>> d_patchP;
  void createPatch(int iteration, const REAL& ghostWidth);
  void updatePatch(int iteration, const REAL& ghostWidth);

  ParticleIDHashMap d_sentParticles;  // sent particles per process
  DEMParticlePArray d_receivedParticles;  // received particles per process
  DEMParticlePArray d_mergedParticles; // merged particles per process

  // stream
  std::ofstream progressInf;
  std::ofstream balancedInf;

};

} // namespace dem

#endif
