#ifndef SMOOTHPARTICLE_HYDRODYNAMICS_H
#define SMOOTHPARTICLE_HYDRODYNAMICS_H

#include <SmoothParticleHydro/SPHContainers.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <Core/Parallel/Patch.h>
#include <InputOutput/Output.h>
#include <fstream>

#ifndef isProc0_macro
  #define isProc0_macro ( sph::SmoothParticleHydro::getMPIRank() == 0 )
  #define proc0cout if( isProc0_macro ) std::cout
  #define proc0cerr if( isProc0_macro ) std::cerr
#endif

namespace sph {

using SPHPatch = dem::Patch<sph::SPHParticlePArray>;
using Output = dem::Output;

class SmoothParticleHydro
{
public:

  SmoothParticleHydro();
  ~SmoothParticleHydro();

  // MPI methods
  inline int getMPIRank() const { return d_mpiRank; }
  inline boost::mpi::communicator getMPIWorld() const { return d_boostWorld; }

  // Accessor methods
  inline const SPHParticlePArray& getAllSPHParticleVec() const
  {
    return d_allSPHParticleVec;
  }
  inline const SPHParticlePArray& getSPHParticleVec() const
  {
    return d_sphParticleVec;
  }
  inline const SPHParticlePArray& getRecvSPHParticleVec() const
  {
    return d_recvSPHParticleVec;
  }

  void clearAllSPHParticleVec() { d_allSPHParticleVec.clear(); }

  void setCommunicator(const boost::mpi::communicator& boostWorldComm,
                       const MPI_Comm& mpiWorldComm, 
                       const MPI_Comm& mpiCartComm,
                       int mpiRank, int mpiSize, int mpiTag,
                       const dem::IntVec& mpiProcs, 
                       const dem::IntVec& mpiCoords)
  {
    d_boostWorld = boostWorldComm;
    d_mpiWorld = mpiWorldComm;
    d_cartComm = mpiCartComm;
    d_mpiRank = mpiRank; 
    d_mpiSize = mpiSize; 
    d_mpiTag = mpiTag;
    d_mpiProcs = mpiProcs;
    d_mpiCoords = mpiCoords;
  }

  inline void setGrid(dem::Box cont) { d_sphGrid = cont; }

  template <int dim>
  void generateSPHParticle(const dem::Box& allContainer,
                           dem::DEMParticlePArray& allDEMParticles);
                          
  template <int dim>
  REAL computeMass(const double& density, const double& length, 
                   const std::size_t& numPts) const;

  template <int dim>
  void createCoords(const dem::Vec& vmin, const dem::Vec& vmax, 
                    const REAL& spaceInterval, 
                    const int& numLayers, 
                    std::vector<REAL>& xCoords, 
                    std::vector<REAL>& yCoords, 
                    std::vector<REAL>& zCoords) const;

  template <int dim>
  void createParticleArray(const REAL& mass,
                           const REAL& density,
                           const std::vector<REAL>& xCoords,
                           const std::vector<REAL>& yCoords,
                           const std::vector<REAL>& zCoords);

  void removeRedundantSPHParticles();

  /*
  // Scatter the sphdynamics particles
  void scatterSPHParticle(const dem::Box& allContainer);

  void findSPHParticleInBox(const dem::Box& container, 
                             const SPHParticlePArray& allParticles,
                             SPHParticlePArray& foundParticles);

  void commuSPHParticle(int iteration,
                         const double& maxDEMParticleRadius);

  void removeSPHParticleOutBox(const dem::Box& container,
                                SPHParticlePArray& d_sphParticleVec);

  bool isBdryProcess();

  void releaseRecvSPHParticle();
  void releaseGatheredSPHParticle();

  void migrateSPHParticle(int iteration);
  void gatherSPHParticle();

  void updateSPHGrid(const SPHParticlePArray& particles);

  void openSPHProgress(std::ofstream& ofs, const std::string& str);
  void printSPHProgress(std::ofstream& ofs, const int iframe) const;
  void printSPHProgressHalf(std::ofstream& ofs, const int iframe) const;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  inline int getNumSPHParticle() const
  {
    return nSPHParticle;
  } // getnSPHParticle - returns number of sph-particles

  // construct Matrix members in
  // d_recvSPHParticleVec, construction here
  // since currently
  // the pointer array in Matrix cannot be transfered well between cpus
  void constructRecvSPHMatrix(); 

  // delete those sph-points that are inside dem particles or vice versa
  // depending on the removeSPHParticles flag (if false DEM particles
  // are removed)
  void removeOverlappingParticles(dem::DEMParticlePArray& particles,
                                  bool removeSPHParticles = true);     

  void removeInsideSPHParticles(const dem::DEMParticlePArray& particles);
  void removeInsideDEMParticles(dem::DEMParticlePArray& particles) const;     

  //-------------------------------------------------------------
  // Output
  //-------------------------------------------------------------
  void createOutputWriter(const std::string& outputFolder, const int& iter); 

  void updateFileNames(const int& iter, const std::string& extension) {
    d_writer->updateFileNames(iter, extension);
  }
  void updateFileNames(const int& iter) {
    d_writer->updateFileNames(iter);
  }
  std::string getSPHParticleFileName() const {
    return d_writer->getSPHParticleFileName();
  }

  // print all particles
  void writeParticlesToFile(int frame) const; 
  // print a subset of particles
  void writeParticlesToFile(SPHParticlePArray& particleVec, int frame) const; 

  void printSPHDomain(const std::string&) const;
  void printRecvSPHDomain(const std::string&) const;
  void printSPHParticle(const std::string& str) const;
  */

private:

  // The output writer pointer
  std::unique_ptr<Output> d_writer;

  // MPI data
  boost::mpi::communicator d_boostWorld;
  MPI_Comm d_mpiWorld, d_cartComm;
  int d_mpiRank, d_mpiSize, d_mpiTag;

  dem::IntVec d_mpiProcs;
  dem::IntVec d_mpiCoords;

  std::unique_ptr<SPHPatch> d_patchP;
  void createPatch(int iteration, const REAL& ghostWidth);
  void updatePatch(int iteration, const REAL& ghostWidth);

  dem::Box d_sphGrid;

  // stream
  std::ofstream sphProgInf;
  std::ofstream sphProgInfHalf;

  // number of SPHParticles in the domain, only for master cpu
  int nSPHParticle; 

  // dimension of the problem, only for master cpu
  int ndim; 

  // only for master cpu
  SPHParticlePArray d_allSPHParticleVec; 

  // coordinates of all the particles in this cpu, local 
  // sph-points in current cpu
  SPHParticlePArray d_sphParticleVec; 
  SPHParticlePArray ghostSPHParticleVec; 

  SPHParticlePArray d_recvSPHParticleVec;
  SPHParticlePArray mergeSPHParticleVec;
};

} // namespace sph

#endif
