#ifndef SMOOTHPARTICLE_HYDRODYNAMICS_H
#define SMOOTHPARTICLE_HYDRODYNAMICS_H

#include <Core/Parallel/Patch.h>
#include <InputOutput/Output.h>
#include <SmoothParticleHydro/SPHContainers.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <fstream>

#ifndef isProc0_macro
#define isProc0_macro (sph::SmoothParticleHydro::getMPIRank() == 0)
#define proc0cout                                                              \
  if (isProc0_macro)                                                           \
  std::cout
#define proc0cerr                                                              \
  if (isProc0_macro)                                                           \
  std::cerr
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

  void setCommunicator(const boost::mpi::communicator& boostWorldComm);
  void copyCommunicator(const boost::mpi::communicator& boostWorldComm,
                        const MPI_Comm& mpiWorldComm,
                        const MPI_Comm& mpiCartComm, int mpiRank, int mpiSize,
                        int mpiTag, const dem::IntVec& mpiProcs,
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

  // Accessor methods
  inline const SPHParticlePArray& getAllSPHParticleVec() const
  {
    return d_allSPHParticleVec;
  }
  inline const SPHParticlePArray& getSPHParticleVec() const
  {
    return d_sphParticleVec;
  }

  inline const SPHParticlePArray& getMergedSPHParticleVec() const
  {
    return d_mergeSPHParticleVec;
  }

  inline void setAllSPHParticleVec(const SPHParticlePArray& particles) 
  {
    d_allSPHParticleVec = particles;
  }

  void clearAllSPHParticleVec() { d_allSPHParticleVec.clear(); }

  // Scatter the sph particles
  void scatterSPHParticle(const dem::Box& allContainer,
                          const REAL& ghostWidth,
                          REAL& bufferLength);

  void findSPHParticleInBox(const dem::Box& container,
                            const SPHParticlePArray& allParticles,
                            SPHParticlePArray& foundParticles);

  void commuSPHParticle(int iteration,
                         const double& ghostWidth);
  /*


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

  // only for master cpu
  SPHParticlePArray d_allSPHParticleVec;

  // Temporary
  SPHParticlePArray d_patchSPHParticleVec;

  // coordinates of all the particles in this cpu, local
  // sph-points in current cpu
  SPHParticlePArray d_sphParticleVec;

  // Particles in patch + ghpst
  SPHParticlePArray d_mergeSPHParticleVec;


};

} // namespace sph

#endif
