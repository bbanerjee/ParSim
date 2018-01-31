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

  inline void setPatchBox(dem::Box cont) { d_sphPatchBox = cont; }

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
  void scatterSPHParticle(const dem::Box& spatialDomain,
                          const REAL& ghostWidth,
                          REAL& bufferLength);

  void commuSPHParticle(int iteration,
                         const double& ghostWidth);

  void migrateSPHParticle(int iteration);

  void gatherSPHParticle();

  template <int dim>
  void updateParticleInteractions(const dem::Box& spatialDomain,
                                  const REAL& bufferWidth,
                                  const REAL& ghostWidth,
                                  const REAL& kernelSize,
                                  const REAL& smoothLength);

  template <int dim>
  void initializeSPHVelocity(const REAL& delT);

  void updateSPHLeapFrogPositionDensity(const REAL& delT);

  void updateSPHLeapFrogVelocity(const REAL& delT);

  //-------------------------------------------------------------
  // Output
  //-------------------------------------------------------------
  void createOutputWriter(const std::string& outputFolder, const int& iter); 

  void updateFileNames(const int& iter, const std::string& extension) {
    d_writer->updateFilenames(iter, extension);
  }
  void updateFileNames(const int& iter) {
    d_writer->updateFilenames(iter);
  }
  std::string getSPHParticleFileName() const {
    return d_writer->getSPHParticleFilename();
  }
  void writeParticlesToFile(int frame, REAL time) const; 
  void printSPHParticle(const char* str) const;

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

  dem::Box d_sphPatchBox;
  SPHPatchGridIndex     d_sphPatchGridIndex;
  SPHPatchGridParticleP d_sphPatchGrid;
  dem::IntVec d_numGridCells;

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

  void findSPHParticleInBox(const dem::Box& domain,
                            const SPHParticlePArray& allParticles,
                            SPHParticlePArray& foundParticles);

  void initializeDensityRateAndAcceleration();

  template <int dim>
  void assignParticlesToPatchGrid(const dem::Box& domain,
                                  const REAL& bufferWidth,
                                  const REAL& ghostWidth,
                                  const REAL& kernelSize);

  template <int dim>
  int getCellIndex(const dem::Vec& cellMinCorner,
                   const REAL& cellWidth,
                   const dem::IntVec& numGridCells,
                   const dem::Vec& pointPosition);

  template <int dim>
  int getIndex(const dem::IntVec& cell, const int& nx,
               const int& ny, const int& nz) const;

  template <int dim>
  std::vector<int> getAdjacentCellIndices(const dem::IntVec& cell,
                                          const dem::IntVec& numGridCells) const;

  template <int dim>
  SPHParticlePArray getParticlesInAdjacentCells(const int& cellIndex,
                                                const dem::IntVec& numGridCells) const;

  template <int dim>
  void computeParticleInteractions(SPHParticleP& sph_part_a,
                                   SPHParticleP& sph_part_b,
                                   const REAL& kernelSize,
                                   const REAL& smoothLength,
                                   const dem::Box& spatialDomain) const;

};

} // namespace sph

#endif
