/*
 * The MIT License
 *
 * Copyright (c) 2017- Parresia Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef __DEM_OUTPUT_H__
#define __DEM_OUTPUT_H__

#include <Core/Geometry/Box.h>
#include <Core/Geometry/OrientedBox.h>
#include <DiscreteElements/DEMContainers.h>
#include <Peridynamics/PeriContainers.h>
#include <SmoothParticleHydro/SPHContainers.h>
#include <DiscreteElements/Gradation.h>

#include <mpi.h>

#include <iostream>
#include <string>

namespace dem {

class Output
{
public:
  friend std::ostream& operator<<(std::ostream& out, const dem::Output& output);

public:
  Output() = delete;
  Output(const std::string& folderName, int iterInterval);
  virtual ~Output();

  void clone(const Output& output);

  virtual void write(int frame, REAL time);

  virtual void setDomain(const Box* domain) {};
  virtual void setPatchBox(const Box* patchBox) {};
  virtual void setParticles(const DEMParticlePArray* particles) {};
  virtual void setParticles(const pd::PeriParticlePArray* particles) {};
  virtual void setParticles(const sph::SPHParticlePArray* particles) {};

  virtual void writeDomain(const Box* domain, REAL time) {};
  virtual void writeDomain(const OrientedBox& domain, REAL time) {};
  virtual void writePatchBoxGrid(const Box* patchBox, REAL time) {};
  virtual void writeParticles(const DEMParticlePArray* particles, 
                              int frame, REAL time) {};
  virtual void writeParticles(const pd::PeriParticlePArray* particles, 
                              int frame, REAL time) {};
  virtual void writeParticles(const sph::SPHParticlePArray* particles, 
                              int frame, REAL time) {};
  virtual void writeSieves(const Gradation* gradation) {};

  void createFilenames();
  void updateFilenames(const int& iteration, const std::string& extension);
  void updateFilenames(const int& iteration);
  std::string getDomainFilename() const { return d_domainFilename; }
  std::string getBoundaryFilename() const { return d_boundaryFilename; }
  std::string getPatchBoxFilename() const { return d_patchBoxFilename; }
  std::string getParticleFilename() const { return d_particleFilename; }
  std::string getPeriParticleFilename() const { return d_periParticleFilename; }
  std::string getSPHParticleFilename() const { return d_sphParticleFilename; }
  std::string getBdryContactFilename() const { return d_bdryContactFilename; }
  std::string getContactFilename() const { return d_contactFilename; }
  void setParticleFilename(const std::string& name) { d_particleFilename = d_outputFolderName + "/" + name; }
  void setPeriParticleFilename(const std::string& name) { d_periParticleFilename = d_outputFolderName + "/" + name; }
  void setSPHParticleFilename(const std::string& name) { d_sphParticleFilename = d_outputFolderName + "/" + name; }
  inline std::string outputFolder() const { return d_outputFolderName; }
  inline int outputIteratonInterval() const { return d_outputIteration; }

  // Set the processor distribution
  void setMPIComm(const MPI_Comm& cartComm) { d_cartComm = cartComm; }

  void setMPIProc(size_t x, size_t y, size_t z)
  {
    d_mpiProcX = x;
    d_mpiProcY = y;
    d_mpiProcZ = z;
  }

protected:

  // Processor distribution
  MPI_Comm d_cartComm;
  std::size_t d_mpiProcX;
  std::size_t d_mpiProcY;
  std::size_t d_mpiProcZ;

  //  Output file name
  std::string d_outputFolderName;
  int d_outputIteration;

  std::string d_domainFilename;
  std::string d_orientedDomainFilename;
  std::string d_boundaryFilename;
  std::string d_patchBoxFilename;
  std::string d_particleFilename;
  std::string d_periParticleFilename;
  std::string d_sphParticleFilename;
  std::string d_bdryContactFilename;
  std::string d_contactFilename;

private:
  Output(const Output&) = delete;
}; // end class

} // end namespace

#endif
