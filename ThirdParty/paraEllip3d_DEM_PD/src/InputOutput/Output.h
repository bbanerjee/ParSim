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
#include <DiscreteElements/DEMContainers.h>
#include <Peridynamics/PeriContainers.h>
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

  virtual void write(int frame);

  virtual void setDomain(const Box* domain) {};
  virtual void setGrid(const Box* grid) {};
  virtual void setParticles(const ParticlePArray* particles) {};
  virtual void setParticles(const pd::PeriParticlePArray* particles) {};

  virtual void writeDomain(const Box* domain) {};
  virtual void writeGrid(const Box* grid) {};
  virtual void writeParticles(const ParticlePArray* particles, int frame) {};
  virtual void writeParticles(const pd::PeriParticlePArray* particles, int frame) {};
  virtual void writeSieves(const Gradation* gradation) {};

  void createFileNames();
  void updateFileNames(const int& iteration, const std::string& extension);
  void updateFileNames(const int& iteration);
  std::string getDomainFileName() const { return d_domainFileName; }
  std::string getBoundaryFileName() const { return d_boundaryFileName; }
  std::string getGridFileName() const { return d_gridFileName; }
  std::string getParticleFileName() const { return d_particleFileName; }
  std::string getPeriParticleFileName() const { return d_periParticleFileName; }
  std::string getBdryContactFileName() const { return d_bdryContactFileName; }
  void setParticleFileName(const std::string& name) { d_particleFileName = d_outputFolderName + "/" + name; }
  void setPeriParticleFileName(const std::string& name) { d_periParticleFileName = d_outputFolderName + "/" + name; }
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

  std::string d_domainFileName;
  std::string d_boundaryFileName;
  std::string d_gridFileName;
  std::string d_particleFileName;
  std::string d_periParticleFileName;
  std::string d_bdryContactFileName;

private:
  Output(const Output&) = delete;
}; // end class

} // end namespace

#endif
