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

#include <InputOutput/Output.h>

#include <exception>
#include <fstream>
#include <unistd.h>
#include <experimental/filesystem>

using namespace dem;

Output::Output(const std::string& fileName, int iterInterval)
{
  d_outputFolderName = fileName;
  d_outputIteration = iterInterval;

  d_cartComm = nullptr;
  d_mpiProcX = 1;
  d_mpiProcY = 1;
  d_mpiProcZ = 1;
}

Output::~Output() = default;

void
Output::clone(const Output& output)
{
  d_outputFolderName = output.d_outputFolderName;
  d_outputIteration = output.d_outputIteration;
  d_mpiProcX = output.d_mpiProcX;
  d_mpiProcY = output.d_mpiProcY;
  d_mpiProcZ = output.d_mpiProcZ;
}

// Create individual file names
void
Output::createFilenames()
{
  std::string folderName = d_outputFolderName;

  std::ostringstream domainOutputFile;
  domainOutputFile << folderName << "/domain_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_domainFilename = domainOutputFile.str();

  std::ostringstream boundaryOutputFile;
  boundaryOutputFile << folderName << "/boundary_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_boundaryFilename = boundaryOutputFile.str();

  std::ostringstream patchOutputFile;
  patchOutputFile << folderName << "/patchGrid_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_patchBoxFilename = patchOutputFile.str();

  std::ostringstream particleOutputFile;
  particleOutputFile << folderName << "/particle_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_particleFilename = particleOutputFile.str();

  std::ostringstream periParticleOutputFile;
  periParticleOutputFile << folderName << "/peri_particle_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_periParticleFilename = periParticleOutputFile.str();

  std::ostringstream sphParticleOutputFile;
  sphParticleOutputFile << folderName << "/sph_particle_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_sphParticleFilename = sphParticleOutputFile.str();

  std::ostringstream bdryContactOutputFile;
  bdryContactOutputFile << folderName << "/bdrycontact_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_bdryContactFilename = bdryContactOutputFile.str();

  //std::cout << "createFilenames: Domain file name = " << d_domainFilename << "\n";
  //std::cout << "createFilenames: Grid file name = " << d_patchFilename << "\n";
  //std::cout << "createFilenames: DEMParticle file name = " << d_particleFilename << "\n";
}

void
Output::updateFilenames(const int& iteration, const std::string& extension) {
  d_outputIteration = iteration;
  createFilenames();
  d_domainFilename += extension;
  d_boundaryFilename += extension;
  d_patchBoxFilename += extension;
  d_particleFilename += extension;
  d_periParticleFilename += extension;
  d_sphParticleFilename += extension;
  d_bdryContactFilename += extension;
}

void
Output::updateFilenames(const int& iteration) {
  d_outputIteration = iteration;
  createFilenames();
}

void
Output::write(int frame)
{
  //std::cout << "Please call the writer routine from the correct derived class.\n";
}

namespace dem {

std::ostream&
operator<<(std::ostream& out, const Output& output)
{
  out.setf(std::ios::floatfield);
  out.precision(6);
  out << "Output file = " << output.d_outputFolderName << std::endl;
  out << "  Output iteration interval = " << output.d_outputIteration
      << std::endl;
  return out;
}

}
