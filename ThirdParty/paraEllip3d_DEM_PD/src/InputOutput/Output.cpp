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
Output::createFileNames()
{
  std::string folderName = d_outputFolderName;

  std::ostringstream domainOutputFile;
  domainOutputFile << folderName << "/domain_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_domainFileName = domainOutputFile.str();

  std::ostringstream boundaryOutputFile;
  boundaryOutputFile << folderName << "/boundary_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_boundaryFileName = boundaryOutputFile.str();

  std::ostringstream gridOutputFile;
  gridOutputFile << folderName << "/grid_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_gridFileName = gridOutputFile.str();

  std::ostringstream particleOutputFile;
  particleOutputFile << folderName << "/particle_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_particleFileName = particleOutputFile.str();

  std::ostringstream bdryContactOutputFile;
  bdryContactOutputFile << folderName << "/bdrycontact_"
    << std::setfill('0') << std::setw(5) << d_outputIteration;
  d_bdryContactFileName = bdryContactOutputFile.str();

  std::cout << "createFileNames: Domain file name = " << d_domainFileName << "\n";
  std::cout << "createFileNames: Grid file name = " << d_gridFileName << "\n";
  std::cout << "createFileNames: Particle file name = " << d_particleFileName << "\n";
}

void
Output::updateFileNames(const int& iteration, const std::string& extension) {
  d_outputIteration = iteration;
  createFileNames();
  d_domainFileName += extension;
  d_boundaryFileName += extension;
  d_gridFileName += extension;
  d_particleFileName += extension;
  d_bdryContactFileName += extension;
}

void
Output::updateFileNames(const int& iteration) {
  d_outputIteration = iteration;
  createFileNames();
}

void
Output::write()
{
  std::cout << "Please call the writer routine from the correct derived class.\n";
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
