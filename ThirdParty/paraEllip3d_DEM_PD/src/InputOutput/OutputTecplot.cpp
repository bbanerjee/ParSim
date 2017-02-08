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

#include <DiscreteElements/Particle.h>
#include <InputOutput/OutputTecplot.h>

#include <fstream>
#include <regex>
#include <unistd.h>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace dem;

OutputTecplot::OutputTecplot(const std::string& fileName, int iterInterval)
  : Output(fileName, iterInterval)
{
  d_domain = nullptr;
  d_grid = nullptr;
  d_particles = nullptr;
  d_cartComm = 0;
}

OutputTecplot::~OutputTecplot()
{
}

void
OutputTecplot::write()
{

  // The domain and the grid have to be set before a write is
  // completed.
  if (!d_domain || !d_grid || !d_particles) {
    std::cout << "**ERROR** Domain and/or Grid and/or Particles have not been "
                 "set.  Nothing "
                 "will be written\n";
    return;
  }

  // Get the file names
  std::ostringstream domainOutputFile, gridOutputFile, particleOutputFile;
  getFileNames(domainOutputFile, gridOutputFile, particleOutputFile);

  // Write files for the domain extents at each timestep
  writeDomain(d_domain, domainOutputFile.str());

  // Write files for the grid representing each processor at each timestep
  writeGrid(d_grid, gridOutputFile.str());

  // Write files for the particle list each timestep
  writeParticles(d_particles, gridOutputFile.str());

  // Increment the output file count
  incrementOutputFileCount();
}

void
OutputTecplot::writeDomain(const Box* domain, const std::string& fileName)
{
  std::ofstream ofs(fileName);
  if (!ofs) {
    debugInf << "stream error: plotBoundary" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1, y1, z1, x2, y2, z2;
  x1 = domain->getMinCorner().getX();
  y1 = domain->getMinCorner().getY();
  z1 = domain->getMinCorner().getZ();
  x2 = domain->getMaxCorner().getX();
  y2 = domain->getMaxCorner().getY();
  z2 = domain->getMaxCorner().getZ();

  ofs << "ZONE N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z1
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y1 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x2 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y2 << std::setw(OWID) << z2
      << std::endl;
  ofs << std::setw(OWID) << x1 << std::setw(OWID) << y1 << std::setw(OWID) << z2
      << std::endl;
  ofs << "1 2 3 4 5 6 7 8" << std::endl;

  ofs.close();
}

void
OutputTecplot::writeGrid(const Box* grid, const std::string& fileName)
{
  std::ofstream ofs(fileName);
  if (!ofs) {
    debugInf << "stream error: plotGrid" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  // Get the number of processes
  int mpiSize = 0;
  MPI_Comm_size(d_cartComm, &mpiSize);

  Vec v1 = grid->getMinCorner();
  Vec v2 = grid->getMaxCorner();
  Vec vspan = v2 - v1;

  ofs << "ZONE N=" << (d_mpiProcX + 1) * (d_mpiProcY + 1) * (d_mpiProcZ + 1)
      << ", E=" << d_mpiProcX * d_mpiProcY * d_mpiProcZ
      << ", DATAPACKING=POINT, ZONETYPE=FEBRICK" << std::endl;

  std::vector<Vec> coords((d_mpiProcX + 1) * (d_mpiProcY + 1) *
                          (d_mpiProcZ + 1));
  std::size_t index = 0;
  for (auto i = 0; i < d_mpiProcX + 1; ++i)
    for (auto j = 0; j < d_mpiProcY + 1; ++j)
      for (auto k = 0; k < d_mpiProcZ + 1; ++k)
        coords[index++] = Vec(v1.getX() + vspan.getX() / d_mpiProcX * i,
                              v1.getY() + vspan.getY() / d_mpiProcY * j,
                              v1.getZ() + vspan.getZ() / d_mpiProcZ * k);

  for (auto i = 0; i < (d_mpiProcX + 1) * (d_mpiProcY + 1) * (d_mpiProcZ + 1);
       ++i)
    ofs << std::setw(OWID) << coords[i].getX() << std::setw(OWID)
        << coords[i].getY() << std::setw(OWID) << coords[i].getZ() << std::endl;

  for (int iRank = 0; iRank < mpiSize; ++iRank) {
    int coords[3];
    MPI_Cart_coords(d_cartComm, iRank, 3, coords);

    int id4 = 1 + coords[0] * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              coords[1] * (d_mpiProcZ + 1) + coords[2];
    int id1 = 1 + (coords[0] + 1) * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              coords[1] * (d_mpiProcZ + 1) + coords[2];
    int id3 = 1 + coords[0] * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              (coords[1] + 1) * (d_mpiProcZ + 1) + coords[2];
    int id2 = 1 + (coords[0] + 1) * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              (coords[1] + 1) * (d_mpiProcZ + 1) + coords[2];

    int id8 = 1 + coords[0] * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              coords[1] * (d_mpiProcZ + 1) + (coords[2] + 1);
    int id5 = 1 + (coords[0] + 1) * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              coords[1] * (d_mpiProcZ + 1) + (coords[2] + 1);
    int id7 = 1 + coords[0] * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              (coords[1] + 1) * (d_mpiProcZ + 1) + (coords[2] + 1);
    int id6 = 1 + (coords[0] + 1) * (d_mpiProcZ + 1) * (d_mpiProcY + 1) +
              (coords[1] + 1) * (d_mpiProcZ + 1) + (coords[2] + 1);

    ofs << std::setw(8) << id1 << std::setw(8) << id2 << std::setw(8) << id3
        << std::setw(8) << id4 << std::setw(8) << id5 << std::setw(8) << id6
        << std::setw(8) << id7 << std::setw(8) << id8 << std::endl;
  }

  ofs.close();
}

void
OutputTecplot::writeParticles(const ParticlePArray* particles,
                              const std::string& fileName)
{
  std::ofstream ofs(fileName);
  if (!ofs) {
    debugInf << "stream error: printParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  ofs << std::setw(OWID) << particles->size() << std::endl;
  ofs << std::setw(OWID) << "id" << std::setw(OWID) << "type" << std::setw(OWID)
      << "radius_a" << std::setw(OWID) << "radius_b" << std::setw(OWID)
      << "radius_c" << std::setw(OWID) << "position_x" << std::setw(OWID)
      << "position_y" << std::setw(OWID) << "position_z" << std::setw(OWID)
      << "axle_a_x" << std::setw(OWID) << "axle_a_y" << std::setw(OWID)
      << "axle_a_z" << std::setw(OWID) << "axle_b_x" << std::setw(OWID)
      << "axle_b_y" << std::setw(OWID) << "axle_b_z" << std::setw(OWID)
      << "axle_c_x" << std::setw(OWID) << "axle_c_y" << std::setw(OWID)
      << "axle_c_z" << std::setw(OWID) << "velocity_x" << std::setw(OWID)
      << "velocity_y" << std::setw(OWID) << "velocity_z" << std::setw(OWID)
      << "omga_x" << std::setw(OWID) << "omga_y" << std::setw(OWID) << "omga_z"
      << std::setw(OWID) << "force_x" << std::setw(OWID) << "force_y"
      << std::setw(OWID) << "force_z" << std::setw(OWID) << "moment_x"
      << std::setw(OWID) << "moment_y" << std::setw(OWID) << "moment_z"
      << std::endl;

  Vec vObj;
  for (const auto& part : *particles) {
    ofs << std::setw(OWID) << part->getId() << std::setw(OWID)
        << part->getType() << std::setw(OWID) << part->getA() << std::setw(OWID)
        << part->getB() << std::setw(OWID) << part->getC();

    vObj = part->getCurrPos();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getCurrDirecA();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getCurrDirecB();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getCurrDirecC();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getCurrVeloc();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getCurrOmga();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getForce();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ();

    vObj = part->getMoment();
    ofs << std::setw(OWID) << vObj.getX() << std::setw(OWID) << vObj.getY()
        << std::setw(OWID) << vObj.getZ() << std::endl;
  }
  ofs.close();
}

void
OutputTecplot::writeSieves(const Gradation* gradation,
                           const std::string& fileName)
{

  std::ofstream ofs(fileName, std::ofstream::app);
  if (!ofs) {
    debugInf << "stream error: printParticle" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  std::size_t sieveNum = gradation->getSieveNum();
  std::vector<REAL> percent = gradation->getPercent();
  std::vector<REAL> size = gradation->getSize();
  ofs << std::endl << std::setw(OWID) << sieveNum << std::endl;
  for (std::size_t i = 0; i < sieveNum; ++i)
    ofs << std::setw(OWID) << percent[i] << std::setw(OWID) << size[i]
        << std::endl;
  ofs << std::endl
      << std::setw(OWID) << gradation->getPtclRatioBA() << std::setw(OWID)
      << gradation->getPtclRatioCA() << std::endl;

  ofs.close();
}

// Get individual file names
void
OutputTecplot::getFileNames(std::ostringstream& domainFileName,
                            std::ostringstream& nodeFileName,
                            std::ostringstream& particleFileName)
{
  std::string name_without_ext = outputFile();
  unsigned int lastIndex = name_without_ext.find_last_of(".");
  name_without_ext = name_without_ext.substr(0, lastIndex);

  if (outputFileCount() == 0) {
    int dircount = 0;
    d_output_dir << "./" << name_without_ext;
    d_output_dir << ".tecplot." << std::setfill('0') << std::setw(3)
                 << dircount;

#if defined(_WIN32)
    _mkdir((d_output_dir.str()).c_str());
#else
    while (opendir((d_output_dir.str()).c_str())) {
      ++dircount;
      d_output_dir.clear();
      d_output_dir.str("");
      d_output_dir << "./" << name_without_ext;
      d_output_dir << ".tecplot." << std::setfill('0') << std::setw(3)
                   << dircount;
    }
    mkdir((d_output_dir.str()).c_str(), 0777);
#endif
  }

  domainFileName << d_output_dir.str() << "/" << name_without_ext << "_d_"
                 << std::setfill('0') << std::setw(5) << outputFileCount();
  nodeFileName << d_output_dir.str() << "/" << name_without_ext << "_n_"
               << std::setfill('0') << std::setw(5) << outputFileCount();
  particleFileName << d_output_dir.str() << "/" << name_without_ext << "_p_"
                   << std::setfill('0') << std::setw(5) << outputFileCount();
}
