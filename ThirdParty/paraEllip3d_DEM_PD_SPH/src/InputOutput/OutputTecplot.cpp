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

#include <DiscreteElements/DEMParticle.h>
#include <Peridynamics/PeriParticle.h>
#include <SmoothParticleHydro/SPHParticle.h>
#include <InputOutput/OutputTecplot.h>

#include <fstream>
#include <regex>
#include <unistd.h>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace dem;

using pd::PeriParticlePArray;
using sph::SPHParticlePArray;

template <typename TArray>
OutputTecplot<TArray>::OutputTecplot(const std::string& folderName, int iterInterval)
  : Output(folderName, iterInterval)
{
  d_domain = nullptr;
  d_patchBox = nullptr;
  d_particles = nullptr;
  d_cartComm = nullptr;

  // Create the basic file names (path + folder + name)
  createFilenames();

}

template <typename TArray>
OutputTecplot<TArray>::~OutputTecplot() = default;

template <typename TArray>
void
OutputTecplot<TArray>::write(int frame, REAL time)
{

  // The domain and the patchGrid have to be set before a write is
  // completed.
  if (!d_domain || !d_patchBox || !d_particles) {
    std::cerr << "**ERROR** Domain and/or Patch Grid and/or Particles have not been "
                 "set.  Nothing "
                 "will be written\n";
    return;
  }

  // Write files for the domain extents at each timestep
  writeDomain(d_domain, time);

  // Write files for the patchGrid representing each processor at each timestep
  writePatchBoxGrid(d_patchBox, time);

  // Write files for the particle list each timestep
  writeParticles(d_particles, frame, time);
}

template <typename TArray>
void
OutputTecplot<TArray>::writeDomain(const Box* domain, REAL time)
{
  // Get the filename
  std::string filename(d_domainFilename);
  filename.append(".dat");

  std::ofstream ofs(filename);
  if (!ofs) {
    debugInf << "stream error: writeBoundaryToFile" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  REAL x1, y1, z1, x2, y2, z2;
  x1 = domain->minCorner().x();
  y1 = domain->minCorner().y();
  z1 = domain->minCorner().z();
  x2 = domain->maxCorner().x();
  y2 = domain->maxCorner().y();
  z2 = domain->maxCorner().z();

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

template <typename TArray>
void
OutputTecplot<TArray>::writePatchBoxGrid(const Box* patchBox, REAL time)
{
  // Get the filename
  std::string filename(d_patchBoxFilename);
  filename.append(".dat");

  std::ofstream ofs(filename);
  if (!ofs) {
    debugInf << "stream error: writePatchGridToFile" << std::endl;
    exit(-1);
  }
  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);

  // Get the number of processes
  int mpiSize = 0;
  MPI_Comm_size(d_cartComm, &mpiSize);

  Vec v1 = patchBox->minCorner();
  Vec v2 = patchBox->maxCorner();
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
        coords[index++] = Vec(v1.x() + vspan.x() / d_mpiProcX * i,
                              v1.y() + vspan.y() / d_mpiProcY * j,
                              v1.z() + vspan.z() / d_mpiProcZ * k);

  for (auto i = 0; i < (d_mpiProcX + 1) * (d_mpiProcY + 1) * (d_mpiProcZ + 1);
       ++i)
    ofs << std::setw(OWID) << coords[i].x() << std::setw(OWID) << coords[i].y()
        << std::setw(OWID) << coords[i].z() << std::endl;

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

template <typename TArray>
void
OutputTecplot<TArray>::writeParticles(const TArray* particles, int frame,
                                      REAL time) {
  std::cout << "**ERROR** Noting to do here. The array of particles is"
            << " not of the correct type\n";
}

template <>
void
OutputTecplot<DEMParticlePArray>::writeParticles(const DEMParticlePArray* particles, 
                                                 int frame, REAL time)
{
  // Get the filename
  std::string filename(d_particleFilename);
  filename.append(".dat");

  //std::cout << "filename = " << filename << "\n";
  std::ofstream ofs(filename);
  if (!ofs) {
    std::cout << "Could not open" << std::endl;
    debugInf << "stream error: printParticlesCSV" << std::endl;
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
        << static_cast<int>(part->getType()) << std::setw(OWID) << part->radiusA() << std::setw(OWID)
        << part->radiusB() << std::setw(OWID) << part->radiusC();

    vObj = part->currentPosition();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentAnglesAxisA();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentAnglesAxisB();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentAnglesAxisC();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentVelocity();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->currentAngularVelocity();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->force();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z();

    vObj = part->moment();
    ofs << std::setw(OWID) << vObj.x() << std::setw(OWID) << vObj.y()
        << std::setw(OWID) << vObj.z() << std::endl;
  }
  ofs.close();
}

template <>
void
OutputTecplot<PeriParticlePArray>::writeParticles(const PeriParticlePArray* particles, 
                                                  int frame, REAL time)
{
  // Get the filename
  std::string filename(d_periParticleFilename);
  filename.append(".dat");

  std::ofstream ofs(filename);
  if (!ofs) {
    debugInf << "stream error: printParticlesCSV" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  if (frame == 0) {
    ofs << "Title = \"DEMParticle Information\"" << std::endl;
    ofs << "VARIABLES = \"X\", \"Y\",\"Z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" "
           "\"Vz\" \"KE\" \"P\" \"Mises\""
        << std::endl;
  }
  ofs << "ZONE T =\" " << frame << "-th Load Step\" " << std::endl;

  // Output the coordinates and the array information
  REAL pressure, vonMisesStress;
  Matrix sigma;
  for (const auto& pt : *particles) {
    sigma = pt->getSigma();
    pressure = sigma(1, 1) + sigma(2, 2) + sigma(3, 3);
    vonMisesStress =
      sqrt(0.5 * ((sigma(1, 1) - sigma(2, 2)) * (sigma(1, 1) - sigma(2, 2)) +
                  (sigma(2, 2) - sigma(3, 3)) * (sigma(2, 2) - sigma(3, 3)) +
                  (sigma(1, 1) - sigma(3, 3)) * (sigma(1, 1) - sigma(3, 3))) +
           3 * (sigma(1, 2) * sigma(1, 2) + sigma(2, 3) * sigma(2, 3) +
                sigma(3, 1) * sigma(3, 1)));
    ofs << std::setw(20)
        << pt->getInitPosition().x() + pt->getDisplacement().x()
        << std::setw(20)
        << pt->getInitPosition().y() + pt->getDisplacement().y()
        << std::setw(20)
        << pt->getInitPosition().z() + pt->getDisplacement().z()
        << std::setw(20) << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20) << pt->getVelocity().x()
        << std::setw(20) << pt->getVelocity().y() << std::setw(20)
        << pt->getVelocity().z() << std::setw(20) << vnormL2(pt->getVelocity())
        << std::setw(20) << pressure << std::setw(20) << vonMisesStress
        << std::endl;
    ofs.flush();
  }

  ofs.close();
}

template <>
void
OutputTecplot<SPHParticlePArray>::writeParticles(const SPHParticlePArray* particles, 
                                                 int frame, REAL time)
{
  // Get the filename
  std::string filename(d_sphParticleFilename);
  filename.append(".dat");

  std::ofstream ofs(filename);
  if (!ofs) {
    debugInf << "stream error: printParticlesCSV" << std::endl;
    exit(-1);
  }

  ofs.setf(std::ios::scientific, std::ios::floatfield);
  ofs.precision(OPREC);
  if (frame == 0) {
    ofs << "Title = \"SPH Particle Information\"" << std::endl;
    ofs << "VARIABLES = \"x\", \"y\",\"z\" \"Ux\" \"Uy\" \"Uz\" \"Vx\" \"Vy\" ";
    ofs << "\"Vz\" \"Pressure\" \"a_x\" \"a_y\" \"a_z\" \"density_dot\" ";
    ofs << "\"density\" " << std::endl; 
  }
  ofs << "ZONE T =\" " << frame << "-th Load Step\" " << std::endl;

  // Output the coordinates and the array information
  for (const auto& pt : *particles) {

    pt->calculatePressure();

    ofs << std::setw(20) << pt->currentPosition().x() << std::setw(20)
        << pt->currentPosition().y() << std::setw(20)
        << pt->currentPosition().z() << std::setw(20)
        << pt->getDisplacement().x() << std::setw(20)
        << pt->getDisplacement().y() << std::setw(20)
        << pt->getDisplacement().z() << std::setw(20)
        << pt->getVelocity().x() << std::setw(20) << pt->getVelocity().y()
        << std::setw(20) << pt->getVelocity().z();

    ofs << std::setw(20) << pt->getPressure();

    ofs << std::setw(20) << pt->accelerationeration().x() << std::setw(20)
        << pt->accelerationeration().y() << std::setw(20)
        << pt->accelerationeration().z() << std::setw(20)
        << pt->densityRate() << std::setw(20)
        << pt->density() << std::endl;

    ofs.flush();
  }

  ofs.close();
}

template <typename TArray>
void
OutputTecplot<TArray>::writeSieves(const Gradation* gradation)
{

  // Get the filename
  std::string filename(d_particleFilename);
  filename.append(".dat");

  std::ofstream ofs(filename, std::ofstream::app);
  if (!ofs) {
    debugInf << "stream error: printParticlesCSV" << std::endl;
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

namespace dem {
  template class OutputTecplot<DEMParticlePArray>;
  template class OutputTecplot<pd::PeriParticlePArray>;
  template class OutputTecplot<sph::SPHParticlePArray>;
}