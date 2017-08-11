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
#include <InputOutput/OutputVTK.h>

#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <fstream>
#include <regex>
#include <unistd.h>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace dem;

using PeriParticlePArray = pd::PeriParticlePArray;
using SPHParticlePArray = sph::SPHParticlePArray;

using vtkHexahedronP = vtkSmartPointer<vtkHexahedron>;
using vtkXMLUnstructuredGridWriterP =
  vtkSmartPointer<vtkXMLUnstructuredGridWriter>;
using vtkDoubleArrayP = vtkSmartPointer<vtkDoubleArray>;

template <typename TArray>
OutputVTK<TArray>::OutputVTK(const std::string& folderName, int iterInterval)
  : Output(folderName, iterInterval)
{
  d_domain = nullptr;
  d_grid = nullptr;
  d_particles = nullptr;
  d_cartComm = nullptr;

  // Create the basic file names (path + folder + name)
  createFileNames();

}

template <typename TArray>
OutputVTK<TArray>::~OutputVTK() = default;

template <typename TArray>
void
OutputVTK<TArray>::write(int frame)
{

  // The domain and the grid have to be set before a write is
  // completed.
  if (!d_domain || !d_grid || !d_particles) {
    std::cerr << "**ERROR** Domain and/or Grid and/or Particles have not been "
                 "set.  Nothing "
                 "will be written\n";
    return;
  }

  // Write files for the domain extents at each timestep
  writeDomain(d_domain);

  // Write files for the grid representing each processor at each timestep
  writeGrid(d_grid);

  // Write files for the particle list each timestep
  writeParticles(d_particles, frame);
}

template <typename TArray>
void
OutputVTK<TArray>::writeDomain(const Box* domain)
{

  // Create a writer
  vtkXMLUnstructuredGridWriterP writer = vtkXMLUnstructuredGridWriterP::New();

  // Get the filename
  std::string fileName(d_domainFileName);
  fileName.append(".").append(writer->GetDefaultFileExtension());
  writer->SetFileName(fileName.c_str());

  // Create a pointer to a VTK Unstructured Grid data set
  vtkUnstructuredGridP dataSet = vtkUnstructuredGridP::New();

  // Set up pointer to point data
  vtkPointsP pts = vtkPointsP::New();
  pts->SetNumberOfPoints(8);

  // Add the time
  double time = 0.0;
  addTimeToVTKDataSet(time, dataSet);

  // Add the domain boundary to the unstructured grid cell data
  addDomainToVTKUnstructuredGrid(domain, pts, dataSet);

  // Set the points
  dataSet->SetPoints(pts);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  writer->SetInput(dataSet);
  writer->SetDataModeToAscii();
  writer->Write();
}

template <typename TArray>
void
OutputVTK<TArray>::writeGrid(const Box* grid)
{

  // Create a writer
  vtkXMLUnstructuredGridWriterP writer = vtkXMLUnstructuredGridWriterP::New();

  // Get the filename
  std::string fileName(d_gridFileName);
  fileName.append(".").append(writer->GetDefaultFileExtension());
  writer->SetFileName(fileName.c_str());

  // Create a pointer to a VTK Unstructured Grid data set
  vtkUnstructuredGridP dataSet = vtkUnstructuredGridP::New();

  // Set up pointer to point data
  vtkPointsP pts = vtkPointsP::New();

  // Count the total number of points to be saved
  int num_pts = (d_mpiProcX + 1) * (d_mpiProcY + 1) * (d_mpiProcZ + 1);
  pts->SetNumberOfPoints(num_pts);

  // Add the time
  double time = 0.0;
  addTimeToVTKDataSet(time, dataSet);

  // Create the individual processor domain extents
  Vec v1 = grid->getMinCorner();
  Vec v2 = grid->getMaxCorner();
  Vec vspan = v2 - v1;
  std::vector<Vec> coords((d_mpiProcX + 1) * (d_mpiProcY + 1) *
                          (d_mpiProcZ + 1));
  std::size_t index = 0;
  for (std::size_t i = 0; i < d_mpiProcX + 1; ++i)
    for (std::size_t j = 0; j < d_mpiProcY + 1; ++j)
      for (std::size_t k = 0; k < d_mpiProcZ + 1; ++k)
        coords[index++] = Vec(v1.x() + vspan.x() / d_mpiProcX * i,
                              v1.y() + vspan.y() / d_mpiProcY * j,
                              v1.z() + vspan.z() / d_mpiProcZ * k);

  // Add the processor boundaries to the unstructured grid cell data
  addProcessorsToVTKUnstructuredGrid(coords, pts, dataSet);

  // Set the points
  dataSet->SetPoints(pts);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  writer->SetInput(dataSet);
  writer->SetDataModeToAscii();
  writer->Write();
}

template <typename TArray>
void
OutputVTK<TArray>::writeParticles(const TArray* particles, int frame) 
{
  std::cout << "**ERROR** Noting to do here. The array of particles is"
            << " not of the correct type\n";
}

template <>
void
OutputVTK<DEMParticlePArray>::writeParticles(const DEMParticlePArray* particles, int frame) 
{
  // Create a writer
  vtkXMLUnstructuredGridWriterP writer = vtkXMLUnstructuredGridWriterP::New();

  // Get the filename
  std::string fileName(d_particleFileName);
  fileName.append(".").append(writer->GetDefaultFileExtension());
  writer->SetFileName(fileName.c_str());
  //std::cout << "writeParticles::DEMParticle file = " << fileName << "\n";

  actuallyWriteParticles(particles, frame, writer);
}

template <>
void
OutputVTK<PeriParticlePArray>::writeParticles(const PeriParticlePArray* particles, 
                                              int frame) 
{
  // Create a writer
  vtkXMLUnstructuredGridWriterP writer = vtkXMLUnstructuredGridWriterP::New();

  // Get the filename
  std::string fileName(d_periParticleFileName);
  fileName.append(".").append(writer->GetDefaultFileExtension());
  writer->SetFileName(fileName.c_str());
  //std::cout << "writeParticles::PeriParticle file = " << fileName << "\n";

  actuallyWriteParticles(particles, frame, writer);
}

template <>
void
OutputVTK<SPHParticlePArray>::writeParticles(const SPHParticlePArray* particles, 
                                             int frame) 
{
  // Create a writer
  vtkXMLUnstructuredGridWriterP writer = vtkXMLUnstructuredGridWriterP::New();

  // Get the filename
  std::string fileName(d_sphParticleFileName);
  fileName.append(".").append(writer->GetDefaultFileExtension());
  writer->SetFileName(fileName.c_str());
  //std::cout << "writeParticles::SPHParticle file = " << fileName << "\n";

  actuallyWriteParticles(particles, frame, writer);
}

template <typename TArray>
void
OutputVTK<TArray>::actuallyWriteParticles(const TArray* particles, int frame,
                                          vtkXMLUnstructuredGridWriterP& writer) 
{
  // Create a pointer to a VTK Unstructured Grid data set
  vtkUnstructuredGridP dataSet = vtkUnstructuredGridP::New();

  // Set up pointer to point data
  vtkPointsP pts = vtkPointsP::New();

  // Count the total number of points
  int num_pts = static_cast<int>(particles->size());
  pts->SetNumberOfPoints(num_pts);

  // Add the time
  double time = 0.0;
  addTimeToVTKDataSet(time, dataSet);

  // Add the particle data to the unstructured grid
  createVTKUnstructuredGrid(particles, pts, dataSet);

  // Set the points
  dataSet->SetPoints(pts);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  writer->SetInput(dataSet);
  writer->SetDataModeToAscii();
  writer->Write();
}

template <typename TArray>
void
OutputVTK<TArray>::addTimeToVTKDataSet(double time, vtkUnstructuredGridP& dataSet)
{
  vtkDoubleArrayP array = vtkDoubleArrayP::New();
  array->SetName("TIME");
  array->SetNumberOfTuples(1);
  array->SetTuple1(0, time);
  dataSet->GetFieldData()->AddArray(array);
}

template <typename TArray>
void
OutputVTK<TArray>::addDomainToVTKUnstructuredGrid(const Box* domain, vtkPointsP& pts,
                                                  vtkUnstructuredGridP& dataSet)
{
  double xmin = static_cast<double>(domain->getMinCorner().x());
  double ymin = static_cast<double>(domain->getMinCorner().y());
  double zmin = static_cast<double>(domain->getMinCorner().z());
  double xmax = static_cast<double>(domain->getMaxCorner().x());
  double ymax = static_cast<double>(domain->getMaxCorner().y());
  double zmax = static_cast<double>(domain->getMaxCorner().z());

  int id = 0;
  pts->SetPoint(id, xmin, ymin, zmin);
  ++id;
  pts->SetPoint(id, xmax, ymin, zmin);
  ++id;
  pts->SetPoint(id, xmax, ymax, zmin);
  ++id;
  pts->SetPoint(id, xmin, ymax, zmin);
  ++id;
  pts->SetPoint(id, xmin, ymin, zmax);
  ++id;
  pts->SetPoint(id, xmax, ymin, zmax);
  ++id;
  pts->SetPoint(id, xmax, ymax, zmax);
  ++id;
  pts->SetPoint(id, xmin, ymax, zmax);

  vtkHexahedronP hex = vtkHexahedronP::New();
  for (int ii = 0; ii < 8; ++ii) {
    hex->GetPointIds()->SetId(ii, ii);
  }
  dataSet->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
}

template <typename TArray>
void
OutputVTK<TArray>::addProcessorsToVTKUnstructuredGrid(const std::vector<Vec>& coords,
                                                      vtkPointsP& pts,
                                                      vtkUnstructuredGridP& dataSet)
{
  // Get the number of processes
  int mpiSize = 0;
  MPI_Comm_size(d_cartComm, &mpiSize);

  // Set the coordinates of the corners of the processor boundaries in
  // physical space
  int id = 0;
  for (const auto& coord : coords) {
    pts->SetPoint(id, coord.x(), coord.y(), coord.z());
    ++id;
  }

  size_t nY = d_mpiProcY;
  size_t nZ = d_mpiProcZ;
  vtkHexahedronP hex = vtkHexahedronP::New();
  for (int rank = 0; rank < mpiSize; ++rank) {

    // Get process coords
    int procCoords[3];
    MPI_Cart_coords(d_cartComm, rank, 3, procCoords);

    // Compute node ids
    int id1 = 1 + (procCoords[0] + 1) * (nZ + 1) * (nY + 1) +
              procCoords[1] * (nZ + 1) + procCoords[2];
    int id3 = 1 + procCoords[0] * (nZ + 1) * (nY + 1) +
              (procCoords[1] + 1) * (nZ + 1) + procCoords[2];
    int id2 = 1 + (procCoords[0] + 1) * (nZ + 1) * (nY + 1) +
              (procCoords[1] + 1) * (nZ + 1) + procCoords[2];
    int id4 = 1 + procCoords[0] * (nZ + 1) * (nY + 1) +
              procCoords[1] * (nZ + 1) + procCoords[2];

    int id5 = 1 + (procCoords[0] + 1) * (nZ + 1) * (nY + 1) +
              procCoords[1] * (nZ + 1) + (procCoords[2] + 1);
    int id7 = 1 + procCoords[0] * (nZ + 1) * (nY + 1) +
              (procCoords[1] + 1) * (nZ + 1) + (procCoords[2] + 1);
    int id6 = 1 + (procCoords[0] + 1) * (nZ + 1) * (nY + 1) +
              (procCoords[1] + 1) * (nZ + 1) + (procCoords[2] + 1);
    int id8 = 1 + procCoords[0] * (nZ + 1) * (nY + 1) +
              procCoords[1] * (nZ + 1) + (procCoords[2] + 1);

    // Set the node ids for this hexahedron
    hex->GetPointIds()->SetId(0, id1 - 1);
    hex->GetPointIds()->SetId(1, id2 - 1);
    hex->GetPointIds()->SetId(2, id3 - 1);
    hex->GetPointIds()->SetId(3, id4 - 1);
    hex->GetPointIds()->SetId(4, id5 - 1);
    hex->GetPointIds()->SetId(5, id6 - 1);
    hex->GetPointIds()->SetId(6, id7 - 1);
    hex->GetPointIds()->SetId(7, id8 - 1);

    dataSet->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
  }
}

template <typename TArray>
void
OutputVTK<TArray>::createVTKUnstructuredGrid(const TArray* particles,
                                             vtkPointsP& pts,
                                             vtkUnstructuredGridP& dataSet)
{
  std::cout << "**ERROR** Noting to do here. The array of particles is"
            << " not of the correct type\n";
}

template <>
void
OutputVTK<DEMParticlePArray>::createVTKUnstructuredGrid(const DEMParticlePArray* particles,
                                                     vtkPointsP& pts,
                                                     vtkUnstructuredGridP& dataSet)
{
  // Set up pointers for material property data
  vtkDoubleArrayP ID = vtkDoubleArrayP::New();
  ID->SetNumberOfComponents(1);
  ID->SetNumberOfTuples(pts->GetNumberOfPoints());
  ID->SetName("ID");

  vtkDoubleArrayP type = vtkDoubleArrayP::New();
  type->SetNumberOfComponents(1);
  type->SetNumberOfTuples(pts->GetNumberOfPoints());
  type->SetName("Type");

  vtkDoubleArrayP radii = vtkDoubleArrayP::New();
  radii->SetNumberOfComponents(3);
  radii->SetNumberOfTuples(pts->GetNumberOfPoints());
  radii->SetName("Radius");

  vtkDoubleArrayP position = vtkDoubleArrayP::New();
  position->SetNumberOfComponents(3);
  position->SetNumberOfTuples(pts->GetNumberOfPoints());
  position->SetName("Position");

  vtkDoubleArrayP axis_a = vtkDoubleArrayP::New();
  axis_a->SetNumberOfComponents(3);
  axis_a->SetNumberOfTuples(pts->GetNumberOfPoints());
  axis_a->SetName("Axis a");

  vtkDoubleArrayP axis_b = vtkDoubleArrayP::New();
  axis_b->SetNumberOfComponents(3);
  axis_b->SetNumberOfTuples(pts->GetNumberOfPoints());
  axis_b->SetName("Axis b");

  vtkDoubleArrayP axis_c = vtkDoubleArrayP::New();
  axis_c->SetNumberOfComponents(3);
  axis_c->SetNumberOfTuples(pts->GetNumberOfPoints());
  axis_c->SetName("Axis c");

  vtkDoubleArrayP velocity = vtkDoubleArrayP::New();
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(pts->GetNumberOfPoints());
  velocity->SetName("Velocity");

  vtkDoubleArrayP omega = vtkDoubleArrayP::New();
  omega->SetNumberOfComponents(3);
  omega->SetNumberOfTuples(pts->GetNumberOfPoints());
  omega->SetName("Omega");

  vtkDoubleArrayP force = vtkDoubleArrayP::New();
  force->SetNumberOfComponents(3);
  force->SetNumberOfTuples(pts->GetNumberOfPoints());
  force->SetName("Force");

  vtkDoubleArrayP moment = vtkDoubleArrayP::New();
  moment->SetNumberOfComponents(3);
  moment->SetNumberOfTuples(pts->GetNumberOfPoints());
  moment->SetName("Moment");

  // Loop through particles
  Vec vObj;
  int id = 0;
  double vec[3];
  for (const auto& particle : *particles) {

    // Position
    vObj = particle->currentPosition();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    pts->SetPoint(id, vec);

    // ID
    ID->InsertValue(id, particle->getId());

    // Type
    type->InsertValue(id, particle->getType());

    // Ellipsoid radii
    vec[0] = particle->getA();
    vec[1] = particle->getB();
    vec[2] = particle->getC();
    radii->InsertTuple(id, vec);

    // Current direction A
    vObj = particle->getCurrDirecA();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    axis_a->InsertTuple(id, vec);

    // Current direction B
    vObj = particle->getCurrDirecB();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    axis_b->InsertTuple(id, vec);

    // Current direction C
    vObj = particle->getCurrDirecC();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    axis_c->InsertTuple(id, vec);

    // Velocity
    vObj = particle->currentVel();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    velocity->InsertTuple(id, vec);

    // Omega
    vObj = particle->currentOmega();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    omega->InsertTuple(id, vec);

    // Force
    vObj = particle->getForce();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    force->InsertTuple(id, vec);

    // Moment
    vObj = particle->getMoment();
    vec[0] = vObj.x();
    vec[1] = vObj.y();
    vec[2] = vObj.z();
    moment->InsertTuple(id, vec);

    ++id;
  }

  // Add points to data set
  dataSet->GetPointData()->AddArray(ID);
  dataSet->GetPointData()->AddArray(type);
  dataSet->GetPointData()->AddArray(radii);
  dataSet->GetPointData()->AddArray(axis_a);
  dataSet->GetPointData()->AddArray(axis_b);
  dataSet->GetPointData()->AddArray(axis_c);
  dataSet->GetPointData()->AddArray(velocity);
  dataSet->GetPointData()->AddArray(omega);
  dataSet->GetPointData()->AddArray(force);
  dataSet->GetPointData()->AddArray(moment);

  // Check point data
  /*
  vtkPointData *pd = dataSet->GetPointData();
  if (pd) {
    //std::cout << " contains point data with " << pd->GetNumberOfArrays() << "
  arrays." << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); i++) {
      //std::cout << "\tArray " << i << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;
    }
  }
  */
}

template <>
void
OutputVTK<PeriParticlePArray>::createVTKUnstructuredGrid(const PeriParticlePArray* particles,
                                                         vtkPointsP& pts,
                                                         vtkUnstructuredGridP& dataSet)
{
  // Set up pointers for material property data
  vtkDoubleArrayP ID = vtkDoubleArrayP::New();
  ID->SetNumberOfComponents(1);
  ID->SetNumberOfTuples(pts->GetNumberOfPoints());
  ID->SetName("ID");

  vtkDoubleArrayP position = vtkDoubleArrayP::New();
  position->SetNumberOfComponents(3);
  position->SetNumberOfTuples(pts->GetNumberOfPoints());
  position->SetName("Position");

  vtkDoubleArrayP displacement = vtkDoubleArrayP::New();
  displacement->SetNumberOfComponents(3);
  displacement->SetNumberOfTuples(pts->GetNumberOfPoints());
  displacement->SetName("Displacement");

  vtkDoubleArrayP velocity = vtkDoubleArrayP::New();
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(pts->GetNumberOfPoints());
  velocity->SetName("Velocity");

  vtkDoubleArrayP kineticEnergy = vtkDoubleArrayP::New();
  kineticEnergy->SetNumberOfComponents(1);
  kineticEnergy->SetNumberOfTuples(pts->GetNumberOfPoints());
  kineticEnergy->SetName("KineticEnergy");

  vtkDoubleArrayP meanStress = vtkDoubleArrayP::New();
  meanStress->SetNumberOfComponents(1);
  meanStress->SetNumberOfTuples(pts->GetNumberOfPoints());
  meanStress->SetName("MeanStress");

  vtkDoubleArrayP vonMisesStress = vtkDoubleArrayP::New();
  vonMisesStress->SetNumberOfComponents(1);
  vonMisesStress->SetNumberOfTuples(pts->GetNumberOfPoints());
  vonMisesStress->SetName("VonMisesStress");

  // Loop through particles
  int id = 0;
  double vec[3];
  for (const auto& particle : *particles) {

    // Get the vector quantities
    Vec pos0 = particle->getInitPosition();
    Vec disp = particle->getDisplacement();
    Vec vel = particle->getVelocity();

    // Compute derived quantities (mean stress and von mises stress)
    Matrix sigma = particle->getSigma();
    REAL sigma_m = (sigma(1,1) + sigma(2,2) + sigma(3,3))/3.0;
    REAL s1s2 = sigma(1,1) - sigma(2,2);
    REAL s2s3 = sigma(2,2) - sigma(3,3);
    REAL s3s1 = sigma(3,3) - sigma(1,1);
    REAL s12Sq = sigma(1,2)*sigma(1,2);
    REAL s23Sq = sigma(2,3)*sigma(2,3);
    REAL s31Sq = sigma(3,1)*sigma(3,1);
    REAL sigma_eq = std::sqrt(0.5*(s1s2*s1s2 + s2s3*s2s3 + s3s1*s3s1 +
                                   3.0*(s12Sq + s23Sq + s31Sq)));

    // Compute kinetic energy
    REAL mass = particle->getMass();
    REAL KE = 0.5 * mass * dot(vel, vel);

    // Position
    Vec pos = pos0 + disp;
    vec[0] = pos.x();
    vec[1] = pos.y();
    vec[2] = pos.z();
    pts->SetPoint(id, vec);

    // ID
    ID->InsertValue(id, particle->getId());

    // Displacement
    vec[0] = disp.x();
    vec[1] = disp.y();
    vec[2] = disp.z();
    displacement->InsertTuple(id, vec);

    // Velocity
    vec[0] = vel.x();
    vec[1] = vel.y();
    vec[2] = vel.z();
    velocity->InsertTuple(id, vec);

    // Kinetic energy
    kineticEnergy->InsertValue(id, KE);

    // Mean stress
    meanStress->InsertValue(id, sigma_m);

    // von Mises stress
    vonMisesStress->InsertValue(id, sigma_eq);

    ++id;
  }

  // Add points to data set
  dataSet->GetPointData()->AddArray(ID);
  dataSet->GetPointData()->AddArray(displacement);
  dataSet->GetPointData()->AddArray(velocity);
  dataSet->GetPointData()->AddArray(kineticEnergy);
  dataSet->GetPointData()->AddArray(meanStress);
  dataSet->GetPointData()->AddArray(vonMisesStress);

  // Check point data
  /*
  vtkPointData *pd = dataSet->GetPointData();
  if (pd) {
    //std::cout << " contains point data with " << pd->GetNumberOfArrays() << "
  arrays." << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); i++) {
      //std::cout << "\tArray " << i << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;
    }
  }
  */
}

template <>
void
OutputVTK<SPHParticlePArray>::createVTKUnstructuredGrid(const SPHParticlePArray* particles,
                                                        vtkPointsP& pts,
                                                        vtkUnstructuredGridP& dataSet)
{
  // Set up pointers for material property data
  vtkDoubleArrayP ID = vtkDoubleArrayP::New();
  ID->SetNumberOfComponents(1);
  ID->SetNumberOfTuples(pts->GetNumberOfPoints());
  ID->SetName("ID");

  vtkDoubleArrayP position = vtkDoubleArrayP::New();
  position->SetNumberOfComponents(3);
  position->SetNumberOfTuples(pts->GetNumberOfPoints());
  position->SetName("Position");

  vtkDoubleArrayP displacement = vtkDoubleArrayP::New();
  displacement->SetNumberOfComponents(3);
  displacement->SetNumberOfTuples(pts->GetNumberOfPoints());
  displacement->SetName("Displacement");

  vtkDoubleArrayP velocity = vtkDoubleArrayP::New();
  velocity->SetNumberOfComponents(3);
  velocity->SetNumberOfTuples(pts->GetNumberOfPoints());
  velocity->SetName("Velocity");

  vtkDoubleArrayP pressure = vtkDoubleArrayP::New();
  pressure->SetNumberOfComponents(1);
  pressure->SetNumberOfTuples(pts->GetNumberOfPoints());
  pressure->SetName("Pressure");

  vtkDoubleArrayP acceleration = vtkDoubleArrayP::New();
  acceleration->SetNumberOfComponents(3);
  acceleration->SetNumberOfTuples(pts->GetNumberOfPoints());
  acceleration->SetName("Acceleration");

  vtkDoubleArrayP densityRate = vtkDoubleArrayP::New();
  densityRate->SetNumberOfComponents(1);
  densityRate->SetNumberOfTuples(pts->GetNumberOfPoints());
  densityRate->SetName("DensityRate");

  vtkDoubleArrayP density = vtkDoubleArrayP::New();
  density->SetNumberOfComponents(1);
  density->SetNumberOfTuples(pts->GetNumberOfPoints());
  density->SetName("Density");

  // Loop through particles
  int id = 0;
  double vec[3];
  for (const auto& particle : *particles) {

    // Compute derived quantities (pressure)
    particle->calculatePressure();

    // Get the vector quantities
    Vec pos  = particle->currentPosition();
    Vec disp = particle->getDisplacement();
    Vec vel = particle->getVelocity();
    Vec acc = particle->getAcceleration();
    REAL press = particle->getPressure();
    REAL rhoDot = particle->getDensityDot();
    REAL rho = particle->getDensity();

    // Position
    vec[0] = pos.x();
    vec[1] = pos.y();
    vec[2] = pos.z();
    pts->SetPoint(id, vec);

    // ID
    ID->InsertValue(id, particle->getId());

    // Displacement
    vec[0] = disp.x();
    vec[1] = disp.y();
    vec[2] = disp.z();
    displacement->InsertTuple(id, vec);

    // Velocity
    vec[0] = vel.x();
    vec[1] = vel.y();
    vec[2] = vel.z();
    velocity->InsertTuple(id, vec);

    // Pressure
    pressure->InsertValue(id, press);

    // Acceleration
    vec[0] = acc.x();
    vec[1] = acc.y();
    vec[2] = acc.z();
    acceleration->InsertTuple(id, vec);

    // Density rate
    densityRate->InsertValue(id, rhoDot);

    // Density
    density->InsertValue(id, rho);

    ++id;
  }

  // Add points to data set
  dataSet->GetPointData()->AddArray(ID);
  dataSet->GetPointData()->AddArray(displacement);
  dataSet->GetPointData()->AddArray(velocity);
  dataSet->GetPointData()->AddArray(pressure);
  dataSet->GetPointData()->AddArray(acceleration);
  dataSet->GetPointData()->AddArray(densityRate);
  dataSet->GetPointData()->AddArray(density);

  // Check point data
  /*
  vtkPointData *pd = dataSet->GetPointData();
  if (pd) {
    //std::cout << " contains point data with " << pd->GetNumberOfArrays() << "
  arrays." << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); i++) {
      //std::cout << "\tArray " << i << " is named "
                << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;
    }
  }
  */
}

namespace dem {
  template class OutputVTK<DEMParticlePArray>;
  template class OutputVTK<pd::PeriParticlePArray>;
  template class OutputVTK<sph::SPHParticlePArray>;
}
