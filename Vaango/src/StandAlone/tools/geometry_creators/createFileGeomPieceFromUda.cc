/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

/*
 *  createFileGeomPieceFromUda.cc
 *
 *     Read a uintah data archive and create a file geometry piece containing
 * the
 *     particle positions and volumes
 *
 *  Written by:
 *   Biswajit Banerjee
 *   September 2013
 *
 */

#include <Core/DataArchive/DataArchive.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Uintah;

// Structure to store particle data
typedef struct
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> volume;
} FileGeomParticleData;

// declarations
void
usage(const std::string& badarg, const std::string& progname);
int
getNumberOfMaterials(DataArchive* da);
void
readFileGeomParticleData(DataArchive* da,
                         int matID,
                         FileGeomParticleData& partData);
void
writeFileGeomParticleData(int matID,
                          FileGeomParticleData& partData,
                          const std::string& output_file);

//-------------------------------------------------------------------------------------
// Usage: createFileGeomPieceFromUda <uda_file> <output_file>
//-------------------------------------------------------------------------------------
void
usage(const std::string& progname)
{
  std::cerr << "Incorrect inputs to program:" << std::endl;
  std::cerr << "Correct usage: " << progname << " <uda file> <output_file>"
            << std::endl;
  std::cerr << "    File will contain [x y z volume] for each particle from the"
            << "       last timestep in the uda file." << std::endl;
  exit(1);
}

//-------------------------------------------------------------------------------------
// Main
//-------------------------------------------------------------------------------------
int
main(int argc, char** argv)
{
  // In correct number of arguments
  if (argc != 3) {
    usage(argv[0]);
  }

  // Parse arguments
  std::string uda_file    = argv[1];
  std::string output_file = argv[2];

  // Read the data archive and print out
  try {
    DataArchive* da = scinew DataArchive(uda_file);

    // Find number of materials
    int numMat = getNumberOfMaterials(da);

    // loop thru all the materials
    for (int matID = 0; matID < numMat; matID++) {
      FileGeomParticleData part_data;
      readFileGeomParticleData(da, matID, part_data);
      writeFileGeomParticleData(matID, part_data, output_file);
    }

    delete da;

  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << endl;
    abort();
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    abort();
  }
}

int
getNumberOfMaterials(DataArchive* da)
{

  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());

  unsigned int timeID = times.size() - 1;
  double time         = times[timeID];
  std::cout << "Time = " << time << endl;

  std::string posVar("p.x");
  int numMat = 0;
  GridP grid = da->queryGrid(timeID);
  for (int levelID = 0; levelID < grid->numLevels(); levelID++) {
    LevelP level = grid->getLevel(levelID);
    for (auto patchIter = level->patchesBegin();
         patchIter != level->patchesEnd();
         patchIter++) {
      const Patch* patch        = *patchIter;
      ConsecutiveRangeSet matls = da->queryMaterials(posVar, patch, timeID);
      int localNumMat           = 0;
      for ([[maybe_unused]] auto matl : matls) {
        ++localNumMat;
      }
      if (localNumMat > numMat) {
        numMat = localNumMat;
      }
    }
  }
  return numMat;
}

//-------------------------------------------------------------------------------
// Read the particle positions and volumes
//-------------------------------------------------------------------------------
void
readFileGeomParticleData(DataArchive* da,
                         int matID,
                         FileGeomParticleData& part_data)
{
  // set defaults for cout
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(8);

  // Check if the particle variables p.position and p.volume are available
  // in the uda
  std::vector<std::string> vars;
  std::vector<int> num_matl;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matl, types);
  ASSERTEQ(vars.size(), types.size());

  bool positionFound = false;
  bool volumeFound   = false;
  std::string posVar("p.x");
  std::string volVar("p.volume");
  for (auto var : vars) {
    if (var == posVar) {
      positionFound = true;
    }
    if (var == volVar) {
      volumeFound = true;
    }
  }
  if (!positionFound) {
    std::cerr << "**ERROR** Variable " << posVar << " not found in the uda.\n";
    exit(1);
  }
  if (!volumeFound) {
    std::cerr << "**ERROR** Variable " << volVar << " not found in the uda.\n";
    exit(1);
  }

  // Now that the variable has been found, get the data for the final timestep
  // from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";
  unsigned int timeID = times.size() - 1;
  double time         = times[timeID];
  std::cout << "Time = " << time << endl;

  // Get grid info for time t
  GridP grid = da->queryGrid(timeID);

  // Loop thru all the levels
  for (int levelID = 0; levelID < grid->numLevels(); levelID++) {

    // Get level info
    LevelP level = grid->getLevel(levelID);

    // Loop thru all the patches
    Level::const_patch_iterator patchIter = level->patchesBegin();
    for (; patchIter != level->patchesEnd(); patchIter++) {

      // Get patch
      const Patch* patch = *patchIter;

      // First find the particle positions for the first material in the uda
      ConsecutiveRangeSet matls = da->queryMaterials(posVar, patch, timeID);
      for (auto matl : matls) {
        if (matl == matID) {

          ParticleVariable<Point> position;
          ParticleVariable<double> volume;
          da->query(position, posVar, matl, patch, timeID);
          da->query(volume, volVar, matl, patch, timeID);
          ParticleSubset* pset = position.getParticleSubset();
          if (pset->numParticles() > 0) {
            for (auto particle : *pset) {
              Point pos = position[particle];
              part_data.x.push_back(pos.x());
              part_data.y.push_back(pos.y());
              part_data.z.push_back(pos.z());
              double vol = volume[particle];
              part_data.volume.push_back(vol);
            }
          }
        }
      }
    } // end of patch loop
  }   // end of level loop
}

//------------------------------------------------------------------------------------
// Write out the data
//------------------------------------------------------------------------------------
void
writeFileGeomParticleData(int matID,
                          FileGeomParticleData& partData,
                          const std::string& output_file)
{
  std::string fileName = output_file + ".mat" + std::to_string(matID);
  std::ofstream file(fileName);
  file.setf(std::ios::scientific, std::ios::floatfield);
  file.precision(8);
  for (unsigned int jj = 0; jj < partData.x.size(); ++jj) {
    file << partData.x[jj] << " " << partData.y[jj] << " " << partData.z[jj]
         << " " << partData.volume[jj];
    file << endl;
  }
  file.close();
  std::cout << "Wrote output file " << output_file << std::endl;
}
