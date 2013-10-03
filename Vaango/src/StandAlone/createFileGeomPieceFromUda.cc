/*
 *  createFileGeomPieceFromUda.cc 
 *
 *     Read a uintah data archive and create a file geometry piece containing the
 *     particle positions and volumes
 *
 *  Written by:
 *   Biswajit Banerjee
 *   September 2013
 *
 */

#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace SCIRun;
using namespace Uintah;

// Structure to store particle data
typedef struct {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> volume;
} FileGeomParticleData;

// declarations
void usage(const std::string& badarg, const std::string& progname);
void readFileGeomParticleData(DataArchive* da, FileGeomParticleData& partData);
void writeFileGeomParticleData(FileGeomParticleData& partData, const std::string& output_file);

//-------------------------------------------------------------------------------------
// Usage: createFileGeomPieceFromUda <uda_file> <output_file>
//-------------------------------------------------------------------------------------
void usage(const std::string& progname)
{
  std::cerr << "Incorrect inputs to program:" << std::endl;
  std::cerr << "Correct usage: " << progname << " <uda file> <output_file>" << std::endl;
  std::cerr << "    File will contain [x y z volume] for each particle from the" 
            << "       last timestep in the uda file." << std::endl; 
  exit(1);
}


//-------------------------------------------------------------------------------------
// Main
//-------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // In correct number of arguments
  if (argc != 3) usage(argv[0]);

  // Parse arguments
  std::string uda_file = argv[1];
  std::string output_file = argv[2];

  // Read the data archive and print out 
  try {
    DataArchive* da = scinew DataArchive(uda_file);

    FileGeomParticleData part_data;
    readFileGeomParticleData(da, part_data);  
    
    writeFileGeomParticleData(part_data, output_file);

    delete da;

  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << endl;
    abort();
  } catch(...){
    std::cerr << "Caught unknown exception\n";
    abort();
  }
}

//-------------------------------------------------------------------------------
// Read the particle positions and volumes
//-------------------------------------------------------------------------------
void readFileGeomParticleData(DataArchive* da, 
                      FileGeomParticleData& part_data)
{
  // set defaults for cout
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(8);

  // Check if the particle variables p.position and p.volume are available
  // in the uda
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());

  bool positionFound = false;
  bool volumeFound = false;
  std::string posVar("p.x");
  std::string volVar("p.volume");
  for(unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == posVar) positionFound = true;
    if (var == volVar) volumeFound = true;
  }
  if (!positionFound) {
    std::cerr << "Variable " << posVar << " not found in the uda.\n"; 
    exit(1);
  }
  if (!volumeFound) {
    std::cerr << "Variable " << volVar << " not found in the uda.\n"; 
    exit(1);
  }

  // Now that the variable has been found, get the data for the final timestep from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";
  unsigned int timeID = times.size() - 1;
  double time = times[timeID];
  std::cout << "Time = " << time << endl;

  // Get grid info for time t
  GridP grid = da->queryGrid(timeID);

  // Loop thru all the levels
  for(int levelID=0; levelID<grid->numLevels(); levelID++) {

    // Get level info
    LevelP level = grid->getLevel(levelID);

    // Loop thru all the patches
    Level::const_patchIterator patchIter = level->patchesBegin(); 
    for(; patchIter != level->patchesEnd(); patchIter++) {

      // Get patch
      const Patch* patch = *patchIter;

      // First find the particle positions for the first material in the uda
      ConsecutiveRangeSet matls = da->queryMaterials(posVar, patch, timeID);
      int matl = *(matls.begin());

      ParticleVariable<Point> position;
      da->query(position, posVar, matl, patch, timeID);
      ParticleSubset* pset = position.getParticleSubset();
      if (pset->numParticles() > 0) {
        ParticleSubset::iterator partIter = pset->begin();
        for (; partIter != pset->end(); partIter++) {
          Point pos = position[*partIter];
          part_data.x.push_back(pos.x());
          part_data.y.push_back(pos.y());
          part_data.z.push_back(pos.z());
        }
      }

      // Next find the particle volumes for the first material in the uda
      matls = da->queryMaterials(volVar, patch, timeID);
      matl = *(matls.begin());

      ParticleVariable<double> volume;
      da->query(volume, volVar, matl, patch, timeID);
      pset = volume.getParticleSubset();
      if (pset->numParticles() > 0) {
        ParticleSubset::iterator partIter = pset->begin();
        for (; partIter != pset->end(); partIter++) {
          double vol = volume[*partIter];
          part_data.volume.push_back(vol);
        }
      }
    } // end of patch loop
  } // end of level loop
}

//------------------------------------------------------------------------------------
// Write out the data
//------------------------------------------------------------------------------------
void writeFileGeomParticleData(FileGeomParticleData& partData, const std::string& output_file)
{
  std::ofstream file(output_file.c_str());
  file.setf(std::ios::scientific, std::ios::floatfield);
  file.precision(8);
  for (unsigned int jj = 0; jj < partData.x.size(); ++jj) {
    file << partData.x[jj] << " " << partData.y[jj] << " " << partData.z[jj] << " " << partData.volume[jj];
    file << endl;
  }
  file.close();
  std::cout << "Wrote output file " << output_file << std::endl;
}
