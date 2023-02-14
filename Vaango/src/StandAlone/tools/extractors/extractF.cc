/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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
 *  extractV.cc: Print out a uintah data archive for particle defGrad data
 *
 *  Written by:
 *   Biswajit Banerjee
 *   July 2004
 *
 */

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Math/Matrix3.h>

// #include <Core/Containers/Array3.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/MinMax.h>
#include <Core/OS/Dir.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdio>
#include <ctime>

using namespace std;
using namespace Uintah;

typedef struct
{
  std::vector<Matrix3> defGrad;
  std::vector<long64> id;
  std::vector<double> time;
  std::vector<int> patch;
  std::vector<int> matl;
} MaterialData;

void
usage(const std::string& badarg, const std::string& progname);

void
printDefGrad(DataArchive* da,
             int matID,
             std::vector<long64>& partID,
             string outFile);

int
main(int argc, char** argv)
{
  string partVar;
  int matID = 0;
  string partIDFile;
  string udaDir;
  string outFile;

  // set defaults for cout
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(8);
  /*
   * Parse arguments
   */
  for (int i = 1; i < argc; i++) {
    string s = argv[i];
    if (s == "-m") {
      string id = argv[++i];
      if (id[0] == '-') {
        usage("-m <material id>", argv[0]);
      }
      matID = atoi(argv[i]);
    } else if (s == "-p") {
      partIDFile = argv[++i];
      if (partIDFile[0] == '-') {
        usage("-p <particle id file>", argv[0]);
      }
    } else if (s == "-uda") {
      udaDir = argv[++i];
      if (udaDir[0] == '-') {
        usage("-uda <archive file>", argv[0]);
      }
    } else if (s == "-o") {
      outFile = argv[++i];
      if (outFile[0] == '-') {
        usage("-o <output file>", argv[0]);
      }
    }
  }
  if (argc != 9) {
    usage("", argv[0]);
  }

  std::cout << "Particle Variable to be extracted = p.deformationGradient\n";
  std::cout << "Material ID to be extracted = " << matID << std::endl;

  // Read the particle ID file
  std::cout << "Particle ID File to be read = " << partIDFile << std::endl;
  std::vector<long64> partID;
  ifstream pidFile(partIDFile.c_str());
  if (!pidFile.is_open()) {
    std::cerr << "Particle ID File " << partIDFile << " not found \n";
    exit(1);
  }
  do {
    long64 id = 0;
    pidFile >> id;
    partID.push_back(id);
  } while (!pidFile.eof());

  std::cout << "  Number of Particle IDs = " << partID.size() << std::endl;
  for (unsigned int ii = 0; ii < partID.size() - 1; ++ii) {
    std::cout << "    p" << (ii + 1) << " = " << partID[ii] << std::endl;
  }

  std::cout << "Output file name = " << outFile << std::endl;
  std::cout << "UDA directory to be read = " << udaDir << std::endl;
  try {
    DataArchive* da = scinew DataArchive(udaDir);

    // Print a particular particle variable
    printDefGrad(da, matID, partID, outFile);
  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << std::endl;
    abort();
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    abort();
  }
}
void
usage(const std::string& badarg, const std::string& progname)
{
  if (badarg != "") {
    std::cerr << "Error parsing argument: " << badarg << std::endl;
  }
  std::cerr << "Usage: " << progname << " -m <material id>"
            << " -p <particle id file>"
            << " -uda <archive file>"
            << " -o <output file>\n\n";
  exit(1);
}

////////////////////////////////////////////////////////////////////////////
//
// Print a particle variable
//
////////////////////////////////////////////////////////////////////////////
void
printDefGrad(DataArchive* da,
             int matID,
             std::vector<long64>& partID,
             string outFile)
{

  // Check if the particle variable is available
  std::vector<std::string> vars;
  std::vector<int> num_matls;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matls, types);
  ASSERTEQ(vars.size(), types.size());
  bool variableFound = false;
  string partVar("p.deformationGradient");
  for (unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == partVar) {
      variableFound = true;
    }
  }
  if (!variableFound) {
    std::cerr << "Variable " << partVar << " not found\n";
    exit(1);
  }

  // Create arrays of material data for each particle ID
  MaterialData* matData = scinew MaterialData[partID.size() - 1];
  // for (unsigned int ii = 0; ii < partID.size() ; ++ii) {
  //   matData[ii] = scinew MaterialData();
  // }

  // Now that the variable has been found, get the data for all
  // available time steps from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";

  // Loop thru all the variables
  for (int v = 0; v < (int)vars.size(); v++) {
    std::string var = vars[v];

    // Find the name of the variable
    if (var == partVar) {

      // Loop thru all time steps
      int startPatch = 1;
      for (unsigned long t = 0; t < times.size(); t++) {
        double time = times[t];
        std::cerr << "t = " << time;
        clock_t start = clock();
        GridP grid    = da->queryGrid(t);

        unsigned int numFound = 0;

        // Loop thru all the levels
        for (int l = 0; l < grid->numLevels(); l++) {
          if (numFound == partID.size() - 1) {
            break;
          }

          LevelP level                    = grid->getLevel(l);
          Level::const_patch_iterator iter = level->patchesBegin();
          int patchIndex                  = 0;

          // Loop thru all the patches
          for (; iter != level->patchesEnd(); iter++) {
            if (numFound == partID.size() - 1) {
              break;
            }

            const Patch* patch = *iter;
            ++patchIndex;
            if (patchIndex < startPatch) {
              continue;
            }

            ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
            ConsecutiveRangeSet::iterator matlIter = matls.begin();

            // loop thru all the materials
            for (; matlIter != matls.end(); matlIter++) {
              if (numFound == partID.size() - 1) {
                break;
              }

              int matl = *matlIter;
              if (matl == matID || matl == matID + 1) {
                ParticleVariable<Matrix3> value;
                da->query(value, var, matl, patch, t);
                ParticleSubset* pset = value.getParticleSubset();
                if (pset->numParticles() > 0) {
                  ParticleVariable<long64> pid;
                  da->query(pid, "p.particleID", matl, patch, t);
                  std::vector<bool> found;
                  for (unsigned int ii = 0; ii < partID.size() - 1; ++ii) {
                    found.push_back(false);
                  }
                  ParticleSubset::iterator iter = pset->begin();
                  for (; iter != pset->end(); iter++) {
                    for (unsigned int ii = 0; ii < partID.size() - 1; ++ii) {
                      if (found[ii]) {
                        continue;
                      }
                      if (partID[ii] != pid[*iter]) {
                        continue;
                      }
                      matData[ii].defGrad.push_back(value[*iter]);
                      matData[ii].id.push_back(pid[*iter]);
                      matData[ii].time.push_back(time);
                      matData[ii].patch.push_back(patchIndex);
                      matData[ii].matl.push_back(matl);
                      found[ii] = true;
                      ++numFound;
                      break;
                    }
                    if (numFound == partID.size() - 1) {
                      break;
                    }
                  }
                  if (numFound > 0 && startPatch == 0) {
                    startPatch = patchIndex;
                  }
                }
              } // end of mat compare if
            }   // end of material loop
          }     // end of patch loop
        }       // end of level loop
        clock_t end      = clock();
        double timetaken = (double)(end - start) / (double)CLOCKS_PER_SEC;
        std::cerr << " CPU Time = " << timetaken << " s"
                  << " found " << numFound << std::endl;
      } // end of time step loop
    }   // end of var compare if
  }     // end of variable loop

  // Create output files for each of the particle IDs
  for (unsigned int ii = 0; ii < partID.size() - 1; ++ii) {
    std::ostringstream name;
    name << outFile << "_p" << setw(2) << setfill('0') << (ii + 1);
    ofstream file(name.str().c_str());
    file.setf(ios::scientific, ios::floatfield);
    file.precision(8);
    std::cout << "Created output file " << name.str() << " for particle ID "
              << partID[ii] << std::endl;
    for (unsigned int jj = 0; jj < matData[ii].time.size(); ++jj) {
      double time     = matData[ii].time[jj];
      int patchIndex  = matData[ii].patch[jj];
      int matl        = matData[ii].matl[jj];
      long64 pid      = matData[ii].id[jj];
      Matrix3 defGrad = matData[ii].defGrad[jj];
      file << time << " " << patchIndex << " " << matl;
      file << " " << pid;
      for (int kk = 0; kk < 3; ++kk) {
        for (int ll = 0; ll < 3; ++ll) {
          file << " " << defGrad(kk, ll);
        }
      }
      file << std::endl;
    }
    file.close();
  }
}
