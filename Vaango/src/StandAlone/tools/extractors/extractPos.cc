/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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
 *  extractPos.cc: 
 *
 *   Input:
 *    1) A selection of particle IDs using selectpart with one timestep only
 *    2) Material id, uda etc.
 *   Output:
 *    1) A file containing the positions of all these particles at
 *       a given time step
 *
 *  Written by:
 *   Biswajit Banerjee
 *   Sep 2015
 *
 */

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Math/Matrix3.h>

#include <Core/Containers/ConsecutiveRangeSet.h>
#include <Core/Containers/Array3.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Math/MinMax.h>
#include <Core/OS/Dir.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace Uintah;
using namespace std;
using namespace Uintah;

typedef struct{
  vector<Point> position;
  vector<long64> id;
  vector<double> time;
  vector<int> patch;
  vector<int> matl;
} MaterialData;

void usage(const std::string& badarg, const std::string& progname);

void printPosition(DataArchive* da, 
                   int matID,
                   unsigned long timestep,
                   vector<long64>& partID,
                   string outFile);

int main(int argc, char** argv)
{
  string partVar;
  int matID = 0;
  string partIDFile;
  string udaDir;
  string outFile;
  unsigned long timeStep = 0;

  // set defaults for cout
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(8);
  /*
   * Parse arguments
   */
  cerr << "Particle Variable to be extracted = p.x\n";
  for(int i=1;i<argc;i++){
    string s=argv[i];
    if (s == "-m") {
      string id = argv[++i];
      if (id[0] == '-')  
        usage("-m <material id>", argv[0]);
      matID = atoi(argv[i]);
    } else if(s == "-p"){
      partIDFile = argv[++i];
      if (partIDFile[0] == '-') 
        usage("-p <particle id file>", argv[0]);
    } else if (s == "-timestep") {
      timeStep = std::stoul(argv[++i]);
    } else if(s == "-uda"){
      udaDir = argv[++i];
      if (udaDir[0] == '-') 
        usage("-uda <archive file>", argv[0]);
    } else if(s == "-o"){
      outFile = argv[++i];
      if (outFile[0] == '-') 
        usage("-o <output file>", argv[0]);
    } 
  }
  cerr << "Number of arguments = " << argc << std::endl;
  if (argc != 10) usage( "", argv[0] );

  cerr << "Material ID to be extracted = " << matID << endl;
  cerr << "Timestep to be extracted = " << timeStep << endl;

  // Read the particle ID file
  cerr << "Particle ID File to be read = " << partIDFile << endl;
  vector<long64> partID;
  ifstream pidFile(partIDFile.c_str());
  if (!pidFile.is_open()) {
    cerr << "Particle ID File " << partIDFile << " not found \n";
    exit(1);
  }
  do {
    long64 id = 0;
    double t = 0.0, x = 0.0, y = 0.0, z = 0.0;
    int patch = 0, mat = 0;
    pidFile >> t >> patch >> mat >> id >> x >> y >> z;
    //pidFile >> id;
    partID.push_back(id);
  } while (!pidFile.eof());
  
  cerr << "  Number of Particle IDs = " << partID.size() << endl;
  for (unsigned int ii = 0; ii < partID.size()-1 ; ++ii) {
    cerr << "    p"<< (ii+1) << " = " << partID[ii] << endl;
  }

  cerr << "Output file name = " << outFile << endl;
  cerr << "UDA directory to be read = " << udaDir << endl;
  try {
    DataArchive* da = scinew DataArchive(udaDir);
    
    // Print a particular particle variable
    printPosition(da, matID, timeStep, partID, outFile);
  } catch (Exception& e) {
    cerr << "Caught exception: " << e.message() << endl;
    abort();
  } catch(...){
    cerr << "Caught unknown exception\n";
    abort();
  }
}
void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") cerr << "Error parsing argument: " << badarg << endl;
  cerr << "Usage: " << progname 
       << " -m <material id>"
       << " -p <particle id file>"
       << " -timestep <timestep #>"
       << " -uda <archive file>"
       << " -o <output file>\n\n";
  exit(1);
}

////////////////////////////////////////////////////////////////////////////
//
// Print a particle variable
//
////////////////////////////////////////////////////////////////////////////
void printPosition(DataArchive* da, 
                   int matID,
                   unsigned long timeStep,
                   vector<long64>& partID,
                   string outFile){

  // Check if the particle variable is available
  vector<string> vars;
  vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());
  bool variableFound = false;
  string partVar("p.x");
  for(unsigned int v=0;v<vars.size();v++){
    std::string var = vars[v];
    if (var == partVar) variableFound = true;
  }
  if (!variableFound) {
    cerr << "Variable " << partVar << " not found\n"; 
    exit(1);
  }

  // Now that the variable has been found, get the data for the 
  // required time step from the data archive
  vector<int> index;
  vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  cerr << "There are " << index.size() << " timesteps:\n";

  // Check that the input timestep exists else quit
  if (timeStep > times.size()-1) {
    std::cerr << "timeStep " << timeStep << " not in range " 
              << "[" << 0 << "," << times.size()-1 << "]" << std::endl;
    usage( "", "extractPos" );
  }

  // Create arrays of material data for each time step
  MaterialData* matData = scinew MaterialData[partID.size()-1];
  //for (unsigned int ii = 0; ii < partID.size() ; ++ii) {
  //  matData[ii] = scinew MaterialData();
  //}

  // Loop thru all the variables 
  for(int v=0;v<(int)vars.size();v++){
    std::string var = vars[v];

    // Find the name of the variable
    if (var == partVar) {

      // Extract the input timestep info
      int startPatch = 1;
      unsigned long t = timeStep;
      double time = times[t];
      cerr << "t = " << time ;
      clock_t start = clock();
      GridP grid = da->queryGrid(t);

      unsigned int numFound = 0;

      // Loop thru all the levels
      for(int l=0;l<grid->numLevels();l++){
        if (numFound == partID.size()-1) break;

        LevelP level = grid->getLevel(l);
        int patchIndex = 0;

        // Loop thru all the patches
        for (auto iter = level->patchesBegin(); iter != level->patchesEnd(); iter++){
          if (numFound == partID.size()-1) break;

          const Patch* patch = *iter;
          ++patchIndex; 
          if (patchIndex < startPatch) continue;

          ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);

          // loop thru all the materials
          for (auto matlIter = matls.begin(); matlIter != matls.end(); matlIter++){
            if (numFound == partID.size()-1) break;

            int matl = *matlIter;
            if (matl == matID) {
              ParticleVariable<Point> position;
              da->query(position, "p.x", matl, patch, t);
              ParticleSubset* pset = position.getParticleSubset();
              if(pset->numParticles() > 0){
                ParticleVariable<long64> pid;
                da->query(pid, "p.particleID", matl, patch, t);
                vector<bool> found;
                for (unsigned int ii = 0; ii < partID.size()-1 ; ++ii) {
                  found.push_back(false);
                }
                for(auto iter = pset->begin();iter != pset->end(); iter++){
                  for (unsigned int ii = 0; ii < partID.size()-1 ; ++ii) {
                    if (found[ii]) continue;
                    if (partID[ii] != pid[*iter]) continue;
                      matData[ii].position.push_back(position[*iter]);
                      matData[ii].id.push_back(pid[*iter]);
                      matData[ii].time.push_back(time);
                      matData[ii].patch.push_back(patchIndex);
                      matData[ii].matl.push_back(matl);
                      //cout << time << " " << patchIndex << " " << matl << " " << pid[*iter]
                      //     << " " << position[*iter].x() << " " << position[*iter].y() 
                      //     << " " << position[*iter].z() << endl;
                      found[ii] = true;
                      ++numFound;
                      break;
                    }
                    if (numFound == partID.size()-1) break;
                  }
                  if (numFound > 0 && startPatch == 0) startPatch = patchIndex;
                }
              } // end of mat compare if
            } // end of material loop
          } // end of patch loop
        } // end of level loop
        clock_t end = clock();
        double timetaken = (double) (end - start)/(double) CLOCKS_PER_SEC;
        cerr << " CPU Time = " << timetaken << " s" << " found " 
             << numFound << endl;
    } // end of var compare if
  } // end of variable loop

  // Write the output to stdout
  for (unsigned int ii = 0; ii < partID.size()-1 ; ++ii) {
    for (unsigned int jj = 0; jj < matData[ii].time.size(); ++jj) {
      double time = matData[ii].time[jj];
      int patchIndex = matData[ii].patch[jj];
      int matl = matData[ii].matl[jj];
      long64 pid = matData[ii].id[jj];
      Point pos = matData[ii].position[jj];
      cout << time << " " << patchIndex << " " << matl << " " << pid
           << " " << pos.x() << " " << pos.y() 
           << " " << pos.z() << endl;
    }
  }
  // Create output files for each of the particle IDs
  /*
  for (unsigned int ii = 0; ii < partID.size()-1 ; ++ii) {
    ostringstream name;
    name << outFile << "_p" << setw(2) << setfill('0') << (ii+1);
    ofstream file(name.str().c_str());
    file.setf(ios::scientific,ios::floatfield);
    file.precision(8);
    cout << "Created output file " << name.str() << " for particle ID "
         << partID[ii] << endl;
    for (unsigned int jj = 0; jj < matData[ii].time.size(); ++jj) {
      double time = matData[ii].time[jj];
      int patchIndex = matData[ii].patch[jj];
      int matl = matData[ii].matl[jj];
      long64 pid = matData[ii].id[jj];
      Vector vel = matData[ii].velocity[jj];
      Point pos = matData[ii].position[jj];
      file << time << " " << patchIndex << " " << matl ;
      file << " " << pid;
      file << " " << vel[0] << " " << vel[1] << " " << vel[2]; 
      file << " " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
    }
    file.close();
  }
  */
}

