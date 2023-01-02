/*
 * The MIT License
 *
 * Copyright (c) 2014-2023 Biswajit Banerjee
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
 *  extractPScalar.cc: Print out a uintah data archive for particle scalar variable data
 *
 *  Written by:
 *   Biswajit Banerjee
 *   June 2016
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
//#include <Core/Containers/Array3.h>
#include <Core/Geometry/Point.h>
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
#include <algorithm>
#include <map>
#include <chrono>


using namespace Uintah;

typedef struct{
  std::vector<int>            patch;
  std::vector<int>            matl;
  std::vector<Uintah::long64> id;
  std::vector<double>         variable;
  std::vector<double>         time;
  std::vector<Uintah::Point>  position;
} MaterialData;

void usage(const std::string& badarg, const std::string& progname);

void printScalarVariable(DataArchive* da, 
                         std::string partVar,
                         int matID,
                         std::vector<long64>& partID,
                         std::string outFile,
                         bool timeFiles);

int main(int argc, char** argv)
{
  std::string partVar;
  int matID = 0;
  std::string particleVariable;
  std::string partIDFile;
  std::string udaDir;
  std::string outFile;
  bool timeFiles = false;

  // set defaults for cout
  std::cout.setf(std::ios::scientific,std::ios::floatfield);
  std::cout.precision(8);
  /*
   * Parse arguments
   */
  for(int i=1;i<argc;i++){
    string s=argv[i];
    if (s == "-m") {
      string id = argv[++i];
      if (id[0] == '-') 
        usage("-m <material id>", argv[0]);
      matID = atoi(argv[i]);
    } else if (s == "-partvar") {
      particleVariable = argv[++i];
      if (particleVariable[0] == '-')
        usage("-partvar <particle variable name>", argv[0]);
    } else if(s == "-p"){
      partIDFile = argv[++i];
      if (partIDFile[0] == '-') 
        usage("-p <particle id file>", argv[0]);
    } else if(s == "-uda"){
      udaDir = argv[++i];
      if (udaDir[0] == '-') 
        usage("-uda <archive file>", argv[0]);
    } else if(s == "-o"){
      outFile = argv[++i];
      if (outFile[0] == '-') 
        usage("-o <output file>", argv[0]);
    } else if (s == "-timefiles") {
      timeFiles = true;
    }
  }
  if (argc < 11) usage( "", argv[0] );

  std::cout << "Particle Variable to be extracted = " << particleVariable << "\n";
  std::cout << "Material ID to be extracted = " << matID << endl;

  // Read the particle ID file
  std::cout << "Particle ID File to be read = " << partIDFile << endl;
  vector<long64> partID;
  std::ifstream pidFile(partIDFile.c_str());
  if (!pidFile.is_open()) {
    std::cerr << "Particle ID File " << partIDFile << " not found \n";
    exit(1);
  }
  do {
    long64 id = 0;
    double t = 0.0, x = 0.0, y = 0.0, z = 0.0;
    int patch = 0, mat = 0;
    pidFile >> t >> patch >> mat >> id >> x >> y >> z;
    //pidFile >> id;
    if (id > 0) {
      partID.push_back(id);
    }
  } while (!pidFile.eof());
  
  std::cout << "  Number of Particle IDs = " << partID.size() << endl;
  for (unsigned int ii = 0; ii < partID.size() ; ++ii) {
    std::cout << "    p"<< (ii+1) << " = " << partID[ii] << endl;
  }

  std::cout << "Output file name = " << outFile << endl;
  std::cout << "UDA directory to be read = " << udaDir << endl;
  try {
    DataArchive* da = scinew DataArchive(udaDir);
    
    // Print a particular particle variable
    printScalarVariable(da, particleVariable, matID, partID, outFile, timeFiles);
  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << endl;
    abort();
  } catch(...){
    std::cerr << "Caught unknown exception\n";
    abort();
  }

  exit(0);
}

void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") std::cerr << "Error parsing argument: " << badarg << endl;
  std::cerr << "Usage: " << progname 
       << " -partvar <scalar particle variable>"
       << " -m <material id> "
       << " -p <particle id file>"
       << " -uda <archive file>"
       << " -o <output file>"
       << " -timesfiles (optional)\n\n";
  exit(1);
}

////////////////////////////////////////////////////////////////////////////
//
// Print a particle variable
//
////////////////////////////////////////////////////////////////////////////
void printScalarVariable(DataArchive* da, 
                         std::string partVar,
                         int matID,
                         std::vector<long64>& partID,
                         std::string outFile,
                         bool timeFiles)
{

  // Check if the particle variable is available
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());
  bool variableFound = false;
  for (auto var : vars) {
    if (var == partVar) variableFound = true;
  }
  if (!variableFound) {
    std::cerr << "Variable " << partVar << " not found\n"; 
    exit(1);
  }

  // Create arrays of material data for each particle ID
  std::map<long64, MaterialData> matData;

  // Now that the variable has been found, get the data for all 
  // available time steps from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";
      
  // Loop thru all the variables 
  for (auto var : vars) {

    // Find the name of the variable
    //std::cout << "Variable = " << var << std::endl;
    if (var == partVar) {

      // Loop thru all time steps 
      for (auto t = 0ul; t < times.size(); t++) {
        double time = times[t];

        //std::cerr << "t = " << time ;
        auto start = std::chrono::system_clock::now();
        GridP grid = da->queryGrid(t);

        unsigned int numFound = 0;

        // Loop thru all the levels
        for (int l = 0; l < grid->numLevels(); l++) {

          //std::cerr << " level = " << l ;
          LevelP level = grid->getLevel(l);
          int patchIndex = 0;

          // Loop thru all the patches
          for(auto iter = level->patchesBegin(); iter != level->patchesEnd(); iter++){

            const Patch* patch = *iter;
            ++patchIndex; 
            //std::cerr << " patch = " << patchIndex ;

            ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);

            // loop thru all the materials
            for (auto matl : matls) {

              if (matl != matID) continue;
              //std::cerr << " matl = " << matl ;

              if (numFound == partID.size()) break;

              ParticleVariable<Point>  position;
              ParticleVariable<long64> pid;
              da->query(position, "p.x",          matl, patch, t);
              da->query(pid,      "p.particleID", matl, patch, t);

              ParticleVariable<double> value;
              da->query(value, var, matl, patch, t);

              ParticleSubset* pset = position.getParticleSubset();
              //std::cerr << " particles = " << pset->numParticles() ;
              if (pset->numParticles() > 0) {

                for (auto particleID : partID) {
                  for(auto particle : *pset) {
                    if (particleID != pid[particle]) {
                      continue;
                    }
                    (matData[particleID].position).push_back(position[particle]);
                    (matData[particleID].variable).push_back(value[particle]);
                    (matData[particleID].id).push_back(pid[particle]);
                    (matData[particleID].time).push_back(time);
                    (matData[particleID].patch).push_back(patchIndex);
                    (matData[particleID].matl).push_back(matl);
                    ++numFound;
                  }
                }
              }
            } // end of material loop
          } // end of patch loop
        } // end of level loop
        auto end = std::chrono::system_clock::now();
        auto timetaken = end - start;
        std::cerr << " CPU Time = " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(timetaken).count() 
                  << " ms" << " found " << numFound << endl;
      } // end of time step loop
    } // end of var compare if
  } // end of variable loop

  
  if (timeFiles) {
    // Create one output file for all of the timesteps
    std::ofstream file(outFile);
    file.setf(std::ios::scientific, std::ios::floatfield);
    file.precision(8);
    std::cout << "Created output file " << outFile << endl;
    //std::cout << "Data size = " << matData.size() << std::endl;
    for (auto jj = 0ul; jj < times.size(); ++jj) {
      double time = times[jj];
      int numFound = 0;
      for (auto particleID : partID) {
        
        if (matData.find(particleID) != matData.end()) {
          if (matData[particleID].time.size() <= times.size()) {
            numFound++;
            auto patchIndex = matData[particleID].patch[jj];
            auto matl = matData[particleID].matl[jj];
            auto pid = matData[particleID].id[jj];
            auto var = matData[particleID].variable[jj];
            auto pos = matData[particleID].position[jj];
            file << time << " " << patchIndex << " " << matl ;
            file << " " << pid;
            file << " " << var;
            file << " " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
          } else {
            std::cerr << "**WARNING** Data not found for all timesteps for particle "
                      << particleID << std::endl;
          }
        }
      }
      //std::cout << "num found = " << numFound << std::endl;
    }
    file.close();
    std::cout << "Closed output file " << outFile << endl;
  } else {
    // Create output files for each of the particle IDs
    for (unsigned int ii = 0; ii < partID.size() ; ++ii) {
      auto particleID = partID[ii];
      if (matData.find(particleID) != matData.end()) {
        ostringstream name;
        name << outFile << "_p" << std::setw(2) << std::setfill('0') << (ii+1);
        std::ofstream file(name.str());
        file.setf(std::ios::scientific,std::ios::floatfield);
        file.precision(8);
        std::cout << "Created output file " << name.str() << " for particle ID "
             << particleID << endl;
        for (unsigned int jj = 0; jj < matData[particleID].time.size(); ++jj) {
          auto time = matData[particleID].time[jj];
          auto patchIndex = matData[particleID].patch[jj];
          auto matl = matData[particleID].matl[jj];
          auto pid = matData[particleID].id[jj];
          auto var = matData[particleID].variable[jj];
          auto pos = matData[particleID].position[jj];
          file << time << " " << patchIndex << " " << matl ;
          file << " " << pid;
          file << " " << var;
          file << " " << pos.x() << " " << pos.y() << " " << pos.z() << endl;
        }
        file.close();
        std::cout << "Closed output file " << outFile << endl;
      }
    }
  }
}

