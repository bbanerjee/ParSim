/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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
 *  selectpart.cc: Get the range of particle ids in a given selection box
 *
 *  Written by:
 *   Biswajit Banerjee
 *   July 2005
 *
 *  Original code works only for the first timestep, and all the materials.
 *  Added -mat, -timesteplow, -timestephigh 
 *     options by Jonah Lee, December 2008
 *
 *  Update:
 *  July 2016:  Allows points and  axis aligned lines, planes, and boxes
 *
 */

#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/SymmMatrix3.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/Variables/ShareAssignParticleVariable.h>
#include <Core/Math/MinMax.h>
#include <Core/Geometry/Point.h>
#include <Core/Grid/Box.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Geometry/Vector.h>
#include <Core/OS/Dir.h>
#include <Core/Containers/Array3.h>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace SCIRun;
using namespace std;
using namespace Uintah;

// declarations
void usage(const std::string& badarg, const std::string& progname);

void printParticleID(DataArchive* da, int mat, 
                     bool pointExists, bool lineExists, bool planeExists, bool boxExists,
                     const Point& pt0, const Point& pt1, 
                     const Point& pt2, const Point& pt3, 
                     bool tslow_set, bool tsup_set, 
                     unsigned long& time_step_lower, unsigned long& time_step_upper);

//borrowed from puda.h
void findTimestep_loopLimits(const bool tslow_set, 
                             const bool tsup_set,
                             const std::vector<double>& times,
                             unsigned long & time_step_lower,
                             unsigned long & time_step_upper);

// From: https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
namespace Vaango {
  class InputParser{

    public:

      InputParser (const int &argc, const char **argv) {
        for (int i=1; i < argc; ++i) {
          d_tokens.push_back(std::string(argv[i]));
        }
        d_currentIter = d_tokens.begin();
      }

      /// @author iain
      const std::string getCmdOption(const std::string &option) {
        auto itr =  std::find(d_tokens.begin(), d_tokens.end(), option);
        if (itr != d_tokens.end() && ++itr != d_tokens.end()) {
          d_currentIter = itr;
          return *itr;
        }
        d_currentIter = d_tokens.begin();
        return std::string("");
      }

      const std::string getNextArg() {
        if (d_currentIter != d_tokens.end() && ++d_currentIter != d_tokens.end()) {
          return *d_currentIter;
        }
        return std::string("");
      }

      /// @author iain
      bool cmdOptionExists(const std::string &option) const{
        return 
          std::find(d_tokens.begin(), d_tokens.end(), option) != d_tokens.end();
      }

    private:
      std::vector <std::string> d_tokens;
      std::vector <std::string>::iterator d_currentIter;

    friend std::ostream& operator<<(std::ostream& out, const InputParser& input) {
      for (auto token: input.d_tokens)  {
        out << token << " ";
      }
      out << endl;
      return out;
    }
  };

}

// Main
int main(const int argc, const char** argv)
{
  /*
   * Default values
   */
  double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  double x1 = 0.5, y1 = 0.5, z1 = 0.5;
  double x2 = 1.0, y2 = 1.0, z2 = 1.0;
  double x3 = 1.0, y3 = 1.0, z3 = 1.0;
  Uintah::Point pt0 = {x0, y0, z0};
  Uintah::Point pt1 = {x1, y1, z1};
  Uintah::Point pt2 = {x2, y2, z2};
  Uintah::Point pt3 = {x3, y3, z3};
  string filebase;

  int mat = -1;
  unsigned long time_step_lower = 0;
  unsigned long time_step_upper = 1;
  unsigned long time_step_inc = 1;
  bool tslow_set = false;
  bool tsup_set = false;

  // set defaults for cout
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(8);

  /*-------------------------------------------------------------------------------
   * Parse arguments
   */

  // Create input parser
  Vaango::InputParser input(argc, argv);
  cerr << input;

  // Check if only one of point, line, box have been given in the input line
  bool pointExists = input.cmdOptionExists("-point");
  bool lineExists = input.cmdOptionExists("-line");
  bool planeExists = input.cmdOptionExists("-plane");
  bool boxExists = input.cmdOptionExists("-box");
  if (!pointExists && !lineExists && !boxExists && !planeExists) {
    std::cerr << "**ERROR** No input search region provided: \n"; 
    std::cerr << "\t Input arguments provided: " <<  input;
    std::cerr << "\t Need one of the following: " << std::endl;
    std::cerr << "\t\t -point x y z " << std::endl;
    std::cerr << "\t\t -line x0 y0 z0 x1 y1 z1 " << std::endl;
    std::cerr << "\t\t -plane x0 y0 z0 x1 y1 z1 x2 y2 z2 " << std::endl;
    std::cerr << "\t\t -box x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3 " << std::endl;
    usage("", argv[0]);
    exit(1);
  }

  if ((pointExists && lineExists) || (pointExists && boxExists) ||
      (pointExists && planeExists) || (lineExists && boxExists) ||
      (lineExists && planeExists) || (boxExists && planeExists)) {
    std::cerr << "**ERROR** More than one input search region provided \n";
    std::cerr << "\t Input arguments provided: " <<  input;
    std::cerr << "\t Need only one of the following: " << std::endl;
    std::cerr << "\t\t -point x y z " << std::endl;
    std::cerr << "\t\t -line x0 y0 z0 x1 y1 z1 " << std::endl;
    std::cerr << "\t\t -plane x0 y0 z0 x1 y1 z1 x2 y2 z2 " << std::endl;
    std::cerr << "\t\t -box x0 y0 z0 x1 y1 z1 x2 y2 z2 x3 y3 z3 " << std::endl;
    usage("", argv[0]);
    exit(1);
  }

  // If a "-point" has been given as input then the next three arguments
  // are the coordinates of the point
  if (pointExists) {
    try {
      x0 = std::stod(input.getCmdOption("-point"));
      y0 = std::stod(input.getNextArg());
      z0 = std::stod(input.getNextArg());
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** Point coordinates must be numeric values." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    pt0.x(x0); pt0.y(y0); pt0.z(z0); 
    std::cerr << "Search point = " << pt0 << std::endl;
  }
  
  // If a "-line" has been given as input then the next six arguments
  // are the coordinates of the start point and end point
  if (lineExists) {
    try {
      x0 = std::stod(input.getCmdOption("-line"));
      y0 = std::stod(input.getNextArg());
      z0 = std::stod(input.getNextArg());
      x1 = std::stod(input.getNextArg());
      y1 = std::stod(input.getNextArg());
      z1 = std::stod(input.getNextArg());
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** Line coordinates must be numeric values." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    pt0.x(x0); pt0.y(y0); pt0.z(z0); 
    pt1.x(x1); pt1.y(y1); pt1.z(z1); 
    if ((pt1 - pt0).length2() < 1.0e-16) {
      std::cerr << "**ERROR** Line has zero length." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    std::cerr << "Search line: Start = " << pt0 <<  " End = " << pt1 << std::endl;
  }

  // If a "-plane" has been given as input then the next six arguments
  // are the coordinates of the start point and two end points
  if (planeExists) {
    try {
      x0 = std::stod(input.getCmdOption("-plane"));
      y0 = std::stod(input.getNextArg());
      z0 = std::stod(input.getNextArg());
      x1 = std::stod(input.getNextArg());
      y1 = std::stod(input.getNextArg());
      z1 = std::stod(input.getNextArg());
      x2 = std::stod(input.getNextArg());
      y2 = std::stod(input.getNextArg());
      z2 = std::stod(input.getNextArg());
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** Plane coordinates must be numeric values." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    pt0.x(x0); pt0.y(y0); pt0.z(z0); 
    pt1.x(x1); pt1.y(y1); pt1.z(z1); 
    pt2.x(x2); pt2.y(y2); pt2.z(z2); 
    Vector crossprod = Cross((pt2 - pt0), (pt1 - pt0));
    if (crossprod.length2() < 1.0e-16) {
      std::cerr << "**ERROR** Plane has zero area." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    std::cerr << "Search plane: Start = " << pt0 << std::endl
              << " End 1 = " << pt1 << " End 2 = " << pt2 << std::endl;
  }

  // If a "-box" has been given as input then the next six arguments
  // are the coordinates of the start point and two end points
  if (boxExists) {
    try {
      x0 = std::stod(input.getCmdOption("-box"));
      y0 = std::stod(input.getNextArg());
      z0 = std::stod(input.getNextArg());
      x1 = std::stod(input.getNextArg());
      y1 = std::stod(input.getNextArg());
      z1 = std::stod(input.getNextArg());
      x2 = std::stod(input.getNextArg());
      y2 = std::stod(input.getNextArg());
      z2 = std::stod(input.getNextArg());
      x3 = std::stod(input.getNextArg());
      y3 = std::stod(input.getNextArg());
      z3 = std::stod(input.getNextArg());
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** Box coordinates must be numeric values." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    pt0.x(x0); pt0.y(y0); pt0.z(z0); 
    pt1.x(x1); pt1.y(y1); pt1.z(z1); 
    pt2.x(x2); pt2.y(y2); pt2.z(z2); 
    pt3.x(x3); pt3.y(y3); pt3.z(z3); 
    Vector crossprod = Cross((pt2 - pt0), (pt1 - pt0));
    if (std::abs(Dot(crossprod, (pt3 - pt0))) < 1.0e-16) {
      std::cerr << "**ERROR** Box has zero volume." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    std::cerr << "Search Box: Start = " << pt0 << std::endl
              << " End 1 = " << pt1 << " End 2 = " << pt2 << std::endl
              << " End 3 = " << pt3 << std::endl;
  }
 
  // Get the material id
  try {
    mat = std::stoi(input.getCmdOption("-mat"));
  } catch (const std::invalid_argument& is) {
    std::cerr << "**ERROR** A material id code must be provided." << std::endl;
    std::cerr << "\t Input arguments provided: " <<  input;
    usage("", argv[0]);
    exit(1);
  }

  // Check the low timestep
  bool timesteplowExists = input.cmdOptionExists("-timesteplow");
  if (timesteplowExists) {
    try {
      time_step_lower = std::stoul(input.getCmdOption("-timesteplow"));
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** A invalid value of timesteplow has been provided." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    tslow_set = true; 
  } 
  
  // Check the high timestep
  bool timestephighExists = input.cmdOptionExists("-timestephigh");
  if (timestephighExists) {
    try {
      time_step_upper = std::stoul(input.getCmdOption("-timestephigh"));
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** A invalid value of timestephigh has been provided." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
    tsup_set = true; 
  } 

  // Check the timestep increment
  bool timestepincExists = input.cmdOptionExists("-timestepinc");
  if (timestepincExists) {
    try {
      time_step_inc = std::stoul(input.getCmdOption("-timestepinc"));
    } catch (const std::invalid_argument& is) {
      std::cerr << "**ERROR** A invalid value of timestepinc has been provided." << std::endl;
      std::cerr << "\t Input arguments provided: " <<  input;
      usage("", argv[0]);
      exit(1);
    }
  } 

  // Get the file name
  filebase = input.getCmdOption("-uda");
  if (filebase == "") {
    std::cerr << "**ERROR** No input uda file name provided." << std::endl;
    std::cerr << "\t Input arguments provided: " <<  input;
    usage("", argv[0]);
    exit(1);
  }
  
  /*-------------------------------------------------------------------------------*/
  try {
    DataArchive* da = scinew DataArchive(filebase);

    // Get the particle IDs
    printParticleID(da, mat, pointExists, lineExists, planeExists, boxExists,
                    pt0, pt1, pt2, pt3, tslow_set, tsup_set, time_step_lower, 
                    time_step_upper);

  } catch (Exception& e) {
    cerr << "Caught exception: " << e.message() << endl;
    exit(1);
  } catch(...){
    cerr << "Caught unknown exception\n";
    exit(1);
  }
}

void usage(const std::string& badarg, const std::string& progname)
{
  if(badarg != "") cerr << "Error parsing argument: " << badarg << endl;

  cerr << "Usage: " << progname << "[options] <archive file>\n\n";
  cerr << "Valid options are:\n";
  cerr << "  exactly one of (required): \n";
  cerr << "    -point <x0> <y0> <z0>\n";
  cerr << "    -line <x0> <y0> <z0> <x1> <y1> <z1>\n";
  cerr << "    -plane <x0> <y0> <z0> <x1> <y1> <z1> <x2> <y2> <z2>\n";
  cerr << "    -box <x0> <y0> <z0> <x1> <y1> <z1> <x2> <y2> <z2> <x3> <y3> <z3>\n";
  cerr << "  -mat <material id>\n";
  cerr << "  -timesteplow <int>  (only outputs timestep from int)\n";
  cerr << "  -timestephigh <int> (only outputs timesteps upto int)\n";
  cerr << "  -uda <filename>\n";
  exit(1);
}

void
findTimestep_loopLimits(bool tslow_set, 
                        bool tsup_set,
                        const vector<double>& times,
                        unsigned long& time_step_lower,
                        unsigned long& time_step_upper)
{
  if( !tslow_set ) {
    time_step_lower = 0;
  }
  else if( time_step_lower >= times.size() ) {
    cerr << "timesteplow must be between 0 and " << times.size()-1 << "\n";
    exit(1);
  }
  if( !tsup_set || time_step_upper > (times.size() - 1)) {
    time_step_upper = times.size() - 1;
  }
  else if( time_step_upper >= times.size() ) {
    cerr << "timestephigh must be between 0 and " << times.size()-1 << "\n";
    exit(1);
  }
}


////////////////////////////////////////////////////////////////////////////
//
// Print particle IDs
//
////////////////////////////////////////////////////////////////////////////
void printParticleID(DataArchive* da, int mat, 
                     bool point, bool line, bool plane, bool box,
                     const Point& pt0, const Point& pt1,
                     const Point& pt2, const Point& pt3,
                     bool tslow_set, bool tsup_set,
                     unsigned long& time_step_lower,
                     unsigned long& time_step_upper)
{
  // Check if the particle variable p.particleId and p.x are available
  std::vector<std::string> vars;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());

  bool variableFound = false;
  for (auto v : vars) {
    if (v == "p.particleID" || v == "p.x") variableFound = true;
  }

  if (!variableFound) {
    cerr << "p.particleID or p.x not found\n"; 
    exit(1);
  }

  // Now that the variable has been found, get the data for the
  // desired timesteps from the archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());

  findTimestep_loopLimits(tslow_set, tsup_set, times, time_step_lower, 
                          time_step_upper); 

  // Loop thru all time steps
  for (auto t = time_step_lower; t <= time_step_upper; t++) {

    double time = times[t];
    GridP grid = da->queryGrid(t);

    // Loop thru all the levels
    for (int l = 0; l < grid->numLevels(); l++) {
     LevelP level = grid->getLevel(l);

     // Loop thru all the patches
     int patchIndex = 0;
     for(auto iter = level->patchesBegin(); iter != level->patchesEnd(); iter++){
       const Patch* patch = *iter;
       ++patchIndex; 

       // Get the cell size
       Vector cellSize = patch->dCell();
       std::cerr << "cellSize = " << cellSize << std::endl;

       // Setup coarse selection box (0.251 is needed because Grid/Box.h is used)
       // Default is that a point will be used for selection
       Box selectionBox(pt0 - 0.25*cellSize, pt0 + 0.251*cellSize);
       if (line) {
         Box coarseBox(pt0 - 0.25*cellSize, pt1 + 0.251*cellSize);
         selectionBox = coarseBox;
       } else if (plane) {
         Vector v1 = pt1 - pt0;         
         Vector v2 = pt2 - pt0;         
         Point pt4 = (v1 + v2).asPoint();
         double minx = std::min(std::min(pt0.x(), pt1.x()),
                                std::min(pt2.x(), pt4.x()));
         double maxx = std::max(std::max(pt0.x(), pt1.x()),
                                std::max(pt2.x(), pt4.x()));
         double miny = std::min(std::min(pt0.y(), pt1.y()),
                                std::min(pt2.y(), pt4.y()));
         double maxy = std::max(std::max(pt0.y(), pt1.y()),
                                std::max(pt2.y(), pt4.y()));
         double minz = std::min(std::min(pt0.z(), pt1.z()),
                                std::min(pt2.z(), pt4.z()));
         double maxz = std::max(std::max(pt0.z(), pt1.z()),
                                std::max(pt2.z(), pt4.z()));
         Point low(minx, miny, minz), high(maxx, maxy, maxz);
         Box coarseBox(low - 0.25*cellSize, high + 0.251*cellSize);
         selectionBox = coarseBox;
       } else if (box) {
         Vector v1 = pt1 - pt0;         
         Vector v2 = pt2 - pt0;         
         Vector v3 = pt3 - pt0;         
         Point pt4 = (v1 + v2 + v3).asPoint();
         Point pt5 = (v1 + v2).asPoint();
         Point pt6 = (v1 + v3).asPoint();
         Point pt7 = (v2 + v3).asPoint();
         double minx = std::min(std::min(std::min(pt0.x(), pt1.x()),
                                         std::min(pt2.x(), pt3.x())),
                                std::min(std::min(pt4.x(), pt5.x()),
                                         std::min(pt6.x(), pt7.x())));
         double maxx = std::max(std::max(std::max(pt0.x(), pt1.x()),
                                         std::max(pt2.x(), pt3.x())),
                                std::max(std::max(pt4.x(), pt5.x()),
                                         std::max(pt6.x(), pt7.x())));
         double miny = std::min(std::min(std::min(pt0.y(), pt1.y()),
                                         std::min(pt2.y(), pt3.y())),
                                std::min(std::min(pt4.y(), pt5.y()),
                                         std::min(pt6.y(), pt7.y())));
         double maxy = std::max(std::max(std::max(pt0.y(), pt1.y()),
                                         std::max(pt2.y(), pt3.y())),
                                std::max(std::max(pt4.y(), pt5.y()),
                                         std::max(pt6.y(), pt7.y())));
         double minz = std::min(std::min(std::min(pt0.z(), pt1.z()),
                                         std::min(pt2.z(), pt3.z())),
                                std::min(std::min(pt4.z(), pt5.z()),
                                         std::min(pt6.z(), pt7.z())));
         double maxz = std::max(std::max(std::max(pt0.z(), pt1.z()),
                                         std::max(pt2.z(), pt3.y())),
                                std::max(std::max(pt4.z(), pt5.z()),
                                         std::max(pt6.z(), pt7.z())));
         Point low(minx, miny, minz), high(maxx, maxy, maxz);
         Box coarseBox(low - 0.25*cellSize, high + 0.251*cellSize);
         selectionBox = coarseBox;
       } 
       cerr << "Selection box = " << selectionBox << endl;
    
       // loop thru all the materials
       auto matls = da->queryMaterials("p.x", patch, t);
       //for(auto matlIter = matls.begin(); matlIter != matls.end(); matlIter++){
       for (auto matl : matls) {
         //int matl = *matlIter;
         if (mat != -1 && matl != mat) continue;

         ParticleVariable<Point> position;
         ParticleVariable<long64> pid;
         da->query(position, "p.x",          matl, patch, t);
         da->query(pid,      "p.particleID", matl, patch, t);

         ParticleSubset* pset = position.getParticleSubset();
         if (pset->numParticles() > 0) {
          //for(auto iter = pset->begin(); iter != pset->end(); iter++){
          // auto pidx = pidx;
          for (auto pidx : *pset) {
	     if (selectionBox.contains(position[pidx])) {
               cout << time << " " << patchIndex << " " << matl ;
               cout << " " << pid[pidx];
               cout << " " << position[pidx](0) 
                    << " " << position[pidx](1)
                    << " " << position[pidx](2) << endl;
	     }
          }
        }
      }
    } // end of patch loop
  } // end of level loop
 } // end of time loop
}

