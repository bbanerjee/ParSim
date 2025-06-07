/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/IntVector.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Primes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <cstring>
#include <fstream>
#include <iostream>

using namespace Uintah;

void
usage(char* prog_name);

/*
particleFileSplitter is used in conjunction with particle geometries given as
a series of points in a file.  Given an ASCII file containing, at a minimum x,y
and z for each particle, particleFileSplitter reads an input file which contains
a "file" geometry piece description, and it also reads in the data in the
"points" file associated with that geometry.  particleFileSplitter then creates
a separate file for each patch, and places in that file those points which lie
within that patch.  These individual files may be ASCII or binary.  Also, the
input file may specify that each of the files is to contain additional per
particle data (e.g. temperature), which of course must also be present in the
original points file.  All of this is done to prevent numerous procs from trying
to access the same file, which is hard on file systems.
*/

void
parseArgs(int argc, char* argv[], string& infile, bool& binmode)
{
  binmode = false;
  infile  = argv[argc - 1];

  if (argc < 2) {
    usage(argv[0]);
  }

  if (string(argv[1]) == "-b") {
    binmode = true;
    argc--;
  }
  if (argc > 2) {
    usage(argv[0]);
  }
}

void
usage(char* prog_name)
{
  std::cout << "Usage: " << prog_name << " [-b] [-B] infile \n";
  exit(1);
}

int
main(int argc, char* argv[])
{
  try {
    // Uintah::Parallel::determineIfRunningUnderMPI( argc, argv );
    Uintah::Parallel::initializeManager(argc, argv);

    string infile;
    bool binmode;

    parseArgs(argc, argv, infile, binmode);

    // Get the problem specification
    ProblemSpecP ups = ProblemSpecReader().readInputFile(infile);

    if (!ups) {
      throw ProblemSetupException(
        "Cannot read problem specification", __FILE__, __LINE__);
    }

    if (ups->getNodeName() != "Uintah_specification") {
      throw ProblemSetupException(
        "Input file is not a Uintah specification", __FILE__, __LINE__);
    }

    const ProcessorGroup* world = Uintah::Parallel::getRootProcessorGroup();

    // Setup the initial grid
    GridP grid = scinew Grid();
    IntVector extraCells(0, 0, 0);

    // save and remove the extra cells before the problem setup
    ProblemSpecP g = ups->findBlock("Grid");
    for (ProblemSpecP levelspec = g->findBlock("Level"); levelspec != 0;
         levelspec              = levelspec->findNextBlock("Level")) {
      for (ProblemSpecP box = levelspec->findBlock("Box"); box != 0;
           box              = box->findNextBlock("Box")) {

        ProblemSpecP cells = box->findBlock("extraCells");
        if (cells != 0) {
          box->get("extraCells", extraCells);
          box->removeChild(cells);
        }
      }
    }

    grid->problemSetup(ups, world, false);

    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);

      IntVector low, high;
      level->findCellIndexRange(low, high);
      IntVector diff = high - low;
      long cells     = diff.x() * diff.y() * diff.z();
      if (cells != level->totalCells()) {
        throw ProblemSetupException(
          "particleFileSplitter can only handle square grids",
          __FILE__,
          __LINE__);
      }

      // Parse the geometry from the UPS
      ProblemSpecP mp  = ups->findBlockWithOutAttribute("MaterialProperties");
      ProblemSpecP mpm = mp->findBlock("MPM");
      for (ProblemSpecP child = mpm->findBlock("material"); child != 0;
           child              = child->findNextBlock("material")) {
        for (ProblemSpecP geom_obj_ps = child->findBlock("geom_object");
             geom_obj_ps != 0;
             geom_obj_ps = geom_obj_ps->findNextBlock("geom_object")) {

          for (ProblemSpecP child = geom_obj_ps->findBlock(); child != 0;
               child              = child->findNextBlock()) {
            std::string go_type = child->getNodeName();

            // Read in points from a file
            if (go_type == "file") {
              string f_name, of_name;
              child->require("name", f_name);
              std::ifstream source(f_name.c_str());
              int ncols = 0;
              // count number of vars, and their sizes
              for (ProblemSpecP varblock = child->findBlock("var"); varblock;
                   varblock              = varblock->findNextBlock("var")) {
                string next_var_name("");
                varblock->get(next_var_name);
                if (next_var_name == "p.volume") {
                  ncols += 1;
                } else if (next_var_name == "p.temperature") {
                  ncols += 1;
                } else if (next_var_name == "p.color") {
                  ncols += 1;
                } else if (next_var_name == "p.externalforce") {
                  ncols += 3;
                } else if (next_var_name == "p.fiberdir") {
                  ncols += 3;
                } else if (next_var_name == "p.velocity") {
                  ncols += 3; // gcd adds
                } else if (next_var_name == "p.rvec1") {
                  ncols += 3; // CPDI, CPTI
                } else if (next_var_name == "p.rvec2") {
                  ncols += 3; // CPDI, CPTI
                } else if (next_var_name == "p.rvec3") {
                  ncols += 3; // CPDI, CPTI
                } else {
                  throw ProblemSetupException("Unexpected field variable of '" +
                                                next_var_name + "'",
                                              __FILE__,
                                              __LINE__);
                }
              }

              std::vector<std::vector<std::pair<Point, std::vector<double>>>> points(
                level->numPatches());
              Point min(1e30, 1e30, 1e30), max(-1e30, -1e30, -1e30);

              double x, y, z;
              while (source >> x >> y >> z) {
                std::vector<double> cols;
                for (int ic = 0; ic < ncols; ic++) {
                  double v;
                  source >> v;
                  cols.push_back(v);
                }

                Point pp(x, y, z);
                if (level->containsPointIncludingExtraCells(pp)) {
                  const Patch* currentpatch =
                    level->selectPatchForCellIndex(level->getCellIndex(pp));
                  int pid = currentpatch->getID();
                  min     = Min(pp, min);
                  max     = Max(pp, max);
                  points[pid].push_back(std::pair<Point, std::vector<double>>(pp, cols));
                }
              }

              source.close();
              for (Level::const_patch_iterator iter = level->patchesBegin();
                   iter != level->patchesEnd();
                   iter++) {
                const Patch* patch = *iter;
                int pid            = patch->getID();

                char fnum[5];
                sprintf(fnum, ".%d", pid);
                of_name = f_name + fnum;

                // ADB: change this to always be 128 bytes, so that we can
                // cleanly read the header off a binary file
                std::ofstream dest(of_name.c_str());
                if (binmode) {
                  double x0 = min.x(), y0 = min.y(), z0 = min.z();
                  double x1 = max.x(), y1 = max.y(), z1 = max.z();
                  dest.write((char*)&x0, sizeof(double));
                  dest.write((char*)&y0, sizeof(double));
                  dest.write((char*)&z0, sizeof(double));
                  dest.write((char*)&x1, sizeof(double));
                  dest.write((char*)&y1, sizeof(double));
                  dest.write((char*)&z1, sizeof(double));
                } else {
                  dest << min.x() << " " << min.y() << " " << min.z() << " "
                       << max.x() << " " << max.y() << " " << max.z() << endl;
                }
                for (int I = 0; I < (int)points[pid].size(); I++) {
                  Point p               = points[pid][I].first;
                  std::vector<double> r = points[pid][I].second;

                  // FIXME: should have way of specifying endiness
                  if (binmode) {
                    double x = p.x();
                    double y = p.y();
                    double z = p.z();
                    dest.write((char*)&x, sizeof(double));
                    dest.write((char*)&y, sizeof(double));
                    dest.write((char*)&z, sizeof(double));
                    for (std::vector<double>::const_iterator rit(r.begin());
                         rit != r.end();
                         rit++) {
                      double v = *rit;
                      dest.write((char*)&v, sizeof(double));
                    }
                  } else {
                    dest << p.x() << " " << p.y() << " " << p.z();
                    for (std::vector<double>::const_iterator rit(r.begin());
                         rit != r.end();
                         rit++) {
                      dest << " " << *rit;
                    }
                    dest << endl;
                  }
                }
                dest.close();
              }
            }
          }
        }
      }
    }

  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << '\n';
    if (e.stackTrace()) {
      std::cerr << "Stack trace: " << e.stackTrace() << '\n';
    }
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
  }
}
