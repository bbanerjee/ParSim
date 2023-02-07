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
 *  puda.cc: Print out a uintah data archive
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   February 2000
 *
 */

#include <StandAlone/tools/puda/puda.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Math/Matrix3.h>

#include <StandAlone/tools/puda/AA_MMS.h>
#include <StandAlone/tools/puda/ER_MMS.h>
#include <StandAlone/tools/puda/GV_MMS.h>
#include <StandAlone/tools/puda/PIC.h>
#include <StandAlone/tools/puda/POL.h>
#include <StandAlone/tools/puda/UniaxialStrain_MMS.h>
#include <StandAlone/tools/puda/asci.h>
#include <StandAlone/tools/puda/jacquie.h>
#include <StandAlone/tools/puda/jim1.h>
#include <StandAlone/tools/puda/jim2.h>
#include <StandAlone/tools/puda/monica1.h>
#include <StandAlone/tools/puda/monica2.h>
#include <StandAlone/tools/puda/rtdata.h>
#include <StandAlone/tools/puda/util.h>
#include <StandAlone/tools/puda/varsummary.h>

// #include <Core/Containers/Array3.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <cmath>
#include <cstdio>

using namespace std;
using namespace Uintah;

/////////////////////////////////////////////////////////////////
// Pre-declarations:
//
void
printParticleVariable(DataArchive* da,
                      string particleVariable,
                      unsigned long time_step_lower,
                      unsigned long time_step_upper,
                      int mat);

/////////////////////////////////////////////////////////////////

void
usage(const std::string& badarg, const std::string& progname)
{
  if (badarg != "") {
    std::cerr << "Error parsing argument: " << badarg << "\n";
  }
  std::cerr << "Usage: " << progname << " [options] <archive file>\n\n";
  std::cerr << "Valid options are:\n";
  std::cerr << "  -h[elp]\n";
  std::cerr << "  -timesteps\n";
  std::cerr << "  -gridstats\n";
  std::cerr << "  -listvariables\n";
  std::cerr << "  -varsummary\n";
  std::cerr << "  -brief               (Makes varsummary print out a subset of "
               "information.)\n";
  std::cerr << "  -jim1\n";
  std::cerr << "  -jim2\n";
  std::cerr << "  -jacquie              (finds burn rate vs pressure)\n";
  cerr
    << "  -monica1             (Finds the maximum pressure in the domain.)\n";
  std::cerr
    << "  -monica2             (Finds the sum of the cell centered kinetic "
       "energy in the domain.)\n";
  std::cerr << "  -AA_MMS_1            (1D periodic bar MMS)\n";
  std::cerr << "  -AA_MMS_2            (3D Axis aligned MMS)\n";
  std::cerr << "  -GV_MMS              (GeneralizedVortex MMS)\n"; // MMS
  std::cerr << "  -ER_MMS              (Expanding Ring MMS)\n";
  std::cerr << "  -US_MMS              (Uniaxial strain MMS)\n";
  std::cerr << "  -partvar <variable name>\n";
  std::cerr << "  -asci\n";
  std::cerr
    << "  -no_extra_cells      (Excludes extra cells when iterating over "
       "cells.\n";
  std::cerr << "                        Default is to include extra cells.)\n";
  std::cerr << "  -cell_stresses\n";
  std::cerr << "  -rtdata <output directory>\n";
  std::cerr << "  -PTvar\n";
  std::cerr << "  -ptonly              (prints out only the point location)\n";
  std::cerr << "  -patch               (outputs patch id with data)\n";
  std::cerr << "  -material            (outputs material number with data)\n";
  std::cerr << "  -NCvar               (double | float | point | vector)\n";
  std::cerr << "  -CCvar               (double | float | point | vector)\n";
  std::cerr << "  -verbose             (prints status of output)\n";
  std::cerr << "  -timesteplow <int>   (only outputs timestep from int)\n";
  std::cerr << "  -timestephigh <int>  (only outputs timesteps upto int)\n";
  std::cerr << "  -matl,mat <int>      (only outputs data for matl)\n";
  std::cerr
    << "  -pic                 (prints particle ids of all particles  in "
       "cell\n";
  std::cerr << "                        <i> <j> <k> [ints] on the specified "
               "timesteps)\n";
  std::cerr
    << "  -pol                 (prints out average of all particles in a "
       "cell over an\n";
  std::cerr
    << "                       entire line on a line of cells and is called "
       "with:\n";
  std::cerr
    << "                       <axis: [x,y,z]> <ortho1> <ortho2> <average; "
       "default=true>\n";
  std::cerr << "                       <stressSplitting; default=false>\n";
  std::cerr << "                       'ortho1' and 'ortho2' inidicate the "
               "coordinates in the plane\n";
  std::cerr << "                       orthogonal to 'axis'.  'average' tells "
               "whether to average\n";
  std::cerr
    << "                       over all particles in the cell, or just to "
       "use the first\n";
  std::cerr
    << "                       particle encountered.  'stressSplitting' "
       "only takes affect\n";
  std::cerr
    << "                       if the particle variable is p.stress, and "
       "splits the stress\n";
  std::cerr
    << "                       into hydrostatic and deviatoric parts.)\n";
  std::cerr << "*NOTE* to use -PTvar or -NVvar -rtdata must be used\n";
  std::cerr << "*NOTE* ptonly, patch, material, timesteplow, timestephigh "
            << "are used in conjuntion with -PTvar.\n\n";

  std::cerr << "USAGE IS NOT FINISHED\n\n";
  exit(1);
}

void
gridstats(DataArchive* da,
          const bool tslow_set,
          const bool tsup_set,
          unsigned long& time_step_lower,
          unsigned long& time_step_upper)
{
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());

  findTimestep_loopLimits(
    tslow_set, tsup_set, times, time_step_lower, time_step_upper);

  for (unsigned long t = time_step_lower; t <= time_step_upper; t++) {
    double time = times[t];
    std::cout << "__________________________________\n";
    std::cout << "Timestep " << t << ": " << time << "\n";
    GridP grid = da->queryGrid(t);
    grid->performConsistencyCheck();
    grid->printStatistics();

    Vector domainLength;
    grid->getLength(domainLength, "minusExtraCells");
    std::cout << "Domain Length:        " << domainLength << "\n";

    BBox box;
    grid->getInteriorSpatialRange(box);
    std::cout << "\nInterior Spatial Range: " << box << "\n";

    grid->getSpatialRange(box);
    std::cout << "Spatial Range:          " << box << "\n\n";

    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);
      std::cout << "Level: index " << level->getIndex() << ", id "
                << level->getID() << "\n";

      IntVector rr(level->getRefinementRatio());
      std::cout << "       refinement ratio: " << rr << "\n";

      BBox lbox;
      level->getInteriorSpatialRange(lbox);
      std::cout << "       Interior Spatial Range: " << lbox << "\n";

      level->getSpatialRange(lbox);
      std::cout << "       Spatial Range:          " << lbox << "\n\n";

      IntVector lo, hi;
      level->findInteriorCellIndexRange(lo, hi);
      std::cout << "Total Number of Cells:" << hi - lo << "\n";
      std::cout << "dx:                   " << level->dCell() << "\n";

      for (Level::const_patchIterator iter = level->patchesBegin();
           iter != level->patchesEnd();
           iter++) {
        const Patch* patch = *iter;
        std::cout << *patch << "\n";
        std::cout << "\t   BC types: x- " << patch->getBCType(Patch::xminus)
                  << ", x+ " << patch->getBCType(Patch::xplus) << ", y- "
                  << patch->getBCType(Patch::yminus) << ", y+ "
                  << patch->getBCType(Patch::yplus) << ", z- "
                  << patch->getBCType(Patch::zminus) << ", z+ "
                  << patch->getBCType(Patch::zplus) << "\n";
      }
    }
  }
} // end gridstats()

int
main(int argc, char** argv)
{
  if (argc <= 1) {
    // Print out the usage and die
    usage("", argv[0]);
  }

  CommandLineFlags clf;

  int mat   = -1; // not part of clf
  int cellx = -1, celly = -1, cellz = -1;
  char axis             = 'n';
  int ortho1            = -1;
  int ortho2            = -1;
  bool doPOLAverage     = true;
  bool doPOLStressSplit = false;

  // set defaults for cout.
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(8);
  /*
   * Parse arguments
   */
  for (int i = 1; i < argc; i++) {
    string s = argv[i];
    if (s == "-timesteps") {
      clf.do_timesteps = true;
    } else if (s == "-no_extra_cells" || s == "-no_extracells" ||
               s == "-no_ExtraCells") {
      clf.use_extra_cells = false;
    } else if (s == "-gridstats" || s == "-gridStats" || s == "-grid_stats") {
      clf.do_gridstats = true;
    } else if (s == "-listvariables" || s == "-listVariables" ||
               s == "-list_variables") {
      clf.do_listvars = true;
    } else if (s == "-varsummary" || s == "-varSummary" ||
               s == "-var_summary") {
      clf.do_varsummary = true;
    } else if (s == "-brief") {
      clf.be_brief = true;
    } else if (s == "-monica1") {
      clf.do_monica1 = true;
    } else if (s == "-monica2") {
      clf.do_monica2 = true;
    } else if (s == "-jacquie") {
      clf.do_jacquie = true;
    } else if (s == "-jim1") {
      clf.do_jim1 = true;
    } else if (s == "-jim2") {
      clf.do_jim2 = true;
    } else if (s == "-pic") {
      clf.do_PIC = true;

      if (i + 3 >= argc) {
        usage("-pic", argv[0]);
        return 0;
      }

      cellx = strtoul(argv[++i], (char**)nullptr, 10);
      celly = strtoul(argv[++i], (char**)nullptr, 10);
      cellz = strtoul(argv[++i], (char**)nullptr, 10);
    } else if (s == "-pol") {
      if (i + 3 >= argc) {
        usage("-pol", argv[0]);
        return 0;
      }

      axis   = *argv[++i];
      ortho1 = strtoul(argv[++i], (char**)nullptr, 10);
      ortho2 = strtoul(argv[++i], (char**)nullptr, 10);

      clf.do_POL = true;

      // check if optional arguments were found
      if (i + 1 < argc) {
        if (string(argv[i + 1]) == "true") {
          doPOLAverage = true;
          i++;
        } else if (string(argv[i + 1]) == "false") {
          doPOLAverage = false;
          i++;
        }
      }

      if (i + 1 < argc) {
        if (string(argv[i + 1]) == "true") {
          doPOLStressSplit = true;
          i++;
        } else if (string(argv[i + 1]) == "false") {
          doPOLStressSplit = false;
          i++;
        }
      }
    } else if (s == "-AA_MMS_1") {
      clf.do_AA_MMS_1 = true;
    } else if (s == "-AA_MMS_2") {
      clf.do_AA_MMS_2 = true;
    } else if (s == "-GV_MMS") { // MMS
      clf.do_GV_MMS = true;
    } else if (s == "-ER_MMS") {
      clf.do_ER_MMS = true;
    } else if (s == "-US_MMS") {
      clf.do_US_MMS = true;
    } else if (s == "-partvar") {
      clf.do_partvar = true;
      if (i + 1 >= argc) {
        usage("-partvar", argv[0]);
        return 0;
      }
      clf.particleVariable = argv[++i];
      if (clf.particleVariable[0] == '-') {
        usage("-partvar <particle variable name>", argv[0]);
      }
    } else if (s == "-asci") {
      clf.do_asci = true;
    } else if (s == "-cell_stresses") {
      clf.do_cell_stresses = true;
    } else if (s == "-rtdata") {
      clf.do_rtdata = true;
      if (++i < argc) {
        s = argv[i];
        if (s[0] == '-') {
          usage("-rtdata", argv[0]);
        }
        clf.raydatadir = s;
      }
    } else if (s == "-NCvar") {
      if (++i < argc) {
        s = argv[i];
        if (s == "double") {
          clf.do_NCvar_double = true;
        } else if (s == "float") {
          clf.do_NCvar_float = true;
        } else if (s == "point") {
          clf.do_NCvar_point = true;
        } else if (s == "vector") {
          clf.do_NCvar_vector = true;
        } else if (s == "matrix3") {
          clf.do_NCvar_matrix3 = true;
        } else {
          usage("-NCvar", argv[0]);
        }
      } else {
        usage("-NCvar", argv[0]);
      }
    } else if (s == "-CCvar") {
      if (++i < argc) {
        s = argv[i];
        if (s == "double") {
          clf.do_CCvar_double = true;
        } else if (s == "float") {
          clf.do_CCvar_float = true;
        } else if (s == "point") {
          clf.do_CCvar_point = true;
        } else if (s == "vector") {
          clf.do_CCvar_vector = true;
        } else if (s == "matrix3") {
          clf.do_CCvar_matrix3 = true;
        } else {
          usage("-CCvar", argv[0]);
        }
      } else {
        usage("-CCvar", argv[0]);
      }
    } else if (s == "-PTvar") {
      clf.do_PTvar = true;
    } else if (s == "-ptonly") {
      clf.do_PTvar_all = false;
    } else if (s == "-patch") {
      clf.do_patch = true;
    } else if (s == "-material" || s == "-matl" || s == "-mat") {
      if (i + 1 >= argc) {
        usage("-mat", argv[0]);
        return 0;
      }
      clf.matl_jim    = strtoul(argv[++i], (char**)nullptr, 10);
      clf.do_material = true;
      mat             = clf.matl_jim;

    } else if (s == "-verbose") {
      clf.do_verbose = true;
    } else if (s == "-timesteplow" || s == "-timeStepLow" ||
               s == "-timestep_low") {
      if (i + 1 >= argc) {
        usage("-timesteplow", argv[0]);
        return 0;
      }
      clf.time_step_lower = strtoul(argv[++i], (char**)nullptr, 10);
      clf.tslow_set       = true;
    } else if (s == "-timestephigh" || s == "-timeStepHigh" ||
               s == "-timestep_high") {
      if (i + 1 >= argc) {
        usage("-timestephigh", argv[0]);
        return 0;
      }
      clf.time_step_upper = strtoul(argv[++i], (char**)nullptr, 10);
      clf.tsup_set        = true;
    } else if (s == "-timestepinc" || s == "-timestepInc" ||
               s == "-timestep_inc") {
      if (i + 1 >= argc) {
        usage("-timestepinc", argv[0]);
        return 0;
      }
      clf.time_step_inc = strtoul(argv[++i], (char**)nullptr, 10);
    } else if ((s == "-help") || (s == "-h")) {
      usage("", argv[0]);
    } else if (clf.filebase == "") {

      if (argv[i][0] == '-') { // File name can't start with a dash.
        usage(s, argv[0]);
      }
      clf.filebase = argv[i];
    } else {
      usage(s, argv[0]);
    }
  }

  if (clf.filebase == "") {
    std::cerr << "No archive file specified\n";
    usage("", argv[0]);
  }

  try {
    DataArchive* da = scinew DataArchive(clf.filebase);
    //__________________________________
    //  LIST TIMESTEPS
    if (clf.do_timesteps) {
      std::vector<int> index;
      std::vector<double> times;
      da->queryTimesteps(index, times);
      ASSERTEQ(index.size(), times.size());
      std::cout << "There are " << index.size() << " timesteps:\n";

      // Please don't change this.  We need 16
      // significant digits for detailed comparative studies. -Todd
      cout.setf(ios::scientific, ios::floatfield);
      cout.precision(16);

      for (int i = 0; i < (int)index.size(); i++) {
        std::cout << index[i] << ": " << times[i] << "\n";
      }
    }
    //__________________________________

    //    DO GRIDSTATS
    if (clf.do_gridstats) {
      gridstats(da,
                clf.tslow_set,
                clf.tsup_set,
                clf.time_step_lower,
                clf.time_step_upper);
    }

    //__________________________________
    //  LIST VARIABLES
    if (clf.do_listvars) {
      std::vector<std::string> vars;
      std::vector<int> num_matl;
      std::vector<const Uintah::TypeDescription*> types;
      da->queryVariables(vars, num_matl, types);
      std::cout << "There are " << vars.size() << " variables:\n";
      for (int i = 0; i < (int)vars.size(); i++) {
        std::cout << vars[i] << ": " << types[i]->getName() << "\n";
      }
    }

    // Print a particular particle variable
    if (clf.do_partvar && !clf.do_POL) {
      std::vector<int> index;
      std::vector<double> times;
      da->queryTimesteps(index, times);
      ASSERTEQ(index.size(), times.size());
      if (!clf.tslow_set) {
        clf.time_step_lower = 0;
      } else if (clf.time_step_lower >= times.size()) {
        std::cerr << "timesteplow must be between 0 and " << times.size() - 1
                  << "\n";
        abort();
      }
      if (!clf.tsup_set) {
        clf.time_step_upper = times.size() - 1;
      } else if (clf.time_step_upper >= times.size()) {
        std::cerr << "timestephigh must be between 0 and " << times.size() - 1
                  << "\n";
        abort();
      }
      printParticleVariable(da,
                            clf.particleVariable,
                            clf.time_step_lower,
                            clf.time_step_upper,
                            mat);
    }
    //______________________________________________________________________
    //              V A R S U M M A R Y   O P T I O N
    if (clf.do_varsummary) {
      varsummary(da, clf, mat);
    }

    if (clf.do_monica1) {
      monica1(da, clf);
    }

    if (clf.do_monica2) {
      monica2(da, clf);
    }

    if (clf.do_jim1) {
      jim1(da, clf);
    }
    if (clf.do_jacquie) {
      jacquie(da, clf);
    }

    if (clf.do_jim2) {
      jim2(da, clf);
    }

    if (clf.do_PIC) {
      PIC(da, clf, cellx, celly, cellz);
    }

    if (clf.do_POL) {
      POL(da, clf, axis, ortho1, ortho2, doPOLAverage, doPOLStressSplit);
    }

    if (clf.do_AA_MMS_1 || clf.do_AA_MMS_2) {
      AA_MMS(da, clf);
    }

    // MMS
    if (clf.do_GV_MMS) {
      GV_MMS(da, clf);
    }
    if (clf.do_ER_MMS) {
      ER_MMS(da, clf);
    }

    if (clf.do_US_MMS) {
      UniaxialStrain_MMS(da, clf);
    }

    if (clf.do_asci) {
      asci(da,
           clf.tslow_set,
           clf.tsup_set,
           clf.time_step_lower,
           clf.time_step_upper);
    }

    //______________________________________________________________________
    //	       DO CELL STRESSES
    if (clf.do_cell_stresses) {
      std::vector<std::string> vars;
      std::vector<int> num_matl;
      std::vector<const Uintah::TypeDescription*> types;
      da->queryVariables(vars, num_matl, types);
      ASSERTEQ(vars.size(), types.size());

      std::cout << "There are " << vars.size() << " variables:\n";
      std::vector<int> index;
      std::vector<double> times;
      da->queryTimesteps(index, times);
      ASSERTEQ(index.size(), times.size());

      std::cout << "There are " << index.size() << " timesteps:\n";

      findTimestep_loopLimits(clf.tslow_set,
                              clf.tsup_set,
                              times,
                              clf.time_step_lower,
                              clf.time_step_upper);

      // obtain the desired timesteps
      unsigned long t = 0, start_time = 0, stop_time = 0;

      std::cout << "Time Step       Value\n";

      for (t = clf.time_step_lower; t <= clf.time_step_upper; t++) {
        double time = times[t];
        std::cout << "    " << t + 1 << "        " << time << "\n";
      }
      std::cout << "\n";
      if (t != (clf.time_step_lower + 1)) {
        std::cout << "Enter start time-step (1 - " << t << "): ";
        cin >> start_time;
        start_time--;
        std::cout << "Enter stop  time-step (1 - " << t << "): ";
        cin >> stop_time;
        stop_time--;
      } else {
        start_time = t - 1;
        stop_time  = t - 1;
      }
      // end of timestep acquisition

      for (t = start_time; t <= stop_time; t++) {

        double time = times[t];
        std::cout << "time = " << time << "\n";
        GridP grid = da->queryGrid(t);
        for (int v = 0; v < (int)vars.size(); v++) {
          std::string var = vars[v];

          // only dumps out data if it is variable g.stressFS
          if (var == "g.stressFS") {
            const Uintah::TypeDescription* td      = types[v];
            const Uintah::TypeDescription* subtype = td->getSubType();
            std::cout << "\tVariable: " << var << ", type " << td->getName()
                      << "\n";
            for (int l = 0; l < grid->numLevels(); l++) {
              LevelP level = grid->getLevel(l);
              for (Level::const_patchIterator iter = level->patchesBegin();
                   iter != level->patchesEnd();
                   iter++) {
                const Patch* patch = *iter;
                std::cout << "\t\tPatch: " << patch->getID() << "\n";
                ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
                // loop over materials
                for (ConsecutiveRangeSet::iterator matlIter = matls.begin();
                     matlIter != matls.end();
                     matlIter++) {
                  int matl = *matlIter;
                  if (mat != -1 && matl != mat) {
                    continue;
                  }

                  // dumps header and variable info to file
                  std::ostringstream fnum, pnum, matnum;
                  string filename;
                  unsigned long timestepnum = t + 1;
                  fnum << setw(4) << setfill('0') << timestepnum;
                  pnum << setw(4) << setfill('0') << patch->getID();
                  matnum << setw(4) << setfill('0') << matl;
                  string partroot("stress.t");
                  string partextp(".p");
                  string partextm(".m");
                  filename = partroot + fnum.str() + partextp + pnum.str() +
                             partextm + matnum.str();
                  ofstream partfile(filename.c_str());
                  partfile << "# x, y, z, st11, st12, st13, st21, st22, st23, "
                              "st31, st32, st33\n";

                  std::cout << "\t\t\tMaterial: " << matl << "\n";
                  switch (td->getType()) {
                    case Uintah::TypeDescription::Type::NCVariable:
                      switch (subtype->getType()) {
                        case Uintah::TypeDescription::Type::Matrix3: {
                          NCVariable<Matrix3> value;
                          da->query(value, var, matl, patch, t);
                          std::cout << "\t\t\t\t" << td->getName() << " over "
                                    << value.getLowIndex() << " to "
                                    << value.getHighIndex() << "\n";
                          IntVector dx(value.getHighIndex() -
                                       value.getLowIndex());
                          if (dx.x() && dx.y() && dx.z()) {
                            NodeIterator iter = patch->getNodeIterator();
                            for (; !iter.done(); iter++) {
                              partfile << (*iter).x() << " " << (*iter).y()
                                       << " " << (*iter).z() << " "
                                       << (value[*iter])(0, 0) << " "
                                       << (value[*iter])(0, 1) << " "
                                       << (value[*iter])(0, 2) << " "
                                       << (value[*iter])(1, 0) << " "
                                       << (value[*iter])(1, 1) << " "
                                       << (value[*iter])(1, 2) << " "
                                       << (value[*iter])(2, 0) << " "
                                       << (value[*iter])(2, 1) << " "
                                       << (value[*iter])(2, 2) << "\n";
                            }
                          }
                        } break;
                        default:
                          std::cerr << "No Matrix3 Subclass available."
                                    << static_cast<std::underlying_type<
                                         TypeDescription::Type>::type>(
                                         subtype->getType())
                                    << "\n";
                          break;
                      }
                      break;
                    default:
                      cerr
                        << "No NC Variables available."
                        << static_cast<
                             std::underlying_type<TypeDescription::Type>::type>(
                             td->getType())
                        << "\n";
                      break;
                  }
                }
              }
            }
          } else {
            std::cout << "No g.stressFS variables avaliable at time " << t
                      << ".\n";
          }
        }
        if (start_time == stop_time) {
          t++;
        }
      }
    } // end do_cell_stresses

    //______________________________________________________________________
    //              DO RTDATA
    if (clf.do_rtdata) {
      rtdata(da, clf);
    }
  } catch (Exception& e) {
    std::cerr << "Caught exception: " << e.message() << "\n";
    abort();
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    abort();
  }
} // end main()

////////////////////////////////////////////////////////////////////////////
//
// Print ParticleVariable
//
void
printParticleVariable(DataArchive* da,
                      string particleVariable,
                      unsigned long time_step_lower,
                      unsigned long time_step_upper,
                      int mat)
{
  // Check if the particle variable is available
  std::vector<std::string> vars;
  std::vector<int> num_matl;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matl, types);
  ASSERTEQ(vars.size(), types.size());
  bool variableFound = false;
  for (unsigned int v = 0; v < vars.size(); v++) {
    std::string var = vars[v];
    if (var == particleVariable) {
      variableFound = true;
    }
  }
  if (!variableFound) {
    std::cerr << "Variable " << particleVariable << " not found\n";
    exit(1);
  }

  // Now that the variable has been found, get the data for all
  // available time steps // from the data archive
  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  // std::cout << "There are " << index.size() << " timesteps:\n";

  bool useParticleID = true;

  // Loop thru all time steps and store the volume and variable (stress/strain)
  for (unsigned long t = time_step_lower; t <= time_step_upper; t++) {
    double time = times[t];
    // std::cout << "Time = " << time << "\n";
    GridP grid = da->queryGrid(t);

    // Loop thru all the levels
    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);

      // Loop thru all the patches
      Level::const_patchIterator iter = level->patchesBegin();
      int patchIndex                  = 0;
      for (; iter != level->patchesEnd(); iter++) {
        const Patch* patch = *iter;
        ++patchIndex;

        // Loop thru all the variables
        for (int v = 0; v < (int)vars.size(); v++) {
          std::string var                        = vars[v];
          const Uintah::TypeDescription* td      = types[v];
          const Uintah::TypeDescription* subtype = td->getSubType();

          // Check if the variable is a ParticleVariable
          if (td->getType() ==
              Uintah::TypeDescription::Type::ParticleVariable) {

            // loop thru all the materials
            ConsecutiveRangeSet matls = da->queryMaterials(var, patch, t);
            ConsecutiveRangeSet::iterator matlIter = matls.begin();
            for (; matlIter != matls.end(); matlIter++) {
              int matl = *matlIter;
              if (mat != -1 && matl != mat) {
                continue;
              }

              // Find the name of the variable
              if (var == particleVariable) {
                // std::cout << "Material: " << matl << "\n";
                switch (subtype->getType()) {
                  case Uintah::TypeDescription::Type::double_type: {
                    ParticleVariable<double> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        if (useParticleID) {
                          std::cout << " " << pid[*iter];
                        }
                        std::cout << " " << value[*iter] << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::float_type: {
                    ParticleVariable<float> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        if (useParticleID) {
                          std::cout << " " << pid[*iter];
                        }
                        std::cout << " " << value[*iter] << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::int_type: {
                    ParticleVariable<int> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        if (useParticleID) {
                          std::cout << " " << pid[*iter];
                        }
                        std::cout << " " << value[*iter] << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Point: {
                    ParticleVariable<Point> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        if (useParticleID) {
                          std::cout << " " << pid[*iter];
                        }
                        std::cout << " " << value[*iter](0) << " "
                                  << value[*iter](1) << " " << value[*iter](2)
                                  << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Vector: {
                    ParticleVariable<Vector> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        if (useParticleID) {
                          std::cout << time << " " << patchIndex << " " << matl;
                        }
                        std::cout << " " << pid[*iter];
                        std::cout << " " << value[*iter][0] << " "
                                  << value[*iter][1] << " " << value[*iter][2]
                                  << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::Matrix3: {
                    ParticleVariable<Matrix3> value;
                    da->query(value, var, matl, patch, t);
                    ParticleVariable<long64> pid;
                    if (useParticleID) {
                      try {
                        // If particleID wasn't saved, just move on...
                        da->query(pid, "p.particleID", matl, patch, t);
                      } catch (Exception& e) {
                        useParticleID = false;
                      }
                    }
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        if (useParticleID) {
                          std::cout << " " << pid[*iter];
                        }
                        for (int ii = 0; ii < 3; ++ii) {
                          for (int jj = 0; jj < 3; ++jj) {
                            std::cout << " " << value[*iter](ii, jj);
                          }
                        }
                        std::cout << "\n";
                      }
                    }
                  } break;
                  case Uintah::TypeDescription::Type::long64_type: {
                    ParticleVariable<long64> value;
                    da->query(value, var, matl, patch, t);
                    ParticleSubset* pset = value.getParticleSubset();
                    if (pset->numParticles() > 0) {
                      ParticleSubset::iterator iter = pset->begin();
                      for (; iter != pset->end(); iter++) {
                        std::cout << time << " " << patchIndex << " " << matl;
                        std::cout << " " << value[*iter] << "\n";
                      }
                    }
                  } break;
                  default:
                    cerr
                      << "Particle Variable of unknown type: "
                      << static_cast<
                           std::underlying_type<TypeDescription::Type>::type>(
                           subtype->getType())
                      << "\n";
                    break;
                }
              } // end of var compare if
            }   // end of material loop
          }     // end of ParticleVariable if
        }       // end of variable loop
      }         // end of patch loop
    }           // end of level loop
  }             // end of time step loop
} // end printParticleVariable()
