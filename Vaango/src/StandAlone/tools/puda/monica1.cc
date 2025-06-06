/*
 * The MIT License
 *
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
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <Core/DataArchive/DataArchive.h>
#include <StandAlone/tools/puda/monica1.h>
#include <StandAlone/tools/puda/util.h>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace Uintah;

using namespace std;

////////////////////////////////////////////////////////////////////////
//              M O N I C A 1   O P T I O N
//   This function will output the maximum pressure achieved during the
//   simulation as well as the maximum pressure at each timestep.  Time-
//   step pressures are output to a file called "maxPressures.dat".

void
Uintah::monica1(DataArchive* da, CommandLineFlags& clf)
{
  std::vector<std::string> vars;
  std::vector<int> num_matl;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matl, types);
  ASSERTEQ(vars.size(), types.size());
  std::cout << "There are " << vars.size() << " variables:\n";
  for (int i = 0; i < (int)vars.size(); i++) {
    std::cout << vars[i] << ": " << types[i]->getName() << std::endl;
  }

  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "There are " << index.size() << " timesteps:\n";
  for (int i = 0; i < (int)index.size(); i++) {
    std::cout << index[i] << ": " << times[i] << std::endl;
  }

  findTimestep_loopLimits(clf.tslow_set,
                          clf.tsup_set,
                          times,
                          clf.time_step_lower,
                          clf.time_step_upper);

  std::ostringstream fnum;
  string filename("maxPressures.dat");
  ofstream outfile(filename.c_str());

  double maxPressure = -9999999; // the max pressure at any time

  for (unsigned long t = clf.time_step_lower; t <= clf.time_step_upper;
       t += clf.time_step_inc) {
    double time = times[t];
    std::cout << "time = " << time << std::endl;
    GridP grid = da->queryGrid(t);

    double pressure = -9999999.0; // the max pressure during the timestep
    LevelP level    = grid->getLevel(grid->numLevels() - 1);
    std::cout << "Level: " << grid->numLevels() - 1 << endl;
    for (Level::const_patch_iterator iter = level->patchesBegin();
         iter != level->patchesEnd();
         iter++) {
      const Patch* patch = *iter;
      int matl           = clf.matl_jim; // material number

      CCVariable<double> press_CC;
      // get all the pressures from the patch
      da->query(press_CC, "press_CC", matl, patch, t);

      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter; // get teh coordinates of the cell

        if (press_CC[c] > pressure) {
          pressure = press_CC[c];
        }

        if (press_CC[c] > maxPressure) {
          maxPressure = press_CC[c];
        }

      } // for cells
    }   // for patches

    std::cout << "Max pressure for timestep was:\t" << pressure << std::endl;

    outfile.precision(15);
    outfile << t << " " << pressure << std::endl;
  }
  std::cout << "Max pressure overall was:\t" << maxPressure << std::endl;
} // end jim2()
