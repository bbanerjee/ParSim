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
#include <StandAlone/tools/puda/monica2.h>
#include <StandAlone/tools/puda/util.h>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace Uintah;

using namespace std;

////////////////////////////////////////////////////////////////////////
//              M O N I C A 2   O P T I O N
// Reads in cell centered mass and velocity, and compiles mean velocity
// magnitude and total kinetic energy for all cells of a specified material
// within the range of timesteps specified

void
Uintah::monica2(DataArchive* da, CommandLineFlags& clf)
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

  // Output File
  std::ostringstream fnum;
  string filename("time_meanvel_KE.dat");
  ofstream outfile(filename.c_str());

  // Timestep Loop
  for (unsigned long t = clf.time_step_lower; t <= clf.time_step_upper;
       t += clf.time_step_inc) {
    double time = times[t];
    std::cout << "time = " << time << std::endl;
    GridP grid = da->queryGrid(t);

    Vector mean_vel(0., 0., 0.);
    double KE         = 0.;
    double total_mass = 0.;
    LevelP level      = grid->getLevel(grid->numLevels() - 1);
    std::cout << "Level: " << grid->numLevels() - 1 << endl;
    for (Level::const_patch_iterator iter = level->patchesBegin();
         iter != level->patchesEnd();
         iter++) {
      const Patch* patch = *iter;
      int matl           = clf.matl_jim;
      //____________________________________________
      //   C E L L C E N T E R E D   V A R I A B L E
      CCVariable<Vector> value_velocity;
      CCVariable<double> value_density;
      CCVariable<double> value_vf;

      da->query(value_velocity, "vel_CC", matl, patch, t);
      da->query(value_density, "rho_CC", matl, patch, t);
      da->query(value_vf, "vol_frac_CC", matl, patch, t);

      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c    = *iter;
        double vel_mag = value_velocity[c].length();
        double mass    = value_density[c] * value_vf[c];
        mean_vel += value_velocity[c] * mass;
        KE += mass * vel_mag * vel_mag;
        total_mass += mass;
      }
      /*
      ParticleVariable<Point> value_pos;
      ParticleVariable<Vector> value_vel;
      ParticleVariable<double> value_mass;
      da->query(value_pos, "p.x",        matl, patch, t);
      da->query(value_vel, "p.velocity", matl, patch, t);
      da->query(value_mass,"p.mass",     matl, patch, t);

      ParticleSubset* pset = value_pos.getParticleSubset();
      if(pset->numParticles() > 0){
        ParticleSubset::iterator piter = pset->begin();
        for(;piter != pset->end(); piter++){
          double vel_mag = value_vel[*piter].length();
          mean_vel+=value_vel[*piter]*value_mass[*piter];
          KE+=value_mass[*piter]*vel_mag*vel_mag;
          total_mass+=value_mass[*piter];
        } // for
        */
    } // for patches
    mean_vel /= total_mass;
    double mean_vel_mag = mean_vel.length();
    KE *= .5;

    outfile.precision(15);
    outfile << time << " " << mean_vel_mag << " " << total_mass << " " << KE
            << std::endl;
  }
} // end monica2()
