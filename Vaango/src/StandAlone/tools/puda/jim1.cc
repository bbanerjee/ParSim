/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <StandAlone/tools/puda/jim1.h>

#include <StandAlone/tools/puda/util.h>

#include <Core/DataArchive/DataArchive.h>

#include <fstream>
#include <iomanip>
#include <vector>

using namespace Uintah;

using namespace std;

////////////////////////////////////////////////////////////////////////
//              J I M 1   O P T I O N
// This currently pulls out particle position, velocity and ID
// and prints that on one line for each particle.  This is useful
// for postprocessing of particle data for particles which move
// across patches.

void
Uintah::jim1(DataArchive* da, CommandLineFlags& clf)
{
  std::vector<std::string> vars;
  std::vector<int> num_matl;
  std::vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, num_matl, types);
  ASSERTEQ(vars.size(), types.size());
  std::cout << "#There are " << vars.size() << " variables:\n";
  for (int i = 0; i < (int)vars.size(); i++) {
    std::cout << vars[i] << ": " << types[i]->getName() << std::endl;
  }

  std::vector<int> index;
  std::vector<double> times;
  da->queryTimesteps(index, times);
  ASSERTEQ(index.size(), times.size());
  std::cout << "#There are " << index.size() << " timesteps:\n";
  //  for( int i = 0; i < (int)index.size(); i++ ) {
  //    std::cout << index[i] << ": " << times[i] << std::endl;
  //  }

  findTimestep_loopLimits(clf.tslow_set,
                          clf.tsup_set,
                          times,
                          clf.time_step_lower,
                          clf.time_step_upper);

  for (unsigned long t = clf.time_step_lower; t <= clf.time_step_upper;
       t += clf.time_step_inc) {
    double time = times[t];
    // cout << "time = " << time << std::endl;
    GridP grid = da->queryGrid(t);
    std::ostringstream fnum;
    string filename;
    fnum << setw(4) << setfill('0') << t / clf.time_step_inc;
    string partroot("partout");
    filename = partroot + fnum.str();
    //    ofstream partfile(filename.c_str());

    double tail_pos = -9.e99;
    double toe_pos  = -9.e99;

    for (int l = 0; l < grid->numLevels(); l++) {
      LevelP level = grid->getLevel(l);
      //      std::cout << "Level: " <<  endl;
      for (Level::const_patchIterator iter = level->patchesBegin();
           iter != level->patchesEnd();
           iter++) {
        const Patch* patch = *iter;
        int matl           = clf.matl_jim;
        //__________________________________
        //   P A R T I C L E   V A R I A B L E
        ParticleVariable<long64> value_pID;
        ParticleVariable<Point> value_pos;
        ParticleVariable<Vector> value_vel;
        ParticleVariable<double> value_vol, value_mas, value_tmp;
        ParticleVariable<Matrix3> value_strs;
        da->query(value_pos, "p.x", matl, patch, t);
        //        da->query(value_vel, "p.velocity",   matl, patch, t);
        //        da->query(value_vol, "p.volume",     matl, patch, t);
        //        da->query(value_mas, "p.mass",       matl, patch, t);
        //        da->query(value_tmp, "p.temperature",matl, patch, t);
        ParticleSubset* pset = value_pos.getParticleSubset();
        if (pset->numParticles() > 0) {
          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            tail_pos = max(value_pos[*iter].y(), tail_pos);
            toe_pos  = max(value_pos[*iter].x(), toe_pos);
          } // for
        }   // if
      }     // for patches
    }       // for levels
    std::cout << time << " " << tail_pos << " " << toe_pos << std::endl;
  }
} // end jim1()
