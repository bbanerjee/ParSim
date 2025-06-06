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

#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/InvalidGrid.h>
#include <StandAlone/tools/puda/AA_MMS.h>
#include <StandAlone/tools/puda/util.h>

#include <fstream>
#include <iomanip>
#include <vector>

using namespace Uintah;

using namespace std;

////////////////////////////////////////////////////////////////////////
//              AA_MMS   O P T I O N
// Compares the exact solution for the 1D or 3D axis aligned MMS problem with
// the MPM solution.  Computes the L-infinity error on the displacement
// at each timestep and reports to the command line.
//
// Reference:
//  M. Steffen, P.C. Wallstedt, J.E. Guilkey, R.M. Kirby, and M. Berzins
//  "Examination and Analysis of Implementation Choices within the Material
//  Point Method (MPM) CMES, vol. 31, no. 2, pp. 107-127, 2008

void
Uintah::AA_MMS(DataArchive* da, CommandLineFlags& clf)
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

  for (unsigned long t = clf.time_step_lower; t <= clf.time_step_upper;
       t += clf.time_step_inc) {
    double time = times[t];
    GridP grid  = da->queryGrid(t);

    //__________________________________
    //  bulletproofing
    IntVector low, high;
    grid->getLevel(0)->findInteriorCellIndexRange(low, high);
    IntVector cellNum = high - low;
    int dir           = -9;
    int num1D_dirs    = 0;
    for (int d = 0; d < 3; d++) {
      if (cellNum[d] == 1) {
        num1D_dirs += 1;
      }
      if (cellNum[d] != 1) {
        dir = d;
      }
    }

    // Is the grid 1D?
    if (clf.do_AA_MMS_1 && num1D_dirs != 2) {
      std::ostringstream warn;
      warn << "\nERROR: You cannot use the 1D MMS solution on a domain thatis "
              "not 1D. "
           << " Number of cells in each direction " << cellNum << " num1D_dirs "
           << num1D_dirs << " \n";
      throw InvalidGrid(warn.str(), __FILE__, __LINE__);
    }
    // Is the grid 3D?
    if (clf.do_AA_MMS_2 && num1D_dirs != 0) {
      std::ostringstream warn;
      warn << "\nERROR: You cannot use the 3D MMS solution on a domain that is "
              "not 3D. "
           << " Number of cells in each direction " << cellNum << " \n";
      throw InvalidGrid(warn.str(), __FILE__, __LINE__);
    }

    //__________________________________
    //  hard coded constants!!!!
    double mu   = 3846.;
    double bulk = 8333.;
    double E    = 9. * bulk * mu / (3. * bulk + mu);
    double rho0 = 1.0;
    double c    = sqrt(E / rho0);
    //    double A0     = 1e-2;                    // << This is normalized
    //    below
    int TotalNumParticles        = 0;
    double max_errorAllLevels    = 0.0;
    double TotalSumError         = 0.0;
    Point worstPosAllLevels      = Point(-9, -9, -9);
    IntVector worstCellAllLevels = IntVector(-9, -9, -9);

    int numLevels = grid->numLevels();
    std::vector<double> LinfLevel(numLevels);
    std::vector<double> L2normLevel(numLevels);
    std::vector<Point> worstPosLevel(numLevels);
    std::vector<IntVector> worstCellLevel(numLevels);
    std::vector<int> numParticles(numLevels);

    //__________________________________
    //  Level loop
    for (int l = 0; l < numLevels; l++) {
      LevelP level = grid->getLevel(l);

      double sumError  = 0.0;
      double max_error = 0;
      numParticles[l]  = 0;
      Point worstPos   = Point(-9, -9, -9);
      IntVector worstCell(-9, -9, -9);

      //      Vector dx = level->dCell();             // you need to normalize
      //      the variable A by the
      double A = 0.05;

      //__________________________________
      // Patch loop
      for (Level::const_patch_iterator iter = level->patchesBegin();
           iter != level->patchesEnd();
           iter++) {

        const Patch* patch = *iter;

        int matl = clf.matl_jim;
        ParticleVariable<Point> value_pos;
        ParticleVariable<Vector> value_disp;

        da->query(value_pos, "p.x", matl, patch, t);
        da->query(value_disp, "p.displacement", matl, patch, t);

        ParticleSubset* pset = value_pos.getParticleSubset();
        numParticles[l] += pset->numParticles();

        //__________________________________
        //  Compute the error.
        if (pset->numParticles() > 0) { // are there particles on this patch

          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {

            Point refx = value_pos[*iter] - value_disp[*iter];
            Vector u_exact(0, 0, 0);

            //__________________________________
            //  Equation 47 of reference
            if (clf.do_AA_MMS_1) {
              double U = A * sin(2 * M_PI * refx(dir)) * cos(M_PI * c * time);
              u_exact[dir] = U;
            }

            //__________________________________
            //  Equation 51 of reference
            if (clf.do_AA_MMS_2) {
              u_exact = A * Vector(sin(M_PI * refx.x()) * sin(c * M_PI * time),
                                   sin(M_PI * refx.y()) *
                                     sin((2. / 3.) * M_PI + c * M_PI * time),
                                   sin(M_PI * refx.z()) *
                                     sin((4. / 3.) * M_PI + c * M_PI * time));
            }
            double error = (u_exact - value_disp[*iter]).length();
            std::cout << refx(dir) << " " << error << std::endl;
            sumError += error * error;

            if (error > max_error) {
              max_error = error;
              worstPos  = value_pos[*iter];
              worstCell = patch->getCellIndex(worstPos);
            }
          } // particle Loop

        } // if
      }   // for patches
      LinfLevel[l]      = max_error;
      worstPosLevel[l]  = worstPos;
      worstCellLevel[l] = worstCell;

      if (sumError != 0) {
        L2normLevel[l] = sqrt(sumError / (double)numParticles[l]);
      } else {
        L2normLevel[l] = 0.0;
      }

      std::cout << "     Level: " << level->getIndex()
                << " L_inf Error: " << LinfLevel[l]
                << ", L2norm: " << L2normLevel[l]
                << " numParticles: " << numParticles[l]
                << " , Worst particle: " << worstPos << ", " << worstCell
                << std::endl;

      TotalSumError += sumError;
      TotalNumParticles += numParticles[l];

      if (max_error > max_errorAllLevels) {
        max_errorAllLevels = max_error;
        worstPosAllLevels  = worstPos;
        worstCellAllLevels = worstCell;
      }
    } // for levels
    double L2norm = sqrt(TotalSumError / (double)TotalNumParticles);

    std::cout << "time: " << time << " , L_inf Error: " << max_errorAllLevels
              << " , L2norm Error: " << L2norm
              << " , Worst particle: " << worstPosAllLevels << " "
              << worstCellAllLevels << std::endl;

    //__________________________________
    // write data to the files (L_norms & L_normsPerLevels)
    FILE* outFile;

    // output level information
    outFile = fopen("L_normsPerLevel", "w");
    fprintf(outFile, "#Time,  Level,   L_inf,    L2norm,    NumParticles\n");
    for (int l = 0; l < numLevels; l++) {
      fprintf(outFile,
              "%16.16le, %i,  %16.16le,  %16.16le  %i\n",
              time,
              l,
              LinfLevel[l],
              L2normLevel[l],
              numParticles[l]);
    }
    fclose(outFile);

    // overall
    outFile = fopen("L_norms", "w");
    fprintf(outFile, "#Time,    L_inf,    L2norm\n");
    fprintf(outFile,
            "%16.16le, %16.16le, %16.16le\n",
            time,
            max_errorAllLevels,
            L2norm);
    fclose(outFile);
  }
} // end AA_MMS()
