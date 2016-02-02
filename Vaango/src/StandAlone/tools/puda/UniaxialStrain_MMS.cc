/*
 * The MIT License
 *
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


#include <StandAlone/tools/puda/UniaxialStrain_MMS.h>
#include <StandAlone/tools/puda/util.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/InvalidGrid.h>

#include <iomanip>
#include <fstream>
#include <vector>

using namespace Uintah;
using namespace SCIRun;
using namespace std;

////////////////////////////////////////////////////////////////////////
//              UniaxialStrain_MMS   O P T I O N
// Compares the exact solution for the 1D uniaxial strain MMS problem with
// the MPM solution.  Computes the L-infinity error on the displacement
// at each timestep and reports to the command line.
//
// Reference:
//  ONR-MURI Report Dec. 30, 2015

void
Uintah::UniaxialStrain_MMS( DataArchive * da, CommandLineFlags & clf )
{
  vector<string> vars;
  vector<const Uintah::TypeDescription*> types;
  da->queryVariables(vars, types);
  ASSERTEQ(vars.size(), types.size());
  
  cout << "There are " << vars.size() << " variables:\n";
  
  for(int i=0;i<(int)vars.size();i++)
    cout << vars[i] << ": " << types[i]->getName() << endl;
      
  vector<int> index;
  vector<double> times;
  
  da->queryTimesteps(index, times);
  
  ASSERTEQ(index.size(), times.size());
  cout << "There are " << index.size() << " timesteps:\n";
  for( int i = 0; i < (int)index.size(); i++ ) {
    cout << index[i] << ": " << times[i] << endl;
  }
      
  findTimestep_loopLimits( clf.tslow_set, clf.tsup_set, times, clf.time_step_lower, clf.time_step_upper);
      
  // Hardcoded density (1.7 gm/cc) and elastic properties (K = 60 MPa, G = 100 MPa)
  double rho0 = 1700.0;
  double kappa = 6.0e7;
  double mu = 1.0e8;
  double lambda = kappa - (mu*2.0)/3.0;
  double cp = std::sqrt((lambda + 2.0*mu)/rho0);

  // Hardcoded amplitude (0.01 m) and frequency (10000 rad/s)
  double alpha = 0.01;
  double omega = 10000.0;

  for(unsigned long t=clf.time_step_lower;t<=clf.time_step_upper;t+=clf.time_step_inc){
    double time = times[t];
    GridP grid = da->queryGrid(t);
    
    
    //__________________________________
    //  bulletproofing
    IntVector low, high;
    grid->getLevel(0)->findInteriorCellIndexRange(low, high);
    IntVector cellNum = high-low;
    int dir = -9;
    int num1D_dirs=0;
    for(int d = 0; d <3; d++){
      if(cellNum[d] == 1){
        num1D_dirs +=1;
      }
      if(cellNum[d] != 1){
        dir = d;
      }
    }
    
    // Is the grid 1D?
    if (clf.do_US_MMS && num1D_dirs != 2){
      ostringstream warn;
      warn << "\nERROR: You cannot use the 1D MMS solution on a domain thatis not 1D. "
           << " Number of cells in each direction " <<cellNum << " num1D_dirs " << num1D_dirs <<" \n";
      throw InvalidGrid(warn.str(), __FILE__, __LINE__);
    
    }

    int    TotalNumParticles  = 0;   
    double max_errorAllLevels = 0.0;
    double TotalSumError      = 0.0;
    Point  worstPosAllLevels  = Point(-9,-9,-9);
    IntVector worstCellAllLevels = IntVector(-9,-9,-9);
    
    int numLevels = grid->numLevels();
    vector<double>    LinfLevel(numLevels);
    vector<double>    L2normLevel(numLevels);
    vector<Point>     worstPosLevel(numLevels);
    vector<IntVector> worstCellLevel(numLevels);
    vector<int>       numParticles(numLevels);
    
    //__________________________________
    //  Level loop
    for(int l=0;l<numLevels;l++){
      LevelP level = grid->getLevel(l);
    
      double sumError  = 0.0;
      double max_error = 0;
      numParticles[l]  = 0;
      Point worstPos   = Point(-9,-9,-9);
      IntVector worstCell(-9,-9,-9);
      
      //__________________________________
      // Patch loop
      for(Level::const_patchIterator iter = level->patchesBegin();
          iter != level->patchesEnd(); iter++){

        const Patch* patch = *iter;
        
        int matl = clf.matl_jim;
        ParticleVariable<Point>  value_pos;
        ParticleVariable<Vector> value_vel;

        da->query(value_pos,  "p.x",           matl, patch, t);
        //da->query(value_vel, "p.displacement",matl, patch, t);
        da->query(value_vel, "p.velocity",     matl, patch, t);
          
        ParticleSubset* pset = value_pos.getParticleSubset(); 
        numParticles[l] += pset->numParticles();
        
        //__________________________________
        //  Compute the error.       
        if(pset->numParticles() > 0){  // are there particles on this patch
        
          for (auto iter = pset->begin(); iter != pset->end(); iter++) {

            Point curx = value_pos[*iter];
    
            // Displacement
            Vector u_exact(0,0,0);
            double U  = alpha*std::cos(omega*(time - curx.x()/cp));
            u_exact[dir] = U;

            // Velocity
            Vector v_exact(0,0,0);
            double V  = -omega*alpha*std::sin(omega*(time - curx.x()/cp));
            v_exact[dir] = V;
          
            double error = (v_exact - value_vel[*iter]).length();
            cout << " x = " << curx(dir) << " v (exact) = " << V << " v (mpm) = " << value_vel[*iter].x()
                 <<  " error = "  << error << endl;
            sumError += error*error;
                
            if (error>max_error){
              max_error = error;
              worstPos  = value_pos[*iter];
              worstCell = patch->getCellIndex(worstPos);
            }
          }  // particle Loop
            
        }  //if
      }  // for patches
      LinfLevel[l]      = max_error;
      worstPosLevel[l]  = worstPos;
      worstCellLevel[l] = worstCell;
      
      if(sumError != 0){
        L2normLevel[l]    = sqrt( sumError/(double)numParticles[l]);
      }else{
        L2normLevel[l]    = 0.0;
      }
      
      cout << "     Level: " << level->getIndex() << " L_inf Error: " << LinfLevel[l] << ", L2norm: " << L2normLevel[l] 
           << " numParticles: " << numParticles[l] << " , Worst particle: " << worstPos << ", " << worstCell << endl;
      
      TotalSumError     += sumError;
      TotalNumParticles += numParticles[l];
      
      if (max_error > max_errorAllLevels) {
        max_errorAllLevels = max_error;
        worstPosAllLevels  = worstPos;
        worstCellAllLevels = worstCell;
      }
    }   // for levels
    double L2norm = sqrt( TotalSumError /(double)TotalNumParticles );
    
    cout << "time: " << time << " , L_inf Error: " << max_errorAllLevels << " , L2norm Error: "<< L2norm << " , Worst particle: " << worstPosAllLevels << " " << worstCellAllLevels << endl;
    
    //__________________________________
    // write data to the files (L_norms & L_normsPerLevels)
    FILE *outFile;
    
    // output level information
    outFile = fopen("L_normsPerLevel","w");
    fprintf(outFile, "#Time,  Level,   L_inf,    L2norm,    NumParticles\n");
    for(int l=0;l<numLevels;l++){
      fprintf(outFile, "%16.16le, %i,  %16.16le,  %16.16le  %i\n", time, l, LinfLevel[l], L2normLevel[l],numParticles[l]);
    }
    fclose(outFile);
    
    // overall 
    outFile = fopen("L_norms","w");
    fprintf(outFile, "#Time,    L_inf,    L2norm\n");
    fprintf(outFile, "%16.16le, %16.16le, %16.16le\n", time, max_errorAllLevels, L2norm);
    fclose(outFile); 
  }
} // end UniaxialStrain_MMS()

