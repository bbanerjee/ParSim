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


#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>


#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Patch.h>

using namespace Uintah;

#include "mpi.h"
#include <testprograms/Regridders/TiledRegridder.h>
#include <testprograms/Regridders/BNRRegridder.h>
#include <testprograms/Regridders/LBNRRegridder.h>
#include <testprograms/Regridders/GBRv1Regridder.h>
#include <testprograms/Regridders/GBRv2Regridder.h>

#include <testprograms/Regridders/common.h>

#define REPEAT1 5
#define REPEAT2 400
void outputTime(double time, int alg);

int main(int argc, char **argv) 
{
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(argc!=5)
  {
    if(rank==0)
    {
      std::cout << "Usage: benchmark patch_size number_of_patches flag_inner_rad(0-1) flag_outter_rad(0-1)\n";
      std::cout << " Command was: ";
      for(int i=0;i<argc;i++)
        std::cout << argv[i] << " ";
      std::cout << std::endl;
    }
    MPI_Finalize();
    return 1;
  }

  //get command line arguements
  IntVector patch_size;
  patch_size[0]=patch_size[1]=patch_size[2]=atoi(argv[1]);
  IntVector num_patches(1,1,1);

  int target_patches=atoi(argv[2]);
  int d=0;
  while(num_patches[0]*num_patches[1]*num_patches[2]<target_patches)
  {
    num_patches[d]*=2;
    d=(d+1)%3;
  }
  //cout << "Num patches: " << num_patches << " total:" << num_patches[0]*num_patches[1]*num_patches[2] << std::endl;
  IntVector cells;
  cells[0]=cells[1]=cells[2]=num_patches[0]*patch_size[0];
  IntVector rr(4,4,4);
#if 0
    patch_size[2]=1;
    num_patches[2]=1;
    cells[2]=1;
    rr[2]=1;
#endif 

  double radin=cells[0]/2.0*atof(argv[3]);
  double radout=cells[0]/2.0*atof(argv[4]);

  Sphere2 s(cells.asVector()/Vector(2,2,2),radin,radout);
  
  //cout << "num patches:" << num_patches << std::endl;
  //cout << "cells: " << cells << std::endl;
  //cout << "rad: " << rad << std::endl;

  //create coarse patch set
  std::vector<Region> patches;
  std::vector<CCVariable<int> * > flags;
  std::vector<IntVector> gflags;
  std::vector<list<IntVector> > lflags;

  int total_patches=num_patches.x()*num_patches.y()*num_patches.z();
  int div=total_patches/num_procs;
  int mod=total_patches%num_procs;

  //cout << "total patches: " << total_patches << std::endl;
  //cout << "div: " << div << " mod: " << mod << std::endl;
  int p=0;
  int p_assigned=0;
  int to_assign=div+int(mod>0);
  int idx=0;
  //cout << "to_assign=" << to_assign << std::endl;

  for(int i=0;i<num_patches.x();i++)
  {
    for(int j=0;j<num_patches.y();j++)
    {
      for(int k=0;k<num_patches.z();k++)
      {
        IntVector low=IntVector(i,j,k)*patch_size;
        IntVector high=low+patch_size;
        if(p==rank)
        {
          patches.push_back(Region(low,high));
        }
        
        p_assigned++;

        if(p_assigned==to_assign)
        { 
          //if(rank==0)
          //  std::cout << p << " assigned: " << to_assign << " patches\n"; 
          p++;
          p_assigned=0;
          to_assign=div+int(mod>p);
        }
        idx++;
      }
    }
  }

  //for(unsigned int i=0;i<patches.size();i++)
  //  std::cout << rank << " patch: " << patches[i] << std::endl;

  //create refinement flags
  flags.resize(patches.size());
  lflags.resize(patches.size());

  //for each local patch
  for(unsigned int patch=0;patch<patches.size();patch++)
  {
    //allocate the variable
    CCVariable<int> *f=new CCVariable<int>();

    //allocate the memory in the variable
    f->allocate(patches[patch].getLow(),patches[patch].getHigh());
    flags[patch]=f;
    
    //determine flag set
    for(CellIterator iter(patches[patch].getLow(),patches[patch].getHigh());!iter.done();iter++)
    {
      if(s.intersects( (*iter).asVector()+Vector(.5,.5,.5)))
      {
        //add to CC variable flags
        (*f)[*iter]=1;
        //add to per patch flags list
        lflags[patch].push_back(*iter);
        //add to per processor flags list
        gflags.push_back(*iter);
      }
      else
        (*f)[*iter]=0;
    }
  }

  ofstream fout;
  
  std::vector<Region> fine_patches,global_patches;

  TiledRegridder tiled(patch_size,rr);
  LBNRRegridder lbnr(.85,rr);
  GBRv1Regridder gbrv1(.85,rr,rank,num_procs);
  GBRv2Regridder gbrv2(.85,rr,rank,num_procs);

  std::cout << setprecision(20);
  clock_t start;
  double time;
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT2;i++)
    tiled.regrid(patches,flags,fine_patches);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT2;
  outputTime(time,0);

#endif
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT2;i++)
    lbnr.regrid(lflags,fine_patches);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT2;
  outputTime(time,1);
#endif
#if 1
  std::vector<IntVector> tmpflags;
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT2;i++)
    makeFlagsList(patches, flags, tmpflags);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT2;
  outputTime(time,2);
#endif
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT2;i++)
    gatherPatches(fine_patches, global_patches);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT2;
  outputTime(time,3);
#endif
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT2;i++)
    splitPatches(fine_patches, global_patches,.25);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT2;
  outputTime(time,4);
#endif
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT1;i++)
    gbrv2.regrid(gflags,fine_patches);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT1;
  outputTime(time,5);
#endif
#if 1
  MPI_Barrier(MPI_COMM_WORLD);
  start=clock();
  for(int i=0;i<REPEAT1;i++)
    gbrv1.regrid(gflags,fine_patches);
  time=(clock()-start)/(double) CLOCKS_PER_SEC / REPEAT1;
  outputTime(time,6);
#endif

  MPI_Finalize();
}

void getTime(double time, double &mint, double &maxt, double &avgt)
{
  Uintah::MPI::Allreduce(&time,&mint,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  Uintah::MPI::Allreduce(&time,&maxt,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  Uintah::MPI::Allreduce(&time,&avgt,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  avgt/=num_procs;
}

void outputTime(double time, int alg)
{
  double mint,maxt,avgt;
  getTime(time,mint,maxt,avgt);

  if(rank==0)
    std::cout << num_procs << " " << alg << " " << avgt << " " << mint << " " << maxt << std::endl; 

}


