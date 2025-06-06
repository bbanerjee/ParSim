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

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Region.h>
#include <testprograms/Regridders/BNRRegridder.h>
#include "mpi.h"
#include <vector>

using namespace Uintah;
using namespace Uintah;

void BOUNDS( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
   Region* inbox=(Region*)(invec);
   Region* outbox=(Region*)(inoutvec);

   outbox->extend(*inbox);
}

class GBRv1Regridder : public BNRRegridder
{
  public:
    GBRv1Regridder(double tol,IntVector rr,int rank, int processors) : BNRRegridder(tol,rr), rank(rank), num_procs(processors)
    {
      MPI_Op_create(BOUNDS, true, &BOUNDS_OP);
    };
    virtual ~GBRv1Regridder() 
    {
      //MPI_Op_free(&BOUNDS_OP);
    };

  private:
    Region computeBounds(std::list<IntVector> &flags);
    void computeHistogram(std::list<IntVector> &flags, const Region &bounds, std::vector<std::vector<unsigned int> > &hist);
    int computeNumFlags(std::list<IntVector> &flags);
    int rank,num_procs;
    MPI_Op BOUNDS_OP;
  
};

Region GBRv1Regridder::computeBounds(std::list<IntVector> &flags)
{
  Region bounds=BNRRegridder::computeBounds(flags);
  
  if(flags.size()==0)
    bounds=Region(IntVector(INT_MAX,INT_MAX,INT_MAX),IntVector(INT_MIN,INT_MIN,INT_MIN));

  //cout << getpid() << " local flags: " << std::endl;
  //for(std::list<IntVector>::iterator iter=flags.begin();iter!=flags.end();iter++)
  //  std::cout << "      " << getpid() << " " << *iter << std::endl;
  //cout << getpid() << " local bounds: " << bounds << std::endl;
 
  Region gbounds;
  //all reduce bounds
  Uintah::MPI::Allreduce(&bounds,&gbounds,sizeof(Region),MPI_BYTE,BOUNDS_OP,MPI_COMM_WORLD);
  //cout << getpid() << " global bounds: " << bounds << std::endl;
  return gbounds;
}
    
int GBRv1Regridder::computeNumFlags(std::list<IntVector> &flags)
{
  unsigned int nflags=flags.size();
  unsigned int ngflags;
  
  //all reduce bounds
  Uintah::MPI::Allreduce(&nflags,&ngflags,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
  //cout << " num local flags: " << nflags << " global: " << ngflags << std::endl;
  return ngflags;
}
    
void GBRv1Regridder::computeHistogram(std::list<IntVector> &flags, const Region &bounds, std::vector<std::vector<unsigned int> > &hist)
{
  std::vector<std::vector<unsigned int> > lhist;
  BNRRegridder::computeHistogram(flags,bounds,lhist);
 
  hist.resize(lhist.size());
  //all reduce hist
  for(size_t d=0;d<hist.size();d++)
  {
    hist[d].resize(lhist[d].size());
    Uintah::MPI::Allreduce(&lhist[d][0],&hist[d][0],lhist[d].size(),MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
  }
}

