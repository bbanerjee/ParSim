/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __CCA_COMPONENTS_SCHEDULERS_RELOCATE_MPI_SCATTER_RECORD_H__
#define __CCA_COMPONENTS_SCHEDULERS_RELOCATE_MPI_SCATTER_RECORD_H__

#include <CCA/Components/Schedulers/RelocateScatterRecord.h>

#include <map>
#include <vector>

namespace Uintah {

using map_type     = std::multimap<std::pair<const Patch*, int>, ScatterRecord*>;
using patches_type = std::vector<const Patch*>;

struct MPIScatterProcessorRecord
{
  patches_type patches;
  void
  sortPatches();
};

using procmap_type = std::map<int, MPIScatterProcessorRecord*>;

struct MPIRecvBuffer
{
  MPIRecvBuffer* next;
  char* databuf;
  int bufsize;
  int numParticles;
  MPIRecvBuffer(char* databuf, int bufsize, int numParticles)
    : next(0)
    , databuf(databuf)
    , bufsize(bufsize)
    , numParticles(numParticles)
  {
  }
};

using recvmap_type = std::map<std::pair<const Patch*, int>, MPIRecvBuffer*>;

class MPIScatterRecords
{
public:

  ~MPIScatterRecords();

  ScatterRecord*
  findOrInsertRecord(const Patch* from,
                     const Patch* to,
                     int matl,
                     int curLevelIndex,
                     ParticleSubset* pset);
  ScatterRecord*
  findRecord(const Patch* from, const Patch* to, int matl, int curLevelIndex);

  void
  addNeighbor(LoadBalancer* lb, const ProcessorGroup* pg, const Patch* to);

  void
  saveRecv(const Patch* to,
           int matl,
           char* databuf,
           int bufsize,
           int numParticles);

  MPIRecvBuffer*
  findRecv(const Patch* to, int matl);

public:

  // map the to patch and matl to the ScatterRecord
  map_type records;
  procmap_type procs;
  recvmap_type recvs;
};

} // namespace Uintah

#endif //__CCA_COMPONENTS_SCHEDULERS_RELOCATE_MPI_SCATTER_RECORD_H__