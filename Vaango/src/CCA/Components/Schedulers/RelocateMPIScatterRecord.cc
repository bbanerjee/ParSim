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

#include <CCA/Components/Schedulers/RelocateMPIScatterRecord.h>

#include <CCA/Components/Schedulers/RelocateScatterRecord.h>

#include <Core/Grid/Patch.h>
#include <Core/Util/FancyAssert.h>

namespace Uintah {

static bool
ComparePatches(const Patch* p1, const Patch* p2)
{
  return p1->getID() < p2->getID();
}

void
MPIScatterProcessorRecord::sortPatches()
{
  sort(patches.begin(), patches.end(), ComparePatches);
}

void
MPIScatterRecords::saveRecv(const Patch* to,
                            int matl,
                            char* databuf,
                            int datasize,
                            int numParticles)
{
  recvmap_type::key_type key(to, matl);
  recvmap_type::iterator iter = recvs.find(key);
  MPIRecvBuffer* record = scinew MPIRecvBuffer(databuf, datasize, numParticles);

  if (iter == recvs.end()) {
    recvs[key] = record;
  } else {
    record->next = iter->second;
    recvs[key]   = record;
  }
}

MPIRecvBuffer*
MPIScatterRecords::findRecv(const Patch* to, int matl)
{
  recvmap_type::iterator iter = recvs.find(std::make_pair(to, matl));
  if (iter == recvs.end()) {
    return 0;
  } else {
    return iter->second;
  }
}

ScatterRecord*
MPIScatterRecords::findOrInsertRecord(const Patch* fromPatch,
                                      const Patch* toPatch,
                                      int matl,
                                      int curLevelIndex,
                                      ParticleSubset* pset)
{
  ASSERT(toPatch != 0);
  IntVector vectorToNeighbor =
    toPatch->getExtraCellLowIndex() - fromPatch->getExtraCellLowIndex();
  const Patch* realToPatch = toPatch->getRealPatch();

  std::pair<map_type::iterator, map_type::iterator> pr =
    records.equal_range(std::make_pair(realToPatch, matl));

  // loop over all scatter records
  // Does this record exist if so return.
  for (; pr.first != pr.second; pr.first++) {
    ScatterRecord* sr = pr.first->second;

    if ((realToPatch == sr->to_patch->getRealPatch()) &&
        (curLevelIndex == sr->level_index) &&
        (vectorToNeighbor == sr->vector_to_neighbor)) {
      return sr;
    }
  }

  //  Create a new record and insert it into
  //  all records
  ScatterRecord* rec =
    scinew ScatterRecord(fromPatch, toPatch, matl, curLevelIndex);
  rec->send_pset = scinew ParticleSubset(0, -1, 0);
  records.insert(map_type::value_type(std::make_pair(realToPatch, matl), rec));
  return rec;
}

ScatterRecord*
MPIScatterRecords::findRecord(const Patch* fromPatch,
                              const Patch* toPatch,
                              int matl,
                              int curLevelIndex)
{
  ASSERT(toPatch != 0);

  IntVector vectorToNeighbor =
    toPatch->getExtraCellLowIndex() - fromPatch->getExtraCellLowIndex();
  const Patch* realToPatch   = toPatch->getRealPatch();
  const Patch* realFromPatch = fromPatch->getRealPatch();

  std::pair<map_type::iterator, map_type::iterator> pr =
    records.equal_range(std::make_pair(realToPatch, matl));

  // loop over all scatter records
  // Does this record exist if so return.
  for (; pr.first != pr.second; pr.first++) {
    ScatterRecord* sr = pr.first->second;

    if ((realToPatch == sr->to_patch->getRealPatch()) &&
        (curLevelIndex == sr->level_index) && (matl == sr->matl) &&
        (vectorToNeighbor == sr->vector_to_neighbor) &&
        (realFromPatch == sr->from_patch->getRealPatch()) &&
        (fromPatch != toPatch)) {
      return sr;
    }
  }

  return 0;
}

void
MPIScatterRecords::addNeighbor(LoadBalancer* load_balancer,
                               const ProcessorGroup* pg,
                               const Patch* neighbor)
{
  neighbor   = neighbor->getRealPatch();
  int toProc = load_balancer->getPatchwiseProcessorAssignment(neighbor);
  ASSERTRANGE(toProc, 0, pg->nRanks());

  procmap_type::iterator iter = procs.find(toProc);

  if (iter == procs.end()) {
    MPIScatterProcessorRecord* pr = scinew MPIScatterProcessorRecord();
    procs[toProc]                 = pr;
    pr->patches.push_back(neighbor);
  } else {
    // This is linear, with the hope that the number of patches per
    // processor will not be huge.
    MPIScatterProcessorRecord* pr = iter->second;
    int i;
    for (i = 0; i < (int)pr->patches.size(); i++) {
      if (pr->patches[i] == neighbor) {
        break;
      }
    }
    if (i == (int)pr->patches.size()) {
      pr->patches.push_back(neighbor);
    }
  }
}

MPIScatterRecords::~MPIScatterRecords()
{
  for (procmap_type::iterator iter = procs.begin(); iter != procs.end();
       iter++) {
    delete iter->second;
  }

  for (map_type::iterator mapiter = records.begin(); mapiter != records.end();
       mapiter++) {
    delete mapiter->second->send_pset;
    delete mapiter->second;
  }

  for (recvmap_type::iterator iter = recvs.begin(); iter != recvs.end();
       iter++) {
    MPIRecvBuffer* p = iter->second;

    while (p) {
      MPIRecvBuffer* next = p->next;
      delete p;
      p = next;
    }
  }
}

} // namespace Uintah
