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

#ifndef __CCA_COMPONENTS_SCHEDULERS_RELOCATE_SCATTER_RECORD_H__
#define __CCA_COMPONENTS_SCHEDULERS_RELOCATE_SCATTER_RECORD_H__

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ParticleSubset.h>

namespace Uintah {

struct ScatterRecord
{
  const Patch* from_patch;
  const Patch* to_patch;
  IntVector vector_to_neighbor;
  int matl;
  int level_index;
  ParticleSubset* send_pset;

  ScatterRecord(const Patch* fromPatch,
                const Patch* toPatch,
                int matl,
                int levelIndex)
    : from_patch(fromPatch)
    , to_patch(toPatch)
    , matl(matl)
    , level_index(levelIndex)
    , send_pset(0)
  {
    ASSERT(fromPatch != 0);
    ASSERT(toPatch != 0);

    vector_to_neighbor =
      toPatch->getExtraCellLowIndex() - fromPatch->getExtraCellLowIndex();
  }

  // Note that when the ScatterRecord going from a real patch to
  // a virtual patch has an equivalent representation going from
  // a virtual patch to a real patch (wrap-around, periodic bound. cond.).
  bool
  equivalent(const ScatterRecord& sr)
  {
    return (to_patch->getRealPatch() == sr.to_patch->getRealPatch()) &&
           (matl == sr.matl) && (vector_to_neighbor == sr.vector_to_neighbor);
  }
};

std::ostream&
operator<<(std::ostream& out, const ScatterRecord& r)
{
  out.setf(std::ios::scientific, std::ios::floatfield);
  out.precision(4);
  out << " Scatter Record, matl: " << r.matl << " Level: " << r.level_index
      << " numParticles " << r.send_pset->numParticles()
      << " (Particle moving from Patch " << r.from_patch->getID()
      << ", to Patch " << r.to_patch->getID() << ")"
      << " vectorToNeighbor " << r.vector_to_neighbor;
  out.setf(std::ios::scientific, std::ios::floatfield);
  return out;
}

} // namespace Uintah

#endif //__CCA_COMPONENTS_SCHEDULERS_RELOCATE_SCATTER_RECORD_H__