/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef PATCH_BVH_NODE_H
#define PATCH_BVH_NODE_H

#include <Core/Grid/PatchBVH/PatchBVHBase.h>
#include <vector>

namespace Uintah {

/**************************************
  CLASS
  PatchBVHNode

  A Bounding Volume Hiearchy for querying patches that are
  within a given range.  This class is a general node of the tree.

  GENERAL INFORMATION

  PatchBVHNode.h

  Justin Luitjens
  Department of Computer Science
  University of Utah

  Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
 ****************************************/

class PatchBVHNode : public PatchBVHBase
{
public:
  PatchBVHNode(std::vector<PatchKeyVal>::iterator begin,
               std::vector<PatchKeyVal>::iterator end);

  ~PatchBVHNode() noexcept(false);

  void
  query(const IntVector& low,
        const IntVector& high,
        std::vector<const Patch*>& patches,
        bool includeExtraCells);

private:
  PatchBVHBase *left_, *right_;
};

} // end namespace Uintah

#endif
