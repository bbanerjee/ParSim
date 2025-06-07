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

#include <Core/Grid/PatchBVH/PatchBVH.h>
#include <Core/Grid/PatchBVH/PatchBVHLeaf.h>
#include <Core/Grid/PatchBVH/PatchBVHNode.h>

namespace Uintah {

PatchBVH::PatchBVH(const std::vector<const Patch*>& patches)
  : root_(nullptr)
{
  if (patches.size() == 0) {
    return;
  }

  for (auto const& patch_p : patches) {
    PatchKeyVal key;

    key.patch = patch_p;
    key.center2 =
      patch_p->getExtraCellLowIndex() + patch_p->getExtraCellHighIndex();

    patches_.push_back(key);
  }

  if (patches_.size() > PatchBVHBase::getLeafSize()) {
    // create a node
    root_ = new PatchBVHNode(patches_.begin(), patches_.end());
  } else {
    // create a leaf
    root_ = new PatchBVHLeaf(patches_.begin(), patches_.end());
  }
}

PatchBVH::PatchBVH(const std::vector<Patch*>& patches)
  : root_(nullptr)
{
  if (patches.size() == 0) {
    return;
  }

  for (auto const& patch_p : patches) {
    PatchKeyVal key;

    key.patch = patch_p;
    key.center2 =
      patch_p->getExtraCellLowIndex() + patch_p->getExtraCellHighIndex();

    patches_.push_back(key);
  }

  if (patches_.size() > PatchBVHBase::getLeafSize()) {
    // create a node
    root_ = new PatchBVHNode(patches_.begin(), patches_.end());
  } else {
    // create a leaf
    root_ = new PatchBVHLeaf(patches_.begin(), patches_.end());
  }
}

PatchBVH::~PatchBVH()
{
  if (root_ != nullptr) {
    delete root_;
  }

  patches_.clear();
}

void
PatchBVH::query(const IntVector& low,
                const IntVector& high,
                std::vector<const Patch*>& patches,
                bool includeExtraCells)
{
  // verify query range is valid
  if (high.x() <= low.x() || high.y() <= low.y() || high.z() <= low.z()) {
    return;
  }

  if (root_ == nullptr) {
    return;
  }

  root_->query(low, high, patches, includeExtraCells);
}

} // end namespace Uintah
