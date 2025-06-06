/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <Core/Geometry/Point.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/BoundaryConditions/DifferenceBCData.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Variables/DifferenceIterator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>
#include <algorithm>
#include <iostream>
#include <set>

namespace {

// Usage: export SCI_DEBUG="DiffBC_dbg:+"
Uintah::Dout bc_dbg{ "DiffBC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Diff BC debug info",
                     false };

} // namespace
namespace Uintah {

DifferenceBCData::DifferenceBCData(std::shared_ptr<BCGeomBase> p1,
                                   std::shared_ptr<BCGeomBase> p2)
  : BCGeomBase()
  , left(p1->clone())
  , right(p2->clone())
{
}

DifferenceBCData::DifferenceBCData(const DifferenceBCData& rhs)
  : BCGeomBase(rhs)
{
  left  = rhs.left->clone();
  right = rhs.right->clone();
}

auto
DifferenceBCData::operator=(const DifferenceBCData& rhs) -> DifferenceBCData&
{
  BCGeomBase::operator=(rhs);

  if (this == &rhs) {
    return *this;
  }

  // Copy the rhs to the lhs
  left  = rhs.left->clone();
  right = rhs.right->clone();

  return *this;
}

DifferenceBCData::~DifferenceBCData() {}

auto
DifferenceBCData::operator==(const BCGeomBase& rhs) const -> bool
{
  const auto* p_rhs = dynamic_cast<const DifferenceBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    return (this->left == p_rhs->left) && (this->right == p_rhs->right);
  }
}

std::shared_ptr<BCGeomBase>
DifferenceBCData::clone()
{
  return std::make_shared<DifferenceBCData>(*this);
}

void
DifferenceBCData::addBCData([[maybe_unused]] BCData& bc)
{
}

void
DifferenceBCData::addBC([[maybe_unused]] BoundCondBaseSP bc)
{
}

void
DifferenceBCData::sudoAddBC(BoundCondBaseSP& bc)
{
  left->sudoAddBC(bc);
}

void
DifferenceBCData::getBCData(BCData& bc) const
{
  left->getBCData(bc);
}

auto
DifferenceBCData::inside(const Point& p) const -> bool
{
  return (left->inside(p) && !right->inside(p));
}

void
DifferenceBCData::print()
{
#if 1
  DOUT(bc_dbg, "Difference Geometry type = " << typeid(this).name());
  DOUT(bc_dbg, "Left");
#endif
  left->print();
#if 1
  DOUT(bc_dbg, "Right");
#endif
  right->print();
}

void
DifferenceBCData::determineIteratorLimits(Patch::FaceType face,
                                          const Patch* patch,
                                          std::vector<Point>& test_pts)
{

#if 0
  std::cout << "DifferenceBC determineIteratorLimits() " << patch->getFaceName(face)<< std::endl;
#endif

  left->determineIteratorLimits(face, patch, test_pts);
  right->determineIteratorLimits(face, patch, test_pts);

  Iterator left_cell, left_node, right_cell, right_node;
  left->getCellFaceIterator(left_cell);
  left->getNodeFaceIterator(left_node);
  right->getCellFaceIterator(right_cell);
  right->getNodeFaceIterator(right_node);

  d_cells = DifferenceIterator(left_cell, right_cell);
  d_nodes = DifferenceIterator(left_node, right_node);

#if 0
#if 0
  std::cout << "DifferenceBC determineIteratorLimits()" << std::endl;
  std::cout << "Doing left determineIteratorLimits()" << std::endl;
#endif
  left->determineIteratorLimits(face,patch,test_pts);
#if 0
  std::cout << "Doing right determineIteratorLimits()" << std::endl;
#endif
  right->determineIteratorLimits(face,patch,test_pts);

#if 0
  std::cout << "Size of boundary = " << boundary.size() << std::endl;
  std::cout << "Size of nboundary = " << nboundary.size() << std::endl;
#endif

  // Need to do the set difference operations for the left and right to get
  // the boundary and nboundary iterators.
  std::vector<IntVector> diff_boundary,   diff_nboundary;
  std::vector<IntVector> *left_boundary,  *right_boundary;
  std::vector<IntVector> *left_nboundary, *right_nboundary;

  left->getBoundaryIterator(left_boundary);
  left->getNBoundaryIterator(left_nboundary);

  right->getBoundaryIterator(right_boundary);
  right->getNBoundaryIterator(right_nboundary);

#if 0
  std::cout << "Size of left_boundary = " << left_boundary->size() << std::endl;
  std::cout << "Size of left_nboundary = " << left_nboundary->size() << std::endl;
  std::cout << "Size of right_boundary = " << right_boundary->size() << std::endl;
  std::cout << "Size of right_nboundary = " << right_nboundary->size() << std::endl;
#endif
  
  for (std::vector<IntVector>::const_iterator it = left_boundary->begin();
       it != left_boundary->end(); ++it) {
    std::vector<IntVector>::const_iterator result = find(right_boundary->begin(),
                                                    right_boundary->end(),*it);
    if (result == right_boundary->end())
      diff_boundary.push_back(*it);
  }

  for (std::vector<IntVector>::const_iterator it = left_nboundary->begin();
       it != left_nboundary->end(); ++it) {
    std::vector<IntVector>::const_iterator result = find(right_nboundary->begin(),
                                                    right_nboundary->end(),*it);
    if (result == right_nboundary->end())
      diff_nboundary.push_back(*it);
  }

  setBoundaryIterator(diff_boundary);
  setNBoundaryIterator(diff_nboundary);

#if 0
  std::cout << "Size of boundary = " << boundary->size() << std::endl;
  std::cout << "Size of nboundary = " << nboundary->size() << std::endl;
#endif

#endif
}

} // namespace Uintah
