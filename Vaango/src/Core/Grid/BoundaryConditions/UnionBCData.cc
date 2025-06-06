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
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/BoundaryConditions/UnionBCData.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/UnionIterator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>
#include <algorithm>
#include <iostream>

namespace {

// Usage: export SCI_DEBUG="UnionBC_dbg:+"
Uintah::Dout bc_dbg{ "UnionBC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Union BC debug info",
                     false };

} // namespace

namespace Uintah {

UnionBCData::UnionBCData()
  : BCGeomBase()
{
}

UnionBCData::~UnionBCData()
{
  child.clear();
}

UnionBCData::UnionBCData(const UnionBCData& mybc)
  : BCGeomBase(mybc)
{
  for (auto& bc : mybc.child) {
    child.push_back(bc->clone());
  }
}

auto
UnionBCData::operator=(const UnionBCData& rhs) -> UnionBCData&
{
  BCGeomBase::operator=(rhs);

  if (this == &rhs) {
    return *this;
  }

  // Delete the lhs
  child.clear();

  // copy the rhs to the lhs
  for (auto& bc : rhs.child) {
    child.push_back(bc->clone());
  }

  return *this;
}

auto
UnionBCData::operator==(const BCGeomBase& rhs) const -> bool
{
  const auto* p_rhs = dynamic_cast<const UnionBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    if (this->child.size() != p_rhs->child.size()) {
      return false;
    }

    return equal(this->child.begin(), this->child.end(), p_rhs->child.begin());
  }
}

auto
UnionBCData::clone() -> std::shared_ptr<BCGeomBase>
{
  return std::make_shared<UnionBCData>(*this);
}

void
UnionBCData::addBCData([[maybe_unused]] BCData& bc)
{
}

void
UnionBCData::addBC([[maybe_unused]] BoundCondBaseSP bc)
{
}

void
UnionBCData::sudoAddBC(BoundCondBaseSP& bc)
{
  for (auto& i : child) {
    i->sudoAddBC(bc); // or add to zero element only?
  }
}

void
UnionBCData::addBCData(std::shared_ptr<BCGeomBase> bc)
{
  child.push_back(bc);
}

void
UnionBCData::getBCData(BCData& bc) const
{
  child[0]->getBCData(bc);
}

auto
UnionBCData::inside(const Point& p) const -> bool
{
  for (auto i : child) {
    if (i->inside(p)) {
      return true;
    }
  }
  return false;
}

void
UnionBCData::print()
{
  DOUT(bc_dbg, "Geometry type = " << typeid(this).name());
  for (auto i : child) {
    i->print();
  }
}

void
UnionBCData::determineIteratorLimits(Patch::FaceType face,
                                     const Patch* patch,
                                     std::vector<Point>& test_pts)
{
#if 0
  std::cout << "UnionBC determineIteratorLimits()" << std::endl;
#endif

  for (auto bc : child) {
    bc->determineIteratorLimits(face, patch, test_pts);
  }

  UnionIterator cells;
  UnionIterator nodes;
  for (auto bc : child) {
    Iterator cell_itr, node_itr;
    bc->getCellFaceIterator(cell_itr);
    bc->getNodeFaceIterator(node_itr);
    Iterator base_ci, base_ni;
    cells = UnionIterator(base_ci, cell_itr);
    nodes = UnionIterator(base_ni, node_itr);
  }
  d_cells = UnionIterator(cells);
  d_nodes = UnionIterator(nodes);

#if 0
  IntVector l,h;
  patch->getFaceCells(face,0,l,h);

  std::vector<IntVector> b,nb;
  std::vector<Point>::iterator pts;
  pts = test_pts.begin();
  for (CellIterator bound(l,h); !bound.done(); bound++,pts++) 
    if (inside(*pts))
      b.push_back(*bound);

  setBoundaryIterator(b);
#if 0
  std::cout << "Size of boundary = " << boundary.size() << std::endl;
#endif
  // Need to determine the boundary iterators for each separate bc.
  for (std::vector<BCGeomBase*>::const_iterator bc = child.begin();  
       bc != child.end(); ++bc) {
    pts = test_pts.begin();
    std::vector<IntVector> boundary_itr;
    for (CellIterator bound(l,h); !bound.done(); bound++, pts++) 
      if ( (*bc)->inside(*pts))
        boundary_itr.push_back(*bound);
#if 0
    std::cout << "Size of boundary_itr = " << boundary_itr.size() << std::endl;
#endif
    (*bc)->setBoundaryIterator(boundary_itr);
  }
    
  IntVector ln,hn;
  patch->getFaceNodes(face,0,ln,hn);
  for (NodeIterator bound(ln,hn);!bound.done();bound++) {
    Point p = patch->getLevel()->getNodePosition(*bound);
    if (inside(p)) 
      nb.push_back(*bound);
  }
  
  setNBoundaryIterator(nb);

#endif
}

} // namespace Uintah