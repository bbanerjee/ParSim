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
#include <Core/Grid/BoundaryConditions/BCGeomBase.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/GridIterator.h>
#include <Core/Grid/Variables/ListOfCellsIterator.h>
#include <Core/Util/DOUT.hpp>
#include <iostream>
#include <vector>


namespace {

// Usage: export SCI_DEBUG="BCGeom_dbg:+"
Uintah::Dout bc_dbg{ "BCGeom_dbg",
                     "Grid_BoundaryConditions",
                     "Grid BC Geom debug info",
                     false };

}

namespace Uintah {
BCGeomBase::BCGeomBase()
{
  d_cells           = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  d_nodes           = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  d_bcname          = "NotSet";
  d_bndtype         = "None";
  d_particleBndSpec = ParticleBndSpec(ParticleBndSpec::NOTSET,
                                      ParticleBndSpec::ELASTIC,
                                      0.0,
                                      0.0);
  d_surfaceArea     = 0.0;
  d_origin          = Uintah::Point(0, 0, 0);
}

BCGeomBase::BCGeomBase(const BCGeomBase& rhs)
{
  d_cells           = rhs.d_cells;
  d_nodes           = rhs.d_nodes;
  d_bcname          = rhs.d_bcname;
  d_bndtype         = rhs.d_bndtype;
  d_particleBndSpec = rhs.d_particleBndSpec;
  d_surfaceArea     = rhs.d_surfaceArea;
  d_origin          = rhs.d_origin;
}

auto
BCGeomBase::operator=(const BCGeomBase& rhs) -> BCGeomBase&
{
  if (this == &rhs) {
    return *this;
  }

  d_cells           = rhs.d_cells;
  d_nodes           = rhs.d_nodes;
  d_bcname          = rhs.d_bcname;
  d_bndtype         = rhs.d_bndtype;
  d_particleBndSpec = rhs.d_particleBndSpec;
  d_surfaceArea     = rhs.d_surfaceArea;
  d_origin          = rhs.d_origin;

  return *this;
}

BCGeomBase::~BCGeomBase() = default;

void
BCGeomBase::getCellFaceIterator(Iterator& b_ptr)
{
  b_ptr = d_cells;
}

void
BCGeomBase::getNodeFaceIterator(Iterator& b_ptr)
{
  b_ptr = d_nodes;
}

void
BCGeomBase::determineIteratorLimits(Patch::FaceType face,
                                    const Patch* patch,
                                    std::vector<Point>& test_pts)
{
#if 0
  std::cout << "BCGeomBase determineIteratorLimits() " << patch->getFaceName(face)<< std::endl;
#endif

  IntVector l, h;
  patch->getFaceCells(face, 0, l, h);
  GridIterator cells(l, h);

  IntVector ln, hn;
  patch->getFaceNodes(face, 0, ln, hn);
  GridIterator nodes(ln, hn);

  Iterator cell_itr(cells), node_itr(nodes);

  auto pts = test_pts.begin();

  std::vector<IntVector> vec_cells;
  for (cell_itr.reset(); !cell_itr.done(); cell_itr++, pts++) {
    if (inside(*pts)) {
      vec_cells.push_back(*cell_itr);
    }
  }

  std::vector<IntVector> vec_nodes;
  for (node_itr.reset(); !node_itr.done(); node_itr++) {
    Point p = patch->getLevel()->getNodePosition(*node_itr);
    if (inside(p)) {
      vec_nodes.push_back(*node_itr);
    }
  }

  if (vec_cells.empty()) {
    d_cells = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  } else {
    ListOfCellsIterator list_cells(vec_cells.size());
    for (const auto& cell : vec_cells) {
      list_cells.add(cell);
    }
    d_cells = list_cells;
  }

  if (vec_nodes.empty()) {
    d_nodes = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  } else {
    ListOfCellsIterator list_nodes(vec_nodes.size());
    for (const auto& node : vec_nodes) {
      list_nodes.add(node);
    }
    d_nodes = list_nodes;
  }
}

void
BCGeomBase::determineInteriorBndIteratorLimits(const Patch::FaceType face,
                                               const Patch* patch)
{
  int pdir        = 0;
  int plusOrMinus = 1;
  switch (face) {
    case Uintah::Patch::xminus:
      plusOrMinus = -1;
      pdir        = 0;
      break;
    case Uintah::Patch::xplus:
      plusOrMinus = 1;
      pdir        = 0;
      break;
    case Uintah::Patch::yminus:
      plusOrMinus = -1;
      pdir        = 1;
      break;
    case Uintah::Patch::yplus:
      plusOrMinus = 1;
      pdir        = 1;
      break;
    case Uintah::Patch::zminus:
      plusOrMinus = -1;
      pdir        = 2;
      break;
    case Uintah::Patch::zplus:
      plusOrMinus = 1;
      pdir        = 2;
      break;
    default:
      break;
  }

  // This will be aligned with the nearest layer of nodes to this geometry
  // object. When constructing this geometry object, Uintah will move its center
  // with the nearest node.
  Point origin = getOrigin();

  // find the nearest cell to the plus or minus side of this interior boundary.
  Point nearestCell = origin + plusOrMinus * patch->dCell() / 2.0;

  Uintah::IntVector nearestCellIdx = patch->getCellIndex(nearestCell);
  // Find how far away from the patch boundary the layer of nearest cells is.
  // The purpose here is to find the layer of cells that are near this interior
  // boundary and figure out which cells straddle the interior boundary.
  Uintah::IntVector offset;
  if (plusOrMinus == -1) { // minus side
    offset = nearestCellIdx - patch->getExtraCellLowIndex();
  } else {
    offset = patch->getExtraCellHighIndex() - nearestCellIdx;
    offset[pdir] -= 1;
  }
  Uintah::IntVector lpts, hpts;
  // Get the layer of cells that are on either the plus or minus side of this
  // interior boundary
  patch->getFaceCells(face, -offset[pdir], lpts, hpts);
  std::vector<Point> test_pts;

  // Now that we have a layer of cells near the interior boundary, if we move
  // 1/2*dx in the plus or minus direction, then we should hit points that could
  // possibly be inside this interior geometry. This will essentially create a
  // layer of points on the same plane that this interior geometry lives on.
  Vector halfdx = patch->dCell() / 2.0;
  for (CellIterator candidatePoints(lpts, hpts); !candidatePoints.done();
       candidatePoints++) {
    Point p = patch->getCellPosition(*candidatePoints);
    p(pdir) -= plusOrMinus *
               halfdx[pdir]; // for a plus face, pick cells on the plus side
    test_pts.push_back(p);
  }

  GridIterator cells(lpts, hpts);
  Iterator cell_itr(cells);
  auto pts = test_pts.begin();

  std::vector<IntVector> vec_cells;

  cell_itr.reset();
  for (; pts != test_pts.end(); ++pts, ++cell_itr) {
    if (patch->containsCell(*cell_itr)) {
      // if ( (*pts)(pdir) == origin(pdir)) {
      if (inside(*pts)) {
        vec_cells.push_back(*cell_itr);
      }
      //}
    }
  }

  if (vec_cells.empty()) {
    d_cells = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  } else {
    ListOfCellsIterator list_cells(vec_cells.size());
    for (const auto& cell : vec_cells) {
      list_cells.add(cell);
    }
    d_cells = list_cells;
  }

  // Now for nodes...
  IntVector ln, hn;
  patch->getFaceNodes(face, -offset[pdir], ln, hn);
  GridIterator nodes(ln, hn);
  Iterator node_itr(nodes);
  std::vector<IntVector> vec_nodes;
  for (node_itr.reset(); !node_itr.done(); node_itr++) {
    Point p = patch->getLevel()->getNodePosition(*node_itr);
    if (inside(p)) {
      vec_nodes.push_back(*node_itr);
    }
  }
  if (vec_nodes.empty()) {
    d_nodes = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  } else {
    ListOfCellsIterator list_nodes(vec_nodes.size());
    for (const auto& node : vec_nodes) {
      list_nodes.add(node);
    }
    d_nodes = list_nodes;
  }

  d_surfaceArea = d_cells.size() * patch->cellArea(face);
}

void
BCGeomBase::printLimits() const
{
  
  std::cout << std::endl;
  std::cout << "d_cells = " << d_cells.begin() << " " << d_cells.end()
            << std::endl;
  std::cout << "d_nodes = " << d_nodes.begin() << " " << d_nodes.end()
            << std::endl;
  std::cout << std::endl;
}

} // namespace Uintah