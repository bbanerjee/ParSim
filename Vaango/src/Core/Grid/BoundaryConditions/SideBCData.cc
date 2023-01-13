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

#include <Core/Grid/BoundaryConditions/SideBCData.h>

#include <Core/Geometry/Point.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/Box.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>
#include <iostream>

namespace {

// Usage: export SCI_DEBUG="BC_dbg:+"
Uintah::Dout bc_dbg{ "BC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Side BC debug info",
                     false };

}

namespace Uintah {

SideBCData::SideBCData()
{
  d_cells  = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  d_nodes  = GridIterator(IntVector(0, 0, 0), IntVector(0, 0, 0));
  d_bcname = "NotSet";
}

SideBCData::~SideBCData() {}

bool
SideBCData::operator==(const BCGeomBase& rhs) const
{
  const SideBCData* p_rhs = dynamic_cast<const SideBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    return true;
  }
}

SideBCData*
SideBCData::clone()
{
  return scinew SideBCData(*this);
}

void
SideBCData::addBCData(BCData& bc)
{
  d_bc = bc;
}

void
SideBCData::addBC(BoundCondBaseP bc)
{
  d_bc.setBCValues(bc);
}

void
SideBCData::sudoAddBC(BoundCondBaseP& bc)
{
  d_bc.setBCValues(bc);
}

void
SideBCData::getBCData(BCData& bc) const
{
  bc = d_bc;
}

bool
SideBCData::inside(const Point& p) const
{
  return true;
}

void
SideBCData::print()
{
  DOUT(bc_dbg, "Geometry type = " << typeid(this).name());
  d_bc.print();
}

void
SideBCData::determineIteratorLimits(Patch::FaceType face,
                                    const Patch* patch,
                                    std::vector<Point>& test_pts)
{
  DOUT(bc_dbg, "SideBC determineIteratorLimits() " << patch->getFaceName(face));

  IntVector l, h;
  patch->getFaceCells(face, 0, l, h);
  d_cells = GridIterator(l, h);

  DOUT(bc_dbg, "d_cells->begin() = " << d_cells.begin() << " d_cells->end() = " 
       << d_cells.end());

  IntVector ln, hn;
  patch->getFaceNodes(face, 0, ln, hn);
  d_nodes = GridIterator(ln, hn);

  DOUT(bc_dbg, "d_nodes->begin() = " << d_nodes.begin() << " d_nodes->end() = " 
       << d_nodes.end());

  //  Iterator iii(d_cells);
  DOUT(bc_dbg, "Iterator output . . . ");
  for (Iterator ii(d_cells); !ii.done(); ii++) {
    DOUT(bc_dbg, ii);
  }
}

} // namespace Uintah
