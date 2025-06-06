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
#include <Core/Grid/BoundaryConditions/RectangleBCData.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>
#include <iostream>

namespace {

// Usage: export SCI_DEBUG="BC_dbg:+"
Uintah::Dout bc_dbg{ "BC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Rectangle BC debug info",
                     false };

}

namespace Uintah {

RectangleBCData::RectangleBCData()
  : BCGeomBase()
{
}

RectangleBCData::RectangleBCData(Point& low, Point& up)
  : BCGeomBase()
  , d_min(low)
  , d_max(up)
{
}

RectangleBCData::~RectangleBCData() = default;

auto
RectangleBCData::operator==(const BCGeomBase& rhs) const -> bool
{
  const auto* p_rhs = dynamic_cast<const RectangleBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    return (this->d_min == p_rhs->d_min) && (this->d_max == p_rhs->d_max);
  }
}

std::shared_ptr<BCGeomBase>
RectangleBCData::clone()
{
  return std::make_shared<RectangleBCData>(*this);
}

void
RectangleBCData::addBCData(BCData& bc)
{
  d_bc = bc;
}

void
RectangleBCData::addBC(BoundCondBaseSP bc)
{
  d_bc.setBCValues(bc);
}

void
RectangleBCData::sudoAddBC(BoundCondBaseSP& bc)
{
  d_bc.setBCValues(bc);
}

void
RectangleBCData::getBCData(BCData& bc) const
{
  bc = d_bc;
}

auto
RectangleBCData::inside(const Point& p) const -> bool
{
  if (d_min.x() == d_max.x()) {
    if (p.y() <= d_max.y() && p.y() >= d_min.y() && p.z() <= d_max.z() &&
        p.z() >= d_min.z()) {
      return true;
    } else {
      return false;
    }
  }

  else if (d_min.y() == d_max.y()) {
    if (p.x() <= d_max.x() && p.x() >= d_min.x() && p.z() <= d_max.z() &&
        p.z() >= d_min.z()) {
      return true;
    } else {
      return false;
    }
  }

  else if (d_min.z() == d_max.z()) {
    if (p.y() <= d_max.y() && p.y() >= d_min.y() && p.x() <= d_max.x() &&
        p.x() >= d_min.x()) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

void
RectangleBCData::print()
{
  DOUT(bc_dbg, "Geometry type = " << typeid(this).name());
  d_bc.print();
}

void
RectangleBCData::determineIteratorLimits(Patch::FaceType face,
                                         const Patch* patch,
                                         std::vector<Point>& test_pts)
{
#if 0
  std::cout << "RectangleBC determineIteratorLimits()" << std::endl;
#endif

  BCGeomBase::determineIteratorLimits(face, patch, test_pts);
}

} // namespace Uintah