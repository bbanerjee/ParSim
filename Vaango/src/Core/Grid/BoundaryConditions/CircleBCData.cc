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

#include <Core/Grid/BoundaryConditions/CircleBCData.h>

#include <Core/Util/DOUT.hpp>
#include <iostream>

namespace {

// Usage: export SCI_DEBUG="CircleBC_dbg:+"
Uintah::Dout bc_dbg{ "CircleBC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Circle BC debug info",
                     false };

}
namespace Uintah {

CircleBCData::CircleBCData()
  : BCGeomBase()
{
}

CircleBCData::CircleBCData(Point& p, double radius)
  : BCGeomBase()
  , d_radius(radius)
  , d_origin(p)
{
}

CircleBCData::~CircleBCData() = default;

auto
CircleBCData::operator==(const BCGeomBase& rhs) const -> bool
{

  const auto* p_rhs = dynamic_cast<const CircleBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    return (this->d_radius == p_rhs->d_radius) &&
           (this->d_origin == p_rhs->d_origin);
  }
}

std::shared_ptr<BCGeomBase>
CircleBCData::clone()
{
  return std::make_shared<CircleBCData>(*this);
}

void
CircleBCData::addBCData(BCData& bc)
{
  d_bc = bc;
}

void
CircleBCData::addBC(BoundCondBaseSP bc)
{
  d_bc.setBCValues(bc);
}

void CircleBCData::sudoAddBC(BoundCondBaseSP& bc) 
{
  d_bc.setBCValues(bc);
}

void
CircleBCData::getBCData(BCData& bc) const
{
  bc = d_bc;
}

auto
CircleBCData::inside(const Point& p) const -> bool
{
  Vector diff = p - d_origin;
  if (diff.length() > d_radius) {
    return false;
  } else {
    return true;
  }
}

void
CircleBCData::print()
{
  DOUT(bc_dbg, "Geometry type = " << typeid(this).name());
  d_bc.print();
}

void
CircleBCData::determineIteratorLimits(Patch::FaceType face,
                                      const Patch* patch,
                                      std::vector<Point>& test_pts)
{
  DOUT(bc_dbg, "Circle determineIteratorLimits() " << patch->getFaceName(face));

  BCGeomBase::determineIteratorLimits(face, patch, test_pts);
}

} // namespace Uintah