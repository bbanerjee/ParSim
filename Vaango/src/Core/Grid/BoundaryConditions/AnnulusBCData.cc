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
#include <Core/Grid/BoundaryConditions/AnnulusBCData.h>
#include <Core/Grid/BoundaryConditions/BoundCondFactory.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/DOUT.hpp>
#include <iostream>

namespace {

Uintah::Dout bc_dbg{ "AnnulusBC_dbg",
                     "Grid_BoundaryConditions",
                     "Grid Annulus BC debug info",
                     false };

}

namespace Uintah {

AnnulusBCData::AnnulusBCData()
  : BCGeomBase()
{
}

AnnulusBCData::AnnulusBCData(Point& p, double inRadius, double outRadius)
  : BCGeomBase()
  , d_innerRadius(inRadius)
  , d_outerRadius(outRadius)
  , d_origin(p)
{
}

AnnulusBCData::~AnnulusBCData() = default;

auto
AnnulusBCData::operator==(const BCGeomBase& rhs) const -> bool
{
  const auto* p_rhs = dynamic_cast<const AnnulusBCData*>(&rhs);

  if (p_rhs == nullptr) {
    return false;
  } else {
    return (this->d_innerRadius == p_rhs->d_innerRadius) &&
           (this->d_outerRadius == p_rhs->d_outerRadius) &&
           (this->d_origin == p_rhs->d_origin);
  }
}

std::shared_ptr<BCGeomBase>
AnnulusBCData::clone()
{
  return std::make_shared<AnnulusBCData>(*this);
}

void
AnnulusBCData::addBCData(BCData& bc)
{
  d_bc = bc;
}

void
AnnulusBCData::addBC(BoundCondBaseSP bc)
{
  d_bc.setBCValues(bc);
}

void
AnnulusBCData::sudoAddBC(BoundCondBaseSP& bc)
{
  d_bc.setBCValues(bc);
}

void
AnnulusBCData::getBCData(BCData& bc) const
{
  bc = d_bc;
}

auto
AnnulusBCData::inside(const Point& p) const -> bool
{
  Vector diff = p - d_origin;

  bool inside_outer  = false;
  bool outside_inner = false;

  if (diff.length() > d_outerRadius) {
    inside_outer = false;
  } else {
    inside_outer = true;
  }

  if (diff.length() > d_innerRadius) {
    outside_inner = true;
  } else {
    outside_inner = false;
  }

  return (inside_outer && outside_inner);
}

void
AnnulusBCData::print()
{
  DOUT(bc_dbg, "Geometry type = " << typeid(this).name());
  d_bc.print();
}

void
AnnulusBCData::determineIteratorLimits(Patch::FaceType face,
                                       const Patch* patch,
                                       std::vector<Point>& test_pts)
{
  DOUT(bc_dbg, "Annulus determineIteratorLimits()");
  BCGeomBase::determineIteratorLimits(face, patch, test_pts);
}

} // namespace Uintah