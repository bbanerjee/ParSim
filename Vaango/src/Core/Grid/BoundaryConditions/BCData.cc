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

#include <Core/Grid/BoundaryConditions/BCData.h>

#include <Core/Grid/BoundaryConditions/BoundCondBase.h>
#include <Core/Util/DOUT.hpp>

#include <algorithm>
#include <iostream>
#include <typeinfo> // for typeid

namespace {

// Usage: export SCI_DEBUG="BCData_dbg:+"
Uintah::Dout BCData_dbg{ "BCData_dbg",
                         "Grid_BoundaryConditions",
                         "Grid BC data debug info",
                         false };

} // namespace

namespace Uintah {

BCData::BCData(const BCData& rhs)
{
  for (const auto& bc : rhs.d_BCData) {
    d_BCData.push_back(bc->clone());
  }
}

auto
BCData::operator=(const BCData& rhs) -> BCData&
{
  if (this == &rhs) {
    return *this;
  }

  // Delete the lhs
  d_BCData.clear();

  // Copy the rhs to the lhs
  for (const auto& bc : rhs.d_BCData) {
    d_BCData.push_back(bc->clone());
  }

  return *this;
}

auto
BCData::operator==(const BCData& rhs) -> bool
{
  if (d_BCData.size() != rhs.d_BCData.size()) {
    return false;
  }

  for (const auto& bc : d_BCData) {
    if (rhs.find(bc->getBCVariable()) == false) {
      return false;
    }
  }

  return true;
}

auto
BCData::operator<(const BCData& rhs) const -> bool
{
  if (d_BCData.size() < rhs.d_BCData.size()) {
    return true;
  } else {
    return false;
  }
}

void
BCData::setBCValues(BoundCondBaseSP bc)
{
  if (!find(bc->getBCVariable())) {
    d_BCData.push_back(bc->clone());
  }
}

auto
BCData::cloneBCValues(const std::string& var_name) const
  -> const BoundCondBaseSP
{
  // The default location for BCs defined for all materials is mat_id = -1.
  // Need to first check the actual mat_id specified.  If this is not found,
  // then will check mat_id = -1 case.  If it isn't found, then return 0.

  for (const auto& bc : d_BCData) {
    if (bc->getBCVariable() == var_name) {
      return bc->clone();
    }
  }
  return nullptr;
}

auto
BCData::getBCValues(const std::string& var_name) const -> const BoundCondBase*
{
  // The default location for BCs defined for all materials is mat_id = -1.
  // Need to first check the actual mat_id specified.  If this is not found,
  // then will check mat_id = -1 case.  If it isn't found, then return 0.
  for (const auto& bc : d_BCData) {
    if (bc->getBCVariable() == var_name) {
      return bc.get();
    }
  }
  return nullptr;
}

auto
BCData::getBCData() const -> const std::vector<BoundCondBaseSP>&
{
  return d_BCData;
}

auto
BCData::find(const std::string& var_name) const -> bool
{
  for (const auto& bc : d_BCData) {
    if (bc->getBCVariable() == var_name) {
      return true;
    }
  }
  return false;
}

auto
BCData::find(const std::string& bc_type, const std::string& bc_variable) const
  -> bool
{
  const BoundCondBase* bc = getBCValues(bc_variable);
  if (bc) {
    if (bc->getBCType() == bc_type) {
      return true;
    }
  }
  return false;
}

void
BCData::combine(BCData& from)
{
  for (const auto& bc : from.d_BCData) {
    setBCValues(bc);
  }
}

void
BCData::print() const
{
  DOUT(BCData_dbg, "size of d_BCData = " << d_BCData.size());
  for (const auto& bc : d_BCData) {
    DOUT(BCData_dbg,
         "BC = " << bc->getBCVariable() << " type = " << bc->getBCType());
  }
}

} // namespace Uintah