/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_BCData_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_BCData_H__

#include <Core/Grid/BoundaryConditions/BoundCondBase.h>
#include <Core/Grid/BoundaryConditions/BoundCondBaseP.h>
#include <Core/Util/Handle.h>

#include <map>
#include <string>
#include <vector>

namespace Uintah {

class BCData
{
public:
  BCData() = default;

  ~BCData() = default;

  BCData(const BCData&);

  auto
  operator=(const BCData&) -> BCData&;

  auto
  operator==(const BCData&) -> bool;

  auto
  operator<(const BCData&) const -> bool;

  void
  setBCValues(BoundCondBaseSP bc);

  /**
   * \brief Clones the boundary conditions associated with this boundary.
   *   This is a deep copy and the user is responsible for deleting the newly
   *   created variable.
   */
  [[nodiscard]] auto
  cloneBCValues(const std::string& type) const -> const BoundCondBaseSP;

  /**
   *  \author Tony Saad
   *  \date   August 2014
   *  \brief Returns a const pointer to the boundary conditions associated with
   * this boundary. This is a lightweight access to the BCValues and does not
   * perform a deep copy.
   */
  [[nodiscard]] auto
  getBCValues(const std::string& type) const -> const BoundCondBase*;

  /**
       *  \author Tony Saad
       *  \date   August 29, 2013
       *  \brief  Returns a reference to the vector of boundary conditions
     (BoundCondBase) attached to this BCData. Recall that BCDataArray stores a
     collection of BCGeometries. Each BCGeomBase has BCData associated with it.
     The BCData stores all the boundary specification via BoundCond and
     BoundCondBase. This boundary specification corresponds to all the spec
     within a given <Face> xml tag in the input file.
       */
  [[nodiscard]] auto
  getBCData() const -> const std::vector<BoundCondBaseSP>&;

  void
  print() const;

  [[nodiscard]] auto
  find(const std::string& type) const -> bool;

  [[nodiscard]] auto
  find(const std::string& bc_type, const std::string& bc_variable) const -> bool;

  void
  combine(BCData& from);

private:
  // The map is for the name of the
  // bc type and then the actual bc data, i.e.
  // "Velocity", VelocityBoundCond
  std::vector<BoundCondBaseSP> d_BCData;
};
} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARYCONDITIONS_BCData_H__
