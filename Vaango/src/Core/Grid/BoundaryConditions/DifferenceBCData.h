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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_DifferenceBCData_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_DifferenceBCData_H__

#include <Core/Grid/BoundaryConditions/BCGeomBase.h>

#include <Core/Geometry/Vector.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <vector>

namespace Uintah {
/*!

\class DifferenceBCData

\ brief Stores the difference of two different boundary condition geometries.

\author John A. Schmidt \n
Department of Mechanical Engineering \n
University of Utah \n
Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n\n

*/

class DifferenceBCData : public BCGeomBase
{
public:
  /// Default constructor
  DifferenceBCData();

  /// Copy constructor
  DifferenceBCData(const DifferenceBCData& rhs);

  /// Assignment operator
  auto
  operator=(const DifferenceBCData& bc) -> DifferenceBCData&;

  /// Constructor with two boundary condition geometries.  The second
  /// argument is subtracted from the first argument.
  DifferenceBCData(std::shared_ptr<BCGeomBase> p1, std::shared_ptr<BCGeomBase> p2);

  /// Destructor
  ~DifferenceBCData() override;

  bool
  operator==(const BCGeomBase&) const override;

  /// Clone the boundary condition geometry -- allocates memory.
  std::shared_ptr<BCGeomBase>
  clone() override;

  /// Get the boundary condition data
  void
  getBCData(BCData& bc) const override;

  /// Add the boundary condition data -- no longer used.
  void
  addBCData(BCData& bc) override;

  /// Add the old boundary condition data -- no longer used.
  void
  addBC(BoundCondBaseSP bc) override;

  void
  sudoAddBC(BoundCondBaseSP& bc) override;

  /// Determine if a point is inside the object.
  [[nodiscard]] bool
  inside(const Point& p) const override;

  /// Print out the boundary condition geometry types.
  void
  print() override;

  /// Determine the cell and node boundary iterators.
  void
  determineIteratorLimits(Patch::FaceType face,
                          const Patch* patch,
                          std::vector<Point>& test_pts) override;

private:
  std::shared_ptr<BCGeomBase> left;
  std::shared_ptr<BCGeomBase> right;

  friend class BCReader;
};

} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARYCONDITIONS_DifferenceBCData_H__
