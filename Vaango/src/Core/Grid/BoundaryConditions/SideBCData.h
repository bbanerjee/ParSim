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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_SideBCData_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_SideBCData_H__

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BCGeomBase.h>
#include <Core/Grid/Variables/GridIterator.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <vector>

namespace Uintah {

/*!

\class SideBCData

\ brief Defines a boundary condition geometry for the entire side of the
domain.

\author John A. Schmidt \n
Department of Mechanical Engineering \n
University of Utah \n
Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n\n

*/

class SideBCData : public BCGeomBase
{

public:
  /// Constructor
  SideBCData();

  /// Assignment Operator
  SideBCData&
  operator=(const SideBCData& bc);

  /// Destructor
  virtual ~SideBCData();

  virtual bool
  operator==(const BCGeomBase&) const;

  /// Clone the boundary condition geometry -- allocates memory.
  SideBCData*
  clone();

  /// Get the boundary condition data
  void
  getBCData(BCData& bc) const;

  /// Add the boundary condition data
  void
  addBCData(BCData& bc);

  /// Add the old boundary condition data -- no longer used.
  void
  addBC(BoundCondBaseP bc);

  /// Add boundary condition within a scheduled task.
  void
  sudoAddBC(BoundCondBaseP& bc);

  /// Determines if a point is inside -- always returns true.
  bool
  inside(const Point& p) const;

  /// Print out the boundary condition geometry type.
  virtual void
  print();

  /// Determine the cell and node centered iterators
  virtual void
  determineIteratorLimits(Patch::FaceType face,
                          const Patch* patch,
                          std::vector<Point>& test_pts);

private:
  BCData d_bc;
};

} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARYCONDITIONS_SideBCData_H__
