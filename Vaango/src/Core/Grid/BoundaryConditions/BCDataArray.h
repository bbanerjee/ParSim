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

#ifndef __CORE_GRID_BOUNDARY_CONDITIONS_BCDataArray_H__
#define __CORE_GRID_BOUNDARY_CONDITIONS_BCDataArray_H__

#include <Core/Geometry/Vector.h>
#include <Core/Grid/BoundaryConditions/BCData.h>
#include <Core/Grid/BoundaryConditions/BCGeomBase.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/Iterator.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <map>
#include <vector>

namespace Uintah {

/*!

  \class BCDataArray

  \brief Supporting class for holding the basic boundary condition geometry
  classes such as side, circle, rectangle, union, and difference.

  The boundary conditions may be applied to \b all materials.  In this case,
  we use the mat_id of -1.  Other materials are stored according to the mat_id
  specified such as 0, 1, 2, etc.

  \author John A. Schmidt \n
  Department of Mechanical Engineering \n
  University of Utah \n
  Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n\n

*/

class BCDataArray
{
public:
  /// Constructor
  BCDataArray() = default;

  /// Destructor
  ~BCDataArray();

  /// Copy constructor
  BCDataArray(const BCDataArray& bc);

  /// Assignment operator
  auto
  operator=(const BCDataArray& bc) -> BCDataArray&;

  /// Make a clone of self.
  auto
  clone() -> std::shared_ptr<BCDataArray>;

  /// Get the boundary condition data for a given material and a given
  /// type for a given child.
  [[nodiscard]] auto
  getBoundCondData(int mat_id, const std::string& type, int ichild) const
    -> const BoundCondBaseSP;

  auto
  checkForBoundCondData(int& mat_id, const std::string& type, int ichild)
    -> bool;

  /// Determine the iterator limits.
  void
  determineIteratorLimits(Patch::FaceType face, const Patch* patch);

  /*
        \author Tony Saad
        \date   September 2014
        \brief  Determine the iterator points associated with interior
     boundaries.
        */
  void
  determineInteriorBndIteratorLimits(Patch::FaceType face, const Patch* patch);

  /// Add boundary condition data
  void
  addBCData(int mat_id, std::shared_ptr<BCGeomBase> bc);

  /// Combine the duplicate BCGeometryTypes into a single BCGeometryType
  void
  combineBCGeometryTypes(int mat_id);
  void
  combineBCGeometryTypes_NEW(int mat_id);

  /// Get the cell centered face iterator for the ith face child on mat_id.
  void
  getCellFaceIterator(int mat_id, Iterator& b_ptr, int ichild) const;

  /// Get the node centered face iterator for the ith face child on mat_id.
  void
  getNodeFaceIterator(int mat_id, Iterator& b_ptr, int ichild) const;

  /// Return the number of children in the vector<BCGeomBase*>.
  [[nodiscard]] auto
  getNumberChildren(int mat_id) const -> int;

  /// Get the ith child.
  [[nodiscard]] auto
  getChild(int mat_id, int ichild) const -> std::shared_ptr<BCGeomBase>;

  /// Print out the various boundary condition geometry types.
  void
  print() const;

  auto
  getBCGeom(int matl_index) -> std::vector<std::shared_ptr<BCGeomBase>>
  {
    return d_BCDataArray[matl_index];
  }

  /// The map is for the mat_id.  -1 is for mat_id = "all", 0, for
  /// mat_id = "0", etc.
  using bcDataArrayType =
    std::map<int, std::vector<std::shared_ptr<BCGeomBase>>>;

private:
  bcDataArrayType d_BCDataArray;

  friend class Patch;
  friend class BoundCondReader;
};

// Used for inserting IntVectors into a set.  The standard < operator
// for IntVectors is too restrictive.

/// Sorts along the x axis
struct ltiv_x
{
  auto
  operator()(const IntVector& i1, const IntVector& i2) -> bool
  {
    if (i2.x() < i1.x()) {
      return false;
    }
    if (i1.y() < i2.y()) {
      return true;
    }
    if (i1.z() < i2.z()) {
      return true;
    }

    return false;
  }
};

/// Sorts along the y axis
struct ltiv_y
{
  auto
  operator()(const IntVector& i1, const IntVector& i2) -> bool
  {
    if (i2.y() < i1.y()) {
      return false;
    }
    if (i1.z() < i2.z()) {
      return true;
    }
    if (i1.x() < i2.x()) {
      return true;
    }

    return false;
  }
};

/// Sorts along the z axis
struct ltiv_z
{
  auto
  operator()(const IntVector& i1, const IntVector& i2) -> bool
  {
    if (i2.z() < i1.z()) {
      return false;
    }
    if (i1.y() < i2.y()) {
      return true;
    }
    if (i1.x() < i2.x()) {
      return true;
    }

    return false;
  }
};

/// A less restrictive < operator rather than the built-in one for
/// IntVectors.
struct ltiv_xyz
{
  auto
  operator()(const IntVector& i1, const IntVector& i2) -> bool
  {
    if (i1.x() < i2.x() && i1.y() < i2.y() && i1.z() < i2.z()) {
      return true;
    } else {
      return false;
    }
  }
};

} // End namespace Uintah

#endif //__CORE_GRID_BOUNDARY_CONDITIONS_BCDataArray_H__
