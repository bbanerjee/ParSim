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

#ifndef __CORE_GRID_BOUNDARYCONDITIONS_BOUNDARYCONDITIONREADER_H__
#define __CORE_GRID_BOUNDARYCONDITIONS_BOUNDARYCONDITIONREADER_H__

#include <Core/Geometry/Point.h>
#include <Core/Grid/BoundaryConditions/BCData.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Grid/BoundaryConditions/BCGeomBase.h>
#include <Core/Grid/Patch.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <map>
#include <string>

namespace Uintah {

/*!

  \class BoundCondReader

  \brief Reads in the boundary conditions and stores them for later processing
  by Patch.cc.

  Reads in the boundary conditions for a given face and a given material.
  Multiple material bc specifications may be combined within a single
  \<Face\> \</Face\>.

  The boundary condition specification consist of two components, a
  geometrical component and a type-value pair component.  The geometrical
  component describes where on a patch's face the boundary condition is
  applied.  The side specification will apply the bc's type-value pair to
  the entire face.  The circle specification will apply the bc to a circular
  region on a face.  The rectangle specification will apply the bc to a
  rectangular region on a face.  There may be multiple specifications of
  circles and rectangles on a give face.  There must be at least one side
  specification for each face.  The type-value pair component are described
  in the BoundCondBase class.

  A significant component of the BoundCondReader class is to take the
  geometrical specification and arrange the various instances of side,
  circle, and rectangle in such a manner as to ensure that the entire face
  has an bc value specified.  To accomplish this, auxillary classes
  UnionBCData, and DifferenceBCData are used to manage the union of multiple
  instances of rectangle and circle bcs for a given side and the difference
  of that collection with the base side bc.  Unions and difference bcs are
  not available to the user.

  Currently, this class supports the old way of specifying the geometric
  component and the new way.  The old way only allowed a side specification.
  The new way allows for circles, rectangles and a side specification
  for a given face.  At some point, the old way infrastructure will be go
  away.


  \author John A. Schmidt \n
  Department of Mechanical Engineering \n
  University of Utah \n
  Center for the Simulation of Accidental Fires and Explosions (C-SAFE) \n\n

*/

class BoundCondReader
{
public:
  /// Constructor
  BoundCondReader() = default;

  /// Destructor
  ~BoundCondReader() = default;

  [[nodiscard]] auto
  is_on_face(const int dir,
             const Point p,
             const std::vector<Point>& points) const -> bool;

  [[nodiscard]] auto
  isPtOnFace(const int dir,
             const int plusMinusFaces,
             const Point pt,
             const std::vector<Point>& grid_LoPts,
             const std::vector<Point>& grid_HiPts) const -> bool;

  void
  whichPatchFace(const std::string fc,
                 Patch::FaceType& face_side,
                 int& plusMinusFaces,
                 int& p_dir) const;

  /// Read in the boundary conditions given a problem specification, ps.
  /// Each face is read in and processed using the function
  /// createBoundaryConditionFace() which indicates whether the boundary
  /// condition should be applied on the side, circle, or rectangle region.
  /// Then the individual boundary conditions such as Pressure, Density,
  /// etc. are processed.  The individual boundary conditions for a given
  /// face may be for several different materials.   Boundary conditions
  /// are then separated out by material ids.  Following separation, the
  /// boundary conditions are then combined (combineBCS()) so that any
  /// circles and rectangles specified on a face are combined into a union.
  /// The union is the subtracted out from the side case using a difference
  /// boundary condition.
  void
  read(ProblemSpecP& ps, const ProblemSpecP& grid_ps, const LevelP level);

  /// Read in the geometric tags: side, circle, and rectangle.  Performs
  /// error checking if the tag is not present or if the circle and rectangle
  /// tags are not specified correctly.
  auto
  createBoundaryConditionFace(ProblemSpecP& ps,
                              const ProblemSpecP& grid_ps,
                              Patch::FaceType& face_side)
    -> std::shared_ptr<BCGeomBase>;

  /*!
  \author  Tony Saad
  \date    September 2014
  \brief   Creates interior boundaries. Currently supported geometries are:
  rectangle, circle, ellipse, and annulus.

  For a given grid resolution, the interior boundary is moved to the closest
  face. If the interior boundary coincides with a cell center, then it is moved
  to the face side (minus/plus) that is specified through the input file. The
  face side (minus/plus) determines which cell iterator is returned. For a minus
  boundary, the cells on the minus side are returned. For a plus boundary, the
  cells on the plus side are returned.

  \todo    Handle unions and differences.
  */
  auto
  createInteriorBndBoundaryConditionFace(ProblemSpecP& ps,
                                         const ProblemSpecP& grid_ps,
                                         Patch::FaceType& face_side,
                                         const LevelP level)
    -> std::shared_ptr<BCGeomBase>;

  /// Combine the boundary conditions for a given face into union and
  /// difference operations for the face.  Multiple circles and rectangles
  /// are stored in a union.  The resultant union is then subtracted from
  /// the side and stored as a difference bc.  This operation only happens
  /// if there are more than one bc specified for a given face.
  void
  combineBCS();

  void
  bulletProofing();

  ///
  auto
  compareBCData(std::shared_ptr<BCGeomBase> b1, std::shared_ptr<BCGeomBase> b2)
    -> bool;

  /// not used
  auto
  getBCDataArray(Patch::FaceType& face) const -> const BCDataArray;

private:
  friend class Level;
  friend class Patch;

  std::map<Patch::FaceType, std::shared_ptr<BCDataArray>> d_BCReaderData;
  std::map<Patch::FaceType, std::shared_ptr<BCDataArray>>
    d_interiorBndBCReaderData;

  void
  readDomainBCs(ProblemSpecP& ps, const ProblemSpecP& grid_ps);
  void
  readInteriorBndBCs(ProblemSpecP& ps,
                     const ProblemSpecP& grid_ps,
                     const LevelP level);

  auto
  createSideBC(const std::map<std::string, std::string>& values,
               Patch::FaceType& face_side) const -> std::shared_ptr<BCGeomBase>;

  auto
  createCircleBC(const std::map<std::string, std::string>& values,
                 const std::vector<Point>& grid_LoPts,
                 const std::vector<Point>& grid_HiPts,
                 Patch::FaceType& face_side) const
    -> std::shared_ptr<BCGeomBase>;

  auto
  createAnnulusBC(const std::map<std::string, std::string>& values,
                  const std::vector<Point>& grid_LoPts,
                  const std::vector<Point>& grid_HiPts,
                  Patch::FaceType& face_side) const
    -> std::shared_ptr<BCGeomBase>;

  auto
  createEllipseBC(const std::map<std::string, std::string>& values,
                  const std::vector<Point>& grid_LoPts,
                  const std::vector<Point>& grid_HiPts,
                  Patch::FaceType& face_side) const
    -> std::shared_ptr<BCGeomBase>;

  auto
  createRectangleBC(const std::map<std::string, std::string>& values,
                    const std::vector<Point>& grid_LoPts,
                    const std::vector<Point>& grid_HiPts,
                    Patch::FaceType& face_side) const
    -> std::shared_ptr<BCGeomBase>;

  auto
  createRectangulusBC(const std::map<std::string, std::string>& values,
                      const std::vector<Point>& grid_LoPts,
                      const std::vector<Point>& grid_HiPts,
                      Patch::FaceType& face_side) const
    -> std::shared_ptr<BCGeomBase>;
};

} // End namespace Uintah

namespace Uintah::BCReaderUtils {

void
print(std::shared_ptr<BCGeomBase> p);

auto
moveToClosestNode(const LevelP level,
                  const int facedir,
                  const int plusMinusFaces,
                  const Point& p0) -> Point;

} // namespace Uintah::BCReaderUtils

#endif //__CORE_GRID_BOUNDARYCONDITIONS_BOUNDARYCONDITIONREADER_H__
