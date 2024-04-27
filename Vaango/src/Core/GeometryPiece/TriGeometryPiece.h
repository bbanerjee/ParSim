/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2014-2023 Biswajit Banerjee
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

#ifndef __VAANGO_CORE_GEOMPIECE_TRI_GEOMETRYPIECE_H__
#define __VAANGO_CORE_GEOMPIECE_TRI_GEOMETRYPIECE_H__

#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/GeometryPiece/UniformGrid.h>
#include <Core/Grid/Box.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Plane.h>
#include <Core/Geometry/Point.h>

#include <memory>
#include <vector>

namespace Uintah {

/**************************************

CLASS
   TriGeometryPiece

   Creates a triangulated surface piece from the xml input file description.

GENERAL INFORMATION

   TriGeometryPiece.h

   John A. Schmidt
   Department of Mechanical Engineering
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)

KEYWORDS
   TriGeometryPiece BoundingBox inside

DESCRIPTION
   Creates a triangulated surface piece from the xml input file description.
   Requires one input: file name (convetion use suffix .dat).
   There are methods for checking if a point is inside the surface
   and also for determining the bounding box for the surface.
   The input form looks like this:
       <tri>
         <file>surface.dat</file>
       </tri>

****************************************/

class TriGeometryPiece : public GeometryPiece
{
public:
  inline static const std::string TYPE_NAME{ "tri" };

public:
  TriGeometryPiece(ProblemSpecP&, const std::string& input_ups_dir = "");
  virtual ~TriGeometryPiece() = default;

  TriGeometryPiece(const TriGeometryPiece&);

  TriGeometryPiece&
  operator=(const TriGeometryPiece&);

  virtual std::string
  getType() const
  {
    return TYPE_NAME;
  }

  virtual GeometryPieceP
  clone() const;

  // Find if a point is inside the triangulated surface
  virtual bool
  inside(const Point& p) const override;

  bool
  inside(const Point& p, bool use_x_crossing) const;

  bool
  inside(const Point& p, int& crossings) const;

  // A relatively newer version of the inside test that uses three
  // nearly orthogonal rays and counts intersections in each of those
  // directions
  bool
  inside(const Point& p, int& crossings, bool allDirections) const;

  virtual Box
  getBoundingBox() const override;

  void
  scale(const double factor);

  double
  surfaceArea() const;

  inline int
  getNumIntersections(const Point& start,
                      const Point& end,
                      double& min_distance)
  {
    int intersections = 0;

    d_grid->countIntersections(start, end, intersections, min_distance);
    return intersections;
  }

  inline int
  getNumIntersections(const Point& start)
  {
    int intersections = 0;
    d_grid->countIntersections(start, intersections);
    return intersections;
  }

  inline std::vector<IntVector>
  getTriangles()
  {
    return d_triangles;
  }

  inline std::vector<Point>
  getPoints()
  {
    return d_points;
  }

private:
  void
  checkInput() const;

  virtual void
  outputHelper(ProblemSpecP& ps) const;

  void
  readPoints(const std::string& file);

  void
  readTriangles(const std::string& file);

  void
  readPLYMesh();

  void
  readOBJMesh();

  void
  readSTLMesh();

  void
  readMeshFromAsciiStlFile(std::ifstream& in);

  void
  readMeshFromBinaryStlFile(std::ifstream in);

  void
  mergeIdenticalVertices();

  void
  scaleTranslateReflect();

  void
  findBoundingBox();

  void
  makePlanes();

  void
  makeTriangleBoxes();

  void
  insideTriangle(Point& p, int i, int& NCS, int& NES) const;

  std::string d_file;
  std::string d_fileType;
  double d_scale_factor;
  Vector d_trans_vector;
  Vector d_reflect_vector;
  IntVector d_axis_sequence;

  Box d_box;
  std::vector<Point> d_points;
  std::vector<IntVector> d_triangles;
  std::vector<Plane> d_planes;
  std::vector<Box> d_boxes;

  std::unique_ptr<UniformGrid> d_grid;
};

} // End namespace Uintah

#endif //__VAANGO_CORE_GEOMPIECE_TRI_GEOMETRYPIECE_H__
