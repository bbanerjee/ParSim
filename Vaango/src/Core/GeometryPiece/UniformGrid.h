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

#ifndef __VAANGO_CORE_GEOMPIECE_UNIFORM_GRID_H__
#define __VAANGO_CORE_GEOMPIECE_UNIFORM_GRID_H__

#include <Core/GeometryPiece/GeometryPiece.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Plane.h>
#include <Core/Geometry/Point.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Variables/Array3.h>

#include <list>
#include <vector>

namespace Uintah {

class Triangle
{
public:
  enum class Coord
  {
    X = 0,
    Y = 1,
    Z = 2
  };

  Triangle(Point& p1, Point& p2, Point& p3);
  Triangle()  = default;
  ~Triangle() = default;

  Point
  centroid();

  Point
  vertex(int i);

  using TriangleList = std::list<Triangle>;
  TriangleList
  makeTriangleList(std::vector<IntVector>& connectivity,
                   std::vector<Point>& coordinates);
  bool
  inside(Point& p) const;

  Plane
  plane() const;

private:
  Point d_points[3];
  Plane d_plane;
};

class LineSeg
{

public:
  enum coord
  {
    X = 0,
    Y = 1,
    Z = 2
  };

  LineSeg(Point& p1, Point& p2);
  LineSeg();
  ~LineSeg();

  Point
  centroid();

  Point
  vertex(int i);

private:
  Point d_points[2];
};

class UniformGrid
{
public:
  UniformGrid(Box& bound_box);
  ~UniformGrid() = default;

  UniformGrid&
  operator=(const UniformGrid&);
  UniformGrid(const UniformGrid&);

  IntVector
  cellID(Point point);

  using TriangleList = std::list<Triangle>;
  void
  buildUniformGrid(TriangleList& polygons);

  /** @brief Assume the ray goes to infinity **/
  void
  countIntersections(const Point& ray, int& crossings);

  void
  countIntersectionsx(const Point& pt, int& crossings);

  void
  countIntersectionsy(const Point& pt, int& crossings);

  void
  countIntersectionsz(const Point& pt, int& crossings);

  /** @brief Let the user specify the second point to define the ray (pt ->
     pt_away). This returns the total number of crossing and the min distance
     from pt **/
  void
  countIntersections(const Point& pt,
                     const Point& pt_away,
                     int& crossings,
                     double& min_distance);

private:
  Array3<TriangleList> d_grid;
  Box d_bound_box;
  Vector d_max_min;
};

} // End namespace Uintah

#endif // __VAANGO_CORE_GEOMPIECE_UNIFORM_GRID_H__
