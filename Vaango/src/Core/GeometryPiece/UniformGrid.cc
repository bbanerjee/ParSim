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

#include <Core/GeometryPiece/UniformGrid.h>
#include <Core/Malloc/Allocator.h>

#include <iostream>
#include <map>

using namespace Uintah;
using TriangleList = std::list<Triangle>;

Triangle::Triangle(Point& p1, Point& p2, Point& p3)
{
  d_points[0] = p1;
  d_points[1] = p2;
  d_points[2] = p3;
  d_plane     = Plane(p1, p2, p3);
}

Point
Triangle::centroid()
{
  Vector cent(0., 0., 0);
  for (int i = 0; i < 3; i++) {
    cent += d_points[i].asVector();
  }

  cent /= 3.;

  return Point(cent.x(), cent.y(), cent.z());
}

Point
Triangle::vertex(int i)
{
  return d_points[i];
}

TriangleList
Triangle::makeTriangleList(std::vector<IntVector>& triangle_nodes,
                           std::vector<Point>& node_coords)
{
  TriangleList tri_list;
  for (auto nodes : triangle_nodes) {
    auto node_1 = nodes[0];
    auto node_2 = nodes[1];
    auto node_3 = nodes[2];
    Triangle triangle(node_coords[node_1],
                      node_coords[node_2],
                      node_coords[node_3]);
    tri_list.emplace_back(node_coords[node_1],
                          node_coords[node_2],
                          node_coords[node_3]);
  }
  /*
  std::vector<IntVector>::const_iterator tri_itr;
  for (tri_itr = tris.begin(); tri_itr != tris.end(); ++tri_itr) {
    Triangle triangle = Triangle(pts[(*tri_itr).x()],pts[(*tri_itr).y()],
                       pts[(*tri_itr).z()]);
    tri_list.push_back(triangle);
  }
  */
  return tri_list;
}

bool
Triangle::inside(Point& pt) const
{
  Vector plane_normal     = d_plane.normal();
  Vector plane_normal_abs = Abs(plane_normal);
  double largest          = plane_normal_abs.maxComponent();

  // WARNING: if dominant_coord is not 1-3, then this code breaks...
  int dominant_coord = -1;
  if (largest == plane_normal_abs.x()) {
    dominant_coord = 1;
  } else if (largest == plane_normal_abs.y()) {
    dominant_coord = 2;
  } else if (largest == plane_normal_abs.z()) {
    dominant_coord = 3;
  }

  int x = static_cast<int>(Triangle::Coord::X);
  int y = static_cast<int>(Triangle::Coord::Y);
  switch (dominant_coord) {
    case 1:
      x = static_cast<int>(Triangle::Coord::Y);
      y = static_cast<int>(Triangle::Coord::Z);
      break;
    case 2:
      x = static_cast<int>(Triangle::Coord::Z);
      y = static_cast<int>(Triangle::Coord::X);
      break;
    case 3:
      x = static_cast<int>(Triangle::Coord::X);
      y = static_cast<int>(Triangle::Coord::Y);
      break;
  }

  double tx = pt(x), ty = pt(y);

  const Point *p1 = &d_points[2], *p2 = d_points;
  int yflag0 = ((*p1)(y) >= ty);

  bool inside = false;

  for (int i = 3; i--;) {
    int yflag1 = ((*p2)(y) >= ty);
    if (yflag0 != yflag1) {
      int xflag0 = ((*p1)(x) >= tx);
      if (xflag0 == ((*p2)(x) >= tx)) {
        if (xflag0) {
          inside = !inside;
        }
      } else {
        double cmp = ((*p2)(x) - ((*p2)(y)-ty) * ((*p1)(x) - (*p2)(x)) /
                                   ((*p1)(y) - (*p2)(y)));
        //	if (cmp >= tx)
        //	if (abs(cmp - tx) >= 0.)
        if (cmp - tx >= -1.e-8) {
          inside = !inside;
        }
      }
    }
    yflag0 = yflag1;
    p1     = p2;
    p2++;
  }
  return bool(inside);
}

Plane
Triangle::plane() const
{
  return d_plane;
}

LineSeg::LineSeg(Point& p1, Point& p2)
{
  d_points[0] = p1;
  d_points[1] = p2;
}

LineSeg::LineSeg() {}

LineSeg::~LineSeg() {}

Point
LineSeg::centroid()
{
  Vector cent(0., 0., 0);
  for (int i = 0; i < 2; i++) {
    cent += d_points[i].asVector();
  }

  cent /= 3.;

  return Point(cent.x(), cent.y(), cent.z());
}

Point
LineSeg::vertex(int i)
{
  return d_points[i];
}

UniformGrid::UniformGrid(Box& bound_box)
{
  const IntVector low(0, 0, 0), hi(10, 10, 10);
  Vector diff =
    Vector(hi.x(), hi.y(), hi.z()) - Vector(low.x(), low.y(), low.z());
  d_bound_box = bound_box;
  d_max_min =
    (bound_box.upper().asVector() - bound_box.lower().asVector()) / diff;
  d_grid.resize(low, hi);
}

UniformGrid::UniformGrid(const UniformGrid& copy)
{
  d_bound_box = copy.d_bound_box;
  d_max_min   = copy.d_max_min;

  d_grid.resize(copy.d_grid.getLowIndex(), copy.d_grid.getHighIndex());

  for (auto gridIter = copy.d_grid.begin(); gridIter != copy.d_grid.end();
       gridIter++) {
    IntVector index = gridIter.getIndex();
    d_grid[index]   = copy.d_grid[index];
  }
}

UniformGrid&
UniformGrid::operator=(const UniformGrid& rhs)
{
  if (this == &rhs) {
    return *this;
  }

  std::cout << "d_grid size = " << d_grid.size() << std::endl;
  if (d_grid.size() != IntVector(0, 0, 0)) {
    // Delete the lhs stuff grid and copy the rhs to it
    for (auto gridIter = d_grid.begin(); gridIter != d_grid.end(); gridIter++) {
      IntVector index = gridIter.getIndex();
      d_grid[index].clear();
    }
  }

  d_grid.resize(rhs.d_grid.getLowIndex(), rhs.d_grid.getHighIndex());

  d_bound_box = rhs.d_bound_box;
  d_max_min   = rhs.d_max_min;

  for (auto gridIter = rhs.d_grid.begin(); gridIter != rhs.d_grid.end();
       gridIter++) {
    IntVector index = gridIter.getIndex();
    d_grid[index]   = rhs.d_grid[index];
  }

  return *this;
}

IntVector
UniformGrid::cellID(Point point)
{
  Vector pt_diff = point.asVector() - (d_bound_box.lower()).asVector();
  Vector id      = pt_diff / d_max_min;
  int i          = (int)floor(id.x());
  int j          = (int)floor(id.y());
  int k          = (int)floor(id.z());
  return IntVector(i, j, k);
}

void
UniformGrid::buildUniformGrid(TriangleList& polygons)
{
  for (auto triangle : polygons) {
    IntVector v0 = cellID(triangle.vertex(0));
    IntVector v1 = cellID(triangle.vertex(1));
    IntVector v2 = cellID(triangle.vertex(2));
#if 0
    if (v0 > d_grid.getHighIndex())
      std::cout << "v0 = " << v0 << std::endl;
    if (v1 > d_grid.getHighIndex())
      std::cout << "v1 = " << v1 << std::endl;
    if (v2 > d_grid.getHighIndex())
      std::cout << "v2 = " << v2 << std::endl;

    std::cout << "Triangle = " << triangle.vertex(0) << " " << triangle.vertex(1) << " " 
	 << triangle.vertex(2) << std::endl;
    std::cout << "v0 " << v0 << " v1 " << v1 << " v2 " << v2 << std::endl;
#endif
    IntVector low = Min(v0, v1);
    low           = Min(low, v2);
    IntVector hi  = Max(v0, v1);
    hi            = Max(hi, v2);
    for (int i = low.x(); i <= hi.x(); i++) {
      for (int j = low.y(); j <= hi.y(); j++) {
        for (int k = low.z(); k <= hi.z(); k++) {
          IntVector id(i, j, k);
          //  std::cout << "Inserting into cellID = " << id << std::endl;
          d_grid[id].push_back(triangle);
        }
      }
    }
  }
}

void
UniformGrid::countIntersections(const Point& pt, int& crossings)
{
  // Make a ray and shoot it in the +x direction
  Vector infinity      = Vector(pt.x() + 1e10, pt.y(), pt.z());
  IntVector test_pt_id = cellID(pt);
  IntVector start      = d_grid.getLowIndex();
  IntVector stop       = d_grid.getHighIndex();

  std::map<double, Triangle> cross_map;
  for (int i = start.x(); i < stop.x(); i++) {
    IntVector curr(i, test_pt_id.y(), test_pt_id.z());
    TriangleList triangles = d_grid[curr];
    for (auto triangle : triangles) {
      Point hit;
      if ((triangle.plane()).Intersect(pt, infinity, hit)) {
        Vector int_ray = hit.asVector() - pt.asVector();
        double cos_angle =
          Dot(infinity, int_ray) / (infinity.length() * int_ray.length());
        if (cos_angle < 0.) {
          continue;
        }
        if (triangle.inside(hit)) {
#if 0
	  std::cout << "Inside_new hit = " << hit << "vertices: " 
	       << triangle.vertex(0) << " " << triangle.vertex(1) <<  " "
	       << triangle.vertex(2) << std::endl;
#endif
          double distance = int_ray.length();
          auto duplicate  = cross_map.find(distance);
          if (duplicate == cross_map.end()) {
            cross_map[distance] = triangle;
            crossings++;
          }
        }
      }
    }
  }
} // count intersections infinite ray in x+ direction

void
UniformGrid::countIntersectionsx(const Point& pt, int& crossings)
{

  // Assuming that the caller intends to check the intersection between pt and
  // infinity:
  Vector infinity      = Vector(pt.x() + 1.e10, pt.y() + 10000., pt.z() + 100.);
  IntVector test_pt_id = cellID(pt);
  IntVector start      = d_grid.getLowIndex();
  IntVector stop       = d_grid.getHighIndex();

  std::map<double, Triangle> cross_map;

  for (int i = start.x(); i < stop.x(); i++) {

    IntVector curr(i, test_pt_id.y(), test_pt_id.z());
    TriangleList triangles = d_grid[curr];

    for (const auto& triangle : triangles) {
      Point hit;
      if (triangle.plane().Intersect(pt, infinity, hit)) {
        Vector int_ray = hit.asVector() - pt.asVector();
        double cos_angle =
          Dot(infinity, int_ray) / (infinity.length() * int_ray.length());

        if (cos_angle < 0.) {
          continue;
        }

        if (triangle.inside(hit)) {
          double distance = int_ray.length();
          auto duplicate  = cross_map.find(distance);

          if (duplicate == cross_map.end()) {
            cross_map[distance] = triangle;
            crossings++;
          }
        }
      }
    }
  }
}

void
UniformGrid::countIntersectionsy(const Point& pt, int& crossings)
{

  // Assuming that the caller intends to check the intersection between pt and
  // infinity:
  Vector infinity      = Vector(pt.x() + 100., pt.y() + 1.e10, pt.z() + 10000.);
  IntVector test_pt_id = cellID(pt);
  IntVector start      = d_grid.getLowIndex();
  IntVector stop       = d_grid.getHighIndex();

  std::map<double, Triangle> cross_map;

  for (int i = start.y(); i < stop.y(); i++) {

    IntVector curr(test_pt_id.x(), i, test_pt_id.z());
    TriangleList triangles = d_grid[curr];

    for (const auto& triangle : triangles) {
      Point hit;
      if (triangle.plane().Intersect(pt, infinity, hit)) {
        Vector int_ray = hit.asVector() - pt.asVector();
        double cos_angle =
          Dot(infinity, int_ray) / (infinity.length() * int_ray.length());

        if (cos_angle < 0.) {
          continue;
        }

        if (triangle.inside(hit)) {
          double distance = int_ray.length();
          auto duplicate  = cross_map.find(distance);

          if (duplicate == cross_map.end()) {
            cross_map[distance] = triangle;
            crossings++;
          }
        }
      }
    }
  }
}

void
UniformGrid::countIntersectionsz(const Point& pt, int& crossings)
{

  // Assuming that the caller intends to check the intersection between pt and
  // infinity:
  Vector infinity      = Vector(pt.x() + 10000., pt.y() + 100., pt.z() + 1.e10);
  IntVector test_pt_id = cellID(pt);
  IntVector start      = d_grid.getLowIndex();
  IntVector stop       = d_grid.getHighIndex();

  std::map<double, Triangle> cross_map;

  for (int i = start.z(); i < stop.z(); i++) {

    IntVector curr(test_pt_id.x(), test_pt_id.y(), i);
    TriangleList triangles = d_grid[curr];

    for (const auto& triangle : triangles) {
      Point hit;
      if (triangle.plane().Intersect(pt, infinity, hit)) {
        Vector int_ray = hit.asVector() - pt.asVector();
        double cos_angle =
          Dot(infinity, int_ray) / (infinity.length() * int_ray.length());

        if (cos_angle < 0.) {
          continue;
        }

        if (triangle.inside(hit)) {
          double distance = int_ray.length();
          auto duplicate  = cross_map.find(distance);

          if (duplicate == cross_map.end()) {
            cross_map[distance] = triangle;
            crossings++;
          }
        }
      }
    }
  }
}

void
UniformGrid::countIntersections(const Point& pt,
                                const Point& pt_away,
                                int& crossings,
                                double& min_distance)
{

  // This method doesn't assume anything about the direction in which the
  // ray is sent out from the original point, thus we must loop over all
  // cells in the bounding box, which makes this approach more expensive.
  // This method does return the max number of intersections.
  // We draw a ray from the origin (pt) to the point away (pt_away) and
  // then attempt to intersect with the triangulated geometry inside of this
  // cell. A crossing is only counted if the length of the intersected ray is
  // less than or equal to the length of the distance between the two points.
  Vector infinity =
    Vector(pt_away.x() - pt.x(), pt_away.y() - pt.y(), pt_away.z() - pt.z());
  IntVector test_pt_id      = cellID(pt);
  IntVector test_pt_away_id = cellID(pt_away);
  IntVector start           = d_grid.getLowIndex();
  IntVector stop            = d_grid.getHighIndex();

  IntVector min = Min(test_pt_id, test_pt_away_id);
  IntVector max = Max(test_pt_id, test_pt_away_id);

  if (min.x() > start.x() && min.x() < stop.x()) {
    start[0] = min[0];
  }
  if (min.y() > start.y() && min.y() < stop.y()) {
    start[1] = min[1];
  }
  if (min.z() > start.z() && min.z() < stop.z()) {
    start[2] = min[2];
  }

  if (max.x() < stop.x() && max.x() > start.x()) {
    stop[0] = max[0];
  }
  if (max.y() < stop.y() && max.y() > start.y()) {
    stop[1] = max[1];
  }
  if (max.z() < stop.z() && max.z() > start.z()) {
    stop[2] = max[2];
  }

  // maybe the bounding boxes don't even overlap
  if (min.x() > stop.x() || max.x() < start.x()) {
    start[0] = stop[0];
  }
  if (min.y() > stop.y() || max.y() < start.y()) {
    start[0] = stop[0];
  }
  if (min.z() > stop.z() || max.z() < start.z()) {
    start[0] = stop[0];
  }

  std::map<double, Triangle> cross_map;

  min_distance = 1.e10;

  for (int i = start.x(); i < stop.x(); i++) {
    for (int j = start.y(); j < stop.y(); j++) {
      for (int k = start.z(); k < stop.z(); k++) {

        IntVector curr(i, j, k);

        TriangleList triangles = d_grid[curr];

        for (const auto& triangle : triangles) {
          Point hit;
          if (triangle.plane().Intersect(pt, infinity, hit)) {
            Vector int_ray = hit.asVector() - pt.asVector();
            double cos_angle =
              Dot(infinity, int_ray) / (infinity.length() * int_ray.length());
            if (cos_angle < 0.) {
              continue;
            }
            if (triangle.inside(hit)) {
              double distance = int_ray.length();
              if (distance <= infinity.length()) {
                auto duplicate = cross_map.find(distance);
                if (duplicate == cross_map.end()) {
                  cross_map[distance] = triangle;
                  crossings++;
                }
                if (distance < min_distance) {
                  min_distance = distance;
                }
              }
            }
          }
        }
      }
    }
  }
}