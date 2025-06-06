#ifndef CORE_GEOMETRY_FACE_EDGE_H
#define CORE_GEOMETRY_FACE_EDGE_H

#include <Core/Math/Vec.h>
#include <iostream>
#include <vector>

namespace dem {

using Point = dem::Vec;
using Vertex = dem::Vec;

struct Edge
{
  std::array<Vertex, 2> vertex;

  Edge(const Vertex& v1, const Vertex& v2)
  {
    vertex[0] = v1;
    vertex[1] = v2;
  }

  // Distance from edge to point
  // Also returns parameter (t) locating the point of projection (p) on the edge
  // such that p = (1 - t) v1 + t v2
  inline
  std::pair<REAL, REAL> distance(const Point& point) const
  {
    Vertex v1v0 = vertex[0] - point;
    Vertex v2v0 = vertex[1] - point;
    Vertex v2v1 = vertex[1] - vertex[0];
    REAL lengthSq = v2v1.lengthSq();

    REAL t = -dot(v1v0, v2v1)/lengthSq;
    REAL dist = cross(-v1v0, -v2v0).length()/std::sqrt(lengthSq);
    return std::make_pair(dist, t);
  }

  friend std::ostream& operator<<(std::ostream& os, const Edge& edge)
  {
    os << "v1: " << edge.vertex[0] << ": v2: " << edge.vertex[1];
    return os;
  }
};

struct Triangle
{
  std::array<Vertex, 3> vertex;

  Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
  {
    vertex[0] = v1;
    vertex[1] = v2;
    vertex[2] = v3;
  }

  inline
  const Vec& v1() const { return vertex[0]; }
  inline
  const Vec& v2() const { return vertex[1]; }
  inline
  const Vec& v3() const { return vertex[2]; }

  std::pair<REAL, REAL> barycentricCoords(const Point& point)
  {
    Vec vp = point - vertex[0];
    Vec v1 = vertex[1] - vertex[0];
    Vec v2 = vertex[2] - vertex[0];
    auto vp1 = dot(vp, v1);
    auto vp2 = dot(vp, v2);
    auto v11 = dot(v1, v1);
    auto v12 = dot(v1, v2);
    auto v22 = dot(v2, v2);
    auto denom = v11 * v22 - v12 * v12;
    auto v = (v11 * vp2 - v12 * vp1) / denom;
    auto u = (vp1 - v * v12) / v11;
    return std::make_pair(u, v);
  }
};

struct Face
{
  enum class Location
  {
    NONE = 0,
    VERTEX = 1,
    EDGE = 2,
    FACE = 3
  };

  std::array<Vertex, 4> vertex;

  Face(const Vertex& v1, const Vertex& v2, const Vertex& v3, const Vertex& v4)
  {
    vertex[0] = v1;
    vertex[1] = v2;
    vertex[2] = v3;
    vertex[3] = v4;
  }

  inline
  const Vec& v1() const { return vertex[0]; }
  inline
  const Vec& v2() const { return vertex[1]; }
  inline
  const Vec& v3() const { return vertex[2]; }
  inline
  const Vec& v4() const { return vertex[3]; }

  inline
  bool isValid()
  {
    REAL area_diag = area();
    if (std::abs(area_diag) < 1.0e-16) return false;

    REAL area_tri1 = 0.5*dot(normal(), 
      cross(vertex[1] - vertex[0], vertex[2] - vertex[0]));
    REAL area_tri2 = 0.5*dot(normal(), 
      cross(vertex[2] - vertex[0], vertex[3] - vertex[0]));
    if (std::abs(area_diag - area_tri1 - area_tri2) > 1.0e-16) return false;
    return true;
  }

  inline
  REAL area()
  {
    return 0.5*dot(normal(), cross(vertex[2] - vertex[0], vertex[3] - vertex[1]));
  }

  inline
  Vec normal() const 
  {
    Vec normal = cross(vertex[1] - vertex[0], vertex[2] - vertex[0]);
    normal.normalizeInPlace();
    return normal;
  }

  // Distance from face to point
  // Also returns parameter (t_j) locating the point of projection (p_j) 
  // on the j-th edge such that p_j = (1 - t_j) v1_j + t_j v2_j
  std::pair<REAL, std::vector<std::pair<REAL, REAL>>> 
    distance(const Point& point) const
  {
    Vec normal = this->normal();
    REAL dist = dot(normal, (point - vertex[0]));
    Point pointOnFace = point - normal*dist;
    //std::cout << "normal = " << normal << " point = " << point 
    //          << " v1 = " << vertex[0] << " point on face = " << pointOnFace
    //          << " dist = " << dist << "\n";

    Edge edge1(vertex[0], vertex[1]);
    std::pair<REAL, REAL> t1d1  = edge1.distance(pointOnFace);

    Edge edge2(vertex[1], vertex[2]);
    std::pair<REAL, REAL> t2d2  = edge2.distance(pointOnFace);

    Edge edge3(vertex[2], vertex[3]);
    std::pair<REAL, REAL> t3d3  = edge3.distance(pointOnFace);

    Edge edge4(vertex[3], vertex[0]);
    std::pair<REAL, REAL> t4d4  = edge4.distance(pointOnFace);

    return std::make_pair(dist,
      std::vector<std::pair<REAL, REAL>>({t1d1, t2d2, t3d3, t4d4}));
  }

  // Distance from face to point
  // Also returns parameters (u1, v1, u2, v2) locating the point of projection (p) 
  // in barycentric coordinates of the two triangles forming the face 
  std::pair<REAL, std::vector<std::pair<REAL, REAL>>> 
    distanceAndBarycentricCoords(const Point& point) const
  {
    Vec normal = this->normal();
    REAL dist = dot(normal, (point - vertex[0]));
    Point pointOnFace = point - normal*dist;
    //std::cout << "normal = " << normal << " point = " << point 
    //          << " v1 = " << vertex[0] << " point on face = " << pointOnFace
    //          << " dist = " << dist << "\n";

    Triangle triangle1(vertex[0], vertex[1], vertex[2]);
    std::pair<REAL, REAL> u1v1  = triangle1.barycentricCoords(pointOnFace);

    Triangle triangle2(vertex[0], vertex[3], vertex[2]);
    std::pair<REAL, REAL> u2v2  = triangle2.barycentricCoords(pointOnFace);

    return std::make_pair(dist, 
      std::vector<std::pair<REAL, REAL>>({u1v1, u2v2}));
  }

  friend std::ostream& operator<<(std::ostream& os, const Face& face)
  {
    os << "v1: " << face.vertex[0] << ": v2: " << face.vertex[1]
       << ": v3: " << face.vertex[2] << ": v4: " << face.vertex[3];
    return os;
  }
};
} // end namespace dem

#endif
