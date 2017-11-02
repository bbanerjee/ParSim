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
  Vertex v1;
  Vertex v2;

  Edge(Vertex v1, Vertex v2)
    : v1(std::move(v1))
    , v2(std::move(v2))
  {
  }

  // Distance from edge to point
  // Also returns parameter (t) locating the point of projection (p) on the edge
  // such that p = (1 - t) v1 + t v2
  std::pair<REAL, REAL> distance(const Point& point)
  {
    Vertex v1v0 = v1 - point;
    Vertex v2v0 = v2 - point;
    Vertex v2v1 = v2 - v1;
    REAL lengthSq = v2v1.lengthSq();

    REAL t = -dot(v1v0, v2v1)/lengthSq;
    REAL dist = cross(-v1v0, -v2v0).length()/std::sqrt(lengthSq);
    return std::make_pair(t, dist);
  }

  friend std::ostream& operator<<(std::ostream& os, const Edge& edge)
  {
    os << "v1: " << edge.v1 << ": v2: " << edge.v2;
    return os;
  }
};

struct Face
{
  Vertex v1;
  Vertex v2;
  Vertex v3;
  Vertex v4;

  Face(Vertex v1, Vertex v2, Vertex v3, Vertex v4)
    : v1(std::move(v1))
    , v2(std::move(v2))
    , v3(std::move(v3))
    , v4(std::move(v4))
  {
  }

  Vec normal() const 
  {
    Vec normal = cross(v2 - v1, v3 - v1);
    normal.normalizeInPlace();
    return normal;
  }

  // Distance from face to point
  // Also returns parameter (t_j) locating the point of projection (p_j) 
  // on the j-th edge such that p_j = (1 - t_j) v1_j + t_j v2_j
  std::pair<std::vector<std::pair<REAL, REAL>>, REAL> 
    distance(const Point& point) const
  {
    Vec normal = this->normal();
    REAL dist = dot(normal, (point - v1));
    Point pointOnFace = point - normal*dist;

    Edge edge1(v1, v2);
    std::pair<REAL, REAL> t1d1  = edge1.distance(pointOnFace);

    Edge edge2(v2, v3);
    std::pair<REAL, REAL> t2d2  = edge2.distance(pointOnFace);

    Edge edge3(v3, v4);
    std::pair<REAL, REAL> t3d3  = edge3.distance(pointOnFace);

    Edge edge4(v4, v1);
    std::pair<REAL, REAL> t4d4  = edge4.distance(pointOnFace);

    return std::make_pair(
      std::vector<std::pair<REAL, REAL>>({t1d1, t2d2, t3d3, t4d4}), dist);
  }

  friend std::ostream& operator<<(std::ostream& os, const Face& face)
  {
    os << "v1: " << face.v1 << ": v2: " << face.v2
       << ": v3: " << face.v3 << ": v4: " << face.v4;
    return os;
  }
};
} // end namespace dem

#endif
