#include <Core/Const/Constants.h>
#include <Core/Geometry/Ellipsoid.h>
#include <Core/Geometry/OrientedBox.h>
#include <Core/Geometry/FaceEdge.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Vec.h>
#include <iostream>

namespace dem {

void 
Ellipsoid::normalize_axes() {
  d_axes[0].normalizeInPlace();
  d_axes[1].normalizeInPlace();
  d_axes[2].normalizeInPlace();
}

OrientedBox
Ellipsoid::getOrientedBoundingBox() const
{
  OrientedBox box(d_center, d_axes[0], d_axes[1], d_axes[2],
                  d_radii[0], d_radii[1], d_radii[2]);
  return box;
}

Matrix3 
Ellipsoid::toUnitSphereTransformationMatrix() const
{
  Matrix3 mat;
  for (auto ii = 0; ii < 3; ++ii) {
    for (auto jj = 0; jj < 3; ++jj) {
      mat(ii, jj) = d_axes[ii][jj]/d_radii[ii];
    }
  }
  return mat;
}

bool 
Ellipsoid::containsPoint(const Vec& point) const
{
  Matrix3 N = toUnitSphereTransformationMatrix();
  Vec pp = N*(point - d_center) ;
  if (pp.length() > 1.0) {
    return false;
  }
  return true;
}

bool
Ellipsoid::intersects(const OrientedBox& box) const
{
  OrientedBox boundingBox = getOrientedBoundingBox();
  if (!boundingBox.intersects(box)) {
    return false;
  }

  // Identify the visible faces
  Vec ww = d_center - box.center();
  std::vector<std::pair<int, float>> visibleFaces;
  for (int ii = 0; ii < 3; ++ii) {
    REAL dd = dot(ww, box.axis(ii));
    if (std::abs(dd) > box.extent(ii)) {
      visibleFaces.push_back(std::make_pair(ii, std::copysign(1, dd)));
    }
  }
  // No face chosen; center of ellipsoid is inside box
  if (visibleFaces.size() < 1) {
    //std::cout << "Visible faces = " << visibleFaces.size() << "\n";
    return true;
  }
                            
  // Compute transformation matrix
  Matrix3 N = toUnitSphereTransformationMatrix();

  // Transform the oriented box into a parallelepiped
  Vec center = N * (-ww);
  std::array<Vec, 3> axes;
  for (auto ii = 0u; ii < 3; ii++) {
    axes[ii] = N * (box.axis(ii) * box.extent(ii));
  }

  // Compute the vertex locations
  std::vector<Vec> vertices; vertices.reserve(8);
  vertices.push_back(center - axes[0] - axes[1] - axes[2]);
  vertices.push_back(center + axes[0] - axes[1] - axes[2]);
  vertices.push_back(center + axes[0] + axes[1] - axes[2]);
  vertices.push_back(center - axes[0] + axes[1] - axes[2]);
  vertices.push_back(center - axes[0] - axes[1] + axes[2]);
  vertices.push_back(center + axes[0] - axes[1] + axes[2]);
  vertices.push_back(center + axes[0] + axes[1] + axes[2]);
  vertices.push_back(center - axes[0] + axes[1] + axes[2]);

  //std::copy(vertices.begin(), vertices.end(),
  //          std::ostream_iterator<Vec>(std::cout, "\n"));

  // Check if the vertices are inside the transformed ellipsoid
  // Note that the transformed ellipsoid is now a sphere centered at (0,0,0)
  for (const auto& vertex : vertices) {
    if (vertex.lengthSq() < 1.0) {
      //std::cout << "Vertex distance = " << vertex.lengthSq() << "\n";
      return true;
    }
  }

  // Set up the edges
  constexpr std::array<std::array<int, 2>, 12> edges = {{
    {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}},
    {{4, 5}}, {{5, 6}}, {{6, 7}}, {{7, 4}},
    {{0, 4}}, {{1, 5}}, {{2, 6}}, {{3, 7}},
    }};
  
  // Check sphere edge intersections
  Point sphereCenter(0, 0, 0);
  for (const auto& edge : edges) {
    int v1 = edge[0];
    int v2 = edge[1];
    Edge ee(vertices[v1], vertices[v2]);
    auto td = ee.distance(sphereCenter);
    if (td.second < 1.0) {
      //std::cout << "Edge distance = " << td.second << "\n";
      if (!(td.first < 0 || td.first > 1)) {
        return true;
      }
    }
  }

  // Set up the faces
  constexpr std::array<std::array<int, 4>, 6> faces = {{
     {{0, 4, 7, 3}}, // x-
     {{1, 2, 6, 5}}, // x+
     {{0, 1, 5, 4}}, // y-
     {{2, 3, 7, 6}}, // y+
     {{0, 3, 2, 1}}, // z-
     {{4, 5, 6, 7}}  // z+
    }};

  for (const auto& face : faces) {
    int v0 = face[0];
    int v1 = face[1];
    int v2 = face[2];
    int v3 = face[3];
    Face ff(vertices[v0], vertices[v1], vertices[v2], vertices[v3]);
    auto td = ff.distance(sphereCenter);
    if (std::abs(td.second) < 1.0) {

      //std::cout << "Face distance = " << td.second << "\n";
      //std::cout << "td = ";
      //for (auto val : td.first) {
      //  std::cout << val.first << " ";
      //}
      //std::cout << "\n";

      if (!(td.first[0].first < 0 || td.first[0].first > 1 ||
            td.first[1].first < 0 || td.first[1].first > 1 ||
            td.first[2].first < 0 || td.first[2].first > 1 ||
            td.first[3].first < 0 || td.first[3].first > 1)) {
        return true;
      }
    }
  }

  return false;
}

std::pair<bool, std::pair<Face::Location, int> >
Ellipsoid::intersects(const Face& face) const
{
  // Compute transformation matrix (ellipsoid to sphere)
  Matrix3 N = toUnitSphereTransformationMatrix();

  // Create transformed face
  Face transformedFace(N * (face.v1() - d_center), 
                       N * (face.v2() - d_center), 
                       N * (face.v3() - d_center), 
                       N * (face.v4() - d_center));

  // Compute distance to face from (0, 0, 0)
  auto face_td = transformedFace.distance(Point(0, 0, 0));

  // If the distance is greater than 1 there is no intersection
  //std::cout << "Face = " << face << "\n";
  //std::cout << "X-Face = " << transformedFace << "\n";
  //std::cout << "Distance from face = " << face_td.second << "\n";
  if (std::abs(face_td.second) > 1.0) {
    return std::make_pair(false, std::make_pair(Face::Location::NONE, 0));
  }

  // Check sphere vertex intersections
  std::array<Vec, 4> vertices = {{
    transformedFace.v1(), transformedFace.v2(), transformedFace.v3(), 
    transformedFace.v4() 
  }};
  int index = 0;
  for (const auto& vertex : vertices) {
    //std::cout << "Distance from vertex: " << vertex << " = " << vertex.lengthSq() << "\n";
    if (vertex.lengthSq() < 1.0) {
      return std::make_pair(true, std::make_pair(Face::Location::VERTEX, index));
    }
    ++index;
  }

  // Check sphere edge intersections
  std::array<Edge, 4> edges = {{
    Edge(transformedFace.v1(), transformedFace.v2()),
    Edge(transformedFace.v2(), transformedFace.v3()),
    Edge(transformedFace.v3(), transformedFace.v4()),
    Edge(transformedFace.v4(), transformedFace.v1())
    }};
  index = 0;
  for (const auto& edge : edges) {
    auto edge_td = edge.distance(Point(0, 0, 0));
    //std::cout << "Distance from Edge = " << edge_td.second << "\n";
    if (edge_td.second < 1.0) {
      if (!(edge_td.first < 0 || edge_td.first > 1)) {
        return std::make_pair(true, std::make_pair(Face::Location::EDGE, index));
      }
    }
    ++index;
  }

  // If none of the edges are intersected
  //std::cout << "Edge intersect locations = " 
  //          << face_td.first[0].first << " "
  //          << face_td.first[1].first << " "
  //          << face_td.first[2].first << " "
  //          << face_td.first[3].first << "\n";
  if (!(face_td.first[0].first < 0 || face_td.first[0].first > 1 ||
        face_td.first[1].first < 0 || face_td.first[1].first > 1 ||
        face_td.first[2].first < 0 || face_td.first[2].first > 1 ||
        face_td.first[3].first < 0 || face_td.first[3].first > 1)) {
    return std::make_pair(true, std::make_pair(Face::Location::FACE, 0));
  }

  return std::make_pair(false, std::make_pair(Face::Location::NONE, 0));
}

void 
Ellipsoid::rotate(REAL angle, const Vec& axis)
{
  Matrix3 rot(angle, axis);
  d_axes[0] = rot*d_axes[0];
  d_axes[1] = rot*d_axes[1];
  d_axes[2] = rot*d_axes[2];
}

void 
Ellipsoid::translate(const Vec& dist)
{
  d_center += dist;
}

std::ostream&
operator<<(std::ostream& os, const Ellipsoid& b)
{
  os << " center = " << b.d_center 
     << " axis_a = " << b.d_axes[0] 
     << " axis_b = " << b.d_axes[1] 
     << " axis_c = " << b.d_axes[2] 
     << " radius_a = " << b.d_radii[0] 
     << " radius_b = " << b.d_radii[1] 
     << " radius_c = " << b.d_radii[2] 
     << "\n";
  return os;
}

} // namespace dem
