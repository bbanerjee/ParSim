/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Parresia research Limited, New Zealand
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

#include <Core/GeometryPiece/TriGeometryPiece.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Plane.h>
#include <Core/GeometryPiece/happly.h>
#include <Core/Grid/Box.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>

namespace Uintah {

using TriangleList = std::list<Triangle>;

TriGeometryPiece::TriGeometryPiece(ProblemSpecP& ps)
{
  d_name = "Unnamed Tri";

  ps->require("file_name_prefix", d_file);
  ps->getWithDefault("file_type", d_fileType, "default");

  ps->getWithDefault("scaling_factor", d_scale_factor, 1.0);
  Vector trans(0.0, 0.0, 0.0);
  ps->getWithDefault("translation_vector", d_trans_vector, trans);
  Vector reflect(1.0, 1.0, 1.0);
  ps->getWithDefault("reflection_vector", d_reflect_vector, reflect);
  IntVector axisSeq(1, 2, 3);
  ps->getWithDefault("axis_sequence", d_axis_sequence, axisSeq);

  std::transform(d_fileType.begin(),
                 d_fileType.end(),
                 d_fileType.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  checkInput();

  if (d_fileType == "default") {
    readPoints(d_file);
    readTriangles(d_file);
  } else if (d_fileType == "ply") {
    readPLYMesh();
  } else if (d_fileType == "obj") {
    readOBJMesh();
  } else {
    readSTLMesh();
  }

  scaleTranslateReflect();
  findBoundingBox();

  makePlanes();
  makeTriangleBoxes();

  // std::cout << "Triangulated surfaces read: \t" <<d_triangles.size() <<endl;
  Triangle tri;
  TriangleList tri_list = tri.makeTriangleList(d_triangles, d_points);
  d_grid                = std::make_unique<UniformGrid>(d_box);
  d_grid->buildUniformGrid(tri_list);
}

void
TriGeometryPiece::checkInput() const
{
  if (d_fileType == "default") {
    std::string f = d_file + ".pts";
    std::ifstream source(f.c_str());
    if (!source) {
      std::ostringstream err;
      err << "\n ERROR: Input file: opening geometry points file (" << f
          << ").\n  The file must be in the same directory the input ups "
             "file.\n"
          << "  Do not enclose the filename in quotation marks\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    source.close();
    f = d_file + ".tri";
    std::ifstream source1(f.c_str());
    if (!source) {
      std::ostringstream err;
      err << "\n ERROR: Input file: opening geometry triangle file (" << f
          << ").\n  The file must be in the same directory as the input ups "
             "file.\n"
          << "  Do not enclose the filename in quotation marks\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    source1.close();
  } else if (d_fileType == "ply") {
    std::string f1 = d_file + ".ply";
    std::string f2 = d_file + ".PLY";
    std::ifstream source1(f1.c_str());
    std::ifstream source2(f2.c_str());
    if (!source1 && !source2) {
      std::ostringstream err;
      err << "\n ERROR: Input file: opening PLY triangulated geometry file ("
          << f1
          << ").\n  The file must be in the same directory the input ups "
             "file.\n"
          << "  Do not enclose the filename in quotation marks\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    source1.close();
    source2.close();
  } else if (d_fileType == "obj") {
    std::string f1 = d_file + ".obj";
    std::string f2 = d_file + ".OBJ";
    std::ifstream source1(f1.c_str());
    std::ifstream source2(f2.c_str());
    if (!source1 && !source2) {
      std::ostringstream err;
      err << "\n ERROR: Input file: opening OBJ triangulated geometry file ("
          << f1
          << ").\n  The file must be in the same directory the input ups "
             "file.\n"
          << "  Do not enclose the filename in quotation marks\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    source1.close();
    source2.close();
  } else if (d_fileType == "stl") {
    std::string f1 = d_file + ".stl";
    std::string f2 = d_file + ".STL";
    std::ifstream source1(f1.c_str());
    std::ifstream source2(f2.c_str());
    if (!source1 && !source2) {
      std::ostringstream err;
      err << "\n ERROR: Input file: opening STL triangulated geometry file ("
          << f1
          << ").\n  The file must be in the same directory the input ups "
             "file.\n"
          << "  Do not enclose the filename in quotation marks\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    source1.close();
    source2.close();
  } else {
    std::ostringstream err;
    err << "\n ERROR: Input file: opening triangulated geometry file of "
           "unknown type ("
        << d_fileType << ").\n  Please select one of 'ply', 'obj', or 'stl'.\n"
        << "  Do not enclose the filename in quotation marks.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  if (d_scale_factor < 0.0) {
    std::ostringstream err;
    err << "\n ERROR: Input file: The geometry scaling factor cannot be less "
           "than zero.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  if (std::abs(d_reflect_vector[0]) != 1 ||
      std::abs(d_reflect_vector[1]) != 1 ||
      std::abs(d_reflect_vector[2]) != 1) {
    std::ostringstream err;
    err << "\n ERROR: Input file: The reflection vector requires the form [1, "
           "-1, 1].\n"
        << " Values other than +1 or -1 are not allowed.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  if (d_axis_sequence != IntVector(1, 2, 3) &&
      d_axis_sequence != IntVector(2, 3, 1) &&
      d_axis_sequence != IntVector(3, 1, 2) &&
      d_axis_sequence != IntVector(1, 3, 2) &&
      d_axis_sequence != IntVector(3, 2, 1) &&
      d_axis_sequence != IntVector(2, 1, 3)) {
    std::ostringstream err;
    err << "\n ERROR: Input file: The axis sequence must be [1,2,3], [2,3,1], "
           "[3,1,2] or\n"
        << " [1,3,2], [3,2,1], [2,1,3].\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

TriGeometryPiece::TriGeometryPiece(const TriGeometryPiece& copy)
{
  d_box       = copy.d_box;
  d_points    = copy.d_points;
  d_triangles = copy.d_triangles;
  d_planes    = copy.d_planes;
  d_boxes     = copy.d_boxes;
  d_grid      = std::make_unique<UniformGrid>(*copy.d_grid);
}

TriGeometryPiece&
TriGeometryPiece::operator=(const TriGeometryPiece& rhs)
{
  if (this == &rhs) {
    return *this;
  }

  // Clean out lhs
  d_points.clear();
  d_triangles.clear();
  d_planes.clear();
  d_boxes.clear();

  // Copy the rhs stuff
  d_box       = rhs.d_box;
  d_points    = rhs.d_points;
  d_triangles = rhs.d_triangles;
  d_planes    = rhs.d_planes;
  d_boxes     = rhs.d_boxes;
  d_grid      = std::make_unique<UniformGrid>(*rhs.d_grid);
  return *this;
}

void
TriGeometryPiece::outputHelper(ProblemSpecP& ps) const
{
  ps->appendElement("file_name_prefix", d_file);
  ps->appendElement("file_type", d_fileType);
  ps->appendElement("scaling_factor", d_scale_factor);
  ps->appendElement("translation_vector", d_trans_vector);
  ps->appendElement("reflection_vector", d_reflect_vector);
  ps->appendElement("axis_sequence", d_axis_sequence);
}

GeometryPieceP
TriGeometryPiece::clone() const
{
  return std::make_shared<TriGeometryPiece>(*this);
}

bool
TriGeometryPiece::inside(const Point& p) const
{
  // use crossings in all three directions
  return inside(p, false);
}

bool
TriGeometryPiece::inside(const Point& p, bool use_x_crossing) const
{
  // Count the number of times a ray from the point p
  // intersects the triangular surface.  If the number
  // of crossings is odd, the point is inside, else it
  // is outside.

  // Check if Point p is outside the bounding box
  if (!(p == Max(p, d_box.lower()) && p == Min(p, d_box.upper()))) {
    return false;
  }

  int cross = 0;
  bool is_inside;

  if (use_x_crossing) {
    is_inside = inside(p, cross);
  } else {
    is_inside = inside(p, cross, true);
  }

  return is_inside;
}

bool
TriGeometryPiece::inside(const Point& p, int& cross) const
{
  // Count the number of times a ray from the point p
  // intersects the triangular surface.  If the number
  // of crossings is odd, the point is inside, else it
  // is outside.

  // This version only tests by casting a ray in the x-direction

  // Check if Point p is outside the bounding box
  if (!(p == Max(p, d_box.lower()) && p == Min(p, d_box.upper()))) {
    return false;
  }

  d_grid->countIntersections(p, cross);
  //  cout << "Point " << p << " has " << cross << " crossings " << std::endl;
  if (cross % 2) {
    return true;
  } else {
    return false;
  }
}

bool
TriGeometryPiece::inside(const Point& p,
                         [[maybe_unused]] int& cross,
                         [[maybe_unused]] bool all_directions) const
{
  // Count the number of times a ray from the point p
  // intersects the triangular surface.  If the number
  // of crossings is odd, the point is inside, else it
  // is outside.

  // This version checks using three roughly orthogonal rays  that are
  // nearly, but not exactly, cast in the 3 ordinal directions.
  // It returns the answer that two or more of those tests agree upon

  // Check if Point p is outside the bounding box
  if (!(p == Max(p, d_box.lower()) && p == Min(p, d_box.upper()))) {
    return false;
  }

  int crossx = 0;
  int crossy = 0;
  int crossz = 0;

  d_grid->countIntersectionsx(p, crossx);
  d_grid->countIntersectionsy(p, crossy);
  d_grid->countIntersectionsz(p, crossz);

  //  cout << "Point " << p << " has " << cross << " crossings " << std::endl;
  if ((crossx % 2 == 1 && crossy % 2 == 1) ||
      (crossx % 2 == 1 && crossz % 2 == 1) ||
      (crossy % 2 == 1 && crossz % 2 == 1)) {
    return true;
  } else {
    return false;
  }
}

Box
TriGeometryPiece::getBoundingBox() const
{
  return d_box;
}

void
TriGeometryPiece::readPoints(const string& file)
{
  std::string f = file + ".pts";
  std::ifstream source(f.c_str());
  if (!source) {
    std::ostringstream warn;
    warn << "\n ERROR: opening geometry pts points file (" << f
         << ").\n  The file must be in the same directory as sus \n"
         << "  Do not enclose the filename in quotation marks\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // Get the axis sequence, e.g. 1, 2, 3 == x, y, z
  //                             1, 3, 2 == x, z, y
  int first        = d_axis_sequence.x() - 1;
  int second       = d_axis_sequence.y() - 1;
  int third        = d_axis_sequence.z() - 1;
  double coords[3] = { 0.0, 0.0, 0.0 };
  double x, y, z;
  while (source >> x >> y >> z) {
    coords[first]  = x;
    coords[second] = y;
    coords[third]  = z;
    // std::cout << " x = " << x << " y = " << y << " z = " << z
    //           << "[" << first << "," << second << "," << third << "]"
    //           << "[" << coords[0] << "," << coords[1] << "," << coords[2]
    //           << "]" << "\n";

    d_points.push_back(Point(coords));
  }

  source.close();
  std::cout << "Read " << d_points.size() << " points from geometry file"
            << "\n";
}

void
TriGeometryPiece::readTriangles(const string& file)
{
  std::string f = file + ".tri";
  std::ifstream source(f.c_str());
  if (!source) {
    std::ostringstream err;
    err << "\n ERROR: opening geometry tri points file (" << f
        << ").\n   The file must be in the same directory as sus"
        << "   Do not enclose the filename in quotation marks\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  int x, y, z;
  while (source >> x >> y >> z) {
    d_triangles.push_back(IntVector(x, y, z));
  }
  source.close();

  std::cout << "Read " << d_triangles.size() << " triangles from geometry file"
            << "\n";
}

void
TriGeometryPiece::readPLYMesh()
{
  // Allow both lowercase and uppercase extensions
  std::string inFile = d_file + ".ply";
  std::ifstream source(inFile.c_str());
  if (!source) {
    inFile = d_file + ".PLY";
  }

  // Read the file with happly
  happly::PLYData plyIn(inFile);
  std::vector<std::array<double, 3>> points = plyIn.getVertexPositions();
  std::vector<std::vector<size_t>> faces    = plyIn.getFaceIndices<size_t>();

  // Save the points
  // Use the axis sequence, e.g. 1, 2, 3 == x, y, z
  //                             1, 3, 2 == x, z, y
  int first        = d_axis_sequence.x() - 1;
  int second       = d_axis_sequence.y() - 1;
  int third        = d_axis_sequence.z() - 1;
  double coords[3] = { 0.0, 0.0, 0.0 };
  for (auto point : points) {
    coords[first]  = point[0];
    coords[second] = point[1];
    coords[third]  = point[2];
    d_points.push_back(Point(coords));
  }
  std::cout << "Read " << d_points.size() << " points from PLY mesh file."
            << "\n";

  // Save the triangles
  for (auto face : faces) {
    d_triangles.push_back(IntVector(face[0], face[1], face[2]));
  }
  std::cout << "Read " << d_triangles.size() << " triangles from PLY mesh file"
            << "\n";
}

void
TriGeometryPiece::readOBJMesh()
{
  // Allow both lowercase and uppercase extensions
  std::string inFile = d_file + ".obj";
  std::ifstream source;
  source.open(inFile.c_str(), std::ifstream::in);
  if (!source.good()) {
    source.close();
    inFile = d_file + ".OBJ";
  }
  source.open(inFile.c_str(), std::ifstream::in);

  // Read the file
  std::string line;
  while (std::getline(source, line)) {
    // erase white spaces from the beginning of line
    line.erase(line.begin(),
               std::find_if(line.begin(),
                            line.end(),
                            std::not1(std::ptr_fun<int, int>(std::isspace))));

    // Ignore empty lines
    if (line.empty()) {
      continue;
    }

    // Skip comment lines
    if (line[0] == '#') {
      continue;
    }

    // Ignore lines containing texture information
    if (line[0] == 'v' && line[1] == 't') {
      continue;
    }

    // Read the vertices
    if (line[0] == 'v') {
      int first        = d_axis_sequence.x() - 1;
      int second       = d_axis_sequence.y() - 1;
      int third        = d_axis_sequence.z() - 1;
      double coords[3] = { 0.0, 0.0, 0.0 };

      std::istringstream data_stream(line);
      char data_id;
      double xcoord = 0.0, ycoord = 0.0, zcoord = 0.0;
      if (!(data_stream >> data_id >> xcoord >> ycoord >> zcoord)) {
        std::ostringstream err;
        err << "**ERROR**Could not read node input data stream from OBJ "
               "file.\n";
        throw ProblemSetupException(err.str(), __FILE__, __LINE__);
      }
      coords[first]  = xcoord;
      coords[second] = ycoord;
      coords[third]  = zcoord;
      d_points.push_back(Point(coords));
    }

    // Read the face connectivity
    if (line[0] == 'f') {
      std::istringstream data_stream4(line);
      std::istringstream data_stream3(line);
      std::istringstream data_stream4_notex(line);
      std::istringstream data_stream3_notex(line);
      char data_id;
      char dummy;
      int node1, node2, node3, node4;
      int tex1, tex2, tex3, tex4;
      std::vector<int> elem;
      if ((data_stream4 >> data_id >> node1 >> dummy >> tex1 >> node2 >>
           dummy >> tex2 >> node3 >> dummy >> tex3 >> node4 >> dummy >> tex4)) {
        // std::cout << "[" << node1 << "," << node2 << "," << node3 << "," <<
        // node4 << "]" << "\n";
        d_triangles.push_back(IntVector(node1, node2, node3));
        d_triangles.push_back(IntVector(node1, node3, node4));
        continue;
      }
      if ((data_stream3 >> data_id >> node1 >> dummy >> tex1 >> node2 >>
           dummy >> tex2 >> node3 >> dummy >> tex3)) {
        // std::cout << "[" << node1 << "," << node2 << "," << node3 << "]" <<
        // "\n";
        d_triangles.push_back(IntVector(node1, node2, node3));
        continue;
      }
      if ((data_stream4_notex >> data_id >> node1 >> node2 >> node3 >> node4)) {
        // std::cout << "[" << node1 << "," << node2 << "," << node3 << "," <<
        // node4 << "]" << "\n";
        d_triangles.push_back(IntVector(node1, node2, node3));
        d_triangles.push_back(IntVector(node1, node3, node4));
        continue;
      }
      if ((data_stream3_notex >> data_id >> node1 >> node2 >> node3)) {
        // std::cout << "[" << node1 << "," << node2 << "," << node3 << "]" <<
        // "\n";
        d_triangles.push_back(IntVector(node1, node2, node3));
        continue;
      } else {
        std::ostringstream err;
        err << "**ERROR**Could not read triangle input data stream from OBJ "
               "file.\n";
        throw ProblemSetupException(err.str(), __FILE__, __LINE__);
      }
    }
  }
  source.close();

  if (d_points.size() == 0 || d_triangles.size() == 0) {
    std::ostringstream err;
    err << "**ERROR**The OBJ file contains either no vertices or no triangles "
           "or both.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

/**
 * From:
 * https://github.com/nmwsharp/geometry-central/blob/master/src/surface/polygon_soup_mesh.cpp
 * May 23, 2020
 */
void
TriGeometryPiece::readSTLMesh()
{
  // Allow both lowercase and uppercase extensions
  std::string inFile = d_file + ".stl";
  std::ifstream source;
  source.open(inFile.c_str(), std::ifstream::in);
  if (!source.good()) {
    source.close();
    inFile = d_file + ".STL";
  }
  source.open(inFile.c_str(), std::ifstream::in);

  // parse stl format
  std::string line;
  std::getline(source, line);
  std::stringstream ss(line);
  std::string token;
  ss >> token;
  if (token == "solid") {
    readMeshFromAsciiStlFile(source);
  } else {
    readMeshFromBinaryStlFile(
      std::ifstream(inFile, std::ios::in | std::ios::binary));
  }
  mergeIdenticalVertices();
}

/**
 * From:
 * https://github.com/nmwsharp/geometry-central/blob/master/src/surface/polygon_soup_mesh.cpp
 * May 23, 2020
 */
// Assumes that first line has already been consumed
void
TriGeometryPiece::readMeshFromAsciiStlFile(std::ifstream& in)
{
  std::string line;
  std::stringstream ss;
  size_t lineNum = 1;

  auto assertToken = [&](const std::string& expected) {
    std::string token;
    ss >> token;
    if (token != expected) {
      std::ostringstream err;
      err << "**ERROR** Failed to parse ASCII stl file."
          << "\n"
          << "Error on line " << lineNum << ". Expected \"" << expected
          << "\" but token \"" << token << "\""
          << "\n"
          << "Full line: \"" << line << "\""
          << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  };

  auto nextLine = [&]() {
    if (!std::getline(in, line)) {
      return false;
    }

    ss = std::stringstream(line);
    lineNum++;
    return true;
  };

  auto startsWithToken = [](const std::string& str, const std::string& prefix) {
    std::stringstream ss(str);
    std::string token;
    ss >> token;
    return token == prefix;
  };

  int first        = d_axis_sequence.x() - 1;
  int second       = d_axis_sequence.y() - 1;
  int third        = d_axis_sequence.z() - 1;
  double coords[3] = { 0.0, 0.0, 0.0 };

  // Parse STL file
  while (nextLine() && !startsWithToken(line, "endsolid")) {
    assertToken("facet");
    assertToken("normal");

    double normal_x, normal_y, normal_z;
    ss >> normal_x >> normal_y >> normal_z;
    Vector normal(normal_x, normal_y, normal_z);

    nextLine();

    assertToken("outer");
    assertToken("loop");

    std::vector<size_t> face;
    while (nextLine() && !startsWithToken(line, "endloop")) {
      assertToken("vertex");

      double position_x, position_y, position_z;
      ss >> position_x >> position_y >> position_z;

      coords[first]  = position_x;
      coords[second] = position_y;
      coords[third]  = position_z;
      d_points.push_back(Point(coords));

      face.push_back(d_points.size() - 1);
    }

    nextLine();
    assertToken("endfacet");

    // Orient face using normal
    Vector faceNormal = Cross(d_points[face[1]] - d_points[face[0]],
                              d_points[face[2]] - d_points[face[0]]);
    if (Dot(faceNormal, normal) < 0) {
      std::reverse(std::begin(face), std::end(face));
    }

    if (face.size() == 4) {
      d_triangles.push_back(IntVector(face[0], face[1], face[2]));
      d_triangles.push_back(IntVector(face[0], face[2], face[3]));
    } else {
      d_triangles.push_back(IntVector(face[0], face[1], face[2]));
    }
  }
}

/**
 * From:
 * https://github.com/nmwsharp/geometry-central/blob/master/src/surface/polygon_soup_mesh.cpp
 * May 23, 2020
 */
void
TriGeometryPiece::readMeshFromBinaryStlFile(std::ifstream in)
{
  auto parseVector3 = [&](std::ifstream& in) {
    char buffer[3 * sizeof(float)];
    in.read(buffer, 3 * sizeof(float));
    float* fVec = (float*)buffer;
    return Vector{ fVec[0], fVec[1], fVec[2] };
  };

  char header[80];
  char nTriangleChars[4];
  in.read(header, 80);
  in.read(nTriangleChars, 4);
  unsigned int* intPtr = (unsigned int*)nTriangleChars;
  size_t nTriangles    = *intPtr;

  int first        = d_axis_sequence.x() - 1;
  int second       = d_axis_sequence.y() - 1;
  int third        = d_axis_sequence.z() - 1;
  double coords[3] = { 0.0, 0.0, 0.0 };

  for (size_t iT = 0; iT < nTriangles; ++iT) {
    Vector normal = parseVector3(in);
    std::vector<size_t> face;
    for (size_t iV = 0; iV < 3; ++iV) {
      Vector point = parseVector3(in);

      coords[first]  = point[0];
      coords[second] = point[1];
      coords[third]  = point[2];
      d_points.push_back(Point(coords));

      face.push_back(d_points.size() - 1);
    }

    // Orient face using normal
    Vector faceNormal = Cross(d_points[face[1]] - d_points[face[0]],
                              d_points[face[2]] - d_points[face[0]]);
    if (Dot(faceNormal, normal) < 0) {
      std::reverse(std::begin(face), std::end(face));
    }

    if (face.size() == 4) {
      d_triangles.push_back(IntVector(face[0], face[1], face[2]));
      d_triangles.push_back(IntVector(face[0], face[2], face[3]));
    } else {
      d_triangles.push_back(IntVector(face[0], face[1], face[2]));
    }
    char dummy[2];
    in.read(dummy, 2);
  }
}

/**
 * From:
 * https://github.com/nmwsharp/geometry-central/blob/master/src/surface/polygon_soup_mesh.cpp
 * May 23, 2020
 */
// Mutate this mesh by merging vertices with identical floating point positions.
// Useful for loading .stl files, which don't contain information about which
// triangle corners meet at vertices.
void
TriGeometryPiece::mergeIdenticalVertices()
{
  std::vector<Point> compressedPositions;
  // Store mapping from original vertex index to merged vertex index
  std::vector<size_t> compressVertex;
  compressVertex.reserve(d_points.size());

  std::unordered_map<Point, size_t> canonicalIndex;

  for (size_t iV = 0; iV < d_points.size(); ++iV) {
    Point v = d_points[iV];
    auto it = canonicalIndex.find(v);

    // Check if vertex exists in map or not
    if (it == canonicalIndex.end()) {
      compressedPositions.push_back(v);
      size_t vecIndex   = compressedPositions.size() - 1;
      canonicalIndex[v] = vecIndex;
      compressVertex.push_back(vecIndex);
    } else {
      size_t vecIndex = it->second;
      compressVertex.push_back(vecIndex);
    }
  }

  d_points = std::move(compressedPositions);

  // Update face indices
  for (auto& triangle : d_triangles) {
    triangle[0] = compressVertex[triangle[0]];
    triangle[1] = compressVertex[triangle[1]];
    triangle[2] = compressVertex[triangle[2]];
  }
}

void
TriGeometryPiece::scaleTranslateReflect()
{
  // Scale, translate, and reflect
  std::for_each(d_points.begin(), d_points.end(), [&](Point& p) {
    p *= d_scale_factor;
    p *= d_reflect_vector;
    p += d_trans_vector;
  });
}

void
TriGeometryPiece::findBoundingBox()
{
  // Find the min and max points so that the bounding box can be determined.
  Point min(1e30, 1e30, 1e30), max(-1e30, -1e30, -1e30);
  for (auto point : d_points) {
    min = Min(point, min);
    max = Max(point, max);
  }
  Vector fudge(1.e-5, 1.e-5, 1.e-5);
  min   = min - fudge;
  max   = max + fudge;
  d_box = Box(min, max);

  std::cout << "Bounding box = " << d_box
            << " scale factor = " << d_scale_factor
            << " translation vector = " << d_trans_vector
            << " reflection vector = " << d_reflect_vector
            << " axis sequence = " << d_axis_sequence << "\n";
}

void
TriGeometryPiece::makePlanes()
{
  for (auto triangle : d_triangles) {
    Point pt[3];
    pt[0] = d_points[triangle.x()];
    pt[1] = d_points[triangle.y()];
    pt[2] = d_points[triangle.z()];
    Plane plane(pt[0], pt[1], pt[2]);
    d_planes.push_back(plane);
  }
}

void
TriGeometryPiece::makeTriangleBoxes()
{
  for (auto triangle : d_triangles) {
    Point pt[3];
    pt[0]     = d_points[triangle.x()];
    pt[1]     = d_points[triangle.y()];
    pt[2]     = d_points[triangle.z()];
    Point min = Min(Min(pt[0], pt[1]), Min(pt[1], pt[2]));
    Point max = Max(Max(pt[0], pt[1]), Max(pt[1], pt[2]));
    Box box(min, max);
    d_boxes.push_back(box);
  }
}

void
TriGeometryPiece::insideTriangle(Point& q, int num, int& NCS, int& NES) const
{
  // Pulled from makemesh.77.c
  //  Now we have to do the pt_in_pgon test to determine if ri is
  //  inside or outside the triangle.  Use Eric Haines idea in
  //  Essential Ray Tracing Algorithms, pg 53.  Don't have to worry
  //  about whether the intersection point is on the edge or vertex,
  //  cause the edge and/or vertex will be defined away.
  //

  // Now we project the ri and the vertices of the triangle onto
  // the dominant coordinate, i.e., the plane's normal largest
  // magnitude.
  //
  Vector plane_normal     = d_planes[num].normal();
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

  if (dominant_coord == -1) {
    std::cout << " dominant coordinate not found " << std::endl;
    throw InternalError("Dominant coordinate not found", __FILE__, __LINE__);
  }
  Point p[3];
  p[0] = d_points[d_triangles[num].x()];
  p[1] = d_points[d_triangles[num].y()];
  p[2] = d_points[d_triangles[num].z()];

  Triangle tri(p[0], p[1], p[2]);
  // bool inside = tri.inside(q);
  //   std::cout << "inside = " << inside << std::endl;

  // Now translate the points that make up the vertices of the triangle.
  Point trans_pt(0., 0., 0.), trans_vt[3];
  trans_vt[0] = Point(0., 0., 0.);
  trans_vt[1] = Point(0., 0., 0.);
  trans_vt[2] = Point(0., 0., 0.);

  if (dominant_coord == 1) {
    trans_pt.x(q.y());
    trans_pt.y(q.z());
    for (int i = 0; i < 3; i++) {
      trans_vt[i].x(p[i].y());
      trans_vt[i].y(p[i].z());
    }
  } else if (dominant_coord == 2) {
    trans_pt.x(q.x());
    trans_pt.y(q.z());
    for (int i = 0; i < 3; i++) {
      trans_vt[i].x(p[i].x());
      trans_vt[i].y(p[i].z());
    }
  } else if (dominant_coord == 3) {
    trans_pt.x(q.x());
    trans_pt.y(q.y());
    for (int i = 0; i < 3; i++) {
      trans_vt[i].x(p[i].x());
      trans_vt[i].y(p[i].y());
    }
  }

  // Now translate the intersecting point to the origin and the vertices
  // as well.
  for (int i = 0; i < 3; i++) {
    trans_vt[i] -= trans_pt.asVector();
  }

  int SH = 0, NSH = 0;
  double out_edge = 0.;

  if (trans_vt[0].y() < 0.0) {
    SH = -1;
  } else {
    SH = 1;
  }

  if (trans_vt[1].y() < 0.0) {
    NSH = -1;
  } else {
    NSH = 1;
  }

  if (SH != NSH) {
    if ((trans_vt[0].x() > 0.0) && (trans_vt[1].x() > 0.0)) {
      NCS += 1;
    } else if ((trans_vt[0].x() > 0.0) || (trans_vt[1].x() > 0.0)) {
      out_edge = (trans_vt[0].x() - trans_vt[0].y() *
                                      (trans_vt[1].x() - trans_vt[0].x()) /
                                      (trans_vt[1].y() - trans_vt[0].y()));
      if (out_edge == 0.0) {
        NES += 1;
        NCS += 1;
      }
      if (out_edge > 0.0) {
        NCS += 1;
      }
    }
    SH = NSH;
  }

  if (trans_vt[2].y() < 0.0) {
    NSH = -1;
  } else {
    NSH = 1;
  }

  if (SH != NSH) {
    if ((trans_vt[1].x() > 0.0) && (trans_vt[2].x() > 0.0)) {
      NCS += 1;
    } else if ((trans_vt[1].x() > 0.0) || (trans_vt[2].x() > 0.0)) {
      out_edge = (trans_vt[1].x() - trans_vt[1].y() *
                                      (trans_vt[2].x() - trans_vt[1].x()) /
                                      (trans_vt[2].y() - trans_vt[1].y()));
      if (out_edge == 0.0) {
        NES += 1;
        NCS += 1;
      }
      if (out_edge > 0.0) {
        NCS += 1;
      }
    }
    SH = NSH;
  }

  if (trans_vt[0].y() < 0.0) {
    NSH = -1;
  } else {
    NSH = 1;
  }

  if (SH != NSH) {
    if ((trans_vt[2].x() > 0.0) && (trans_vt[0].x() > 0.0)) {
      NCS += 1;
    }

    else if ((trans_vt[2].x() > 0.0) || (trans_vt[0].x() > 0.0)) {
      out_edge = (trans_vt[2].x() - trans_vt[2].y() *
                                      (trans_vt[0].x() - trans_vt[2].x()) /
                                      (trans_vt[0].y() - trans_vt[2].y()));
      if (out_edge == 0.0) {
        NES += 1;
        NCS += 1;
      }
      if (out_edge > 0.0) {
        NCS += 1;
      }
    }
    SH = NSH;
  }
}

void
TriGeometryPiece::scale(const double factor)
{
  Vector origin(0., 0., 0.);
  for (auto point : d_points) {
    origin = origin + point.asVector();
  }
  origin = origin / (static_cast<double>(d_points.size()));
  for (auto& point : d_points) {
    point = factor * (point - origin) + origin;
  }
}

double
TriGeometryPiece::surfaceArea() const
{
  double surfaceArea = 0.;
  for (auto triangle : d_triangles) {
    Point pt[3];
    pt[0] = d_points[triangle.x()];
    pt[1] = d_points[triangle.y()];
    pt[2] = d_points[triangle.z()];
    Vector v[2];
    v[0] = pt[0].asVector() - pt[1].asVector();
    v[1] = pt[2].asVector() - pt[1].asVector();

    Vector area = Cross(v[0], v[1]);
    surfaceArea += .5 * area.length();
  }
  return surfaceArea;
}

} // end namespace Uintah