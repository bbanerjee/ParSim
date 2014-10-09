/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <Core/Crack.h>
#include <Core/Node.h>
#include <Core/Bond.h>
#include <Core/Exception.h>
#include <Types/Types.h>
#include <Geometry/Point3D.h>
#include <InputOutput/ProblemSpecUtil.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Containers/StringUtil.h>

#include <del_interface.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

using namespace Matiti;

Crack::Crack()
  : d_factor(1.0)
{
}

Crack::~Crack()
{
}

void Crack::initialize(const Uintah::ProblemSpecP& ps)
{
  // Check for the <BoundaryPoints> block
  Uintah::ProblemSpecP crack_ps = ps->findBlock("BoundaryPoints");
  if (crack_ps) {

    // Get the crack boundary points and create linestring
    Uintah::ProblemSpecP node = crack_ps->findBlock("point");
    if (node) { 

      Matiti_ProblemSpecUtil::readBoundary(node, d_boundary);

    } else {
      // Read points from a separate file containing x, y, z
      std::string crack_input_file;
      Uintah::ProblemSpecP file = crack_ps->get("crack_boundary_file", crack_input_file);
      std::cout << "Input crack file = " << crack_input_file << std::endl; 
      d_factor = 1.0;
      crack_ps->get("coord_scaling_factor", d_factor);
      readCrackFile(crack_input_file); 
    }

    // Triangulate the crack and save elements
    triangulate();
  }
}

// Read the crack input file
void
Crack::readCrackFile(const std::string& fileName) 
{
  // Try to open file
  std::ifstream file(fileName);
  if (!file.is_open()) {
    std::string out = "Could not open crack input file " + fileName + " for reading \n";
    throw Exception(out, __FILE__, __LINE__);
  }

  // Read file
  int counter = 0;
  std::string line;
  while (std::getline(file, line)) {

    // Ignore empty lines
    if (line.empty()) continue;

    // erase white spaces from the beginning of line
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
         std::not1(std::ptr_fun<int, int>(std::isspace))));
    
    // Skip comment lines
    if (line[0] == '#') continue;

    // Read the data
    std::istringstream data_stream(line);
    double xcoord, ycoord, zcoord;
    if (!(data_stream >> xcoord >> ycoord >> zcoord)) {
      throw Exception("Could not read crack input data stream", __FILE__, __LINE__);
    }
    xcoord *= d_factor;
    ycoord *= d_factor;
    zcoord *= d_factor;

    //std::cout << "x" << counter << "=" << xcoord
    //          << " y" << counter << "=" << ycoord
    //          << " z" << counter << "=" << zcoord << std::endl;

    // Save the data
    ++counter;
    //d_boundary.push_back(boost::geometry::make<Point3D>(xcoord, ycoord, zcoord));  
    d_boundary.addVertex(Point3D(xcoord, ycoord, zcoord));  
  }
  if (counter < 3) {
    std::ostringstream out;
    out << "**ERROR** A crack boundary cannot have less than two points" << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  } 
  std::cout << "Completed reading crack input file: num crack boundary points = " << d_boundary.numVertices() << std::endl;
}

void 
Crack::triangulate()
{
  // **WARNING** A number of checks are needed at this stage.  Skipping and assuming 
  // that all points lie on the xy plane.

  // Check whether the polygon is on the xz or yz plane
  //   For simplicity just check the first three points
  int counter = 0;
  double xx[3], yy[3];
  for (auto iter = d_boundary.begin(); iter != d_boundary.end(); ++iter) {
    xx[counter] = (*iter).x();
    yy[counter] = (*iter).y();
    //std::cout << "x" << counter << "=" << (*iter).x() << ":" << xx[counter]
    //          << " y" << counter << "=" << (*iter).y() << ":" << yy[counter]
    //          << " z" << counter << "=" << (*iter).z() << std::endl;
    ++counter;
    if (counter > 2) break;
  }
  
  // Points are on yz plane
  bool yz_plane = false;
  if ((xx[0] == xx[1]) && (xx[0] == xx[2])) yz_plane = true;

  // Points are on xz plane
  bool xz_plane = false;
  if ((yy[0] == yy[1]) && (yy[0] == yy[2])) xz_plane = true;

  //std::cout << "yz_plane = " << std::boolalpha << yz_plane << " xz_plane = " << xz_plane << std::endl;

  if (yz_plane && xz_plane) {
    std::ostringstream out;
    out << "**ERROR** Crack cannot be both on xz and yz planes." << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Project points so that they lie on xy plane
  tpp::Delaunay::Point tempP;
  std::vector<tpp::Delaunay::Point> vec;
  for (auto iter = d_boundary.begin(); iter != d_boundary.end(); ++iter) {
    if (yz_plane) {
      tempP[0] = (*iter).y();
      tempP[1] = (*iter).z();
    } else if (xz_plane) {
      tempP[0] = (*iter).x();
      tempP[1] = (*iter).z();
    } else {
      tempP[0] = (*iter).x();
      tempP[1] = (*iter).y();
    }
    vec.push_back(tempP);
    //std::cout << "tempP[0] = " << tempP[0] << " tempP[1] = " << tempP[1] << std::endl;
  }

  // Triangulate
  tpp::Delaunay delobject(vec);
  int flag = delobject.Triangulate();
  if (flag) {
    std::ostringstream out;
    out << "Could not triangulate crack " << flag << std::endl;
    out << *this << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Save triangles
  for (tpp::Delaunay::fIterator iter = delobject.fbegin(); iter != delobject.fend(); ++iter) {
    d_origin.emplace_back(delobject.Org(iter));
    d_destination.emplace_back(delobject.Dest(iter));
    d_apex.emplace_back(delobject.Apex(iter));
  }

  std::cout << "Crack triangulation complete: number of triangles = " << d_origin.size() << std::endl;
}

void
Crack::breakBonds(const NodeP node, BondPArray& family) const
{
  // Get node location
  const Point3D& seg_start = node->position();  

  if (!(family.size() > 0)) {
    std::ostringstream out;
    out << "**ERROR** Number of initial bonds is zero for node " << node->getID();
    throw Exception(out.str(), __FILE__, __LINE__);
  }
   //std::cout << "Node = " << node->getID() << " Num bonds before = " << family.size();
  // Loop through triangles
  auto o_iter = d_origin.begin();
  auto d_iter = d_destination.begin();
  auto a_iter = d_apex.begin();
  for (; o_iter != d_origin.end(); ++o_iter, ++d_iter, ++a_iter) {

    // Get the three vertices of the triangle
    const Point3D& orig = d_boundary[*o_iter];
    const Point3D& dest = d_boundary[*d_iter];
    const Point3D& apex = d_boundary[*a_iter];

    // Get family node location
    // Get intersection of segment with triangle and remove if true
    auto lambda_func = 
        [&](const BondP& bond)
        {
          const Point3D& seg_end = bond->second()->position();
          return Crack::intersectSegmentWithTriangle(seg_start, seg_end, orig, dest, apex); 
        };
    family.erase(std::remove_if(family.begin(), family.end(), lambda_func), family.end());

  } // end triangle loop
  //std::cout << " after = " << family.size() << std::endl;
}

// Algorithm from: http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle%28%29
bool 
Crack::intersectSegmentWithTriangle(const Point3D& seg_start, const Point3D& seg_end,
                                    const Point3D& tri_orig, const Point3D& tri_dest, 
                                    const Point3D& tri_apex) const
{
  // Triangle edge vectors
  Vector3D uu(tri_orig, tri_dest);
  Vector3D vv(tri_orig, tri_apex);
  Vector3D normal = uu.cross(vv);

  // Check if triangle is degenerate
  if (normal == Vector3D(0.0,0.0,0.0)) return false;

  // Ray direction
  Vector3D segment(seg_start, seg_end);
  Vector3D vec0(tri_orig, seg_start);
  double normal_dot_vec0 = -normal.dot(vec0);
  double normal_dot_seg = normal.dot(segment);

  // Check is segment is parallel to plane of triangle (**WARNING** Hardcoded tolerance)
  if (std::abs(normal_dot_seg) < 1.0e-8) return false;

  // Check if segment intersects the plane of the triangle
  double tt = normal_dot_vec0/normal_dot_seg;
  if (tt < 0.0 || tt > 1.0) return false;

  // Find point of intersection of segment and plane
  Point3D intersect = seg_start + segment*tt;

  // Find if the intersection point is inside the triangle
  double u_dot_u = uu.dot(uu);
  double u_dot_v = uu.dot(vv);
  double v_dot_v = vv.dot(vv);
  Vector3D ww(tri_orig, intersect);
  double w_dot_u = ww.dot(uu);
  double w_dot_v = ww.dot(vv);
  double denom = u_dot_v*u_dot_v - u_dot_u*v_dot_v;
  if (std::abs(denom) < 1.0e-8) return false;

  // Barycentric coordinates
  double ss = (u_dot_v*w_dot_v - v_dot_v*w_dot_u)/denom;
  if (ss < 0.0 || ss > 1.0) return false;

  double rr = (u_dot_v*w_dot_u - u_dot_u*w_dot_v)/denom;
  if (rr < 0.0 || (rr+ss) > 1.0) return false;

  return true;  
}

namespace Matiti {

  std::ostream& operator<<(std::ostream& out, const Crack& crack)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Crack geometry points:" << std::endl;
    for (auto iter = crack.d_boundary.begin(); iter != crack.d_boundary.end(); ++iter) {
      out << *iter ;
    }
    out << std::endl;
    // out << boost::geometry::dsv(crack.d_boundary) << std::endl;
    int num_elem = crack.d_origin.size();
    for (int ii = 0; ii < num_elem; ii++) {
      out << "Element " << ii << " :[" << crack.d_origin[ii] << ", "
              << crack.d_destination[ii] << ", " << crack.d_apex[ii] << "]" << std::endl;
    }
    return out;
  }
}
