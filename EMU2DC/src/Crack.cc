#include <Crack.h>
#include <Node.h>
#include <Exception.h>
#include <Types.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Containers/StringUtil.h>

#include <del_interface.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

using namespace Emu2DC;

Crack::Crack()
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
    Uintah::Vector point(0.0, 0.0, 0.0); 
    Uintah::ProblemSpecP node = crack_ps->findBlock("point");
    if (node) { 

      // Read points from the input file
      parseVector(node->getNodeValue(), point);

      // Counter for crack boundary points
      int counter = 1;

      //d_boundary.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
      d_boundary.addVertex(Point3D(point[0], point[1], point[2]));  

      while ((node = node->findNextBlock("point"))) {
        parseVector(node->getNodeValue(), point);

        ++counter;

        //d_boundary.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
        d_boundary.addVertex(Point3D(point[0], point[1], point[2]));  
      }
      if (counter < 3) {
	std::ostringstream out;
        out << "**ERROR** A crack boundary cannot have less than three points" << std::endl;
        throw Exception(out.str(), __FILE__, __LINE__);
      } 
    } else {
      // Read points from a separate file containing x, y, z
      std::string crack_input_file;
      Uintah::ProblemSpecP file = crack_ps->get("crack_boundary_file", crack_input_file);
      std::cout << "Input crack file = " << crack_input_file << std::endl; 
      readCrackFile(crack_input_file); 
    }

    // Triangulate the crack and save elements
    triangulate();
  }
}

// Bit of code to parse a Uintah::Vector input
void 
Crack::parseVector(const std::string& stringValue, Uintah::Vector& value)
{
  // Parse out the [num,num,num]
  // Now pull apart the stringValue
  std::string::size_type i1 = stringValue.find("[");
  std::string::size_type i2 = stringValue.find_first_of(",");
  std::string::size_type i3 = stringValue.find_last_of(",");
  std::string::size_type i4 = stringValue.find("]");

  std::string x_val(stringValue,i1+1,i2-i1-1);
  std::string y_val(stringValue,i2+1,i3-i2-1);
  std::string z_val(stringValue,i3+1,i4-i3-1);
    
  checkForInputError( x_val );
  checkForInputError( y_val );
  checkForInputError( z_val );
    
  value.x(std::stod(x_val));
  value.y(std::stod(y_val));
  value.z(std::stod(z_val));
}

void
Crack::checkForInputError(const std::string & stringValue)
{
  std::string validChars(" -+.0123456789eE");
  std::string::size_type  pos = stringValue.find_first_not_of(validChars);
  if (pos != std::string::npos){
    std::ostringstream warn;
    warn << "Bad Float string: Found '"<< stringValue[pos]
         << "' inside of \""<< stringValue << "\" at position " << pos << "\n";
    throw Exception(warn.str(), __FILE__, __LINE__);
  }
  //__________________________________
  // check for two or more "."
  std::string::size_type p1 = stringValue.find_first_of(".");
  std::string::size_type p2 = stringValue.find_last_of(".");
  if (p1 != p2){
    std::ostringstream warn;
    warn << "Input file error: I found two (..) "
         << "inside of "<< stringValue << "\n";
    throw Exception(warn.str(), __FILE__, __LINE__);
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

}

void
Crack::breakBonds(const NodeP& node, NodePArray& family) const
{
  // Get node location
  Array3 node_pos = node->position();  
  Point3D seg_start(node_pos[0], node_pos[1], node_pos[2]);  

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
        [&](const NodeP& fam_node)
        {
          Array3 fam_pos = fam_node->position();
          Point3D seg_end(fam_pos[0], fam_pos[1], fam_pos[2]);  
          return Crack::intersectSegmentWithTriangle(seg_start, seg_end, orig, dest, apex); 
        };
    family.erase(std::remove_if(family.begin(), family.end(), lambda_func), family.end());

  } // end triangle loop
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

namespace Emu2DC {

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
