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
      //d_boundary.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
      d_boundary.addVertex(Point3D(point[0], point[1], point[2]));  
      while ((node = node->findNextBlock("point"))) {
        parseVector(node->getNodeValue(), point);
        //d_boundary.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
        d_boundary.addVertex(Point3D(point[0], point[1], point[2]));  
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

void
Crack::breakBonds(NodeP& node, const NodePArray& family) const
{
  // Get node location
  Array3 node_pos = node->position();  
  Point3D seg_start(node_pos[0], node_pos[1], node_pos[2]);  

  // Loop through node family
  for (auto iter = family.begin(); iter != family.end(); ++iter) {

    // Get family node location
    Array3 fam_pos = (*iter)->position();
    Point3D seg_end(fam_pos[0], fam_pos[1], fam_pos[2]);  

    // Loop through triangles
    auto o_iter = d_origin.begin();
    auto d_iter = d_destination.begin();
    auto a_iter = d_apex.begin();
    for (; o_iter != d_origin.end(); ++o_iter, ++d_iter, ++a_iter) {

      // Get the three vertices of the triangle
      const Point3D& node1 = d_boundary[*o_iter];
      const Point3D& node2 = d_boundary[*d_iter];
      const Point3D& node3 = d_boundary[*a_iter];
    } // end triangle loop
  } // end family loop
}

bool 
Crack::intersectSegmentWithTriangle(const Point3D& start, const Point3D& end,
                                    const Point3D& orig, const Point3D& dest, const Point3D& apex) const
{
  // Triangle edge vectors
  return false;  
}

void 
Crack::triangulate()
{
  // **WARNING** A number of checks are needed at this stage.  Skipping and assuming 
  // that all points lie on the xy plane.
  tpp::Delaunay::Point tempP;
  std::vector<tpp::Delaunay::Point> vec;
  for (auto iter = d_boundary.begin(); iter != d_boundary.end(); ++iter) {
    tempP[0] = (*iter).x();
    tempP[1] = (*iter).y();
    vec.push_back(tempP);
    //std::cout << "tempP[0] = " << tempP[0] << " tempP[1] = " << tempP[1] << std::endl;
  }

  // Triangulate
  tpp::Delaunay delobject(vec);
  int flag = delobject.Triangulate();
  if (flag) {
    std::ostringstream out;
    out << "Could not triangulate crack " << flag << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }

  // Save triangles
  for (tpp::Delaunay::fIterator iter = delobject.fbegin(); iter != delobject.fend(); ++iter) {
    d_origin.emplace_back(delobject.Org(iter));
    d_destination.emplace_back(delobject.Dest(iter));
    d_apex.emplace_back(delobject.Apex(iter));
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
    //d_boundary.push_back(boost::geometry::make<Point3D>(xcoord, ycoord, zcoord));  
    d_boundary.addVertex(Point3D(xcoord, ycoord, zcoord));  
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
