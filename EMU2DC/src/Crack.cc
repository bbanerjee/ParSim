#include <Crack.h>
#include <Exception.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Emu2DC;

Crack::Crack()
{
}

Crack::~Crack()
{
}

void Crack::initialize(const Uintah::ProblemSpecP& ps)
{
  // Check for the <LineString> block
  Uintah::ProblemSpecP crack_ps = ps->findBlock("LineString");
  if (!crack_ps) return;

  // Get the crack boundary points and create linestring
  Uintah::Vector point(0.0, 0.0, 0.0); 
  Uintah::ProblemSpecP node = crack_ps->findBlock("point");
  parseVector(node->getNodeValue(), point);
  d_lineString.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
  int ii =1;
  while ((node = node->findNextBlock("point"))) {
    parseVector(node->getNodeValue(), point);
    d_lineString.push_back(boost::geometry::make<Point3D>(point[0], point[1], point[2]));  
  }
}

void 
Crack::triangulate()
{
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
    out << boost::geometry::dsv(crack.d_lineString) << std::endl;
    return out;
  }
}
