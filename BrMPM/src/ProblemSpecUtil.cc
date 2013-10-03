#include <ProblemSpecUtil.h>
#include <Geometry/Point3D.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Exception.h>
#include <iostream>

namespace Matiti_ProblemSpecUtil
{
  // A routine to read in the boundary of a two-dimensional region
  void
  readBoundary(Uintah::ProblemSpecP& ps, Matiti::Polygon3D& boundary)
  {
    // Read points from the input file
    SCIRun::Vector point(0.0, 0.0, 0.0); 
    Matiti_ProblemSpecUtil::parseVector(ps->getNodeValue(), point);

    // Counter for area boundary points
    int counter = 1;

    boundary.addVertex(Matiti::Point3D(point[0], point[1], point[2]));  

    while ((ps = ps->findNextBlock("point"))) {
      Matiti_ProblemSpecUtil::parseVector(ps->getNodeValue(), point);

      ++counter;

      boundary.addVertex(Matiti::Point3D(point[0], point[1], point[2]));  
    }
    if (counter < 3) {
      std::ostringstream out;
      out << "**ERROR** A area boundary cannot have less than three points" << std::endl;
      throw Matiti::Exception(out.str(), __FILE__, __LINE__);
    } 
  }

  // Bit of code to parse a Uintah::Vector input
  void 
  parseVector(const std::string& stringValue, SCIRun::Vector& value)
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
  checkForInputError(const std::string & stringValue)
  {
    std::string validChars(" -+.0123456789eE");
    std::string::size_type  pos = stringValue.find_first_not_of(validChars);
    if (pos != std::string::npos){
      std::ostringstream warn;
      warn << "Bad Float string: Found '"<< stringValue[pos]
           << "' inside of \""<< stringValue << "\" at position " << pos << "\n";
      throw Matiti::Exception(warn.str(), __FILE__, __LINE__);
    }
    //__________________________________
    // check for two or more "."
    std::string::size_type p1 = stringValue.find_first_of(".");
    std::string::size_type p2 = stringValue.find_last_of(".");
    if (p1 != p2){
      std::ostringstream warn;
      warn << "Input file error: I found two (..) "
           << "inside of "<< stringValue << "\n";
      throw Matiti::Exception(warn.str(), __FILE__, __LINE__);
    }
  }
}
