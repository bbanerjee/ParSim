#ifndef __EMU2DC_CRACK_H__
#define __EMU2DC_CRACK_H__

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp> 

#include <string>

namespace bg = boost::geometry;
typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D;
typedef bg::model::linestring<Point3D> LineString;

namespace Emu2DC {

  
 class Crack {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Crack& crack);

  public:

   Crack();
   virtual ~Crack();

   void initialize(const Uintah::ProblemSpecP& ps);

   void triangulate();

  private:

   LineString d_lineString;

   // Bit of code to parse a Uintah::Vector input
   void parseVector(const std::string& stringValue, SCIRun::Vector& value);
   void checkForInputError(const std::string& stringValue);
    
 }; // end class

}; // end namespace

#endif

