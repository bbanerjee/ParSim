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

  protected:

   // Read the crack input file
   void readCrackFile(const std::string& fileName);
 
   // Triangulate the crack
   void triangulate();

   // Bit of code to parse a Uintah::Vector input
   void parseVector(const std::string& stringValue, SCIRun::Vector& value);
   void checkForInputError(const std::string& stringValue);

  private:

   // Boundary nodes that describe the crack. The nodes are numbered starting from zero.
   LineString d_boundary;

   // Triangle element connectivity.  Each triangle has three nodes: the origin, the
   // destination, and the apex.  
   std::vector<int> d_origin;
   std::vector<int> d_destination;
   std::vector<int> d_apex;

 }; // end class

}; // end namespace

#endif

