#ifndef __EMU2DC_CRACK_H__
#define __EMU2DC_CRACK_H__

#include <NodeP.h>
#include <NodePArray.h>

#include <Geometry/Point3D.h>
#include <Geometry/Vector3D.h>
#include <Geometry/Polygon3D.h>

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Geometry/Vector.h>

#include <string>

// #include <boost/geometry.hpp>
// #include <boost/geometry/geometries/linestring.hpp> 

// namespace bg = boost::geometry;
// typedef bg::model::point<double, 3, bg::cs::cartesian> Point3D;
// typedef bg::model::linestring<Point3D> Polygon3D;

namespace Emu2DC {

  
 class Crack {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Emu2DC::Crack& crack);

  public:

   Crack();
   virtual ~Crack();

   void initialize(const Uintah::ProblemSpecP& ps);

   void breakBonds(const NodeP& node, NodePArray& family) const;

  protected:

   // Read the crack input file
   void readCrackFile(const std::string& fileName);
 
   // Triangulate the crack (only xy cracks for now)
   void triangulate();

   // Check if a bond intersects a triangle
   bool intersectSegmentWithTriangle(const Point3D& start, const Point3D& end,
                                     const Point3D& orig, const Point3D& dest, 
                                     const Point3D& apex) const;

   // Bit of code to parse a Uintah::Vector input
   void parseVector(const std::string& stringValue, SCIRun::Vector& value);
   void checkForInputError(const std::string& stringValue);

  private:

   // Boundary nodes that describe the crack. The nodes are numbered starting from zero.
   Polygon3D d_boundary;

   // Triangle element connectivity.  Each triangle has three nodes: the origin, the
   // destination, and the apex.  
   std::vector<int> d_origin;
   std::vector<int> d_destination;
   std::vector<int> d_apex;

 }; // end class

}; // end namespace

#endif

