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

#ifndef __MATITI_CRACK_H__
#define __MATITI_CRACK_H__

#include <Pointers/NodeP.h>
#include <Containers/BondPArray.h>

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

namespace Matiti {

  
 class Crack {

  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Crack& crack);

  public:

   Crack();
   virtual ~Crack();

   void initialize(const Uintah::ProblemSpecP& ps);

   void breakBonds(const NodeP node, BondPArray& family) const;

  protected:

   // Read the crack input file
   void readCrackFile(const std::string& fileName);
 
   // Triangulate the crack (only xy cracks for now)
   void triangulate();

   // Check if a bond intersects a triangle
   bool intersectSegmentWithTriangle(const Point3D& start, const Point3D& end,
                                     const Point3D& orig, const Point3D& dest, 
                                     const Point3D& apex) const;

  private:

   // Boundary nodes that describe the crack. The nodes are numbered starting from zero.
   Polygon3D d_boundary;
   double d_factor; // Coordinate scaling factor

   // Triangle element connectivity.  Each triangle has three nodes: the origin, the
   // destination, and the apex.  
   std::vector<int> d_origin;
   std::vector<int> d_destination;
   std::vector<int> d_apex;

 }; // end class

}; // end namespace

#endif

