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

#ifndef MATITI_ELEMENT_H
#define MATITI_ELEMENT_H

#include <Pointers/NodeP.h>
#include <Containers/NodePArray.h>
#include <Geometry/Point3D.h>

namespace Matiti {
    
  class Element {
      
  public:

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Element& elem);

  public:

    Element();
    Element(const int& id, const NodePArray& nodes);
    ~Element();

    void initialize(const int id, const NodePArray& nodes);

    inline int id() const {return d_id;}
    const NodePArray& nodes() const {return d_nodes;}
    inline int numNodes() const {return d_nodes.size();}

    void computeGeometry2D(double& area, double& xlength, double& ylength) const;

    void computeVolume();
    void computeVolume2D();
    void computeVolume3D();

    double volume() const {return d_volume;}

  protected:

    int d_id;
    double d_volume;
    NodePArray d_nodes; 

  private:

    double computeVolumeTetrahedron(const Point3D& p0, const Point3D& p1, const Point3D& p2,
                                    const Point3D& p3) const;

    // Prevent copy construction and operator=
    Element(const Element& element);
    Element& operator=(const Element& element);

  }; // end class

} // end namespace

#endif


