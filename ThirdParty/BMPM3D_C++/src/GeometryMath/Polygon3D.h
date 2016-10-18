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

#ifndef __MATITI_POLYGON3D_H__
#define __MATITI_POLYGON3D_H__

#include <GeometryMath/Point3D.h>

#include <vector>

namespace BrMPM {

  class Polygon3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Polygon3D& poly);

  public:

    Polygon3D();
    Polygon3D(const std::vector<Point3D>& poly);
    Polygon3D(const Polygon3D& poly);
    ~Polygon3D();

    bool operator==(const Polygon3D& poly) const;

    unsigned int numVertices() const;
    void addVertex(const Point3D& pt);
    Polygon3D& operator+=(const Point3D& pt);
    
    const Point3D& vertex(const int& index) const;
    const Point3D& operator[](const int& index) const;

    std::vector<Point3D>::iterator begin() {return d_vertices.begin();}
    std::vector<Point3D>::iterator end() {return d_vertices.end();}

    std::vector<Point3D>::const_iterator begin() const {return d_vertices.begin();}
    std::vector<Point3D>::const_iterator end() const {return d_vertices.end();}

  private:

    std::vector<Point3D> d_vertices;

    // Prevent copying (for now)
    Polygon3D& operator=(const Polygon3D& poly);

  }; // end class Polygon3D
 
} // End namespace 

#endif //ifndef __MATITI_POLYGON3D_H__
