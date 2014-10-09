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

#ifndef __MATITI_BOX3D_H__
#define __MATITI_BOX3D_H__

#include <Geometry/Point3D.h>

#include <iostream>

namespace Matiti {

  class Box3D 
  {
  public:

    friend std::ostream& operator<<(std::ostream& os, const Box3D& box);

  public:

    Box3D();
    Box3D(const Box3D& box);
    Box3D(const Point3D& lower, const Point3D& upper);

    ~Box3D();

    Box3D& operator=(const Box3D& box);

    bool overlaps(const Box3D& box, double epsilon=1.0e-6) const;

    bool contains(const Point3D& pt) const;

    bool isDegenerate();

    Point3D lower() const;
    Point3D upper() const;

  private:

    void correctBoundingBox();

  private:

    Point3D d_lower;
    Point3D d_upper;

  }; // end class Box3D
 
} // End namespace 

#endif //ifndef __MATITI_BOX3D_H__
