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

#ifndef __MATITI_DOMAIN_H__
#define __MATITI_DOMAIN_H__

#include <Types/Types.h>
#include <Pointers/BodySP.h>
#include <Geometry/Point3D.h>
#include <Containers/VelocityBCSPArray.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Geometry/Vector.h>
#include <iostream>

namespace Matiti {

  class Domain {

  public:  

    friend std::ostream& operator<<(std::ostream& out, const Matiti::Domain& domain);

  public:  

    Domain() ;
    virtual ~Domain();

    void clone(const Domain& domain);

    Domain(const Point3D& lower, const Point3D& upper);

    Domain(const Point3D& lower, const Point3D& upper, const IntArray3& numCells);
    
    Domain(const Point3D& lower, const Point3D& upper, const Uintah::Vector& cellSize);

    virtual void initialize(const Uintah::ProblemSpecP& ps);

    const Point3D& lower() const;
    const Point3D& upper() const;
    const double& xrange() const;
    const double& yrange() const;
    const double& zrange() const;
    const IntArray3& numCells() const;
    const Uintah::Vector& cellSize() const;
    const double totalCells() const;

    void findCellIndex(const Point3D& point,
                       IntArray3& cell) const;
    void findCellIndex(const long64& cell_key,
                       IntArray3& cell) const;

    bool inside(const Point3D& point) const;

    void applyVelocityBC(BodySP body) const;

    bool intersection(const Point3D& point, const Vector3D& ray,
                      Point3D& hitPoint) const;

  private:

    Point3D d_lower;
    Point3D d_upper;

    double d_xrange;
    double d_yrange;
    double d_zrange;

    IntArray3 d_num_cells;
    Uintah::Vector d_cell_size;

    VelocityBCSPArray d_vel_BC;

    // Don't allow copy
    Domain(const Domain& dom);
    Domain& operator=(const Domain& dom);

  };  // end class
}  // end namespace
#endif
