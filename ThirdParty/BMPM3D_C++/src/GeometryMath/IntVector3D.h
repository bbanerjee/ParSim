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

/*
 * IntIntVector3D.h
 *
 *  Created on: 22/10/2013
 *      Author: banerjee
 */

#ifndef INTVECTOR3D_H_
#define INTVECTOR3D_H_

#include <Types.h>

#include <iostream>
#include <limits>

namespace BrMPM {

  class IntVector3D {

  public:

    friend std::ostream& operator<<(std::ostream& os, const BrMPM::IntVector3D& p);

  public:

    IntVector3D(): d_vec({{std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::max(),
                        std::numeric_limits<int>::max()}}) {}
    IntVector3D(int val) : d_vec({{val, val, val}}) {}
    IntVector3D(int x, int y, int z): d_vec({{x, y, z}}) {}
    IntVector3D(const IntVector3D& vec);

    ~IntVector3D() {}

    //  **WARNING** Not checking index.  Do checking outside.
    int operator[](int index) const {return d_vec[index];}
    int& operator[](int index) {return d_vec[index];}

    void set(const int& val) {d_vec[0] = val; d_vec[1] = val; d_vec[2] = val;}
    void x(const int xx) {d_vec[0] = xx;}
    int x() const {return d_vec[0];}
    void y(const int yy) {d_vec[1] = yy;}
    int y() const {return d_vec[1];}
    void z(const int zz) {d_vec[2] = zz;}
    int z() const {return d_vec[2];}

    bool operator==(const IntVector3D& vec) const;
    bool operator!=(const IntVector3D& vec) const;

    void operator=(const IntVector3D& vec);

    IntVector3D operator*(const int fac) const;
    IntVector3D operator*(const IntVector3D& vec) const;

    IntVector3D operator+(const IntVector3D& vec) const;
    IntVector3D operator-(const IntVector3D& vec) const;

    IntVector3D& operator*=(const int);
    IntVector3D& operator+=(const IntVector3D& vec);
    IntVector3D& operator-=(const IntVector3D& vec);

    void reset() {d_vec[0] = 0; d_vec[1] = 0; d_vec[2] = 0;}

    int max() const;
    int min() const;

  private:

    IntArray3 d_vec;

  }; // end class

} // end namespace

namespace BrMPM {

 IntVector3D min(const IntVector3D& v1, const IntVector3D& v2);
 IntVector3D max(const IntVector3D& v1, const IntVector3D& v2);

} /* namespace BrMPM */
#endif /* INTVECTOR3D_H_ */
