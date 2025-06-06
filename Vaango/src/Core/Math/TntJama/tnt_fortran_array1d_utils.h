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
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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
*
* Template Numerical Toolkit (TNT)
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*/

#ifndef TNT_FORTRAN_ARRAY1D_UTILS_H
#define TNT_FORTRAN_ARRAY1D_UTILS_H

#include <iostream>

namespace TNT {

/**
        Write an array to a character outstream.  Output format is one that can
        be read back in via the in-stream operator: one integer
        denoting the array dimension (n), followed by n elements,
        one per line.

*/
template <class T>
std::ostream&
operator<<(std::ostream& s, const Fortran_Array1D<T>& A)
{
  int N = A.dim1();

  s << N << "\n";
  for (int j = 1; j <= N; j++) {
    s << A(j) << "\n";
  }
  s << "\n";

  return s;
}

/**
        Read an array from a character stream.  Input format
        is one integer, denoting the dimension (n), followed
        by n whitespace-separated elments.  Newlines are ignored

        <p>
        Note: the array being read into references new memory
        storage. If the intent is to fill an existing conformant
        array, use <code> cin >> B;  A.inject(B) ); </code>
        instead or read the elements in one-a-time by hand.

        @param s the charater to read from (typically <code>std::in</code>)
        @param A the array to read into.
*/
template <class T>
std::istream&
operator>>(std::istream& s, Fortran_Array1D<T>& A)
{
  int N;
  s >> N;

  Fortran_Array1D<T> B(N);
  for (int i = 1; i <= N; i++)
    s >> B(i);
  A = B;
  return s;
}

template <class T>
Fortran_Array1D<T>
operator+(const Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() != n)
    return Fortran_Array1D<T>();

  else {
    Fortran_Array1D<T> C(n);

    for (int i = 1; i <= n; i++) {
      C(i) = A(i) + B(i);
    }
    return C;
  }
}

template <class T>
Fortran_Array1D<T>
operator-(const Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() != n)
    return Fortran_Array1D<T>();

  else {
    Fortran_Array1D<T> C(n);

    for (int i = 1; i <= n; i++) {
      C(i) = A(i) - B(i);
    }
    return C;
  }
}

template <class T>
Fortran_Array1D<T> operator*(const Fortran_Array1D<T>& A,
                             const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() != n)
    return Fortran_Array1D<T>();

  else {
    Fortran_Array1D<T> C(n);

    for (int i = 1; i <= n; i++) {
      C(i) = A(i) * B(i);
    }
    return C;
  }
}

template <class T>
Fortran_Array1D<T>
operator/(const Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() != n)
    return Fortran_Array1D<T>();

  else {
    Fortran_Array1D<T> C(n);

    for (int i = 1; i <= n; i++) {
      C(i) = A(i) / B(i);
    }
    return C;
  }
}

template <class T>
Fortran_Array1D<T>&
operator+=(Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() == n) {
    for (int i = 1; i <= n; i++) {
      A(i) += B(i);
    }
  }
  return A;
}

template <class T>
Fortran_Array1D<T>&
operator-=(Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() == n) {
    for (int i = 1; i <= n; i++) {
      A(i) -= B(i);
    }
  }
  return A;
}

template <class T>
Fortran_Array1D<T>&
operator*=(Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() == n) {
    for (int i = 1; i <= n; i++) {
      A(i) *= B(i);
    }
  }
  return A;
}

template <class T>
Fortran_Array1D<T>&
operator/=(Fortran_Array1D<T>& A, const Fortran_Array1D<T>& B)
{
  int n = A.dim1();

  if (B.dim1() == n) {
    for (int i = 1; i <= n; i++) {
      A(i) /= B(i);
    }
  }
  return A;
}

} // namespace TNT

#endif
