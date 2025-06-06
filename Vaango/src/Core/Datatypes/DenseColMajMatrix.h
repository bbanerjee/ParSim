/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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
 *  DenseColMajMatrix.h:  DenseColMaj matrices
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   October 1994
 *
 */

#ifndef __VAANGO_CORE_DATATYPE_DENSE_COL_MAJOR_MATRIX_H__
#define __VAANGO_CORE_DATATYPE_DENSE_COL_MAJOR_MATRIX_H__ 1

#include <Core/Datatypes/Matrix.h>

#include <Core/Math/MiscMath.h>

namespace Uintah {

class DenseColMajMatrix : public Matrix
{
  double* dataptr_;

public:
  //! Constructors
  DenseColMajMatrix();
  DenseColMajMatrix(int r, int c);
  DenseColMajMatrix(const DenseColMajMatrix&);

  //! Destructor
  virtual ~DenseColMajMatrix();

  //! Public member functions
  virtual DenseColMajMatrix*
  clone();

  DenseColMajMatrix&
  operator=(const DenseColMajMatrix&);

  virtual DenseMatrix*
  dense();

  virtual SparseRowMatrix*
  sparse();

  virtual ColumnMatrix*
  column();

  virtual DenseColMajMatrix*
  dense_col_maj();

  virtual double*
  get_data_pointer();

  virtual size_t
  get_data_size();

  //! slow setters/getter for polymorphic operations
  virtual void
  zero();

  virtual double
  get(int r, int c) const;

  virtual void
  put(int r, int c, double val);

  virtual void
  add(int r, int c, double val);

  virtual void
  getRowNonzerosNoCopy(int r,
                       int& size,
                       int& stride,
                       int*& cols,
                       double*& vals);

  virtual DenseColMajMatrix*
  transpose();

  virtual void
  mult(const ColumnMatrix& x,
       ColumnMatrix& b,
       int& flops,
       int& memrefs,
       int beg   = -1,
       int end   = -1,
       int spVec = 0) const;

  virtual void
  mult_transpose(const ColumnMatrix& x,
                 ColumnMatrix& b,
                 int& flops,
                 int& memrefs,
                 int beg   = -1,
                 int end   = -1,
                 int spVec = 0) const;

  virtual MatrixHandle
  submatrix(int r1, int c1, int r2, int c2);

  double
  sumOfCol(int);

  double
  sumOfRow(int);

  //! fast accessors
  inline double&
  iget(int r, int c)
  {
    return dataptr_[c * nrows_ + r];
  };

  //! fast accessors
  inline const double&
  iget(int r, int c) const
  {
    return dataptr_[c * nrows_ + r];
  };

  //! Throws an assertion if not square
  double
  determinant();

  static DenseColMajMatrix*
  identity(int size);

  virtual void
  print() const;

  virtual void
  print(std::ostream&) const;
};

} // End namespace Uintah

#endif //__VAANGO_CORE_DATATYPE_DENSE_COL_MAJOR_MATRIX_H__ 1
