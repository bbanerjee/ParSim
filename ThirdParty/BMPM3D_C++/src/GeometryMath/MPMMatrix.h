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

#ifndef __MPMMATRIX_H__
#define __MPMMATRIX_H__

#include <vector>
#include <utility>

namespace BrMPM {

  template<class T, int rows, int cols>
  class MPMMatrix
  {
    public:
      MPMMatrix();
      MPMMatrix(const T& initialValue);
          
      const T& get(int row, int column) const;
      T& get(int row, int column);
          
      void set(int row, int column, const T& value);
          
      inline int numRows() const {return d_num_rows;}
      inline int numColumns() const {return d_num_columns;}

    private:
      int d_num_rows;
      int d_num_columns;
      std::vector<std::vector<T> > d_data;

  }; // end class

  static const int DIM = 3;  //dimension
  static const int SHAPESIZE = 27;  // Max allowed search area (GIMP)

  typedef MPMMatrix<double, 1, DIM>  MatrixVec;
  typedef std::vector<MatrixVec> ArrayMatrixVec;
  typedef MPMMatrix<double, DIM, DIM> Matrix;
  typedef std::vector<Matrix> ArrayMatrix;

  typedef MPMMatrix<int, 1, SHAPESIZE>  IntMatrixVecShape;
  typedef std::vector<IntMatrixVecShape> ArrayIntMatrixVecShape;

  typedef MPMMatrix<double, 1, SHAPESIZE>  MatrixVecShape;
  typedef std::vector<MatrixVecShape> ArrayMatrixVecShape;
  typedef MPMMatrix<double, SHAPESIZE, DIM> MatrixShape;
  typedef std::vector<MatrixShape> ArrayMatrixShape;

  typedef std::pair<int, ArrayMatrixVec> MatArrayMatrixVec;
  typedef std::pair<int, ArrayMatrix> MatArrayMatrix;
  typedef std::pair<int, ArrayIntMatrixVecShape> MatArrayIntMatrixVecShape;
  typedef std::pair<int, ArrayMatrixVecShape> MatArrayMatrixVecShape;
  typedef std::pair<int, ArrayMatrixShape> MatArrayMatrixShape;

} // end namespace
              





































#endif
