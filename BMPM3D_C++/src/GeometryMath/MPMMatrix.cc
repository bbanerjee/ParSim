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

#include <MPMMatrix.h>

using namespace BrMPM;

template<class T, int rows, int cols>
MPMMatrix<T, rows, cols>::MPMMatrix()
    : d_num_rows(rows), d_num_columns(cols)
{
}

template<class T, int rows, int cols>
MPMMatrix<T, rows, cols>::MPMMatrix(const T& initialValue)
 : d_num_rows(rows), d_num_columns(cols)
{
   d_data.resize(d_num_rows*d_num_columns, initialValue);
}


template<class T, int rows, int cols>
const T& MPMMatrix<T, rows, cols>::get(int row, int column) const
{
  return d_data[row][column];
}


template<class T, int rows, int cols>
T& MPMMatrix<T, rows, cols>::get(int row, int column)
{
  return d_data[row][column];
}

template<class T, int rows, int cols>
void MPMMatrix<T, rows, cols>::set(int row, int column, const T& value)
{
  d_data[row][column] = value;
}













