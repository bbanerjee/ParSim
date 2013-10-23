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













