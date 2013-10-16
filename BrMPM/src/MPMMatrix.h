#ifndef __MPMMATRIX_H__
#define __MPMMATRIX_H__

#include <vector>

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

  typedef MPMMatrix<int, 27, 1>  NodeIndexVector;
  typedef MPMMatrix<double, 27, 1>  NodeWeightVector;
  typedef MPMMatrix<double, 27, 3> ShapeGradientMatrix;

} // end namespace
              





































#endif
