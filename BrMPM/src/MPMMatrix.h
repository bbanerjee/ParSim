#ifndef __MPMMATRIX_H__
#define __MPMMATRIX_H__

#include <vector>

namespace MPM {

  template<class T, int numRows, int numColumns>
  class MPMMatrix
  {
    public:
      MPMMatrix();
      //    MPMMatrix(int numRows, int numColumns);
      MPMMatrix(const T& initialValue);
      //    MPMMatrix(int numRows, int numColumns, const T& initialValue);
          
      const T& get(int row, int column) const;
      T& get(int row, int column);
          
      void set(int row, int column, const T& value);
          
      inline int numRows() const {return d_num_rows;}
      inline int numColumns() const {return d_num_columns;}

    private:
      int d_num_rows;
      int d_num_columns;
      std::vector<T> d_vector;
  };
}
              





































#endif
