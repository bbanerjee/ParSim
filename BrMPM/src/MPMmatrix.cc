#include <MPMmatrix.h>


namespace MPM;

//template<class T>
template<class T, int numRows, int numColumns>
MPMmatrix<T>::MPMmatrix()
// : d_num_rows(1), d_num_columns(1);
 : d_num_rows(numRows), d_num_columns(numColumns);
{
}


/*
//template<class T>
template<class T, int numRows, int numColumns>
MPMmatrix<T>::MPMmatrix(int numRows, int numColumns)
 : d_num_rows(numRows), d_num_columns(numColumns);
{
   d_vector.resize(numRows*numColumns);
}  */


//template<class T>
template<class T, int numRows, int numColumns>
//MPMmatrix<T>::MPMmatrix(int numRows, int numColumns, const T& initialValue)
MPMmatrix<T>::MPMmatrix(const T& initialValue)
 : d_num_rows(numRows), d_num_columns(numColumns);
{
   d_vector.resize(numRows*numColumns, initialValue);
}


//template<class T>
template<class T, int numRows, int numColumns>
 const T& MPMmatrix<T>::get(int row, int column) const
{
     return d_vector[row*d_num_columns+column);
}


//template<class T>
template<class T, int numRows, int numColumns>
 T& MPMmatrix<T>::get(int row, int column)
{
     return d_vector[row*d_num_columns+column);
}



//template<class T>
template<class T, int numRows, int numColumns>
void MPMmatrix<T>::set(int row, int column, const T& value)
{
     d_vector[row*d_num_columns+column) = value;
}
















