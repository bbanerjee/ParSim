
#include <Core/Types/realtypes.h>
#include <boost/mpi.hpp>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

#ifndef MATRIX_H
#define MATRIX_H

// three ways to initialize/ assign value to a Matrix object:
// 1, Matrix A; A.appendRow(std::vector<REAL>);	this is not a convenient way,
// 2, Matrix A(int i, int j); (=> Aij = 0) then we can use A(i,j)=REAL to assign
// value
// each element of A
// 3, Matrix A; A = ["1,2,2;2,3,4"] ([] is overloaded which will return a
// Matrix), this is plan

namespace dem {

class Matrix
{
public:
  REAL* value; // 2-dimension, values are stored by row first, namely
               // Matrix(i,j) is value[i*num_col+j]
  int num_row; // number of rows
  int num_col; // number of columns

  Matrix()
    : num_row(0)
    , num_col(0)
  {
    value = NULL;
  }                 // constructor
  Matrix(int, int); // constructor, way 2
  Matrix(const Matrix&);
  ~Matrix();

  //	void calcDimensions();
  //	REAL* getCol(int i) const;	// get ith colume, 1<= i <=num_col
  //	REAL* getRow(int i) const;	// get ith row
  //	void appendRow(REAL*);
  //	void appendCol(REAL*);
  void
  clear(); // clear all value of the Matrix, assign num_row and num_col to be 0
           //	void getInvs();	// return the inverse
           //	void getTrans();	// return the transpose
  REAL getNorm(); // return the norm of Matrix, at present it can only calculate
                  // norm of a vector
                  // void add(REAL);	// row vector
  void LU(Matrix& L, Matrix& U) const; // LU decomposition, November 24
  REAL getVal(int, int)
    const; // access to the element at ith row and jth column: 1<= i <= num_rol

  Matrix& operator+=(const Matrix& A)
  {
    (*this) = (*this) + A;
    return *this;
  }
  Matrix& operator+=(REAL k)
  {
    (*this) = (*this) + k;
    return *this;
  }
  Matrix& operator-=(const Matrix& A)
  {
    (*this) = (*this) - A;
    return *this;
  }
  Matrix& operator-=(REAL k)
  {
    (*this) = (*this) - k;
    return *this;
  }
  Matrix& operator*=(const Matrix& A)
  {
    (*this) = (*this) * A;
    return *this;
  }
  Matrix& operator*=(REAL k)
  {
    (*this) = (*this) * k;
    return *this;
  }
  Matrix& operator=(const Matrix& A); // overload = for Matrix

  // access to the element at ith row and jth column: 1<= i <= num_rol
  REAL& operator()( int i, int j); 
  REAL operator()( int i, int j) const; 
  std::string print() const; // used for debug

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  { // need to think about how to serialize the pointer array
    //      ar & value;
    ar& num_row;
    ar& num_col;
  }

  // non-member functions
  friend Matrix operator+(const Matrix& A,
                          const Matrix& B);     // overload + for two Matrix
  friend Matrix operator+(const Matrix&, REAL); // overload + for two Matrix
  friend Matrix operator+(REAL,
                          const Matrix&); // overload + for scalar and Matrix
  friend Matrix operator-(const Matrix& A,
                          const Matrix& B); // overload - for Matrix and scalar
  friend Matrix operator-(REAL,
                          const Matrix&); // overload - for scalar and Matrix
  friend Matrix operator-(const Matrix&,
                          REAL); // overload - for Matrix and scalar
  friend Matrix operator*(const Matrix& A,
                          const Matrix& B); // overload * for two Matrix
  friend Matrix operator*(const Matrix& A,
                          REAL k); // overload * for Matrix and scalar
  friend Matrix operator*(
    REAL k, const Matrix& A); // overload * for a scalar and a Matrix
  // friend Matrix operator / (Matrix A, Matrix B);
  friend Matrix operator/(const Matrix& A,
                          REAL k); // overload / for Matrix and a scalar
  friend Matrix operator%(Matrix& A,
                          Matrix& B); // November 15, leftdivision, "\"

  // friend Matrix operator [] (stringstream);	// overload [] for assignment
  friend Matrix expm(const Matrix&); // calculate the exponential of each
                                     // element in the Matrix. Not used here

  friend std::ostream& operator<<(std::ostream& os, const Matrix& mat);
}; // end Matrix

// October 31, 2013

bool isnan_mat(Matrix&);
int size(const Matrix&, int); // the same size function in matlab
Matrix ones(int, int);        // the same ones function in matlab
Matrix zeros(int, int);       // the same zeros function in matlab
REAL max(Matrix&);            // the same max function in matlab
REAL min(Matrix&);            // the same min function in matlab
// Matrix abs(Matrix &);	// the same abs function in matlab
int length(const Matrix&); // the same length function in matlab
REAL norm(Matrix&);        // the same norm function in matlab
REAL det(Matrix&);         // the same det function in matlab
REAL trace(Matrix&);
Matrix linspace(REAL, REAL, int num);
Matrix inv(const Matrix&);
Matrix trans(Matrix&);

void MatrixEqnSolver(Matrix& dsol, Matrix& K, Matrix& F);

} // end of namespace dem

#endif
