/*
 * Matrix3D.cc
 *
 *  Created on: 16/10/2013
 *      Author: banerjee
 */

#include <GeometryMath/Matrix3D.h>

using namespace BrMPM;

// Default Constructor
// Initialization to 0.0
Matrix3D::Matrix3D()
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] = 0.0;
    }
  }
}

// Constructor
// With initialization to a single value
Matrix3D::Matrix3D(double value)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] = value;
    }
  }
}

// Constructor
// With full initialization
Matrix3D::Matrix3D(double a00,double a01,double a02,
                          double a10,double a11,double a12,
                          double a20,double a21,double a22)
{
  mat3[0][0] = a00; mat3[0][1] = a01; mat3[0][2] = a02;
  mat3[1][0] = a10; mat3[1][1] = a11; mat3[1][2] = a12;
  mat3[2][0] = a20; mat3[2][1] = a21; mat3[2][2] = a22;
}

// Copy Constructor
Matrix3D::Matrix3D(const Matrix3D& copy)
{
  mat3[0][0] = copy.mat3[0][0]; mat3[0][1] = copy.mat3[0][1]; mat3[0][2] = copy.mat3[0][2];
  mat3[1][0] = copy.mat3[1][0]; mat3[1][1] = copy.mat3[1][1]; mat3[1][2] = copy.mat3[1][2];
  mat3[2][0] = copy.mat3[2][0]; mat3[2][1] = copy.mat3[2][1]; mat3[2][2] = copy.mat3[2][2];
}

// Create a matrix using dyadic multiplication M = a*b^T
// Dyadic multiplication of two vectors a and b
Matrix3D::Matrix3D(const Vector3D& a, const Vector3D& b)
{
  mat3[0][0] = a[0]*b[0]; mat3[0][1] = a[0]*b[1]; mat3[0][2] = a[0]*b[2];
  mat3[1][0] = a[1]*b[0]; mat3[1][1] = a[1]*b[1]; mat3[1][2] = a[1]*b[2];
  mat3[2][0] = a[2]*b[0]; mat3[2][1] = a[2]*b[1]; mat3[2][2] = a[2]*b[2];
}

Matrix3D::~Matrix3D()
{
}

// Assign the Matrix3 the value components
void Matrix3D::set(int i, int j, double value)
{
  mat3[i][j] = value;
}

// Return the inverse of a 3x3 matrix
// This looks ugly but it works -- just for 3x3
Matrix3D Matrix3D::Inverse() const
{
  double det;
  Matrix3D inv_matrix(0.0);

  det = this->Determinant();
  if ( det == 0.0 )
  {
    std::cerr << "Singular matrix in matrix inverse..." << std::endl;
    exit(1);
  }
  else
  {
    inv_matrix(0,0) = (*this)(1,1)*(*this)(2,2) - (*this)(1,2)*(*this)(2,1);
    inv_matrix(0,1) = -(*this)(0,1)*(*this)(2,2) + (*this)(2,1)*(*this)(0,2);
    inv_matrix(0,2) = (*this)(0,1)*(*this)(1,2) - (*this)(1,1)*(*this)(0,2);
    inv_matrix(1,0) = -(*this)(1,0)*(*this)(2,2) + (*this)(2,0)*(*this)(1,2);
    inv_matrix(1,1) = (*this)(0,0)*(*this)(2,2) - (*this)(0,2)*(*this)(2,0);
    inv_matrix(1,2) = -(*this)(0,0)*(*this)(1,2) + (*this)(1,0)*(*this)(0,2);
    inv_matrix(2,0) = (*this)(1,0)*(*this)(2,1) - (*this)(2,0)*(*this)(1,1);
    inv_matrix(2,1) = -(*this)(0,0)*(*this)(2,1) + (*this)(2,0)*(*this)(0,1);
    inv_matrix(2,2) = (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);

    inv_matrix = inv_matrix/det;

  }
  return inv_matrix;

} //end Inverse()

// A recursive Taylor series expansion (USE WITH CARE)
// **WARNING** Expansion may not be convergent in which case use
// eigenvalue expansion (not implemented)
// Based on Ortiz, Radovitzsky, Repetto (2001)
Matrix3D Matrix3D::Exponential(int num_terms) const
{
  Matrix3D exp(0.0);

  // The k = 0 term
  Matrix3D exp_k_term; exp_k_term.Identity();
  exp += exp_k_term;

  for (int kk = 0; kk < num_terms; ++kk) {
    exp_k_term = (exp_k_term*(*this))*(1.0/(double)(kk+1));
    exp += exp_k_term;
  }
  return exp;
}

// A recursive Taylor series expansion (USE WITH CARE)
// **WARNING** Expansion may not be convergent in which case use
// eigenvalue expansion (not implemented)
// Based on Ortiz, Radovitzsky, Repetto (2001)
Matrix3D Matrix3D::Logarithm(int num_terms) const
{
  Matrix3D log(0.0);
  Matrix3D One; One.Identity();

  // The k = 0 term
  Matrix3D log_0_term(0.0), log_k_term(0.0);
  log_0_term = *this - One;
  log_k_term = log_0_term;
  log += log_k_term;

  for (int ii = 1; ii <= num_terms; ++ii) {
    log_k_term = (log_k_term*log_0_term)*((double)ii/(double)(ii+1));
    log += log_k_term;
  }
  return log;
}

// Return the trace of a 3x3 matrix
double Matrix3D::Trace() const
{
  double trace = 0.0;
  for (int i = 0; i< 3; i++) {
    trace += mat3[i][i];
  }
  return trace;
}

// Return the norm of a 3x3 matrix
double Matrix3D::NormSquared() const
{
  double norm = 0.0;
  for (int i = 0; i< 3; i++) {
    for(int j=0;j<3;j++){
      norm += mat3[i][j]*mat3[i][j];
    }
  }
  return norm;
}

double Matrix3D::Norm() const
{
  return sqrt(NormSquared());
}

// Return the contraction of this matrix with another
double Matrix3D::Contract(const Matrix3D& mat) const
{
  double contract = 0.0;
  for (int i = 0; i< 3; i++) {
    for(int j=0;j<3;j++){
      contract += mat3[i][j]*(mat.mat3[i][j]);
    }
  }
  return contract;
}

double Matrix3D::MaxAbsElem() const
{
  double max = 0;
  double absval;
  for (int i = 0; i< 3; i++) {
    for(int j=0;j<3;j++){
      absval = fabs(mat3[i][j]);
      if (absval > max) max = absval;
    }
  }
  return max;
}

// Return the transpose of a 3x3 matrix
Matrix3D Matrix3D::Transpose() const
{
  return Matrix3D(mat3[0][0],mat3[1][0],mat3[2][0],
      mat3[0][1],mat3[1][1],mat3[2][1],
      mat3[0][2],mat3[1][2],mat3[2][2]);
}

// Set a matrix3 to the identity
void Matrix3D::Identity()
{
  mat3[0][0] = mat3[1][1] = mat3[2][2] = 1.0;
  mat3[0][1] = mat3[0][2] = mat3[1][0] = 0.0;
  mat3[1][2] = mat3[2][0] = mat3[2][1] = 0.0;
}

// Copy value from right hand side of assignment
void Matrix3D::operator=(const Matrix3D &m3)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] = m3(i,j);
    }
  }
}

// Assign the Matrix3D the value components
void Matrix3D::set(double value)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] = value;
    }
  }
}

// Multiply each component of the Matrix3D by the value
void Matrix3D::operator*=(const double value)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] *= value;
    }
  }
}

// += operator
void Matrix3D::operator+=(const Matrix3D &m3)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] += m3(i,j);
    }
  }
}

// -= operator
void Matrix3D::operator-=(const Matrix3D &m3)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] -= m3(i,j);
    }
  }
}

// Divide each component of the Matrix3D by the value
void Matrix3D::operator/=(const double value)
{
  double ivalue = 1./value;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      mat3[i][j] *= ivalue;
    }
  }
}

//   Multiply a Matrix3D by a constant
Matrix3D Matrix3D::operator*(const double value) const
{
  return Matrix3D(mat3[0][0]*value,mat3[0][1]*value,mat3[0][2]*value,
      mat3[1][0]*value,mat3[1][1]*value,mat3[1][2]*value,
      mat3[2][0]*value,mat3[2][1]*value,mat3[2][2]*value);
}

//   Multiply a Matrix3D by a Matrix3
Matrix3D Matrix3D::operator*(const Matrix3D &m3) const
{
  return Matrix3D(mat3[0][0]*m3(0,0)+mat3[0][1]*m3(1,0)+mat3[0][2]*m3(2,0),
      mat3[0][0]*m3(0,1)+mat3[0][1]*m3(1,1)+mat3[0][2]*m3(2,1),
      mat3[0][0]*m3(0,2)+mat3[0][1]*m3(1,2)+mat3[0][2]*m3(2,2),

      mat3[1][0]*m3(0,0)+mat3[1][1]*m3(1,0)+mat3[1][2]*m3(2,0),
      mat3[1][0]*m3(0,1)+mat3[1][1]*m3(1,1)+mat3[1][2]*m3(2,1),
      mat3[1][0]*m3(0,2)+mat3[1][1]*m3(1,2)+mat3[1][2]*m3(2,2),

      mat3[2][0]*m3(0,0)+mat3[2][1]*m3(1,0)+mat3[2][2]*m3(2,0),
      mat3[2][0]*m3(0,1)+mat3[2][1]*m3(1,1)+mat3[2][2]*m3(2,1),
      mat3[2][0]*m3(0,2)+mat3[2][1]*m3(1,2)+mat3[2][2]*m3(2,2));
}

Vector3D Matrix3D::operator*(const Vector3D& V) const
{
  return Vector3D(mat3[0][0] * V[0] + mat3[0][1] * V[1] + mat3[0][2] * V[2],
      mat3[1][0] * V[0] + mat3[1][1] * V[1] + mat3[1][2] * V[2],
      mat3[2][0] * V[0] + mat3[2][1] * V[1] + mat3[2][2] * V[2]);
}

//   Add a Matrix3D to a Matrix3
Matrix3D Matrix3D::operator+(const Matrix3D &m3) const
{
  return Matrix3D(mat3[0][0] + m3(0,0),mat3[0][1] + m3(0,1),mat3[0][2] + m3(0,2),
      mat3[1][0] + m3(1,0),mat3[1][1] + m3(1,1),mat3[1][2] + m3(1,2),
      mat3[2][0] + m3(2,0),mat3[2][1] + m3(2,1),mat3[2][2] + m3(2,2));
}

//   Subtract a Matrix3D from a Matrix3
Matrix3D Matrix3D::operator-(const Matrix3D &m3) const
{
  return Matrix3D(mat3[0][0] - m3(0,0),mat3[0][1] - m3(0,1),mat3[0][2] - m3(0,2),
      mat3[1][0] - m3(1,0),mat3[1][1] - m3(1,1),mat3[1][2] - m3(1,2),
      mat3[2][0] - m3(2,0),mat3[2][1] - m3(2,1),mat3[2][2] - m3(2,2));
}

bool Matrix3D::operator==(const Matrix3D &m3) const
{
  return
      mat3[0][0] == m3(0,0) && mat3[0][1] == m3(0,1) && mat3[0][2] == m3(0,2) &&
      mat3[1][0] == m3(1,0) && mat3[1][1] == m3(1,1) && mat3[1][2] == m3(1,2) &&
      mat3[2][0] == m3(2,0) && mat3[2][1] == m3(2,1) && mat3[2][2] == m3(2,2);
}

//   Divide a Matrix3D by a constant
Matrix3D Matrix3D::operator/(const double value) const
{
  double ivalue = 1.0/value;
  return Matrix3D(mat3[0][0]*ivalue,mat3[0][1]*ivalue,mat3[0][2]*ivalue,
      mat3[1][0]*ivalue,mat3[1][1]*ivalue,mat3[1][2]*ivalue,
      mat3[2][0]*ivalue,mat3[2][1]*ivalue,mat3[2][2]*ivalue);
}

// Return the determinant of a 3x3 matrix
double Matrix3D::Determinant() const
{
  double temp = 0.0;
  temp= mat3[0][0]*mat3[1][1]*mat3[2][2] +
      mat3[0][1]*mat3[1][2]*mat3[2][0] +
      mat3[0][2]*mat3[1][0]*mat3[2][1] -
      mat3[0][2]*mat3[1][1]*mat3[2][0] -
      mat3[0][1]*mat3[1][0]*mat3[2][2] -
      mat3[0][0]*mat3[1][2]*mat3[2][1];
  return temp;
}

bool Matrix3D::Orthogonal() const
{
  // Identity Matrix
  Matrix3D I;  I.Identity();

  // Calculate Matrix*(Matrix.Transpose)
  Matrix3D R(*this);
  Matrix3D RRT = R*(R.Transpose());

  // Check that RRT = I
  if (RRT != I) return false;
  return true;
}

bool Matrix3D::solveCramer(Vector3D& b, Vector3D& x) const
{
  double det = this->Determinant();
  if (det == 0.0) return false;
  double odet = 1.0/det;
  x[0] = odet*(b[0]*(mat3[1][1]*mat3[2][2] - mat3[2][1]*mat3[1][2])
      - mat3[0][1]*(b[1]*mat3[2][2] - b[2]*mat3[1][2])
      + mat3[0][2]*(b[1]*mat3[2][1] - b[2]*mat3[1][1]));

  x[1] = odet*(mat3[0][0]*(b[1]*mat3[2][2] - b[2]*mat3[1][2])
      - b[0]*(mat3[1][0]*mat3[2][2] - mat3[2][0]*mat3[1][2])
      + mat3[0][2]*(mat3[1][0]*b[2]  - mat3[2][0]*b[1]));

  x[2] = odet*(mat3[0][0]*(mat3[1][1]*b[2]  - mat3[2][1]*b[1])
      - mat3[0][1]*(mat3[1][0]*b[2]  - mat3[2][0]*b[1])
      + b[0]*(mat3[1][0]*mat3[2][1] - mat3[2][0]*mat3[1][1]));
  return true;
}

// Access the i,j component
double Matrix3D::operator()(int i, int j) const
{
  return mat3[i][j];
}

// Access the i,j component
double &Matrix3D::operator()(int i, int j)
{
  return mat3[i][j];
}

namespace BrMPM {

  Matrix3D operator*( double c, const Matrix3D &m3 )
  {
    return Matrix3D(m3 * c);
  }

  // This is backwards: if the vector comes first then it should
  // multiply the matrix columnwise instead of rowwise.  For now,
  // I won't fix it because changing it may break other code. -- witzel
  Vector3D operator*(const Vector3D& v, const Matrix3D& m3)
  {
    // Right multiply a Vector3D by a Matrix3

    double x = v.x()*m3(0,0)+v.y()*m3(0,1)+v.z()*m3(0,2);
    double y = v.x()*m3(1,0)+v.y()*m3(1,1)+v.z()*m3(1,2);
    double z = v.x()*m3(2,0)+v.y()*m3(2,1)+v.z()*m3(2,2);

    return Vector3D(x, y, z);
  }

  Matrix3D dyadicProduct(const Vector3D& v1, const Vector3D& v2)
  {
    return Matrix3D(v1, v2);
  }
}

