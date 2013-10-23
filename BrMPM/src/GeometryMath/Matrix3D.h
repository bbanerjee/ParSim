/*
 * Matrix3D.h
 *
 *  Created on: 16/10/2013
 *      Author: banerjee
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include <GeometryMath/Vector3D.h>

namespace BrMPM {

  class Matrix3D
  {
  public:
    Matrix3D();
    Matrix3D(double value);
    Matrix3D(double v00, double v01, double v02,
        double v10, double v11, double v12,
        double v20, double v21, double v22);
    Matrix3D(const Matrix3D&);
    Matrix3D(const Vector3D& v1, const Vector3D& v2);

    ~Matrix3D();

    void set(double val);
    void set(int i, int j, double val);
    void operator = (const Matrix3D &m3);
    double operator() (int i,int j) const;
    double & operator() (int i,int j);
    Matrix3D operator * (const double value) const;
    Matrix3D operator / (const double value) const;
    void operator += (const Matrix3D &m3);
    void operator -= (const Matrix3D &m3);
    Matrix3D operator + (const Matrix3D &m3) const;
    Matrix3D operator * (const Matrix3D &m3) const;
    Vector3D operator * (const Vector3D& V) const;
    Matrix3D operator - (const Matrix3D &m3) const;
    bool operator==(const Matrix3D &m3) const;
    bool operator!=(const Matrix3D &m3) const
              { return !(*this == m3); }
    void operator *= (const double value);
    void operator /= (const double value);
    double Determinant() const;
    bool Orthogonal() const;
    bool solveCramer(Vector3D& rhs, Vector3D& x) const;
    Matrix3D Inverse() const;
    double Trace() const;
    double Norm() const;
    double NormSquared() const;
    double Contract(const Matrix3D& mat) const;
    double MaxAbsElem() const;
    void Identity();
    Matrix3D Transpose() const;
    Matrix3D Exponential(int num_terms) const;
    Matrix3D Logarithm(int num_terms) const;

  private:
    double mat3[3][3];
  };

  Matrix3D operator*(double c, const Matrix3D &m3);
  Vector3D operator*(const Vector3D& v, const Matrix3D& m3);
  Matrix3D dyadicProduct(const Vector3D& v1, const Vector3D& v2);

} /* namespace BrMPM */
#endif /* MATRIX3D_H_ */
