/*
 * Matrix3D.h
 *
 *  Created on: 16/10/2013
 *      Author: banerjee
 */

#ifndef MATRIX3D_H_
#define MATRIX3D_H_

#include <Geometry/Vector3D.h>

namespace BrMPM {

  class Matrix3D
  {
  public:
    inline Matrix3D();
    inline Matrix3D(double value);
    inline Matrix3D(double v00, double v01, double v02,
        double v10, double v11, double v12,
        double v20, double v21, double v22);
    inline Matrix3D(const Matrix3D&);
    inline Matrix3D(const Vector3D& v1, const Vector3D& v2);

    inline ~Matrix3D();

    inline void set(double val);
    void set(int i, int j, double val);
    inline void operator = (const Matrix3D &m3);
    inline double operator() (int i,int j) const;
    inline double & operator() (int i,int j);
    inline Matrix3D operator * (const double value) const;
    inline Matrix3D operator / (const double value) const;
    inline void operator += (const Matrix3D &m3);
    inline void operator -= (const Matrix3D &m3);
    inline Matrix3D operator + (const Matrix3D &m3) const;
    inline Matrix3D operator * (const Matrix3D &m3) const;
    inline Vector3D operator * (const Vector3D& V) const;
    inline Matrix3D operator - (const Matrix3D &m3) const;
    inline bool operator==(const Matrix3D &m3) const;
    inline bool operator!=(const Matrix3D &m3) const
              { return !(*this == m3); }
    inline void operator *= (const double value);
    inline void operator /= (const double value);
    inline double Determinant() const;
    inline bool Orthogonal() const;
    inline bool solveCramer(Vector3D& rhs, Vector3D& x) const;
    Matrix3D Inverse() const;
    inline double Trace() const;
    inline double Norm() const;
    inline double NormSquared() const;
    inline double Contract(const Matrix3D& mat) const;
    inline double MaxAbsElem() const;
    inline void Identity();
    inline Matrix3D Transpose() const;
    Matrix3D Exponential(int num_terms) const;
    Matrix3D Logarithm(int num_terms) const;

  private:
    double mat3[3][3];
  };

  inline Matrix3D operator*(double c, const Matrix3D &m3);
  inline Vector3D operator*(const Vector3D& v, const Matrix3D& m3);

} /* namespace BrMPM */
#endif /* MATRIX3D_H_ */
