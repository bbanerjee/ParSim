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

#ifndef __SYMM_MATRIX6_H__
#define __SYMM_MATRIX6_H__

#include <Eigen/Dense>
#include <Core/Math/Matrix3.h>

#include <vector>

namespace Uintah {

  using Vector6d = Eigen::Matrix<double, 6, 1>;
  using Vector21d = Eigen::Matrix<double, 21, 1>;
  using Matrix6d = Eigen::Matrix<double, 6, 6>;

  /*! \class SymmMatrix6
    \brief Special matrix operations for a symmetric 6x6 matrix.  
           This is a wrapper for Eigen::MatrixXd.
   
  */

  class SymmMatrix6 {

  public:
    // constructors
    inline SymmMatrix6();
    inline SymmMatrix6(const SymmMatrix6& mat);
    inline SymmMatrix6(const Matrix6d& mat);

    /*! Size of matrix = 6x6; need to store n(n+1)/2 = 6*7/2 = 21 elements 
        The std:vector of values stores the sequence 
          M_{11}, M_{12}, M_{13}, M_{14}, M_{15}, M_{16}
          M_{22}, M_{23}, M_{24}, M_{25}, M_{26},
          M_{33}, M_{34}, M_{35}, M_{36},
          M_{44}, M_{45}, M_{46},
          M_{55}, M_{56}, 
          M_{66} 
    */
    inline SymmMatrix6(const std::vector<double>& vals);
    inline SymmMatrix6(const Vector21d& vec);

    // destructor
    virtual inline ~SymmMatrix6();

    /*! Access operators */
    inline void operator=(const SymmMatrix6& mat);
    inline void operator=(const Matrix6d& mat);
    inline void operator=(const Vector21d& vec);
    inline double operator() (int i, int j) const;
    inline double& operator() (int i, int j);

    /*! Matrix-vector product */
    Matrix3 operator*(const Matrix3& mat) const;

    /*! Compute Identity matrix */
    inline void Identity();

    /*! Compute trace of the matrix */
    inline double Trace() const;

    /*! Compute the inverse of a SymmMatrix6 */
    void inverse(SymmMatrix6& inv);

    /*! Calculate eigenvalues */
    void eigen(Vector6d& eval, Matrix6d& evec);

  private:

    const static int d_size = 6;
    Matrix6d d_mat6;

  };

  inline SymmMatrix6::SymmMatrix6()
  {
    // Initialization to 0.0
    for(int ii = 0; ii < d_size; ii++) {
      for(int jj = 0; jj < d_size; jj++) {
        d_mat6(ii, jj) = 0.0;
      }
    }
  }

  inline SymmMatrix6::SymmMatrix6(const SymmMatrix6& copy)
  {
    d_mat6 = copy.d_mat6;
  }

  inline SymmMatrix6::SymmMatrix6(const Matrix6d& mat)
  {
    d_mat6 = mat; 
  }

  inline SymmMatrix6::SymmMatrix6(const std::vector<double>& vals)
  {
    int count = 0;
    for (int ii = 0; ii < d_size; ++ii) {
      for (int jj = ii; jj < d_size; jj++) {
        d_mat6(jj,ii) = vals[count++];
        if (ii != jj) {
          d_mat6(ii,jj) = d_mat6(jj,ii);
        }
      }
    }
  }

  inline SymmMatrix6::SymmMatrix6(const Vector21d& vec)
  {
    int count = 0;
    for (int ii = 0; ii < d_size; ++ii) {
      for (int jj = ii; jj < d_size; jj++) {
        d_mat6(jj,ii) = vec(count++);
        if (ii != jj) {
          d_mat6(ii,jj) = d_mat6(jj,ii);
        }
      }
    }
  }

  inline SymmMatrix6::~SymmMatrix6()
  {
  } 

  inline void SymmMatrix6::operator=(const SymmMatrix6& mat)
  {
    d_mat6 = mat.d_mat6;
  }

  inline void SymmMatrix6::operator=(const Matrix6d& mat)
  {
    d_mat6 = mat;
  }

  inline void SymmMatrix6::operator=(const Vector21d& vec)
  {
    int count = 0;
    for (int ii = 0; ii < d_size; ++ii) {
      for (int jj = ii; jj < d_size; jj++) {
        d_mat6(jj,ii) = vec(count++);
        if (ii != jj) {
          d_mat6(ii,jj) = d_mat6(jj,ii);
        }
      }
    }
  }

  inline double SymmMatrix6::operator()(int i, int j) const
  {
    return d_mat6(i,j); 
  }

  inline double& SymmMatrix6::operator()(int i, int j)
  {
    return d_mat6(i,j); 
  }

  inline void SymmMatrix6::Identity()
  {
    for (int ii = 0; ii < d_size; ++ii) {
      for (int jj = ii; jj < d_size; jj++) {
        if (ii == jj) {
          d_mat6(ii,jj) = 1.0;
        } else {
          d_mat6(ii,jj) = 0.0;
        }
      }
    }
  }

  inline double SymmMatrix6::Trace() const
  {
    double trace = 0.0;
    for (int ii = 0; ii < d_size; ++ii) {
      trace += d_mat6(ii, ii);
    }
    return trace;
  }

} // End namespace Uintah

#endif  // __SYMM_MATRIX6_H__

