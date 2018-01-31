/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#ifndef __ELLIP3D_CORE_MATH_MATRIX2_H__
#define __ELLIP3D_CORE_MATH_MATRIX2_H__

#include <Core/Math/Vec.h>

#include <cmath>

#include <iosfwd>
#include <vector>
#include <string>

namespace dem {

  class Matrix2 {

  private:
    double mat2[2][2];

  public:
    inline Matrix2();
    inline Matrix2(double value);
    inline Matrix2(double v00, double v01,
                   double v10, double v11);
    inline Matrix2(const Matrix2&);

    // assign a value to all components of the Matrix2
    inline void set(double val);

    // assign a value to components of the Matrix2
     void set(int i, int j, double val);

    // assignment operator
    inline void operator = (const Matrix2 &m2);

    // access operator
    inline double operator() (int i,int j) const;
    inline double & operator() (int i,int j);

    // multiply Matrix2 by a constant: Mat * 3  
    inline Matrix2 operator * (const double value) const;
    // multiply constant by Matrix2:   3 * Mat
     friend Matrix2 operator * (double c, const Matrix2 &m2);

    // divide Matrix2 by a constant
    inline Matrix2 operator / (const double value) const;
 
    // modify by adding right hand side
    inline void operator += (const Matrix2 &m2);

    // modify by subtracting right hand side
    inline void operator -= (const Matrix2 &m2);

    // add two Matrix2s
    inline Matrix2 operator + (const Matrix2 &m2) const;

    // multiply two Matrix2s
    inline Matrix2 operator * (const Matrix2 &m2) const;

    // subtract two Matrix2s
    inline Matrix2 operator - (const Matrix2 &m2) const;

    // compare two Matrix2s
    inline bool operator==(const Matrix2 &m2) const;
    inline bool operator!=(const Matrix2 &m2) const
      { return !(*this == m2); }

    // multiply by constant
    inline void operator *= (const double value);

    // divide by constant
    inline void operator /= (const double value);

    //Determinant
    inline double Determinant() const;

    //Trace
    inline double Trace() const;

    //Norm, sqrt(M:M)
    inline double Norm() const;

    //NormSquared, M:M
    inline double NormSquared() const;

    //Contraction of two tensors
    inline double Contract(const Matrix2& mat) const;

    //Maximum element absolute value
    inline double MaxAbsElem() const;
  
    //Identity
    inline void Identity();

    //Transpose
    inline Matrix2 Transpose() const;

    // Output in nice format
    void prettyPrint(std::ostream &out_file) const;


  private:

     friend std::ostream & operator << (std::ostream &out_file, const Matrix2 &m2);

  }; // end class Matrix2

  inline double Matrix2::Trace() const {
      double trace = 0.0;
      for (int i = 0; i < 2; i++) {
        trace += mat2[i][i];
      }
      return trace;
  }

  inline Matrix2::Matrix2()
    {
      // Default Constructor
      // Initialization to 0.0
      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] = 0.0;
        }
      }
    }

  inline Matrix2::Matrix2(double value)
    {
      // Constructor
      // With initialization to a single value

      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] = value;
        }
      }
    }


  inline Matrix2::Matrix2(double a00,double a01,
                          double a10,double a11)
    {
      // Constructor
      // With full initialization

      mat2[0][0] = a00; mat2[0][1] = a01; 
      mat2[1][0] = a10; mat2[1][1] = a11;

    }

  inline Matrix2::Matrix2(const Matrix2& copy)
    {
      // Copy Constructor

      mat2[0][0] = copy.mat2[0][0]; mat2[0][1] = copy.mat2[0][1];
      mat2[1][0] = copy.mat2[1][0]; mat2[1][1] = copy.mat2[1][1];

    }

  inline double Matrix2::NormSquared() const
    {
      // Return the norm of a 2x2 matrix

      double norm = 0.0;

      for (int i = 0; i< 2; i++) {
        for(int j=0;j<2;j++){
          norm += mat2[i][j]*mat2[i][j];
        }
      }
      return norm;
    }

  inline double Matrix2::Norm() const
    {
      return sqrt(NormSquared());
    }

  inline double Matrix2::Contract(const Matrix2& mat) const
    {
      // Return the contraction of this matrix with another 

      double contract = 0.0;

      for (int i = 0; i< 2; i++) {
        for(int j=0;j<2;j++){
          contract += mat2[i][j]*(mat.mat2[i][j]);
        }
      }
      return contract;
    }

  inline double Matrix2::MaxAbsElem() const
    {
      double max = 0;
      double absval;
      for (int i = 0; i< 2; i++) {
        for(int j=0;j<2;j++){
          absval = fabs(mat2[i][j]);
          if (absval > max) max = absval;
        }
      }
      return max;
    }

  inline Matrix2 Matrix2::Transpose() const
    {
      // Return the transpose of a 2x2 matrix

      return Matrix2(mat2[0][0],mat2[1][0],
                     mat2[0][1],mat2[1][1]);

    }

  inline void Matrix2::Identity()
    {
      // Set a matrix2 to the identity

      mat2[0][0] = mat2[1][1] = 1.0;
      mat2[0][1] = mat2[1][0] = 0.0;

    }

  inline void Matrix2::operator = (const Matrix2 &m2)
    {
      // Copy value from right hand side of assignment

      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] = m2(i,j);
        }
      }

    }

  inline void Matrix2::set(double value)
    {
      // Assign the Matrix2 the value components
      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] = value;
        }
      }
    }

  inline void Matrix2::operator *= (const double value)
    {
      // Multiply each component of the Matrix2 by the value

      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] *= value;
        }
      }
    }

  inline void Matrix2::operator += (const Matrix2 &m2)
    {
      // += operator 

      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] += m2(i,j);
        }
      }

    }

  inline void Matrix2::operator -= (const Matrix2 &m2)
    {
      // -= operator 

      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] -= m2(i,j);
        }
      }

    }

  inline void Matrix2::operator /= (const double value)
    {
      // Divide each component of the Matrix2 by the value

      assert(value != 0.);
      double ivalue = 1./value;
      for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
          mat2[i][j] *= ivalue;
        }
      }
    }

  inline Matrix2 Matrix2::operator * (const double value) const
    {
      //   Multiply a Matrix2 by a constant

      return Matrix2(mat2[0][0]*value,mat2[0][1]*value,
                     mat2[1][0]*value,mat2[1][1]*value);
    }

  inline Matrix2 operator * ( double c, const Matrix2 &m2 )
    {
      return m2 * c;
    }

  inline Matrix2 Matrix2::operator * (const Matrix2 &m2) const
    {
      //   Multiply a Matrix2 by a Matrix2

      return Matrix2(mat2[0][0]*m2(0,0)+mat2[0][1]*m2(1,0),
                     mat2[0][0]*m2(0,1)+mat2[0][1]*m2(1,1),

                     mat2[1][0]*m2(0,0)+mat2[1][1]*m2(1,0),
                     mat2[1][0]*m2(0,1)+mat2[1][1]*m2(1,1));

    }


  inline Matrix2 Matrix2::operator + (const Matrix2 &m2) const
    {
      //   Add a Matrix2 to a Matrix2

      return Matrix2(mat2[0][0] + m2(0,0),mat2[0][1] + m2(0,1),
                     mat2[1][0] + m2(1,0),mat2[1][1] + m2(1,1));
    }

  inline Matrix2 Matrix2::operator - (const Matrix2 &m2) const
    {
      //   Subtract a Matrix2 from a Matrix2

      return Matrix2(mat2[0][0] - m2(0,0),mat2[0][1] - m2(0,1),
                     mat2[1][0] - m2(1,0),mat2[1][1] - m2(1,1));
    }

  inline bool Matrix2::operator==(const Matrix2 &m2) const
    {
      return
        mat2[0][0] == m2(0,0) && mat2[0][1] == m2(0,1) && 
        mat2[1][0] == m2(1,0) && mat2[1][1] == m2(1,1) ;
    }

  inline Matrix2 Matrix2::operator / (const double value) const
    {
      //   Divide a Matrix2 by a constant

      assert(value != 0.);
      double ivalue = 1.0/value;

      return Matrix2(mat2[0][0]*ivalue,mat2[0][1]*ivalue,
                     mat2[1][0]*ivalue,mat2[1][1]*ivalue);
    }

  inline double Matrix2::Determinant() const
    {
      // Return the determinant of a 2x2 matrix

      double temp = 0.0;

      temp= mat2[0][0]*mat2[1][1] - mat2[0][1]*mat2[1][0];

      // return result

      return temp;
    }

  inline double Matrix2::operator () (int i, int j) const
    {
      // Access the i,j component
      return mat2[i][j];
    }

  inline double &Matrix2::operator () (int i, int j)
    {
      // Access the i,j component
      return mat2[i][j];
    }

} // End namespace dem


#endif  // __ELLIP3D_CORE_MATH_MATRIX2_H__

