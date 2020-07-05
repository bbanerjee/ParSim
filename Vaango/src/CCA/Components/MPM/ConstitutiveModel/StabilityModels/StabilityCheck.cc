/*
 * The MIT License
 *
 * Copyright (c) 1997-20(0,1) The University of Utah
 * Copyright (c) 20(0,2)-20(0,3) Callaghan Innovation, New Zealand
 * Copyright (c) 20(0,4)-2020 Parresia Research Limited, New Zealand
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

#include "StabilityCheck.h"

using namespace Uintah;
using namespace Vaango;

Matrix3 
StabilityCheck::elasticAcousticTensor(const Vaango::Tensor::Matrix6Mandel& C_e,
                                      const Vector& n) const
{
  // Get Cijkl values from Mandel values
  Vaango::Tensor::Matrix6Mandel C(C_e);
  C.block<3,3>(0,3) *= Vaango::Tensor::one_sqrt_two;
  C.block<3,3>(3,0) *= Vaango::Tensor::one_sqrt_two;
  C.block<3,3>(3,3) *= Vaango::Tensor::half;
  
  Matrix3 nCn;
  nCn(0,0) = n[0]*(C(0,0)*n[0] + C(0,5)*n[1] + C(0,4)*n[2]) + 
             n[1]*(C(0,5)*n[0] + C(5,5)*n[1] + C(4,5)*n[2]) + 
             n[2]*(C(0,4)*n[0] + C(4,5)*n[1] + C(4,4)*n[2]);
  nCn(0,1) = n[0]*(C(0,5)*n[0] + C(0,1)*n[1] + C(0,3)*n[2]) + 
             n[1]*(C(5,5)*n[0] + C(1,5)*n[1] + C(3,5)*n[2]) + 
             n[2]*(C(4,5)*n[0] + C(1,4)*n[1] + C(3,4)*n[2]);
  nCn(0,2) = n[0]*(C(0,4)*n[0] + C(0,3)*n[1] + C(0,2)*n[2]) + 
             n[1]*(C(4,5)*n[0] + C(3,5)*n[1] + C(2,5)*n[2]) + 
             n[2]*(C(4,4)*n[0] + C(3,4)*n[1] + C(2,4)*n[2]);
  //nCn(1,0) = n[0]*(C(0,5)*n[0] + C(5,5)*n[1] + C(4,5)*n[2]) + 
  //           n[1]*(C(0,1)*n[0] + C(1,5)*n[1] + C(1,4)*n[2]) + 
  //           n[2]*(C(0,3)*n[0] + C(3,5)*n[1] + C(3,4)*n[2]);
  nCn(1,0) = nCn(0,1);
  nCn(1,1) = n[0]*(C(5,5)*n[0] + C(1,5)*n[1] + C(3,5)*n[2]) + 
             n[1]*(C(1,5)*n[0] + C(1,1)*n[1] + C(1,3)*n[2]) + 
             n[2]*(C(3,5)*n[0] + C(1,3)*n[1] + C(3,3)*n[2]);
  nCn(1,2) = n[0]*(C(4,5)*n[0] + C(3,5)*n[1] + C(2,5)*n[2]) + 
             n[1]*(C(1,4)*n[0] + C(1,3)*n[1] + C(1,2)*n[2]) + 
             n[2]*(C(3,4)*n[0] + C(3,3)*n[1] + C(2,3)*n[2]);
  //nCn(2,0) = n[0]*(C(0,4)*n[0] + C(4,5)*n[1] + C(4,4)*n[2]) + 
  //           n[1]*(C(0,3)*n[0] + C(3,5)*n[1] + C(3,4)*n[2]) + 
  //           n[2]*(C(0,2)*n[0] + C(2,5)*n[1] + C(2,4)*n[2]);
  nCn(2,0) = nCn(0,2);
  //nCn(2,1) = n[0]*(C(4,5)*n[0] + C(1,4)*n[1] + C(3,4)*n[2]) + 
  //           n[1]*(C(3,5)*n[0] + C(1,3)*n[1] + C(3,3)*n[2]) + 
  //           n[2]*(C(2,5)*n[0] + C(1,2)*n[1] + C(2,3)*n[2]);
  nCn(2,1) = nCn(1,2);
  nCn(2,2) = n[0]*(C(4,4)*n[0] + C(3,4)*n[1] + C(2,4)*n[2]) + 
             n[1]*(C(3,4)*n[0] + C(3,3)*n[1] + C(2,3)*n[2]) + 
             n[2]*(C(2,4)*n[0] + C(2,3)*n[1] + C(2,2)*n[2]);

  return nCn;
}

Matrix3 
StabilityCheck::elasticPlasticAcousticTensor(const Vaango::Tensor::Matrix6Mandel& C_e,
                                             const Vaango::Tensor::Vector6Mandel& P_vec,
                                             const Vaango::Tensor::Vector6Mandel& N_vec,
                                             double H,
                                             const Vector& n) const
{
  // Compute C:N
  auto CN_vec = C_e * N_vec;

  // Compute P:N + H
  double PN_H = P_vec.transpose() * N_vec + H;

  // Convert back into Matrix3
  auto CN = Vaango::Tensor::constructMatrix3(CN_vec);
  auto P = Vaango::Tensor::constructMatrix3(P_vec);

  // Compute n . P
  Vector nP(n[0] * P(0,0) + n[1] * P(1,0) + n[2] * P(2,0),
            n[0] * P(0,1) + n[1] * P(1,1) + n[2] * P(2,1),
            n[0] * P(0,2) + n[1] * P(1,2) + n[2] * P(2,2));

  // Compute CN . n
  Vector CNn(CN(0,0) * n[0] + CN(0,1) * n[1] + CN(0,2) * n[2],
             CN(1,0) * n[0] + CN(1,1) * n[1] + CN(1,2) * n[2],
             CN(2,0) * n[0] + CN(2,1) * n[1] + CN(2,2) * n[2]);

  // Compute (n . P) otimes (CN . n)
  Matrix3 nCpn(nP, CNn);
  nCpn /= PN_H;

  // Compute elastic acoustic tensor
  auto nCn = elasticAcousticTensor(C_e, n);
  
  return nCn - nCpn; 
}

std::tuple<Matrix3, double, Matrix3>
StabilityCheck::computeAandJ(const Vaango::Tensor::Matrix6Mandel& C_e,
                             const Vaango::Tensor::Vector6Mandel& P_vec,
                             const Vaango::Tensor::Vector6Mandel& N_vec,
                             double H,
                             const Vector& n) const
{
  // Compute C:N
  auto CN_vec = C_e * N_vec;

  // Compute P:N + H
  double PN_H = P_vec.transpose() * N_vec + H;

  // Convert back into Matrix3
  auto CN = Vaango::Tensor::constructMatrix3(CN_vec);
  auto P = Vaango::Tensor::constructMatrix3(P_vec);

  // Compute n . P
  Vector nP(n[0] * P(0,0) + n[1] * P(1,0) + n[2] * P(2,0),
            n[0] * P(0,1) + n[1] * P(1,1) + n[2] * P(2,1),
            n[0] * P(0,2) + n[1] * P(1,2) + n[2] * P(2,2));

  // Compute CN . n
  Vector CNn(CN(0,0) * n[0] + CN(0,1) * n[1] + CN(0,2) * n[2],
             CN(1,0) * n[0] + CN(1,1) * n[1] + CN(1,2) * n[2],
             CN(2,0) * n[0] + CN(2,1) * n[1] + CN(2,2) * n[2]);

  // Compute (n . P) otimes (CN . n)
  Matrix3 nCpn(nP, CNn);
  nCpn /= PN_H;

  // Compute elastic acoustic tensor
  auto nCn = elasticAcousticTensor(C_e, n);

  // Compute elastic-plastic acoustic tensor
  auto A = nCn - nCpn; 
  auto detA = A.Determinant();
  auto Ainv = A.Inverse();
  auto B = Ainv.Transpose();

  // Compute elastic-plastic tangent
  auto P_CN = Vaango::Tensor::constructMatrix6Mandel(P_vec, CN_vec);
  auto C = C_e - P_CN / PN_H;

  // Compute J = det(A) * C_mkln A^{-1}_lk
  double J11 = B(0,0)*C(0,0) + B(0,1)*C(0,5) + B(0,2)*C(0,4) + 
               B(1,0)*C(5,0) + B(1,1)*C(5,5) + B(1,2)*C(5,4) + 
               B(2,0)*C(4,0) + B(2,1)*C(4,5) + B(2,2)*C(4,4);
  double J12 = B(0,0)*C(5,0) + B(0,1)*C(5,5) + B(0,2)*C(5,4) + 
               B(1,0)*C(1,0) + B(1,1)*C(1,5) + B(1,2)*C(1,4) + 
               B(2,0)*C(3,0) + B(2,1)*C(3,5) + B(2,2)*C(3,4);
  double J13 = B(0,0)*C(4,0) + B(0,1)*C(4,5) + B(0,2)*C(4,4) + 
               B(1,0)*C(3,0) + B(1,1)*C(3,5) + B(1,2)*C(3,4) + 
               B(2,0)*C(2,0) + B(2,1)*C(2,5) + B(2,2)*C(2,4);
  double J21 = B(0,0)*C(0,5) + B(0,1)*C(0,1) + B(0,2)*C(0,3) + 
               B(1,0)*C(5,5) + B(1,1)*C(5,1) + B(1,2)*C(5,3) + 
               B(2,0)*C(4,5) + B(2,1)*C(4,1) + B(2,2)*C(4,3);
  double J22 = B(0,0)*C(5,5) + B(0,1)*C(5,1) + B(0,2)*C(5,3) + 
               B(1,0)*C(1,5) + B(1,1)*C(1,1) + B(1,2)*C(1,3) + 
               B(2,0)*C(3,5) + B(2,1)*C(3,1) + B(2,2)*C(3,3);
  double J23 = B(0,0)*C(4,5) + B(0,1)*C(4,1) + B(0,2)*C(4,3) + 
               B(1,0)*C(3,5) + B(1,1)*C(3,1) + B(1,2)*C(3,3) + 
               B(2,0)*C(2,5) + B(2,1)*C(2,1) + B(2,2)*C(2,3);
  double J31 = B(0,0)*C(0,4) + B(0,1)*C(0,3) + B(0,2)*C(0,2) + 
               B(1,0)*C(5,4) + B(1,1)*C(5,3) + B(1,2)*C(5,2) + 
               B(2,0)*C(4,4) + B(2,1)*C(4,3) + B(2,2)*C(4,2);
  double J32 = B(0,0)*C(5,4) + B(0,1)*C(5,3) + B(0,2)*C(5,2) + 
               B(1,0)*C(1,4) + B(1,1)*C(1,3) + B(1,2)*C(1,2) + 
               B(2,0)*C(3,4) + B(2,1)*C(3,3) + B(2,2)*C(3,2);
  double J33 = B(0,0)*C(4,4) + B(0,1)*C(4,3) + B(0,2)*C(4,2) + 
               B(1,0)*C(3,4) + B(1,1)*C(3,3) + B(1,2)*C(3,2) + 
               B(2,0)*C(2,4) + B(2,1)*C(2,3) + B(2,2)*C(2,2);
  Matrix3 J(J11, J12, J13, J21, J22, J23, J31, J32, J33);
  J *= detA;
  
  return std::make_tuple(A, detA, J);
}
