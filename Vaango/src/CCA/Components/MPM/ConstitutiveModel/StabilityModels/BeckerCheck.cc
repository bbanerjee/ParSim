/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include "BeckerCheck.h"
#include <Core/Math/SymmMatrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Eigen/Dense>

#include <cmath>
#include <vector>

using namespace Uintah;

BeckerCheck::BeckerCheck(ProblemSpecP&)
{
}

BeckerCheck::BeckerCheck(const BeckerCheck*)
{
}

BeckerCheck::~BeckerCheck() = default;

void
BeckerCheck::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP stability_ps = ps->appendChild("stability_check");
  stability_ps->setAttribute("type", "becker");
}

bool
BeckerCheck::checkStability(const Matrix3& stress, const Matrix3&,
                            const TangentModulusTensor& M, Vector&)
{
  // Find the magnitudes and directions of the principal stresses
  // sorted in descending order of eigenvalues
  SymmMatrix3 sigma(stress);
  Vector sig(0.0, 0.0, 0.0);
  Matrix3 evec;
  sigma.eigen(sig, evec);
  // std::cout << "stress = \n";
  // std::cout << stress << "\n";
  // std::cout << "Eigenvalues : " << sig << "\n";
  // std::cout << "Eigenvectors : " << evec << "\n";

  // Get components of M in eigenbasis
  Eigen::Matrix3d Q;
  Q << evec(0,0), evec(0,1), evec(0,2),
       evec(1,0), evec(1,1), evec(1,2),
       evec(2,0), evec(2,1), evec(2,2);
  Q = Q.transpose();
  TangentModulusTensor M_prime;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          for (int m = 0; m < 3; ++m) {
            for (int n = 0; n < 3; ++n) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  M_prime(i, j, k, l) += Q(i, m) * Q(j, n) * Q(k, p) * Q(l, q) * M(m, n, p, q);
                } } } } } } } }

  // Compute quadratic terms
  double A = M_prime(2, 0, 2, 0) * (-sig[2] + 2.0 * M_prime(2, 2, 2, 2));
  double C = M_prime(2, 0, 2, 0) * (-sig[0] + 2.0 * M_prime(0, 0, 0, 0));
  double B =
    M_prime(2, 0, 2, 0) *
      (sig[2] - 2.0 * M_prime(0, 0, 2, 2) + sig[0] - 2.0 * M_prime(2, 2, 0, 0)) +
    sig[0] * (M_prime(2, 2, 0, 0) - M_prime(2, 2, 2, 2)) +
    sig[2] * (M_prime(0, 0, 2, 2) - M_prime(0, 0, 0, 0)) +
    2.0 * (-M_prime(0, 0, 2, 2) * M_prime(2, 2, 0, 0) + M_prime(0, 0, 0, 0) * M_prime(2, 2, 2, 2));

  // Solve the quadric
  // Substitute x^2 by y and solve the quadratic
  double B2_4AC = B * B - 4.0 * A * C;
  if (B2_4AC < 0.0) {
    // No real roots - no bifurcation
    return false;
  } else {
    ASSERT(!(A == 0));
    double yplus = (-B + sqrt(B2_4AC)) / (2.0 * A);
    double yminus = (-B - sqrt(B2_4AC)) / (2.0 * A);
    if (yplus < 0.0 && yminus < 0.0) {
      // No real roots - no bifurcation
      return false;
    } else {
      if (yplus < 0.0 || yminus < 0.0) {
        // Two real roots -  bifurcation ? (parabolic)
        return false;
      }
    }
  }
  // Four real roots -  bifurcation
  return true;
}

bool 
BeckerCheck::checkStability(const Matrix3& stress,
                            const Matrix3& /*deformRate*/,
                            const Vaango::Tensor::Matrix6Mandel& C_e,
                            const Vaango::Tensor::Vector6Mandel& P_vec,
                            const Vaango::Tensor::Vector6Mandel& N_vec,
                            double H,
                            Vector& /*direction*/)
{
  // Find the magnitudes and directions of the principal stresses
  // sorted in ascending order of eigenvalues
  Eigen::Matrix3d sigma;
  sigma << stress(0,0), stress(0,1), stress(0,2),
           stress(1,0), stress(1,1), stress(1,2),
           stress(2,0), stress(2,1), stress(2,2);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(sigma);
  Eigen::Vector3d eval = eigen_solver.eigenvalues();
  Eigen::Matrix3d evec = eigen_solver.eigenvectors();
  double sig_max = eval[2];
  double sig_min = eval[0];
  // std::cout << "stress = \n";
  // std::cout << stress << "\n";
  // std::cout << "Eigenvalues : " << eval << "\n";
  // std::cout << "Eigenvectors : " << evec << "\n";

  // Compute tangent modulus
  auto CN_vec = C_e * N_vec;
  double PN_H = P_vec.transpose() * N_vec + H;
  auto P_CN = Vaango::Tensor::constructMatrix6Mandel(P_vec, CN_vec);
  auto C_ep = C_e - P_CN / PN_H;

  // Get components of C_ep in eigenbasis
  auto Qhat = coordTransformMatrix(evec);
  auto Chat_ep = Qhat * C_ep * Qhat.transpose();

  // Convert components of C_ep to M_{ijkl} form
  Vaango::Tensor::Matrix6Mandel M(Chat_ep);
  M.block<3,3>(0,3) *= Vaango::Tensor::one_sqrt_two;
  M.block<3,3>(3,0) *= Vaango::Tensor::one_sqrt_two;
  M.block<3,3>(3,3) *= Vaango::Tensor::half;

  // Set up quadric equation
  double A = M(5, 5) * (-sig_min + 2.0 * M(2, 2));
  double C = M(5, 5) * (-sig_max + 2.0 * M(0, 0));
  double B =
    M(5, 5) *
      (sig_min - 2.0 * M(0, 2) + sig_max - 2.0 * M(2, 0)) +
    sig_max * (M(2, 0) - M(2, 2)) +
    sig_min * (M(0, 2) - M(0, 0)) +
    2.0 * (-M(0, 2) * M(2, 0) + M(0, 0) * M(2, 2));

  // Solve the quadric
  // Substitute x^2 by y and solve the quadratic
  double B2_4AC = B * B - 4.0 * A * C;
  if (B2_4AC < 0.0) {
    // No real roots - no bifurcation
    return false;
  } else {
    ASSERT(!(A == 0));
    double yplus = (-B + sqrt(B2_4AC)) / (2.0 * A);
    double yminus = (-B - sqrt(B2_4AC)) / (2.0 * A);
    if (yplus < 0.0 && yminus < 0.0) {
      // No real roots - no bifurcation
      return false;
    } else {
      if (yplus < 0.0 || yminus < 0.0) {
        // Two real roots -  bifurcation ? (parabolic)
        return false;
      }
    }
  }
  // Four real roots -  bifurcation
  return true;
}
