/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#include "DruckerBeckerCheck.h"
#include <Core/Math/SymmMatrix3.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Eigen/Dense>

#include <cmath>
#include <vector>

using namespace Uintah;
using namespace std;

DruckerBeckerCheck::DruckerBeckerCheck(ProblemSpecP&)
{
}

DruckerBeckerCheck::DruckerBeckerCheck(const DruckerBeckerCheck*)
{
}

DruckerBeckerCheck::~DruckerBeckerCheck() = default;

void
DruckerBeckerCheck::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP stability_ps = ps->appendChild("stability_check");
  stability_ps->setAttribute("type", "drucker_becker");
}

bool
DruckerBeckerCheck::checkStability(const Matrix3& stress,
                                   const Matrix3& deformRate,
                                   const TangentModulusTensor& C_ep, Vector&)
{
  // Do the Drucker stability check
  Matrix3 stressRate(0.0);
  C_ep.contract(deformRate, stressRate);
  double val = stressRate.Contract(deformRate);
  if (val <= 0.0)
    return true; // Bifurcation

  // Do the Becker check
  // Find the magnitudes and directions of the principal stresses
  SymmMatrix3 sigma(stress);
  Vector sig(0.0, 0.0, 0.0);
  Matrix3 evec;
  sigma.eigen(sig, evec);

  // Get components of C in eigenbasis
  Eigen::Matrix3d Q;
  Q << evec(0,0), evec(0,1), evec(0,2),
       evec(1,0), evec(1,1), evec(1,2),
       evec(2,0), evec(2,1), evec(2,2);
  Q = Q.transpose();
  TangentModulusTensor C;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          for (int m = 0; m < 3; ++m) {
            for (int n = 0; n < 3; ++n) {
              for (int p = 0; p < 3; ++p) {
                for (int q = 0; q < 3; ++q) {
                  C(i, j, k, l) += Q(i, m) * Q(j, n) * Q(k, p) * Q(l, q) * C_ep(m, n, p, q);
                } } } } } } } }

  // Compute quadratic terms
  double a1 = C(2,2,2,2) * C(0,2,2,0);
  double a2 = C(0,0,0,0) * C(2,2,2,2) - C(0,0,2,2) * C(2,2,0,0) * C(2,1,1,2) -
              (C(0,0,2,2) + C(2,2,0,0)) * C(2,1,1,2) * C(3,1,1,3) - 
              (C(2,1,2,1) - 1.0) * C(3,1,1,3) * C(1,3,1,3);
  double a3 = C(0,0,0,0) * C(3,1,1,3) - C(0,0,2,2) * C(2,2,0,0) * C(1,2,2,1) - 
              (C(0,0,2,2) + C(2,2,0,0)) * C(3,1,1,3) * C(1,2,2,1);

  // Solve the quadric
  // Substitute x^2 by y and solve the quadratic
  double B2_4AC = a2 * a2 - 4.0 * a1 * a3;
  if (B2_4AC < 0.0) {
    // No real roots - no bifurcation
    return false;
  } else {
    ASSERT(!(a1 == 0));
    double yplus = (-a2 + sqrt(B2_4AC)) / (2.0 * a1);
    double yminus = (-a2 - sqrt(B2_4AC)) / (2.0 * a1);
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
DruckerBeckerCheck::checkStability(const Matrix3& cauchyStress,
                                   const Matrix3& deformRate,
                                   const Vaango::Tensor::Matrix6Mandel& C_e,
                                   const Vaango::Tensor::Vector6Mandel& P_vec,
                                   const Vaango::Tensor::Vector6Mandel& N_vec,
                                   double H,
                                   Vector& direction)
{
  // Compute elastic-plastic tangent modulus
  auto CN_vec = C_e * N_vec;
  double PN_H = P_vec.transpose() * N_vec + H;
  auto P_CN = Vaango::Tensor::constructMatrix6Mandel(P_vec, CN_vec);
  auto C_ep = C_e - P_CN / PN_H;

  // Calculate the stress rate
  auto d_vec = Vaango::Tensor::constructVector6Mandel(deformRate);
  Vaango::Tensor::Vector6Mandel stressRate = C_ep * d_vec;

  double val = stressRate.transpose() * d_vec;
  if (val <= 0.0)
    return true; // Bifurcation

  // Find the magnitudes and directions of the principal stresses
  // sorted in ascending order of eigenvalues
  Eigen::Matrix3d sigma;
  sigma << cauchyStress(0,0), cauchyStress(0,1), cauchyStress(0,2),
           cauchyStress(1,0), cauchyStress(1,1), cauchyStress(1,2),
           cauchyStress(2,0), cauchyStress(2,1), cauchyStress(2,2);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(sigma);
  Eigen::Matrix3d evec = eigen_solver.eigenvectors();
  // std::cout << "Stress = \n" << sigma << "\n";
  // std::cout << "Eigenvectors : " << evec << "\n";

  // Get components of C_ep in eigenbasis
  auto Qhat = coordTransformMatrix(evec);
  auto Chat_ep = Qhat * C_ep * Qhat.transpose();

  // Convert components of C_ep to Voigt M_{ijkl} form
  Vaango::Tensor::Matrix6Mandel C(Chat_ep);
  C.block<3,3>(0,3) *= Vaango::Tensor::one_sqrt_two;
  C.block<3,3>(3,0) *= Vaango::Tensor::one_sqrt_two;
  C.block<3,3>(3,3) *= Vaango::Tensor::half;

  // Set up quadric equation
  double a1 = C(2,2) * C(4,4);
  double a2 = C(0,0) * C(2,2) - C(0,2) * C(2,0) * C(3,3) -
              (C(0,2) + C(2,0)) * C(3,3) * C(4,4) - 
              (C(3,3) - 1.0) * C(4,4) * C(4,4);
  double a3 = C(0,0) * C(4,4) - C(0,2) * C(2,0) * C(5,5) - 
              (C(0,2) + C(2,0)) * C(4,4) * C(5,5);

  // Solve the quadric
  // Substitute x^2 by y and solve the quadratic
  double B2_4AC = a2 * a2 - 4.0 * a1 * a3;
  if (B2_4AC < 0.0) {
    // Some imaginary roots - no bifurcation
    return false;
  } else {
    ASSERT(!(a1 == 0));
    double yplus = (-a2 + sqrt(B2_4AC)) / (2.0 * a1);
    double yminus = (-a2 - sqrt(B2_4AC)) / (2.0 * a1);
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
