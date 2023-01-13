/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include "DruckerCheck.h"
#include <Core/ProblemSpec/ProblemSpec.h>
#include <cmath>
#include <vector>

using namespace Uintah;


DruckerCheck::DruckerCheck(ProblemSpecP&)
{
}

DruckerCheck::DruckerCheck(const DruckerCheck*)
{
}

DruckerCheck::~DruckerCheck() = default;

void
DruckerCheck::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP stability_ps = ps->appendChild("stability_check");
  stability_ps->setAttribute("type", "drucker");
}

bool
DruckerCheck::checkStability(const Matrix3&, const Matrix3& deformRate,
                             const TangentModulusTensor& Cep, Vector&)
{
  // Calculate the stress rate
  Matrix3 stressRate(0.0);
  Cep.contract(deformRate, stressRate);

  // cout << "Deform Rate = \n" << deformRate << endl;
  // cout << "Cep = \n" << Cep ;
  // cout << "Stress Rate = \n" << stressRate << endl;

  double val = stressRate.Contract(deformRate);
  // cout << "val = " << val << endl << endl;
  if (val > 0.0)
    return false;
  return true;
}

bool 
DruckerCheck::checkStability(const Matrix3& cauchyStress,
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

  // cout << "Deform Rate = \n" << d_vec << endl;
  // cout << "Cep = \n" << Cep ;
  // cout << "Stress Rate = \n" << stressRate << endl;

  double val = stressRate.transpose() * d_vec;
  // cout << "val = " << val << endl << endl;
  if (val > 0.0)
    return false;
  return true;
}
