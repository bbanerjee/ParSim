/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Solvers_CGSolverParams_h
#define Packages_Uintah_CCA_Components_Solvers_CGSolverParams_h

#include <CCA/Ports/SolverInterface.h>

namespace Uintah {

class CGSolverParams : public SolverParameters
{
public:
  double tolerance{ 1.0e-8 };
  double initial_tolerance{ 1.0e-15 };
  int maxiterations;

  enum class Norm
  {
    L1,
    L2,
    LInfinity
  };

  Norm norm;

  enum class Criteria
  {
    Absolute,
    Relative
  };

  Criteria criteria;

  CGSolverParams()
    : norm(Norm::L2)
    , criteria(Criteria::Relative)
  {
  }

  ~CGSolverParams() = default;
};


} // end namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_CGSolverParams_h
