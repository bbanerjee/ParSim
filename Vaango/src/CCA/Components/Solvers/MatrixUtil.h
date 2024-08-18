/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef Packages_Uintah_CCA_Components_Solvers_MatrixUtil_h
#define Packages_Uintah_CCA_Components_Solvers_MatrixUtil_h

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/Stencil4.h>
#include <Core/Grid/Variables/Stencil7.h>

namespace Uintah {
class SFCXTypes
{
public:
  using matrix_type           = constSFCXVariable<Stencil7>;
  using symmetric_matrix_type = constSFCXVariable<Stencil4>;
  using const_double_type     = constSFCXVariable<double>;
  using double_type           = SFCXVariable<double>;
};

class SFCYTypes
{
public:
  using matrix_type           = constSFCYVariable<Stencil7>;
  using symmetric_matrix_type = constSFCYVariable<Stencil4>;
  using const_double_type     = constSFCYVariable<double>;
  using double_type           = SFCYVariable<double>;
};

class SFCZTypes
{
public:
  using matrix_type           = constSFCZVariable<Stencil7>;
  using symmetric_matrix_type = constSFCZVariable<Stencil4>;
  using const_double_type     = constSFCZVariable<double>;
  using double_type           = SFCZVariable<double>;
};

class CCTypes
{
public:
  using matrix_type           = constCCVariable<Stencil7>;
  using symmetric_matrix_type = constCCVariable<Stencil4>;
  using const_double_type     = constCCVariable<double>;
  using double_type           = CCVariable<double>;
};

class NCTypes
{
public:
  using matrix_type           = constNCVariable<Stencil7>;
  using symmetric_matrix_type = constNCVariable<Stencil4>;
  using const_double_type     = constNCVariable<double>;
  using double_type           = NCVariable<double>;
};
} // namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_MatrixUtil_h
