/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#ifndef _BB_SHEARMODULUSMODELFACTORY_H_
#define _BB_SHEARMODULUSMODELFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>

#include <memory>

namespace Vaango {

// Forward declarations
class ShearModulusModel;
class MPMEquationOfState;

/*! \class ShearModulusModelFactory
 *  \brief Creates instances of Shear Modulus Models
 *  \author  Biswajit Banerjee,
 *  \author  C-SAFE and Department of Mechanical Engineering,
 *  \author  University of Utah.
 */

class ShearModulusModelFactory
{

public:
  //! Create a shear modulus model from the input file problem specification.
  static std::unique_ptr<ShearModulusModel>
  create(Uintah::ProblemSpecP& ps);
  static std::unique_ptr<ShearModulusModel>
  create(Uintah::ProblemSpecP& ps, MPMEquationOfState* eos);
  static std::unique_ptr<ShearModulusModel>
  createCopy(const ShearModulusModel* yc);
};
} // namespace Vaango

#endif /* _SHEARMODULUSMODELFACTORY_H_ */
