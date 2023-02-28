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

#ifndef _INTERNAL_VARIABLE_MODEL_FACTORY_H_
#define _INTERNAL_VARIABLE_MODEL_FACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>

namespace Vaango {

class InternalVariableModel;
class ElasticModuliModel;
class ShearModulusModel;

class InternalVariableModelFactory
{
public:
  // this function has a switch for all known mat_types
  static std::unique_ptr<InternalVariableModel>
  create(Uintah::ProblemSpecP& ps);

  static std::unique_ptr<InternalVariableModel>
  create(Uintah::ProblemSpecP& ps, ElasticModuliModel* elastic);

  static std::unique_ptr<InternalVariableModel>
  create(Uintah::ProblemSpecP& ps, ShearModulusModel* shear);

  static std::unique_ptr<InternalVariableModel>
  createCopy(const InternalVariableModel* pm);
};
} // namespace Vaango

#endif /* _INTERNAL_VARIABLE_MODEL_FACTORY_H_ */
