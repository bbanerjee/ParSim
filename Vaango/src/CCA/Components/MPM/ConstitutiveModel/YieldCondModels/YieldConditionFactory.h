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

#ifndef _MPM_CONSTITUTIVEMODELS_MODELS_YIELDCONDITIONFACTORY_H_
#define _MPM_CONSTITUTIVEMODELS_MODELS_YIELDCONDITIONFACTORY_H_

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {
// Forward declarations
class FlowStressModel;
} // namespace Uintah

namespace Vaango {

// Forward declarations
class YieldCondition;
class IntVar_Metal;
class IntVar_BorjaPressure;
class IntVar_TabularCap;

/*! \class YieldConditionFactory
 *  \brief Creates instances of Yield Conditions
 *  \author  Biswajit Banerjee,
 *  \author  C-SAFE and Department of Mechanical Engineering,
 *  \author  University of Utah.
 *  \warning Currently implemented yield conditions:
 *           von Mises, Gurson-Tvergaard-Needleman, Rousselier
 */

class YieldConditionFactory
{

public:
  //! Create a yield condition from the input file problem specification.
  static std::unique_ptr<YieldCondition>
  create(Uintah::ProblemSpecP& ps);

  static std::unique_ptr<YieldCondition>
  create(Uintah::ProblemSpecP& ps,
         IntVar_Metal* intvar,
         const Uintah::FlowStressModel* flow);

  static std::unique_ptr<YieldCondition>
  create(Uintah::ProblemSpecP& ps, IntVar_BorjaPressure* intvar);

  static std::unique_ptr<YieldCondition>
  create(Uintah::ProblemSpecP& ps, IntVar_TabularCap* intvar);

  static std::unique_ptr<YieldCondition>
  createCopy(const YieldCondition* yc);
};
} // namespace Vaango

#endif /* _MPM_CONSTITUTIVEMODELS_MODELS_YIELDCONDITIONFACTORY_H_ */
