/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __CCA_COMPONENTS_MPM_CORE_CONTACT_DEFS_H__
#define __CCA_COMPONENTS_MPM_CORE_CONTACT_DEFS_H__

#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ComputeSet.h>

namespace Uintah {

using constNCint         = constNCVariable<int>;
using NCint              = NCVariable<int>;
using constNCdouble      = constNCVariable<double>;
using NCdouble           = NCVariable<double>;
using constNCPoint       = constNCVariable<Point>;
using NCPoint            = NCVariable<Point>;
using constNCVector      = constNCVariable<Vector>;
using NCVector           = NCVariable<Vector>;
using constNCdoubleArray = std::vector<constNCdouble>;
using NCdoubleArray      = std::vector<NCdouble>;
using constNCVectorArray = std::vector<constNCVector>;
using NCVectorArray      = std::vector<NCVector>;

} // End namespace Uintah

#endif // __CCA_COMPONENTS_MPM_CORE_CONTACT_DEFS_H__
