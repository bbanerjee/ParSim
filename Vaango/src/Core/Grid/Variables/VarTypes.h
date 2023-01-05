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

#ifndef __CORE_GRID_VARIABLES_VARTYPES_H__
#define __CORE_GRID_VARIABLES_VARTYPES_H__

#include <Core/Disclosure/TypeUtils.h>
#include <Core/Grid/Variables/ReductionVariable.h>
#include <Core/Grid/Variables/Reductions.h>
#include <Core/Grid/Variables/SoleVariable.h>

namespace Uintah {

// System vars related to the application
const std::string timeStep_name("timeStep");
const std::string simTime_name("simulationTime");
const std::string delT_name("delT");

const std::string outputInterval_name("outputInterval");
const std::string checkpointInterval_name("checkpointInterval");
const std::string outputTimeStep_name("outputTimeStep");
const std::string outputPreviousTimeStep_name("outputPreviousTimeStep");
const std::string checkpointTimeStep_name("checkpointTimeStep");
const std::string checkpointPreviousTimeStep_name("checkpointPreviousTimeStep");
const std::string recomputeTimeStep_name("recomputeTimeStep");
const std::string abortTimeStep_name("abortTimeStep");
const std::string endSimulation_name("endSimulation");

using timeStep_vartype = SoleVariable<unsigned int>;
using simTime_vartype  = SoleVariable<double>;

using delt_vartype = ReductionVariable<double, Reductions::Min<double>>;

using max_vartype = ReductionVariable<double, Reductions::Max<double>>;

using min_vartype = ReductionVariable<double, Reductions::Min<double>>;

using sum_vartype = ReductionVariable<double, Reductions::Sum<double>>;

using bool_and_vartype = ReductionVariable<bool, Reductions::And<bool>>;

using bool_or_vartype = ReductionVariable<bool, Reductions::Or<bool>>;

using minvec_vartype = ReductionVariable<Vector, Reductions::Min<Vector>>;

using maxvec_vartype = ReductionVariable<Vector, Reductions::Max<Vector>>;

using sumvec_vartype = ReductionVariable<Vector, Reductions::Sum<Vector>>;

using sumlong_vartype = ReductionVariable<long64, Reductions::Sum<long64>>;

using sumlonglong_vartype =
    ReductionVariable<long64, Reductions::Sum<long long>>;

}  // End namespace Uintah

#endif  //__CORE_GRID_VARIABLES_VARTYPES_H__
