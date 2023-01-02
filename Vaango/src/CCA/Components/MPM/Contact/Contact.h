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

#ifndef __CONTACT_H__
#define __CONTACT_H__

#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <cmath>

namespace Uintah {

using constNCint    = constNCVariable<int>;
using NCint         = NCVariable<int>;
using constNCdouble = constNCVariable<double>;
using NCdouble      = NCVariable<double>;
using constNCPoint  = constNCVariable<Point>;
using NCPoint       = NCVariable<Point>;
using constNCVector = constNCVariable<Vector>;
using NCVector      = NCVariable<Vector>;
using constNCdoubleArray =  std::vector<constNCdouble>;
using NCdoubleArray =  std::vector<NCdouble>;
using constNCVectorArray =  std::vector<constNCVector>;
using NCVectorArray =  std::vector<NCVector>;

class DataWarehouse;
class MPMLabel;
class MPMFlags;
class ProcessorGroup;
class Patch;
class VarLabel;
class Task;

class Contact : public UintahParallelComponent
{
public:
  // Constructor
  Contact(const ProcessorGroup* myworld, MPMLabel* Mlb, MPMFlags* MFlag,
          ProblemSpecP ps);
  virtual ~Contact();

  virtual void outputProblemSpec(ProblemSpecP& ps) = 0;

  // Basic contact methods
  virtual void exchangeMomentum(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw, DataWarehouse* new_dw,
                                const VarLabel* label) = 0;

  virtual void addComputesAndRequires(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls,
                                      const VarLabel* label) = 0;

  inline bool needNormals() const { return d_needNormals;}
  inline bool useLogisticRegression() const { return d_useLogisticRegression;}
  inline int  oneOrTwoStep() const {return d_oneOrTwoStep;}
    
protected:
  MPMLabel* lb;
  MPMFlags* flag;

  ContactMaterialSpec d_matls;

  bool d_needNormals;
  bool d_useLogisticRegression;
  int  d_oneOrTwoStep;
};

inline bool
compare(double num1, double num2)
{
  double EPSILON = 1.e-14;

  return (fabs(num1 - num2) <= EPSILON);
}

} // End namespace Uintah

#endif // __CONTACT_H__
