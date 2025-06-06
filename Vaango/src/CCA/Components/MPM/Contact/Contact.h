/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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
#include <CCA/Components/MPM/Core/ContactDefs.h>

#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <cmath>

namespace Uintah {

class DataWarehouse;
class MPMLabel;
class MPMFlags;
class ProcessorGroup;
class Patch;
class VarLabel;
class Task;
class Output;

class Contact
{
public:
  // Constructor
  Contact(const ProcessorGroup* myworld,
          const MaterialManagerP& mat_manager,
          const MPMLabel* Mlb,
          const MPMFlags* MFlag,
          ProblemSpecP& ps);

  virtual ~Contact() = default;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) = 0;

  // Basic contact methods
  virtual void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) = 0;

  virtual void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) = 0;

  inline bool
  needNormals() const
  {
    return d_need_normals;
  }

  inline bool
  useLogisticRegression() const
  {
    return d_use_logistic_regression;
  }

  inline int
  oneOrTwoStep() const
  {
    return d_one_or_two_step;
  }

  // Enable setting material attributes (isRigid, needsNormals, etc)
  // based on the chosen contact model
  virtual void
  setContactMaterialAttributes() = 0;

protected:
  Output* d_output{ nullptr };

  MaterialManagerP d_mat_manager{ nullptr };
  const MPMLabel* d_mpm_labels{ nullptr };
  const MPMFlags* d_mpm_flags{ nullptr };
  ContactMaterialSpec d_matls;

  bool d_need_normals{ false };
  bool d_use_logistic_regression{ false };

  int d_one_or_two_step{ 2 };
  int d_exclude_material{ -999 };
  int d_num_ghost_particles{ 1 };
  int d_num_ghost_nodes{ 1 };

  double d_vol_const{ 0.0 };
};

inline bool
compare(double num1, double num2)
{
  double EPSILON = 1.e-14;

  return (std::abs(num1 - num2) <= EPSILON);
}

} // End namespace Uintah

#endif // __CONTACT_H__
