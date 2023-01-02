/*
 * The MIT License
 *
 * Copyright (c) 1997-2014 The University of Utah
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

#ifndef UINTAH_RD_JGSCALARDIFFUSION_H
#define UINTAH_RD_JGSCALARDIFFUSION_H

#include <CCA/Components/MPM/ReactionDiffusion/ScalarDiffusionModel.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class Task;
class MPMFlags;
class MPMLabel;
class MPMMaterial;
class DataWarehouse;
class ProcessorGroup;

class JGConcentrationDiffusion : public ScalarDiffusionModel
{
public:
  JGConcentrationDiffusion(ProblemSpecP& ps,
                           MaterialManagerP& matManager,
                           MPMFlags* Mflag,
                           std::string diff_type);
  ~JGConcentrationDiffusion() = default;

  JGConcentrationDiffusion(const JGConcentrationDiffusion&) = delete;
  JGConcentrationDiffusion&
  operator=(const JGConcentrationDiffusion&) = delete;

  virtual void
  scheduleComputeFlux(Task* task,
                      const MPMMaterial* matl,
                      const PatchSet* patch) const;

  virtual void
  computeFlux(const Patch* patch,
              const MPMMaterial* matl,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw);

};

} // end namespace Uintah
#endif
