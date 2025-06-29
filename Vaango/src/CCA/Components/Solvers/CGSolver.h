/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef Packages_Uintah_CCA_Components_Solvers_CGSolver_h
#define Packages_Uintah_CCA_Components_Solvers_CGSolver_h

#include <CCA/Components/Solvers/SolverCommon.h>
#include <CCA/Components/Solvers/CGSolverParams.h>

#include <memory>

namespace Uintah {

class CGSolver : public SolverCommon
{

public:
  CGSolver(const ProcessorGroup* myworld);
  virtual ~CGSolver() = default;

  virtual void
  readParameters(ProblemSpecP& params, const std::string& name);

  virtual SolverParameters*
  getParameters()
  {
    return d_params.get();
  }

  virtual void
  scheduleSolve(const LevelP& level,
                SchedulerP& sched,
                const MaterialSet* matls,
                const VarLabel* A,
                Task::WhichDW which_A_dw,
                const VarLabel* x,
                bool modifies_x,
                const VarLabel* b,
                Task::WhichDW which_b_dw,
                const VarLabel* guess,
                Task::WhichDW which_guess_dw,
                bool isFirstSolve = true);

  virtual std::string
  getName();

  // CGSolver does not require initialization... but we need an empty
  // routine to satisfy inheritance.
  virtual void
  scheduleInitialize([[maybe_unused]] const LevelP& level,
                     [[maybe_unused]] SchedulerP& sched,
                     [[maybe_unused]] const MaterialSet* matls)
  {
  }

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level,
                            [[maybe_unused]] SchedulerP& sched,
                            [[maybe_unused]] const MaterialSet* matls)
  {
  }

private:
  std::unique_ptr<CGSolverParams> d_params{ nullptr };
};

} // end namespace Uintah

#endif // Packages_Uintah_CCA_Components_Solvers_CGSolver_h
