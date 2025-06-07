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

#include <CCA/Components/Solvers/CGSolver.h>
#include <CCA/Components/Solvers/CGStencil7.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/Stencil7.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Timers/Timers.hpp>

#include <iomanip>

using namespace Uintah;
//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "CGSOLVER_DOING_COUT:+"

static DebugStream cout_doing("CGSOLVER_DOING_COUT", false);

namespace Uintah {

CGSolver::CGSolver(const ProcessorGroup* myworld)
  : SolverCommon(myworld)
{
  d_params = std::make_unique<CGSolverParams>();
}

void
CGSolver::readParameters(ProblemSpecP& params_ps, const string& varname)
{
  if (params_ps) {
    for (ProblemSpecP param_ps = params_ps->findBlock("Parameters");
         param_ps != nullptr;
         param_ps = param_ps->findNextBlock("Parameters")) {
      std::string variable;
      if (param_ps->getAttribute("variable", variable) && variable != varname) {
        continue;
      }
      param_ps->get("initial_tolerance", d_params->initial_tolerance);
      param_ps->get("tolerance", d_params->tolerance);
      param_ps->getWithDefault("maxiterations", d_params->maxiterations, 75);

      std::string norm;
      if (param_ps->get("norm", norm)) {
        if (norm == "L1" || norm == "l1") {
          d_params->norm = CGSolverParams::Norm::L1;
        } else if (norm == "L2" || norm == "l2") {
          d_params->norm = CGSolverParams::Norm::L2;
        } else if (norm == "LInfinity" || norm == "linfinity") {
          d_params->norm = CGSolverParams::Norm::LInfinity;
        } else {
          throw ProblemSetupException("Unknown norm type: ",
                                      __FILE__,
                                      __LINE__);
        }
      }
      std::string criteria;
      if (param_ps->get("criteria", criteria)) {
        if (criteria == "Absolute" || criteria == "absolute") {
          d_params->criteria = CGSolverParams::Criteria::Absolute;
        } else if (criteria == "Relative" || criteria == "relative") {
          d_params->criteria = CGSolverParams::Criteria::Relative;
        } else {
          throw ProblemSetupException("Unknown criteria: ", __FILE__, __LINE__);
        }
      }
    }
  }

  if (d_params->norm == CGSolverParams::Norm::L2) {
    d_params->tolerance *= d_params->tolerance;
  }
}

void
CGSolver::scheduleSolve(const LevelP& level,
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
                        bool isFirstSolve)
{
  Task* task;
  // The extra handle arg ensures that the stencil7 object will get freed
  // when the task gets freed.  The downside is that the refcount gets
  // tweaked everytime solve is called.

  TypeDescription::Type domtype = A->typeDescription()->getType();
  ASSERTEQ(
    static_cast<std::underlying_type<TypeDescription::Type>::type>(domtype),
    static_cast<std::underlying_type<TypeDescription::Type>::type>(
      x->typeDescription()->getType()));
  ASSERTEQ(
    static_cast<std::underlying_type<TypeDescription::Type>::type>(domtype),
    static_cast<std::underlying_type<TypeDescription::Type>::type>(
      b->typeDescription()->getType()));

  Ghost::GhostType Around;

  switch (domtype) {
    case TypeDescription::Type::SFCXVariable: {
      Around = Ghost::AroundFaces;
      CGStencil7<SFCXTypes>* that =
        scinew CGStencil7<SFCXTypes>(sched.get_rep(),
                                     d_myworld,
                                     level.get_rep(),
                                     matls,
                                     Around,
                                     A,
                                     which_A_dw,
                                     x,
                                     modifies_x,
                                     b,
                                     which_b_dw,
                                     guess,
                                     which_guess_dw,
                                     d_params.get());
      Handle<CGStencil7<SFCXTypes>> handle = that;
      task = scinew Task("CGSolver::Matrix solve(SFCX)",
                         that,
                         &CGStencil7<SFCXTypes>::solve,
                         handle);
    } break;
    case TypeDescription::Type::SFCYVariable: {
      Around = Ghost::AroundFaces;
      CGStencil7<SFCYTypes>* that =
        scinew CGStencil7<SFCYTypes>(sched.get_rep(),
                                     d_myworld,
                                     level.get_rep(),
                                     matls,
                                     Around,
                                     A,
                                     which_A_dw,
                                     x,
                                     modifies_x,
                                     b,
                                     which_b_dw,
                                     guess,
                                     which_guess_dw,
                                     d_params.get());
      Handle<CGStencil7<SFCYTypes>> handle = that;
      task = scinew Task("CGSolver::Matrix solve(SFCY)",
                         that,
                         &CGStencil7<SFCYTypes>::solve,
                         handle);
    } break;
    case TypeDescription::Type::SFCZVariable: {
      Around = Ghost::AroundFaces;
      CGStencil7<SFCZTypes>* that =
        scinew CGStencil7<SFCZTypes>(sched.get_rep(),
                                     d_myworld,
                                     level.get_rep(),
                                     matls,
                                     Around,
                                     A,
                                     which_A_dw,
                                     x,
                                     modifies_x,
                                     b,
                                     which_b_dw,
                                     guess,
                                     which_guess_dw,
                                     d_params.get());
      Handle<CGStencil7<SFCZTypes>> handle = that;
      task = scinew Task("CGSolver::Matrix solve(SFCZ)",
                         that,
                         &CGStencil7<SFCZTypes>::solve,
                         handle);
    } break;
    case TypeDescription::Type::CCVariable: {
      Around                    = Ghost::AroundCells;
      CGStencil7<CCTypes>* that = scinew CGStencil7<CCTypes>(sched.get_rep(),
                                                             d_myworld,
                                                             level.get_rep(),
                                                             matls,
                                                             Around,
                                                             A,
                                                             which_A_dw,
                                                             x,
                                                             modifies_x,
                                                             b,
                                                             which_b_dw,
                                                             guess,
                                                             which_guess_dw,
                                                             d_params.get());
      Handle<CGStencil7<CCTypes>> handle = that;
      task = scinew Task("CGSolver::Matrix solve(CC)",
                         that,
                         &CGStencil7<CCTypes>::solve,
                         handle);
    } break;
    case TypeDescription::Type::NCVariable: {
      Around                    = Ghost::AroundNodes;
      CGStencil7<NCTypes>* that = scinew CGStencil7<NCTypes>(sched.get_rep(),
                                                             d_myworld,
                                                             level.get_rep(),
                                                             matls,
                                                             Around,
                                                             A,
                                                             which_A_dw,
                                                             x,
                                                             modifies_x,
                                                             b,
                                                             which_b_dw,
                                                             guess,
                                                             which_guess_dw,
                                                             d_params.get());
      Handle<CGStencil7<NCTypes>> handle = that;
      task = scinew Task("CGSolver::Matrix solve(NC)",
                         that,
                         &CGStencil7<NCTypes>::solve,
                         handle);
    } break;
    default:
      throw InternalError("Unknown variable type in scheduleSolve",
                          __FILE__,
                          __LINE__);
  }

  task->needs(which_A_dw, A, Ghost::None, 0);
  if (guess) {
    task->needs(which_guess_dw, guess, Around, 1);
  }
  if (modifies_x) {
    task->modifies(x);
  } else {
    task->computes(x);
  }

  task->needs(which_b_dw, b, Ghost::None, 0);
  task->hasSubScheduler();

  if (d_params->getRecomputeTimestepOnFailure()) {
    task->computes(VarLabel::find(abortTimestep_name));
    task->computes(VarLabel::find(recomputeTimestep_name));
  }

  LoadBalancer* lb                = sched->getLoadBalancer();
  const PatchSet* perproc_patches = lb->getPerProcessorPatchSet(level);
  sched->addTask(task, perproc_patches, matls);
}

std::string
CGSolver::getName()
{
  return "CGSolver";
}

} // end namespace Uintah
