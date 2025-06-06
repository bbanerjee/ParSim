/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <CCA/Components/Solvers/CGSolverParams.h>
#include <CCA/Components/Solvers/MatrixUtil.h>
#include <CCA/Components/Solvers/SolverUtils.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Util/DebugStream.h>
#include <Core/Util/Timers/Timers.hpp>

namespace Uintah {

//  To turn on normal output
//  setenv SCI_DEBUG "CGSTENCIL7_COUT:+"
inline Uintah::DebugStream cout_doing("CGSTENCIL7_COUT", false);

template<class GridVarType>
class CGStencil7 : public RefCounted
{
public:
  CGStencil7(Scheduler* sched,
             const ProcessorGroup* world,
             const Level* level,
             const MaterialSet* matlset,
             Ghost::GhostType Around,
             const VarLabel* A,
             Task::WhichDW which_A_dw,
             const VarLabel* x,
             bool modifies_x,
             const VarLabel* b,
             Task::WhichDW which_b_dw,
             const VarLabel* guess,
             Task::WhichDW which_guess_dw,
             const CGSolverParams* params)
    : sched(sched)
    , world(world)
    , level(level)
    , matlset(matlset)
    , Around(Around)
    , A_label(A)
    , which_A_dw(which_A_dw)
    , X_label(x)
    , B_label(b)
    , which_b_dw(which_b_dw)
    , guess_label(guess)
    , which_guess_dw(which_guess_dw)
    , params(params)
    , modifies_x(modifies_x)
  {
    switch (which_A_dw) {
      case Task::OldDW:
        parent_which_A_dw = Task::ParentOldDW;
        break;
      case Task::NewDW:
        parent_which_A_dw = Task::ParentNewDW;
        break;
      default:
        throw ProblemSetupException("Unknown data warehouse for A matrix",
                                    __FILE__,
                                    __LINE__);
    }
    switch (which_b_dw) {
      case Task::OldDW:
        parent_which_b_dw = Task::ParentOldDW;
        break;
      case Task::NewDW:
        parent_which_b_dw = Task::ParentNewDW;
        break;
      default:
        throw ProblemSetupException("Unknown data warehouse for b rhs",
                                    __FILE__,
                                    __LINE__);
    }
    switch (which_guess_dw) {
      case Task::OldDW:
        parent_which_guess_dw = Task::ParentOldDW;
        break;
      case Task::NewDW:
        parent_which_guess_dw = Task::ParentNewDW;
        break;
      default:
        throw ProblemSetupException("Unknown data warehouse for initial guess",
                                    __FILE__,
                                    __LINE__);
    }
    using double_type =  typename GridVarType::double_type;
    R_label =
      VarLabel::create(A->getName() + " R", double_type::getTypeDescription());
    D_label =
      VarLabel::create(A->getName() + " D", double_type::getTypeDescription());
    Q_label =
      VarLabel::create(A->getName() + " Q", double_type::getTypeDescription());
    d_label =
      VarLabel::create(A->getName() + " d", sum_vartype::getTypeDescription());
    aden_label = VarLabel::create(A->getName() + " aden",
                                  sum_vartype::getTypeDescription());
    diag_label = VarLabel::create(A->getName() + " inverse diagonal",
                                  double_type::getTypeDescription());

    tolerance_label =
      VarLabel::create("tolerance", sum_vartype::getTypeDescription());

    VarLabel* tmp_flop_label =
      VarLabel::create(A->getName() + " flops",
                       sumlong_vartype::getTypeDescription());
    tmp_flop_label->isReductionTask(false);
    flop_label = tmp_flop_label;

    VarLabel* tmp_memref_label =
      VarLabel::create(A->getName() + " memrefs",
                       sumlong_vartype::getTypeDescription());
    tmp_memref_label->isReductionTask(false);
    memref_label = tmp_memref_label;

    switch (params->norm) {
      case CGSolverParams::Norm::L1:
        err_label = VarLabel::create(A->getName() + " err",
                                     sum_vartype::getTypeDescription());
        break;
      case CGSolverParams::Norm::L2:
        err_label = d_label; // Uses d
        break;
      case CGSolverParams::Norm::LInfinity:
        err_label = VarLabel::create(A->getName() + " err",
                                     max_vartype::getTypeDescription());
        break;
      default:
        break;
    }
  }

  virtual ~CGStencil7()
  {
    VarLabel::destroy(R_label);
    VarLabel::destroy(D_label);
    VarLabel::destroy(Q_label);
    VarLabel::destroy(d_label);
    VarLabel::destroy(diag_label);
    VarLabel::destroy(flop_label);
    VarLabel::destroy(memref_label);
    VarLabel::destroy(tolerance_label);

    if (err_label != d_label) {
      VarLabel::destroy(err_label);
    }
    VarLabel::destroy(aden_label);
  }
  //______________________________________________________________________
  //
  void
  step1(const ProcessorGroup*,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse* old_dw,
        DataWarehouse* new_dw)
  {
    DataWarehouse* A_dw = new_dw->getOtherDataWarehouse(parent_which_A_dw);
    // Step 1 - requires A(parent), D(old, 1 ghost) computes aden(new)
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      for (int m = 0; m < matls->size(); m++) {
        int matl = matls->get(m);

        typename GridVarType::double_type Q;
        new_dw->allocateAndPut(Q, Q_label, matl, patch);

        typename GridVarType::matrix_type A;
        A_dw->get(A, A_label, matl, patch, Ghost::None, 0);

        typename GridVarType::const_double_type D;
        old_dw->get(D, D_label, matl, patch, Around, 1);

        typedef typename GridVarType::double_type double_type;
        Patch::VariableBasis basis = Patch::translateTypeToBasis(
          double_type::getTypeDescription()->getType(),
          true);

        IntVector l, h;
        if (params->getSolveOnExtraCells()) {
          l = patch->getExtraLowIndex(basis, IntVector(0, 0, 0));
          h = patch->getExtraHighIndex(basis, IntVector(0, 0, 0));
        } else {
          l = patch->getLowIndex(basis);
          h = patch->getHighIndex(basis);
        }

        CellIterator iter(l, h);

        IntVector ll(l);
        IntVector hh(h);
        ll -=
          IntVector(patch->getBCType(Patch::xminus) == Patch::Neighbor ? 1 : 0,
                    patch->getBCType(Patch::yminus) == Patch::Neighbor ? 1 : 0,
                    patch->getBCType(Patch::zminus) == Patch::Neighbor ? 1 : 0);

        hh +=
          IntVector(patch->getBCType(Patch::xplus) == Patch::Neighbor ? 1 : 0,
                    patch->getBCType(Patch::yplus) == Patch::Neighbor ? 1 : 0,
                    patch->getBCType(Patch::zplus) == Patch::Neighbor ? 1 : 0);
        hh -= IntVector(1, 1, 1);

        // Q = A*D
        long64 flops   = 0;
        long64 memrefs = 0;
        // Must be qualified with :: for the IBM xlC compiler.
        double aden;
        Uintah::Solver::Mult(Q, A, D, iter, ll, hh, flops, memrefs, aden);
        new_dw->put(sum_vartype(aden), aden_label);

        new_dw->put(sumlong_vartype(flops), flop_label);
        new_dw->put(sumlong_vartype(memrefs), memref_label);
      }
    }
  }
  //______________________________________________________________________
  //
  void
  step2(const ProcessorGroup*,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse* old_dw,
        DataWarehouse* new_dw)
  {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      if (cout_doing.active()) {
        cout_doing << "CGSolver::step2 on patch" << patch->getID() << std::endl;
      }
      for (int m = 0; m < matls->size(); m++) {
        int matl = matls->get(m);
        typedef typename GridVarType::double_type double_type;
        Patch::VariableBasis basis = Patch::translateTypeToBasis(
          double_type::getTypeDescription()->getType(),
          true);

        IntVector l, h;
        if (params->getSolveOnExtraCells()) {
          l = patch->getExtraLowIndex(basis, IntVector(0, 0, 0));
          h = patch->getExtraHighIndex(basis, IntVector(0, 0, 0));
        } else {
          l = patch->getLowIndex(basis);
          h = patch->getHighIndex(basis);
        }
        CellIterator iter(l, h);

        // Step 2 - requires d(old), aden(new) D(old), X(old) R(old)  computes
        // X, R, Q, d
        typename GridVarType::const_double_type D;
        old_dw->get(D, D_label, matl, patch, Ghost::None, 0);

        typename GridVarType::const_double_type X, R, diagonal;
        old_dw->get(X, X_label, matl, patch, Ghost::None, 0);
        old_dw->get(R, R_label, matl, patch, Ghost::None, 0);
        old_dw->get(diagonal, diag_label, matl, patch, Ghost::None, 0);

        typename GridVarType::double_type Xnew, Rnew;
        new_dw->allocateAndPut(Xnew, X_label, matl, patch, Ghost::None, 0);
        new_dw->allocateAndPut(Rnew, R_label, matl, patch, Ghost::None, 0);

        typename GridVarType::double_type Q;
        new_dw->getModifiable(Q, Q_label, matl, patch);

        sum_vartype aden;
        new_dw->get(aden, aden_label);

        sum_vartype d;
        old_dw->get(d, d_label);

        long64 flops   = 0;
        long64 memrefs = 0;
        double a       = d / aden;

        // X = a*D+X
        Uintah::Solver::ScMult_Add(Xnew, a, D, X, iter, flops, memrefs);
        // R = -a*Q+R
        Uintah::Solver::ScMult_Add(Rnew, -a, Q, R, iter, flops, memrefs);

        // Simple Preconditioning...
        Uintah::Solver::Mult(Q, Rnew, diagonal, iter, flops, memrefs);

        // Calculate coefficient bk and direction vectors p and pp
        double dnew = Uintah::Solver::Dot(Q, Rnew, iter, flops, memrefs);

        // Calculate error term

        switch (params->norm) {
          case CGSolverParams::Norm::L1: {
            double err = Uintah::Solver::L1(Q, iter, flops, memrefs);
            new_dw->put(sum_vartype(err), err_label);
          } break;
          case CGSolverParams::Norm::L2:
            // Nothing...
            break;
          case CGSolverParams::Norm::LInfinity: {
            double err = Uintah::Solver::LInf(Q, iter, flops, memrefs);
            new_dw->put(max_vartype(err), err_label);
          } break;
          default:
            break;
        }
        new_dw->put(sum_vartype(dnew), d_label);
        new_dw->put(sumlong_vartype(flops), flop_label);
        new_dw->put(sumlong_vartype(memrefs), memref_label);
      }
    }
    new_dw->transferFrom(old_dw, diag_label, patches, matls);
  }
  //______________________________________________________________________
  //
  void
  step3(const ProcessorGroup*,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse* old_dw,
        DataWarehouse* new_dw)
  {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      if (cout_doing.active()) {
        cout_doing << "CGSolver::step3 on patch" << patch->getID() << std::endl;
      }
      for (int m = 0; m < matls->size(); m++) {
        int matl = matls->get(m);
        typedef typename GridVarType::double_type double_type;
        Patch::VariableBasis basis = Patch::translateTypeToBasis(
          double_type::getTypeDescription()->getType(),
          true);

        IntVector l, h;
        if (params->getSolveOnExtraCells()) {
          l = patch->getExtraLowIndex(basis, IntVector(0, 0, 0));
          h = patch->getExtraHighIndex(basis, IntVector(0, 0, 0));
        } else {
          l = patch->getLowIndex(basis);
          h = patch->getHighIndex(basis);
        }
        CellIterator iter(l, h);

        sum_vartype dnew, dold;
        old_dw->get(dold, d_label);
        new_dw->get(dnew, d_label);
        typename GridVarType::const_double_type Q;
        new_dw->get(Q, Q_label, matl, patch, Ghost::None, 0);

        typename GridVarType::const_double_type D;
        old_dw->get(D, D_label, matl, patch, Ghost::None, 0);

        // Step 3 - requires D(old), Q(new), d(new), d(old), computes D
        double b = dnew / dold;

        // D = b*D+Q
        typename GridVarType::double_type Dnew;
        new_dw->allocateAndPut(Dnew, D_label, matl, patch, Ghost::None, 0);
        long64 flops   = 0;
        long64 memrefs = 0;
        Uintah::Solver::ScMult_Add(Dnew, b, D, Q, iter, flops, memrefs);
        new_dw->put(sumlong_vartype(flops), flop_label);
        new_dw->put(sumlong_vartype(memrefs), memref_label);
      }
    }
  }
  //______________________________________________________________________
  void
  setup(const ProcessorGroup*,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse*,
        DataWarehouse* new_dw)
  {
    DataWarehouse* A_dw = new_dw->getOtherDataWarehouse(parent_which_A_dw);
    DataWarehouse* b_dw = new_dw->getOtherDataWarehouse(parent_which_b_dw);
    DataWarehouse* guess_dw =
      new_dw->getOtherDataWarehouse(parent_which_guess_dw);
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      if (cout_doing.active()) {
        cout_doing << "CGSolver::setup on patch " << patch->getID() << std::endl;
      }

      for (int m = 0; m < matls->size(); m++) {
        int matl = matls->get(m);
        typedef typename GridVarType::double_type double_type;
        Patch::VariableBasis basis = Patch::translateTypeToBasis(
          double_type::getTypeDescription()->getType(),
          true);

        IntVector l, h;
        if (params->getSolveOnExtraCells()) {
          l = patch->getExtraLowIndex(basis, IntVector(0, 0, 0));
          h = patch->getExtraHighIndex(basis, IntVector(0, 0, 0));
        } else {
          l = patch->getLowIndex(basis);
          h = patch->getHighIndex(basis);
        }
        CellIterator iter(l, h);

        typename GridVarType::double_type R, Xnew, diagonal;
        new_dw->allocateAndPut(R, R_label, matl, patch);
        new_dw->allocateAndPut(Xnew, X_label, matl, patch);
        new_dw->allocateAndPut(diagonal, diag_label, matl, patch);
        typename GridVarType::const_double_type B;
        typename GridVarType::matrix_type A;
        b_dw->get(B, B_label, matl, patch, Ghost::None, 0);
        A_dw->get(A, A_label, matl, patch, Ghost::None, 0);

        long64 flops   = 0;
        long64 memrefs = 0;
        if (guess_label) {
          typename GridVarType::const_double_type X;
          guess_dw->get(X, guess_label, matl, patch, Around, 1);

          // R = A*X
          IntVector ll(l);
          IntVector hh(h);
          ll -= IntVector(
            patch->getBCType(Patch::xminus) == Patch::Neighbor ? 1 : 0,
            patch->getBCType(Patch::yminus) == Patch::Neighbor ? 1 : 0,
            patch->getBCType(Patch::zminus) == Patch::Neighbor ? 1 : 0);

          hh += IntVector(
            patch->getBCType(Patch::xplus) == Patch::Neighbor ? 1 : 0,
            patch->getBCType(Patch::yplus) == Patch::Neighbor ? 1 : 0,
            patch->getBCType(Patch::zplus) == Patch::Neighbor ? 1 : 0);
          hh -= IntVector(1, 1, 1);

          Uintah::Solver::Mult(R, A, X, iter, ll, hh, flops, memrefs);

          // R = B-R
          Uintah::Solver::Sub(R, B, R, iter, flops, memrefs);
          Xnew.copy(X, iter.begin(), iter.end());
        } else {
          R.copy(B);
          Xnew.initialize(0);
        }

        // D = R/Ap
        typename GridVarType::double_type D;
        new_dw->allocateAndPut(D, D_label, matl, patch);

        Uintah::Solver::InverseDiagonal(diagonal, A, iter, flops, memrefs);
        Uintah::Solver::Mult(D, R, diagonal, iter, flops, memrefs);

        double dnew = Uintah::Solver::Dot(R, D, iter, flops, memrefs);
        new_dw->put(sum_vartype(dnew), d_label);
        new_dw->put(sum_vartype(params->tolerance), tolerance_label);

        // Calculate error term
        double residualNormalization = params->getResidualNormalizationFactor();

        switch (params->norm) {
          case CGSolverParams::Norm::L1: {
            double err = Uintah::Solver::L1(R, iter, flops, memrefs);
            err /= residualNormalization;
            new_dw->put(sum_vartype(err), err_label);
          } break;
          case CGSolverParams::Norm::L2:
            // Nothing...
            break;
          case CGSolverParams::Norm::LInfinity: {
            double err = Uintah::Solver::LInf(R, iter, flops, memrefs);
            err /= residualNormalization;
            new_dw->put(max_vartype(err), err_label);
          } break;
          default:
            break;
        }

        new_dw->put(sumlong_vartype(flops), flop_label);
        new_dw->put(sumlong_vartype(memrefs), memref_label);
      }
    }
  }

  //______________________________________________________________________
  void
  solve(const ProcessorGroup* pg,
        const PatchSubset* patches,
        const MaterialSubset* matls,
        DataWarehouse* old_dw,
        DataWarehouse* new_dw,
        Handle<CGStencil7<GridVarType>>)
  {
    if (cout_doing.active()) {
      cout_doing << "CGSolver::solve" << std::endl;
    }

    Timers::Simple timer;
    timer.start();

    SchedulerP subsched = sched->createSubScheduler();
    DataWarehouse::ScrubMode old_dw_scrubmode =
      old_dw->setScrubbing(DataWarehouse::ScrubNone);
    DataWarehouse::ScrubMode new_dw_scrubmode =
      new_dw->setScrubbing(DataWarehouse::ScrubNone);
    subsched->initialize(3, 1);
    subsched->setParentDWs(old_dw, new_dw);
    subsched->clearMappings();
    subsched->mapDataWarehouse(Task::ParentOldDW, 0);
    subsched->mapDataWarehouse(Task::ParentNewDW, 1);
    subsched->mapDataWarehouse(Task::OldDW, 2);
    subsched->mapDataWarehouse(Task::NewDW, 3);

    GridP grid = level->getGrid();
    IntVector l, h;
    level->findCellIndexRange(l, h);

    int niter = 0;

    subsched->advanceDataWarehouse(grid);

    //__________________________________
    // Schedule the setup
    if (cout_doing.active()) {
      cout_doing << "CGSolver::schedule setup" << std::endl;
    }
    Task* task =
      scinew Task("CGSolver:setup", this, &CGStencil7<GridVarType>::setup);
    task->needs(parent_which_b_dw, B_label, Ghost::None, 0);
    task->needs(parent_which_A_dw, A_label, Ghost::None, 0);

    if (guess_label) {
      task->needs(parent_which_guess_dw, guess_label, Around, 1);
    }

    task->computes(memref_label);
    task->computes(R_label);
    task->computes(X_label);
    task->computes(D_label);
    task->computes(d_label);
    task->computes(tolerance_label);
    task->computes(diag_label);
    if (params->norm != CGSolverParams::Norm::L2) {
      task->computes(err_label);
    }
    task->computes(flop_label);
    subsched->addTask(task, level->eachPatch(), matlset);

    subsched->compile();

    DataWarehouse* subNewDW = subsched->get_dw(3);
    subNewDW->setScrubbing(DataWarehouse::ScrubNone);
    subsched->execute(); // execute CGSolver:setup and the reduction tasks

    //__________________________________
    // At this point the tolerance_label and err_label have
    // been computed by CGSolver:setup and they have been reduced
    double tolerance = 0;
    sum_vartype tol;
    subNewDW->get(tol, tolerance_label);
    tolerance = tol;
    double e  = 0;

    switch (params->norm) {
      case CGSolverParams::Norm::L1:
      case CGSolverParams::Norm::L2: {
        sum_vartype err;
        subNewDW->get(err, err_label);
        e = err;
      } break;
      case CGSolverParams::Norm::LInfinity: {
        max_vartype err;
        subNewDW->get(err, err_label);
        e = err;
      } break;
      default:
        break;
    }
    double err0 = e;
    sumlong_vartype f;
    subNewDW->get(f, flop_label);

    long64 flops = f;
    subNewDW->get(f, memref_label);
    long64 memrefs = f;

    //__________________________________
    if (!(e < params->initial_tolerance)) {
      subsched->initialize(3, 1);
      subsched->setParentDWs(old_dw, new_dw);
      subsched->clearMappings();
      subsched->mapDataWarehouse(Task::ParentOldDW, 0);
      subsched->mapDataWarehouse(Task::ParentNewDW, 1);
      subsched->mapDataWarehouse(Task::OldDW, 2);
      subsched->mapDataWarehouse(Task::NewDW, 3);

      //__________________________________
      // Step 1 - requires A(parent), D(old, 1 ghost) computes aden(new)
      if (cout_doing.active()) {
        cout_doing << "CGSolver::schedule Step 1" << std::endl;
      }
      task =
        scinew Task("CGSolver:step1", this, &CGStencil7<GridVarType>::step1);
      task->needs(parent_which_A_dw, A_label, Ghost::None, 0);
      task->needs(Task::OldDW, D_label, Around, 1);
      task->computes(aden_label);
      task->computes(Q_label);
      task->computes(flop_label);
      task->computes(memref_label);
      subsched->addTask(task, level->eachPatch(), matlset);

      //__________________________________
      // schedule
      // Step 2 - requires d(old), aden(new) D(old), X(old) R(old)  computes X,
      // R, Q, d
      if (cout_doing.active()) {
        cout_doing << "CGSolver::schedule Step 2" << std::endl;
      }
      task =
        scinew Task("CGSolver:step2", this, &CGStencil7<GridVarType>::step2);
      task->needs(Task::OldDW, d_label);
      task->needs(Task::NewDW, aden_label);
      task->needs(Task::OldDW, D_label, Ghost::None, 0);
      task->needs(Task::OldDW, X_label, Ghost::None, 0);
      task->needs(Task::OldDW, R_label, Ghost::None, 0);
      task->needs(Task::OldDW, diag_label, Ghost::None, 0);
      task->computes(X_label);
      task->computes(R_label);
      task->modifies(Q_label);
      task->computes(d_label);
      task->computes(diag_label);
      task->computes(flop_label);
      task->modifies(memref_label);

      if (params->norm != CGSolverParams::Norm::L2) {
        task->computes(err_label);
      }
      subsched->addTask(task, level->eachPatch(), matlset);

      //__________________________________
      // schedule
      // Step 3 - requires D(old), Q(new), d(new), d(old), computes D
      if (cout_doing.active()) {
        cout_doing << "CGSolver::schedule Step 3" << std::endl;
      }
      task =
        scinew Task("CGSolver:step3", this, &CGStencil7<GridVarType>::step3);
      task->needs(Task::OldDW, D_label, Ghost::None, 0);
      task->needs(Task::NewDW, Q_label, Ghost::None, 0);
      task->needs(Task::NewDW, d_label);
      task->needs(Task::OldDW, d_label);
      task->computes(D_label);
      task->computes(flop_label);
      task->modifies(memref_label);
      subsched->addTask(task, level->eachPatch(), matlset);
      subsched->compile();

      //__________________________________
      //  Main iteration
      while (niter < params->maxiterations && !(e < tolerance)) {
        niter++;
        subsched->advanceDataWarehouse(grid);
        DataWarehouse* subOldDW = subsched->get_dw(2);
        DataWarehouse* subNewDW = subsched->get_dw(3);

        subOldDW->setScrubbing(DataWarehouse::ScrubComplete);
        subNewDW->setScrubbing(DataWarehouse::ScrubNonPermanent);

        subsched->execute();

        //__________________________________
        switch (params->norm) {
          case CGSolverParams::Norm::L1:
          case CGSolverParams::Norm::L2: {
            sum_vartype err;
            subNewDW->get(err, err_label);
            e = err;
          } break;
          case CGSolverParams::Norm::LInfinity: {
            max_vartype err;
            subNewDW->get(err, err_label);
            e = err;
          } break;
          default:
            break;
        }
        if (params->criteria == CGSolverParams::Criteria::Relative) {
          e /= err0;
        }
        sumlong_vartype f;
        subNewDW->get(f, flop_label);
        flops += f;
        subNewDW->get(f, memref_label);
        memrefs += f;
      }
    }

    //__________________________________
    //  Pull the solution out of subsched new DW and put it into our X
    if (modifies_x) {
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        for (int m = 0; m < matls->size(); m++) {
          int matl = matls->get(m);
          typedef typename GridVarType::double_type double_type;
          Patch::VariableBasis basis = Patch::translateTypeToBasis(
            double_type::getTypeDescription()->getType(),
            true);

          IntVector l, h;
          if (params->getSolveOnExtraCells()) {
            l = patch->getExtraLowIndex(basis, IntVector(0, 0, 0));
            h = patch->getExtraHighIndex(basis, IntVector(0, 0, 0));
          } else {
            l = patch->getLowIndex(basis);
            h = patch->getHighIndex(basis);
          }
          CellIterator iter(l, h);

          typename GridVarType::double_type Xnew;
          typename GridVarType::const_double_type X;
          new_dw->getModifiable(Xnew, X_label, matl, patch);

          subsched->get_dw(3)->get(X, X_label, matl, patch, Ghost::None, 0);
          Xnew.copy(X, l, h);
        }
      }
    } else {
      new_dw->transferFrom(subsched->get_dw(3), X_label, patches, matls);
    }

    // Restore the scrubbing mode
    old_dw->setScrubbing(old_dw_scrubmode);
    new_dw->setScrubbing(new_dw_scrubmode);

    double dt      = timer().seconds();
    double mflops  = (double(flops) * 1.e-6) / dt;
    double memrate = (double(memrefs) * 1.e-9) / dt;

    if (niter < params->maxiterations) {
      proc0cout << "Solve of " << X_label->getName() << " on level "
                << level->getIndex() << " completed in " << dt << " seconds ("
                << niter << " iterations, " << e << " residual, " << mflops
                << " MFLOPS, " << memrate << " GB/sec)\n";
    } else if (params->getRecomputeTimestepOnFailure()) {
      proc0cout << "CGSolver not converging, requesting smaller time step\n";
      proc0cout << "    niters:   " << niter << "\n"
                << "    residual: " << e << std::endl;

      new_dw->put(bool_or_vartype(true), VarLabel::find(abortTimestep_name));
      new_dw->put(bool_or_vartype(true),
                  VarLabel::find(recomputeTimestep_name));
    } else {
      throw ConvergenceFailure("CGSolve variable: " + X_label->getName(),
                               niter,
                               e,
                               tolerance,
                               __FILE__,
                               __LINE__);
    }
  }
  //______________________________________________________________________
  //
private:
  Scheduler* sched;
  const ProcessorGroup* world;
  const Level* level;
  const MaterialSet* matlset;
  Ghost::GhostType Around;
  const VarLabel* A_label;
  Task::WhichDW which_A_dw, parent_which_A_dw;
  const VarLabel* X_label;
  const VarLabel* B_label;
  Task::WhichDW which_b_dw, parent_which_b_dw;
  const VarLabel* guess_label;
  Task::WhichDW which_guess_dw, parent_which_guess_dw;

  const VarLabel* R_label;
  const VarLabel* D_label;
  const VarLabel* diag_label;
  const VarLabel* Q_label;
  const VarLabel* d_label;
  const VarLabel* err_label;
  const VarLabel* aden_label;
  const VarLabel* flop_label;
  const VarLabel* memref_label;
  const VarLabel* tolerance_label;

  const CGSolverParams* params;
  bool modifies_x;
};
} // namespace Uintah
