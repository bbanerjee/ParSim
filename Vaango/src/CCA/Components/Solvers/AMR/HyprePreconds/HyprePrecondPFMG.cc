/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

//--------------------------------------------------------------------------
// File: HyprePrecondPFMG.cc
// 
// Hypre PFMG (geometric multigrid #2) preconditioner.
//--------------------------------------------------------------------------

#include <CCA/Components/Solvers/AMR/HyprePreconds/HyprePrecondPFMG.h>
#include <CCA/Components/Solvers/AMR/HypreDriver.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <CCA/Components/Solvers/AMR/HypreSolverParams.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Util/DebugStream.h>

using namespace Uintah;
//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "SOLVER_DOING_COUT:+"

static DebugStream cout_doing("SOLVER_DOING_COUT", false);

Priorities
HyprePrecondPFMG::initPriority(void)
  //___________________________________________________________________
  // Function HyprePrecondPFMG::initPriority~
  // Set the Hypre interfaces that PFMG can work with. Currently, only
  // the Struct interface is supported here. The vector of interfaces
  // is sorted by descending priority.
  //___________________________________________________________________
{
  Priorities priority;
  priority.push_back(HypreStruct);
  return priority;
}

void
HyprePrecondPFMG::setup(void)
  //___________________________________________________________________
  // Function HyprePrecondPFMG::setup~
  // Set up the preconditioner object. After this function call, a
  // Hypre solver can use the preconditioner.
  //___________________________________________________________________
{
  const HypreDriver* driver = _solver->getDriver();
  const HypreSolverParams* params = driver->getParams();
  const HypreInterface& interface = driver->getInterface();
  if (interface == HypreStruct) {
    HYPRE_StructSolver precond_solver;
    HYPRE_StructPFMGCreate(driver->getPG()->getComm(), &precond_solver);
    HYPRE_StructPFMGSetMaxIter(precond_solver, 1);
    HYPRE_StructPFMGSetTol(precond_solver, 0.0);
    HYPRE_StructPFMGSetZeroGuess(precond_solver);
    /* weighted Jacobi = 1; red-black GS = 2 */
    HYPRE_StructPFMGSetRelaxType(precond_solver, 1);
    HYPRE_StructPFMGSetNumPreRelax(precond_solver, params->nPre);
    HYPRE_StructPFMGSetNumPostRelax(precond_solver, params->nPost);
    HYPRE_StructPFMGSetSkipRelax(precond_solver, params->skip);
    HYPRE_StructPFMGSetLogging(precond_solver, 0);
    _precond = (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSolve;
    _pcsetup = (HYPRE_PtrToSolverFcn)HYPRE_StructPFMGSetup;
    _precond_solver = (HYPRE_Solver) precond_solver;
  }
}

HyprePrecondPFMG::~HyprePrecondPFMG(void)
  //___________________________________________________________________
  // HyprePrecondPFMG destructor~
  // Destroy the Hypre objects associated with the preconditioner.
  //___________________________________________________________________
{
  const HypreDriver* driver = _solver->getDriver();
  const HypreInterface& interface = driver->getInterface();
  if (interface == HypreStruct) {
    HYPRE_StructPFMGDestroy((HYPRE_StructSolver) _precond_solver);
  }
}
