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

#ifndef __VAANGO_CCA_Ports_SolverInterace_h__
#define __VAANGO_CCA_Ports_SolverInterace_h__

#include <Core/Parallel/UintahParallelPort.h>

#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SchedulerP.h>

#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ReductionVariable.h>
#include <Core/Grid/Variables/Reductions.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <string>

namespace Uintah {

class UintahParallelComponent;
class VarLabel;

class SolverParameters
{
public:
  SolverParameters()          = default;
  virtual ~SolverParameters() = default;

  void
  setSolveOnExtraCells(bool s)
  {
    d_solve_on_extra_cells = s;
  }

  bool
  getSolveOnExtraCells() const
  {
    return d_solve_on_extra_cells;
  }

  void
  setUseStencil4(bool s)
  {
    d_use_stencil_4 = s;
  }

  bool
  getUseStencil4() const
  {
    return d_use_stencil_4;
  }

  void
  setSymmetric(bool s)
  {
    d_symmetric = s;
  }

  bool
  getSymmetric() const
  {
    return d_symmetric;
  }

  void
  setResidualNormalizationFactor(double s)
  {
    d_residual_normalization_factor = s;
  }

  double
  getResidualNormalizationFactor() const
  {
    return d_residual_normalization_factor;
  }

  // If convergence fails call for the timestep to be restarted.
  void
  setRecomputeTimestepOnFailure(bool s)
  {
    d_recomputable_timestep = s;
  }

  bool
  getRecomputeTimestepOnFailure() const
  {
    return d_recomputable_timestep;
  }

  // Used for outputting A, X & B to files
  void
  setOutputFileName(std::string s)
  {
    d_output_file_name = s;
  }

  void
  getOutputFileName(std::vector<std::string>& fname) const
  {
    fname.push_back("A" + d_output_file_name);
    fname.push_back("b" + d_output_file_name);
    fname.push_back("x" + d_output_file_name);
  }

  void
  setSetupFrequency(const int freq)
  {
    d_setup_frequency = freq;
  }

  int
  getSetupFrequency() const
  {
    return d_setup_frequency;
  }
  void
  setUpdateCoefFrequency(const int freq)
  {
    d_update_coef_frequency = freq;
  }

  int
  getUpdateCoefFrequency() const
  {
    return d_update_coef_frequency;
  }

  // oldDW or ParentOldDW
  void
  setWhichOldDW(Task::WhichDW dw)
  {
    ASSERT(((dw == Task::OldDW) || (dw == Task::ParentOldDW)));
    d_which_old_dw = dw;
  }

  Task::WhichDW
  getWhichOldDW() const
  {
    return d_which_old_dw;
  }

private:
  bool d_use_stencil_4{ false };
  bool d_symmetric{ true };
  bool d_solve_on_extra_cells{ false };
  double d_residual_normalization_factor{ 1.0d };
  bool d_recomputable_timestep{ false };
  std::string d_output_file_name;
  int d_setup_frequency{ 1 };       // delete matrix and recreate it and update
                                    // coefficients. Needed if Stencil changes.
  int d_update_coef_frequency{ 1 }; // do not modify matrix stencil/sparsity -
                                    // only change values of coefficients
  Task::WhichDW d_which_old_dw{
    Task::OldDW
  }; // DataWarehouse either old_dw or parent_old_dw
};

class SolverInterface : public UintahParallelPort
{
public:
  SolverInterface() = default;

  virtual ~SolverInterface()
  {
    for (size_t i = 0; i < d_var_labels.size(); ++i) {
      VarLabel::destroy(d_var_labels[i]);
    }
  }

  /*! Disallow copy and move */
  SolverInterface(const SolverInterface&) = delete;
  SolverInterface(SolverInterface&&)      = delete;

  SolverInterface&
  operator=(const SolverInterface&) = delete;
  SolverInterface&
  operator=(SolverInterface&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp) = 0;

  virtual void
  getComponents() = 0;

  virtual void
  releaseComponents() = 0;

  virtual void
  readParameters(ProblemSpecP& params, const std::string& name) = 0;

  virtual SolverParameters*
  getParameters() = 0;

  virtual void
  scheduleInitialize(const LevelP& level,
                     SchedulerP& sched,
                     const MaterialSet* matls) = 0;

  virtual void
  scheduleRestartInitialize(const LevelP& level,
                            SchedulerP& sched,
                            const MaterialSet* matls) = 0;

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
                bool is_first_solve = true) = 0;

  virtual std::string
  getName() = 0;

  /**
   \brief Enforces solvability condition on periodic problems or in domains
   where boundary conditions on the Poisson system are zero Neumann (dp/dn = 0).
   \param bLabel Varlabel of the Poisson system right hand side (RHS). The RHS
   MUST live in the newDW (i.e. be modifiable). The remaining parameters take
   the standard form of other Uintah tasks. \param rkStage: In a multistep
   integration scheme, Uintah is incapable of dealing with multiple reductions
   on the same variable. Hence the need for a stage number (e.g. rkStage) to
   create unique varlabels
   */
  template<typename FieldT>
  void
  scheduleEnforceSolvability(const LevelP& level,
                             SchedulerP& sched,
                             const MaterialSet* matls,
                             const VarLabel* bLabel,
                             const int rkStage);

  /**
   \brief Set a reference pressure in the domain. The user picks a reference
   cell (cRef) and a desired reference value pRef. The pressure in the reference
   cell is pCell = p[cRef]. Then, pCell + dp = pRef or dp = pRef - pCell The
   value dp is the pressure difference that needs to be added to the pressure
   solution at ALL points in the domain to adjust for the reference pressure.
   */
  template<typename FieldT>
  void
  scheduleSetReferenceValue(const LevelP& level,
                            SchedulerP& sched,
                            const MaterialSet* matls,
                            const VarLabel* xLabel,
                            const int rkStage,
                            const IntVector refCell,
                            const double refValue);

private:

  std::vector<VarLabel*> d_var_labels;

  /**
   \brief The user picks a reference cell (cRef) and a
   desired reference value pRef. The pressure in the reference cell is pCell =
   p[cRef]. Then, pCell + dp = pRef or dp = pRef - pCell This task computes dp
   and broadcasts it across all processors. We do this using a reduction
   variable because Uintah doesn't provide us with a nice interface for doing a
   broadcast.
   */
  template<typename FieldT>
  void
  findRefValueDiff(const Uintah::ProcessorGroup*,
                   const Uintah::PatchSubset* patches,
                   const Uintah::MaterialSubset* materials,
                   Uintah::DataWarehouse* old_dw,
                   Uintah::DataWarehouse* new_dw,
                   const VarLabel* xLabel,
                   VarLabel* refValueLabel,
                   const IntVector refCell,
                   const double refValue);

  /**
   \brief Computes the volume integral of the RHS of the Poisson equation: 1/V *
   int(rhs*dV) Since Uintah deals with uniform structured grids, the above
   equation can be simplified, discretely, to: 1/n * sum(rhs) where n is the
   total number of cell sin the domain.
   */
  template<typename FieldT>
  void
  computeRHSIntegral(const Uintah::ProcessorGroup*,
                     const Uintah::PatchSubset* patches,
                     const Uintah::MaterialSubset* materials,
                     Uintah::DataWarehouse* old_dw,
                     Uintah::DataWarehouse* new_dw,
                     const VarLabel* bLabel,
                     VarLabel* rhsIntegralLabel);

  /**
   \brief This task adds dp (see above) to the pressure at all points in the
   domain to reflect the reference value specified by the user/developer.
   */
  template<typename FieldT>
  void
  setRefValue(const Uintah::ProcessorGroup*,
              const Uintah::PatchSubset* patches,
              const Uintah::MaterialSubset* materials,
              Uintah::DataWarehouse* old_dw,
              Uintah::DataWarehouse* new_dw,
              const VarLabel* xLabel,
              VarLabel* refValueLabel);

  /**
   \brief Modifies the RHS of the Poisson equation to satisfy the solvability
   condition on periodic problems.
   */
  template<typename FieldT>
  void
  enforceSolvability(const Uintah::ProcessorGroup*,
                     const Uintah::PatchSubset* patches,
                     const Uintah::MaterialSubset* materials,
                     Uintah::DataWarehouse* old_dw,
                     Uintah::DataWarehouse* new_dw,
                     const VarLabel* bLabel,
                     VarLabel* rhsIntegralLabel);
};
}

#endif //__VAANGO_CCA_Ports_SolverInterace_h__
