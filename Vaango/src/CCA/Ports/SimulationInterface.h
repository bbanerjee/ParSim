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

#ifndef VAANGO_CCA_PORTS_SimulationInterface_H
#define VAANGO_CCA_PORTS_SimulationInterface_H

#include <Core/Parallel/UintahParallelPort.h>

#include <CCA/Ports/SchedulerP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/InfoMapper.h>

namespace Uintah {

/***********************************************************************************
 * CLASS
 *   SimulationInterface
 * GENERAL INFORMATION
 *   SimulationInterface.h
 *     Steven G. Parker
 *     Department of Computer Science
 *     University of Utah
 *     Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
 * KEYWORDS
 *   Simulation_Interface
 * DESCRIPTION
 * WARNING
 ************************************************************************************/

class UintahParallelComponent;

class DataWarehouse;
class Regridder;
class Output;
class VarLabel;

using ValidateFlag = unsigned char;

class SimulationInterface : public UintahParallelPort
{

  // NOTE: No physics simulation component should be a friend class.
  // They should access the data via the DataWarehouse.
  friend class SimulationController;
  friend class AMRSimulationController;
  friend class DataArchiver;

  friend class SchedulerCommon;
  friend class DynamicMPIScheduler;
  friend class MPIScheduler;
  friend class UnifiedScheduler;
  friend class DetailedTasks;

  friend class LoadBalancersCommon;
  friend class DynamicLoadBalancer;
  friend class ParticleLoadBalancer;

  friend class RegridderCommon;

  friend class Switcher;

public:
  SimulationInterface();
  virtual ~SimulationInterface();

  SimulationInterface(const SimulationInterface&) = delete;
  SimulationInterface(SimulationInterface&&)      = delete;
  SimulationInterface&
  operator=(const SimulationInterface&) = delete;
  SimulationInterface&
  operator=(SimulationInterface&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp) = 0;
  virtual void
  getComponents() = 0;
  virtual void
  releaseComponents() = 0;

  virtual Scheduler*
  getScheduler() = 0;
  virtual Regridder*
  getRegridder() = 0;
  virtual Output*
  getOutput() = 0;

  // Top level problem set up called by vaango.
  virtual void
  problemSetup(const ProblemSpecP& prob_spec) = 0;

  virtual void
  problemSetupDeltaT(const ProblemSpecP& prob_spec) = 0;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "") = 0;

  virtual void
  preGridProblemSetup(const ProblemSpecP& params, GridP& grid) = 0;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) = 0;

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& scheduler) = 0;

  // On a restart, schedule an initialization task
  virtual void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& scheduler) = 0;

  // Used by the switcher
  virtual void
  setupForSwitching() = 0;

  // Get the task graph the application wants to execute. Returns an
  // index into the scheduler's list of task graphs.
  virtual void
  setTaskGraphIndex(int index) = 0;
  virtual int
  getTaskGraphIndex() = 0;

  //////////
  // restartInitialize() is called once and only once if and when a simulation
  // is restarted. This allows the simulation component to handle
  // initializations that are necessary when a simulation is restarted.
  //
  virtual void restartInitialize()
  {

  }

  virtual void
  restartInitialize(const ProcessorGroup*, const PatchSubset*, const MaterialSubset*, DataWarehouse*, DataWarehouse*)
  {
  }

  // Schedule the initial switching.
  virtual void
  scheduleSwitchInitialization(const LevelP& level, SchedulerP& sched) = 0;

  virtual void
  scheduleSwitchTest(const LevelP& level, SchedulerP& scheduler) = 0;

  // Schedule the actual time step advancement tasks.
  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP& scheduler) = 0;

  // this is for wrapping up a timestep when it can't be done in
  // scheduleTimeAdvance.
  virtual void
  scheduleFinalizeTimestep(const LevelP& level, SchedulerP& scheduler) = 0;

  // Optionally schedule analysis tasks.
  virtual void
  scheduleAnalysis(const LevelP& level, SchedulerP& scheduler) = 0;

  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP& scheduler) = 0;

  // Reduce the system wide values such as the next delta T.
  virtual void
  scheduleReduceSystemVars(const GridP& grid,
                           const PatchSet* perProcPatchSet,
                           SchedulerP& scheduler) = 0;

  // Schedule the initialization of system values such at the time step.
  virtual void
  scheduleInitializeSystemVars(const GridP& grid,
                               const PatchSet* perProcPatchSet,
                               SchedulerP& scheduler) = 0;

  // Schedule the updating of system values such at the time step.
  virtual void
  scheduleUpdateSystemVars(const GridP& grid,
                           const PatchSet* perProcPatchSet,
                           SchedulerP& scheduler) = 0;

  virtual void
  scheduleRefine(const PatchSet* patches, SchedulerP& scheduler) = 0;

  virtual void
  scheduleRefineInterface(const LevelP& fineLevel,
                          SchedulerP& scheduler,
                          bool needCoarseOld,
                          bool needCoarseNew) = 0;

  virtual void
  scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& scheduler) = 0;

  /// Schedule to mark flags for AMR regridding
  virtual void
  scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched) = 0;

  /// Schedule to mark initial flags for AMR regridding
  virtual void
  scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                               SchedulerP& sched) = 0;

  // Used to get the progress ratio of an AMR regridding subcycle.
  virtual double
  getSubCycleProgress(DataWarehouse* fineNewDW) = 0;

  virtual void
  switchInitialize([[maybe_unused]] const LevelP& level,
                   [[maybe_unused]] SchedulerP&)
  {
  }

  // Redo a timestep if current time advance is not converging.
  // Returned time is the new dt to use.
  virtual void
  recomputeDelT() = 0;
  virtual double
  recomputeDelT(const double delT) = 0;

  // Updates the tiem step and the delta T.
  virtual void
  prepareForNextTimestep() = 0;

  // Ask the component if it needs to be recompiled
  virtual bool
  needRecompile(const GridP& grid) = 0;

  // Instruct component to add a new material
  // virtual void
  // addMaterial(const ProblemSpecP& params, GridP& grid) = 0;

  // Labels for access value in the data warehouse.
  virtual const VarLabel*
  getTimestepLabel() const = 0;
  virtual const VarLabel*
  getSimTimeLabel() const = 0;
  virtual const VarLabel*
  getDelTLabel() const = 0;

  virtual void
  setAMR(bool val) = 0;
  virtual bool
  isAMR() const = 0;

  virtual void
  setLockstepAMR(bool val) = 0;
  virtual bool
  isLockstepAMR() const = 0;

  virtual void
  setDynamicRegridding(bool val) = 0;
  virtual bool
  isDynamicRegridding() const = 0;

  // Boolean for vars changed by the in-situ.
  virtual void
  haveModifiedVars(bool val) = 0;
  virtual bool
  haveModifiedVars() const = 0;

  // For restarting.
  virtual bool
  isRestartTimestep() const = 0;
  virtual void
  setRestartTimestep(bool val) = 0;

  // For regridding.
  virtual bool
  isRegridTimestep() const = 0;
  virtual void
  setRegridTimestep(bool val) = 0;
  virtual int
  getLastRegridTimestep() = 0;
  virtual bool
  wasRegridLastTimestep() const = 0;

  // Some applications can set reduction variables
  virtual void
  addReductionVariable(std::string name,
                       const TypeDescription* varType,
                       bool varActive = false) = 0;
  virtual unsigned int
  numReductionVariable() const = 0;
  virtual void
  activateReductionVariable(std::string name, bool val) = 0;
  virtual bool
  activeReductionVariable(std::string name) const = 0;
  virtual bool
  isBenignReductionVariable(std::string name) = 0;
  virtual void
  setReductionVariable(DataWarehouse* new_dw, std::string name, bool val) = 0;
  virtual void
  setReductionVariable(DataWarehouse* new_dw, std::string name, double val) = 0;
  virtual void
  overrideReductionVariable(DataWarehouse* new_dw,
                            std::string name,
                            bool val) = 0;
  virtual void
  overrideReductionVariable(DataWarehouse* new_dw,
                            std::string name,
                            double val) = 0;

  // Get application specific reduction values all cast to doubles.
  virtual double
  getReductionVariable(std::string name) const = 0;
  virtual double
  getReductionVariable(unsigned int index) const = 0;
  virtual std::string
  getReductionVariableName(unsigned int index) const = 0;
  virtual unsigned int
  getReductionVariableCount(unsigned int index) const = 0;
  virtual bool
  overriddenReductionVariable(unsigned int index) const = 0;

  // Access methods for member classes.
  virtual MaterialManagerP
  getMaterialManagerP() const = 0;

  // Simulation statistics
  enum SimulationStatsEnum
  {
    DummyEnum = 999
  };

  virtual ReductionInfoMapper<SimulationStatsEnum, double>&
  getSimulationStats() = 0;

  virtual void
  resetSimulationStats(double val) = 0;

  virtual void
  reduceSimulationStats(bool allReduce, const ProcessorGroup* myWorld) = 0;

  virtual void
  setDelTOverrideRestart(double val) = 0;
  virtual double
  getDelTOverrideRestart() const = 0;

  virtual void
  setDelTInitialMax(double val) = 0;
  virtual double
  getDelTInitialMax() const = 0;

  virtual void
  setDelTInitialRange(double val) = 0;
  virtual double
  getDelTInitialRange() const = 0;

  virtual void
  setDelTMultiplier(double val) = 0;
  virtual double
  getDelTMultiplier() const = 0;

  virtual void
  setDelTMaxIncrease(double val) = 0;
  virtual double
  getDelTMaxIncrease() const = 0;

  virtual void
  setDelTMin(double val) = 0;
  virtual double
  getDelTMin() const = 0;

  virtual void
  setDelTMax(double val) = 0;
  virtual double
  getDelTMax() const = 0;

  virtual void
  setSimTimeEndAtMax(bool val) = 0;
  virtual bool
  getSimTimeEndAtMax() const = 0;

  virtual void
  setSimTimeMax(double val) = 0;
  virtual double
  getSimTimeMax() const = 0;

  virtual void
  setSimTimeClampToOutput(bool val) = 0;
  virtual bool
  getSimTimeClampToOutput() const = 0;

  virtual void
  setTimestepsMax(int val) = 0;
  virtual int
  getTimestepsMax() const = 0;

  virtual void
  setWallTimeMax(double val) = 0;
  virtual double
  getWallTimeMax() const = 0;

private:
  // Flag for outputting or checkpointing if the next delta is invalid
  virtual void
  setOutputIfInvalidNextDelT(ValidateFlag flag) = 0;
  virtual ValidateFlag
  getOutputIfInvalidNextDelT() const = 0;

  virtual void
  setCheckpointIfInvalidNextDelT(ValidateFlag flag) = 0;
  virtual ValidateFlag
  getCheckpointIfInvalidNextDelT() const = 0;

  // Delta T methods
  virtual void
  setDelT(double delT) = 0;
  virtual double
  getDelT() const = 0;
  virtual void
  setDelTForAllLevels(SchedulerP& scheduler,
                      const GridP& grid,
                      const int totalFine) = 0;

  virtual void
  setNextDelT(double delT, bool restart = false) = 0;
  virtual double
  getNextDelT() const = 0;

  virtual void
  setSimTime(double simTime) = 0;
  virtual double
  getSimTime() const = 0;

  // Returns the integer time step index of the simulation.  All
  // simulations start with a time step number of 0.  This value is
  // incremented by one for each simulation time step processed.

  // The 'set' function should only be called by the
  // SimulationController at the beginning of a simulation.  The
  // 'increment' function is called by the SimulationController at
  // the beginning of each time step.
  virtual void
  setTimestep(int timeStep) = 0;
  virtual void
  incrementTimestep() = 0;
  virtual int
  getTimestep() const = 0;

  virtual bool
  isLastTimestep(double walltime) = 0;
  virtual bool
  maybeLastTimestep(double walltime) const = 0;
};
} // End namespace Uintah

#endif
