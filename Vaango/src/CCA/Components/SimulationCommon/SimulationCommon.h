/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2021-2023 Biswajit Banerjee
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

#ifndef __CCA_COMPONENTS_SIMULATIONCOMMON_H__
#define __CCA_COMPONENTS_SIMULATIONCOMMON_H__

#include <CCA/Ports/SimulationInterface.h>
#include <Core/Parallel/UintahParallelComponent.h>

#include <CCA/Components/SimulationCommon/SimulationReductionVariable.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/SchedulerP.h>
#include <CCA/Ports/SolverInterface.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/OS/Dir.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/Handle.h>

namespace Uintah {

class DataWarehouse;
class Output;
class Regridder;

class SimulationCommon
  : public UintahParallelComponent
  , public SimulationInterface
{

  friend class Switcher;
  friend class PostProcessUda;

public:
  SimulationCommon(const ProcessorGroup* myworld,
                   const MaterialManagerP materialManager);

  virtual ~SimulationCommon();

  SimulationCommon(const SimulationCommon&) = delete;
  SimulationCommon&
  operator=(const SimulationCommon&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp);
  virtual void
  getComponents();
  virtual void
  releaseComponents();

  virtual Scheduler*
  getScheduler()
  {
    return d_scheduler;
  }
  virtual Regridder*
  getRegridder()
  {
    return d_regridder;
  }
  virtual Output*
  getOutput()
  {
    return d_output;
  }

  // Top level problem set up called by sus.
  virtual void
  problemSetup(const ProblemSpecP& prob_spec);

  virtual void
  problemSetupDeltaT(const ProblemSpecP& prob_spec);

  // Top level problem set up called by simulation controller.
  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "") = 0;

  // Called to add missing grid based UPS specs.
  virtual void
  preGridProblemSetup([[maybe_unused]] const ProblemSpecP& params,
                      [[maybe_unused]] GridP& grid){};

  // Used to write parts of the problem spec.
  virtual void
  outputProblemSpec(ProblemSpecP& ps) = 0;

  // Schedule the inital setup of the problem.
  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP& scheduler) = 0;

  // On a restart schedule an initialization task.
  virtual void
  scheduleRestartInitialize(const LevelP& level, SchedulerP& scheduler) = 0;

  // Used by the switcher
  virtual void
  setupForSwitching()
  {
  }

  // Get the task graph the application wants to execute. Returns an
  // index into the scheduler's list of task graphs.
  virtual void
  setTaskGraphIndex(int index)
  {
    d_taskGraphIndex = index;
  }
  virtual int
  getTaskGraphIndex()
  {
    return d_taskGraphIndex;
  }

  // Schedule the inital switching.
  virtual void
  scheduleSwitchInitialization([[maybe_unused]] const LevelP& level,
                               [[maybe_unused]] SchedulerP& sched)
  {
  }

  virtual void
  scheduleSwitchTest([[maybe_unused]] const LevelP& level,
                     [[maybe_unused]] SchedulerP& scheduler){};

  // Schedule the actual time step advencement tasks.
  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP& scheduler);

  // Optionally schedule tasks that can not be done in scheduleTimeAdvance.
  virtual void
  scheduleFinalizeTimestep([[maybe_unused]] const LevelP& level,
                           [[maybe_unused]] SchedulerP& scheduler)
  {
  }

  // Optionally schedule analysis tasks.
  virtual void
  scheduleAnalysis([[maybe_unused]] const LevelP& level,
                   [[maybe_unused]] SchedulerP& scheduler)
  {
  }

  // Optionally schedule a task that determines the next delt T value.
  virtual void
  scheduleComputeStableTimestep(const LevelP& level, SchedulerP& scheduler) = 0;

  // Reduce the system wide values such as the next delta T.
  virtual void
  scheduleReduceSystemVars(const GridP& grid,
                           const PatchSet* perProcPatchSet,
                           SchedulerP& scheduler);

  void
  reduceSystemVars(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw);

  // An optional call for the application to check their reduction vars.
  virtual void
  checkReductionVars([[maybe_unused]] const ProcessorGroup* pg,
                     [[maybe_unused]] const PatchSubset* patches,
                     [[maybe_unused]] const MaterialSubset* matls,
                     [[maybe_unused]] DataWarehouse* old_dw,
                     [[maybe_unused]] DataWarehouse* new_dw){};

  // Schedule the initialization of system values such at the time step.
  virtual void
  scheduleInitializeSystemVars(const GridP& grid,
                               const PatchSet* perProcPatchSet,
                               SchedulerP& scheduler);

  void
  initializeSystemVars(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* /*matls*/,
                       DataWarehouse* /*old_dw*/,
                       DataWarehouse* new_dw);

  // Schedule the updating of system values such at the time step.
  virtual void
  scheduleUpdateSystemVars(const GridP& grid,
                           const PatchSet* perProcPatchSet,
                           SchedulerP& scheduler);

  virtual void
  updateSystemVars(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* /*matls*/,
                   DataWarehouse* /*old_dw*/,
                   DataWarehouse* new_dw);

  // Methods used for scheduling AMR regridding.
  virtual void
  scheduleRefine(const PatchSet* patches, SchedulerP& scheduler);

  virtual void
  scheduleRefineInterface(const LevelP& fineLevel,
                          SchedulerP& scheduler,
                          bool needCoarseOld,
                          bool needCoarseNew);

  virtual void
  scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& scheduler);

  // Schedule to mark flags for AMR regridding
  virtual void
  scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& scheduler);

  // Schedule to mark initial flags for AMR regridding
  virtual void
  scheduleInitialErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched);

  // Used to get the progress ratio of an AMR regridding subcycle.
  virtual double
  getSubCycleProgress(DataWarehouse* fineNewDW);

  // Recompute a time step if current time advance is not
  // converging.  The returned time is the new delta T.
  virtual void
  recomputeDelT();
  virtual double
  recomputeDelT(const double delT);

  // Updates the time step and the delta T.
  virtual void
  prepareForNextTimestep();

  // Asks the application if it needs to be recompiled.
  virtual bool
  needRecompile([[maybe_unused]] const GridP& grid)
  {
    return false;
  };

  // Labels for access value in the data warehouse.
  virtual const VarLabel*
  getTimestepLabel() const
  {
    return d_timeStepLabel;
  }
  virtual const VarLabel*
  getSimTimeLabel() const
  {
    return d_simulationTimeLabel;
  }
  virtual const VarLabel*
  getDelTLabel() const
  {
    return d_delTLabel;
  }

  //////////
  virtual void
  setAMR(bool val)
  {
    d_AMR = val;
  }
  virtual bool
  isAMR() const
  {
    return d_AMR;
  }

  virtual void
  setLockstepAMR(bool val)
  {
    d_lockstepAMR = val;
  }
  virtual bool
  isLockstepAMR() const
  {
    return d_lockstepAMR;
  }

  virtual void
  setDynamicRegridding(bool val)
  {
    d_dynamicRegridding = val;
  }
  virtual bool
  isDynamicRegridding() const
  {
    return d_dynamicRegridding;
  }

  // Boolean for vars chanegd by the in-situ.
  virtual void
  haveModifiedVars(bool val)
  {
    d_haveModifiedVars = val;
  }
  virtual bool
  haveModifiedVars() const
  {
    return d_haveModifiedVars;
  }

  // For restarting.
  virtual bool
  isRestartTimestep() const
  {
    return d_isRestartTimestep;
  }
  virtual void
  setRestartTimestep(bool val)
  {
    d_isRestartTimestep = val;
  }

  // For regridding.
  virtual bool
  isRegridTimestep() const
  {
    return d_isRegridTimestep;
  }
  virtual void
  setRegridTimestep(bool val)
  {
    d_isRegridTimestep = val;

    if (d_isRegridTimestep) {
      d_lastRegridTimestep = d_timeStep;
    }
  }
  virtual int
  getLastRegridTimestep()
  {
    return d_lastRegridTimestep;
  }

  virtual bool
  wasRegridLastTimestep() const
  {
    return ((d_timeStep - d_lastRegridTimestep - 1) == 0);
  }

  // Some applications can set reduction variables
  virtual unsigned int
  numReductionVariable() const
  {
    return d_simReductionVars.size();
  }

  virtual void
  addReductionVariable(std::string name,
                       const TypeDescription* varType,
                       bool varActive = false)
  {
    d_simReductionVars[name] =
      std::make_unique<SimulationReductionVariable>(name, varType, varActive);
  }

  virtual void
  activateReductionVariable(std::string name, bool val)
  {
    d_simReductionVars[name]->setActive(val);
  }
  virtual bool
  activeReductionVariable(std::string name) const
  {
    if (d_simReductionVars.find(name) != d_simReductionVars.end()) {
      return d_simReductionVars.find(name)->second->getActive();
    } else {
      return false;
    }
  }

  virtual bool
  isBenignReductionVariable(std::string name)
  {
    return d_simReductionVars[name]->isBenignValue();
  }
  virtual void
  overrideReductionVariable(DataWarehouse* new_dw, std::string name, bool val)
  {
    d_simReductionVars[name]->setValue(new_dw, val);
  }
  virtual void
  overrideReductionVariable(DataWarehouse* new_dw, std::string name, double val)
  {
    d_simReductionVars[name]->setValue(new_dw, val);
  }
  virtual void
  setReductionVariable(DataWarehouse* new_dw, std::string name, bool val)
  {
    d_simReductionVars[name]->setValue(new_dw, val);
  }
  virtual void
  setReductionVariable(DataWarehouse* new_dw, std::string name, double val)
  {
    d_simReductionVars[name]->setValue(new_dw, val);
  }
  // Get application specific reduction values all cast to doubles.
  virtual double
  getReductionVariable(std::string name) const
  {
    if (d_simReductionVars.find(name) != d_simReductionVars.end()) {
      return d_simReductionVars.at(name)->getValue();
    } else {
      return 0;
    }
  }

  virtual double
  getReductionVariable(unsigned int index) const
  {
    for (const auto& var : d_simReductionVars) {
      if (index == 0) {
        return var.second->getValue();
      } else {
        --index;
      }
    }

    return 0;
  }

  virtual std::string
  getReductionVariableName(unsigned int index) const
  {
    for (const auto& var : d_simReductionVars) {
      if (index == 0) {
        return var.second->getName();
      } else {
        --index;
      }
    }

    return "Unknown";
  }

  virtual unsigned int
  getReductionVariableCount(unsigned int index) const
  {
    for (const auto& var : d_simReductionVars) {
      if (index == 0) {
        return var.second->getCount();
      } else {
        --index;
      }
    }

    return 0;
  }

  virtual bool
  overriddenReductionVariable(unsigned int index) const
  {
    for (const auto& var : d_simReductionVars) {
      if (index == 0) {
        return var.second->overridden();
      } else {
        --index;
      }
    }

    return false;
  }

  // Access methods for member classes.
  virtual MaterialManagerP
  getMaterialManagerP() const
  {
    return d_materialManager;
  }

  virtual ReductionInfoMapper<SimulationStatsEnum, double>&
  getSimulationStats()
  {
    return d_simulation_stats;
  };

  virtual void
  resetSimulationStats(double val)
  {
    d_simulation_stats.reset(val);
  };

  virtual void
  reduceSimulationStats(bool allReduce, const ProcessorGroup* myWorld)
  {
    d_simulation_stats.reduce(allReduce, myWorld);
  };

public:
  virtual void
  setDelTOverrideRestart(double val)
  {
    d_delTOverrideRestart = val;
  }
  virtual double
  getDelTOverrideRestart() const
  {
    return d_delTOverrideRestart;
  }

  virtual void
  setDelTInitialMax(double val)
  {
    d_delTInitialMax = val;
  }
  virtual double
  getDelTInitialMax() const
  {
    return d_delTInitialMax;
  }

  virtual void
  setDelTInitialRange(double val)
  {
    d_delTInitialRange = val;
  }
  virtual double
  getDelTInitialRange() const
  {
    return d_delTInitialRange;
  }

  virtual void
  setDelTMultiplier(double val)
  {
    d_delTMultiplier = val;
  }
  virtual double
  getDelTMultiplier() const
  {
    return d_delTMultiplier;
  }

  virtual void
  setDelTMaxIncrease(double val)
  {
    d_delTMaxIncrease = val;
  }
  virtual double
  getDelTMaxIncrease() const
  {
    return d_delTMaxIncrease;
  }

  virtual void
  setDelTMin(double val)
  {
    d_delTMin = val;
  }
  virtual double
  getDelTMin() const
  {
    return d_delTMin;
  }

  virtual void
  setDelTMax(double val)
  {
    d_delTMax = val;
  }
  virtual double
  getDelTMax() const
  {
    return d_delTMax;
  }

  virtual void
  setSimTimeEndAtMax(bool val)
  {
    d_simTimeEndAtMax = val;
  }
  virtual bool
  getSimTimeEndAtMax() const
  {
    return d_simTimeEndAtMax;
  }

  virtual void
  setSimTimeMax(double val)
  {
    d_simTimeMax = val;
  }
  virtual double
  getSimTimeMax() const
  {
    return d_simTimeMax;
  }

  virtual void
  setSimTimeClampToOutput(bool val)
  {
    d_simTimeClampToOutput = val;
  }
  virtual bool
  getSimTimeClampToOutput() const
  {
    return d_simTimeClampToOutput;
  }

  virtual void
  setTimestepsMax(int val)
  {
    d_timeStepsMax = val;
  }
  virtual int
  getTimestepsMax() const
  {
    return d_timeStepsMax;
  }

  virtual void
  setWallTimeMax(double val)
  {
    d_wallTimeMax = val;
  }
  virtual double
  getWallTimeMax() const
  {
    return d_wallTimeMax;
  }

private:
  // The classes are private because only the top level application
  // should be changing them. This only really matters when there are
  // applications built upon multiple applications. The children
  // applications will not have valid values. They should ALWAYS get
  // the values via the data warehouse.

  // Flag for outputting or checkpointing if the next delta is invalid
  virtual void
  setOutputIfInvalidNextDelT(ValidateFlag flag)
  {
    d_outputIfInvalidNextDelTFlag = flag;
  }
  virtual ValidateFlag
  getOutputIfInvalidNextDelT() const
  {
    return d_outputIfInvalidNextDelTFlag;
  }

  virtual void
  setCheckpointIfInvalidNextDelT(ValidateFlag flag)
  {
    d_checkpointIfInvalidNextDelTFlag = flag;
  }
  virtual ValidateFlag
  getCheckpointIfInvalidNextDelT() const
  {
    return d_checkpointIfInvalidNextDelTFlag;
  }

  //////////
  virtual void
  setDelT(double delT)
  {
    d_delT = delT;
  }
  virtual double
  getDelT() const
  {
    return d_delT;
  }
  virtual void
  setDelTForAllLevels(SchedulerP& scheduler,
                      const GridP& grid,
                      const int totalFine);

  virtual void
  setNextDelT(double delT, bool restart = false);
  virtual double
  getNextDelT() const
  {
    return d_delTNext;
  }
  virtual ValidateFlag
  validateNextDelT(double& delTNext, unsigned int level);

  //////////
  virtual void
  setSimTime(double simTime);
  virtual double
  getSimTime() const
  {
    return d_simTime;
  };

  // Returns the integer time step index of the simulation.  All
  // simulations start with a time step number of 0.  This value is
  // incremented by one before a time step is processed.  The 'set'
  // function should only be called by the SimulationController at the
  // beginning of a simulation.  The 'increment' function is called by
  // the SimulationController at the beginning of each time step.
  virtual void
  setTimestep(int timeStep);
  virtual void
  incrementTimestep();
  virtual int
  getTimestep() const
  {
    return d_timeStep;
  }

  virtual bool
  isLastTimestep(double walltime);
  virtual bool
  maybeLastTimestep(double walltime) const;

protected:
  Scheduler* d_scheduler{ nullptr };
  LoadBalancer* d_loadBalancer{ nullptr };
  SolverInterface* d_solver{ nullptr };
  Regridder* d_regridder{ nullptr };
  Output* d_output{ nullptr };

  // Use a map to store the reduction variables.
  std::map<std::string, std::unique_ptr<SimulationReductionVariable>>
    d_simReductionVars;

  enum VALIDATE_ENUM // unsigned char
  {
    DELTA_T_MAX_INCREASE = 0x01,
    DELTA_T_MIN          = 0x02,
    DELTA_T_MAX          = 0x04,
    DELTA_T_INITIAL_MAX  = 0x08,

    CLAMP_TIME_TO_OUTPUT     = 0x10,
    CLAMP_TIME_TO_CHECKPOINT = 0x20,
    CLAMP_TIME_TO_MAX        = 0x40
  };

private:
  bool d_AMR{ false };
  bool d_lockstepAMR{ false };

  bool d_dynamicRegridding{ false };

  bool d_isRestartTimestep{ false };

  bool d_isRegridTimestep{ false };
  // While it may not have been a "re"-grid, the original grid is
  // created on time step 0.
  int d_lastRegridTimestep{ 0 };

  bool d_haveModifiedVars{ false };

  const VarLabel* d_timeStepLabel;
  const VarLabel* d_simulationTimeLabel;
  const VarLabel* d_delTLabel;

  // Some applications may use multiple task graphs.
  int d_taskGraphIndex{ 0 };

  // The simulation runs to either the maximum number of time steps
  // (timeStepsMax) or the maximum simulation time (simTimeMax), which
  // ever comes first. If the "max_Timestep" is not specified in the .ups
  // file, then it is set to zero.
  double d_delT{ 0.0 };
  double d_delTNext{ 0.0 };

  double d_delTOverrideRestart{ 0 }; // Override the restart delta T value
  double d_delTInitialMax{ 0 };      // Maximum initial delta T
  double d_delTInitialRange{
    0
  }; // Simulation time range for the initial delta T

  double d_delTMin{ 0 };          // Minimum delta T
  double d_delTMax{ 0 };          // Maximum delta T
  double d_delTMultiplier{ 1.0 }; // Multiple for increasing delta T
  double d_delTMaxIncrease{ 0 };  // Maximum delta T increase.

  double d_simTime{ 0.0 };  // Current sim time
  double d_simTimeMax{ 0 }; // Maximum simulation time
  bool d_simTimeEndAtMax{
    false
  }; // End the simulation at exactly this sim time.
  bool d_simTimeClampToOutput{
    false
  }; // Clamp the simulation time to the next output or checkpoint

  int d_timeStep{ 0 };     // Current time step
  int d_timeStepsMax{ 0 }; // Maximum number of time steps to run.

  double d_wallTimeMax{ 0 }; // Maximum wall time.

  ValidateFlag d_outputIfInvalidNextDelTFlag{ 0 };
  ValidateFlag d_checkpointIfInvalidNextDelTFlag{ 0 };

protected:
  MaterialManagerP d_materialManager{ nullptr };

  ReductionInfoMapper<SimulationStatsEnum, double> d_simulation_stats;
};

} // End namespace Uintah

#endif //__CCA_COMPONENTS_SIMULATIONCOMMON_H__
