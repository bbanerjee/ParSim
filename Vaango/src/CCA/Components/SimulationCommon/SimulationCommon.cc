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

#include <CCA/Components/SimulationCommon/SimulationCommon.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SolverInterface.h>
#include <Core/DataArchive/DataArchive.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>

using namespace Uintah;

namespace {

Dout g_deltaT_warn_initial(
    "DeltaTWarnInitial",
    "SimulationCommon",
    "Warn if the next delta T is greater than the initial maximum",
    true);
Dout g_deltaT_warn_increase("DeltaTWarnIncrease",
                            "SimulationCommon",
                            "Warn if the next delta T is increases more than a "
                            "fraction of the previous",
                            true);
Dout g_deltaT_warn_minimum("DeltaTWarnMinimum",
                           "SimulationCommon",
                           "Warn if the next delta T is less than the minimum",
                           true);
Dout g_deltaT_warn_maximum(
    "DeltaTWarnMaximum",
    "SimulationCommon",
    "Warn if the next delta T is greater than the maximum",
    true);
Dout g_deltaT_warn_clamp(
    "DeltaTWarnClamp",
    "SimulationCommon",
    "Warn if the next delta T is clamped for output, checkpoint, or max time",
    true);

Dout g_deltaT_prevalidate(
    "DeltaTPreValidate",
    "SimulationCommon",
    "Before reducing validate the next delta T w/warnings for each rank ",
    false);
Dout g_deltaT_prevalidate_sum("DeltaTPreValidateSum",
                              "SimulationCommon",
                              "Before reducing validate the next delta T "
                              "w/summary warning over all ranks ",
                              false);

}  // namespace

SimulationCommon::SimulationCommon(const ProcessorGroup* myworld,
                                   const MaterialManagerP materialManager)
    : UintahParallelComponent(myworld), d_materialManager(materialManager) {
  // There should only be one MaterialManager. If there is a single
  // application the ComponentFactory will pass in a null pointer
  // which will trigger the MaterialManager to be created.

  // If there are multiple applications the Switcher (which is an
  // application) will create the MaterialManager and then pass that
  // to the child applications.

  // If there are combined applications (aka MPMICE) it will create
  // the MaterialManager and then pass that to the child applications.

  if (d_materialManager == nullptr) {
    d_materialManager = std::make_shared<MaterialManager>();
  }

  //__________________________________
  //  These variables can be modified by an application.

  // Time Step
  d_timeStepLabel =
      VarLabel::create(timeStep_name, timeStep_vartype::getTypeDescription());

  // Simulation Time
  d_simulationTimeLabel =
      VarLabel::create(simTime_name, simTime_vartype::getTypeDescription());

  // delta t
  VarLabel* nonconstDelT =
      VarLabel::create(delT_name, delt_vartype::getTypeDescription());
  nonconstDelT->schedReductionTask(false);
  d_delTLabel = nonconstDelT;

  d_simulation_stats.calculateRankMinimum(true);
  d_simulation_stats.calculateRankStdDev(true);

  // Reduction vars local to the application.

  // An application can request that an output or checkpoint been done
  // immediately.

  // output time step
  d_simReductionVars[outputTimeStep_name] =
      std::make_unique<SimulationReductionVariable>(
          outputTimeStep_name, bool_or_vartype::getTypeDescription());

  // checkpoint time step
  d_simReductionVars[checkpointTimeStep_name] =
      std::make_unique<SimulationReductionVariable>(
          checkpointTimeStep_name, bool_or_vartype::getTypeDescription());

  // An application may adjust the output interval or the checkpoint
  // interval during a simulation.  For example in deflagration ->
  // detonation simulations (Models/HEChem/DDT1.cc

  // output interval
  d_simReductionVars[outputInterval_name] =
      std::make_unique<SimulationReductionVariable>(
          outputInterval_name, min_vartype::getTypeDescription());

  // checkpoint interval
  d_simReductionVars[checkpointInterval_name] =
      std::make_unique<SimulationReductionVariable>(
          checkpointInterval_name, min_vartype::getTypeDescription());

  // An application may also request that the time step be recomputed,
  // aborted or the simulation end early.

  // Recompute the time step
  d_simReductionVars[recomputeTimeStep_name] =
      std::make_unique<SimulationReductionVariable>(
          recomputeTimeStep_name, bool_or_vartype::getTypeDescription());

  // Abort the time step
  d_simReductionVars[abortTimeStep_name] =
      std::make_unique<SimulationReductionVariable>(
          abortTimeStep_name, bool_or_vartype::getTypeDescription());

  // End the simulation
  d_simReductionVars[endSimulation_name] =
      std::make_unique<SimulationReductionVariable>(
          endSimulation_name, bool_or_vartype::getTypeDescription());
}

SimulationCommon::~SimulationCommon() {
  releaseComponents();

  VarLabel::destroy(d_timeStepLabel);
  VarLabel::destroy(d_simulationTimeLabel);
  VarLabel::destroy(d_delTLabel);

  d_simReductionVars.clear();

  // No need to delete the material manager as it is refcounted
  d_materialManager = nullptr;
}

void
SimulationCommon::setComponents(UintahParallelComponent* comp) {
  SimulationCommon* parent = dynamic_cast<SimulationCommon*>(comp);

  attachPort("scheduler", parent->d_scheduler);
  attachPort("load balancer", parent->d_loadBalancer);
  attachPort("solver", parent->d_solver);
  attachPort("regridder", parent->d_regridder);
  attachPort("output", parent->d_output);

  getComponents();
}

void
SimulationCommon::getComponents() {
  d_scheduler = dynamic_cast<Scheduler*>(getPort("scheduler"));

  if (!d_scheduler) {
    throw InternalError(
        "dynamic_cast of 'd_scheduler' failed!", __FILE__, __LINE__);
  }

  d_loadBalancer = dynamic_cast<LoadBalancer*>(getPort("load balancer"));

  if (!d_loadBalancer) {
    throw InternalError(
        "dynamic_cast of 'd_loadBalancer' failed!", __FILE__, __LINE__);
  }

  d_solver = dynamic_cast<SolverInterface*>(getPort("solver"));

  if (!d_solver) {
    throw InternalError(
        "dynamic_cast of 'd_solver' failed!", __FILE__, __LINE__);
  }

  d_regridder = dynamic_cast<Regridder*>(getPort("regridder"));

  if (isDynamicRegridding() && !d_regridder) {
    throw InternalError(
        "dynamic_cast of 'd_regridder' failed!", __FILE__, __LINE__);
  }

  d_output = dynamic_cast<Output*>(getPort("output"));

  if (!d_output) {
    throw InternalError(
        "dynamic_cast of 'd_output' failed!", __FILE__, __LINE__);
  }
}

void
SimulationCommon::releaseComponents() {
  releasePort("scheduler");
  releasePort("load balancer");
  releasePort("solver");
  releasePort("regridder");
  releasePort("output");

  d_scheduler    = nullptr;
  d_loadBalancer = nullptr;
  d_solver       = nullptr;
  d_regridder    = nullptr;
  d_output       = nullptr;
}

void
SimulationCommon::problemSetup(const ProblemSpecP& prob_spec) {
  // Check for an AMR attribute with the grid.
  ProblemSpecP grid_ps = prob_spec->findBlock("Grid");

  if (grid_ps) {
    grid_ps->getAttribute("doAMR", d_AMR);
    d_dynamicRegridding = d_AMR;
  }

  // If the AMR block is defined default to turning AMR on.
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");

  if (amr_ps) {
    d_AMR = true;

    std::string type;
    amr_ps->getAttribute("type", type);

    d_dynamicRegridding = (type.empty() || type == std::string("Dynamic"));

    amr_ps->get("useLockStep", d_lockstepAMR);
  }

  // Get the common time stepping specs
  ProblemSpecP time_ps = prob_spec->findBlock("Time");

  if (!time_ps) {
    throw ProblemSetupException(
        "ERROR SimulationTime \n"
        "Can not find the <Time> block.",
        __FILE__,
        __LINE__);
  }

  // Sim time limits

  // Initial simulation time - will be written to the data warehouse
  // when SimulationController::timeStateSetup() is called.
  time_ps->require("initTime", d_simTime);

  // Maximum simulation time
  time_ps->require("maxTime", d_simTimeMax);

  // End the simulation at exactly the maximum simulation time
  if (!time_ps->get("end_at_max_time_exactly", d_simTimeEndAtMax)) {
    d_simTimeEndAtMax = false;
  }

  // Output time
  if (!time_ps->get("clamp_time_to_output", d_simTimeClampToOutput)) {
    d_simTimeClampToOutput = false;
  }

  // Time step limit
  if (!time_ps->get("max_Timesteps", d_timeStepsMax)) {
    d_timeStepsMax = 0;
  }

  // Wall time limit
  if (!time_ps->get("max_wall_time", d_wallTimeMax)) {
    d_wallTimeMax = 0;
  }

  // Delta T values - also used by the Switcher.
  problemSetupDeltaT(prob_spec);
}

void
SimulationCommon::problemSetupDeltaT(const ProblemSpecP& prob_spec) {
  ProblemSpecP time_ps = prob_spec->findBlock("Time");

  if (!time_ps) {
    throw ProblemSetupException(
        "ERROR SimulationTime \n"
        "Can not find the <Time> block.",
        __FILE__,
        __LINE__);
  }

  // Delta T limits
  ProblemSpecP tmp_ps;
  std::string flag;

  d_outputIfInvalidNextDelTFlag     = 0;
  d_checkpointIfInvalidNextDelTFlag = 0;

  // When restarting use this delta T value
  if (!time_ps->get("override_restart_delt", d_delTOverrideRestart)) {
    d_delTOverrideRestart = 0.0;
  }

  // Multiply the next delta T value by this value
  time_ps->require("timestep_multiplier", d_delTMultiplier);

  // The maximum delta T can increase as a percent over the previous value
  if (!time_ps->get("max_delt_increase", d_delTMaxIncrease)) {
    d_delTMaxIncrease = 0;
  } else  // Can optionally output and/or checkpoint if exceeded
  {
    tmp_ps = time_ps->findBlock("max_delt_increase");
    tmp_ps->getAttribute("output", flag);
    if (flag == std::string("true")) {
      d_outputIfInvalidNextDelTFlag |= DELTA_T_MAX_INCREASE;
    }

    tmp_ps->getAttribute("checkpoint", flag);
    if (!flag.empty() && flag == std::string("true")) {
      d_checkpointIfInvalidNextDelTFlag |= DELTA_T_MAX_INCREASE;
    }
  }

  // The maximum delta T for the initial simulation time from
  // initial_delt_range
  if (!time_ps->get("delt_init", d_delTInitialMax)) {
    d_delTInitialMax = 0;
  } else  // Can optionally output and/or checkpoint if exceeded
  {
    tmp_ps = time_ps->findBlock("delt_init");
    tmp_ps->getAttribute("output", flag);
    if (flag == std::string("true")) {
      d_outputIfInvalidNextDelTFlag |= DELTA_T_INITIAL_MAX;
    }

    tmp_ps->getAttribute("checkpoint", flag);
    if (!flag.empty() && flag == std::string("true")) {
      d_checkpointIfInvalidNextDelTFlag |= DELTA_T_INITIAL_MAX;
    }
  }

  // The maximum simulation time which to enforce the delt_init
  if (!time_ps->get("initial_delt_range", d_delTInitialRange)) {
    d_delTInitialRange = 0;
  }

  // The minimum delta T value
  time_ps->require("delt_min", d_delTMin);
  // Can optionally output and/or checkpoint if exceeded
  tmp_ps = time_ps->findBlock("delt_min");
  tmp_ps->getAttribute("output", flag);
  if (flag == std::string("true")) {
    d_outputIfInvalidNextDelTFlag |= DELTA_T_MIN;
  }

  tmp_ps->getAttribute("checkpoint", flag);
  if (flag == std::string("true")) {
    d_checkpointIfInvalidNextDelTFlag |= DELTA_T_MIN;
  }

  // The maximum delta T value
  time_ps->require("delt_max", d_delTMax);
  // Can optionally output and/or checkpoint if exceeded
  tmp_ps = time_ps->findBlock("delt_max");
  tmp_ps->getAttribute("output", flag);
  if (flag == std::string("true")) {
    d_outputIfInvalidNextDelTFlag |= DELTA_T_MAX;
  }

  tmp_ps->getAttribute("checkpoint", flag);
  if (!flag.empty() && flag == std::string("true")) {
    d_checkpointIfInvalidNextDelTFlag |= DELTA_T_MAX;
  }
}

void
SimulationCommon::scheduleTimeAdvance(const LevelP&, SchedulerP&) {
  throw InternalError(
      "scheduleTimeAdvance is not implemented for this application",
      __FILE__,
      __LINE__);
}

void
SimulationCommon::scheduleReduceSystemVars(const GridP& grid,
                                           const PatchSet* perProcPatchSet,
                                           SchedulerP& scheduler) {
  // Reduce the system vars which are on a per patch basis to a per
  // rank basis.
  Task* task = scinew Task("SimulationCommon::reduceSystemVars",
                           this,
                           &SimulationCommon::reduceSystemVars);

  task->setType(Task::OncePerProc);
  task->usesMPI(true);

  // coarsen delT task requires that delT is computed on every level,
  // even if no tasks are run on that level.  I think this is a bug.
  // --Todd
  // TODO: Look into this - APH 02/27/17

  // Coarsen delT computes the global delT variable
  task->computes(d_delTLabel);

  for (int i = 0; i < grid->numLevels(); i++) {
    task->requires(Task::NewDW, d_delTLabel, grid->getLevel(i).get_rep());
  }

  // These are the application reduction variables. An application may
  // also request that the time step be recomputed, aborted, and/or the
  // simulation end early.

  // Check for a task computing the reduction variable, if found add
  // in a requires and activate the variable it will be tested.
  for (auto& var : d_simReductionVars) {
    const VarLabel* label = var.second->getLabel();

    if (scheduler->getComputedVars().find(label) !=
        scheduler->getComputedVars().end()) {
      activateReductionVariable(var.first, true);

      task->requires(Task::NewDW, label);
      task->computes(label);
    }
  }

  // These two reduction vars may be set by the application via a
  // compute in which case a requires is needed (done above). Or if
  // the flag is set by SimulationCommon, no requires is needed but
  // the reduction var needs to be active for the reduction and
  // subsequent test.
  if (d_outputIfInvalidNextDelTFlag) {
    activateReductionVariable(outputTimeStep_name, true);
  }

  if (d_checkpointIfInvalidNextDelTFlag) {
    activateReductionVariable(checkpointTimeStep_name, true);
  }

  // The above three tasks are on a per proc basis any rank can make
  // the request because it is a either benign or a set value.
  scheduler->addTask(task, perProcPatchSet, d_materialManager->allMaterials());
}

//______________________________________________________________________
//
void
SimulationCommon::reduceSystemVars(const ProcessorGroup* pg,
                                   const PatchSubset* patches,
                                   const MaterialSubset* matls,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw) {
  ValidateFlag validDelT = 0;

  // The goal of this task is to line up the delT across all levels.
  // If the coarse delT already exists (the one without an associated
  // level), then the application is not doing AMR.
  Patch* patch = nullptr;

  if (patches->size() != 0 && !new_dw->exists(d_delTLabel, -1, patch)) {
    // Start with the time step multiplier.
    double multiplier = d_delTMultiplier;

    const GridP grid = patches->get(0)->getLevel()->getGrid();

    for (int l = 0; l < grid->numLevels(); l++) {
      const LevelP level = grid->getLevel(l);

      if (l > 0 && !d_lockstepAMR) {
        multiplier *= level->getRefinementRatioMaxDim();
      }

      if (new_dw->exists(d_delTLabel, -1, *level->patchesBegin())) {
        delt_vartype delTvar;
        new_dw->get(delTvar, d_delTLabel, level.get_rep());

        // Adjust the local next delT by the multiplier
        d_delTNext = delTvar * multiplier;

        // Valiadate before the reduction. This assures that there will
        // not be any possible round off error for the next delta T.
        if (g_deltaT_prevalidate || g_deltaT_prevalidate_sum) {
          validDelT = validateNextDelT(d_delTNext, l);
        }

        new_dw->put(delt_vartype(d_delTNext), d_delTLabel);
      }

      // What should happen if there is no delta T???
    }
  }

  if (d_myworld->nRanks() > 1) {
    new_dw->reduceMPI(d_delTLabel, 0, 0, -1);
  }

  // Get the reduced next delta T
  delt_vartype delTvar;
  new_dw->get(delTvar, d_delTLabel);
  d_delTNext = delTvar;

  // Validate after the reduction. NOTE: Because each rank will
  // independently modify delta T the resulting values may be
  // different due to round off.
  if (!g_deltaT_prevalidate && !g_deltaT_prevalidate_sum) {
    // Validate and put the value into the warehouse if it changed.
    if ((validDelT = validateNextDelT(d_delTNext, -1))) {
      new_dw->override(delt_vartype(d_delTNext), d_delTLabel);
    }
  }

  // If delta T has been changed and if requested, for that change
  // output or checkpoint. Must be done before the reduction call.
  if (validDelT & d_outputIfInvalidNextDelTFlag) {
    setReductionVariable(new_dw, outputTimeStep_name, true);
  }

  if (validDelT & d_checkpointIfInvalidNextDelTFlag) {
    setReductionVariable(new_dw, checkpointTimeStep_name, true);
  }

  // Reduce the application specific reduction variables. If no value
  // was computed on an MPI rank, a benign value will be set. If the
  // reduction result is also a benign value, that means no MPI rank
  // wants to change the value and it will be ignored.
  for (auto& var : d_simReductionVars) {
    var.second->reduce(new_dw);
  }

  // When checking a reduction var, if it is not a benign value then
  // it was set at some point by at least one rank. Which is the only
  // time the value should be use.

  // Specific handling for reduction vars that need the grid.
  if (patches->size() != 0) {
    const GridP grid = patches->get(0)->getLevel()->getGrid();

    if (!isBenignReductionVariable(outputTimeStep_name)) {
      d_output->setOutputTimeStep(true, grid);
    }

    if (!isBenignReductionVariable(checkpointTimeStep_name)) {
      d_output->setCheckpointTimeStep(true, grid);
    }
  }

  // Specific handling for other reduction vars.
  if (!isBenignReductionVariable(outputInterval_name)) {
    d_output->setOutputInterval(getReductionVariable(outputInterval_name));
  }

  if (!isBenignReductionVariable(checkpointInterval_name)) {
    d_output->setCheckpointInterval(
        getReductionVariable(checkpointInterval_name));
  }

  checkReductionVars(pg, patches, matls, old_dw, new_dw);

}  // end reduceSysVar()

//______________________________________________________________________
//
void
SimulationCommon::scheduleInitializeSystemVars(const GridP& grid,
                                               const PatchSet* perProcPatchSet,
                                               SchedulerP& scheduler) {
  // Initialize the system vars which are on a per rank basis.
  Task* task = scinew Task("SimulationCommon::initializeSystemVars",
                           this,
                           &SimulationCommon::initializeSystemVars);

  task->setType(Task::OncePerProc);

  task->computes(d_timeStepLabel);
  task->computes(d_simulationTimeLabel);

  // treatAsOld copyData noScrub notCopyData noCheckpoint
  scheduler->overrideVariableBehavior(
      d_timeStepLabel->getName(), false, false, false, true, true);
  scheduler->overrideVariableBehavior(
      d_simulationTimeLabel->getName(), false, false, false, true, true);

  scheduler->addTask(task, perProcPatchSet, d_materialManager->allMaterials());
}

//______________________________________________________________________
//
void
SimulationCommon::initializeSystemVars(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* /*matls*/,
                                       DataWarehouse* /*old_dw*/,
                                       DataWarehouse* new_dw) {
  // Initialize the time step.
  new_dw->put(timeStep_vartype(d_timeStep), d_timeStepLabel);

  // Initialize the simulation time.
  new_dw->put(simTime_vartype(d_simTime), d_simulationTimeLabel);
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleUpdateSystemVars(const GridP& grid,
                                           const PatchSet* perProcPatchSet,
                                           SchedulerP& scheduler) {
  // Update the system vars which are on a per rank basis.
  Task* task = scinew Task("SimulationCommon::updateSystemVars",
                           this,
                           &SimulationCommon::updateSystemVars);

  task->setType(Task::OncePerProc);

  task->computes(d_timeStepLabel);
  task->computes(d_simulationTimeLabel);

  // treatAsOld copyData noScrub notCopyData noCheckpoint
  scheduler->overrideVariableBehavior(
      d_timeStepLabel->getName(), false, false, false, true, true);
  scheduler->overrideVariableBehavior(
      d_simulationTimeLabel->getName(), false, false, false, true, true);

  scheduler->addTask(task, perProcPatchSet, d_materialManager->allMaterials());
}

//______________________________________________________________________
//
void
SimulationCommon::updateSystemVars(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* /*matls*/,
                                   DataWarehouse* /*old_dw*/,
                                   DataWarehouse* new_dw) {
  // If recomputing a time step do not update the time step or the simulation
  // time.
  if (!getReductionVariable(recomputeTimeStep_name)) {
    // Store the time step so it can be incremented at the top of the
    // time step where it is over written.
    new_dw->put(timeStep_vartype(d_timeStep), d_timeStepLabel);

    // Update the simulation time.
    d_simTime += d_delT;

    // Question - before putting the value into the warehouse should
    // it be broadcasted to assure it is sync'd acrosss all ranks?
    // Uintah::MPI::Bcast( &d_simTime, 1, MPI_DOUBLE, 0, d_myworld->getComm() );

    new_dw->put(simTime_vartype(d_simTime), d_simulationTimeLabel);
  }
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleRefine(const PatchSet*, SchedulerP&) {
  throw InternalError("scheduleRefine not implemented for this application\n",
                      __FILE__,
                      __LINE__);
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleRefineInterface(const LevelP&,
                                          SchedulerP&,
                                          bool,
                                          bool) {
  throw InternalError(
      "scheduleRefineInterface is not implemented for this application\n",
      __FILE__,
      __LINE__);
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleCoarsen(const LevelP&, SchedulerP&) {
  throw InternalError(
      "scheduleCoarsen is not implemented for this application\n",
      __FILE__,
      __LINE__);
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleErrorEstimate(const LevelP&, SchedulerP&) {
  throw InternalError(
      "scheduleErrorEstimate is not implemented for this application",
      __FILE__,
      __LINE__);
}

//______________________________________________________________________
//
void
SimulationCommon::scheduleInitialErrorEstimate(const LevelP& /*coarseLevel*/,
                                               SchedulerP& /*sched*/) {
  throw InternalError(
      "scheduleInitialErrorEstimate is not implemented for this application",
      __FILE__,
      __LINE__);
}

//______________________________________________________________________
//
double
SimulationCommon::getSubCycleProgress(DataWarehouse* fineDW) {
  // DWs are always created in order of time.
  int fineID      = fineDW->getID();
  int coarseNewID = fineDW->getOtherDataWarehouse(Task::CoarseNewDW)->getID();

  // Need to do this check, on init timestep, old DW is nullptr, and
  // getOtherDW will throw exception.
  if (fineID == coarseNewID) {
    return 1.0;
  }

  int coarseOldID = fineDW->getOtherDataWarehouse(Task::CoarseOldDW)->getID();

  return ((double)fineID - coarseOldID) / (coarseNewID - coarseOldID);
}

//______________________________________________________________________
//
void
SimulationCommon::recomputeDelT() {
  // Get the new delT from the actual application.

  // Call the actual application method which if defined overrides the
  // virtual default method for the delta T
  double new_delT = recomputeDelT(d_delT);

  proc0cout << "WARNING Recomputng time step " << d_timeStep << " "
            << "and sim time " << d_simTime << " "
            << ", changing delT from " << d_delT << " to " << new_delT
            << std::endl;

  // Bulletproofing
  if (new_delT < d_delTMin || new_delT <= 0) {
    std::ostringstream warn;
    warn << "The new delT (" << new_delT << ") is either less than "
         << "minDelT (" << d_delTMin << ") or equal to 0";
    throw InternalError(warn.str(), __FILE__, __LINE__);
  }

  // When recomputing the delT, rank 0 determines the value and
  // sends it to all other ranks.
  Uintah::MPI::Bcast(&new_delT, 1, MPI_DOUBLE, 0, d_myworld->getComm());

  d_delT = new_delT;
}

//______________________________________________________________________
//
double
SimulationCommon::recomputeDelT(const double delT) {
  throw InternalError("recomputeDelT is not implemented for this application",
                      __FILE__,
                      __LINE__);
}

//______________________________________________________________________
//
void
SimulationCommon::prepareForNextTimeStep() {
  // Increment (by one) the current time step number so components know
  // what time step they are on and get the delta T that will be used.
  incrementTimeStep();

  // Get the delta that will be used for the time step.
  delt_vartype delt_var;
  d_scheduler->getLastDW()->get(delt_var, d_delTLabel);
  d_delT = delt_var;

  // Clear the time step based reduction variables.
  for (auto& var : d_simReductionVars) {
    var.second->reset();
  }
}

//______________________________________________________________________
//
void
SimulationCommon::setDelTForAllLevels(SchedulerP& scheduler,
                                      const GridP& grid,
                                      const int totalFine) {
  // Adjust the delT for each level and store it in all applicable dws.
  double delT_fine = d_delT;
  int skip         = totalFine;

  for (int i = 0; i < grid->numLevels(); ++i) {
    const Level* level = grid->getLevel(i).get_rep();

    if (isAMR() && i != 0 && !isLockstepAMR()) {
      int trr = level->getRefinementRatioMaxDim();
      delT_fine /= trr;
      skip /= trr;
    }

    for (int idw = 0; idw < totalFine; idw += skip) {
      DataWarehouse* dw = scheduler->get_dw(idw);
      dw->override(delt_vartype(delT_fine), d_delTLabel, level);

      // In a similar fashion write the time step and simulation time
      // to all DWs when running AMR grids.
      dw->override(timeStep_vartype(d_timeStep), d_timeStepLabel);

      dw->override(simTime_vartype(d_simTime), d_simulationTimeLabel);
    }
  }

  // Override for the global level as well (only matters on dw 0)
  DataWarehouse* oldDW = scheduler->get_dw(0);
  oldDW->override(delt_vartype(delT_fine), d_delTLabel);
}

//______________________________________________________________________
//
// This method is called only at restart - see
// SimulationController::timeStateSetup()  or by the in situ - see
// visit_DeltaTVariableCallback().

void
SimulationCommon::setNextDelT(double delT, bool restart) {
  // Restart - Check to see if the user has set a restart delT.
  if (restart && d_delTOverrideRestart) {
    proc0cout << "Overriding restart delT " << d_delT << " with "
              << d_delTOverrideRestart << "\n";

    d_delTNext = d_delTOverrideRestart;

    d_scheduler->getLastDW()->override(delt_vartype(d_delTNext), d_delTLabel);
  }

  // Restart - Otherwise get the next delta T from the archive.
  else if (restart && d_scheduler->getLastDW()->exists(d_delTLabel)) {
    delt_vartype delt_var;
    d_scheduler->getLastDW()->get(delt_var, d_delTLabel);
    d_delTNext = delt_var;
  }

  // All else fails use the delta T passed in. If from a restart it
  // would be the value used at the last time step.
  else {
    d_delTNext = delT;
    d_scheduler->getLastDW()->override(delt_vartype(d_delTNext), d_delTLabel);
  }
}

//______________________________________________________________________
//
ValidateFlag
SimulationCommon::validateNextDelT(double& delTNext, unsigned int level) {
  // NOTE: This check is performed BEFORE the simulation time is
  // updated. As such, being that the time step has completed, the
  // actual simulation time is the current simulation time plus the
  // current delta T.

  // The invalid flag is a bitwise XOR for the local rank. Each bit
  // represents which threshold was exceeded. It is reduced on all
  // ranks if the pre-validating min flag is true (and the
  // pre-validating maxflag is false).
  std::ostringstream header, message;

  header << "WARNING ";

  // For the pre-validate report the rank and level.
  if (g_deltaT_prevalidate)
    header << "Rank-" << d_myworld->myRank() << " for level " << level << " ";

  header << "at time step " << d_timeStep << " and sim time "
         << d_simTime + d_delT << " : ";

  ValidateFlag invalid = 0;

  // Check to see if the next delT was increased too much over the
  // current delT
  double delt_tmp = (1.0 + d_delTMaxIncrease) * d_delT;

  if (d_delTMaxIncrease > 0 && delt_tmp > 0 && delTNext > delt_tmp) {
    invalid |= DELTA_T_MAX_INCREASE;

    if (g_deltaT_warn_increase) {
      if (!message.str().empty()) message << std::endl;

      message << header.str() << "lowering the next delT from " << delTNext
              << " to the maximum: " << delt_tmp
              << " (maximum increase permitted is " << 1.0 + d_delTMaxIncrease
              << "x the previous)";
    }

    delTNext = delt_tmp;
  }

  // Check to see if the next delT is below the minimum delt
  if (d_delTMin > 0 && delTNext < d_delTMin) {
    invalid |= DELTA_T_MIN;

    if (g_deltaT_warn_minimum) {
      if (!message.str().empty()) {
        message << std::endl;
      }

      message << header.str() << "raising the next delT from " << delTNext
              << " to the minimum: " << d_delTMin;
    }

    delTNext = d_delTMin;
  }

  // Check to see if the next delT exceeds the maximum delt
  if (d_delTMax > 0 && delTNext > d_delTMax) {
    invalid |= DELTA_T_MAX;

    if (g_deltaT_warn_maximum) {
      if (!message.str().empty()) {
        message << std::endl;
      }

      message << header.str() << "lowering the next delT from " << delTNext
              << " to the maximum: " << d_delTMax;
    }

    delTNext = d_delTMax;
  }

  // Check to see if the next delT exceeds the maximum initial delt
  // This check shoud be last because it is for the initial time steps.
  if (d_delTInitialMax > 0 && d_simTime + d_delT <= d_delTInitialRange &&
      delTNext > d_delTInitialMax) {
    invalid |= DELTA_T_INITIAL_MAX;

    if (g_deltaT_warn_initial) {
      if (!message.str().empty()) {
        message << std::endl;
      }

      message << header.str() << "for the initial time up to "
              << d_delTInitialRange << " lowering the next delT from "
              << delTNext << " to the maximum: " << d_delTInitialMax;
    }

    delTNext = d_delTInitialMax;
  }

  // Perform last so to possibly not cause other checks and warnings
  // as these checks may reduce the delta T to be smaller than the
  // minimums. Unless requested, no warning is issued as there is no
  // problem with the next delta T.

  // Adjust the next delT to clamp the simulation time to the requested
  // output and/or checkpoint times.
  if (d_simTimeClampToOutput) {
    // Adjust the next delta T to clamp the simulation time to the
    // output time.
    double nextOutput = d_output->getNextOutputTime();
    if (!d_output->isOutputTimeStep() && nextOutput != 0 &&
        d_simTime + d_delT + delTNext > nextOutput) {
      invalid |= CLAMP_TIME_TO_OUTPUT;

      if (g_deltaT_warn_clamp) {
        if (!message.str().empty()) message << std::endl;

        message << header.str() << "lowering the next delT from " << delTNext
                << " to " << nextOutput - (d_simTime + d_delT)
                << " to line up with the next output time: " << nextOutput;
      }

      delTNext = nextOutput - (d_simTime + d_delT);
    }

    // Adjust the next delta T to clamp the simulation time to the
    // checkpoint time.
    double nextCheckpoint = d_output->getNextCheckpointTime();
    if (!d_output->isCheckpointTimeStep() && nextCheckpoint != 0 &&
        d_simTime + d_delT + delTNext > nextCheckpoint) {
      invalid |= CLAMP_TIME_TO_CHECKPOINT;

      if (g_deltaT_warn_clamp) {
        if (!message.str().empty()) {
          message << std::endl;
        }

        message << header.str() << "lowering the next delT from " << delTNext
                << " to " << nextCheckpoint - (d_simTime + d_delT)
                << " to line up with the next checkpoint time: "
                << nextCheckpoint;
      }

      delTNext = nextCheckpoint - (d_simTime + d_delT);
    }
  }

  // Adjust delta T so to end at the max simulation time.
  if (d_simTimeEndAtMax && d_simTime + d_delT + delTNext > d_simTimeMax) {
    invalid |= CLAMP_TIME_TO_MAX;

    if (g_deltaT_warn_clamp) {
      if (!message.str().empty()) message << std::endl;

      message << header.str() << "lowering the next delT from " << delTNext
              << " to " << d_simTimeMax - (d_simTime + d_delT)
              << " to line up with the maximum simulation time of "
              << d_simTimeMax;
    }

    delTNext = d_simTimeMax - (d_simTime + d_delT);
  }

  // Check for a message which indicates that delta T was adjusted and
  // the user wants to be warned (see the g_deltaT_major_warnings and
  // g_deltaT_minor_warnings flags).
  if (!message.str().empty()) {
    // The pre-validate flag is true but not the pre-validate sum flag
    // report for all ranks or if no pre-validating flags are set
    // then a post-validate (default) so reprort for rank 0 only.
    if ((g_deltaT_prevalidate && !g_deltaT_prevalidate_sum) ||
        (!g_deltaT_prevalidate && !g_deltaT_prevalidate_sum &&
         d_myworld->myRank() == 0)) {
      DOUT(true, message.str());
    }
  }

  // Report if pre-validating sum flag is set and only for level zero
  // as it is a summary.
  if (g_deltaT_prevalidate_sum && level == 0) {
    // Gather all of the bits where the threshold was exceeded.
    ValidateFlag invalidAll;

    Uintah::MPI::Reduce(&invalid,
                        &invalidAll,
                        1,
                        MPI_UNSIGNED_CHAR,
                        MPI_BOR,
                        0,
                        d_myworld->getComm());

    // Only report the summary on rank 0. One line for each instance
    // where the threshold was exceeded.
    if (d_myworld->myRank() == 0) {
      std::ostringstream header;
      header << "WARNING "
             << "at time step " << d_timeStep << " "
             << "and sim time " << d_simTime + d_delT << " : ";

      std::ostringstream message;

      // Report the warnings
      if (g_deltaT_warn_increase && (invalidAll & DELTA_T_MAX_INCREASE)) {
        if (!message.str().empty()) {
          message << std::endl;
        }

        message << header.str()
                << "for one or more ranks the next delta T was lowered."
                << " The maximum increase permitted is "
                << 1.0 + d_delTMaxIncrease << "x the previous";
      }

      if (g_deltaT_warn_minimum && (invalidAll & DELTA_T_MIN)) {
        if (!message.str().empty()) message << std::endl;

        message << header.str() << "for one or more ranks the next delta T was "
                << "raised to the minimum: " << d_delTMin;
      }

      if (g_deltaT_warn_maximum && (invalidAll & DELTA_T_MAX)) {
        if (!message.str().empty()) {
          message << std::endl;
        }

        message << header.str() << "for one or more ranks the next delta T was "
                << "lowered to the maximum: " << d_delTMax;
      }

      if (g_deltaT_warn_initial && (invalidAll & DELTA_T_INITIAL_MAX)) {
        if (!message.str().empty()) {
          message << std::endl;
        }

        message << header.str() << "for one or more ranks "
                << "for the initial time up to " << d_delTInitialRange
                << " the next delT was lowered to " << d_delTInitialMax;
      }

      if (g_deltaT_warn_clamp) {
        if (invalidAll & CLAMP_TIME_TO_OUTPUT) {
          if (!message.str().empty()) {
            message << std::endl;
          }

          message << header.str()
                  << "for one or more ranks the next delta T was "
                  << "lowered to line up with the next output time: "
                  << d_output->getNextOutputTime();
        }

        if (invalidAll & CLAMP_TIME_TO_CHECKPOINT) {
          if (!message.str().empty()) {
            message << std::endl;
          }

          message << header.str()
                  << "for one or more ranks the next delta T was "
                  << "lowered to line up with the next checkpoint time: "
                  << d_output->getNextCheckpointTime();
        }

        if (invalidAll & CLAMP_TIME_TO_MAX) {
          if (!message.str().empty()) message << std::endl;

          message << header.str()
                  << "for one or more ranks the next delta T was "
                  << "lowered to line up with the maximum simulation time of "
                  << d_simTimeMax;
        }
      }

      // Finally, output the summary.
      if (!message.str().empty()) {
        DOUT(true, message.str());
      }
    }
  }

  return invalid;
}

//______________________________________________________________________
//
// Determines if the time step is the last one.
bool
SimulationCommon::isLastTimeStep(double walltime) {
  if (getReductionVariable(endSimulation_name)) {
    return true;
  }

  if (getReductionVariable(abortTimeStep_name)) {
    return true;
  }

  if (d_simTimeMax > 0 && d_simTime >= d_simTimeMax) {
    return true;
  }

  if (d_timeStepsMax > 0 && d_timeStep >= d_timeStepsMax) {
    return true;
  }

  if (d_wallTimeMax > 0) {
    // When using the wall clock time, rank 0 determines the time and
    // sends it to all other ranks.
    Uintah::MPI::Bcast(&walltime, 1, MPI_DOUBLE, 0, d_myworld->getComm());

    if (walltime >= d_wallTimeMax) {
      return true;
    }
  }

  return false;
}

//______________________________________________________________________
//
// Determines if the time step may be the last one. The simulation
// time, d_delt, and the time step are known. The only real unknown is
// the wall time for the simulation calculation. The best guess is
// based on the ExpMovingAverage of the previous time steps.
//
// MaybeLast should be called before any time step work is done.

bool
SimulationCommon::maybeLastTimeStep(double walltime) const {
  if (d_simTimeMax > 0 && d_simTime + d_delT >= d_simTimeMax) {
    return true;
  }

  if (d_timeStepsMax > 0 && d_timeStep + 1 >= d_timeStepsMax) {
    return true;
  }

  if (d_wallTimeMax > 0) {
    // When using the wall clock time, rank 0 determines the time and
    // sends it to all other ranks.
    Uintah::MPI::Bcast(&walltime, 1, MPI_DOUBLE, 0, d_myworld->getComm());

    if (walltime >= d_wallTimeMax) {
      return true;
    }
  }

  return false;
}

//______________________________________________________________________
//
// This method is called only at restart or initialization -
// see SimulationController::timeStateSetup().

void
SimulationCommon::setTimeStep(int timeStep) {
  d_timeStep = timeStep;

  // Write the time step to the initial DW so apps can get to it when
  // scheduling.
  d_scheduler->getLastDW()->override(timeStep_vartype(d_timeStep),
                                     d_timeStepLabel);
}

//______________________________________________________________________
//
void
SimulationCommon::incrementTimeStep() {
  ++d_timeStep;

  // Write the new time to the new data warehouse as the scheduler has
  // not yet advanced to the next data warehouse - see
  // SchedulerCommon::advanceDataWarehouse()
  d_scheduler->getLastDW()->override(timeStep_vartype(d_timeStep),
                                     d_timeStepLabel);
}

//______________________________________________________________________
//
// This method is called only at restart or initialization -
// see SimulationController::timeStateSetup().

void
SimulationCommon::setSimTime(double simTime) {
  d_simTime = simTime;

  // Write the time step to the initial DW so apps can get to it when
  // scheduling.
  d_scheduler->getLastDW()->override(simTime_vartype(d_simTime),
                                     d_simulationTimeLabel);
}
