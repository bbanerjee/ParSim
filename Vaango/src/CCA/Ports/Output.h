/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __CCA_PORTS_OUTPUT_H__
#define __CCA_PORTS_OUTPUT_H__

#include <Core/Parallel/UintahParallelPort.h>

#include <CCA/Components/Schedulers/RuntimeStatsEnum.h>
#include <CCA/Ports/SchedulerP.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/OS/Dir.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Util/InfoMapper.h>

#include <string>

namespace Uintah {

class UintahParallelComponent;
class ProcessorGroup;

class Patch;
class VarLabel;

class Output : public UintahParallelPort
{
public:
  Output()          = default;
  virtual ~Output() = default;

  // Disallow copy and move
  Output(const Output&) = delete;
  Output(Output&&)      = delete;

  Output&
  operator=(const Output&) = delete;
  Output&
  operator=(Output&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents(UintahParallelComponent* comp) = 0;

  virtual void
  getComponents() = 0;

  virtual void
  releaseComponents() = 0;

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               const MaterialManagerP& material_manager) = 0;

  virtual void
  initializeOutput(const ProblemSpecP& params, const GridP& grid) = 0;

  // Call this when restarting from a checkpoint after calling
  // problemSetup.
  virtual void
  restartSetup(Dir& restartFromDir,
               int startTimestep,
               int timestep,
               double time,
               bool fromScratch,
               bool removeOldDir) = 0;

  // set timeinfoFlags
  virtual void
  postProcessUdaSetup(Dir& fromDir) = 0;

  virtual bool
  needRecompile(const GridP& grid) = 0;

  virtual void
  recompile(const GridP& grid) = 0;

  // Call this after all other tasks have been added to the scheduler
  virtual void
  finalizeTimestep(const GridP& grid,
                   SchedulerP& scheduler,
                   bool recompile  = false,
                   int addMaterial = 0) = 0;

  // schedule all output tasks
  virtual void
  sched_allOutputTasks(const GridP& grid,
                       SchedulerP& scheduler,
                       bool recompile = false) = 0;

  // Call this after the timestep has been executed.
  virtual void
  findNext_OutputCheckPointTimestep(bool restart, const GridP& grid) = 0;

  //! Called after a time step recompute where delta t is adjusted
  //! to make sure an output and/or checkpoint time step is needed.
  virtual void
  recompute_OutputCheckPointTimestep() = 0;

  // update or write to the xml files
  virtual void
  writeto_xml_files(const GridP& grid) = 0;

  virtual const std::string
  getOutputLocation() const = 0;
  virtual bool
  doesOutputDirExist() const = 0;

  virtual void
  setScrubSavedVariables(bool val) = 0;

  // Get the time/timestep the next output will occur
  virtual double
  getNextOutputTime() const = 0;
  virtual int
  getNextOutputTimestep() const = 0;

  // Pushes output back by one time step.
  virtual void
  postponeNextOutputTimestep() = 0;

  // Get the time/timestep the next checkpoint will occur
  virtual double
  getNextCheckpointTime() const = 0;
  virtual int
  getNextCheckpointTimestep() const = 0;
  virtual int
  getNextCheckpointWallTime() const = 0; // integer - seconds

  // Returns true if data will be output this timestep
  virtual void
  setOutputTimestep(bool val, const GridP& grid) = 0;
  virtual bool
  isOutputTimestep() const = 0;

  // Returns true if data will be checkpointed this timestep
  virtual void
  setCheckpointTimestep(bool val, const GridP& grid) = 0;
  virtual bool
  isCheckpointTimestep() const = 0;

  // Returns true if the label is being saved
  virtual bool
  isLabelSaved(const std::string& label) const = 0;

  // output interval
  virtual void
  setOutputInterval(double inv) = 0;
  virtual double
  getOutputInterval() const = 0;
  virtual void
  setOutputTimestepInterval(int inv) = 0;
  virtual int
  getOutputTimestepInterval() const = 0;

  // checkpoint interval
  virtual void
  setCheckpointInterval(double inv) = 0;
  virtual double
  getCheckpointInterval() const = 0;
  virtual void
  setCheckpointTimestepInterval(int inv) = 0;
  virtual int
  getCheckpointTimestepInterval() const = 0;
  virtual void
  setCheckpointWallTimeInterval(int inv) = 0;
  virtual int
  getCheckpointWallTimeInterval() const = 0;

  // Returns true if the UPS file has specified to save the UDA using PIDX
  // format.
  virtual bool
  savingAsPIDX() const = 0;

  // Instructs the output source (DataArchivers) on which format to use when
  // saving data.
  virtual void
  setSaveAsUDA() = 0;
  virtual void
  setSaveAsPIDX() = 0;

  virtual void
  maybeLastTimestep(bool val) = 0;
  virtual bool
  maybeLastTimestep() = 0;

  virtual void
  setElapsedWallTime(double val) = 0;
  virtual double
  getElapsedWallTime() const = 0;

  virtual void
  setCheckpointCycle(int val) = 0;
  virtual double
  getCheckpointCycle() const = 0;

  virtual void
  setUseLocalFileSystems(bool val) = 0;
  virtual bool
  getUseLocalFileSystems() const = 0;

  virtual void
  setSwitchState(bool val) = 0;
  virtual bool
  getSwitchState() const = 0;

  // Get the directory of the current time step for outputting info.
  virtual const std::string&
  getLastTimestepOutputLocation() const = 0;

  virtual void
  setRuntimeStats(
    ReductionInfoMapper<RuntimeStatsEnum, double>* runtimeStats) = 0;

  // Returns trus if an output or checkpoint exists for the time step
  virtual bool
  outputTimestepExists(unsigned int ts) = 0;
  virtual bool
  checkpointTimestepExists(unsigned int ts) = 0;
};

} // End namespace Uintah

#endif //__CCA_PORTS_OUTPUT_H__
