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

#ifndef __CCA_PORTS_MODEL_INTERFACE_H__
#define __CCA_PORTS_MODEL_INTERFACE_H__

#include <Core/Parallel/UintahParallelComponent.h>

#include <CCA/Ports/SchedulerP.h>

#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class SimulationInterface;
class Regridder;
class Output;

class DataWarehouse;
class Material;
class ProcessorGroup;
class VarLabel;

class ModelInterface : public UintahParallelComponent
{
public:
  ModelInterface(const ProcessorGroup* my_world, MaterialManagerP mat_manager);
  virtual ~ModelInterface() = default;

  // Disallow copy and move
  ModelInterface(const ModelInterface&) = delete;
  ModelInterface(ModelInterface&&)      = delete;
  ModelInterface&
  operator=(const ModelInterface&) = delete;
  ModelInterface&
  operator=(ModelInterface&&) = delete;

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents([[maybe_unused]] UintahParallelComponent* comp) override
  {
  }

  virtual void
  setComponents(SimulationInterface* comp);

  virtual void
  getComponents() override;

  virtual void
  releaseComponents() override;

  virtual void
  problemSetup(GridP& grid, const bool is_restart) = 0;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) = 0;

  virtual void
  scheduleInitialize(SchedulerP& scheduler, const LevelP& level) = 0;

  virtual void
  scheduleRestartInitialize(SchedulerP& scheduler, const LevelP& level) = 0;

  virtual void
  scheduleComputeStableTimestep(SchedulerP& sched, const LevelP& level) = 0;

  virtual void
  scheduleRefine([[maybe_unused]] const PatchSet* patches,
                 [[maybe_unused]] SchedulerP& scheduler)
  {
  }

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
  setDynamicRegridding(bool val)
  {
    d_dynamic_regridding = val;
  }
  virtual bool
  isDynamicRegridding() const
  {
    return d_dynamic_regridding;
  }

protected:
  SimulationInterface* d_simulator{ nullptr };
  Scheduler* d_scheduler{ nullptr };
  Regridder* d_regridder{ nullptr };
  Output* d_output{ nullptr };

  MaterialManagerP d_materialManager{ nullptr };

  bool d_AMR{ false };
  bool d_dynamic_regridding{ false };
};

} // End namespace Uintah

#endif //__CCA_PORTS_MODEL_INTERFACE_H__
