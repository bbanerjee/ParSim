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

#ifndef VAANGO_CCA_COMPONENTS_SIMULATIONCONTROLLER_AMRSIMULATIONCONTROLLER_H
#define VAANGO_CCA_COMPONENTS_SIMULATIONCONTROLLER_AMRSIMULATIONCONTROLLER_H

#include <CCA/Components/SimulationController/SimulationController.h>

namespace Uintah {

//! Controls the execution of an AMR Simulation
class AMRSimulationController : public SimulationController
{
public:
  AMRSimulationController(const ProcessorGroup* myworld,
                          ProblemSpecP pspec,
                          const std::string& input_ups_dir = "");

  ~AMRSimulationController() override = default;

  // eliminate copy, assignment and move
  AMRSimulationController(const AMRSimulationController&) = delete;
  AMRSimulationController(AMRSimulationController&&)      = delete;

  auto
  operator=(const AMRSimulationController&)
    -> AMRSimulationController& = delete;
  auto
  operator=(AMRSimulationController&&) -> AMRSimulationController& = delete;

  void
  run() override;

protected:
  //! Set up, compile, and execute initial timestep
  void
  doInitialTimestep();

  //! Execute a time step
  void
  executeTimestep(int totalFine);

  //! If doing AMR do the regridding
  auto
  doRegridding(bool initialTimestep) -> bool;

  void
  compileTaskGraph(int totalFine);

  auto
  doRegridding(GridP& grid, bool initialTimestep) -> bool;

  //! recursively schedule refinement, coarsening, and time advances for
  //! finer levels - compensating for time refinement.  Builds one taskgraph
  void
  subCycleCompile(int startDW, int dwStride, int numLevel, int step);

  //! recursively executes taskgraphs, as several were executed.  Similar to
  //! subCycleCompile, except that this executes the recursive taskgraphs, and
  //! compile builds one taskgraph (to exsecute once) recursively.
  void
  subCycleExecute(int startDW, int dwStride, int numLevel, bool rootCycle);

  void
  scheduleComputeStableTimestep();

protected:
  // Optional flag for scrubbing, defaulted to true.
  bool m_scrub_datawarehouse{ true };

  // Barrier timers used when running and regridding.
  Timers::Simple m_barrier_timer;
  double m_barrier_times[5];
};

} // End namespace Uintah

#endif
