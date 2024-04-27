/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#ifndef __CCA_COMPONENTS_POSTPROCESSUDA_POSTPROCESSUDA_H__
#define __CCA_COMPONENTS_POSTPROCESSUDA_POSTPROCESSUDA_H__

#include <CCA/Components/SimulationCommon/SimulationCommon.h>

#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>
#include <CCA/Components/PostProcessUda/Module.h>

#include <memory>
#include <vector>

namespace Uintah {
class LoadBalancer;
class Module;

class PostProcessUda : public SimulationCommon
{

public:
  PostProcessUda(const ProcessorGroup* myworld,
                 const MaterialManagerP materialManager,
                 const std::string& udaDir);

  virtual ~PostProcessUda();

  virtual void
  problemSetup(const ProblemSpecP& params,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               const std::string& input_ups_dir = "");

  virtual void
  scheduleInitialize(const LevelP& level, SchedulerP&);

  virtual void
  scheduleRestartInitialize([[maybe_unused]] const LevelP& level,
                            [[maybe_unused]] SchedulerP&){};

  virtual void
  restartInitialize()
  {
  }

  virtual void
  scheduleComputeStableTimestep(const LevelP&, SchedulerP&);

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&);

  virtual bool
  needRecompile(const GridP& grid);

  virtual void
  scheduleFinalizeTimestep([[maybe_unused]] const LevelP& level,
                           [[maybe_unused]] SchedulerP&){};

  // stubs
  virtual void
  scheduleInitialErrorEstimate(const LevelP&, SchedulerP&){};
  virtual void
  scheduleCoarsen(const LevelP&, SchedulerP&){};
  virtual void
  scheduleRefine(const PatchSet*, SchedulerP&){};
  virtual void
  scheduleRefineInterface(const LevelP&, SchedulerP&, bool, bool){};

  GridP
  getGrid(const GridP& currentGrid);

  virtual void
  outputProblemSpec(Uintah::ProblemSpecP&)
  {
  }

  //______________________________________________________________________
  //
private:
  PostProcessUda(const PostProcessUda&) = delete;
  PostProcessUda(PostProcessUda&&)      = delete;
  PostProcessUda&
  operator=(const PostProcessUda&) = delete;
  PostProcessUda&
  operator=(PostProcessUda&&) = delete;

  void
  computeDelT(const ProcessorGroup*,
              const PatchSubset* patches,
              const MaterialSubset* matls,
              DataWarehouse* /*old_dw*/,
              DataWarehouse* new_dw);

  void
  sched_readDataArchive(const LevelP& level, SchedulerP& sched);

  void
  readDataArchive(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* /*old_dw*/,
                  DataWarehouse* new_dw);

  void
  doAnalysis(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse* old_dw,
             DataWarehouse* new_dw);

  std::string d_udaDir;

  std::vector<int> d_udaTimesteps;
  std::vector<int> d_numMatls;
  std::vector<double> d_udaTimes;
  std::vector<double> d_udaDelT;
  std::vector<VarLabel*> d_udaSavedLabels;

  DataArchive* d_dataArchive = nullptr;
  int d_simTimestep          = 0;

  std::vector<Module*> d_Modules; // postProcess modules
  std::vector<std::unique_ptr<AnalysisModule>>
    d_analysisModules; // OnTheFly modules
};
} // End namespace Uintah

#endif //__CCA_COMPONENTS_POSTPROCESSUDA_POSTPROCESSUDA_H__
