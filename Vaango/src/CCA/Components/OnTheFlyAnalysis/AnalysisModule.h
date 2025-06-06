/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#ifndef Packages_Uintah_CCA_Ports_AnalysisModule_h
#define Packages_Uintah_CCA_Ports_AnalysisModule_h

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/SchedulerP.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class SimulationInterface;
class Output;
class Scheduler;

class Material;
class VarLabel;

class AnalysisModule : public UintahParallelComponent
{
public:
  AnalysisModule(const ProcessorGroup* myworld,
                 const MaterialManagerP& materialManager,
                 const ProblemSpecP& module_spec);

  virtual ~AnalysisModule();

  enum objectType
  {
    plane = 0,
    line  = 1,
    box   = 2,
    none  = -9
  };

  // Methods for managing the components attached via the ports.
  virtual void
  setComponents([[maybe_unused]] UintahParallelComponent* comp){};
  virtual void
  setComponents(SimulationInterface* comp);
  virtual void
  getComponents();
  virtual void
  releaseComponents();

  virtual void
  problemSetup(const ProblemSpecP& prob_spec,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               std::vector<std::vector<const VarLabel*>>& PState,
               std::vector<std::vector<const VarLabel*>>& PState_preReloc) = 0;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) = 0;

  virtual void
  scheduleInitialize(SchedulerP& sched, const LevelP& level) = 0;

  virtual void
  scheduleRestartInitialize(SchedulerP& sched, const LevelP& level) = 0;

  virtual void
  scheduleDoAnalysis(SchedulerP& sched, const LevelP& level) = 0;

  virtual void
  scheduleDoAnalysis_preReloc(SchedulerP& sched, const LevelP& level) = 0;

  int
  createDirectory(mode_t mode,
                  const std::string& rootPath,
                  const std::string& path);

  void
  bulletProofing_LinesPlanes(const objectType obj,
                             const GridP& grid,
                             const std::string message,
                             const Point start,
                             const Point end);

  //__________________________________
  //  time related
  struct timeVars
  {
    int timeStep;
    double prevAnlysTime{ -9e9 };
    double nextAnlysTime{ -9e9 };
    double now{ -9e9 };
    bool isItTime{ false };
  };

  void
  sched_TimeVars(Task* t,
                 const LevelP& level,
                 const VarLabel* prev_AnlysTimeLabel,
                 const bool addComputes);

  bool
  getTimeVars(DataWarehouse* old_dw,
              const Level* level,
              const VarLabel* prev_AnlysTime,
              timeVars& tv);

  void
  putTimeVars(DataWarehouse* new_dw,
              const VarLabel* prev_AnlysTimeLabel,
              timeVars tv);

  bool
  isItTime(DataWarehouse* old_dw,
           const Level* level,
           const VarLabel* prev_AnlysTime);

  std::string
  getName(objectType type)
  {
    switch (type) {
      case plane:
        return "plane";
      case line:
        return "line";
      case box:
        return "box";
      default:
        throw InternalError("Invalid type", __FILE__, __LINE__);
    }
  };

  //__________________________________
  //  Variables
protected:
  SimulationInterface* m_simulator{ nullptr };
  Output* d_output{ nullptr };
  Scheduler* m_scheduler{ nullptr };

  MaterialManagerP d_materialManager{ nullptr };
  ProblemSpecP m_module_spec{ nullptr };

  const VarLabel* m_timeStepLabel{ nullptr };
  const VarLabel* m_simulationTimeLabel{ nullptr };
  const VarLabel* m_delTLabel{ nullptr };

  double m_analysisFreq; // analysis frequency: units 1/sec
  double d_startTime{ 0.0 };
  double d_stopTime{ DBL_MAX };

  Ghost::GhostType m_gn = Ghost::None;
  MaterialSubset* m_zeroMatl{ nullptr };
  MaterialSet* m_zeroMatlSet{ nullptr };

private:
  std::set<std::string>
    m_DirExists; // keep track when a directory has been created
};
}

#endif
