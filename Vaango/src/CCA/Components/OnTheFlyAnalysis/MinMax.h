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

#ifndef Packages_Uintah_CCA_Components_ontheflyAnalysis_MinMax_h
#define Packages_Uintah_CCA_Components_ontheflyAnalysis_MinMax_h
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModule.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Output.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/GridIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <map>
#include <vector>

namespace Uintah {

/**************************************
GENERAL INFORMATION

   MinMax.h

   This computes the minimum and maximum values in the
   computational domain.

   Todd Harman
   Department of Mechanical Engineering
   University of Utah
****************************************/
class MinMax : public AnalysisModule
{
public:
  MinMax(const ProcessorGroup* myworld,
         const MaterialManagerP& materialManager,
         const ProblemSpecP& module_spec);
  MinMax() = default;

  virtual ~MinMax();

  virtual void
  problemSetup(const ProblemSpecP& prob_spec,
               const ProblemSpecP& restart_prob_spec,
               GridP& grid,
               std::vector<std::vector<const VarLabel*>>& PState,
               std::vector<std::vector<const VarLabel*>>& PState_preReloc);

  virtual void
  outputProblemSpec([[maybe_unused]] ProblemSpecP& ps){};

  virtual void
  scheduleInitialize(SchedulerP& sched, const LevelP& level);

  virtual void
  scheduleRestartInitialize(SchedulerP& sched, const LevelP& level);

  virtual void
  scheduleDoAnalysis(SchedulerP& sched, const LevelP& level);

  virtual void
  scheduleDoAnalysis_preReloc([[maybe_unused]] SchedulerP& sched, [[maybe_unused]] const LevelP& level){};

private:
  bool
  isRightLevel(const int myLevel, const int L_indx, const Level* level);

  void
  initialize(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset*,
             DataWarehouse*,
             DataWarehouse* new_dw);

  void
  computeMinMax(const ProcessorGroup* pg,
                const PatchSubset* patches,
                const MaterialSubset*,
                DataWarehouse* old_dw,
                DataWarehouse* new_dw);

  void
  doAnalysis(const ProcessorGroup* pg,
             const PatchSubset* patches,
             const MaterialSubset*,
             DataWarehouse*,
             DataWarehouse* new_dw);

  void
  createFile(const std::string& filename,
             const std::string& levelIndex,
             FILE*& fp);

  template<class Tvar, class Ttype>
  void
  findMinMax(DataWarehouse* new_dw,
             const VarLabel* varLabel,
             const int indx,
             const Patch* patch,
             GridIterator iter);

  // general labels
  class MinMaxLabel
  {
  public:
    VarLabel* lastCompTimeLabel;
    VarLabel* fileVarsStructLabel;
  };

  MinMaxLabel* d_lb;

  //__________________________________
  // global constants always begin with "d_"

  struct varProperties
  {
    VarLabel* label;
    VarLabel* reductionMinLabel;
    VarLabel* reductionMaxLabel;
    int matl;
    int level;
    MaterialSubset* matSubSet;
  };

  std::vector<varProperties> d_analyzeVars;

  const Material* d_matl{ nullptr };
  MaterialSet* d_matl_set{ nullptr };
  std::vector<bool> d_firstPassThrough;
};
} // namespace Uintah

#endif
