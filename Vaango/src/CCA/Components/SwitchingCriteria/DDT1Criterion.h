/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#ifndef Packages_Uintah_CCA_Components_Switching_DDT1_h
#define Packages_Uintah_CCA_Components_Switching_DDT1_h

#include <CCA/Ports/SwitchingCriteria.h>

#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class ProcessorGroup;
class DataWarehouse;

class ICELabel;
class MPMLabel;
class MPMICELabel;

class DDT1Criterion : public SwitchingCriteria
{
public:
  DDT1Criterion(ProblemSpecP& ps);
  virtual ~DDT1Criterion() = default;

  virtual void
  problemSetup(const ProblemSpecP& ps,
               const ProblemSpecP& restart_prob_spec,
               MaterialManagerP& materialManager);

  virtual void
  scheduleSwitchTest(const LevelP& level, SchedulerP& sched);

  void
  switchTest(const ProcessorGroup*,
             const PatchSubset* patches,
             const MaterialSubset* matls,
             DataWarehouse*,
             DataWarehouse*);

private:
  unsigned int d_material;
  double d_temperature;
  double d_BP; // Number of Particles at Boundary

  MaterialManagerP d_materialManager;
  std::unique_ptr<MPMLabel> d_mpm_labels;
  std::unique_ptr<MPMICELabel> d_mpm_ice_labels;
  std::unique_ptr<ICELabel> d_ice_labels;
#define d_SMALL_NUM 1e-100
};
} // End namespace Uintah

#endif
