/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Parresia Research Limited
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

#ifndef __VAANGO_CCA_COMPONENTS_MPM_UPDATE_STRESS_LAST_H__
#define __VAANGO_CCA_COMPONENTS_MPM_UPDATE_STRESS_LAST_H__

#include <CCA/Components/MPM/SerialMPM.h>

namespace Uintah {

class MPM_UpdateStressLast : public SerialMPM
{
public:
  MPM_UpdateStressLast(const ProcessorGroup* myworld,
                       const MaterialManagerP& mat_manager);
  virtual ~MPM_UpdateStressLast() = default;

  MPM_UpdateStressLast(const MPM_UpdateStressLast&) = delete;
  MPM_UpdateStressLast(MPM_UpdateStressLast&&)      = delete;
  MPM_UpdateStressLast&
  operator=(const MPM_UpdateStressLast&) = delete;
  MPM_UpdateStressLast&
  operator=(MPM_UpdateStressLast&&) = delete;

  virtual void
  scheduleTimeAdvance(const LevelP& level, SchedulerP&) override;

  virtual void
  scheduleInterpolateToParticlesAndUpdate(SchedulerP&,
                                          const PatchSet*,
                                          const MaterialSet*) override;

  virtual void
  interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw) override;
};

} // end namespace Uintah

#endif
