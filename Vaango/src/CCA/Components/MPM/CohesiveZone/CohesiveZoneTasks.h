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

#ifndef __COHESIVE_ZONE_TASKS_H_
#define __COHESIVE_ZONE_TASKS_H_

#include <CCA/Ports/Scheduler.h>
#include <Core/GeometryPiece/GeometryPiece.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <map>
#include <vector>

namespace Uintah {
typedef int particleIndex;
typedef int particleId;

class GeometryObject;
class Patch;
class DataWarehouse;
class MPMFlags;
class CZMaterial;
class MPMLabel;
class CZLabel;
class ParticleSubset;
class VarLabel;

class CohesiveZoneTasks final
{
public:
  CohesiveZoneTasks(const ProblemSpecP& ps,
                    const MaterialManagerP& mat_manager,
                    const MPMLabel* mpm_labels,
                    const CZLabel* cz_labels,
                    const MPMFlags* mpm_flags);

  ~CohesiveZoneTasks() = default;

  void
  cohesiveZoneProblemSetup(const ProblemSpecP& prob_spec);

  void
  outputProblemSpec(ProblemSpecP& ps);

  void
  scheduleInitialize(const LevelP& level, SchedulerP& sched);

  void
  scheduleUpdate(SchedulerP& sched,
                 const PatchSet* patches,
                 const MaterialSet* matls);

  void
  scheduleAddCohesiveZoneForces(SchedulerP&,
                                const PatchSet*,
                                const MaterialSubset*,
                                const MaterialSubset*,
                                const MaterialSet*);

  void
  addCohesiveZoneForces(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw);

  void
  scheduleUpdateCohesiveZones(SchedulerP&,
                              const PatchSet*,
                              const MaterialSubset*,
                              const MaterialSubset*,
                              const MaterialSet*);

  void
  updateCohesiveZones(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw);

  using constVarLabel2DArray = std::vector<std::vector<const VarLabel*>>;
  void
  scheduleParticleRelocation(MaterialSubset* new_mss,
                             constVarLabel2DArray& old_labels,
                             constVarLabel2DArray& new_labels);

  constVarLabel2DArray d_cohesiveZoneState;
  constVarLabel2DArray d_cohesiveZoneState_preReloc;

protected:
  const MPMLabel* d_mpm_labels;
  const CZLabel* d_cz_labels;
  const MPMFlags* d_mpm_flags;
  MaterialManagerP d_mat_manager;
  int d_num_ghost_nodes;
  int d_num_ghost_particles;
};

} // End of namespace Uintah

#endif // __COHESIVE_ZONE_TASKS_H_
