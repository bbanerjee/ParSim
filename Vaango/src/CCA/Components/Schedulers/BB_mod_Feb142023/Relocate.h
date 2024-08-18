/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __CCA_COMPONENTS_SCHEDULERS_RELOCATE_H__
#define __CCA_COMPONENTS_SCHEDULERS_RELOCATE_H__

#include <Core/Grid/LevelP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ComputeSet.h>

#include <Core/Parallel/UintahMPI.h>

#include <vector>

namespace Uintah {
class DataWarehouse;
class LoadBalancer;
class ProcessorGroup;
class Scheduler;
class VarLabel;
class MPIScatterRecords;

class Relocate
{
  using VarLabel2DVector = std::vector<std::vector<const VarLabel*>>;

public:
  Relocate() = default;
  virtual ~Relocate();

  void
  scheduleParticleRelocation(Scheduler*,
                             const ProcessorGroup* pg,
                             LoadBalancer* lb,
                             const LevelP& level,
                             const VarLabel* old_posLabel,
                             const VarLabel2DVector& old_labels,
                             const VarLabel* new_posLabel,
                             const VarLabel2DVector& new_labels,
                             const VarLabel* particleIDLabel,
                             const MaterialSet* matls);

  // Schedule particle relocation without the need to provide pre-relocation
  // labels. Warning: This is experimental and has not been fully tested yet.
  // Use with caution (tsaad).
  void
  scheduleParticleRelocation(Scheduler*,
                             const ProcessorGroup* pg,
                             LoadBalancer* lb,
                             const LevelP& level,
                             const VarLabel* posLabel,
                             const VarLabel2DVector& otherLabels,
                             const MaterialSet* matls);

  const MaterialSet*
  getMaterialSet() const
  {
    return d_reloc_matls;
  }

private:

  // Callback function for particle relocation that doesn't use pre-Relocation
  // variables.
  void
  relocateParticlesModifies(const ProcessorGroup* proc,
                            const PatchSubset* patches,
                            const MaterialSubset* matls,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw,
                            const Level* coarsestLevelwithParticles);

  void
  relocateParticles(const ProcessorGroup* proc,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw,
                    const Level* coarsestLevelwithParticles);

  void
  exchangeParticles(const ProcessorGroup* proc,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw,
                    MPIScatterRecords* scatter_records,
                    int total_reloc[3]);

  void
  findNeighboringPatches(const Patch* patch,
                         const Level* level,
                         const bool findFiner,
                         const bool findCoarser,
                         Patch::selectType& AllNeighborPatches);

  void
  finalizeCommunication();

private:

  // varlabels created for the modifies version of relocation
  std::vector<const Uintah::VarLabel*> d_destroy_me;

  const VarLabel* d_reloc_old_pos_label{ nullptr };
  VarLabel2DVector d_reloc_old_labels;
  const VarLabel* d_reloc_new_pos_label{ nullptr };
  VarLabel2DVector d_reloc_new_labels;
  const VarLabel* d_particle_id_label{ nullptr };
  const MaterialSet* d_reloc_matls{ nullptr };
  LoadBalancer* d_load_balancer{ nullptr };
  std::vector<char*> d_recv_buffers;
  std::vector<char*> d_send_buffers;
  std::vector<MPI_Request> d_send_requests;

};
} // End namespace Uintah

#endif // __CCA_COMPONENTS_SCHEDULERS_RELOCATE_H__
