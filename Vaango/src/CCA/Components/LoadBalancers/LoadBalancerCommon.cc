/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>

#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/DataArchive/DataArchive.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DebugStream.h>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/NotFinished.h>

#include <sci_values.h>
#include <sstream>

namespace {

Uintah::Dout g_lb_dbg("LoadBal",
                      "LoadBalancerCommon",
                      "general info on LB patch assignment",
                      false);

Uintah::Dout g_neighborhood_dbg("LoadBal_dbg1",
                                "LoadBalancerCommon",
                                "report processor neighborhood contents",
                                false);

Uintah::Dout g_neighborhood_size_dbg(
  "LoadBal_dbg2",
  "LoadBalancerCommon",
  "report patch neighborhood sizes, local & distal",
  false);

Uintah::Dout g_patch_assignment("LoadBal_dbg3",
                                "LoadBalancerCommon",
                                "report per-process patch assignment",
                                false);

}

namespace Uintah {

DebugStream g_profile_stats("ProfileStats", "LoadBalancerCommon", "", false);
DebugStream g_profile_stats2("ProfileStats2", "LoadBalancerCommon", "", false);

// If defined, the space-filling curve will be computed in parallel,
// this may not be a good idea because the time to compute the
// space-filling curve is so small that it might not parallelize well.
#define SFC_PARALLEL

LoadBalancerCommon::LoadBalancerCommon(const ProcessorGroup* myworld)
  : UintahParallelComponent(myworld)
  , d_sfc(myworld)
  , stats("LBStats", "LoadBalancerCommon", "", false)
  , times("LBTimes", "LoadBalancerCommon", "", false)
  , lbout("LBOut", "LoadBalancerCommon", "", false)
{
}

void
LoadBalancerCommon::getComponents()
{
  d_scheduler = dynamic_cast<Scheduler*>(getPort("scheduler"));
  if (!d_scheduler) {
    throw InternalError("dynamic_cast of 'd_scheduler' failed!",
                        __FILE__,
                        __LINE__);
  }

  d_simulator = dynamic_cast<SimulationInterface*>(getPort("simulator"));
  if (!d_simulator) {
    throw InternalError("dynamic_cast of 'd_simulator' failed!",
                        __FILE__,
                        __LINE__);
  }
}

void
LoadBalancerCommon::releaseComponents()
{
  releasePort("scheduler");
  releasePort("simulator");

  d_scheduler = nullptr;
  d_simulator = nullptr;

  d_mat_manager = nullptr;
}

void
LoadBalancerCommon::assignResources(DetailedTasks& graph)
{
  int nTasks = graph.numTasks();

  DOUTR(g_lb_dbg, " Assigning Tasks to Resources! (" << nTasks << " tasks)");

  for (int i = 0; i < nTasks; i++) {
    DetailedTask* task = graph.getTask(i);

    const PatchSubset* patches = task->getPatches();
    if (patches && patches->size() > 0 &&
        task->getTask()->getType() != Task::OncePerProc &&
        task->getTask()->getType() != Task::Hypre) {
      const Patch* patch = patches->get(0);

      int idx = getPatchwiseProcessorAssignment(patch);
      ASSERTRANGE(idx, 0, d_myworld->nRanks());

      if (task->getTask()->getType() == Task::Output) {
        task->assignResource(getOutputProc(patch));
      } else {
        task->assignResource(idx);
      }

      DOUTR(g_lb_dbg,
            " Task " << *(task->getTask()) << " put on resource " << idx);

#if SCI_ASSERTION_LEVEL > 0
      std::ostringstream ostr;
      ostr << patch->getID() << ':' << idx;

      for (int i = 1; i < patches->size(); i++) {
        const Patch* p = patches->get(i);
        int pidx       = getPatchwiseProcessorAssignment(p);
        ostr << ' ' << p->getID() << ';' << pidx;
        ASSERTRANGE(pidx, 0, d_myworld->nRanks());

        if (pidx != idx && task->getTask()->getType() != Task::Output) {
          DOUTR(true,
                " WARNING: inconsistent task ("
                  << task->getTask()->getName() << ") assignment (" << pidx
                  << ", " << idx << ") in LoadBalancerCommon");
        }
      }
#endif
    } else {
      auto task_p    = task->getTask();
      auto task_type = task_p->getType();
      if (task_p->isReductionTask()) {
        task->assignResource(d_myworld->myRank());

        DOUTR(g_lb_dbg,
              "  Resource (for no patch task) "
                << *task_p << " is : " << d_myworld->myRank());

      } else if (task_type == Task::InitialSend) {
        // Already assigned, do nothing
        ASSERT(task->getAssignedResourceIndex() != -1);
      } else if (task_type == Task::OncePerProc || task_type == Task::Hypre) {

        // patch-less task, not execute-once, set to run on all procs
        // once per patch subset (empty or not)
        // at least one example is the multi-level (impAMRICE) pressureSolve
        for (auto neighbor : d_local_neighbor_processes) {
          if (patches == task_p->getPatchSet()->getSubset(neighbor)) {
            task->assignResource(neighbor);
            DOUTR(g_lb_dbg,
                  " OncePerProc Task " << *task_p << " put on resource "
                                       << neighbor);
          }
        }
      } else {
        task->assignResource(0);
        DOUTR(g_lb_dbg,
              " Unknown-type Task " << *task_p << " put on resource " << 0);
      }
    }

  } // end nTasks loop
}

int
LoadBalancerCommon::getPatchwiseProcessorAssignment(const Patch* patch)
{
  // If on a copy-data timestep and we ask about an old patch, that could cause
  // problems.
  if (d_scheduler->isCopyDataTimestep() &&
      patch->getRealPatch()->getID() < d_assignment_base_patch) {
    return -patch->getID();
  }

  ASSERTRANGE(patch->getRealPatch()->getID(),
              d_assignment_base_patch,
              d_assignment_base_patch + (int)d_processor_assignment.size());
  int proc = d_processor_assignment[patch->getRealPatch()->getGridIndex()];

  ASSERTRANGE(proc, 0, d_myworld->nRanks());
  return proc;
}

int
LoadBalancerCommon::getOldProcessorAssignment(const Patch* patch)
{
  // On an initial-regrid-timestep, this will get called from createNeighborhood
  // and can have a patch with a higher index than we have.
  if ((int)patch->getRealPatch()->getID() < d_old_assignment_base_patch ||
      patch->getRealPatch()->getID() >=
        d_old_assignment_base_patch + (int)d_old_assignment.size()) {
    return -9999;
  }

  if (patch->getGridIndex() >= (int)d_old_assignment.size()) {
    return -999;
  }

  int proc = d_old_assignment[patch->getRealPatch()->getGridIndex()];
  ASSERTRANGE(proc, 0, d_myworld->nRanks());
  return proc;
}

void
LoadBalancerCommon::useSpaceFillingCurve(const LevelP& level, int* order)
{
  std::vector<DistributedIndex> indices; // output
  std::vector<double> positions;

  IntVector min_patch_size(INT_MAX, INT_MAX, INT_MAX);

  // get the overall range in all dimensions from all patches
  IntVector high(INT_MIN, INT_MIN, INT_MIN);
  IntVector low(INT_MAX, INT_MAX, INT_MAX);

#ifdef SFC_PARALLEL
  // store how many patches each patch has originally
  std::vector<int> originalPatchCount(d_myworld->nRanks(), 0);
#endif
  for (auto iter = level->patchesBegin(); iter != level->patchesEnd(); iter++) {
    const Patch* patch = *iter;

    // calculate patchset bounds
    high = Max(high, patch->getCellHighIndex());
    low  = Min(low, patch->getCellLowIndex());

    // calculate minimum patch size
    IntVector size = patch->getCellHighIndex() - patch->getCellLowIndex();
    min_patch_size = std::min(min_patch_size, size);

    // create positions vector

#ifdef SFC_PARALLEL
    // place in long longs to avoid overflows with large numbers of patches and
    // processors
    long long pindex      = patch->getLevelIndex();
    long long num_patches = d_myworld->nRanks();
    long long proc = (pindex * num_patches) / (long long)level->numPatches();

    ASSERTRANGE(proc, 0, d_myworld->nRanks());
    if (d_myworld->myRank() == (int)proc) {
      Vector point =
        (patch->getCellLowIndex() + patch->getCellHighIndex()).asVector() / 2.0;
      for (int d = 0; d < d_numDims; d++) {
        positions.push_back(point[d_activeDims[d]]);
      }
    }
    originalPatchCount[proc]++;
#else
    Vector point =
      (patch->getCellLowIndex() + patch->getCellHighIndex()).asVector() / 2.0;
    for (int d = 0; d < dim; d++) {
      positions.push_back(point[dimensions[d]]);
    }
#endif
  }

#ifdef SFC_PARALLEL
  // compute patch starting locations
  std::vector<int> originalPatchStart(d_myworld->nRanks(), 0);
  for (int p = 1; p < d_myworld->nRanks(); p++) {
    originalPatchStart[p] =
      originalPatchStart[p - 1] + originalPatchCount[p - 1];
  }
#endif

  // Patchset dimensions
  IntVector range = high - low;

  // Center of patchset
  Vector center = (high + low).asVector() / 2.0;

  double r[3]     = { (double)range[d_activeDims[0]],
                      (double)range[d_activeDims[1]],
                      (double)range[d_activeDims[2]] };
  double c[3]     = { (double)center[d_activeDims[0]],
                      (double)center[d_activeDims[1]],
                      (double)center[d_activeDims[2]] };
  double delta[3] = { (double)min_patch_size[d_activeDims[0]],
                      (double)min_patch_size[d_activeDims[1]],
                      (double)min_patch_size[d_activeDims[2]] };

  // Create SFC
  d_sfc.SetDimensions(r);
  d_sfc.SetCenter(c);
  d_sfc.SetRefinementsByDelta(delta);
  d_sfc.SetLocations(&positions);
  d_sfc.SetOutputVector(&indices);

#ifdef SFC_PARALLEL
  d_sfc.SetLocalSize(originalPatchCount[d_myworld->myRank()]);
  d_sfc.GenerateCurve();
#else
  d_sfc.SetLocalSize(level->numPatches());
  d_sfc.GenerateCurve(SERIAL);
#endif

#ifdef SFC_PARALLEL
  if (d_myworld->nRanks() > 1) {
    std::vector<int> recvcounts(d_myworld->nRanks(), 0);
    std::vector<int> displs(d_myworld->nRanks(), 0);

    for (unsigned i = 0; i < recvcounts.size(); i++) {
      displs[i] = originalPatchStart[i] * sizeof(DistributedIndex);
      if (displs[i] < 0) {
        throw InternalError("Displacments < 0", __FILE__, __LINE__);
      }
      recvcounts[i] = originalPatchCount[i] * sizeof(DistributedIndex);
      if (recvcounts[i] < 0) {
        throw InternalError("Recvcounts < 0", __FILE__, __LINE__);
      }
    }

    std::vector<DistributedIndex> rbuf(level->numPatches());

    // Gather curve
    Uintah::MPI::Allgatherv(&indices[0],
                            recvcounts[d_myworld->myRank()],
                            MPI_BYTE,
                            &rbuf[0],
                            &recvcounts[0],
                            &displs[0],
                            MPI_BYTE,
                            d_myworld->getComm());

    indices.swap(rbuf);
  }

  // Convert distributed indices to normal indices.
  for (unsigned int i = 0; i < indices.size(); i++) {
    DistributedIndex di = indices[i];
    order[i]            = originalPatchStart[di.p] + di.i;
  }
#else
  // Write order array
  for (unsigned int i = 0; i < indices.size(); i++) {
    order[i] = indices[i].i;
  }
#endif

#if 0
  std::cout << "SFC order: ";
  for (int i = 0; i < level->numPatches(); i++) {
    std::cout << order[i] << " ";
  }
  std::cout << std::endl;
#endif
#if 0
  if(d_myworld->myRank()==0) {
    std::cout << "Warning checking SFC correctness\n";
  }
  for (int i = 0; i < level->numPatches(); i++) {
    for (int j = i+1; j < level->numPatches(); j++) {
      if (order[i] == order[j]) 
      {
        std::cout << "Rank:" << d_myworld->myRank() <<  ":   ALERT!!!!!! index done twice: index " << i << " has the same value as index " << j << " " << order[i] << std::endl;
        throw InternalError("SFC unsuccessful", __FILE__, __LINE__);
      }
    }
  }
#endif
}

void
LoadBalancerCommon::restartInitialize(DataArchive* archive,
                                      const int time_index,
                                      [[maybe_unused]] const string& ts_url,
                                      const GridP& grid)
{
  // Here we need to grab the uda data to reassign patch data to the processor
  // that will get the data.
  int num_patches          = 0;
  const Patch* first_patch = *(grid->getLevel(0)->patchesBegin());
  int startingID           = first_patch->getID();
  int prevNumProcs         = 0;

  for (int l = 0; l < grid->numLevels(); l++) {
    const LevelP& level = grid->getLevel(l);
    num_patches += level->numPatches();
  }

  d_processor_assignment.resize(num_patches);
  d_assignment_base_patch = startingID;
  for (int& i : d_processor_assignment) {
    i = -1;
  }

  if (archive->queryPatchwiseProcessor(first_patch, time_index) != -1) {
    // for uda 1.1 - if proc is saved with the patches
    for (int l = 0; l < grid->numLevels(); l++) {
      const LevelP& level = grid->getLevel(l);
      for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
           iter++) {
        d_processor_assignment[(*iter)->getID() - startingID] =
          archive->queryPatchwiseProcessor(*iter, time_index) %
          d_myworld->nRanks();
      }
    }
  } // end queryPatchwiseProcessor
  else {
    // Before uda 1.1 - DELETED THIS CODE - we don't support pre 1.1 UDAs any
    // more.
    throw InternalError(
      "LoadBalancerCommon::restartInitialize() - UDA too old...",
      __FILE__,
      __LINE__);
  }
  for (unsigned i = 0; i < d_processor_assignment.size(); i++) {
    if (d_processor_assignment[i] == -1) {
      std::cout << "index " << i << " == -1\n";
    }
    ASSERT(d_processor_assignment[i] != -1);
  }

  d_old_assignment            = d_processor_assignment;
  d_old_assignment_base_patch = d_assignment_base_patch;

  if (prevNumProcs != d_myworld->nRanks() || d_output_Nth_proc > 1) {
    if (d_myworld->myRank() == 0) {
      DOUTR(g_lb_dbg,
            "  Original run had " << prevNumProcs << ", this has "
                                  << d_myworld->nRanks());
    }
    d_check_after_restart = true;
  }

  if (d_myworld->myRank() == 0) {
    DOUTR(g_lb_dbg, " check after restart: " << d_check_after_restart);
#if 0
    int startPatch = (int)(*grid->getLevel(0)->patchesBegin())->getID();
    std::ostringstream msg;
    for (unsigned i = 0; i < d_processor_assignment.size(); i++) {
      msg << d_myworld->myRank() << " patch " << i << " (real "
          << i + startPatch << ") -> proc " << d_processor_assignment[i]
          << " (old " << d_old_assignment[i] << ") - "
          << d_processor_assignment.size() << ' ' << d_old_assignment.size()
          << "\n";
    }
    DOUTR(true, message.str();)
#endif
  }
} // end restartInitialize()

bool
LoadBalancerCommon::possiblyDynamicallyReallocate(const GridP& grid, int state)
{
  if (state != LoadBalancer::CHECK_LB) {
    // Have it create a new patch set, and have the DLB version call this.
    // This is a good place to do it, as it is automatically called when the
    // grid changes.
    d_level_perproc_patchsets.clear();
    d_output_patchsets.clear();
    d_grid_perproc_patchsets = createPerProcessorPatchSet(grid);

    for (int i = 0; i < grid->numLevels(); i++) {
      d_level_perproc_patchsets.push_back(
        createPerProcessorPatchSet(grid->getLevel(i)));
      d_output_patchsets.push_back(createOutputPatchSet(grid->getLevel(i)));
    }
  }
  return false;
}

// Creates a PatchSet containing PatchSubsets for each processor for a single
// level.
const PatchSet*
LoadBalancerCommon::createPerProcessorPatchSet(const LevelP& level)
{
  auto patches = scinew PatchSet();
  patches->createEmptySubsets(d_myworld->nRanks());

  for (auto iter = level->patchesBegin(); iter != level->patchesEnd(); iter++) {
    const Patch* patch = *iter;
    int proc           = getPatchwiseProcessorAssignment(patch);
    ASSERTRANGE(proc, 0, d_myworld->nRanks());
    PatchSubset* subset = patches->getSubset(proc);
    subset->add(patch);
  }
  patches->sortSubsets();
  return patches;
}

// Creates a PatchSet containing PatchSubsets for each processor for an
// entire grid.
const PatchSet*
LoadBalancerCommon::createPerProcessorPatchSet(const GridP& grid)
{
  auto patches = scinew PatchSet();
  patches->createEmptySubsets(d_myworld->nRanks());

  for (int i = 0; i < grid->numLevels(); i++) {
    const LevelP level = grid->getLevel(i);

    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         iter++) {
      const Patch* patch = *iter;
      int proc           = getPatchwiseProcessorAssignment(patch);
      ASSERTRANGE(proc, 0, d_myworld->nRanks());
      PatchSubset* subset = patches->getSubset(proc);
      subset->add(patch);

      // DEBUG: report patch level assignment
      if (g_patch_assignment) {
        std::ostringstream mesg;
        mesg << "  LoadBal:createPerProcessorPatchSet Patch: " << patch->getID()
             << " is on level: " << patch->getLevel()->getIndex();
        DOUTR(true, mesg.str());
      }
    }
  }
  patches->sortSubsets();

  // DEBUG: report per-proc patch assignment
  if (g_patch_assignment) {
    const PatchSubset* my_patches = patches->getSubset(d_myworld->myRank());
    std::ostringstream mesg;
    mesg << "      LoadBal:assigned patches: {";
    for (auto p = 0; p < my_patches->size(); p++) {
      mesg << ((p == 0 || p == my_patches->size()) ? " " : ", ")
           << my_patches->get(p)->getID();
    }
    mesg << " }";
    DOUTR(true, mesg.str());
  }

  return patches;
}

const PatchSet*
LoadBalancerCommon::createOutputPatchSet(const LevelP& level)
{
  if (d_output_Nth_proc == 1) {
    // assume the perProcessor set on the level was created first
    return d_level_perproc_patchsets[level->getIndex()].get_rep();
  } else {
    auto patches = scinew PatchSet();
    patches->createEmptySubsets(d_myworld->nRanks());

    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         iter++) {
      const Patch* patch = *iter;
      int proc =
        (static_cast<long long>(getPatchwiseProcessorAssignment(patch)) /
         static_cast<long long>(d_output_Nth_proc)) *
        d_output_Nth_proc;
      ASSERTRANGE(proc, 0, d_myworld->nRanks());
      PatchSubset* subset = patches->getSubset(proc);
      subset->add(patch);
    }
    patches->sortSubsets();
    return patches;
  }
}

void
LoadBalancerCommon::createNeighborhood(const GridP& grid,
                                       const GridP& oldGrid,
                                       bool hasDistalReqs)
{
  int my_rank = d_myworld->myRank();

  d_local_neighbor_patches.clear();
  d_local_neighbor_processes.clear();
  d_distal_neighbor_patches.clear();
  d_distal_neighbor_processes.clear();

  // this processor should always be in all neighborhoods
  d_local_neighbor_processes.insert(my_rank);
  d_distal_neighbor_processes.insert(my_rank);

  // go through all patches on all levels, and if the patch-wise
  // processor assignment equals the current processor, then store the
  // patch's neighbors in the load balancer array
  for (int l = 0; l < grid->numLevels(); l++) {
    LevelP level = grid->getLevel(l);

    for (auto iter = level->patchesBegin(); iter != level->patchesEnd();
         iter++) {
      const Patch* patch = *iter;

      // we need to check both where the patch is and where
      // it used to be (in the case of a dynamic reallocation)
      int proc    = getPatchwiseProcessorAssignment(patch);
      int oldproc = getOldProcessorAssignment(patch);

      IntVector low(
        patch->getExtraLowIndex(Patch::CellBased, IntVector(0, 0, 0)));
      IntVector high(
        patch->getExtraHighIndex(Patch::CellBased, IntVector(0, 0, 0)));

      // we also need to see if the output processor for patch is this proc,
      // in case it wouldn't otherwise have been in the neighborhood
      int outputproc = (static_cast<long long>(proc) /
                        static_cast<long long>(d_output_Nth_proc)) *
                       d_output_Nth_proc;

      if (proc == my_rank || oldproc == my_rank || outputproc == my_rank) {

        // add owning processors (local)
        const int maxLocalGhost = d_scheduler->getMaxGhost();
        IntVector localGhost(maxLocalGhost, maxLocalGhost, maxLocalGhost);
        addPatchesAndProcsToNeighborhood(level.get_rep(),
                                         low - localGhost,
                                         high + localGhost,
                                         d_local_neighbor_patches,
                                         d_local_neighbor_processes);

        // add owning processors (distal)
        if (hasDistalReqs) {
          const int maxDistalGhost = d_scheduler->getMaxDistalGhost();
          IntVector distalGhost(maxDistalGhost, maxDistalGhost, maxDistalGhost);
          addPatchesAndProcsToNeighborhood(level.get_rep(),
                                           low - distalGhost,
                                           high + distalGhost,
                                           d_distal_neighbor_patches,
                                           d_distal_neighbor_processes);
        }

        if (d_scheduler->isCopyDataTimestep() && proc == my_rank) {
          if (oldGrid->numLevels() > l) {
            // on copy data timestep we need old patches that line up with
            // this proc's patches, get the other way around at the end
            Patch::selectType oldPatches;
            const LevelP& oldLevel = oldGrid->getLevel(l);
            oldLevel->selectPatches(patch->getExtraCellLowIndex() - localGhost,
                                    patch->getExtraCellHighIndex() + localGhost,
                                    oldPatches);

            // add owning processors (they are the old owners)
            for (auto& old_patch_p : oldPatches) {
              d_local_neighbor_patches.insert(old_patch_p->getRealPatch());
              int nproc = getPatchwiseProcessorAssignment(old_patch_p);
              if (nproc >= 0) {
                d_local_neighbor_processes.insert(nproc);
              }
              int oproc = getOldProcessorAssignment(old_patch_p);
              if (oproc >= 0) {
                d_local_neighbor_processes.insert(oproc);
              }
            }
          }
        }

        // add amr stuff - so the patch will know about coarsening and
        // refining
        if (l > 0 &&
            (proc == my_rank ||
             (oldproc == my_rank && !d_scheduler->isCopyDataTimestep()))) {
          LevelP coarseLevel = level;

          // get the max level offset and max ghost cells to consider for
          // neighborhood creation
          int maxLevelOffset      = d_scheduler->getMaxLevelOffset();
          const int maxLocalGhost = d_scheduler->getMaxGhost();
          IntVector localGhost(maxLocalGhost, maxLocalGhost, maxLocalGhost);

          for (int offset = 1;
               offset <= maxLevelOffset && coarseLevel->hasCoarserLevel();
               ++offset) {

            // add owning processors (local)
            localGhost  = localGhost * coarseLevel->getRefinementRatio();
            coarseLevel = coarseLevel->getCoarserLevel();

            addPatchesAndProcsToNeighborhood(
              coarseLevel.get_rep(),
              level->mapCellToCoarser(low - localGhost, offset),
              level->mapCellToCoarser(high + localGhost, offset),
              d_local_neighbor_patches,
              d_local_neighbor_processes);

            // add owning processors (distal)
            if (hasDistalReqs) {
              const int maxDistalGhost = d_scheduler->getMaxDistalGhost();
              IntVector distalGhost(maxDistalGhost,
                                    maxDistalGhost,
                                    maxDistalGhost);
              distalGhost = distalGhost * coarseLevel->getRefinementRatio();

              addPatchesAndProcsToNeighborhood(
                coarseLevel.get_rep(),
                level->mapCellToCoarser(low - distalGhost, offset),
                level->mapCellToCoarser(high + distalGhost, offset),
                d_distal_neighbor_patches,
                d_distal_neighbor_processes);
            }
          }
        }

        // Second look up a single level (finer)
        if (l < grid->numLevels() - 1 &&
            (proc == my_rank ||
             (oldproc == my_rank && !d_scheduler->isCopyDataTimestep()))) {

          const LevelP& fineLevel = level->getFinerLevel();
          Patch::selectType finePatches;
          fineLevel->selectPatches(level->mapCellToFiner(low - localGhost),
                                   level->mapCellToFiner(high + localGhost),
                                   finePatches);
          for (auto& fine_patch_p : finePatches) { // add owning processors
            d_local_neighbor_patches.insert(fine_patch_p->getRealPatch());
            int nproc = getPatchwiseProcessorAssignment(fine_patch_p);
            if (nproc >= 0) {
              d_local_neighbor_processes.insert(nproc);
            }
            int oproc = getOldProcessorAssignment(fine_patch_p);
            if (oproc >= 0) {
              d_local_neighbor_processes.insert(oproc);
            }
          }
        }
      }
    }
  }

  if (d_scheduler->isCopyDataTimestep()) {
    // Regrid timestep postprocess
    // 1)- go through the old grid and
    //     find which patches used to be on this proc
    for (int l = 0; l < oldGrid->numLevels(); l++) {
      if (grid->numLevels() <= l) {
        continue;
      }

      // NOTE: all other components use uniform ghost cells across levels,
      // global and non-uniform halos are specific cases.
      LevelP oldLevel = oldGrid->getLevel(l);
      LevelP newLevel = grid->getLevel(l);

      for (auto iter = oldLevel->patchesBegin(); iter != oldLevel->patchesEnd();
           iter++) {
        const Patch* oldPatch = *iter;

        // we need to check both where the patch is and where
        // it used to be (in the case of a dynamic reallocation)
        int oldproc = getOldProcessorAssignment(oldPatch);

        if (oldproc == my_rank) {
          const int maxLocalGhost = d_scheduler->getMaxGhost();
          IntVector ghost(maxLocalGhost, maxLocalGhost, maxLocalGhost);
          Patch::selectType neighbor_patches;
          newLevel->selectPatches(oldPatch->getExtraCellLowIndex() - ghost,
                                  oldPatch->getExtraCellHighIndex() + ghost,
                                  neighbor_patches);
          d_local_neighbor_patches.insert(oldPatch);

          int nproc = getPatchwiseProcessorAssignment(oldPatch);
          if (nproc >= 0) {
            d_local_neighbor_processes.insert(nproc);
          }

          int oproc = getOldProcessorAssignment(oldPatch);
          if (oproc >= 0) {
            d_local_neighbor_processes.insert(oproc);
          }

          for (auto& neighbor_patch_p : neighbor_patches) {
            d_local_neighbor_patches.insert(neighbor_patch_p->getRealPatch());

            int nproc = getPatchwiseProcessorAssignment(neighbor_patch_p);
            if (nproc >= 0) {
              d_local_neighbor_processes.insert(nproc);
            }

            int oproc = getOldProcessorAssignment(neighbor_patch_p);
            if (oproc >= 0) {
              d_local_neighbor_processes.insert(oproc);
            }
          }
        }
      }
    }
  }

  if (g_neighborhood_dbg) {
    std::ostringstream message;

    message << "  LoadBal: Neighborhood contains procs: ";
    for (const auto& process : d_local_neighbor_processes) {
      message << process << ", ";
    }

    DOUTR(true, message.str());

    DOUTR(true, "  LoadBal: Neighborhood contains: ");
    for (const auto& patch : d_local_neighbor_patches) {
      DOUTR(true,
            "  LoadBal:    patch: " << patch->getID() << " from proc "
                                    << getPatchwiseProcessorAssignment(patch));
    }
  }

  if (g_neighborhood_size_dbg) {
    DOUTR(true,
          "  LoadBal:     m_neighbors size:      "
            << std::setw(4) << d_local_neighbor_patches.size()
            << " m_neighbor_processes size:        " << std::setw(4)
            << d_local_neighbor_processes.size());
    DOUTR(true,
          "  LoadBal:     m_distal_neighbors size: "
            << std::setw(4) << d_distal_neighbor_patches.size()
            << " m_distal_neighbor_processes size: " << std::setw(4)
            << d_distal_neighbor_processes.size());
  }
} // end createNeighborhood()

void
LoadBalancerCommon::addPatchesAndProcsToNeighborhood(
  const Level* const level,
  const IntVector& low,
  const IntVector& high,
  std::unordered_set<const Patch*>& neighbors,
  std::unordered_set<int>& processors)
{
  // each call to level->selectPatches must be done with an empty patch set
  // or otherwise it will conflict with the sorted order of the cached patches
  Patch::selectType neighborPatches;
  level->selectPatches(low, high, neighborPatches);

  for (auto i = 0u; i < neighborPatches.size(); ++i) {
    neighbors.insert(neighborPatches[i]->getRealPatch());
    int nproc = getPatchwiseProcessorAssignment(neighborPatches[i]);
    if (nproc >= 0) {
      processors.insert(nproc);
    }
    int oproc = getOldProcessorAssignment(neighborPatches[i]);
    if (oproc >= 0) {
      processors.insert(oproc);
    }
  }
}

bool
LoadBalancerCommon::inNeighborhood(const PatchSubset* patches,
                                   bool hasDistalReqs)
{
  // accept a subset with no patches as being inNeighborhood.
  if (patches->size() == 0) {
    return true;
  }

  bool found    = false;
  int patch_idx = 0;

  while (!found && patch_idx < patches->size()) {
    const Patch* patch = patches->get(patch_idx);
    if (hasDistalReqs) {
      found = (d_distal_neighbor_patches.find(patch) !=
               d_distal_neighbor_patches.end());
    } else {
      found = (d_local_neighbor_patches.find(patch) !=
               d_local_neighbor_patches.end());
    }
    patch_idx++;
  }
  return found;
}

bool
LoadBalancerCommon::inNeighborhood(const Patch* patch, bool hasDistalReqs)
{
  if (hasDistalReqs) {
    return d_distal_neighbor_patches.find(patch) !=
           d_distal_neighbor_patches.end();
  }
  return d_local_neighbor_patches.find(patch) != d_local_neighbor_patches.end();
}

void
LoadBalancerCommon::problemSetup(ProblemSpecP& pspec,
                                 [[maybe_unused]] GridP& grid,
                                 const MaterialManagerP& mat_manager)
{
  d_mat_manager = mat_manager;

  ProblemSpecP p    = pspec->findBlock("LoadBalancer");
  d_output_Nth_proc = 1;

  if (p != nullptr) {
    p->getWithDefault("outputNthProc", d_output_Nth_proc, 1);
  }
}

void
LoadBalancerCommon::setDimensionality(bool x, bool y, bool z)
{
  d_numDims = 0;

  int currentDim = 0;
  bool args[3]   = { x, y, z };

  for (int i = 0; i < 3; i++) {
    if (args[i]) {
      d_numDims++;
      d_activeDims[currentDim] = i;

      ++currentDim;
    }
  }
}

// Cost profiling functions
void
LoadBalancerCommon::addContribution([[maybe_unused]] DetailedTask* task, [[maybe_unused]] double cost)
{
  static bool warned = false;
  if (!warned) {
    proc0cout
      << "Warning: addContribution not implemented for LoadBalancerCommon.\n";
    warned = true;
  }
}

// Finalize the contributions (updates the weight, should be called once per
// timestep):
void
LoadBalancerCommon::finalizeContributions([[maybe_unused]] const GridP& currentGrid)
{
  static bool warned = false;
  if (!warned) {
    proc0cout << "Warning: finalizeContributions not implemented for "
                 "LoadBalancerCommon.\n";
    warned = true;
  }
}

// Initializes the regions in the new level that are not in the old level.
void
LoadBalancerCommon::initializeWeights([[maybe_unused]] const Grid* oldgrid, [[maybe_unused]] const Grid* newgrid)
{
  static bool warned = false;
  if (!warned) {
    proc0cout
      << "Warning: initializeWeights not implemented for LoadBalancerCommon.\n";
    warned = true;
  }
}

// Resets the profiler counters to zero
void
LoadBalancerCommon::resetCostForecaster()
{
  static bool warned = false;
  if (!warned) {
    proc0cout << "Warning: resetCostForecaster not implemented for "
                 "LoadBalancerCommon.\n";
    warned = true;
  }
}

} // namespace Uintah