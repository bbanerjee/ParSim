/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/Schedulers/OnDemandDataWarehouse.h>

#include <CCA/Components/Schedulers/DependencyException.h>
#include <CCA/Components/Schedulers/DetailedTasks.h>
#include <CCA/Components/Schedulers/MPIScheduler.h>
#include <CCA/Components/Schedulers/SchedulerCommon.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>

#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/UnknownVariable.h>

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/PSPatchMatlGhost.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <Core/Malloc/Allocator.h>

#include <Core/OS/ProcessInfo.h>

#include <Core/Parallel/BufferInfo.h>
#include <Core/Parallel/CrowdMonitor.h>
#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/FancyAssert.h>
#include <Core/Util/ProgressiveWarning.h>

#ifdef HAVE_CUDA
#include <CCA/Components/Schedulers/GPUGridVariableInfo.h>
#include <Core/Grid/Variables/GPUStencil7.h>
#include <Core/Util/DebugStream.h>
#endif

#include <climits>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// we want a particle message to have a unique tag per patch/matl/batch/dest.
// we only have 32K message tags, so this will have to do.
//   We need this because the possibility exists (particularly with DLB) of
//   two messages with the same tag being sent from the same processor.  Even
//   if these messages are sent to different processors, they can get crossed in
//   the mail or one can overwrite the other.
#define PARTICLESET_TAG 0x4000 | batch->messageTag

namespace {

// Tags for each CrowdMonitor
struct varDB_tag
{};
struct levelDB_tag
{};
struct psetDB_tag
{};
struct addsetDB_tag
{};
struct delsetDB_tag
{};
struct task_acced_send_statetag
{};

using varDB_monitor    = Uintah::CrowdMonitor<varDB_tag>;
using levelDB_monitor  = Uintah::CrowdMonitor<levelDB_tag>;
using psetDB_monitor   = Uintah::CrowdMonitor<psetDB_tag>;
using addsetDB_monitor = Uintah::CrowdMonitor<addsetDB_tag>;
using delsetDB_monitor = Uintah::CrowdMonitor<delsetDB_tag>;
using task_acced_send_statemonitor =
  Uintah::CrowdMonitor<task_acced_send_statetag>;

Dout g_foreign_dbg("ForeignVariables",
                   "OnDemandDataWarehouse",
                   "report when foreign variable is added to DW",
                   false);
Dout g_dw_get_put_dbg("OnDemandDW",
                      "OnDemandDataWarehouse",
                      "report general dbg info for OnDemandDW",
                      false);
Dout g_particles_dbg("DWParticleExchanges",
                     "OnDemandDataWarehouse",
                     "report MPI particle exchanges (sends/recvs)",
                     false);
Dout g_check_accesses("DWCheckTaskAccess",
                      "OnDemandDataWarehouse",
                      "report on task DW access checking (DBG-only)",
                      false);
Dout g_warnings_dbg("DWWarnings",
                    "OnDemandDataWarehouse",
                    "report DW GridVar progressive warnings",
                    false);

Uintah::MasterLock g_running_tasks_lock{};

}

namespace Uintah {

OnDemandDataWarehouse::OnDemandDataWarehouse(const ProcessorGroup* myworld,
                                             Scheduler* scheduler,
                                             int generation,
                                             const GridP& grid,
                                             bool isInitializationDW /*=false*/)
  : DataWarehouse(myworld, scheduler, generation)
  , d_grid{ grid }
  , d_is_initialization_DW{ isInitializationDW }
{
  doReserve();

#ifdef HAVE_CUDA
  if (Uintah::Parallel::usingDevice()) {
    int numDevices;
    cudaError_t retVal;
    CUDA_RT_SAFE_CALL(retVal = cudaGetDeviceCount(&numDevices));

    for (int i = 0; i < numDevices; i++) {
      // those gpuDWs should only live host side.
      // Ideally these don't need to be created at all as a separate
      // datawarehouse, but could be contained within this datawarehouse
      GPUDataWarehouse* gpuDW = (GPUDataWarehouse*)malloc(
        sizeof(GPUDataWarehouse) -
        sizeof(GPUDataWarehouse::dataItem) * MAX_VARDB_ITEMS);

      std::ostringstream out;
      out << "Host-side GPU DW";

      gpuDW->init(i, out.str());
      gpuDW->setDebug(gpudbg.active());
      d_gpuDWs.push_back(gpuDW);
    }
  }
#endif
}

OnDemandDataWarehouse::~OnDemandDataWarehouse()
{
  clear();
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::clear()
{
  {
    psetDB_monitor psetDB_lock{ Uintah::CrowdMonitor<psetDB_tag>::WRITER };

    for (auto iter = d_pset_DB.begin(); iter != d_pset_DB.end(); iter++) {
      if (iter->second->removeReference()) {
        delete iter->second;
      }
    }

    for (auto iter = d_delset_DB.begin(); iter != d_delset_DB.end(); iter++) {
      if (iter->second->removeReference()) {
        delete iter->second;
      }
    }

    for (auto iter = d_addset_DB.begin(); iter != d_addset_DB.end(); iter++) {
      for (auto pvar_itr = iter->second->begin();
           pvar_itr != iter->second->end();
           pvar_itr++) {
        delete pvar_itr->second;
      }
      delete iter->second;
    }
  }

  d_var_DB.clear();
  d_level_DB.clear();
  d_running_tasks.clear();

#ifdef HAVE_CUDA
  if (Uintah::Parallel::usingDevice()) {
    // clear out the host side GPU Datawarehouses.  This does NOT touch the task
    // DWs.
    for (size_t i = 0; i < d_gpuDWs.size(); i++) {
      d_gpuDWs[i]->clear();
      d_gpuDWs[i]->cleanup();
      free(d_gpuDWs[i]) d_gpuDWs[i] = nullptr;
    }
  }
#endif
}

bool
OnDemandDataWarehouse::isFinalized() const
{
  return d_finalized;
}

void
OnDemandDataWarehouse::finalize()
{
  d_var_DB.cleanForeign();
  d_finalized = true;
}

void
OnDemandDataWarehouse::unfinalize()
{
  // this is for processes that need to make small modifications to the DW
  // after it has been finalized.
  d_finalized = false;
}

void
OnDemandDataWarehouse::refinalize()
{
  d_finalized = true;
}

void
OnDemandDataWarehouse::put(Variable* var,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch)
{
  union
  {
    ReductionVariableBase* reduction;
    SoleVariableBase* sole;
    PerPatchBase* perpatch;
    ParticleVariableBase* particle;
    GridVariableBase* gv;
  } castVar;

  if ((castVar.reduction = dynamic_cast<ReductionVariableBase*>(var)) !=
      nullptr) {
    put(*castVar.reduction, label, patch ? patch->getLevel() : 0, matlIndex);
  } else if ((castVar.sole = dynamic_cast<SoleVariableBase*>(var)) != nullptr) {
    put(*castVar.sole, label, patch ? patch->getLevel() : 0, matlIndex);
  } else if ((castVar.sole = dynamic_cast<PerPatchBase*>(var)) != nullptr) {
    put(*castVar.perpatch, label, matlIndex, patch);
  } else if ((castVar.particle = dynamic_cast<ParticleVariableBase*>(var)) !=
             nullptr) {
    put(*castVar.particle, label);
  } else if ((castVar.gv = dynamic_cast<GridVariableBase*>(var)) != nullptr) {
    put(*castVar.gv, label, matlIndex, patch);
  } else {
    SCI_THROW(InternalError("Unknown Variable type", __FILE__, __LINE__));
  }
}

void
OnDemandDataWarehouse::copyKeyDB(KeyDatabase<Patch>& varkeyDB,
                                 KeyDatabase<Level>& levelkeyDB)
{
  d_var_key_DB.merge(varkeyDB);
  d_level_key_DB.merge(levelkeyDB);
}

void
OnDemandDataWarehouse::doReserve()
{
  d_var_DB.doReserve(&d_var_key_DB);
  d_level_DB.doReserve(&d_level_key_DB);
}

void
OnDemandDataWarehouse::get(ReductionVariableBase& var,
                           const VarLabel* label,
                           const Level* level,
                           int matlIndex /*= -1*/)
{
  checkGetAccess(label, matlIndex, nullptr);

  if (!d_level_DB.exists(label, matlIndex, level)) {
    std::string levelIndx = level ? to_string(level->getIndex()) : "nullptr";

    DOUTR(true,
          "get(ReductionVariableBase) failed in dw: "
            << this << ", level: " << levelIndx
            << ", matlIndex: " << matlIndex);

    d_level_DB.print(d_myworld->myRank());
    d_level_key_DB.print(d_myworld->myRank());

    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              level,
                              matlIndex,
                              "on reduction",
                              __FILE__,
                              __LINE__));
  }

  d_level_DB.get(label, matlIndex, level, var);
}

template<typename DataType, typename ReductionVarType>
std::map<int, DataType>
OnDemandDataWarehouse::getReductionVar(const VarLabel* label,
                                       const MaterialSubset* matls)
{
  std::map<int, DataType> reductionVars;
  for (auto m = 0; m < matls->size(); ++m) {
    int dwi = matls->get(m);

    ReductionVarType reductionVar;

    get(reductionVar, label, nullptr, dwi);
    reductionVars[dwi] = reductionVar;
  }

  return reductionVars;
}

template<class DataType>
std::map<int, DataType>
OnDemandDataWarehouse::get_sum_vartypeT(const VarLabel* label,
                                        const MaterialSubset* matls)
{
  using ReductionVarType =
    ReductionVariable<DataType, Reductions::Sum<DataType>>;
  return getReductionVar<DataType, ReductionVarType>(label, matls);
}

std::map<int, double>
OnDemandDataWarehouse::get_sum_vartypeD(const VarLabel* label,
                                        const MaterialSubset* matls)
{
  return get_sum_vartypeT<double>(label, matls);
}

std::map<int, Vector>
OnDemandDataWarehouse::get_sum_vartypeV(const VarLabel* label,
                                        const MaterialSubset* matls)
{
  return get_sum_vartypeT<Vector>(label, matls);
}

void
OnDemandDataWarehouse::get(SoleVariableBase& var,
                           const VarLabel* label,
                           const Level* level,
                           int matlIndex /*= -1*/)
{
  checkGetAccess(label, matlIndex, 0);

  if (!d_level_DB.exists(label, matlIndex, level)) {
    std::string levelIndx = level ? to_string(level->getIndex()) : "nullptr";

    DOUTR(true,
          "get(SoleVariableBase) failed in dw: "
            << this << ", level: " << levelIndx
            << ", matlIndex: " << matlIndex);

    d_level_DB.print(d_myworld->myRank());
    d_level_key_DB.print(d_myworld->myRank());

    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              level,
                              matlIndex,
                              "on sole",
                              __FILE__,
                              __LINE__));
  }

  d_level_DB.get(label, matlIndex, level, var);
}

bool
OnDemandDataWarehouse::exists(const VarLabel* label,
                              int matlIndex,
                              const Patch* patch) const
{
  if (patch && d_var_DB.exists(label, matlIndex, patch)) {
    return true;
  }

  // level-independent reduction vars can be stored with a null level
  if (d_level_DB.exists(label,
                        matlIndex,
                        patch ? patch->getLevel() : nullptr)) {
    return true;
  }

  return false;
}

bool
OnDemandDataWarehouse::exists(const VarLabel* label,
                              int matlIndex,
                              const Level* level) const
{
  if (level && d_level_DB.exists(label, matlIndex, level)) {
    return true;
  }

  return false;
}

bool
OnDemandDataWarehouse::exists(const VarLabel* label) const
{
  levelDB_monitor levelDB_lock{ Uintah::CrowdMonitor<levelDB_tag>::READER };

  // level-independent reduction vars can be stored with a null level
  if (d_level_DB.exists(label, -1, nullptr)) {
    return true;
  } else {
    return false;
  }
}

ReductionVariableBase*
OnDemandDataWarehouse::getReductionVariable(const VarLabel* label,
                                            int matlIndex,
                                            const Level* level) const
{
  if (d_level_DB.exists(label, matlIndex, level)) {
    ReductionVariableBase* var = dynamic_cast<ReductionVariableBase*>(
      d_level_DB.get(label, matlIndex, level));
    return var;
  } else {
    return nullptr;
  }
}

#ifdef HAVE_CUDA

void
OnDemandDataWarehouse::uintahSetCudaDevice(int deviceNum)
{
  //  CUDA_RT_SAFE_CALL( cudaSetDevice(deviceNum) );
}

int
OnDemandDataWarehouse::getNumDevices()
{
  int numDevices = 0;
  cudaError_t retVal;

  if (Uintah::Parallel::usingDevice()) {
    numDevices = 1;
  }

  // if multiple devices are desired, use this:
  CUDA_RT_SAFE_CALL(retVal = cudaGetDeviceCount(&numDevices));

  return numDevices;
}

size_t
OnDemandDataWarehouse::getTypeDescriptionSize(const TypeDescription::Type& type)
{
  switch (type) {
    case TypeDescription::Type::double_type: {
      return sizeof(double);
      break;
    }
    case TypeDescription::Type::float_type: {
      return sizeof(float);
      break;
    }
    case TypeDescription::Type::int_type: {
      return sizeof(int);
      break;
    }
    case TypeDescription::Type::Stencil7: {
      return sizeof(Stencil7);
      break;
    }
    default: {
      SCI_THROW(InternalError("OnDemandDataWarehouse::getTypeDescriptionSize "
                              "unsupported GPU Variable base type: " +
                                type,
                              __FILE__,
                              __LINE__));
    }
  }
}

GPUGridVariableBase*
OnDemandDataWarehouse::createGPUGridVariable(const TypeDescription::Type& type)
{
  // Note: For C++11, these should return a unique_ptr.
  GPUGridVariableBase* device_var = nullptr;
  switch (type) {
    case TypeDescription::Type::double_type: {
      device_var = new GPUGridVariable<double>();
      break;
    }
    case TypeDescription::Type::float_type: {
      device_var = new GPUGridVariable<float>();
      break;
    }
    case TypeDescription::Type::int_type: {
      device_var = new GPUGridVariable<int>();
      break;
    }
    case TypeDescription::Type::Stencil7: {
      device_var = new GPUGridVariable<GPUStencil7>();
      break;
    }
    default: {
      SCI_THROW(InternalError(
        "createGPUGridVariable, unsupported GPUGridVariable type: ",
        __FILE__,
        __LINE__));
    }
  }
  return device_var;
}

GPUPerPatchBase*
OnDemandDataWarehouse::createGPUPerPatch(const TypeDescription::Type& type)
{
  GPUPerPatchBase* device_var = nullptr;

  switch (type) {
    case TypeDescription::Type::double_type: {
      device_var = new GPUPerPatch<double>();
      break;
    }
    case TypeDescription::Type::float_type: {
      device_var = new GPUPerPatch<float>();
      break;
    }
    case TypeDescription::Type::int_type: {
      device_var = new GPUPerPatch<int>();
      break;
    }
    case TypeDescription::Type::Stencil7: {
      device_var = new GPUPerPatch<GPUStencil7>();
      break;
    }
    default: {
      SCI_THROW(
        InternalError("createGPUPerPatch, unsupported GPUPerPatch type: ",
                      __FILE__,
                      __LINE__));
    }
  }

  return device_var;
}

GPUReductionVariableBase*
OnDemandDataWarehouse::createGPUReductionVariable(
  const TypeDescription::Type& type)
{
  GPUReductionVariableBase* device_var = nullptr;

  switch (type) {
    case TypeDescription::Type:
    double_type : {
      device_var = new GPUReductionVariable<double>();
      break;
    }
    case TypeDescription::Type::float_type: {
      device_var = new GPUReductionVariable<float>();
      break;
    }
    case TypeDescription::Type::int_type: {
      device_var = new GPUReductionVariable<int>();
      break;
    }
    case TypeDescription::Type::Stencil7: {
      device_var = new GPUReductionVariable<GPUStencil7>();
      break;
    }
    default: {
      SCI_THROW(InternalError(
        "createGPUReductionVariable, unsupported GPUReductionVariable type: ",
        __FILE__,
        __LINE__));
    }
  }

  return device_var;
}

#endif

void
OnDemandDataWarehouse::exchangeParticleQuantities(DetailedTasks* dts,
                                                  LoadBalancer* lb,
                                                  const VarLabel* pos_var,
                                                  int iteration)
{
  // If this DW is being used for a timestep restart, then it has already done
  // this...
  if (d_exchange_particle_quantities == false) {
    return;
  }

  d_exchange_particle_quantities = false;

  ParticleExchangeVar& recvs = dts->getParticleRecvs();
  ParticleExchangeVar& sends = dts->getParticleSends();

  // need to be sized here, otherwise you risk reallocating the array after a
  // send/recv has been posted
  vector<vector<int>> senddata(sends.size()), recvdata(recvs.size());

  vector<MPI_Request> sendrequests, recvrequests;

  int data_index = 0;
  for (auto&& [proc_index, recv_set] : recvs) {
    if (recv_set.size() > 0) {
      recvdata[data_index].resize(recv_set.size());

      if (g_particles_dbg.active()) {
        std::stringstream mesg;
        mesg << d_myworld->myRank() << " Posting PARTICLES receives for "
             << recv_set.size() << " subsets from proc " << proc_index
             << " index " << data_index;
        DOUT(true, mesg.str());
      }

      MPI_Request req;
      MPI_Irecv(&(recvdata[data_index][0]),
                recv_set.size(),
                MPI_INT,
                proc_index,
                16666,
                d_myworld->getComm(),
                &req);
      recvrequests.push_back(req);
      data_index++;
    }
  }

  data_index = 0;
  for (auto&& [proc_index, send_set] : sends) {
    if (send_set.size() > 0) {
      vector<int>& data = senddata[data_index];
      data.resize(send_set.size());
      int i = 0;
      for (auto siter = send_set.begin(); siter != send_set.end();
           siter++, i++) {
        const PSPatchMatlGhostRange& pmg = *siter;
        if ((pmg.dwid_ == DetailedDep::FirstIteration && iteration > 0) ||
            (pmg.dwid_ == DetailedDep::SubsequentIterations &&
             iteration == 0)) {
          // not used
          data[i] = -2;
        } else if (pmg.dwid_ == DetailedDep::FirstIteration && iteration == 0 &&
                   lb->getOldProcessorAssignment(pmg.patch_) == proc_index) {
          // signify that the recving proc already has this data.  Only use for
          // the FirstIteration after a LB send -1 rather than force the recving
          // end above to iterate through its set
          data[i] = -1;
        } else {
          if (!d_var_DB.exists(pos_var, pmg.matl_, pmg.patch_)) {
            std::cout << d_myworld->myRank() << "  Naughty: patch "
                      << pmg.patch_->getID() << " matl " << pmg.matl_ << " id "
                      << pmg.dwid_ << endl;
            SCI_THROW(UnknownVariable(pos_var->getName(),
                                      getID(),
                                      pmg.patch_,
                                      pmg.matl_,
                                      "in exchangeParticleQuantities",
                                      __FILE__,
                                      __LINE__));
          }

          // Make sure sendset is unique...
          ASSERT(!d_send_state.find_sendset(iter->first,
                                            pmg.patch_,
                                            pmg.matl_,
                                            pmg.low_,
                                            pmg.high_,
                                            d_generation));

          ParticleSubset* p_sendset = scinew ParticleSubset(0,
                                                            pmg.matl_,
                                                            pmg.patch_,
                                                            pmg.low_,
                                                            pmg.high_);
          constParticleVariable<Point> pos;
          get(pos, pos_var, pmg.matl_, pmg.patch_);
          ParticleSubset* pset = pos.getParticleSubset();
          for (auto& particle : *pset) {
            if (Patch::containsIndex(pmg.low_,
                                     pmg.high_,
                                     pmg.patch_->getCellIndex(pos[particle]))) {
              sendset->addParticle(particle);
            }
          }
          d_send_state.add_sendset(p_sendset,
                                   proc_index,
                                   pmg.patch_,
                                   pmg.matl_,
                                   pmg.low_,
                                   pmg.high_,
                                   d_generation);
          data[i] = p_sendset->numParticles();
        }

        if (g_particles_dbg) {
          std::ostringstream mesg;
          mesg << d_myworld->myRank() << " Sending PARTICLES to proc "
               << iter->first << ": patch " << pmg.patch_->getID() << " matl "
               << pmg.matl_ << " low " << pmg.low_ << " high " << pmg.high_
               << " index " << i << ": " << senddata[data_index][i]
               << " particles";
          DOUT(true, mesg.str());
        }
      }

      DOUT(g_particles_dbg,
           d_myworld->myRank()
             << " Sending PARTICLES: " << send_set.size() << " subsets to proc "
             << proc_index << " index " << data_index);

      MPI_Request req;
      MPI_Isend(&(senddata[data_index][0]),
                s.size(),
                MPI_INT,
                proc_index,
                16666,
                d_myworld->getComm(),
                &req);
      sendrequests.push_back(req);
      data_index++;
    }
  }

  Uintah::MPI::Waitall(recvrequests.size(),
                       &recvrequests[0],
                       MPI_STATUSES_IGNORE);
  Uintah::MPI::Waitall(sendrequests.size(),
                       &sendrequests[0],
                       MPI_STATUSES_IGNORE);

  // create particle subsets from recvs
  data_index = 0;
  for (auto&& [proc_index, recv_set] : recvs) {
    if (recv_set.size() > 0) {
      vector<int>& data = recvdata[data_index];
      int i             = 0;
      for (auto riter = recv_set.begin(); riter != recv_set.end();
           riter++, i++) {
        const PSPatchMatlGhostRange& pmg = *riter;

        if (g_particles_dbg) {
          std::ostringstream mesg;
          mesg << d_myworld->myRank() << " Recving PARTICLES from proc "
               << proc_index << ": patch " << pmg.patch_->getID() << " matl "
               << pmg.matl_ << " low " << pmg.low_ << " high " << pmg.high_
               << ": " << data[i];

          DOUT(true, mesg.str());
        }

        if (data[i] == -2) {
          continue;
        }
        if (data[i] == -1) {
          ASSERT(pmg.dwid_ == DetailedDep::FirstIteration && iteration == 0 &&
                 haveParticleSubset(pmg.matl_, pmg.patch_));
          continue;
        }

        int& foreign_particles =
          d_foreign_particle_quantities[std::make_pair(pmg.matl_, pmg.patch_)];
        ParticleSubset* subset = createParticleSubset(data[i],
                                                      pmg.matl_,
                                                      pmg.patch_,
                                                      pmg.low_,
                                                      pmg.high_);

        // make room for other multiple subsets pointing into one variable -
        // additional subsets referenced at the index above the last index of
        // the previous subset
        if (data[i] > 0 && foreign_particles > 0) {

          DOUTR(g_particles_dbg,
                "  adjusting particles by " << foreign_particles);

          for (ParticleSubset::iterator iter = subset->begin();
               iter != subset->end();
               iter++) {
            *iter = *iter + foreign_particles;
          }
        }
        foreign_particles += data[i];

        DOUTR(g_particles_dbg,
              "  Setting foreign particles of patch "
                << pmg.patch_->getID() << " matl " << pmg.matl_
                << " foreign_particles" << foreign_particles);
      }
      data_index++;
    }
  }
}

void
OnDemandDataWarehouse::sendMPI(DependencyBatch* batch,
                               const VarLabel* pos_var,
                               BufferInfo& buffer,
                               OnDemandDataWarehouse* old_dw,
                               const DetailedDep* dep,
                               LoadBalancer* lb)
{
  if (dep->isNonDataDependency()) {
    // A non-data dependency -- send an empty message.
    // This would be used, for example, when a task is to modify data that
    // was previously required with ghost-cells.
    // buffer.add(0, 0, MPI_INT, false);
    return;
  }

  const VarLabel* label = dep->req->var;
  const Patch* patch    = dep->from_patch;
  int matlIndex         = dep->matl;

  switch (label->typeDescription()->getType()) {

    case TypeDescription::Type::ParticleVariable: {
      IntVector low  = dep->low;
      IntVector high = dep->high;

      if (!d_var_DB.exists(label, matlIndex, patch)) {
        SCI_THROW(UnknownVariable(label->getName(),
                                  getID(),
                                  patch,
                                  matlIndex,
                                  "in sendMPI",
                                  __FILE__,
                                  __LINE__));
      }

      ParticleVariableBase* var = dynamic_cast<ParticleVariableBase*>(
        d_var_DB.get(label, matlIndex, patch));

      int dest = batch->toTasks.front()->getAssignedResourceIndex();
      ASSERTRANGE(dest, 0, d_myworld->nRanks());

      ParticleSubset* sendset = 0;
      // first check to see if the receiving proc alrady has the (old) data
      // if data is relocating (of a regrid or re-load-balance), then the other
      // proc may already have it (since in most cases particle data comes from
      // the old dw) if lb is non-null, that means the particle data is on the
      // old dw
      if (lb && lb->getOldProcessorAssignment(patch) == dest) {
        if (this == old_dw) {
          // We don't need to know how many particles there are OR send any
          // particle data...
          return;
        }
        ASSERT(old_dw->haveParticleSubset(matlIndex, patch));
        sendset = old_dw->getParticleSubset(matlIndex, patch);

      } else {
        sendset = old_dw->d_send_state.find_sendset(dest,
                                                    patch,
                                                    matlIndex,
                                                    low,
                                                    high,
                                                    old_dw->d_generation);
      }

      // New dw send.  The NewDW doesn't yet know (on the first time) about this
      // subset if it is on a different processor.  Go ahead and calculate it,
      // but there is no need to send it, since the other proc already knows
      // about it.
      if (!sendset) {
        ASSERT(old_dw != this);

        ParticleSubset* pset = var->getParticleSubset();
        sendset = scinew ParticleSubset(0, matlIndex, patch, low, high);

        constParticleVariable<Point> pos;
        old_dw->get(pos, pos_var, pset);
        for (auto* particle : *pset) {
          if (Patch::containsIndex(low,
                                   high,
                                   patch->getCellIndex(pos[particle]))) {
            sendset->addParticle(particle);
          }
        }

        old_dw->d_send_state.add_sendset(sendset,
                                         dest,
                                         patch,
                                         matlIndex,
                                         low,
                                         high,
                                         old_dw->d_generation);
        DOUTR(g_particles_dbg,
              "  NO SENDSET: posVarLabel: "
                << pos_var->getName() << " Patch: " << patch->getID()
                << " matl " << matlIndex << " " << low << " " << high
                << " dest: " << dest << "\n    " << *sendset);

        old_dw->d_send_state.print();
      }

      ASSERT(sendset);
      if (sendset->numParticles() > 0) {
        var->getMPIBuffer(buffer, sendset);
        buffer.addSendlist(var->getRefCounted());
        buffer.addSendlist(var->getParticleSubset());
      }
    } break;

    case TypeDescription::Type::NCVariable:
    case TypeDescription::Type::CCVariable:
    case TypeDescription::Type::SFCXVariable:
    case TypeDescription::Type::SFCYVariable:
    case TypeDescription::Type::SFCZVariable: {

      if (!d_var_DB.exists(label, matlIndex, patch)) {
        DOUT(true,
             d_myworld->myRank() << "  Needed by " << *dep << " on task "
                                 << *dep->to_tasks.front());
        SCI_THROW(UnknownVariable(label->getName(),
                                  getID(),
                                  patch,
                                  matlIndex,
                                  "in sendMPI",
                                  __FILE__,
                                  __LINE__));
      }
      GridVariableBase* var;
      var =
        dynamic_cast<GridVariableBase*>(d_var_DB.get(label, matlIndex, patch));
      var->getMPIBuffer(buffer, dep->low, dep->high);
      buffer.addSendlist(var->getRefCounted());
    } break;

    case TypeDescription::Type::PetPatch:
    case TypeDescription::Type::ReductionVariable:
    case TypeDescription::Type::SoleVariable:
    default:
      SCI_THROW(InternalError("sendMPI not implemented for " +
                                label->getFullName(matlIndex, patch),
                              __FILE__,
                              __LINE__));
  } // end switch( label->getType() );
}

void
OnDemandDataWarehouse::recvMPI(DependencyBatch* batch,
                               BufferInfo& buffer,
                               OnDemandDataWarehouse* old_dw,
                               const DetailedDep* dep,
                               LoadBalancer* lb)
{
  if (dep->isNonDataDependency()) {
    // A non-data dependency -- send an empty message.
    // This would be used, for example, for dependencies between a modifying
    // task and a task the requires the data before it is to be modified.
    // buffer.add(0, 0, MPI_INT, false);
    return;
  }

  const VarLabel* label = dep->req->var;
  const Patch* patch    = dep->from_patch;
  int matlIndex         = dep->matl;
  int my_rank           = d_myworld->myRank();

  switch (label->typeDescription()->getType()) {
    case TypeDescription::Type::ParticleVariable: {
      IntVector low         = dep->low;
      IntVector high        = dep->high;
      bool whole_patch_pset = false;

      // First, get the particle set.  We should already have it
      //      if(!old_dw->haveParticleSubset(matlIndex, patch, gt, ngc)){

      // if we already have a subset for the entire patch, there's little point
      // in getting another one (and if we did, it would cause synchronization
      // problems - see comment in sendMPI)
      ParticleSubset* recvset = nullptr;
      if (lb && (lb->getOldProcessorAssignment(patch) == my_rank ||
                 lb->getPatchwiseProcessorAssignment(patch) == my_rank)) {
        // first part of the conditional means "we used to own the ghost data so
        // use the same particles" second part means "we were just assigned to
        // this patch and need to receive the whole thing" we will never get
        // here if they are both true, as mpi wouldn't need to be scheduled
        ASSERT(old_dw->haveParticleSubset(matlIndex, patch));
        recvset          = old_dw->getParticleSubset(matlIndex, patch);
        whole_patch_pset = true;
      } else {
        recvset = old_dw->getParticleSubset(matlIndex, patch, low, high);
      }
      ASSERT(recvset);

      ParticleVariableBase* var = nullptr;
      if (d_var_DB.exists(label, matlIndex, patch)) {

        var = dynamic_cast<ParticleVariableBase*>(
          d_var_DB.get(label, matlIndex, patch));
        ASSERT(var->isForeign());

      } else {

        var = dynamic_cast<ParticleVariableBase*>(
          label->typeDescription()->createInstance());
        ASSERT(var != 0);
        var->setForeign();

        // set the foreign before the allocate
        // (allocate CAN take multiple P Subsets, but only if it's foreign)
        if (whole_patch_pset) {
          var->allocate(recvset);
        } else {
          // don't give this a pset as it could be a conatiner for several
          int allocated_particles =
            old_dw
              ->d_foreign_particle_quantities[std::make_pair(matlIndex, patch)];
          var->allocate(allocated_particles);
        }
        d_var_DB
          .put(label, matlIndex, patch, var, d_scheduler->copyTimestep(), true);
      }

      if (recvset->numParticles() > 0 &&
          !(lb && lb->getOldProcessorAssignment(patch) == my_rank &&
            this == old_dw)) {
        var->getMPIBuffer(buffer, recvset);
      }

    } break;

    case TypeDescription::Type::NCVariable:

    case TypeDescription::Type::CCVariable:

    case TypeDescription::Type::SFCXVariable:

    case TypeDescription::Type::SFCYVariable:

    case TypeDescription::Type::SFCZVariable: {

      // Allocate the variable
      GridVariableBase* var = dynamic_cast<GridVariableBase*>(
        label->typeDescription()->createInstance());
      var->allocate(dep->low, dep->high);

      // Set the var as foreign
      var->setForeign();
      var->setInvalid();

      // Add the var to the dependency batch and set it as invalid.
      // The variable is now invalid because there is outstanding MPI pointing
      // to the variable.
      batch->addVar(var);
      IntVector low, high, size;
      var->getSizes(low, high, size);

      DOUT(
        g_foreign_dbg,
        "Rank-" << Parallel::getMPIRank()
                << "  adding foreign var: " << std::setw(10) << *label
                << "  patch: " << patch->getID() << "  matl: " << matlIndex
                << "  level: " << patch->getLevel()->getIndex()
                << "  from proc: " << lb->getPatchwiseProcessorAssignment(patch)
                << "  low: " << low << "  high: " << high << " sizes: " << size
                << "  num ghost cells: " << dep->m_req->m_num_ghost_cells);

      d_var_DB.putForeign(
        label,
        matlIndex,
        patch,
        var,
        d_scheduler->copyTimestep()); // put new var in data warehouse

      var->getMPIBuffer(buffer, dep->low, dep->high);
    }

    break;

    case TypeDescription::PerPatch:
    case TypeDescription::ReductionVariable:
    case TypeDescription::SoleVariable:
    default:
      SCI_THROW(InternalError("recvMPI not implemented for " +
                                label->getFullName(matlIndex, patch),
                              __FILE__,
                              __LINE__));
  } // end switch( label->getType() );

} // end recvMPI()

void
OnDemandDataWarehouse::reduceMPI(const VarLabel* label,
                                 const Level* level,
                                 const MaterialSubset* inmatls,
                                 const int nComm)
{
  const MaterialSubset* matls = inmatls;
  if (!matls) {
    std::unique_ptr<MaterialSubset> tmpmatls =
      std::make_unique<MaterialSubset>();
    tmpmatls->add(-1);
    matls = std::move(tmpmatls);
  }

  // Count the number of data elements in the reduction array
  int nmatls            = matls->size();
  int count             = 0;
  MPI_Op op             = MPI_OP_NULL;
  MPI_Datatype datatype = MPI_DATATYPE_NULL;

  for (int m = 0; m < nmatls; m++) {
    int matlIndex = matls->get(m);

    ReductionVariableBase* var;

    nt levelIndx = level ? level->getIndex() : -1;
    DOUTR(g_mpi_dbg,
          " DW:reduceMPI label: "
            << label->getName() << " matlIndex " << matlIndex
            << " level: " << levelIndx
            << " exists: " << m_level_DB.exists(label, matlIndex, level));

    if (d_level_DB.exists(label, matlIndex, level)) {

      var = dynamic_cast<ReductionVariableBase*>(
        d_level_DB.get(label, matlIndex, level));

    } else {

      //  Create and initialize the variable if it doesn't exist
      var = dynamic_cast<ReductionVariableBase*>(
        label->typeDescription()->createInstance());
      var->setBenignValue();

      d_level_DB
        .put(label, matlIndex, level, var, d_scheduler->copyTimestep(), true);
    }

    int sendcount;
    MPI_Datatype senddatatype = MPI_DATATYPE_NULL;
    MPI_Op sendop             = MPI_OP_NULL;
    var->getMPIInfo(sendcount, senddatatype, sendop);

    if (m == 0) {
      op       = sendop;
      datatype = senddatatype;
    } else {
      ASSERTEQ(op, sendop);
      ASSERTEQ(datatype, senddatatype);
    }
    count += sendcount;
  }

  int packsize;
  MPI_Pack_size(count, datatype, d_myworld->getGlobalComm(nComm), &packsize);
  vector<char> sendbuf(packsize);

  int packindex = 0;
  for (int m = 0; m < nmatls; m++) {
    int matlIndex = matls->get(m);

    ReductionVariableBase* var;
    try {
      var = dynamic_cast<ReductionVariableBase*>(
        d_level_DB.get(label, matlIndex, level));
    } catch (const UnknownVariable& e) {
      SCI_THROW(UnknownVariable(label->getName(),
                                getID(),
                                level,
                                matlIndex,
                                "on reduceMPI(pass 2)",
                                __FILE__,
                                __LINE__));
    }
    var->getMPIData(sendbuf, packindex);
  }

  vector<char> recvbuf(packsize);

  DOUTR(g_mpi_dbg,
        " allreduce, name " << label->getName() << " sendbuf.size() "
                            << sendbuf.size() << " level "
                            << (level ? level->getID() : -1));

  int error = MPI_Allreduce(&sendbuf[0],
                            &recvbuf[0],
                            count,
                            datatype,
                            op,
                            d_myworld->getGlobalComm(nComm));

  DOUTR(g_mpi_dbg,
        " allreduce, done " << label->getName() << " recvbuf.size() "
                            << recvbuf.size() << " level "
                            << (level ? level->getID() : -1));

  if (error) {
    DOUT(true, "reduceMPI: Uintah::MPI::Allreduce error: " << error);
    SCI_THROW(InternalError("reduceMPI: MPI error", __FILE__, __LINE__));
  }

  int unpackindex = 0;
  for (int m = 0; m < nmatls; m++) {
    int matlIndex = matls->get(m);

    ReductionVariableBase* var;
    try {
      var = dynamic_cast<ReductionVariableBase*>(
        d_level_DB.get(label, matlIndex, level));
    } catch (const UnknownVariable& e) {
      SCI_THROW(UnknownVariable(label->getName(),
                                getID(),
                                level,
                                matlIndex,
                                "on reduceMPI(pass 2)",
                                __FILE__,
                                __LINE__));
    }
    var->putMPIData(recvbuf, unpackindex);
  }
}

void
OnDemandDataWarehouse::put(const ReductionVariableBase& var,
                           const VarLabel* label,
                           const Level* level,
                           int matlIndex /* = -1 */)
{
  ASSERT(!d_finalized);

  checkPutAccess(label, matlIndex, nullptr,
                 false /* it actually may be replaced, but it doesn't need
                          to explicitly modify with multiple reduces in the
                          task graph */);

  //  Put it in the database
  bool init = (d_scheduler->copyTimestep()) ||
              !(d_level_DB.exists(label, matlIndex, level));
  d_level_DB.putReduce(label, matlIndex, level, var.clone(), init);
}

template<typename DataType, typename ReductionVarType>
void
OnDemandDataWarehouse::putReductionVar(std::map<int, DataType>& reductionVars,
                                       const VarLabel* label,
                                       const MaterialSubset* matls)
{
  for (auto m = 0; m < matls->size(); ++m) {
    int dwi = matls->get(m);

    using reduction_type = ReductionVariable<DataType, ReductionVarType>;
    put(reduction_type(reductionVars[dwi]), label, nullptr, dwi);
  }
}

template<class DataType>
void
OnDemandDataWarehouse::put_sum_vartypeT(std::map<int, DataType>& reductionVars,
                                        const VarLabel* label,
                                        const MaterialSubset* matls)
{
  using SumType = Reductions::Sum<DataType>;
  putReductionVar<DataType, SumType>(reductionVars, label, matls);
}

void
OnDemandDataWarehouse::put_sum_vartype(std::map<int, Vector>& reductionVars,
                                       const VarLabel* label,
                                       const MaterialSubset* matls)
{
  put_sum_vartypeT<Vector>(reductionVars, label, matls);
}

void
OnDemandDataWarehouse::put_sum_vartype(std::map<int, double>& reductionVars,
                                       const VarLabel* label,
                                       const MaterialSubset* matls)
{
  put_sum_vartypeT<double>(reductionVars, label, matls);
}

void
OnDemandDataWarehouse::override(const ReductionVariableBase& var,
                                const VarLabel* label,
                                const Level* level,
                                int matlIndex /*=-1*/)
{
  checkPutAccess(label, matlIndex, 0, true);

  printDebuggingPutInfo(label, matlIndex, level, __LINE__);

  // Put it in the database, replace whatever may already be there
  d_level_DB.put(label, matlIndex, level, var.clone(), true, true);
}

void
OnDemandDataWarehouse::override(const SoleVariableBase& var,
                                const VarLabel* label,
                                const Level* level,
                                int matlIndex /*=-1*/)
{

  checkPutAccess(label, matlIndex, 0, true);

  printDebuggingPutInfo(label, matlIndex, level, __LINE__);

  // Put it in the database, replace whatever may already be there
  d_level_DB.put(label,
                 matlIndex,
                 level,
                 var.clone(),
                 d_scheduler->copyTimestep(),
                 true);
}

void
OnDemandDataWarehouse::put(const SoleVariableBase& var,
                           const VarLabel* label,
                           const Level* level,
                           int matlIndex /* = -1 */)
{
  ASSERT(!d_finalized);

  checkPutAccess(label, matlIndex, 0,
                 false /* it actually may be replaced, but it doesn't need
                          to explicitly modify with multiple soles in the
                          task graph */);

  printDebuggingPutInfo(label, matlIndex, level, __LINE__);

  // Put it in the database
  if (!d_level_DB.exists(label, matlIndex, level)) {
    d_level_DB.put(label,
                   matlIndex,
                   level,
                   var.clone(),
                   d_scheduler->copyTimestep(),
                   false);
  }
}

ParticleSubset*
OnDemandDataWarehouse::createParticleSubset(particleIndex numParticles,
                                            int matlIndex,
                                            const Patch* patch,
                                            IntVector low /* = (0,0,0) */,
                                            IntVector high /* = (0,0,0) */)
{
  if (low == high && high == IntVector(0, 0, 0)) {
    low  = patch->getExtraCellLowIndex();
    high = patch->getExtraCellHighIndex();
  }

  if (dbg.active()) {
    dbg << d_myworld->myRank() << " DW ID " << getID()
        << " createParticleSubset: MI: " << matlIndex
        << " P: " << patch->getID() << " (" << low << ", " << high
        << ") size: " << numParticles << "\n";
  }

  ASSERT(!patch->isVirtual());

  ParticleSubset* psubset =
    scinew ParticleSubset(numParticles, matlIndex, patch, low, high);

  insertPSetRecord(d_pset_DB, patch, low, high, matlIndex, psubset);

  return psubset;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::saveParticleSubset(ParticleSubset* psubset,
                                          int matlIndex,
                                          const Patch* patch,
                                          IntVector low /* = (0,0,0) */,
                                          IntVector high /* = (0,0,0) */)
{
  ASSERTEQ(psubset->getPatch(), patch);
  ASSERTEQ(psubset->getMatlIndex(), matlIndex);
  ASSERT(!patch->isVirtual());

  if (low == high && high == IntVector(0, 0, 0)) {
    low  = patch->getExtraCellLowIndex();
    high = patch->getExtraCellHighIndex();
  }

  if (dbg.active()) {
    dbg << d_myworld->myRank() << " DW ID " << getID()
        << " saveParticleSubset: MI: " << matlIndex << " P: " << patch->getID()
        << " (" << low << ", " << high << ") size: " << psubset->numParticles()
        << "\n";
  }

  insertPSetRecord(d_pset_DB, patch, low, high, matlIndex, psubset);
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::printParticleSubsets()
{
  std::cout << "----------------------------------------------\n";
  std::cout << "-- Particle Subsets: \n\n";
  psetDBType::iterator iter;
  std::cout << d_myworld->myRank() << " Availabel psets on DW " << d_generation
            << ":\n";
  for (iter = d_pset_DB.begin(); iter != d_pset_DB.end(); iter++) {
    std::cout << d_myworld->myRank() << " " << *(iter->second) << endl;
  }
  std::cout << "----------------------------------------------\n";
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::insertPSetRecord(psetDBType& subsetDB,
                                        const Patch* patch,
                                        IntVector low,
                                        IntVector high,
                                        int matlIndex,
                                        ParticleSubset* psubset)
{
  psubset->setLow(low);
  psubset->setHigh(high);

#if SCI_ASSERTION_LEVEL >= 1
  ParticleSubset* subset =
    queryPSetDB(subsetDB, patch, matlIndex, low, high, 0, true);
  if (subset != 0) {
    if (d_myworld->myRank() == 0) {
      std::cout << d_myworld->myRank() << "  Duplicate: " << patch->getID()
                << " matl " << matlIndex << " " << low << " " << high << endl;
      printParticleSubsets();
    }
    SCI_THROW(
      InternalError("tried to create a particle subset that already exists",
                    __FILE__,
                    __LINE__));
  }
#endif

  d_plock.writeLock();
  psetDBType::key_type key(patch->getRealPatch(), matlIndex, getID());
  subsetDB.insert(
    std::pair<psetDBType::key_type, ParticleSubset*>(key, psubset));
  psubset->addReference();
  d_plock.writeUnlock();
}

//______________________________________________________________________
ParticleSubset*
OnDemandDataWarehouse::queryPSetDB(psetDBType& subsetDB,
                                   const Patch* patch,
                                   int matlIndex,
                                   IntVector low,
                                   IntVector high,
                                   const VarLabel* pos_var,
                                   bool exact)
{
  ParticleSubset* subset = 0;

  psetDBType::key_type key(patch->getRealPatch(), matlIndex, getID());
  int best_volume   = INT_MAX;
  int target_volume = Region::getVolume(low, high);

  d_plock.readLock();
  std::pair<psetDBType::const_iterator, psetDBType::const_iterator> ret =
    subsetDB.equal_range(key);

  // search multimap for best subset
  // std::cout << "Patch = " << patch << " matlIndex = " << matlIndex
  //      << " DW getID() = " << getID() << " ret = " << &(ret.first)
  //      << std::endl;
  for (auto iter = ret.first; iter != ret.second; ++iter) {

    ParticleSubset* ss = iter->second;
    IntVector sslow    = ss->getLow();
    IntVector sshigh   = ss->getHigh();
    int vol            = Region::getVolume(sslow, sshigh);

    // check if volume is better than current best
    // std::cout << " Volume = " << vol << " Best volume = " << best_volume
    //           << " Target vol = " << target_volume << std::endl;

    if (vol < best_volume) {
      // std::cout << " vol < best_volume" << std::endl;
      // intersect ranges
      // std::cout << " outside: intersect ranges: low = " << low
      //      << " sslow = " << sslow << " high = " << high
      //      << " sshigh = " << sshigh << std::endl;
      if (low.x() >= sslow.x() && low.y() >= sslow.y() &&
          low.z() >= sslow.z() && sshigh.x() >= high.x() &&
          sshigh.y() >= high.y() && sshigh.z() >= high.z()) {
        // std::cout << " inside: intersect ranges: low = " << low
        //    << " sslow = " << sslow << " high = " << high
        //    << " sshigh = " << sshigh << std::endl;

        // take this range
        subset      = ss;
        best_volume = vol;

        // short circuit out if we have already found the best possible solution
        if (best_volume == target_volume) {
          break;
        }
      }
    }
  }
  // std::cout << "2: exact = " << exact << " best_volume = " << best_volume
  //           << " target volume = " << target_volume << std::endl;
  d_plock.readUnlock();

  if (exact && best_volume != target_volume) {
    // std::cout << "queryPSetDB: exact = " << exact << " best_volume = " <<
    // best_volume
    //           << " target volume = " << target_volume << "\n";
    return 0;
  }

  // if we don't need to filter or we already have an exact match just return
  // this subset
  if (pos_var == 0 || best_volume == target_volume) {
    // std::cout << "queryPSetDB: pos_var = " << pos_var << " best_volume = " <<
    // best_volume
    //           << " target volume = " << target_volume << "\n";
    // std::cout << "queryPSetDB: particle subset size = " <<
    // subset->numParticles() << "\n";
    return subset;
  }

  // otherwise filter out particles that are not within this range
  constParticleVariable<Point> pos;

  ASSERT(subset != 0);
  get(pos, pos_var, subset);

  ParticleSubset* newsubset =
    scinew ParticleSubset(0, matlIndex, patch->getRealPatch(), low, high);

  for (ParticleSubset::iterator iter = subset->begin(); iter != subset->end();
       iter++) {
    particleIndex idx = *iter;
    if (Patch::containsIndex(low, high, patch->getCellIndex(pos[idx]))) {
      // std::cout << "queryPSetDB: particle = " << idx
      //           << " pos = " << pos[idx]
      //           << " low = " << low << " high = " << high
      //           << " cell = " << patch->getCellIndex(pos[idx]) << "\n";

      newsubset->addParticle(idx);
    }
  }

  // save subset for future queries
  d_plock.writeLock();
  subsetDB.insert(
    std::pair<psetDBType::key_type, ParticleSubset*>(key, newsubset));
  newsubset->addReference();
  d_plock.writeUnlock();

  return newsubset;
}

//______________________________________________________________________
//
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(int matlIndex, const Patch* patch)
{
  return getParticleSubset(matlIndex,
                           patch,
                           patch->getExtraCellLowIndex(),
                           patch->getExtraCellHighIndex());
}

//______________________________________________________________________
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(int matlIndex,
                                         const Patch* patch,
                                         IntVector low,
                                         IntVector high)
{
  const Patch* realPatch = (patch != 0) ? patch->getRealPatch() : 0;
  ParticleSubset* subset = 0;

  subset = queryPSetDB(d_pset_DB, realPatch, matlIndex, low, high, 0);

  // bulletproofing
  if (!subset) {
    printParticleSubsets();
    ostringstream s;
    s << "ParticleSubset, (low: " << low << ", high: " << high << " DWID "
      << getID() << ')';
    SCI_THROW(UnknownVariable(s.str().c_str(),
                              getID(),
                              realPatch,
                              matlIndex,
                              "Cannot find particle set on patch",
                              __FILE__,
                              __LINE__));
  }
  return subset;
}

//______________________________________________________________________
//
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(int matlIndex,
                                         const Patch* patch,
                                         IntVector low,
                                         IntVector high,
                                         const VarLabel* pos_var)
{
  const Patch* realPatch = (patch != 0) ? patch->getRealPatch() : 0;
  ParticleSubset* subset = 0;

  subset = queryPSetDB(d_pset_DB, realPatch, matlIndex, low, high, pos_var);

  // bulletproofing
  if (!subset) {
    printParticleSubsets();
    ostringstream s;
    s << "ParticleSubset, (low: " << low << ", high: " << high << " DWID "
      << getID() << ')';
    SCI_THROW(UnknownVariable(s.str().c_str(),
                              getID(),
                              realPatch,
                              matlIndex,
                              "Cannot find particle set on patch",
                              __FILE__,
                              __LINE__));
  }
  return subset;
}

//______________________________________________________________________
//
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(int matlIndex,
                                         const Patch* patch,
                                         Ghost::GhostType gtype,
                                         int numGhostCells,
                                         const VarLabel* pos_var)
{

  IntVector lowIndex, highIndex;
  patch->computeVariableExtents(Patch::CellBased,
                                pos_var->getBoundaryLayer(),
                                gtype,
                                numGhostCells,
                                lowIndex,
                                highIndex);

  if (gtype == Ghost::None || (lowIndex == patch->getExtraCellLowIndex() &&
                               highIndex == patch->getExtraCellHighIndex())) {
    return getParticleSubset(matlIndex, patch);
  }

  return getParticleSubset(matlIndex, lowIndex, highIndex, patch, pos_var);
}

//______________________________________________________________________
//
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(
  int matlIndex,
  IntVector lowIndex,
  IntVector highIndex,
  const Patch* relPatch,
  const VarLabel* pos_var,
  const Level*
    oldLevel) // level is ONLY used when querying from an old grid, otherwise
              // the level will be determined from the patch
{

  // relPatch can be nullptr if trying to get a particle subset for an arbitrary
  // spot on the level
  Patch::selectType neighbors;

  ASSERT(relPatch !=
         0); // you should pass in the patch on which the task was called on
  const Level* level = relPatch->getLevel();

  // compute intersection between query range and patch
  IntVector low  = Min(lowIndex, relPatch->getExtraCellLowIndex());
  IntVector high = Max(highIndex, relPatch->getExtraCellHighIndex());

  // if the user passed in the old level then query its patches
  if (oldLevel != 0) {
    oldLevel->selectPatches(
      lowIndex,
      highIndex,
      neighbors); // find all intersecting patches with the range
  }
  // if the query range is larger than the patch
  else if (low != relPatch->getExtraCellLowIndex() ||
           high != relPatch->getExtraCellHighIndex()) {
    level->selectPatches(
      lowIndex,
      highIndex,
      neighbors); // find all intersecting patches with the range
  } else {
    // just add this patch, do not query the whole level
    neighbors.push_back(relPatch);
  }

  particleIndex totalParticles = 0;
  std::vector<ParticleVariableBase*> neighborvars;
  std::vector<ParticleSubset*> subsets;
  std::vector<const Patch*> vneighbors;

  for (int i = 0; i < neighbors.size(); i++) {
    const Patch* neighbor     = neighbors[i];
    const Patch* realNeighbor = neighbor->getRealPatch();
    if (neighbor) {
      IntVector newLow;
      IntVector newHigh;

      if (level->getIndex() == 0) {
        newLow  = Max(lowIndex, neighbor->getExtraCellLowIndex());
        newHigh = Min(highIndex, neighbor->getExtraCellHighIndex());
      } else {
        // if in a copy-data timestep, only grab extra cells if on domain
        // boundary
        newLow =
          Max(lowIndex, neighbor->getLowIndexWithDomainLayer(Patch::CellBased));
        newHigh = Min(highIndex,
                      neighbor->getHighIndexWithDomainLayer(Patch::CellBased));
      }

      if (neighbor->isVirtual()) {
        // rather than offsetting each point of pos_var's data,
        // just adjust the box to compare it with.
        IntVector cellOffset = neighbor->getVirtualOffset();
        newLow -= cellOffset;
        newHigh -= cellOffset;
      }
      ParticleSubset* pset;

      if (relPatch->getLevel() == level && relPatch != neighbor) {
        relPatch->cullIntersection(Patch::CellBased,
                                   IntVector(0, 0, 0),
                                   realNeighbor,
                                   newLow,
                                   newHigh);
        if (newLow == newHigh) {
          continue;
        }
      }

      // get the particle subset for this patch
      pset = getParticleSubset(matlIndex, neighbor, newLow, newHigh, pos_var);

      // add subset to our current list
      totalParticles += pset->numParticles();
      subsets.push_back(pset);
      vneighbors.push_back(neighbors[i]);
    }
  }

  // create a new subset
  ParticleSubset* newsubset = scinew ParticleSubset(totalParticles,
                                                    matlIndex,
                                                    relPatch,
                                                    lowIndex,
                                                    highIndex,
                                                    vneighbors,
                                                    subsets);
  return newsubset;
}

/* Create a particle subset for a subset of a patch and its
   neighboring patches defined by a local lowIndex and a local highIndex.
   If the particles are contained outside the current patch, use the
   numGhostCells to get the outside particles */
ParticleSubset*
OnDemandDataWarehouse::getParticleSubset(int matlIndex,
                                         const Patch* patch,
                                         IntVector localLowIndex,
                                         IntVector localHighIndex,
                                         Ghost::GhostType ghostType,
                                         int numGhostCells,
                                         const VarLabel* posVar)
{
  // First find the total extents for the particles
  IntVector lowIndex, highIndex;
  patch->computeVariableExtents(Patch::CellBased,
                                posVar->getBoundaryLayer(),
                                ghostType,
                                numGhostCells,
                                lowIndex,
                                highIndex);

  // Check that the local lowIndex and local highIndex are contained in the
  // variable extents for the patch + extra cells
  if (localLowIndex.x() < lowIndex.x() || localLowIndex.y() < lowIndex.y() ||
      localLowIndex.z() < lowIndex.z()) {
    std::ostringstream msg_str;
    msg_str << "getParticleSubset: Cannot get variable in the local range ("
            << localLowIndex << " , " << localHighIndex << ")" << std::endl;
    msg_str << "Variable's window is (" << lowIndex << " - " << highIndex << ")"
            << std::endl;
    SCI_THROW(InternalError(msg_str.str(), __FILE__, __LINE__));
  }

  if (localHighIndex.x() > highIndex.x() ||
      localHighIndex.y() > highIndex.y() ||
      localHighIndex.z() > highIndex.z()) {
    std::ostringstream msg_str;
    msg_str << "getParticleSubset: Cannot get variable in the local range ("
            << localLowIndex << " , " << localHighIndex << ")" << std::endl;
    msg_str << "Variable's window is (" << lowIndex << " - " << highIndex << ")"
            << std::endl;
    SCI_THROW(InternalError(msg_str.str(), __FILE__, __LINE__));
  }

  return getParticleSubset(matlIndex,
                           localLowIndex,
                           localHighIndex,
                           patch,
                           posVar);
}

// Get the particle index values of a set of ParticleIDs
void
OnDemandDataWarehouse::getParticleIndex(
  ParticleSubset* pset,
  const VarLabel* partIDLabel,
  const std::vector<long64>& partIDList,
  std::vector<particleIndex>& partIndexList)
{
  // **TODO**
  // 1) Check that the label is indeed a particleIDLabel

  // Read the particleIDs for this ParticleSubset
  constParticleVariable<long64> partIDSubset;
  get(partIDSubset, partIDLabel, pset);

  // Loop thru part IDs in the input list
  for (unsigned int ii = 0; ii < partIDList.size(); ii++) {
    long64 currentID = partIDList[ii];

    // Loop through particle subset
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
         ++iter) {
      particleIndex idx = *iter;
      if (partIDSubset[idx] == currentID) {
        partIndexList.push_back(idx);
        break;
      }
    }
  }

  // Check that partIDList and partIndexList have the same size
  if (partIDList.size() != partIndexList.size()) {
    throw InternalError("Input particle ID list does not have the same size as "
                        "the output particle index list.",
                        __FILE__,
                        __LINE__);
  }
}

// Create a map between the long64 particleIDs and the particle indices in a
// ParticleSubset
void
OnDemandDataWarehouse::createParticleIDMap(ParticleSubset* pset,
                                           const VarLabel* partIDLabel,
                                           ParticleIDMap& partIDMap)
{
  // **TODO**
  // 1) Check that the label is indeed a particleIDLabel

  // Read the particleIDs for this ParticleSubset
  constParticleVariable<long64> pParticleID;
  get(pParticleID, partIDLabel, pset);

  // Loop through particle subset
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
       ++iter) {
    particleIndex idx = *iter;
    partIDMap.insert(std::pair<long64, int>(pParticleID[idx], idx));
  }
}

// Get the particle index value of a ParticleID after the partIDMap
// has been created
void
OnDemandDataWarehouse::getParticleIndex(const ParticleIDMap& partIDMap,
                                        const long64& pParticleID,
                                        particleIndex& pParticleIndex)
{
  if (partIDMap.empty()) {
    throw InternalError(
      "The map from particle IDs to the index in the particle subset is empty.",
      __FILE__,
      __LINE__);
  }

  auto mapPair = partIDMap.find(pParticleID);
  if (mapPair == partIDMap.end()) {
    throw InternalError(
      "Could not find the input particle ID in the ID->index map..",
      __FILE__,
      __LINE__);
  }
  pParticleIndex = mapPair->second;
}

//______________________________________________________________________
//
ParticleSubset*
OnDemandDataWarehouse::getDeleteSubset(int matlIndex, const Patch* patch)
{

  const Patch* realPatch = (patch != 0) ? patch->getRealPatch() : 0;
  ParticleSubset* subset = queryPSetDB(d_delset_DB,
                                       realPatch,
                                       matlIndex,
                                       patch->getExtraCellLowIndex(),
                                       patch->getExtraCellHighIndex(),
                                       0);

  if (subset == 0) {
    SCI_THROW(UnknownVariable("DeleteSet",
                              getID(),
                              realPatch,
                              matlIndex,
                              "Cannot find delete set on patch",
                              __FILE__,
                              __LINE__));
  }
  return subset;
}

ParticleLabelVariableMap*
OnDemandDataWarehouse::getNewParticleState(int matlIndex, const Patch* patch)
{
  d_pslock.readLock();
  const Patch* realPatch = (patch != 0) ? patch->getRealPatch() : 0;
  psetAddDBType::key_type key(matlIndex, realPatch);
  psetAddDBType::iterator iter = d_addset_DB.find(key);
  if (iter == d_addset_DB.end()) {
    d_pslock.readUnlock();
    return 0;
  }
  d_pslock.readUnlock();
  return iter->second;
}

//______________________________________________________________________
//
bool
OnDemandDataWarehouse::haveParticleSubset(int matlIndex,
                                          const Patch* patch,
                                          IntVector low /* = (0,0,0) */,
                                          IntVector high /* = (0,0,0) */,
                                          bool exact /*=false*/)
{
  if (low == high && high == IntVector(0, 0, 0)) {
    low  = patch->getExtraCellLowIndex();
    high = patch->getExtraCellHighIndex();
  }
  const Patch* realPatch = patch->getRealPatch();
  // query subset
  ParticleSubset* subset =
    queryPSetDB(d_pset_DB, realPatch, matlIndex, low, high, 0);

  // if no subset was returned there are no suitable subsets
  if (subset == 0) {
    return false;
  }

  if (exact) { // check if the user wanted an exact match
    return subset->getLow() == low && subset->getHigh() == high;
  } else {
    return true;
  }
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::get(constParticleVariableBase& constVar,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch)
{

  checkGetAccess(label, matlIndex, patch);

  if (!d_var_DB.exists(label, matlIndex, patch)) {
    print();
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              patch,
                              matlIndex,
                              "",
                              __FILE__,
                              __LINE__));
  }
  constVar =
    *dynamic_cast<ParticleVariableBase*>(d_var_DB.get(label, matlIndex, patch));
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::get(constParticleVariableBase& constVar,
                           const VarLabel* label,
                           ParticleSubset* pset)
{
  int matlIndex      = pset->getMatlIndex();
  const Patch* patch = pset->getPatch();

  // pset center patch and neighbor patch are not in same level
  // (probably on an AMR copy data timestep)
  if ((pset->getNeighbors().size() == 0) ||
      (pset->getNeighbors().front()->getLevel() == patch->getLevel() &&
       pset->getLow() == patch->getExtraCellLowIndex() &&
       pset->getHigh() == patch->getExtraCellHighIndex())) {
    get(constVar, label, matlIndex, patch);
  } else {
    checkGetAccess(label, matlIndex, patch);
    ParticleVariableBase* var = constVar.cloneType();

    const vector<const Patch*>& neighborPatches = pset->getNeighbors();
    const vector<ParticleSubset*>& neighbor_subsets =
      pset->getNeighborSubsets();

    vector<ParticleVariableBase*> neighborvars(neighborPatches.size());

    for (size_t i = 0; i < neighborPatches.size(); i++) {
      const Patch* neighborPatch = neighborPatches[i];

      if (!d_var_DB.exists(label, matlIndex, neighborPatches[i])) {
        SCI_THROW(
          UnknownVariable(label->getName(),
                          getID(),
                          neighborPatch,
                          matlIndex,
                          neighborPatch == patch ? "on patch" : "on neighbor",
                          __FILE__,
                          __LINE__));
      }

      neighborvars[i] = var->cloneType();

      d_var_DB.get(label, matlIndex, neighborPatch, *neighborvars[i]);
    }

    // Note that when the neighbors are virtual patches (i.e. periodic
    // boundaries), then if var is a ParticleVariable<Point>, the points
    // of neighbors will be translated by its virtualOffset.

    var->gather(pset, neighbor_subsets, neighborvars, neighborPatches);

    constVar = *var;

    for (size_t i = 0; i < neighborPatches.size(); i++) {
      delete neighborvars[i];
    }
    delete var;
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::getModifiable(ParticleVariableBase& var,
                                     const VarLabel* label,
                                     ParticleSubset* pset)
{
  int matlIndex      = pset->getMatlIndex();
  const Patch* patch = pset->getPatch();
  checkModifyAccess(label, matlIndex, patch);

  if (pset->getLow() == patch->getExtraCellLowIndex() &&
      pset->getHigh() == patch->getExtraCellHighIndex()) {
    if (!d_var_DB.exists(label, matlIndex, patch)) {
      SCI_THROW(UnknownVariable(label->getName(),
                                getID(),
                                patch,
                                matlIndex,
                                "",
                                __FILE__,
                                __LINE__));
    }
    d_var_DB.get(label, matlIndex, patch, var);
  } else {
    SCI_THROW(InternalError("getModifiable (Particle Variable (" +
                              label->getName() +
                              ") ).  The particleSubset low/high index does "
                              "not match the patch low/high indices",
                            __FILE__,
                            __LINE__));
  }
}
//______________________________________________________________________
//
ParticleVariableBase*
OnDemandDataWarehouse::getParticleVariable(const VarLabel* label,
                                           ParticleSubset* pset)
{
  int matlIndex      = pset->getMatlIndex();
  const Patch* patch = pset->getPatch();

  if (pset->getLow() == patch->getExtraCellLowIndex() &&
      pset->getHigh() == patch->getExtraCellHighIndex()) {
    return getParticleVariable(label, matlIndex, patch);
  } else {
    SCI_THROW(InternalError("getParticleVariable (Particle Variable (" +
                              label->getName() +
                              ") ).  The particleSubset low/high index does "
                              "not match the patch low/high indices",
                            __FILE__,
                            __LINE__));
  }
}
//______________________________________________________________________
//
ParticleVariableBase*
OnDemandDataWarehouse::getParticleVariable(const VarLabel* label,
                                           int matlIndex,
                                           const Patch* patch)
{
  ParticleVariableBase* var = 0;

  // in case the it's a virtual patch -- only deal with real patches
  if (patch != 0) {
    patch = patch->getRealPatch();
  }

  checkModifyAccess(label, matlIndex, patch);

  if (!d_var_DB.exists(label, matlIndex, patch)) {
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              patch,
                              matlIndex,
                              "",
                              __FILE__,
                              __LINE__));
  }
  var =
    dynamic_cast<ParticleVariableBase*>(d_var_DB.get(label, matlIndex, patch));

  return var;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::allocateTemporary(ParticleVariableBase& var,
                                         ParticleSubset* pset)
{
  var.allocate(pset);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::allocateAndPut(ParticleVariableBase& var,
                                      const VarLabel* label,
                                      ParticleSubset* pset)
{

  int matlIndex      = pset->getMatlIndex();
  const Patch* patch = pset->getPatch();

  // Error checking
  if (d_var_DB.exists(label, matlIndex, patch)) {
    SCI_THROW(
      InternalError("Particle variable already exists: " + label->getName(),
                    __FILE__,
                    __LINE__));
  }

  var.allocate(pset);
  put(var, label);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::put(ParticleVariableBase& var,
                           const VarLabel* label,
                           bool replace /*= false*/)
{
  ASSERT(!d_finalized);

  ParticleSubset* pset = var.getParticleSubset();

  const Patch* patch = pset->getPatch();

  if (pset->getLow() != patch->getExtraCellLowIndex() ||
      pset->getHigh() != patch->getExtraCellHighIndex()) {
    SCI_THROW(InternalError(" put(Particle Variable (" + label->getName() +
                              ") ).  The particleSubset low/high index does "
                              "not match the patch low/high indices",
                            __FILE__,
                            __LINE__));
  }

  int matlIndex = pset->getMatlIndex();

  checkPutAccess(label, matlIndex, patch, replace);

  // Put it in the database
  printDebuggingPutInfo(label, matlIndex, patch, __LINE__);
  d_var_DB.put(label,
               matlIndex,
               patch,
               var.clone(),
               d_scheduler->copyTimestep(),
               replace);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::copyOut(ParticleVariableBase& var,
                               const VarLabel* label,
                               ParticleSubset* pset)
{
  constParticleVariableBase* constVar = var.cloneConstType();
  this->get(*constVar, label, pset);
  var.copyData(&constVar->getBaseRep());
  delete constVar;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::getCopy(ParticleVariableBase& var,
                               const VarLabel* label,
                               ParticleSubset* pset)
{
  constParticleVariableBase* constVar = var.cloneConstType();
  this->get(*constVar, label, pset);
  var.allocate(pset);
  var.copyData(&constVar->getBaseRep());
  delete constVar;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::get(constGridVariableBase& constVar,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch,
                           Ghost::GhostType gtype,
                           int numGhostCells)
{
  GridVariableBase* var = constVar.cloneType();

  checkGetAccess(label, matlIndex, patch, gtype, numGhostCells);
  getGridVar(*var, label, matlIndex, patch, gtype, numGhostCells);

  constVar = *var;
  delete var;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::getModifiable(GridVariableBase& var,
                                     const VarLabel* label,
                                     int matlIndex,
                                     const Patch* patch,
                                     Ghost::GhostType gtype,
                                     int numGhostCells)
{
  // checkModifyAccess(label, matlIndex, patch);
  getGridVar(var, label, matlIndex, patch, gtype, numGhostCells);
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::allocateTemporary(GridVariableBase& var,
                                         const Patch* patch,
                                         Ghost::GhostType gtype,
                                         int numGhostCells)
{
  IntVector boundaryLayer(0, 0, 0); // Is this right?

  IntVector lowIndex, highIndex;
  IntVector lowOffset, highOffset;
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(var.virtualGetTypeDescription()->getType(),
                                false);
  Patch::getGhostOffsets(var.virtualGetTypeDescription()->getType(),
                         gtype,
                         numGhostCells,
                         lowOffset,
                         highOffset);

  patch->computeExtents(basis,
                        boundaryLayer,
                        lowOffset,
                        highOffset,
                        lowIndex,
                        highIndex);

  var.allocate(lowIndex, highIndex);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::allocateAndPut(GridVariableBase& var,
                                      const VarLabel* label,
                                      int matlIndex,
                                      const Patch* patch,
                                      Ghost::GhostType gtype,
                                      int numGhostCells)
{
  if (d_finalized) {
    std::cerr << "  DW " << getID() << " finalized!\n";
  }
  ASSERT(!d_finalized);

  // Note: almost the entire function is write locked in order to prevent dual
  // allocations in a multi-threaded environment.  Whichever patch in a
  // super patch group gets here first, does the allocating for the entire
  // super patch group.

#if 0
  if (!hasRunningTask()) {
    SCI_THROW(InternalError("OnDemandDataWarehouse::AllocateAndPutGridVar can only be used when the dw has a running task associated with it.", __FILE__, __LINE__));
  }
#endif

  checkPutAccess(label, matlIndex, patch, false);
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);

  IntVector lowIndex, highIndex;
  IntVector lowOffset, highOffset;
  Patch::getGhostOffsets(var.virtualGetTypeDescription()->getType(),
                         gtype,
                         numGhostCells,
                         lowOffset,
                         highOffset);
  patch->computeExtents(basis,
                        label->getBoundaryLayer(),
                        lowOffset,
                        highOffset,
                        lowIndex,
                        highIndex);

  if (!d_combineMemory) {
    bool exists = d_var_DB.exists(label, matlIndex, patch);
    if (exists) {
      // it had been allocated and put as part of the superpatch of
      // another patch
      d_var_DB.get(label, matlIndex, patch, var);

      // The var's window should be the size of the patch or smaller than it.
      ASSERTEQ(Min(var.getLow(), lowIndex), lowIndex);
      ASSERTEQ(Max(var.getHigh(), highIndex), highIndex);

      // this is just a tricky way to uninitialize var
      Variable* tmpVar = dynamic_cast<Variable*>(var.cloneType());
      var.copyPointer(*tmpVar);
      delete tmpVar;
    }
    // allocate the memory
    var.allocate(lowIndex, highIndex);

    printDebuggingPutInfo(label, matlIndex, patch, __LINE__);

    // put the variable in the database
    d_var_DB.put(label,
                 matlIndex,
                 patch,
                 var.clone(),
                 d_scheduler->copyTimestep(),
                 true);

  } else {
    d_lock.writeLock();
    bool exists = d_var_DB.exists(label, matlIndex, patch);
    if (exists) {

      // it had been allocated and put as part of the superpatch of
      // another patch
      d_var_DB.get(label, matlIndex, patch, var);

      // The var's window should be the size of the patch or smaller than it.
      ASSERTEQ(Min(var.getLow(), lowIndex), lowIndex);
      ASSERTEQ(Max(var.getHigh(), highIndex), highIndex);

      if (var.getLow() !=
            patch->getExtraLowIndex(basis, label->getBoundaryLayer()) ||
          var.getHigh() !=
            patch->getExtraHighIndex(basis, label->getBoundaryLayer()) ||
          var.getBasePointer() == 0 /* place holder for ghost patch */) {
        // It wasn't allocated as part of another patch's superpatch;
        // it existed as ghost patch of another patch.. so we have no
        // choice but to blow it away and replace it.
        d_var_DB
          .put(label, matlIndex, patch, 0, d_scheduler->copyTimestep(), true);

        // this is just a tricky way to uninitialize var
        Variable* tmpVar = dynamic_cast<Variable*>(var.cloneType());
        var.copyPointer(*tmpVar);
        delete tmpVar;
      } else {
        // It was allocated and put as part of the superpatch of another patch
        d_lock.writeUnlock();
        var.rewindow(lowIndex, highIndex);
        return; // got it -- done
      }
    }

    IntVector superLowIndex, superHighIndex;
    // requiredSuper[Low/High]'s don't take numGhostCells into consideration
    // -- just includes ghosts that will be required by later tasks.
    IntVector requiredSuperLow, requiredSuperHigh;

    const vector<const Patch*>* superPatchGroup =
      d_scheduler->getSuperPatchExtents(label,
                                        matlIndex,
                                        patch,
                                        gtype,
                                        numGhostCells,
                                        requiredSuperLow,
                                        requiredSuperHigh,
                                        superLowIndex,
                                        superHighIndex);

    ASSERT(superPatchGroup != 0);

    var.allocate(superLowIndex, superHighIndex);

#if SCI_ASSERTION_LEVEL >= 3

    // check for dead portions of a variable (variable space that isn't covered
    // by any patch). This will happen with L-shaped patch configs and ngc >
    // extra cells. find all dead space and mark it with a bogus value.

    if (1) { // numGhostCells > ec) { (numGhostCells is 0, query it from the
             // superLowIndex...
      std::deque<Box> b1, b2, difference;
      b1.push_back(
        Box(Point(superLowIndex(0), superLowIndex(1), superLowIndex(2)),
            Point(superHighIndex(0), superHighIndex(1), superHighIndex(2))));
      for (size_t i = 0; i < (*superPatchGroup).size(); i++) {
        const Patch* p = (*superPatchGroup)[i];
        IntVector low  = p->getExtraLowIndex(basis, label->getBoundaryLayer());
        IntVector high = p->getExtraHighIndex(basis, label->getBoundaryLayer());
        b2.push_back(
          Box(Point(low(0), low(1), low(2)), Point(high(0), high(1), high(2))));
      }
      difference = Box::difference(b1, b2);

#if 0
      if (difference.size() > 0) {
        std::cout << "Box difference: " << superLowIndex << " " << superHighIndex << " with patches " << endl;
        for (size_t i = 0; i < (*superPatchGroup).size(); i++) {
          const Patch* p = (*superPatchGroup)[i];
          std::cout << p->getExtraLowIndex(basis, label->getBoundaryLayer()) << " " << p->getExtraHighIndex(basis, label->getBoundaryLayer()) << endl;
        }

        for (size_t i = 0; i < difference.size(); i++) {
          std::cout << difference[i].lower() << " " << difference[i].upper() << endl;
        }
      }
#endif
      // get more efficient way of doing this...
      for (size_t i = 0; i < difference.size(); i++) {
        Box b = difference[i];
        IntVector low((int)b.lower()(0), (int)b.lower()(1), (int)b.lower()(2));
        IntVector high((int)b.upper()(0), (int)b.upper()(1), (int)b.upper()(2));
        if (GridVariable<double>* typedVar =
              dynamic_cast<GridVariable<double>*>(&var)) {
          for (CellIterator iter(low, high); !iter.done(); iter++) {
            (*typedVar)[*iter] = -5.555555e256;
          }
        } else if (GridVariable<Vector>* typedVar =
                     dynamic_cast<GridVariable<Vector>*>(&var)) {
          for (CellIterator iter(low, high); !iter.done(); iter++) {
            (*typedVar)[*iter] = -5.555555e256;
          }
        }
      }
    }
#endif

    Patch::selectType encompassedPatches;
    if (requiredSuperLow == lowIndex && requiredSuperHigh == highIndex) {
      // only encompassing the patch currently being allocated
      encompassedPatches.push_back(patch);
    } else {
      // Use requiredSuperLow/High instead of superLowIndex/superHighIndex
      // so we don't put the var for patches in the datawarehouse that won't be
      // required (this is important for scrubbing).
      patch->getLevel()->selectPatches(requiredSuperLow,
                                       requiredSuperHigh,
                                       encompassedPatches);
    }

    // Make a set of the non ghost patches that
    // has quicker lookup than the vector.
    std::set<const Patch*> nonGhostPatches;
    for (size_t i = 0; i < superPatchGroup->size(); ++i) {
      nonGhostPatches.insert((*superPatchGroup)[i]);
    }

    Patch::selectType::iterator iter = encompassedPatches.begin();
    for (; iter != encompassedPatches.end(); ++iter) {
      const Patch* patchGroupMember = *iter;
      GridVariableBase* clone       = var.clone();
      IntVector groupMemberLowIndex =
        patchGroupMember->getExtraLowIndex(basis, label->getBoundaryLayer());
      IntVector groupMemberHighIndex =
        patchGroupMember->getExtraHighIndex(basis, label->getBoundaryLayer());
      IntVector enclosedLowIndex  = Max(groupMemberLowIndex, superLowIndex);
      IntVector enclosedHighIndex = Min(groupMemberHighIndex, superHighIndex);

      clone->rewindow(enclosedLowIndex, enclosedHighIndex);
      if (patchGroupMember == patch) {
        // this was checked already
        exists = false;
      } else {
        exists = d_var_DB.exists(label, matlIndex, patchGroupMember);
      }
      if (patchGroupMember->isVirtual()) {
        // Virtual patches can only be ghost patches.
        ASSERT(nonGhostPatches.find(patchGroupMember) == nonGhostPatches.end());
        clone->offsetGrid(IntVector(0, 0, 0) -
                          patchGroupMember->getVirtualOffset());
        enclosedLowIndex  = clone->getLow();
        enclosedHighIndex = clone->getHigh();
        patchGroupMember  = patchGroupMember->getRealPatch();
        IntVector dummy;
        if (d_scheduler->getSuperPatchExtents(label,
                                              matlIndex,
                                              patchGroupMember,
                                              gtype,
                                              numGhostCells,
                                              dummy,
                                              dummy,
                                              dummy,
                                              dummy) != 0) {
          // The virtual patch refers to a real patch in which the label
          // is computed locally, so don't overwrite the local copy.
          delete clone;
          continue;
        }
      }
      if (exists) {
        // variable section already exists in this patchGroupMember
        // (which is assumed to be a ghost patch)
        // so check if one is enclosed in the other.

        GridVariableBase* existingGhostVar = dynamic_cast<GridVariableBase*>(
          d_var_DB.get(label, matlIndex, patchGroupMember));
        IntVector existingLow  = existingGhostVar->getLow();
        IntVector existingHigh = existingGhostVar->getHigh();
        IntVector minLow       = Min(existingLow, enclosedLowIndex);
        IntVector maxHigh      = Max(existingHigh, enclosedHighIndex);

        if (existingGhostVar->isForeign()) {
          // data already being received, so don't replace it
          delete clone;
        } else if (minLow == enclosedLowIndex && maxHigh == enclosedHighIndex) {
          // this new ghost variable section encloses the old one,
          // so replace the old one
          printDebuggingPutInfo(label, matlIndex, patchGroupMember, __LINE__);
          d_var_DB.put(label,
                       matlIndex,
                       patchGroupMember,
                       clone,
                       d_scheduler->copyTimestep(),
                       true);
        } else {
          // Either the old ghost variable section encloses this new one
          // (so leave it), or neither encloses the other (so just forget
          // about it -- it'll allocate extra space for it when receiving
          // the ghost data in recvMPIGridVar if nothing else).
          delete clone;
        }
      } else {
        // it didn't exist before -- add it
        printDebuggingPutInfo(label, matlIndex, patchGroupMember, __LINE__);
        d_var_DB.put(label,
                     matlIndex,
                     patchGroupMember,
                     clone,
                     d_scheduler->copyTimestep(),
                     false);
      }
    }
  }
  d_lock.writeUnlock();
  var.rewindow(lowIndex, highIndex);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::copyOut(GridVariableBase& var,
                               const VarLabel* label,
                               int matlIndex,
                               const Patch* patch,
                               Ghost::GhostType gtype,
                               int numGhostCells)
{
  GridVariableBase* tmpVar = var.cloneType();
  getGridVar(*tmpVar, label, matlIndex, patch, gtype, numGhostCells);
  var.copyData(tmpVar);
  delete tmpVar;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::getCopy(GridVariableBase& var,
                               const VarLabel* label,
                               int matlIndex,
                               const Patch* patch,
                               Ghost::GhostType gtype,
                               int numGhostCells)
{
  GridVariableBase* tmpVar = var.cloneType();
  getGridVar(*tmpVar, label, matlIndex, patch, gtype, numGhostCells);
  var.allocate(tmpVar);
  var.copyData(tmpVar);
  delete tmpVar;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::put(GridVariableBase& var,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch,
                           bool replace /*= false*/)
{
  ASSERT(!d_finalized);
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);
  ASSERTEQ(
    basis,
    Patch::translateTypeToBasis(var.virtualGetTypeDescription()->getType(),
                                true));

  checkPutAccess(label, matlIndex, patch, replace);

#if DAV_DEBUG
  std::cerr << "Putting: " << *label << " MI: " << matlIndex
            << " patch: " << *patch << " into DW: " << d_generation << "\n";
#endif
  // Put it in the database
  IntVector low  = patch->getExtraLowIndex(basis, label->getBoundaryLayer());
  IntVector high = patch->getExtraHighIndex(basis, label->getBoundaryLayer());
  if (Min(var.getLow(), low) != var.getLow() ||
      Max(var.getHigh(), high) != var.getHigh()) {
    ostringstream msg_str;
    msg_str << "put: Variable's window (" << var.getLow() << " - "
            << var.getHigh() << ") must encompass patches extent (" << low
            << " - " << high;
    SCI_THROW(InternalError(msg_str.str(), __FILE__, __LINE__));
  }
  USE_IF_ASSERTS_ON(bool no_realloc =) var.rewindow(low, high);
  // error would have been thrown above if the any reallocation would be
  // needed
  ASSERT(no_realloc);
  d_var_DB.put(label,
               matlIndex,
               patch,
               var.clone(),
               d_scheduler->copyTimestep(),
               true);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::get(PerPatchBase& var,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch)
{
  checkGetAccess(label, matlIndex, patch);
  if (!d_var_DB.exists(label, matlIndex, patch)) {
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              patch,
                              matlIndex,
                              "perpatch data",
                              __FILE__,
                              __LINE__));
  }
  d_var_DB.get(label, matlIndex, patch, var);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::put(PerPatchBase& var,
                           const VarLabel* label,
                           int matlIndex,
                           const Patch* patch,
                           bool replace /*= false*/)
{
  ASSERT(!d_finalized);
  checkPutAccess(label, matlIndex, patch, replace);

  // Put it in the database
  d_var_DB.put(label,
               matlIndex,
               patch,
               var.clone(),
               d_scheduler->copyTimestep(),
               true);
}

//______________________________________________________________________
// This returns a constGridVariable for *ALL* patches on a level.
// This method is essentially identical to "getRegion" except the call to
// level->selectPatches( ) has been replaced by level->allPactches()
// For grids containing a large number of patches selectPatches() is very slow
// This assumes that the variable is not in the DWDatabase<Level>  d_level_DB;
//______________________________________________________________________
void
OnDemandDataWarehouse::getLevel(constGridVariableBase& constGridVar,
                                const VarLabel* label,
                                int matlIndex,
                                const Level* level)
{

  IntVector level_lowIndex, level_highIndex;
  level->findCellIndexRange(level_lowIndex,
                            level_highIndex); // including extra cells

  GridVariableBase* gridVar = constGridVar.cloneType();
  gridVar->allocate(level_lowIndex, level_highIndex);
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);

  std::vector<const Patch*> missing_patches; // for bulletproofing

  //__________________________________
  // define the patches for the entire level
  const PatchSet* myPatchesSet = level->allPatches();
  std::vector<const Patch*> patches(level->numPatches());
  for (int m = 0; m < myPatchesSet->size(); m++) {
    const PatchSubset* myPatches = myPatchesSet->getSubset(m);

    for (int p = 0; p < myPatches->size(); p++) {
      patches[p] = myPatches->get(p);
    }
  }

  int totalCells = 0;

  for (size_t i = 0; i < patches.size(); i++) {
    const Patch* patch = patches[i];

    std::vector<Variable*> varlist;
    d_var_DB.getlist(label, matlIndex, patch, varlist);
    GridVariableBase* this_var = nullptr;

    //__________________________________
    //  is this variable on this patch?
    for (std::vector<Variable*>::iterator rit = varlist.begin();; ++rit) {
      if (rit == varlist.end()) {
        this_var = nullptr;
        break;
      }

      // verify that the variable is valid
      this_var = dynamic_cast<GridVariableBase*>(*rit);

      if ((this_var != nullptr) && this_var->isValid()) {
        break;
      }
    }

    // just like a "missing patch": got data on this patch, but it either
    // corresponds to a different region or is incomplete"
    if (this_var == nullptr) {
      missing_patches.push_back(patch->getRealPatch());
      continue;
    }

    GridVariableBase* tmpVar = gridVar->cloneType();
    tmpVar->copyPointer(*this_var);

    // if patch is virtual, it is probably a boundary layer/extra cell that has
    // been requested (from AMR)
    if (patch->isVirtual()) {
      tmpVar->offset(patch->getVirtualOffset());
    }

    //__________________________________
    //  copy this patch's data
    IntVector lo = patch->getExtraLowIndex(basis, label->getBoundaryLayer());
    IntVector hi = patch->getExtraHighIndex(basis, label->getBoundaryLayer());

    try {
      gridVar->copyPatch(tmpVar, lo, hi);
    } catch (InternalError& e) {
      std::cout << " getLevel(" << label->getName()
                << "): copyPatch Bad range.\n"
                << " actual patch lo:" << lo << " hi:" << hi
                << " variable range: " << tmpVar->getLow() << " "
                << tmpVar->getHigh() << std::endl;
      throw e;
    }

    delete tmpVar;
    IntVector diff(hi - lo);
    totalCells += diff.x() * diff.y() * diff.z();
  } // patches loop

  [[maybe_unused]] long totalLevelCells = level->totalCells();

#ifdef BULLETPROOFING_FOR_CUBIC_DOMAINS
  //__________________________________
  //  This is not a valid check on non-cubic domains
  if (totalLevelCells != totalCells && missing_patches.size() > 0) {
    std::cout << d_myworld->myRank() << "  Unknown Variable " << *label
              << ", matl " << matlIndex << ", L-" << level->getIndex()
              << ", for patch(es): ";

    for (size_t i = 0; i < missing_patches.size(); i++) {
      std::cout << *missing_patches[i] << " ";
    }
    std::cout << " copied cells: " << totalCells
              << " requested cells: " << totalLevelCells << std::endl;
    std::cout << "  *** If the computational domain is non-cubic, (L-shaped or "
                 "necked down)  you can remove this bulletproofing"
              << std::endl;
    throw InternalError("Missing patches in getRegion", __FILE__, __LINE__);
  }
#endif

  //__________________________________
  //  Diagnostics
  if (dbg.active()) {
    cerrLock.lock();
    dbg << d_myworld->myRank() << "getLevel:  Variable " << *label << ", matl "
        << matlIndex << ", L-" << level->getIndex() << std::endl;
    cerrLock.unlock();
  }

  ASSERT(totalLevelCells <= totalCells);

  constGridVar = *dynamic_cast<GridVariableBase*>(gridVar);
  delete gridVar;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::getRegion(constGridVariableBase& constVar,
                                 const VarLabel* label,
                                 int matlIndex,
                                 const Level* level,
                                 const IntVector& low,
                                 const IntVector& high,
                                 bool useBoundaryCells /*=true*/)
{

  GridVariableBase* var = constVar.cloneType();
  var->allocate(low, high);
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);

  IntVector adjustment = IntVector(1, 1, 1);
  if (basis == Patch::XFaceBased) {
    adjustment = IntVector(1, 0, 0);
  } else if (basis == Patch::YFaceBased) {
    adjustment = IntVector(0, 1, 0);
  } else if (basis == Patch::ZFaceBased) {
    adjustment = IntVector(0, 0, 1);
  }

  Patch::selectType patches;

  // if in AMR and one node intersects from another patch and that patch is
  // missing ignore the error instead of throwing an exception (should only be
  // for node-based vars)
  vector<const Patch*> missing_patches;

  // make sure we grab all the patches, sometimes we might call only with an
  // extra cell region, which selectPatches doesn't detect
  IntVector tmpLow(low - adjustment);
  IntVector tmpHigh(high + adjustment);
  level->selectPatches(tmpLow, tmpHigh, patches);

  int totalCells = 0;
  for (int i = 0; i < patches.size(); i++) {
    const Patch* patch = patches[i];
    IntVector l, h;

    // the caller should determine whether or not he wants extra cells.
    // It will matter in AMR cases with corner-aligned patches
    if (useBoundaryCells) {
      l = Max(patch->getExtraLowIndex(basis, label->getBoundaryLayer()), low);
      h = Min(patch->getExtraHighIndex(basis, label->getBoundaryLayer()), high);
    } else {
      l = Max(patch->getLowIndex(basis), low);
      h = Min(patch->getHighIndex(basis), high);
    }

    if (l.x() >= h.x() || l.y() >= h.y() || l.z() >= h.z()) {
      continue;
    }
    vector<Variable*> varlist;
    d_var_DB.getlist(label, matlIndex, patch, varlist);
    GridVariableBase* v = nullptr;

    for (auto rit = varlist.rbegin();; ++rit) {
      if (rit == varlist.rend()) {
        v = nullptr;
        break;
      }
      v = dynamic_cast<GridVariableBase*>(*rit);
      // verify that the variable is valid and matches the dependencies
      // requirements.
      if ((v != nullptr) && v->isValid() &&
          Min(l, v->getLow()) == v->getLow() &&
          Max(h, v->getHigh()) == v->getHigh()) { // find a completed region
        break;
      }
    }

    // just like a "missing patch": got data on this patch, but it either
    // corresponds to a different region or is incomplete"
    if (v == nullptr) {
      missing_patches.push_back(patch->getRealPatch());
      continue;
    }

    GridVariableBase* tmpVar = var->cloneType();
    tmpVar->copyPointer(*v);

    if (patch->isVirtual()) {
      // if patch is virtual, it is probable a boundary layer/extra cell that
      // has been requested (from AMR) let Bryan know if this doesn't work.  We
      // need to adjust the source but not the dest by the virtual offset
      tmpVar->offset(patch->getVirtualOffset());
    }
    try {
      var->copyPatch(tmpVar, l, h);
    } catch (InternalError& e) {
      std::cout << " Bad range: " << low << " " << high
                << ", patch intersection: " << l << " " << h << " actual patch "
                << patch->getLowIndex(basis) << " "
                << patch->getHighIndex(basis)
                << " var range: " << tmpVar->getLow() << " "
                << tmpVar->getHigh() << endl;
      throw e;
    }
    delete tmpVar;
    IntVector diff(h - l);
    totalCells += diff.x() * diff.y() * diff.z();
  } // patches loop

  IntVector diff(high - low);

#ifdef BULLETPROOFING_FOR_CUBIC_DOMAINS
  if (diff.x() * diff.y() * diff.z() > totalCells &&
      missing_patches.size() > 0) {
    std::cout << d_myworld->myRank() << "  Unknown Variable " << *label
              << ", matl " << matlIndex << ", L-" << level->getIndex()
              << ", for patch(es): ";
    for (size_t i = 0; i < missing_patches.size(); i++) {
      std::cout << *missing_patches[i] << " ";
    }
    std::cout << endl << " Original region: " << low << " " << high << endl;
    std::cout << " copied cells: " << totalCells
              << " requested cells: " << diff.x() * diff.y() * diff.z() << endl;
    throw InternalError("Missing patches in getRegion", __FILE__, __LINE__);
  }
#endif
  if (dbg.active()) {
    cerrLock.lock();
    dbg << d_myworld->myRank() << "  Variable " << *label << ", matl "
        << matlIndex << ", L-" << level->getIndex() << " For region: " << low
        << " " << high << "  has been gotten" << endl;
    cerrLock.unlock();
  }

  ASSERT(diff.x() * diff.y() * diff.z() <= totalCells);

  constVar = *dynamic_cast<GridVariableBase*>(var);
  delete var;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::getRegion(GridVariableBase& var,
                                 const VarLabel* label,
                                 int matlIndex,
                                 const Level* level,
                                 const IntVector& low,
                                 const IntVector& high,
                                 bool useBoundaryCells /*=true*/)
{

  var.allocate(low, high);
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);

  IntVector adjustment = IntVector(1, 1, 1);
  if (basis == Patch::XFaceBased) {
    adjustment = IntVector(1, 0, 0);
  } else if (basis == Patch::YFaceBased) {
    adjustment = IntVector(0, 1, 0);
  } else if (basis == Patch::ZFaceBased) {
    adjustment = IntVector(0, 0, 1);
  }

  Patch::selectType patches;

  // if in AMR and one node intersects from another patch and that patch is
  // missing ignore the error instead of throwing an exception (should only be
  // for node-based vars)
  std::vector<const Patch*> missing_patches;

  // make sure we grab all the patches, sometimes we might call only with an
  // extra cell region, which selectPatches doesn't detect
  IntVector tmpLow(low - adjustment);
  IntVector tmpHigh(high + adjustment);
  level->selectPatches(tmpLow, tmpHigh, patches);

  int totalCells = 0;
  for (int i = 0; i < patches.size(); i++) {
    const Patch* patch = patches[i];
    IntVector l, h;

    // the caller should determine whether or not he wants extra cells.
    // It will matter in AMR cases with corner-aligned patches
    if (useBoundaryCells) {
      l = Max(patch->getExtraLowIndex(basis, label->getBoundaryLayer()), low);
      h = Min(patch->getExtraHighIndex(basis, label->getBoundaryLayer()), high);
    } else {
      l = Max(patch->getLowIndex(basis), low);
      h = Min(patch->getHighIndex(basis), high);
    }

    if (l.x() >= h.x() || l.y() >= h.y() || l.z() >= h.z()) {
      continue;
    }
    std::vector<Variable*> varlist;
    d_var_DB.getlist(label, matlIndex, patch, varlist);
    GridVariableBase* v = nullptr;

    for (auto rit = varlist.begin();; ++rit) {
      if (rit == varlist.end()) {
        v = nullptr;
        break;
      }
      v = dynamic_cast<GridVariableBase*>(*rit);
      // verify that the variable is valid and matches the dependencies
      // requirements.
      if ((v != nullptr) && v->isValid() &&
          Min(l, v->getLow()) == v->getLow() &&
          Max(h, v->getHigh()) == v->getHigh()) { // find a completed region
        break;
      }
    }

    // just like a "missing patch": got data on this patch, but it either
    // corresponds to a different region or is incomplete"
    if (v == nullptr) {
      missing_patches.push_back(patch->getRealPatch());
      continue;
    }

    GridVariableBase* tmpVar = var.cloneType();
    tmpVar->copyPointer(*v);

    if (patch->isVirtual()) {
      // if patch is virtual, it is probable a boundary layer/extra cell that
      // has been requested (from AMR) let Bryan know if this doesn't work.  We
      // need to adjust the source but not the dest by the virtual offset
      tmpVar->offset(patch->getVirtualOffset());
    }
    try {
      var.copyPatch(tmpVar, l, h);
    } catch (InternalError& e) {
      std::cout << " Bad range: " << low << " " << high
                << ", patch intersection: " << l << " " << h << " actual patch "
                << patch->getLowIndex(basis) << " "
                << patch->getHighIndex(basis)
                << " var range: " << tmpVar->getLow() << " "
                << tmpVar->getHigh() << std::endl;
      throw e;
    }
    delete tmpVar;
    IntVector diff(h - l);
    totalCells += diff.x() * diff.y() * diff.z();
  } // patches loop

  IntVector diff(high - low);

#ifdef BULLETPROOFING_FOR_CUBIC_DOMAINS
  //__________________________________
  //  Not valid for non-cubic domains
  if (diff.x() * diff.y() * diff.z() > totalCells &&
      missing_patches.size() > 0) {
    std::cout << d_myworld->myRank() << "  Unknown Variable " << *label
              << ", matl " << matlIndex << ", L-" << level->getIndex()
              << ", for patch(es): ";
    for (size_t i = 0; i < missing_patches.size(); i++) {
      std::cout << *missing_patches[i] << " ";
    }
    std::cout << std::endl
              << " Original region: " << low << " " << high << std::endl;
    std::cout << " copied cells: " << totalCells
              << " requested cells: " << diff.x() * diff.y() * diff.z()
              << std::endl;
    throw InternalError("Missing patches in getRegion", __FILE__, __LINE__);
  }
#endif
  if (dbg.active()) {
    coutLock.lock();
    dbg << d_myworld->myRank() << "  Variable " << *label << ", matl "
        << matlIndex << ", L-" << level->getIndex() << " For region: " << low
        << " " << high << "  has been gotten" << std::endl;
    coutLock.unlock();
  }

  ASSERT(diff.x() * diff.y() * diff.z() <= totalCells);

  //  constVar = *dynamic_cast<GridVariableBase*>(var);
  //  delete var;
}

//______________________________________________________________________
//
void
OnDemandDataWarehouse::emit(OutputContext& oc,
                            const VarLabel* label,
                            int matlIndex,
                            const Patch* patch)
{
  checkGetAccess(label, matlIndex, patch);

  Variable* var = nullptr;
  IntVector l, h;
  if (patch) {
    // Save with the boundary layer, otherwise restarting from the DataArchive
    // won't work.
    patch->computeVariableExtents(label->typeDescription()->getType(),
                                  label->getBoundaryLayer(),
                                  Ghost::None,
                                  0,
                                  l,
                                  h);
    switch (label->typeDescription()->getType()) {
      case TypeDescription::Type::NCVariable:
      case TypeDescription::Type::CCVariable:
      case TypeDescription::Type::SFCXVariable:
      case TypeDescription::Type::SFCYVariable:
      case TypeDescription::Type::SFCZVariable:
        // get list
        {
          vector<Variable*> varlist;
          d_var_DB.getlist(label, matlIndex, patch, varlist);

          GridVariableBase* v = nullptr;
          for (auto rit = varlist.begin();; ++rit) {
            if (rit == varlist.end()) {
              v = nullptr;
              break;
            }
            v = dynamic_cast<GridVariableBase*>(*rit);
            // verify that the variable is valid and matches the dependencies
            // requirements.
            if (v && v->isValid() && Min(l, v->getLow()) == v->getLow() &&
                Max(h, v->getHigh()) ==
                  v->getHigh()) { // find a completed region
              break;
            }
          }
          var = v;
        }
        break;
      case TypeDescription::Type::ParticleVariable:
        var = d_var_DB.get(label, matlIndex, patch);
        break;
      default:
        var = d_var_DB.get(label, matlIndex, patch);
    }
  } else {
    l = h = IntVector(-1, -1, -1);

    const Level* level = patch ? patch->getLevel() : 0;
    if (d_level_DB.exists(label, matlIndex, level)) {
      var = d_level_DB.get(label, matlIndex, level);
    }
  }

  if (var == nullptr) {
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              patch,
                              matlIndex,
                              "on emit",
                              __FILE__,
                              __LINE__));
  }
  var->emit(oc, l, h, label->getCompressionMode());
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::print(ostream& intout,
                             const VarLabel* label,
                             const Level* level,
                             int matlIndex /* = -1 */)
{

  try {
    checkGetAccess(label, matlIndex, 0);
    ReductionVariableBase* var = dynamic_cast<ReductionVariableBase*>(
      d_level_DB.get(label, matlIndex, level));
    var->print(intout);
  } catch (const UnknownVariable& e) {
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              level,
                              matlIndex,
                              "on emit reduction",
                              __FILE__,
                              __LINE__));
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::deleteParticles(ParticleSubset* delset)
{
  int matlIndex          = delset->getMatlIndex();
  Patch* patch           = (Patch*)delset->getPatch();
  const Patch* realPatch = (patch != 0) ? patch->getRealPatch() : 0;

  d_pslock.writeLock();
  psetDBType::key_type key(realPatch, matlIndex, getID());
  psetDBType::iterator iter = d_delset_DB.find(key);
  ParticleSubset* currentDelset;
  if (iter != d_delset_DB.end()) { // update existing delset
    //    SCI_THROW(InternalError("deleteParticles called twice for patch",
    //    __FILE__, __LINE__));
    // Concatenate the delsets into the delset that already exists in the DB.
    currentDelset = iter->second;
    for (ParticleSubset::iterator d = delset->begin(); d != delset->end();
         d++) {
      currentDelset->addParticle(*d);
    }

    d_delset_DB.erase(key);
    d_delset_DB.insert(
      std::pair<psetDBType::key_type, ParticleSubset*>(key, currentDelset));

    delete delset;

  } else {
    d_delset_DB.insert(
      std::pair<psetDBType::key_type, ParticleSubset*>(key, delset));
    delset->addReference();
  }
  d_pslock.writeUnlock();
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::addParticles(const Patch* patch,
                                    int matlIndex,
                                    ParticleLabelVariableMap* addedState)
{
  d_pslock.writeLock();
  psetAddDBType::key_type key(matlIndex, patch);
  psetAddDBType::iterator iter = d_addset_DB.find(key);
  if (iter != d_addset_DB.end()) {
    // SCI_THROW(InternalError("addParticles called twice for patch", __FILE__,
    // __LINE__));
    std::cerr << "addParticles called twice for patch" << endl;
  }

  else {
    d_addset_DB[key] = addedState;
  }

  d_pslock.writeUnlock();
}
//______________________________________________________________________
//
int
OnDemandDataWarehouse::decrementScrubCount(const VarLabel* var,
                                           int matlIndex,
                                           const Patch* patch)
{

  int count = 0;
  switch (var->typeDescription()->getType()) {
    case TypeDescription::Type::NCVariable:
    case TypeDescription::Type::CCVariable:
    case TypeDescription::Type::SFCXVariable:
    case TypeDescription::Type::SFCYVariable:
    case TypeDescription::Type::SFCZVariable:
    case TypeDescription::Type::PerPatch:
      // try {
      count = d_var_DB.decrementScrubCount(var, matlIndex, patch);
      //}
      // catch (AssertionFailed& e) {
      // std::cout << d_myworld->myRank() << " DW " << getID() << " caught
      // exception.\n"; throw e;
      //}
      break;
    case TypeDescription::Type::ParticleVariable:
      count = d_var_DB.decrementScrubCount(var, matlIndex, patch);
      break;
    case TypeDescription::Type::SoleVariable:
      SCI_THROW(InternalError("decrementScrubCount called for sole variable: " +
                                var->getName(),
                              __FILE__,
                              __LINE__));
    case TypeDescription::Type::ReductionVariable:
      SCI_THROW(InternalError(
        "decrementScrubCount called for reduction variable: " + var->getName(),
        __FILE__,
        __LINE__));
    default:
      SCI_THROW(InternalError(
        "decrementScrubCount for variable of unknown type: " + var->getName(),
        __FILE__,
        __LINE__));
  }
  return count;
}
//______________________________________________________________________
//
DataWarehouse::ScrubMode
OnDemandDataWarehouse::setScrubbing(ScrubMode scrubMode)
{
  ScrubMode oldmode = d_scrubMode;
  d_scrubMode       = scrubMode;
  return oldmode;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::setScrubCount(const VarLabel* var,
                                     int matlIndex,
                                     const Patch* patch,
                                     int count)
{
  switch (var->typeDescription()->getType()) {
    case TypeDescription::Type::NCVariable:
    case TypeDescription::Type::CCVariable:
    case TypeDescription::Type::SFCXVariable:
    case TypeDescription::Type::SFCYVariable:
    case TypeDescription::Type::SFCZVariable:
    case TypeDescription::Type::PerPatch:
    case TypeDescription::Type::ParticleVariable:
      d_var_DB.setScrubCount(var, matlIndex, patch, count);
      break;
    case TypeDescription::Type::SoleVariable:
      SCI_THROW(InternalError("setScrubCount called for sole variable: " +
                                var->getName(),
                              __FILE__,
                              __LINE__));
    case TypeDescription::Type::ReductionVariable:
      // Reductions are not scrubbed
      SCI_THROW(InternalError("setScrubCount called for reduction variable: " +
                                var->getName(),
                              __FILE__,
                              __LINE__));
    default:
      SCI_THROW(InternalError("setScrubCount for variable of unknown type: " +
                                var->getName(),
                              __FILE__,
                              __LINE__));
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::scrub(const VarLabel* var,
                             int matlIndex,
                             const Patch* patch)
{
  switch (var->typeDescription()->getType()) {
    case TypeDescription::Type::NCVariable:
    case TypeDescription::Type::CCVariable:
    case TypeDescription::Type::SFCXVariable:
    case TypeDescription::Type::SFCYVariable:
    case TypeDescription::Type::SFCZVariable:
    case TypeDescription::Type::PerPatch:
    case TypeDescription::Type::ParticleVariable:
      d_var_DB.scrub(var, matlIndex, patch);
      break;
    case TypeDescription::Type::SoleVariable:
      SCI_THROW(
        InternalError("scrub called for sole variable: " + var->getName(),
                      __FILE__,
                      __LINE__));
    case TypeDescription::Type::ReductionVariable:
      // Reductions are not scrubbed
      SCI_THROW(
        InternalError("scrub called for reduction variable: " + var->getName(),
                      __FILE__,
                      __LINE__));
    default:
      SCI_THROW(
        InternalError("scrub for variable of unknown type: " + var->getName(),
                      __FILE__,
                      __LINE__));
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::initializeScrubs(
  int dwid,
  const FastHashTable<ScrubItem>* scrubcounts,
  bool add)
{
  d_var_DB.initializeScrubs(dwid, scrubcounts, add);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::getGridVar(GridVariableBase& var,
                                  const VarLabel* label,
                                  int matlIndex,
                                  const Patch* patch,
                                  Ghost::GhostType gtype,
                                  int numGhostCells)
{
  Patch::VariableBasis basis =
    Patch::translateTypeToBasis(label->typeDescription()->getType(), false);
  ASSERTEQ(
    basis,
    Patch::translateTypeToBasis(var.virtualGetTypeDescription()->getType(),
                                true));

  if (!d_var_DB.exists(label, matlIndex, patch)) {
    // print();
    std::cout << d_myworld->myRank() << " unable to find variable '"
              << label->getName() << " on patch: " << patch->getID()
              << " matl: " << matlIndex << endl;
    // WAIT_FOR_DEBUGGER();
    SCI_THROW(UnknownVariable(label->getName(),
                              getID(),
                              patch,
                              matlIndex,
                              "",
                              __FILE__,
                              __LINE__));
  }
  if (patch->isVirtual()) {
    d_var_DB.get(label, matlIndex, patch->getRealPatch(), var);
    var.offsetGrid(patch->getVirtualOffset());
  } else {
    d_var_DB.get(label, matlIndex, patch, var);
  }

  IntVector low  = patch->getExtraLowIndex(basis, label->getBoundaryLayer());
  IntVector high = patch->getExtraHighIndex(basis, label->getBoundaryLayer());

  if (gtype == Ghost::None) {
    if (numGhostCells != 0) {
      SCI_THROW(InternalError("Ghost cells specified with type: None!\n",
                              __FILE__,
                              __LINE__));
    }
    // if this assertion fails, then it is having problems getting the
    // correct window of the data.
    USE_IF_ASSERTS_ON(bool no_realloc =) var.rewindow(low, high);
    ASSERT(no_realloc);
  } else {
    IntVector dn = high - low;
    long total   = dn.x() * dn.y() * dn.z();

    Patch::selectType neighbors;
    IntVector lowIndex, highIndex;
    patch->computeVariableExtents(basis,
                                  label->getBoundaryLayer(),
                                  gtype,
                                  numGhostCells,
                                  lowIndex,
                                  highIndex);

    if (numGhostCells > 0) {
      patch->getLevel()->selectPatches(lowIndex, highIndex, neighbors);
    } else {
      neighbors.push_back(patch);
    }
    IntVector oldLow = var.getLow(), oldHigh = var.getHigh();
    if (!var.rewindow(lowIndex, highIndex)) {
      // reallocation needed
      // Ignore this if this is the initialization dw in its old state.
      // The reason for this is that during initialization it doesn't
      // know what ghost cells will be required of it for the next timestep.
      // (This will be an issue whenever the taskgraph changes to require
      // more ghost cells from the old datawarehouse).
      static bool warned = false;
      bool ignore        = d_isInitializationDW && d_finalized;
      if (!ignore && !warned) {
        // warned = true;
        // static ProgressiveWarning rw("Warning: Reallocation needed for ghost
        // region you requested.\nThis means the data you get back will be a
        // copy of what's in the DW", 100); if (rw.invoke()) {
        //  print out this message if the ProgressiveWarning does
        /*ostringstream errmsg;
          errmsg << d_myworld->myRank() << " This occurrence for " <<
          label->getName(); if (patch) errmsg << " on patch " << patch->getID();
          errmsg << " for material " << matlIndex;

          errmsg << ".  Old range: " << oldLow << " " << oldHigh << " - new
          range " << lowIndex << " " << highIndex << " NGC " << numGhostCells;
          warn << errmsg.str() << '\n';*/

        //}
      }
    }

    for (int i = 0; i < neighbors.size(); i++) {
      const Patch* neighbor = neighbors[i];
      if (neighbor && (neighbor != patch)) {
        IntVector low =
          Max(neighbor->getExtraLowIndex(basis, label->getBoundaryLayer()),
              lowIndex);
        IntVector high =
          Min(neighbor->getExtraHighIndex(basis, label->getBoundaryLayer()),
              highIndex);

        if (patch->getLevel()->getIndex() > 0 && patch != neighbor) {
          patch->cullIntersection(basis,
                                  label->getBoundaryLayer(),
                                  neighbor,
                                  low,
                                  high);
        }
        if (low == high) {
          continue;
        }

        if (!d_var_DB.exists(label, matlIndex, neighbor)) {
          SCI_THROW(
            UnknownVariable(label->getName(),
                            getID(),
                            neighbor,
                            matlIndex,
                            neighbor == patch ? "on patch" : "on neighbor",
                            __FILE__,
                            __LINE__));
        }

        vector<Variable*> varlist;
        d_var_DB.getlist(label, matlIndex, neighbor, varlist);
        GridVariableBase* v = nullptr;

        for (auto it = varlist.begin();; ++it) {
          if (it == varlist.end()) {
            v = nullptr;
            break;
          }

          v = dynamic_cast<GridVariableBase*>(*it);
          // verify that the variable is valid and matches the depedencies
          // requirements
          if ((v != nullptr) && (v->isValid())) {
            if (neighbor->isVirtual()) {
              if (Min(v->getLow(), low - neighbor->getVirtualOffset()) ==
                    v->getLow() &&
                  Max(v->getHigh(), high - neighbor->getVirtualOffset()) ==
                    v->getHigh()) {
                break;
              }
            } else {
              if (Min(v->getLow(), low) == v->getLow() &&
                  Max(v->getHigh(), high) == v->getHigh()) {
                break;
              }
            }
          }
        } // end for vars
        if (v == nullptr) {
          // std::cout << d_myworld->myRank()  << " cannot copy var " << *label
          // << " from patch " << neighbor->getID()
          // << " " << low << " " << high <<  ", DW has " << srcvar->getLow() <<
          // " " << srcvar->getHigh() << endl;
          SCI_THROW(
            UnknownVariable(label->getName(),
                            getID(),
                            neighbor,
                            matlIndex,
                            neighbor == patch ? "on patch" : "on neighbor",
                            __FILE__,
                            __LINE__));
        }

        GridVariableBase* srcvar = var.cloneType();
        srcvar->copyPointer(*v);

        if (neighbor->isVirtual()) {
          srcvar->offsetGrid(neighbor->getVirtualOffset());
        }

        try {
          var.copyPatch(srcvar, low, high);
        } catch (InternalError& e) {
          std::cout << " Bad range: " << low << " " << high
                    << " source var range: " << srcvar->getLow() << " "
                    << srcvar->getHigh() << endl;
          throw e;
        }
        delete srcvar;
        dn = high - low;
        total += dn.x() * dn.y() * dn.z();
      } // end if neigbor
    }   // end for neigbours

    // dn = highIndex - lowIndex;
    // long wanted = dn.x()*dn.y()*dn.z();
    // ASSERTEQ(wanted, total);
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::transferFrom(DataWarehouse* from,
                                    const VarLabel* var,
                                    const PatchSubset* patches,
                                    const MaterialSubset* matls,
                                    bool replace /*=false*/,
                                    const PatchSubset* newPatches /*=0*/)
{
  OnDemandDataWarehouse* fromDW = dynamic_cast<OnDemandDataWarehouse*>(from);
  ASSERT(fromDW != 0);
  ASSERT(!d_finalized);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch     = patches->get(p);
    const Patch* copyPatch = (newPatches ? newPatches->get(p) : patch);
    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      checkPutAccess(var, matl, patch, replace);
      switch (var->typeDescription()->getType()) {
        case TypeDescription::Type::NCVariable:
        case TypeDescription::Type::CCVariable:
        case TypeDescription::Type::SFCXVariable:
        case TypeDescription::Type::SFCYVariable:
        case TypeDescription::Type::SFCZVariable: {
          if (!fromDW->d_var_DB.exists(var, matl, patch)) {
            SCI_THROW(UnknownVariable(var->getName(),
                                      fromDW->getID(),
                                      patch,
                                      matl,
                                      "in transferFrom",
                                      __FILE__,
                                      __LINE__));
          }
          GridVariableBase* v = dynamic_cast<GridVariableBase*>(
                                  fromDW->d_var_DB.get(var, matl, patch))
                                  ->clone();
          d_var_DB
            .put(var, matl, copyPatch, v, d_scheduler->copyTimestep(), replace);
        } break;
        case TypeDescription::Type::ParticleVariable: {
          if (!fromDW->d_var_DB.exists(var, matl, patch)) {
            SCI_THROW(UnknownVariable(var->getName(),
                                      getID(),
                                      patch,
                                      matl,
                                      "in transferFrom",
                                      __FILE__,
                                      __LINE__));
          }

          // or else the readLock in haveParticleSubset will hang
          ParticleSubset* subset;
          if (!haveParticleSubset(matl, copyPatch)) {
            ParticleSubset* oldsubset = fromDW->getParticleSubset(matl, patch);
            subset =
              createParticleSubset(oldsubset->numParticles(), matl, copyPatch);
          } else {
            subset = getParticleSubset(matl, copyPatch);
          }
          ParticleVariableBase* v = dynamic_cast<ParticleVariableBase*>(
            fromDW->d_var_DB.get(var, matl, patch));
          if (patch == copyPatch) {
            d_var_DB.put(var,
                         matl,
                         copyPatch,
                         v->clone(),
                         d_scheduler->copyTimestep(),
                         replace);
          } else {
            ParticleVariableBase* newv = v->cloneType();
            newv->copyPointer(*v);
            newv->setParticleSubset(subset);
            d_var_DB.put(var,
                         matl,
                         copyPatch,
                         newv,
                         d_scheduler->copyTimestep(),
                         replace);
          }
        } break;
        case TypeDescription::Type::PerPatch: {
          if (!fromDW->d_var_DB.exists(var, matl, patch)) {
            SCI_THROW(UnknownVariable(var->getName(),
                                      getID(),
                                      patch,
                                      matl,
                                      "in transferFrom",
                                      __FILE__,
                                      __LINE__));
          }
          PerPatchBase* v =
            dynamic_cast<PerPatchBase*>(fromDW->d_var_DB.get(var, matl, patch));
          d_var_DB.put(var,
                       matl,
                       copyPatch,
                       v->clone(),
                       d_scheduler->copyTimestep(),
                       replace);
        } break;
        case TypeDescription::Type::ReductionVariable:
          SCI_THROW(
            InternalError("transferFrom doesn't work for reduction variable: " +
                            var->getName(),
                          __FILE__,
                          __LINE__));
          break;
        case TypeDescription::Type::SoleVariable:
          SCI_THROW(InternalError(
            "transferFrom doesn't work for sole variable: " + var->getName(),
            __FILE__,
            __LINE__));
          break;
        default:
          SCI_THROW(InternalError("Unknown variable type in transferFrom: " +
                                    var->getName(),
                                  __FILE__,
                                  __LINE__));
      }
    }
  }
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::logMemoryUse(ostream& out,
                                    unsigned long& total,
                                    const std::string& tag)
{
  int dwid = d_generation;
  d_var_DB.logMemoryUse(out, total, tag, dwid);

  // Log the psets.
  for (psetDBType::iterator iter = d_pset_DB.begin(); iter != d_pset_DB.end();
       iter++) {
    ParticleSubset* pset = iter->second;
    ostringstream elems;
    elems << pset->numParticles();
    logMemory(out,
              total,
              tag,
              "particles",
              "ParticleSubset",
              pset->getPatch(),
              pset->getMatlIndex(),
              elems.str(),
              pset->numParticles() * sizeof(particleIndex),
              pset->getPointer(),
              dwid);
  }
}
//______________________________________________________________________
//
inline void
OnDemandDataWarehouse::checkGetAccess(const VarLabel* label,
                                      int matlIndex,
                                      const Patch* patch,
                                      Ghost::GhostType gtype,
                                      int numGhostCells)
{
#if 1
#if SCI_ASSERTION_LEVEL >= 1
  std::list<RunningTaskInfo>* runningTasks = getRunningTasksInfo();

  if (runningTasks != 0) {
    for (std::list<RunningTaskInfo>::iterator iter = runningTasks->begin();
         iter != runningTasks->end();
         iter++) {
      RunningTaskInfo& runningTaskInfo = *iter;

      //   RunningTaskInfo& runningTaskInfo = runningTasks->back();
      const Task* runningTask = runningTaskInfo.d_task;
      if (runningTask == 0) {
        // don't check if done outside of any task (i.e. SimulationController)
        return;
      }

      IntVector lowOffset, highOffset;
      Patch::getGhostOffsets(label->typeDescription()->getType(),
                             gtype,
                             numGhostCells,
                             lowOffset,
                             highOffset);

      VarAccessMap& runningTaskAccesses = runningTaskInfo.d_accesses;

      map<VarLabelMatl<Patch>, AccessInfo>::iterator findIter;
      findIter =
        runningTaskAccesses.find(VarLabelMatl<Patch>(label, matlIndex, patch));

      if (!hasGetAccess(runningTask,
                        label,
                        matlIndex,
                        patch,
                        lowOffset,
                        highOffset,
                        &runningTaskInfo) &&
          !hasPutAccess(runningTask, label, matlIndex, patch, true) &&
          !hasPutAccess(runningTask, label, matlIndex, patch, false)) {

        // If it was accessed by the current task already, then it should
        // have get access (i.e. if you put it in, you should be able to get it
        // right back out).
        if (findIter != runningTaskAccesses.end() &&
            lowOffset == IntVector(0, 0, 0) &&
            highOffset == IntVector(0, 0, 0)) {
          // allow non ghost cell get if any access (get, put, or modify) is
          // allowed
          // std::cout << "allowing non-ghost cell access\n";
          return;
        }

        if (runningTask == 0 ||
            !(string(runningTask->getName()) == "Relocate::relocateParticles" ||
              string(runningTask->getName()) ==
                "SchedulerCommon::copyDataToNewGrid")) {
          string has;
          switch (getWhichDW(&runningTaskInfo)) {
            case Task::NewDW:
              has = "Task::NewDW";
              break;
            case Task::OldDW:
              has = "Task::OldDW";
              break;
            case Task::ParentNewDW:
              has = "Task::ParentNewDW";
              break;
            case Task::ParentOldDW:
              has = "Task::ParentOldDW";
              break;
            default:
              has = "UnknownDW";
          }
          has += " datawarehouse get";

          if (numGhostCells > 0) {
            ostringstream ghost_str;
            ghost_str << " for " << numGhostCells << " layer";

            if (numGhostCells > 1) {
              ghost_str << "s";
            }
            ghost_str << " of ghosts around " << Ghost::getGhostTypeName(gtype);
            has += ghost_str.str();
          }
          string needs = "task requires";
#if 1
          SCI_THROW(DependencyException(runningTask,
                                        label,
                                        matlIndex,
                                        patch,
                                        has,
                                        needs,
                                        __FILE__,
                                        __LINE__));
#else
          if (d_myworld->myRank() == 0) {
            std::cout << DependencyException::makeMessage(runningTask,
                                                          label,
                                                          matlIndex,
                                                          patch,
                                                          has,
                                                          needs)
                      << endl;
            // WAIT_FOR_DEBUGGER();
          }
#endif
        }
      } else {
        // access granted
        if (findIter == runningTaskAccesses.end()) {
          AccessInfo& accessInfo =
            runningTaskAccesses[VarLabelMatl<Patch>(label, matlIndex, patch)];
          accessInfo.accessType = GetAccess;
          accessInfo.encompassOffsets(lowOffset, highOffset);

          int ID = 0;
          if (patch) {
            ID = patch->getID();
          }
          string varname = "noname";
          if (label) {
            varname = label->getName();
          }
          if (dbg.active()) {
            cerrLock.lock();
            dbg << d_myworld->myRank()
                << " Task running is: " << runningTask->getName();
            dbg << std::left;
            dbg.width(10);
            dbg << "\t" << varname;
            dbg << std::left;
            dbg.width(10);
            dbg << " \t on patch " << ID << " and matl: " << matlIndex
                << " has been gotten\n";
            cerrLock.unlock();
          }
        } else {
          findIter->second.encompassOffsets(lowOffset, highOffset);
        }
      }
    }
  } // running task loop
#endif
#endif
}
//______________________________________________________________________
//
inline void
OnDemandDataWarehouse::checkPutAccess(const VarLabel* label,
                                      int matlIndex,
                                      const Patch* patch,
                                      bool replace)
{
#if 1
#if SCI_ASSERTION_LEVEL >= 1
  std::list<RunningTaskInfo>* runningTasks = getRunningTasksInfo();
  if (runningTasks != 0) {
    for (std::list<RunningTaskInfo>::iterator iter = runningTasks->begin();
         iter != runningTasks->end();
         iter++) {
      RunningTaskInfo& runningTaskInfo = *iter;
      const Task* runningTask          = runningTaskInfo.d_task;

      if (runningTask == 0) {
        return; // don't check if outside of any task (i.e.
                // SimulationController)
      }

      VarAccessMap& runningTaskAccesses = runningTaskInfo.d_accesses;

      if (!hasPutAccess(runningTask, label, matlIndex, patch, replace)) {
        if (string(runningTask->getName()) != "Relocate::relocateParticles") {
          string has, needs;
          switch (getWhichDW(&runningTaskInfo)) {
            case Task::NewDW:
              has = "Task::NewDW";
              break;
            case Task::OldDW:
              has = "Task::OldDW";
              break;
            case Task::ParentNewDW:
              has = "Task::ParentNewDW";
              break;
            case Task::ParentOldDW:
              has = "Task::ParentOldDW";
              break;
            default:
              has = "UnknownDW";
          }
          if (replace) {
            has += " datawarehouse put";
            needs = "task computes(replace)";
          } else {
            has += " datawarehouse put";
            needs = "task computes";
          }
#if 1
          SCI_THROW(DependencyException(runningTask,
                                        label,
                                        matlIndex,
                                        patch,
                                        has,
                                        needs,
                                        __FILE__,
                                        __LINE__));
#else
          if (d_myworld->myRank() == 0) {
            std::cout << DependencyException::makeMessage(runningTask,
                                                          label,
                                                          matlIndex,
                                                          patch,
                                                          has,
                                                          needs)
                      << endl;
          }
          // WAIT_FOR_DEBUGGER();
#endif
        }
      } else {
        runningTaskAccesses[VarLabelMatl<Patch>(label, matlIndex, patch)]
          .accessType = replace ? ModifyAccess : PutAccess;
      }
    }
  }
#endif
#endif
}

inline void
OnDemandDataWarehouse::checkModifyAccess(const VarLabel* label,
                                         int matlIndex,
                                         const Patch* patch)
{
  checkPutAccess(label, matlIndex, patch, true);
}

//______________________________________________________________________
//
inline Task::WhichDW
OnDemandDataWarehouse::getWhichDW(RunningTaskInfo* info)
{
  if (this == OnDemandDataWarehouse::getOtherDataWarehouse(Task::NewDW, info)) {
    return Task::NewDW;
  }
  if (this == OnDemandDataWarehouse::getOtherDataWarehouse(Task::OldDW, info)) {
    return Task::OldDW;
  }
  if (this ==
      OnDemandDataWarehouse::getOtherDataWarehouse(Task::ParentNewDW, info)) {
    return Task::ParentNewDW;
  }
  if (this ==
      OnDemandDataWarehouse::getOtherDataWarehouse(Task::ParentOldDW, info)) {
    return Task::ParentOldDW;
  }
  throw InternalError("Unknown DW\n", __FILE__, __LINE__);
}
//______________________________________________________________________
//
inline bool
OnDemandDataWarehouse::hasGetAccess(const Task* runningTask,
                                    const VarLabel* label,
                                    int matlIndex,
                                    const Patch* patch,
                                    IntVector lowOffset,
                                    IntVector highOffset,
                                    RunningTaskInfo* info)
{
  return runningTask->hasRequires(label,
                                  matlIndex,
                                  patch,
                                  lowOffset,
                                  highOffset,
                                  getWhichDW(info));
}
//______________________________________________________________________
//
inline bool
OnDemandDataWarehouse::hasPutAccess(const Task* runningTask,
                                    const VarLabel* label,
                                    int matlIndex,
                                    const Patch* patch,
                                    bool replace)
{
  return runningTask->hasComputes(label, matlIndex, patch);
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::pushRunningTask(const Task* task,
                                       vector<OnDemandDataWarehouseP>* dws)
{
  ASSERT(task);
  d_runningTasks[Thread::self()->myid()].push_back(RunningTaskInfo(task, dws));
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::popRunningTask()
{
  d_runningTasks[Thread::self()->myid()].pop_back();
}
//______________________________________________________________________
//
inline std::list<OnDemandDataWarehouse::RunningTaskInfo>*
OnDemandDataWarehouse::getRunningTasksInfo()
{
  if (d_runningTasks[Thread::self()->myid()].empty()) {
    return 0;
  } else {
    return &d_runningTasks[Thread::self()->myid()];
  }
}
//______________________________________________________________________
//
inline bool
OnDemandDataWarehouse::hasRunningTask()
{
  if (d_runningTasks[Thread::self()->myid()].empty()) {
    return false;
  } else {
    return true;
  }
}
//______________________________________________________________________
//
inline OnDemandDataWarehouse::RunningTaskInfo*
OnDemandDataWarehouse::getCurrentTaskInfo()
{
  if (d_runningTasks[Thread::self()->myid()].empty()) {
    return 0;
  } else {
    return &d_runningTasks[Thread::self()->myid()].back();
  }
}
//______________________________________________________________________
//
DataWarehouse*
OnDemandDataWarehouse::getOtherDataWarehouse(Task::WhichDW dw,
                                             RunningTaskInfo* info)
{
  int dwindex           = info->d_task->mapDataWarehouse(dw);
  DataWarehouse* result = (*info->dws)[dwindex].get_rep();
  return result;
}
//______________________________________________________________________
//
DataWarehouse*
OnDemandDataWarehouse::getOtherDataWarehouse(Task::WhichDW dw)
{
  RunningTaskInfo* info = getCurrentTaskInfo();
  int dwindex           = info->d_task->mapDataWarehouse(dw);
  DataWarehouse* result = (*info->dws)[dwindex].get_rep();
  return result;
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::checkTasksAccesses(const PatchSubset* /*patches*/,
                                          const MaterialSubset* /*matls*/)
{
#if 0
#if SCI_ASSERTION_LEVEL >= 1

  d_lock.readLock();
  
  RunningTaskInfo* currentTaskInfo = getCurrentTaskInfo();
  ASSERT(currentTaskInfo != 0);
  const Task* currentTask = currentTaskInfo->d_task;
  ASSERT(currentTask != 0);  
  
  if (isFinalized()) {
    checkAccesses(currentTaskInfo, currentTask->getRequires(), GetAccess,
                  patches, matls);
  }
  else {
    checkAccesses(currentTaskInfo, currentTask->getRequires(), GetAccess,
                  patches, matls);
    checkAccesses(currentTaskInfo, currentTask->getComputes(), PutAccess,
                  patches, matls);
    checkAccesses(currentTaskInfo, currentTask->getModifies(), ModifyAccess,
                  patches, matls);
  }

  d_lock.readUnlock();

#endif
#endif
}
//______________________________________________________________________
//
void
OnDemandDataWarehouse::checkAccesses(RunningTaskInfo* currentTaskInfo,
                                     const Task::Dependency* dep,
                                     AccessType accessType,
                                     const PatchSubset* domainPatches,
                                     const MaterialSubset* domainMatls)
{
  ASSERT(currentTaskInfo != 0);
  const Task* currentTask = currentTaskInfo->d_task;
  if (currentTask->isReductionTask()) {
    return; // no need to check reduction tasks.
  }

  VarAccessMap& currentTaskAccesses = currentTaskInfo->d_accesses;

  Handle<PatchSubset> default_patches  = scinew PatchSubset();
  Handle<MaterialSubset> default_matls = scinew MaterialSubset();
  default_patches->add(0);
  default_matls->add(-1);

  for (; dep != 0; dep = dep->next) {
#if 0
    if ((isFinalized() && dep->dw == Task::NewDW) ||
        (!isFinalized() && dep->dw == Task::OldDW))
      continue;
#endif

    const VarLabel* label = dep->var;
    IntVector lowOffset, highOffset;
    Patch::getGhostOffsets(label->typeDescription()->getType(),
                           dep->gtype,
                           dep->numGhostCells,
                           lowOffset,
                           highOffset);

    constHandle<PatchSubset> patches =
      dep->getPatchesUnderDomain(domainPatches);
    constHandle<MaterialSubset> matls =
      dep->getMaterialsUnderDomain(domainMatls);

    if (label->typeDescription() &&
        label->typeDescription()->isReductionVariable()) {
      patches = default_patches.get_rep();
    } else if (patches == 0) {
      patches = default_patches.get_rep();
    }
    if (matls == 0) {
      matls = default_matls.get_rep();
    }

    if (currentTask->getName() == "Relocate::relocateParticles") {
      continue;
    }

    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);

      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);

        VarLabelMatl<Patch> key(label, matl, patch);
        map<VarLabelMatl<Patch>, AccessInfo>::iterator find_iter;
        find_iter = currentTaskAccesses.find(key);
        if (find_iter == currentTaskAccesses.end() ||
            (*find_iter).second.accessType != accessType) {
          if ((*find_iter).second.accessType == ModifyAccess &&
              accessType == GetAccess) { // If you require with ghost cells
            continue;                    // and modify, it can get into this
          }                              // situation.

#if 1
          // THIS OLD HACK PERHAPS CAN GO AWAY
          if (lowOffset == IntVector(0, 0, 0) &&
              highOffset == IntVector(0, 0, 0)) {
            // In checkGetAccess(), this case does not record the fact
            // that the var was accessed, so don't throw exception here.
            continue;
          }
#endif
          if (find_iter == currentTaskAccesses.end()) {
            std::cout << "Error: did not find " << label->getName() << "\n";
            std::cout << "Mtl: " << m << ", Patch: " << *patch << "\n";
          } else {
            std::cout << "Error: accessType is not GetAccess for "
                      << label->getName() << "\n";
          }
          std::cout << "For Task:\n";
          currentTask->displayAll(std::cout);

          // Makes request that is never followed through.
          string has, needs;
          if (accessType == GetAccess) {
            has = "task requires";
            if (isFinalized()) {
              needs = "get from the old datawarehouse";
            } else {
              needs = "get from the new datawarehouse";
            }
          } else if (accessType == PutAccess) {
            has   = "task computes";
            needs = "datawarehouse put";
          } else {
            has   = "task modifies";
            needs = "datawarehouse modify";
          }
          // WAIT_FOR_DEBUGGER();
          SCI_THROW(DependencyException(currentTask,
                                        label,
                                        matl,
                                        patch,
                                        has,
                                        needs,
                                        __FILE__,
                                        __LINE__));
        }
        else if (((*find_iter).second.lowOffset != lowOffset ||
                  (*find_iter).second.highOffset != highOffset) &&
                 accessType != ModifyAccess /* Can == ModifyAccess when you require with
                                               ghost cells and modify */ ) {
          // Makes request for ghost cells that are never gotten.
          AccessInfo accessInfo = (*find_iter).second;
          ASSERT(accessType == GetAccess);

          // Assert that the request was greater than what was asked for
          // because the other cases (where it asked for more than the request)
          // should have been caught in checkGetAccess().
          ASSERT(Max((*find_iter).second.lowOffset, lowOffset) == lowOffset);
          ASSERT(Max((*find_iter).second.highOffset, highOffset) == highOffset);

          string has, needs;
          has = "task requires";
          ostringstream ghost_str;
          ghost_str << " requesting " << dep->numGhostCells << " layer";
          if (dep->numGhostCells > 1) {
            ghost_str << "s";
          }
          ghost_str << " of ghosts around "
                    << Ghost::getGhostTypeName(dep->gtype);
          has += ghost_str.str();

          if (isFinalized()) {
            needs = "get from the old datawarehouse";
          } else {
            needs = "get from the new datawarehouse";
          }
          needs += " that includes these ghosts";

          // WAIT_FOR_DEBUGGER();
          SCI_THROW(DependencyException(currentTask,
                                        label,
                                        matl,
                                        patch,
                                        has,
                                        needs,
                                        __FILE__,
                                        __LINE__));
        }
      }
    }
  }
}

//______________________________________________________________________
//
// For timestep abort/restart
bool
OnDemandDataWarehouse::timestepAborted()
{
  return aborted;
}
//__________________________________
//
bool
OnDemandDataWarehouse::timestepRestarted()
{
  return restart;
}
//__________________________________
//
void
OnDemandDataWarehouse::abortTimestep()
{
  // BJW - timestep aborting does not work in MPI - disabling until we get
  // fixed.
  if (d_myworld->nRanks() == 0) {
    aborted = true;
  }
}
//__________________________________
//
void
OnDemandDataWarehouse::restartTimestep()
{
  restart = true;
}
//__________________________________
//
void
OnDemandDataWarehouse::getVarLabelMatlLevelTriples(
  vector<VarLabelMatl<Level>>& vars) const
{
  d_level_DB.getVarLabelMatlTriples(vars);
}

void
OnDemandDataWarehouse::print()
{
  std::cout << d_myworld->myRank() << " VARIABLES in DW " << getID() << "\n"
            << d_myworld->myRank() << " Variable Patch Material\n"
            << "  -----------------------\n";
  d_var_DB.print(std::cout, d_myworld->myRank());
  d_level_DB.print(std::cout, d_myworld->myRank());
}

//______________________________________________________________________
//  print debugging information
void
OnDemandDataWarehouse::printDebuggingPutInfo(const VarLabel* label,
                                             int matlIndex,
                                             const Patch* patch,
                                             int line)
{
  if (dbg.active()) {
    cerrLock.lock();
    int L_indx = patch->getLevel()->getIndex();
    dbg << d_myworld->myRank() << " Putting (line: " << line << ") ";
    dbg << std::left;
    dbg.width(20);
    dbg << *label << " MI: " << matlIndex << " L-" << L_indx << " " << *patch
        << " \tinto DW: " << d_generation << "\n";
    cerrLock.unlock();
  }
}

} // namespace Uintah