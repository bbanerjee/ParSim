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

#include <Core/Util/StringUtil.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/FancyAssert.h>

#include <set>

namespace Uintah {

void
Task::initialize()
{
  m_comp_head = nullptr;
  m_comp_tail = nullptr;
  m_req_head = nullptr;
  m_req_tail = nullptr;
  m_mod_head = nullptr;
  m_mod_tail = nullptr;
  m_patch_set = nullptr;
  m_matl_set = nullptr;

  m_uses_mpi = false;
  m_uses_threads = false;
  m_uses_device = false;
  m_subpatch_capable = false;
  m_has_subscheduler = false;

  for (int i = 0; i < TotalDWs; i++) {
    m_dwmap[i] = Task::InvalidDW;
  }

  m_sorted_order = -1;
  m_phase = -1;
  m_comm = -1;

  // The 0th level has a max ghost cell of zero.  Other levels are left
  // uninitialized.
  m_max_ghost_cells[0] = 0;
  m_max_level_offset = 0;
}

Task::~Task()
{
  delete m_action;

  Dependency* dep = m_req_head;
  while (dep) {
    Dependency* next = dep->next;
    delete dep;
    dep = next;
  }

  dep = m_comp_head;
  while (dep) {
    Dependency* next = dep->next;
    delete dep;
    dep = next;
  }

  dep = m_mod_head;
  while (dep) {
    Dependency* next = dep->next;
    delete dep;
    dep = next;
  }

  if (m_matl_set && m_matl_set->removeReference()) {
    delete m_matl_set;
  }

  if (m_patch_set && m_patch_set->removeReference()) {
    delete m_patch_set;
  }

  // easier to periodically delete this than to force a call to a cleanup
  // function, and probably not very expensive.
  if (s_global_material_subset && s_global_material_subset->removeReference()) {
    delete s_global_material_subset;
  }

  s_global_material_subset = nullptr;
}

void
Task::setSets(const PatchSet* ps, const MaterialSet* ms)
{
  // NOTE: the outer [patch/matl]Set checks are related to temporal scheduling,
  // e.g. more then 1 regular task graph
  //
  // This is called from TaskGraph::addTask() in which a single task may be
  // added to more than 1 Normal task graphs. In this case, first time here,
  // m_path/matl_set will be nullptr and subsequent visits will be the same
  // pointer as ps and ms respectively. Without these checks, the refCount gets
  // artificially inflated and ComputeSubsets (Patch/Matl)) are not deleted,
  // resulting in a mem leak. APH, 06/08/17
  if (m_patch_set == nullptr) {
    m_patch_set = ps;
    if (m_patch_set) {
      m_patch_set->addReference();
    }
  }
  if (m_matl_set == nullptr) {
    m_matl_set = ms;
    if (m_matl_set) {
      m_matl_set->addReference();
    }
  }
}

const MaterialSubset*
Task::getGlobalMatlSubset()
{
  if (s_global_material_subset == nullptr) {
    s_global_material_subset = scinew MaterialSubset();
    s_global_material_subset->add(-1);
    s_global_material_subset->addReference();
  }
  return s_global_material_subset;
}

void
Task::usesMPI(bool state)
{
  m_uses_mpi = state;
}

void
Task::hasSubScheduler(bool state)
{
  m_has_subscheduler = state;
}

void
Task::usesThreads(bool state)
{
  m_uses_threads = state;
}

void
Task::usesDevice(bool state, int maxStreamsPerTask)
{
  m_uses_device = state;
  m_max_streams_per_task = maxStreamsPerTask;
}

/* clang-format off */
void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const PatchSubset* patches,
               PatchDomainSpec patches_dom,
               int level_offset,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  if (matls == 0 && var->typeDescription()->isReductionVariable()) {
    // default material for a reduction variable is the global
    // material (-1)
    matls = getGlobalMatlSubset();
    matls_dom = OutOfDomain;
  } else if (matls != 0 && matls->size() == 0) {
    return; // no materials, no dependency
  }

  Dependency* dep = scinew Dependency(Requires, this, dw, var, whichTG, patches,
                                      matls, patches_dom, matls_dom, gtype,
                                      numGhostCells, level_offset);

  if (level_offset > m_max_level_offset) {
    m_max_level_offset = level_offset;
  }

  dep->next = 0;
  if (m_req_tail) {
    m_req_tail->next = dep;
  } else {
    m_req_head = dep;
  }
  m_req_tail = dep;

  if (dw == OldDW) {
    m_requires_old_dw.insert(std::make_pair(var, dep));
  } else {
    m_requires.insert(std::make_pair(var, dep));
  }
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const PatchSubset* patches,
               PatchDomainSpec patches_dom,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  int offset = 0;
  if (patches_dom == CoarseLevel || patches_dom == FineLevel) {
    offset = 1;
  }

  requires(dw, var, patches, patches_dom, offset, matls, matls_dom, gtype,
           numGhostCells, whichTG);
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  requires(dw, var, patches, ThisLevel, matls, NormalDomain, gtype, numGhostCells,
           whichTG);
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  requires(dw, var, nullptr, ThisLevel, nullptr, NormalDomain, gtype,
           numGhostCells, whichTG);
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const MaterialSubset* matls,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG)
{
  requires(
    dw, var, nullptr, ThisLevel, matls, NormalDomain, gtype, numGhostCells, whichTG);
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  requires(dw, var, nullptr, ThisLevel, matls, matls_dom, gtype, numGhostCells,
           whichTG);
}

void 
Task::requires(WhichDW dw,
               const VarLabel* var,
               const PatchSubset* patches,
               Ghost::GhostType gtype,
               int numGhostCells,
               SearchTG whichTG) 
{
  requires(dw, var, patches, ThisLevel, nullptr, NormalDomain, gtype,
           numGhostCells, whichTG);
}

void
Task::requires(WhichDW dw,
               const VarLabel* var,
               const PatchSubset* patches,
               const MaterialSubset* matls) 
{
  TypeDescription::Type vartype = var->typeDescription()->getType();
  if (vartype == TypeDescription::Type::SoleVariable) {
    requires(dw, var, (const Level*)0, matls);
  } else if (vartype == TypeDescription::Type::PerPatch) {
    requires(dw, var, patches, ThisLevel, matls, NormalDomain, Ghost::None, 0);
  } else {
    SCI_THROW(InternalError(
              "Requires should specify ghost type or level for this variable",
              __FILE__, __LINE__));
  }
}

void Task::requires(WhichDW dw,
                    const VarLabel* var,
                    const MaterialSubset* matls,
                    SearchTG whichTG) 
{
  TypeDescription::Type vartype = var->typeDescription()->getType();

  if (vartype == TypeDescription::Type::ReductionVariable) {
    requires(dw, var, (const Level*) nullptr, matls, NormalDomain, whichTG);
  } else if (vartype == TypeDescription::Type::SoleVariable) {
    requires(dw, var, (const Level*) nullptr, matls);
  } else if (vartype == TypeDescription::Type::PerPatch) {
    requires(dw, var, nullptr, ThisLevel, matls, NormalDomain,
            Ghost::None, 0, whichTG);
  } else {
    SCI_THROW(InternalError("Requires should specify ghost type for this variable", 
                            __FILE__, __LINE__));
  }
}


//__________________________________
void Task::requires(WhichDW dw,
                    const VarLabel* var,
                    const Level* level,
                    const MaterialSubset* matls,
                    MaterialDomainSpec matls_dom,
                    SearchTG whichTG) 
{
  TypeDescription::Type vartype = var->typeDescription()->getType();
  if ((vartype == TypeDescription::Type::ReductionVariable ||
       vartype == TypeDescription::Type::SoleVariable)) {

    if (matls == nullptr) {
      // default material for a reduction variable is the global
      // material (-1)
      matls = getGlobalMatlSubset();
      matls_dom = OutOfDomain;
    } else if (matls->size() == 0) {
      return; // no materials, no dependency
    }

    Dependency* dep = scinew Dependency(
        Requires, this, dw, var, whichTG, level, matls, matls_dom);
    dep->next = nullptr;

    if (m_req_tail) {
      m_req_tail->next = dep;
    } else {
      m_req_head = dep;
    }
    m_req_tail = dep;

    if (dw == OldDW) {
      m_requires_old_dw.insert(std::make_pair(var, dep));
    } else {
      m_requires.insert(std::make_pair(var, dep));
    }
  } else {
    SCI_THROW(InternalError("Requires should specify ghost type for this variable", 
                            __FILE__, __LINE__));
  }
}

/* clang-format on */
//__________________________________
void
Task::computes(const VarLabel* var,
               const PatchSubset* patches,
               PatchDomainSpec patches_dom,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom)
{
  if (var->typeDescription()->isReductionVariable()) {
    if (matls == nullptr) {
      // default material for a reduction variable is the global material
      // (-1)
      matls = getGlobalMatlSubset();
      matls_dom = OutOfDomain;
    }
    ASSERT(patches == nullptr);
  }

  Dependency* dep = scinew Dependency(Computes,
                                      this,
                                      NewDW,
                                      var,
                                      SearchTG::NewTG,
                                      patches,
                                      matls,
                                      patches_dom,
                                      matls_dom);
  dep->next = nullptr;

  if (m_comp_tail) {
    m_comp_tail->next = dep;
  } else {
    m_comp_head = dep;
  }
  m_comp_tail = dep;

  m_computes.insert(std::make_pair(var, dep));
}

void
Task::computes(const VarLabel* var,
               const PatchSubset* patches,
               const MaterialSubset* matls)
{
  TypeDescription::Type vartype = var->typeDescription()->getType();
  if (vartype == TypeDescription::Type::ReductionVariable ||
      vartype == TypeDescription::Type::SoleVariable) {
    computes(var, (const Level*)nullptr, matls);
  } else {
    computes(var, patches, ThisLevel, matls, NormalDomain);
  }
}

void
Task::computes(const VarLabel* var, const MaterialSubset* matls)
{
  computes(var, nullptr, ThisLevel, matls, NormalDomain);
}

void
Task::computes(const VarLabel* var,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom)
{
  computes(var, nullptr, ThisLevel, matls, matls_dom);
}

void
Task::computes(const VarLabel* var,
               const PatchSubset* patches,
               PatchDomainSpec patches_dom)
{
  computes(var, patches, patches_dom, nullptr, NormalDomain);
}

void
Task::computes(const VarLabel* var,
               const Level* level,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom)
{
  TypeDescription::Type vartype = var->typeDescription()->getType();
  if (vartype == TypeDescription::Type::ReductionVariable ||
      vartype == TypeDescription::Type::SoleVariable) {

    if (matls == 0) {
      // default material for a reduction variable is the global material (-1)
      matls = getGlobalMatlSubset();
      matls_dom = OutOfDomain;
    } else if (matls->size() == 0) {
      throw InternalError(
        "Computes of an empty material set!", __FILE__, __LINE__);
    }

    Dependency* dep = scinew Dependency(
      Computes, this, NewDW, var, SearchTG::NewTG, level, matls, matls_dom);
    dep->next = nullptr;

    if (m_comp_tail) {
      m_comp_tail->next = dep;
    } else {
      m_comp_head = dep;
    }
    m_comp_tail = dep;

    m_computes.insert(std::make_pair(var, dep));

  } else {
    SCI_THROW(
      InternalError("Computes should only be used for reduction variable",
                    __FILE__,
                    __LINE__));
  }
}

void
Task::computesWithScratchGhost(const VarLabel* var,
                               const MaterialSubset* matls,
                               MaterialDomainSpec matls_dom,
                               Ghost::GhostType gtype,
                               int numGhostCells,
                               SearchTG whichTG)
{
  if (var->typeDescription()->isReductionVariable()) {
    SCI_THROW(InternalError(
      "ComputeswithScratchGhost should not be used for reduction variable",
      __FILE__,
      __LINE__));
  }

  Dependency* dep = scinew Dependency(Computes,
                                      this,
                                      NewDW,
                                      var,
                                      whichTG,
                                      nullptr,
                                      matls,
                                      ThisLevel,
                                      matls_dom,
                                      gtype,
                                      numGhostCells);
  dep->next = nullptr;

  if (m_comp_tail) {
    m_comp_tail->next = dep;
  } else {
    m_comp_head = dep;
  }

  m_comp_tail = dep;

  m_computes.insert(std::make_pair(var, dep));
}

void
Task::modifiesWithScratchGhost(const VarLabel* var,
                               const PatchSubset* patches,
                               PatchDomainSpec patches_dom,
                               const MaterialSubset* matls,
                               MaterialDomainSpec matls_dom,
                               Ghost::GhostType gtype,
                               int numGhostCells,
                               SearchTG whichTG)
{
  this->requires(
    NewDW, var, patches, patches_dom, matls, matls_dom, gtype, numGhostCells);
  this->modifies(var, patches, patches_dom, matls, matls_dom);
}

//__________________________________
void
Task::modifies(const VarLabel* var,
               const PatchSubset* patches,
               PatchDomainSpec patches_dom,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom,
               SearchTG whichTG /* = NewTG */)
{
  if (matls == nullptr && var->typeDescription()->isReductionVariable()) {
    // default material for a reduction variable is the global material (-1)
    matls = getGlobalMatlSubset();
    matls_dom = OutOfDomain;
    ASSERT(patches == nullptr);
  }

  Dependency* dep = scinew Dependency(Modifies,
                                      this,
                                      NewDW,
                                      var,
                                      whichTG,
                                      patches,
                                      matls,
                                      patches_dom,
                                      matls_dom);
  dep->next = 0;
  if (m_mod_tail) {
    m_mod_tail->next = dep;
  } else {
    m_mod_head = dep;
  }
  m_mod_tail = dep;

  m_requires.insert(std::make_pair(var, dep));
  m_computes.insert(std::make_pair(var, dep));
  m_modifies.insert(std::make_pair(var, dep));
}

void
Task::modifies(const VarLabel* var,
               const Level* level,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_domain,
               SearchTG whichTG)
{
  const TypeDescription* vartype = var->typeDescription();

  if (matls == nullptr && vartype->isReductionVariable()) {
    // default material for a reduction variable is the global material (-1)
    matls = getGlobalMatlSubset();
    matls_domain = OutOfDomain;
  }

  if (!vartype->isReductionVariable()) {
    SCI_THROW(InternalError(
      "modifies with level should only be used for reduction variable",
      __FILE__,
      __LINE__));
  }

  Dependency* dep = scinew Dependency(
    Modifies, this, NewDW, var, whichTG, level, matls, matls_domain);
  dep->next = nullptr;
  if (m_mod_tail) {
    m_mod_tail->next = dep;
  } else {
    m_mod_head = dep;
  }
  m_mod_tail = dep;

  m_requires.insert(std::make_pair(var, dep));
  m_computes.insert(std::make_pair(var, dep));
  m_modifies.insert(std::make_pair(var, dep));
}

void
Task::modifies(const VarLabel* var,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               SearchTG whichTG)
{
  modifies(var, patches, ThisLevel, matls, NormalDomain, whichTG);
}

void
Task::modifies(const VarLabel* var, SearchTG whichTG)
{
  modifies(var, nullptr, ThisLevel, nullptr, NormalDomain, whichTG);
}

void
Task::modifies(const VarLabel* var,
               const MaterialSubset* matls,
               SearchTG whichTG)
{
  modifies(var, nullptr, ThisLevel, matls, NormalDomain, whichTG);
}

void
Task::modifies(const VarLabel* var,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom,
               SearchTG whichTG)
{
  modifies(var, nullptr, ThisLevel, matls, matls_dom, whichTG);
}

bool
Task::hasComputes(const VarLabel* var, int matlIndex, const Patch* patch) const
{
  return isInDepMap(m_computes, var, matlIndex, patch);
}

bool
Task::hasRequires(const VarLabel* var,
                  int matlIndex,
                  const Patch* patch,
                  IntVector lowOffset,
                  IntVector highOffset,
                  WhichDW dw) const
{
  DepMap depMap = m_requires;

  if (dw == OldDW) {
    depMap = m_requires_old_dw;
  }

  Dependency* dep = isInDepMap(depMap, var, matlIndex, patch);

  if (dep) {
    // make sure we are within the allowed ghost cell limit
    IntVector allowableLowOffset, allowableHighOffset;

    Patch::getGhostOffsets(var->typeDescription()->getType(),
                           dep->gtype,
                           dep->num_ghost_cells,
                           allowableLowOffset,
                           allowableHighOffset);

    return ((Max(allowableLowOffset, lowOffset) == allowableLowOffset) &&
            (Max(allowableHighOffset, highOffset) == allowableHighOffset));
  }
  return false;
}

bool
Task::hasDistalRequires() const
{
  for (auto dep = m_req_head; dep != nullptr; dep = dep->next) {
    if (dep->num_ghost_cells >= MAX_HALO_DEPTH) {
      return true;
    }
  }
  return false;
}

bool
Task::hasModifies(const VarLabel* var, int matlIndex, const Patch* patch) const
{
  return isInDepMap(m_modifies, var, matlIndex, patch);
}

Task::Dependency*
Task::isInDepMap(const DepMap& depMap,
                 const VarLabel* var,
                 int matlIndex,
                 const Patch* patch) const
{
  DepMap::const_iterator found_iter = depMap.find(var);

  // loop over dependency map and search for the right dependency
  while (found_iter != depMap.end() && (*found_iter).first->equals(var)) {
    Dependency* dep = (*found_iter).second;
    const PatchSubset* patches = dep->patches;
    const MaterialSubset* matls = dep->matls;

    bool hasPatches = false, hasMatls = false;

    if (patches == nullptr) {
      // if patches==0 then the requirement for patches
      // is satisfied
      hasPatches = true;
    } else {
      if (dep->patches_dom == Task::CoarseLevel) {
        // check that the level of the patches matches
        // the coarse level
        hasPatches =
          getLevel(getPatchSet())->getRelativeLevel(-dep->level_offset) ==
          getLevel(patches);
      } else if (dep->patches_dom == Task::FineLevel) {
        // check that the level of the patches
        // matches the fine level
        hasPatches =
          getLevel(getPatchSet())->getRelativeLevel(dep->level_offset) ==
          getLevel(patches);
      } else {
        // check that the patches subset contain the requested patch
        hasPatches = patches->contains(patch);
      }
    }

    if (matls == nullptr) {
      // if matls==0 then the requierment for matls is satisfied
      hasMatls = true;
    } else {
      // check thta the malts subset contains the matlIndex
      hasMatls = matls->contains(matlIndex);
    }

    if (hasMatls && hasPatches) {
      // if this dependency contains both the
      // matls and patches return the dependency
      return dep;
    }
    found_iter++;
  }
  return nullptr;
}

//----------------------------------------------------------
Task::Dependency::Dependency(DepType deptype,
                             Task* task,
                             WhichDW whichdw,
                             const VarLabel* var,
                             SearchTG whichTG,
                             const PatchSubset* patches,
                             const MaterialSubset* matls,
                             PatchDomainSpec patches_dom,
                             MaterialDomainSpec matls_dom,
                             Ghost::GhostType gtype,
                             int numGhostCells,
                             int level_offset)

  : dep_type(deptype)
  , task(task)
  , var(var)
  , look_in_old_tg(whichTG == SearchTG::OldTG)
  , patches(patches)
  , matls(matls)
  , patches_dom(patches_dom)
  , matls_dom(matls_dom)
  , gtype(gtype)
  , whichdw(whichdw)
  , num_ghost_cells(numGhostCells)
  , level_offset(level_offset)
{
  if (var) {
    var->addReference();
  }

  if (patches) {
    patches->addReference();
  }

  if (matls) {
    matls->addReference();
  }
}

Task::Dependency::Dependency(DepType deptype,
                             Task* task,
                             WhichDW whichdw,
                             const VarLabel* var,
                             SearchTG whichTG,
                             const Level* reductionLevel,
                             const MaterialSubset* matls,
                             MaterialDomainSpec matls_dom)

  : dep_type(deptype)
  , task(task)
  , var(var)
  , look_in_old_tg(whichTG == SearchTG::OldTG)
  , matls(matls)
  , reduction_level(reductionLevel)
  , matls_dom(matls_dom)
  , gtype(Ghost::None)
  , whichdw(whichdw)
{
  if (var) {
    var->addReference();
  }

  if (matls) {
    matls->addReference();
  }
}

Task::Dependency::~Dependency()
{
  VarLabel::destroy(var); // just remove the ref

  if (patches && patches->removeReference()) {
    delete patches;
  }

  if (matls && matls->removeReference()) {
    delete matls;
  }
}

constHandle<PatchSubset>
Task::Dependency::getPatchesUnderDomain(const PatchSubset* domainPatches) const
{
  switch (patches_dom) {
    case Task::ThisLevel:
    case Task::OtherGridDomain: // use the same patches, we'll figure out
                                // where it corresponds on the other grid
      return PatchSubset::intersection(patches, domainPatches);
    case Task::CoarseLevel:
    case Task::FineLevel:
      return getOtherLevelPatchSubset(
        patches_dom, level_offset, patches, domainPatches, num_ghost_cells);
    default:
      SCI_THROW(
        InternalError(string("Unknown patch domain ") + " type " +
                        Uintah::to_string(static_cast<int>(patches_dom)),
                      __FILE__,
                      __LINE__));
  }
}

constHandle<MaterialSubset>
Task::Dependency::getMaterialsUnderDomain(
  const MaterialSubset* domainMaterials) const
{
  switch (matls_dom) {
    case Task::NormalDomain:
      return MaterialSubset::intersection(matls, domainMaterials);
    case Task::OutOfDomain:
      return matls;
    default:
      SCI_THROW(InternalError(string("Unknown matl domain ") + " type " +
                                Uintah::to_string(static_cast<int>(matls_dom)),
                              __FILE__,
                              __LINE__));
  }
}

constHandle<PatchSubset>
Task::Dependency::getOtherLevelPatchSubset(Task::PatchDomainSpec dom,
                                           int level_offset,
                                           const PatchSubset* subset,
                                           const PatchSubset* domainSubset,
                                           int ngc)
{
  constHandle<PatchSubset> myLevelSubset =
    PatchSubset::intersection(subset, domainSubset);

  int levelOffset = 0;
  switch (dom) {
    case Task::CoarseLevel:
      levelOffset = -level_offset;
      break;
    case Task::FineLevel:
      levelOffset = level_offset;
      break;
    default:
      SCI_THROW(InternalError("Unhandled DomainSpec in "
                              "Task::Dependency::getOtherLevelComputeSubset",
                              __FILE__,
                              __LINE__));
  }

  std::set<const Patch*, Patch::Compare> patches;
  for (int p = 0; p < myLevelSubset->size(); p++) {
    const Patch* patch = myLevelSubset->get(p);
    Patch::selectType somePatches;
    patch->getOtherLevelPatches(levelOffset, somePatches, ngc);
    patches.insert(somePatches.begin(), somePatches.end());
  }

  return constHandle<PatchSubset>(
    scinew PatchSubset(patches.begin(), patches.end()));
}

void
Task::doit(DetailedTask* dtask,
           CallBackEvent event,
           const ProcessorGroup* pg,
           const PatchSubset* patches,
           const MaterialSubset* matls,
           std::vector<DataWarehouseP>& dws,
           void* oldTaskGpuDW,
           void* newTaskGpuDW,
           void* stream,
           int deviceID)
{
  DataWarehouse* fromDW = mapDataWarehouse(Task::OldDW, dws);
  DataWarehouse* toDW = mapDataWarehouse(Task::NewDW, dws);
  if (m_action) {
    m_action->doit(dtask,
                   event,
                   pg,
                   patches,
                   matls,
                   fromDW,
                   toDW,
                   oldTaskGpuDW,
                   newTaskGpuDW,
                   stream,
                   deviceID);
  }
}

void
Task::display(std::ostream& out) const
{
  out << Parallel::getMPIRank() << " " << getName() << " (" << m_tasktype
      << "): [";
  if (m_uses_device) {
    out << ": GPU task,";
  }

  out << " (" << m_tasktype << ")";

  if (m_tasktype == Task::Normal && m_patch_set != nullptr) {
    out << ", Level " << getLevel(m_patch_set)->getIndex();
  }

  if (m_matl_set == nullptr) {
    out << ", No-Matl-Set";
  } else {
    out << ", " << *m_matl_set;
  }
  out << ", DWs: ";
  for (int i = 0; i < TotalDWs; i++) {
    if (i != 0)
      out << ", ";
    out << m_dwmap[i];
  }
  if (m_patch_set == nullptr) {
    out << ", No-Patch-Set";
  } else {
    out << ", " << *m_patch_set;
  }
}

std::ostream&
operator<<(std::ostream& out, const Uintah::Task::Dependency& dep)
{
  out << "[";
  out << std::left;
  out.width(20);
  out << *(dep.var) << ", ";

  // reduction variable
  if (dep.var->typeDescription()->isReductionVariable()) {
    if (dep.reduction_level) {
      out << " reduction Level: " << dep.reduction_level->getIndex();
    } else {
      out << " Global level";
    }
  } else {
    // all other variables:
    if (dep.patches) {
      out << " Level: " << getLevel(dep.patches)->getIndex();
      out << " Patches: ";
      for (int i = 0; i < dep.patches->size(); i++) {
        if (i > 0) {
          out << ",";
        }
        out << dep.patches->get(i)->getID();
      }
    } else if (dep.reduction_level) {
      out << " reduction Level: " << dep.reduction_level->getIndex();
    } else if (dep.patches_dom) {
      switch (dep.patches_dom) {
        case Task::CoarseLevel:
          out << "coarseLevel";
          break;
        case Task::FineLevel:
          out << "fineLevel";
          break;
        case Task::OtherGridDomain:
          out << "OtherGridDomain";
          break;
        case Task::ThisLevel:
          out << "ThisLevel";
          break;
        default:
          break;
      }
    } else {
      out << "all Patches";
    }
  }

  out << ", MI: ";
  if (dep.matls) {
    for (int i = 0; i < dep.matls->size(); i++) {
      if (i > 0)
        out << ",";
      out << dep.matls->get(i);
    }
  } else {
    out << "none";
  }
  out << ", ";
  switch (dep.whichdw) {
    case Task::OldDW:
      out << "OldDW";
      break;
    case Task::NewDW:
      out << "NewDW";
      break;
    case Task::CoarseOldDW:
      out << "CoarseOldDW";
      break;
    case Task::CoarseNewDW:
      out << "CoarseNewDW";
      break;
    case Task::ParentOldDW:
      out << "ParentOldDW";
      break;
    case Task::ParentNewDW:
      out << "ParentNewDW";
      break;
    default:
      out << "Unknown DW!";
      break;
  }
  out << " (mapped to dw index " << dep.task->mapDataWarehouse(dep.whichdw)
      << ")";
  out << ", matl domain:";
  switch (dep.matls_dom) {
    case Task::NormalDomain:
      out << "normal, ";
      break;
    case Task::OutOfDomain:
      out << "OutOfDomain, ";
      break;
    default:
      out << "Unknown, ";
      break;
  }

  switch (dep.gtype) {
    case Ghost::None:
      out << "Ghost::None";
      break;
    case Ghost::AroundNodes:
      out << "Ghost::AroundNodes";
      break;
    case Ghost::AroundCells:
      out << "Ghost::AroundCells";
      break;
    case Ghost::AroundFacesX:
      out << "Ghost::AroundFacesX";
      break;
    case Ghost::AroundFacesY:
      out << "Ghost::AroundFacesY";
      break;
    case Ghost::AroundFacesZ:
      out << "Ghost::AroundFacesZ";
      break;
    case Ghost::AroundFaces:
      out << "Ghost::AroundFaces";
      break;
    default:
      out << "Unknown ghost type";
      break;
  }
  if (dep.gtype != Ghost::None)
    out << ":" << dep.num_ghost_cells;

  out << "]";
  return out;
}

ostream&
operator<<(std::ostream& out, const Task& task)
{
  task.display(out);
  return out;
}

ostream&
operator<<(std::ostream& out, const Task::TaskType& tt)
{
  switch (tt) {
    case Task::Normal:
      out << "Normal";
      break;
    case Task::OutputGlobalVars:
      out << "OutputGlobalVars";
      break;
    case Task::Reduction:
      out << "Reduction";
      break;
    case Task::InitialSend:
      out << "InitialSend";
      break;
    case Task::Output:
      out << "Output";
      break;
    case Task::OncePerProc:
      out << "OncePerProc";
      break;
    case Task::Spatial:
      out << "Spatial";
      break;
    case Task::Hypre:
      out << "Hypre";
      break;
  }
  return out;
}

void
Task::displayAll_DOUT(Uintah::Dout& dbg) const
{
  if (dbg.active()) {
    std::ostringstream message;
    displayAll(message);
    DOUT(dbg, message.str());
  }
}

void
Task::displayAll(std::ostream& out) const
{
  display(out);
  out << '\n';
  for (Task::Dependency* req = m_req_head; req != nullptr; req = req->next) {
    out << Parallel::getMPIRank() << "  requires: " << *req << '\n';
  }
  for (Task::Dependency* comp = m_comp_head; comp != nullptr; comp = comp->next) {
    out << Parallel::getMPIRank() << "  computes: " << *comp << '\n';
  }
  for (Task::Dependency* mod = m_mod_head; mod != nullptr; mod = mod->next) {
    out << Parallel::getMPIRank() << "  modifies: " << *mod << '\n';
  }
}

void
Task::setMapping(int dwmap[TotalDWs])
{
  for (int i = 0; i < TotalDWs; i++) {
    this->m_dwmap[i] = dwmap[i];
  }
}

int
Task::mapDataWarehouse(WhichDW dw) const
{
  ASSERTRANGE(dw, 0, Task::TotalDWs);
  return m_dwmap[dw];
}

//__________________________________
DataWarehouse*
Task::mapDataWarehouse(WhichDW dw, std::vector<DataWarehouseP>& dws) const
{
  ASSERTRANGE(dw, 0, Task::TotalDWs);
  if (m_dwmap[dw] == Task::NoDW) {
    return 0;
  } else {
    ASSERTRANGE(m_dwmap[dw], 0, (int)dws.size());
    return dws[m_dwmap[dw]].get_rep();
  }
}

} // namespace Uintah