/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef VAANGO_CCA_COMPONENTS_SCHEDULERS_DWDatabase_H
#define VAANGO_CCA_COMPONENTS_SCHEDULERS_DWDatabase_H

#include <CCA/Components/Schedulers/MemoryLog.h>

#include <Core/Grid/UnknownVariable.h>
#include <Core/Grid/Variables/ReductionVariableBase.h>
#include <Core/Grid/Variables/ScrubItem.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarLabelMatl.h>

#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/Parallel.h>

#include <Core/Containers/FastHashTable.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Malloc/Allocator.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/FancyAssert.h>

#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

Uintah::MasterLock g_keyDB_lock{};
Uintah::MasterLock g_mvars_lock{};

}

namespace Uintah {

template<class DomainType>
class KeyDatabase
{

  template<class T>
  friend class DWDatabase;

public:
  KeyDatabase(){};

  ~KeyDatabase() noexcept(false) {};

  void
  clear();

  void
  insert(const VarLabel* label, int matlIndex, const DomainType* dom);

  int
  lookup(const VarLabel* label, int matlIndex, const DomainType* dom);

  void
  merge(const KeyDatabase<DomainType>& newDB);

  void
  print(const int rank) const;

private:
  using keyDBtype = std::unordered_map<VarLabelMatl<DomainType>, int>;
  keyDBtype m_keys;

  int m_key_count{ 0 };
};

template<class DomainType>
class DWDatabase
{
public:
  DWDatabase() = default;
  ~DWDatabase();

  // eliminate copy, assignment and move
  DWDatabase(const DWDatabase&) = delete;
  DWDatabase&
  operator=(const DWDatabase&) = delete;
  DWDatabase(DWDatabase&&)     = delete;
  DWDatabase&
  operator=(DWDatabase&&) = delete;

  void
  clear();

  void
  doReserve(KeyDatabase<DomainType>* keydb);

  bool
  exists(const VarLabel* label, int matlIndex, const DomainType* dom) const;

  void
  put(const VarLabel* label,
      int matlindex,
      const DomainType* dom,
      Variable* var,
      bool init,
      bool replace);

  void
  putReduce(const VarLabel* label,
            int matlindex,
            const DomainType* dom,
            ReductionVariableBase* var,
            bool init);

  void
  putForeign(const VarLabel* label,
             int matlindex,
             const DomainType* dom,
             Variable* var,
             bool init);

  void
  get(const VarLabel* label,
      int matlindex,
      const DomainType* dom,
      Variable& var) const;

  void
  getlist(const VarLabel* label,
          int matlIndex,
          const DomainType* dom,
          std::vector<Variable*>& varlist) const;

  inline Variable*
  get(const VarLabel* label, int matlindex, const DomainType* dom) const;

  void
  print(int rank) const;

  void
  cleanForeign();

  // Scrub counter manipulator functions -- when the scrub count goes to
  // zero, the data is scrubbed.  Return remaining count

  // How Scrubbing works:  If a variable is in the OldDW, at the beginning of
  // the timestep initializeScrubs will be called for each of those variables.
  // For each variable computed or copied via MPI, setScrubCount will be called
  // on it, based on the scrubCountTable in DetailedTasks.  Then, when the
  // variable is used, decrementScrubCount is called on it and if the count
  // reaches zero, it is scrubbed.
  int
  decrementScrubCount(const VarLabel* label,
                      int matlindex,
                      const DomainType* dom);

  void
  setScrubCount(const VarLabel* label,
                int matlindex,
                const DomainType* dom,
                int count);

  void
  scrub(const VarLabel* label, int matlindex, const DomainType* dom);

  // add means increment the scrub count instead of setting it.  This is for
  // when a DW can act as a CoarseOldDW as well as an OldDW
  void
  initializeScrubs(int dwid,
                   const FastHashTable<ScrubItem>* scrubcounts,
                   bool add);

  void
  logMemoryUse(std::ostream& out,
               unsigned long& total,
               const std::string& tag,
               int dwid);

  void
  getVarLabelMatlTriples(std::vector<VarLabelMatl<DomainType>>& vars) const;

private:
  struct DataItem
  {
    DataItem() = default;

    ~DataItem() noexcept(false)
    {
      if (next) {
        delete next;
      }
      ASSERT(var);
      delete var;
    }

    Variable* var;
    struct DataItem* next;
  };

  DataItem*
  getDataItem(const VarLabel* label,
              int matlindex,
              const DomainType* dom) const;

private:
  KeyDatabase<DomainType>* m_keyDB{ nullptr };

  using varDBtype   = std::vector<DataItem*>;
  using scrubDBtype = std::vector<int>;

  varDBtype m_vars{};
  scrubDBtype m_scrubs{};
};

template<class DomainType>
int
KeyDatabase<DomainType>::lookup(const VarLabel* label,
                                int matlIndex,
                                const DomainType* dom)
{
  VarLabelMatl<DomainType> v(label, matlIndex, getRealDomain(dom));
  typename keyDBtype::const_iterator const_iter = m_keys.find(v);
  if (const_iter == m_keys.end()) {
    return -1;
  } else {
    return const_iter->second;
  }
}

template<class DomainType>
void
KeyDatabase<DomainType>::merge(const KeyDatabase<DomainType>& newDB)
{
  for (auto const_keyiter = newDB.m_keys.cbegin();
       const_keyiter != newDB.m_keys.cend();
       ++const_keyiter) {
    auto const_db_iter = m_keys.find(const_keyiter->first);
    if (const_db_iter == m_keys.end()) {
      m_keys.insert(
        std::pair<VarLabelMatl<DomainType>, int>(const_keyiter->first,
                                                 m_key_count++));
    }
  }
}

template<class DomainType>
void
KeyDatabase<DomainType>::insert(const VarLabel* label,
                                int matlIndex,
                                const DomainType* dom)
{
  VarLabelMatl<DomainType> v(label, matlIndex, getRealDomain(dom));
  typename keyDBtype::const_iterator const_iter = m_keys.find(v);
  if (const_iter == m_keys.end()) {
    m_keys.insert(std::pair<VarLabelMatl<DomainType>, int>(v, m_key_count++));
  }
}

template<class DomainType>
void
KeyDatabase<DomainType>::clear()
{
  m_keys.clear();
  m_key_count = 0;
}

template<class DomainType>
void
KeyDatabase<DomainType>::print(const int rank) const
{
  DOUT(true,
       "Rank-" << rank << " __________________________________KeyDatabase ")
  for (auto keyiter = m_keys.begin(); keyiter != m_keys.end(); ++keyiter) {
    const VarLabelMatl<DomainType>& vlm = keyiter->first;
    const DomainType* dom               = vlm.m_domain;
    if (dom) {
      DOUT(true,
           "Rank-" << rank << " Name: " << vlm.m_label->getName()
                   << "  domain: " << *dom << "  matl:" << vlm.m_matl_index);
    } else {
      DOUT(true,
           "Rank-" << rank << " Name: " << vlm.m_label->getName()
                   << "  domain: nullptr  matl: " << vlm.m_matl_index);
    }
  }
  DOUT(true, "Rank-" << rank << " __________________________________")
}

template<class DomainType>
DWDatabase<DomainType>::~DWDatabase()
{
  clear();
}

template<class DomainType>
void
DWDatabase<DomainType>::clear()
{
  for (auto iter = m_vars.begin(); iter != m_vars.end(); ++iter) {
    if (*iter) {
      delete *iter;
    }
    *iter = nullptr;
  }
  m_vars.clear();
}

template<class DomainType>
void
DWDatabase<DomainType>::cleanForeign()
{
  for (auto iter = m_vars.begin(); iter != m_vars.end(); ++iter) {
    if (*iter && (*iter)->var->isForeign()) {
      delete (*iter);
      (*iter) = nullptr;
    }
  }
}

template<class DomainType>
int
DWDatabase<DomainType>::decrementScrubCount(const VarLabel* label,
                                            int matlIndex,
                                            const DomainType* dom)
{
  // Dav's conjectures on how this works:
  //   setScrubCount is called the first time with "count" set to some X.
  //   This X represents the number of tasks that will use the var.  Later,
  //   after a task has used the var, it will call decrementScrubCount
  //   If scrubCount then is equal to 0, the var is scrubbed.

  ASSERT(matlIndex >= -1);

  std::lock_guard<MasterLock> decrement_scrub_count_lock(g_keyDB_lock);

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    return 0;
  }
  if (!m_vars[idx]) {
    return 0;
  }

  int rt = __sync_sub_and_fetch(&(m_scrubs[idx]), 1);
  if (rt == 0) {
    delete m_vars[idx];
    m_vars[idx] = nullptr;
  }

  return rt;
}

template<class DomainType>
void
DWDatabase<DomainType>::setScrubCount(const VarLabel* label,
                                      int matlIndex,
                                      const DomainType* dom,
                                      int count)
{
  std::lock_guard<MasterLock> set_scrub_count_lock(g_keyDB_lock);

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    SCI_THROW(UnknownVariable(label->getName(),
                              -99,
                              dom,
                              matlIndex,
                              "DWDatabase::setScrubCount",
                              __FILE__,
                              __LINE__));
  }
  m_scrubs[idx] = count;
}

template<class DomainType>
void
DWDatabase<DomainType>::scrub(const VarLabel* label,
                              int matlIndex,
                              const DomainType* dom)
{
  ASSERT(matlIndex >= -1);

  std::lock_guard<MasterLock> scrub_lock(g_keyDB_lock);

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx != -1 && m_vars[idx]) {
    delete m_vars[idx];
    m_vars[idx] = nullptr;
  }
}

template<class DomainType>
void
DWDatabase<DomainType>::initializeScrubs(
  int dwid,
  const FastHashTable<ScrubItem>* scrubcounts,
  bool add)
{
  // loop over each variable, probing the scrubcount map. Set the
  // scrubcount appropriately.  if the variable has no entry in
  // the scrubcount map, delete it
  for (auto keyiter = m_keyDB->m_keys.begin();
       keyiter != m_keyDB->m_keys.end();) {
    if (m_vars[keyiter->second]) {
      VarLabelMatl<DomainType> vlm = keyiter->first;
      // See if it is in the scrubcounts map.
      ScrubItem key(vlm.m_label, vlm.m_matl_index, vlm.m_domain, dwid);
      ScrubItem* result = scrubcounts->lookup(&key);
      if (!result && !add) {
        delete m_vars[keyiter->second];
        m_vars[keyiter->second] = nullptr;
      } else {
        if (result) {
          if (add) {
            __sync_add_and_fetch(&(m_scrubs[keyiter->second]), result->count);
          } else {
            if (!__sync_bool_compare_and_swap(&(m_scrubs[keyiter->second]),
                                              0,
                                              result->count)) {
              SCI_THROW(InternalError("initializing non-zero scrub counter",
                                      __FILE__,
                                      __LINE__));
            }
          }
        }
        ++keyiter;
      }
    } else {
      ++keyiter;
    }
  }
}

template<class DomainType>
void
DWDatabase<DomainType>::doReserve(KeyDatabase<DomainType>* keydb)
{
  m_keyDB = keydb;
  m_vars.resize(m_keyDB->m_key_count + 1,
                (DataItem*)nullptr);
  m_scrubs.resize(m_keyDB->m_key_count + 1, 0);
}

template<class DomainType>
bool
DWDatabase<DomainType>::exists(const VarLabel* label,
                               int matlIndex,
                               const DomainType* dom) const
{
  std::lock_guard<MasterLock> exists_lock(g_keyDB_lock);

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    return false;
  }
  if (m_vars[idx] == nullptr) {
    return false;
  }
  return true;
}

template<class DomainType>
void
DWDatabase<DomainType>::put(const VarLabel* label,
                            int matlIndex,
                            const DomainType* dom,
                            Variable* var,
                            bool init,
                            bool replace)
{
  ASSERT(matlIndex >= -1);

  std::lock_guard<MasterLock> put_lock(g_keyDB_lock);

  if (init) {
    m_keyDB->insert(label, matlIndex, dom);
    this->doReserve(m_keyDB);
  }

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    SCI_THROW(UnknownVariable(label->getName(),
                              -1,
                              dom,
                              matlIndex,
                              "DWDatabase::put",
                              __FILE__,
                              __LINE__));
  }

  if (m_vars[idx]) {
    if (m_vars[idx]->next) {
      SCI_THROW(
        InternalError("More than one vars on this label", __FILE__, __LINE__));
    }
    if (!replace) {
      SCI_THROW(InternalError("Put replacing old vars", __FILE__, __LINE__));
    }
    ASSERT(m_vars[idx]->var != var);
    delete m_vars[idx];
  }

  DataItem* newdi = new DataItem();
  newdi->var      = var;
  m_vars[idx]     = newdi;
}

//______________________________________________________________________
//
template<class DomainType>
void
DWDatabase<DomainType>::putReduce(const VarLabel* label,
                                  int matlIndex,
                                  const DomainType* dom,
                                  ReductionVariableBase* var,
                                  bool init)
{
  ASSERT(matlIndex >= -1);

  if (init) {
    m_keyDB->insert(label, matlIndex, dom);
    this->doReserve(m_keyDB);
  }
  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    SCI_THROW(UnknownVariable(label->getName(),
                              -1,
                              dom,
                              matlIndex,
                              "check task computes",
                              __FILE__,
                              __LINE__));
  }

  DataItem* newdi = new DataItem();
  newdi->var      = var;
  do {
    DataItem* olddi = __sync_lock_test_and_set(&m_vars[idx], 0);
    if (olddi == nullptr) {
      olddi = newdi;
    } else {
      ReductionVariableBase* oldvar =
        dynamic_cast<ReductionVariableBase*>(olddi->var);
      ReductionVariableBase* newvar =
        dynamic_cast<ReductionVariableBase*>(newdi->var);
      oldvar->reduce(*newvar);
      delete newdi;
    }
    newdi = __sync_lock_test_and_set(&m_vars[idx], olddi);
  } while (newdi != nullptr);
}

template<class DomainType>
void
DWDatabase<DomainType>::putForeign(const VarLabel* label,
                                   int matlIndex,
                                   const DomainType* dom,
                                   Variable* var,
                                   bool init)
{
  ASSERT(matlIndex >= -1);

  std::lock_guard<MasterLock> put_foreign_lock(g_keyDB_lock);

  if (init) {
    m_keyDB->insert(label, matlIndex, dom);
    this->doReserve(m_keyDB);
  }

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    SCI_THROW(UnknownVariable(label->getName(),
                              -1,
                              dom,
                              matlIndex,
                              "DWDatabase::putForeign",
                              __FILE__,
                              __LINE__));
  }

  DataItem* newdi = new DataItem();
  newdi->var    = var;
  do {
    newdi->next = m_vars[idx];
  } while (!__sync_bool_compare_and_swap(&m_vars[idx],
                                         newdi->next,
                                         newdi)); // vars[iter->second] = newdi;
}

template<class DomainType>
typename DWDatabase<DomainType>::DataItem*
DWDatabase<DomainType>::getDataItem(const VarLabel* label,
                                    int matlIndex,
                                    const DomainType* dom) const
{
  ASSERT(matlIndex >= -1);

  int idx = m_keyDB->lookup(label, matlIndex, dom);
  if (idx == -1) {
    SCI_THROW(UnknownVariable(label->getName(),
                              -99,
                              dom,
                              matlIndex,
                              "DWDatabase::getDataItem",
                              __FILE__,
                              __LINE__));
  }
  return m_vars[idx];
}

template<class DomainType>
inline Variable*
DWDatabase<DomainType>::get(const VarLabel* label,
                            int matlIndex,
                            const DomainType* dom) const
{
  std::lock_guard<MasterLock> get_lock(g_keyDB_lock);

  const DataItem* dataItem = getDataItem(label, matlIndex, dom);
  ASSERT(dataItem != nullptr);         // should have thrown an exception before
  ASSERT(dataItem->next == nullptr); // should call getlist()
  return dataItem->var;
}

template<class DomainType>
void
DWDatabase<DomainType>::get(const VarLabel* label,
                            int matlIndex,
                            const DomainType* dom,
                            Variable& var) const
{
  Variable* tmp = get(label, matlIndex, dom);
  var.copyPointer(*tmp);
}

template<class DomainType>
void
DWDatabase<DomainType>::getlist(const VarLabel* label,
                                int matlIndex,
                                const DomainType* dom,
                                std::vector<Variable*>& varlist) const
{
  std::lock_guard<MasterLock> get_list_lock(g_keyDB_lock);

  for (DataItem* dataItem = getDataItem(label, matlIndex, dom);
       dataItem != nullptr;
       dataItem = dataItem->next) {
    varlist.push_back(dataItem->var);
  }
}

template<class DomainType>
void
DWDatabase<DomainType>::print(int rank) const
{
  DOUT(true,
       "Rank-" << rank << " __________________________________DWDatabase ")

  for (auto keyiter = m_keyDB->m_keys.begin(); keyiter != m_keyDB->m_keys.end();
       ++keyiter) {
    if (m_vars[keyiter->second]) {
      const VarLabelMatl<DomainType>& vlm = keyiter->first;
      const DomainType* dom               = vlm.m_domain;

      if (dom) {
        DOUT(true,
             "Rank-" << rank << " Name: " << vlm.m_label->getName()
                     << "  domain: " << *dom << "  matl:" << vlm.m_matl_index);
      } else {
        DOUT(true,
             "Rank-" << rank << " Name: " << vlm.m_label->getName()
                     << "  domain: nullptr  matl: " << vlm.m_matl_index);
      }
    }
  }
  DOUT(true, "Rank-" << rank << " __________________________________")
}

template<class DomainType>
void
DWDatabase<DomainType>::logMemoryUse(std::ostream& out,
                                     unsigned long& total,
                                     const std::string& tag,
                                     int dwid)
{
  for (auto keyiter = m_keyDB->m_keys.begin(); keyiter != m_keyDB->m_keys.end();
       ++keyiter) {
    if (m_vars[keyiter->second]) {
      Variable* var                = m_vars[keyiter->second]->var;
      VarLabelMatl<DomainType> vlm = keyiter->first;
      const VarLabel* label        = vlm.m_label;
      std::string elems;
      unsigned long totsize;
      void* ptr;
      var->getSizeInfo(elems, totsize, ptr);
      const TypeDescription* td = label->typeDescription();

      logMemory(out,
                total,
                tag,
                label->getName(),
                (td ? td->getName() : "-"),
                vlm.m_domain,
                vlm.m_matl_index,
                elems,
                totsize,
                ptr,
                dwid);
    }
  }
}

template<class DomainType>
void
DWDatabase<DomainType>::getVarLabelMatlTriples(
  std::vector<VarLabelMatl<DomainType>>& v) const
{
  std::lock_guard<MasterLock> get_var_label_mat_triples_lock(
    g_keyDB_lock);

  for (auto keyiter = m_keyDB->m_keys.begin(); keyiter != m_keyDB->m_keys.end();
       ++keyiter) {
    const VarLabelMatl<DomainType>& vlm = keyiter->first;
    if (m_vars[keyiter->second]) {
      v.push_back(vlm);
    }
  }
}

} // End namespace Uintah

//______________________________________________________________________
//
// Custom hash function for VarLabelMatl
//
// Specialize std::hash structure and inject into the std namespace so that
// VarLabelMatl<DomainType> can be used as a key in std::unordered_map
//
// NOTE: this is legit, and likely the easiest way to get this done due to
// templates (APH - 10/25/18).
//
// It is allowed to add template specializations for any standard library class
// template to the std namespace only if the declaration depends on at least one
// program-defined type and the specialization satisfies all requirements for
// the original template
// (https://en.cppreference.com/w/cpp/language/extending_std).
//
namespace std {

template<class DomainType>
struct hash<Uintah::VarLabelMatl<DomainType>>
{
  size_t
  operator()(const Uintah::VarLabelMatl<DomainType>& v) const
  {
    size_t h  = 0u;
    char* str = const_cast<char*>(v.m_label->getName().data());
    while (int c = *str++) {
      h = h * 7 + c;
    }
    return ((((size_t)v.m_label) << (sizeof(size_t) / 2) ^
             ((size_t)v.m_label) >> (sizeof(size_t) / 2)) ^
            (size_t)v.m_domain ^ (size_t)v.m_matl_index);
  }
};

} // end namespace std

#endif
