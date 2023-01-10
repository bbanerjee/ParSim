/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#ifndef __CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_COMPTABLE_H__
#define __CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_COMPTABLE_H__

#include <Core/Grid/Task.h>

namespace Uintah {

class DetailedTask;
class DetailedTasks;
class Patch;

class CompTable
{
  struct Data
  {
    unsigned int
    string_hash(const char* p)
    {
      unsigned int sum = 0;
      while (*p) {
        sum = sum * 7 + (unsigned char)*p++;
      }
      return sum;
    }

    Data(DetailedTask* dtask,
         Task::Dependency* comp,
         const Patch* patch,
         int matl)
      : m_dtask(dtask)
      , m_comp(comp)
      , m_patch(patch)
      , m_matl(matl)
    {
      m_hash =
        (unsigned int)(((unsigned int)comp->mapDataWarehouse() << 3) ^
                       (string_hash(comp->m_var->getName().c_str())) ^ matl);
      if (patch) {
        m_hash ^= (unsigned int)(patch->getID() << 4);
      }
    }

    ~Data() {}

    bool
    operator==(const Data& c)
    {
      return m_matl == c.m_matl && m_patch == c.m_patch &&
             m_comp->m_reduction_level == c.m_comp->m_reduction_level &&
             m_comp->mapDataWarehouse() == c.m_comp->mapDataWarehouse() &&
             m_comp->m_var->equals(c.m_comp->m_var);
    }

    Data* m_next{ nullptr };
    DetailedTask* m_dtask{ nullptr };
    Task::Dependency* m_comp{ nullptr };
    const Patch* m_patch{ nullptr };
    int m_matl{};
    unsigned int m_hash{};
  };

  FastHashTable<Data> m_data{};

  void
  insert(Data* data);

public:
  CompTable() = default;

  ~CompTable() = default;

  void
  remembercomp(DetailedTask* dtask,
               Task::Dependency* comp,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               const ProcessorGroup* pg);

  bool
  findcomp(Task::Dependency* req,
           const Patch* patch,
           int matlIndex,
           DetailedTask*& dtask,
           Task::Dependency*& comp,
           const ProcessorGroup* pg);

  bool
  findReductionComps(Task::Dependency* req,
                     const Patch* patch,
                     int matlIndex,
                     std::vector<DetailedTask*>& dt,
                     const ProcessorGroup* pg);

private:
  void
  remembercomp(Data* newData, const ProcessorGroup* pg);
}; // class CompTable

#endif //__CCA_COMPONENTS_SCHEDULERS_TASKGRAPH_COMPTABLE_H__