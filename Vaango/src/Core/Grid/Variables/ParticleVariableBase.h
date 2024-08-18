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

#ifndef __CORE_GRID_VARIABLES_ParticleVariableBase_H__
#define __CORE_GRID_VARIABLES_ParticleVariableBase_H__

#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/Variable.h>
#include <Core/Grid/Variables/constVariable.h>

#include <Core/Util/Assert.h>

#include <vector>

namespace Uintah {

class BufferInfo;
class OutputContext;
class ParticleSubset;
class Patch;
class ProcessorGroup;
class TypeDescription;

using constParticleVariableBase = constVariableBase<ParticleVariableBase>;

class ParticleVariableBase : public Variable
{
public:
  virtual ~ParticleVariableBase();

  virtual ParticleVariableBase*
  clone() = 0;

  virtual ParticleVariableBase*
  cloneSubset(ParticleSubset*) = 0;

  // Make a new default object of the base class.
  virtual ParticleVariableBase*
  cloneType() const = 0;

  virtual constParticleVariableBase*
  cloneConstType() const = 0;

  // not something we normally do, but helps in AMR when we copy
  // data from one patch to another where they are the same boundaries
  // instead of copying all the data
  void
  setParticleSubset(ParticleSubset* pset);

  virtual void
  copyData(const ParticleVariableBase* src) = 0;

  virtual void
  allocate(const Patch*,
           const Uintah::IntVector& boundary)
    override = 0; // will throw an InternalError

  virtual void
  allocate(ParticleSubset*) = 0;

  virtual void
  allocate(int totalParticles) = 0;

  virtual void
  gather(ParticleSubset* dest,
         const std::vector<ParticleSubset*>& subsets,
         const std::vector<ParticleVariableBase*>& srcs,
         particleIndex extra = 0) = 0;

  virtual void
  gather(ParticleSubset* dest,
         const std::vector<ParticleSubset*>& subsets,
         const std::vector<ParticleVariableBase*>& srcs,
         const std::vector<const Patch*>& srcPatches,
         particleIndex extra = 0) = 0;

  virtual void
  unpackMPI(void* buf,
            int bufsize,
            int* bufpos,
            const ProcessorGroup* pg,
            ParticleSubset* pset) = 0;

  virtual void
  packMPI(void* buf,
          int bufsize,
          int* bufpos,
          const ProcessorGroup* pg,
          ParticleSubset* pset) = 0;

  virtual void
  packMPI(void* buf,
          int bufsize,
          int* bufpos,
          const ProcessorGroup* pg,
          ParticleSubset* pset,
          const Patch* forPatch) = 0;

  virtual void
  packsizeMPI(int* bufpos, const ProcessorGroup* pg, ParticleSubset* pset) = 0;

  virtual int
  size() = 0;

  virtual size_t
  getDataSize() const override = 0;

  virtual bool
  copyOut(void* dst) const override = 0;

  ParticleSubset*
  getParticleSubset() const
  {
    ASSERT(!isForeign());
    return d_pset;
  }

  virtual void*
  getBasePointer() const = 0;

  void
  getMPIBuffer(BufferInfo& buffer, ParticleSubset* sendset);

  virtual const TypeDescription*
  virtualGetTypeDescription() const override = 0;

  virtual RefCounted*
  getRefCounted() override = 0;

  virtual void
  getSizeInfo(std::string& elems,
              unsigned long& totsize,
              void*& ptr) const override = 0;

protected:
  ParticleVariableBase(const ParticleVariableBase&);

  ParticleVariableBase(ParticleSubset* pset);

  ParticleVariableBase&
  operator=(const ParticleVariableBase&);

  ParticleSubset* d_pset;

private:
};

} // End namespace Uintah

#endif // __CORE_GRID_VARIABLES_ParticleVariableBase_H__
