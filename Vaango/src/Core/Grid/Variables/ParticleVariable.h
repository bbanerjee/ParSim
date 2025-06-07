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

#ifndef __CORE_GRID_VARIABLES_PARTICLEVARIABLE_H__
#define __CORE_GRID_VARIABLES_PARTICLEVARIABLE_H__

#include <Core/Grid/Variables/ParticleVariableBase.h>

#include <CCA/Ports/InputContext.h>
#include <CCA/Ports/OutputContext.h>
#include <CCA/Ports/PIDXOutputContext.h>

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/MPMIntVarTypes.h>
#include <Core/Grid/Variables/ParticleData.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/constGridVariable.h>

#include <Core/Malloc/Allocator.h>

#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/Assert.h>
#include <Core/Util/Endian.h>
#include <Core/Util/FancyAssert.h>

#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>

#include <cstring>
#include <iostream>

namespace Uintah {

class ProcessorGroup;
class TypeDescription;

template<class T>
class ParticleVariable : public ParticleVariableBase
{
  friend class constVariable<ParticleVariableBase,
                             ParticleVariable<T>,
                             T,
                             particleIndex>;

public:
  ParticleVariable();
  virtual ~ParticleVariable();
  ParticleVariable(ParticleSubset* pset);
  ParticleVariable(ParticleData<T>*, ParticleSubset* pset);

  static const TypeDescription*
  getTypeDescription();

  void
  resync()
  {
    d_pdata->resize(getParticleSubset()->numParticles());
  }

  virtual ParticleVariableBase*
  clone() override;

  virtual const ParticleVariableBase*
  clone() const;

  virtual ParticleVariableBase*
  cloneSubset(ParticleSubset*) override;

  virtual const ParticleVariableBase*
  cloneSubset(ParticleSubset*) const;

  virtual ParticleVariableBase*
  cloneType() const override
  {
    return scinew ParticleVariable<T>();
  }

  virtual constParticleVariableBase*
  cloneConstType() const override
  {
    return scinew constVariable<ParticleVariableBase,
                                ParticleVariable<T>,
                                T,
                                particleIndex>();
  }

  void
  copyData(const ParticleVariable<T>& src);

  virtual void
  copyData(const ParticleVariableBase* src) override
  {
    copyData(castFromBase(src));
  }

  inline T&
  operator[](particleIndex idx)
  {
    ASSERTRANGE(idx, 0, (particleIndex)d_pdata->size);
    return d_pdata->data[idx];
  }

  //////////
  // Insert Documentation Here:
  inline const T&
  operator[](particleIndex idx) const
  {
    ASSERTRANGE(idx, 0, (particleIndex)d_pdata->size);
    return d_pdata->data[idx];
  }

  virtual void
  copyPointer(ParticleVariable<T>&);

  virtual void
  copyPointer(Variable&) override;

  virtual void
  allocate(ParticleSubset*) override;

  virtual void
  allocate(int totalParticles) override;

  virtual void
  allocate(const Patch*, const IntVector& /*boundary*/) override
  {
    SCI_THROW(
      InternalError("Should not call ParticleVariable<T>::allocate(const "
                    "Patch*), use allocate(ParticleSubset*) instead.",
                    __FILE__,
                    __LINE__));
  }

  virtual int
  size() override
  {
    return d_pdata->size;
  }

  // specialized for T=Point
  virtual void
  gather(ParticleSubset* dest,
         const std::vector<ParticleSubset*>& subsets,
         const std::vector<ParticleVariableBase*>& srcs,
         const std::vector<const Patch*>& /*srcPatches*/,
         particleIndex extra = 0) override;

  virtual void
  gather(ParticleSubset* dest,
         const std::vector<ParticleSubset*>& subsets,
         const std::vector<ParticleVariableBase*>& srcs,
         particleIndex extra = 0) override;

  virtual void
  unpackMPI(void* buf,
            int bufsize,
            int* bufpos,
            const ProcessorGroup* pg,
            ParticleSubset* pset) override;

  virtual void
  packMPI(void* buf,
          int bufsize,
          int* bufpos,
          const ProcessorGroup* pg,
          ParticleSubset* pset) override;

  // specialized for T=Point
  virtual void
  packMPI(void* buf,
          int bufsize,
          int* bufpos,
          const ProcessorGroup* pg,
          ParticleSubset* pset,
          const Patch* /*forPatch*/) override;

  virtual void
  packsizeMPI(int* bufpos,
              const ProcessorGroup* pg,
              ParticleSubset* pset) override;

  virtual void
  emitNormal(std::ostream& out,
             const IntVector& l,
             const IntVector& h,
             ProblemSpecP varnode,
             bool outputDoubleAsFloat) override;

#ifdef HAVE_PIDX
  virtual void
  emitPIDX(
    PIDXOutputContext& oc,
    unsigned char* buffer,
    const IntVector& /* l */,
    const IntVector& /* h */,
    const size_t pidx_bufferSize); // buffer size used for bullet proofing.
#endif

  virtual void
  readNormal(std::istream& in, bool swapBytes) override;

  virtual void*
  getBasePointer() const override;

  virtual const TypeDescription*
  virtualGetTypeDescription() const override;

  virtual RefCounted*
  getRefCounted() override
  {
    return d_pdata;
  }

  virtual void
  getSizeInfo(std::string& elems,
              unsigned long& totsize,
              void*& ptr) const override
  {
    std::ostringstream str;
    str << getParticleSubset()->numParticles();
    elems   = str.str();
    totsize = getParticleSubset()->numParticles() * sizeof(T);
    ptr     = getBasePointer();
  }

  virtual size_t
  getDataSize() const override
  {
    return getParticleSubset()->numParticles() * sizeof(T);
  }

  virtual bool
  copyOut(void* dst) const override
  {
    void* src       = (void*)this->getBasePointer();
    size_t numBytes = getDataSize();
    void* retVal    = std::memcpy(dst, src, numBytes);
    return (retVal == dst) ? true : false;
  }

protected:
  static TypeDescription* td;
  ParticleVariable(const ParticleVariable<T>&);

  ParticleVariable<T>&
  operator=(const ParticleVariable<T>&);

private:
  ParticleData<T>* d_pdata;
  Vector offset_; // only used when T is Point

  static const ParticleVariable<T>&
  castFromBase(const ParticleVariableBase* srcptr);

  static TypeDescription::Register registerMe;

  static Variable*
  maker();
};

template<class T>
TypeDescription* ParticleVariable<T>::td = 0;

template<class T>
TypeDescription::Register ParticleVariable<T>::registerMe(getTypeDescription());

template<class T>
const TypeDescription*
ParticleVariable<T>::getTypeDescription()
{
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::ParticleVariable,
                                "ParticleVariable",
                                &maker,
                                fun_getTypeDescription((T*)0));
  }
  return td;
}

template<class T>
Variable*
ParticleVariable<T>::maker()
{
  return scinew ParticleVariable<T>();
}

template<class T>
ParticleVariable<T>::ParticleVariable()
  : ParticleVariableBase(0)
  , d_pdata(0)
{
}

template<class T>
ParticleVariable<T>::~ParticleVariable()
{
  if (d_pdata && d_pdata->removeReference()) {
    delete d_pdata;
  }
}

template<class T>
ParticleVariable<T>::ParticleVariable(ParticleSubset* pset)
  : ParticleVariableBase(pset)
{
  d_pdata = scinew ParticleData<T>(pset->numParticles());
  d_pdata->addReference();
}

template<class T>
void
ParticleVariable<T>::allocate(int totalParticles)
{
  ASSERT(isForeign());
  ASSERT(d_pset == 0);

  // this is a pset-less storage as it could have several.  Should be used for
  // foreign data only.  To iterate over particles in this pset, use gather
  d_pdata = scinew ParticleData<T>(totalParticles);
  d_pdata->addReference();
}

template<class T>
void
ParticleVariable<T>::allocate(ParticleSubset* pset)
{
  if (d_pdata && d_pdata->removeReference()) {
    delete d_pdata;
  }
  if (d_pset && d_pset->removeReference()) {
    delete d_pset;
  }

  d_pset = pset;
  d_pset->addReference();
  d_pdata = scinew ParticleData<T>(pset->numParticles());
  d_pdata->addReference();
}

template<class T>
ParticleVariableBase*
ParticleVariable<T>::clone()
{
  return scinew ParticleVariable<T>(*this);
}

template<class T>
const ParticleVariableBase*
ParticleVariable<T>::clone() const
{
  return scinew ParticleVariable<T>(*this);
}

template<class T>
ParticleVariableBase*
ParticleVariable<T>::cloneSubset(ParticleSubset* pset)
{
  return scinew ParticleVariable<T>(d_pdata, pset);
}

template<class T>
const ParticleVariableBase*
ParticleVariable<T>::cloneSubset(ParticleSubset* pset) const
{
  return scinew ParticleVariable<T>(d_pdata, pset);
}

template<class T>
const ParticleVariable<T>&
ParticleVariable<T>::castFromBase(const ParticleVariableBase* srcptr)
{
  const ParticleVariable<T>* c =
    dynamic_cast<const ParticleVariable<T>*>(srcptr);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in Particle variable", __FILE__, __LINE__));
  }
  return *c;
}

template<class T>
void
ParticleVariable<T>::copyData(const ParticleVariable<T>& src)
{
  ASSERT(*d_pset == *src.d_pset);
  *d_pdata = *src.d_pdata;
}

template<class T>
ParticleVariable<T>::ParticleVariable(ParticleData<T>* pdata,
                                      ParticleSubset* pset)
  : ParticleVariableBase(pset)
  , d_pdata(pdata)
{
  if (d_pdata) {
    d_pdata->addReference();
  }
}

template<class T>
ParticleVariable<T>::ParticleVariable(const ParticleVariable<T>& copy)
  : ParticleVariableBase(copy)
  , d_pdata(copy.d_pdata)
{
  if (d_pdata) {
    d_pdata->addReference();
  }
}

template<class T>
void
ParticleVariable<T>::copyPointer(ParticleVariable<T>& copy)
{
  if (this != &copy) {
    ParticleVariableBase::operator=(copy);
    if (d_pdata && d_pdata->removeReference()) {
      delete d_pdata;
    }
    d_pdata = copy.d_pdata;
    if (d_pdata) {
      d_pdata->addReference();
    }
  }
}

template<class T>
void
ParticleVariable<T>::copyPointer(Variable& copy)
{
  ParticleVariable<T>* c = dynamic_cast<ParticleVariable<T>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in particle variable", __FILE__, __LINE__));
  }
  copyPointer(*c);
}

// specialization for T=Point
template<>
void
ParticleVariable<Point>::gather(ParticleSubset* pset,
                                const std::vector<ParticleSubset*>& subsets,
                                const std::vector<ParticleVariableBase*>& srcs,
                                const std::vector<const Patch*>& srcPatches,
                                particleIndex extra);

template<class T>
void
ParticleVariable<T>::gather(ParticleSubset* pset,
                            const std::vector<ParticleSubset*>& subsets,
                            const std::vector<ParticleVariableBase*>& srcs,
                            const std::vector<const Patch*>& /*srcPatches*/,
                            particleIndex extra)
{
  gather(pset, subsets, srcs, extra);
}

template<class T>
void
ParticleVariable<T>::gather(ParticleSubset* pset,
                            const std::vector<ParticleSubset*>& subsets,
                            const std::vector<ParticleVariableBase*>& srcs,
                            particleIndex extra)
{
  if (d_pdata && d_pdata->removeReference()) {
    delete d_pdata;
  }
  if (d_pset && d_pset->removeReference()) {
    delete d_pset;
  }
  d_pset = pset;
  pset->addReference();
  d_pdata = scinew ParticleData<T>(pset->numParticles());
  d_pdata->addReference();
  ASSERTEQ(subsets.size(), srcs.size());
  ParticleSubset::iterator dstiter = pset->begin();
  for (int i = 0; i < (int)subsets.size(); i++) {
    ParticleVariable<T>* srcptr = dynamic_cast<ParticleVariable<T>*>(srcs[i]);
    if (!srcptr) {
      SCI_THROW(TypeMismatchException(
        "Type mismatch in ParticleVariable::gather", __FILE__, __LINE__));
    }
    ParticleVariable<T>& src = *srcptr;
    ParticleSubset* subset   = subsets[i];
    for (ParticleSubset::iterator srciter = subset->begin();
         srciter != subset->end();
         srciter++) {
      (*this)[*dstiter] = src[*srciter];
      dstiter++;
    }
  }
  ASSERT(dstiter + extra == pset->end());
}

template<class T>
void*
ParticleVariable<T>::getBasePointer() const
{
  return &d_pdata->data[0];
}

template<class T>
const TypeDescription*
ParticleVariable<T>::virtualGetTypeDescription() const
{
  return getTypeDescription();
}

template<class T>
void
ParticleVariable<T>::unpackMPI(void* buf,
                               int bufsize,
                               int* bufpos,
                               const ProcessorGroup* pg,
                               ParticleSubset* pset)
{
  // This should be fixed for variable sized types!
  const TypeDescription* td = getTypeDescription()->getSubType();
  if (td->isFlat()) {
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
         iter++) {
      Uintah::MPI::Unpack(buf,
                          bufsize,
                          bufpos,
                          &d_pdata->data[*iter],
                          1,
                          td->getMPIType(),
                          pg->getComm());
    }
  } else {
    SCI_THROW(InternalError("packMPI not finished\n", __FILE__, __LINE__));
  }
}

// specialized for T=Point
template<>
void
ParticleVariable<Point>::packMPI(void* buf,
                                 int bufsize,
                                 int* bufpos,
                                 const ProcessorGroup* pg,
                                 ParticleSubset* pset,
                                 const Patch* forPatch);
template<class T>
void
ParticleVariable<T>::packMPI(void* buf,
                             int bufsize,
                             int* bufpos,
                             const ProcessorGroup* pg,
                             ParticleSubset* pset,
                             const Patch* /*forPatch*/)
{
  packMPI(buf, bufsize, bufpos, pg, pset);
}

template<class T>
void
ParticleVariable<T>::packMPI(void* buf,
                             int bufsize,
                             int* bufpos,
                             const ProcessorGroup* pg,
                             ParticleSubset* pset)
{
  // This should be fixed for variable sized types!
  const TypeDescription* td = getTypeDescription()->getSubType();
  if (td->isFlat()) {
    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
         iter++) {
      Uintah::MPI::Pack(&d_pdata->data[*iter],
                        1,
                        td->getMPIType(),
                        buf,
                        bufsize,
                        bufpos,
                        pg->getComm());
    }
  } else {
    SCI_THROW(InternalError("packMPI not finished\n", __FILE__, __LINE__));
  }
}

template<class T>
void
ParticleVariable<T>::packsizeMPI(int* bufpos,
                                 const ProcessorGroup* pg,
                                 ParticleSubset* pset)
{
  // This should be fixed for variable sized types!
  const TypeDescription* td = getTypeDescription()->getSubType();
  int n                     = pset->numParticles();
  if (td->isFlat()) {
    int size;
    Uintah::MPI::Pack_size(n, td->getMPIType(), pg->getComm(), &size);
    (*bufpos) += size;
  } else {
    SCI_THROW(InternalError("packsizeMPI not finished\n", __FILE__, __LINE__));
  }
}

// Specialized in ParticleVariable_special.cc
template<>
void
ParticleVariable<double>::emitNormal(std::ostream& out,
                                     const IntVector&,
                                     const IntVector&,
                                     ProblemSpecP varnode,
                                     bool outputDoubleAsFloat);

template<class T>
void
ParticleVariable<T>::emitNormal(std::ostream& out,
                                const IntVector&,
                                const IntVector&,
                                ProblemSpecP varnode,
                                bool /*outputDoubleAsFloat*/)
{
  const TypeDescription* td = fun_getTypeDescription((T*)0);

  if (varnode->findBlock("numParticles") == 0) {
    varnode->appendElement("numParticles", d_pset->numParticles());
  }
  if (!td->isFlat()) {
    SCI_THROW(InternalError(
      "Cannot yet write non-flat objects!\n", __FILE__, __LINE__));
  } else {
    // This could be optimized...
    ParticleSubset::iterator iter = d_pset->begin();
    while (iter != d_pset->end()) {
      particleIndex start = *iter;
      iter++;
      particleIndex end = start + 1;
      while (iter != d_pset->end() && *iter == end) {
        end++;
        iter++;
      }
      ssize_t size = (ssize_t)(sizeof(T) * (end - start));
      out.write((char*)&(*this)[start], size);
    }
  }
}

#ifdef HAVE_PIDX
template<class T>
void
ParticleVariable<T>::emitPIDX(
  PIDXOutputContext& oc,
  unsigned char* buffer,
  const IntVector& /* l */,
  const IntVector& /* h */,
  const size_t pidx_bufferSize) // buffer size used for bullet proofing.
{
  const TypeDescription* td = fun_getTypeDescription((T*)nullptr);

  // int numParticles = d_pset->numParticles();

  if (!td->isFlat()) { // Not certain what this is for?
    SCI_THROW(InternalError(
      "Cannot yet write non-flat objects!\n", __FILE__, __LINE__));
  } else {

    // FIXME: Add in assertion that pidx_bufferSize == "calculate the
    // d_pdata->size * size of element"

    memcpy(buffer, d_pdata->data, pidx_bufferSize);
  }
}
#endif

template<class T>
void
ParticleVariable<T>::readNormal(std::istream& in, bool swapBytes)
{
  const TypeDescription* td = fun_getTypeDescription((T*)0);
  if (!td->isFlat()) {
    SCI_THROW(
      InternalError("Cannot yet read non-flat objects!\n", __FILE__, __LINE__));
  } else {
    // This could be optimized...
    ParticleSubset::iterator iter = d_pset->begin();
    while (iter != d_pset->end()) {
      particleIndex start = *iter;
      iter++;
      particleIndex end = start + 1;
      while (iter != d_pset->end() && *iter == end) {
        end++;
        iter++;
      }
      ssize_t size = (ssize_t)(sizeof(T) * (end - start));
      in.read((char*)&(*this)[start], size);
      if (swapBytes) {
        for (particleIndex idx = start; idx != end; idx++) {
          Uintah::swapbytes((*this)[idx]);
        }
      }
    }
  }
}

template<class T>
class constParticleVariable
  : public constVariable<ParticleVariableBase,
                         ParticleVariable<T>,
                         T,
                         particleIndex>
{
public:
  constParticleVariable()
    : constVariable<ParticleVariableBase,
                    ParticleVariable<T>,
                    T,
                    particleIndex>()
  {
  }

  constParticleVariable(const ParticleVariable<T>& copy)
    : constVariable<ParticleVariableBase,
                    ParticleVariable<T>,
                    T,
                    particleIndex>(copy)
  {
  }

  ParticleSubset*
  getParticleSubset() const
  {
    return this->rep_.getParticleSubset();
  }
};

} // End namespace Uintah

#ifdef __PGI
#include <Core/Grid/Variables/ParticleVariable_special.cc>
#endif

#endif //__CORE_GRID_VARIABLES_PARTICLEVARIABLE_H__
