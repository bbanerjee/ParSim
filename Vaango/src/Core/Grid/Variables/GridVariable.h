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

#ifndef __CORE_GRID_VARIABLES_GRIDVARIABLE_H__
#define __CORE_GRID_VARIABLES_GRIDVARIABLE_H__

#include <Core/Grid/Variables/Array3.h>
#include <Core/Grid/Variables/GridVariableBase.h>

#include <CCA/Ports/InputContext.h>
#include <CCA/Ports/OutputContext.h>
#include <CCA/Ports/PIDXOutputContext.h>

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>

#include <Core/Geometry/Vector.h>
#include <Core/Malloc/Allocator.h>

#include <cstring>

namespace Uintah {

class TypeDescription;

template<class T>
class GridVariable
  : public GridVariableBase
  , public Array3<T>
{
public:
  GridVariable()
    : GridVariableBase()
  {
  }
  virtual ~GridVariable() = default;

  inline void
  copyPointer(GridVariable<T>& copy)
  {
    Array3<T>::copyPointer(copy);
  }

  virtual void
  copyPointer(Variable&) override;

  virtual bool
  rewindow(const IntVector& low, const IntVector& high) override
  {
    return Array3<T>::rewindow(low, high);
  }

  virtual void
  offset(const IntVector& offset) override
  {
    Array3<T>::offset(offset);
  }

  // offset the indexing into the array (useful when getting virtual
  // patch data -- i.e. for periodic boundary conditions)
  virtual void
  offsetGrid(const IntVector& offset) override
  {
    Array3<T>::offset(offset);
  }

  static const GridVariable<T>&
  castFromBase(const GridVariableBase* srcptr);

  virtual void
  allocate(const IntVector& lowIndex, const IntVector& highIndex) override;

  void
  copyPatch(const GridVariable<T>& src,
            const IntVector& lowIndex,
            const IntVector& highIndex);

  virtual void
  copyPatch(const GridVariableBase* src,
            const IntVector& lowIndex,
            const IntVector& highIndex) override
  {
    copyPatch(castFromBase(src), lowIndex, highIndex);
  }

  void
  copyData(const GridVariable<T>& src)
  {
    copyPatch(src, src.getLowIndex(), src.getHighIndex());
  }

  virtual void
  copyData(const GridVariableBase* src) override
  {
    copyPatch(src, src->getLow(), src->getHigh());
  }

  virtual void*
  getBasePointer() const override
  {
    return (void*)this->getPointer();
  }

  virtual void
  getSizes(IntVector& low, IntVector& high, IntVector& siz) const override;

  virtual void
  getSizes(IntVector& low,
           IntVector& high,
           IntVector& dataLow,
           IntVector& siz,
           IntVector& strides) const override;

  virtual void
  getSizeInfo(std::string& elems,
              unsigned long& totsize,
              void*& ptr) const override
  {
    IntVector siz = this->size();
    std::ostringstream str;
    str << siz.x() << "x" << siz.y() << "x" << siz.z();
    elems   = str.str();
    totsize = siz.x() * siz.y() * siz.z() * sizeof(T);
    ptr     = (void*)this->getPointer();
  }

  virtual size_t
  getDataSize() const override
  {
    IntVector siz = this->size();
    return siz.x() * siz.y() * siz.z() * sizeof(T);
  }

  virtual bool
  copyOut(void* dst) const override
  {
    void* src       = (void*)this->getPointer();
    size_t numBytes = getDataSize();
    void* retVal    = std::memcpy(dst, src, numBytes);
    return (retVal == dst) ? true : false;
  }

  virtual IntVector
  getLow() const override
  {
    return this->getLowIndex();
  }

  virtual IntVector
  getHigh() const override
  {
    return this->getHighIndex();
  }

#if HAVE_PIDX
  virtual void
  emitPIDX(PIDXOutputContext& pc,
           unsigned char* pidx_buffer,
           const IntVector& l,
           const IntVector& h,
           const size_t pidx_bufferSize)
  {
    // This seems inefficient  -Todd
    ProblemSpecP dummy;
    bool outputDoubleAsFloat = pc.isOutputDoubleAsFloat();

    // read the Array3 variable into tmpStream
    std::ostringstream tmpStream;

    emitNormal(tmpStream, l, h, dummy, outputDoubleAsFloat);

    // Create a string from the ostringstream
    std::string tmpString = tmpStream.str();

    // create a c-string
    const char* writebuffer = tmpString.c_str();
    size_t uda_bufferSize   = tmpString.size();

    // copy the write buffer into the pidx_buffer
    if (uda_bufferSize == pidx_bufferSize) {
      memcpy(pidx_buffer, writebuffer, uda_bufferSize);
    } else {
      std::cout << "uda_bufferSize: " << uda_bufferSize
                << ", pidx_bufferSize: " << pidx_bufferSize << "\n";

      throw InternalError("ERROR: Variable::emitPIDX() error reading in buffer",
                          __FILE__,
                          __LINE__);
    }
  }
#endif

  virtual void
  emitNormal(std::ostream& out,
             const IntVector& l,
             const IntVector& h,
             ProblemSpecP /*varnode*/,
             bool outputDoubleAsFloat) override
  {
    const TypeDescription* td = fun_getTypeDescription((T*)nullptr);
    if (td->isFlat()) {
      Array3<T>::write(out, l, h, outputDoubleAsFloat);
    } else {
      SCI_THROW(InternalError(
        "Cannot yet write non-flat objects!\n", __FILE__, __LINE__));
    }
  }

  virtual void
  readNormal(std::istream& in, bool swapBytes) override
  {
    const TypeDescription* td = fun_getTypeDescription((T*)0);
    if (td->isFlat()) {
      Array3<T>::read(in, swapBytes);
    } else {
      SCI_THROW(InternalError(
        "Cannot yet read non-flat objects!\n", __FILE__, __LINE__));
    }
  }

  virtual RefCounted*
  getRefCounted() override
  {
    return this->getWindow();
  }

protected:
  GridVariable(const GridVariable<T>& copy)
    : GridVariableBase()
    , Array3<T>(copy)
  {
  }

private:
  GridVariable(Array3Window<T>* window)
    : Array3<T>(window)
  {
  }

  GridVariable<T>&
  operator=(const GridVariable<T>&) = delete;
};

template<class T>
void
GridVariable<T>::copyPointer(Variable& copy)
{
  GridVariable<T>* c = dynamic_cast<GridVariable<T>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in Grid variable", __FILE__, __LINE__));
  }
  copyPointer(*c);
}

template<class T>
const GridVariable<T>&
GridVariable<T>::castFromBase(const GridVariableBase* srcptr)
{
  const GridVariable<T>* c = dynamic_cast<const GridVariable<T>*>(srcptr);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in CC variable", __FILE__, __LINE__));
  }
  return *c;
}

template<class T>
void
GridVariable<T>::allocate(const IntVector& lowIndex, const IntVector& highIndex)
{
  if (this->getWindow()) {
    SCI_THROW(InternalError("Allocating a Gridvariable that "
                            "is apparently already allocated!",
                            __FILE__,
                            __LINE__));
  }
  this->resize(lowIndex, highIndex);
}

template<class T>
void
GridVariable<T>::copyPatch(const GridVariable<T>& src,
                           const IntVector& lowIndex,
                           const IntVector& highIndex)
{
  if (this->getWindow()->getData() == src.getWindow()->getData() &&
      this->getWindow()->getOffset() == src.getWindow()->getOffset()) {
    // No copy needed
    return;
  }

#if 0
    for(int i=lowIndex.x();i<highIndex.x();i++)
      for(int j=lowIndex.y();j<highIndex.y();j++)
        for(int k=lowIndex.z();k<highIndex.z();k++)
          (*this)[IntVector(i, j, k)] = src[IntVector(i,j,k)];
#endif
  this->copy(src, lowIndex, highIndex);
}

template<class T>
void
GridVariable<T>::getSizes(IntVector& low, IntVector& high, IntVector& siz) const
{
  low  = this->getLowIndex();
  high = this->getHighIndex();
  siz  = this->size();
}

template<class T>
void
GridVariable<T>::getSizes(IntVector& low,
                          IntVector& high,
                          IntVector& dataLow,
                          IntVector& siz,
                          IntVector& strides) const
{
  low     = this->getLowIndex();
  high    = this->getHighIndex();
  dataLow = this->getWindow()->getOffset();
  siz     = this->size();
  strides = IntVector(sizeof(T),
                      (int)(sizeof(T) * siz.x()),
                      (int)(sizeof(T) * siz.y() * siz.x()));
}

} // end namespace Uintah

#endif //__CORE_GRID_VARIABLES_GRIDVARIABLE_H__
