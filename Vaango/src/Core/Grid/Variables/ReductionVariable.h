/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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

#ifndef __CORE_GRID_VARIABLES_REDUCTION_VARIABLE_H__
#define __CORE_GRID_VARIABLES_REDUCTION_VARIABLE_H__

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>
#include <Core/Grid/Variables/DataItem.h>
#include <Core/Grid/Variables/ReductionVariableBase.h>
#include <Core/Grid/Variables/Reductions.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/Endian.h>

#include <cstring>
#include <iosfwd>
#include <iostream>
#include <memory>
#include <sstream>

namespace Uintah {

template<class T, class Op>
class ReductionVariable : public ReductionVariableBase
{
public:
  inline ReductionVariable()
    : ReductionVariableBase()
    , d_value(std::make_shared<T>())
  {
  }

  inline ReductionVariable(T value)
    : ReductionVariableBase()
    , d_value(std::make_shared<T>(value))
  {
  }

  inline ReductionVariable(const ReductionVariable<T, Op>& copy)
    : ReductionVariableBase()
    , d_value(copy.d_value)
  {
  }

  virtual void
  copyPointer(Variable&);

  virtual ~ReductionVariable(){};

  virtual const TypeDescription*
  virtualGetTypeDescription() const
  {
    return getTypeDescription();
  }

  static const TypeDescription*
  getTypeDescription();

  inline operator T() const { return *d_value; }
  inline T&
  get()
  {
    return *d_value;
  }
  inline const T&
  get() const
  {
    return *d_value;
  }

  void
  setData(const T& val)
  {
    d_value = std::make_shared<T>(val);
  };

  virtual ReductionVariableBase*
  clone() const
  {
    return scinew ReductionVariable<T, Op>(*this);
  }

private:
  ReductionVariable<T, Op>&
  operator=(const ReductionVariable<T, Op>& copy)
  {
    d_value = copy.d_value;
    return *this;
  };

public:
  virtual void
  getSizeInfo(std::string& elems, unsigned long& totsize, void*& ptr) const
  {
    elems   = "1";
    totsize = sizeof(T);
    ptr     = getBasePointer();
  }

  virtual size_t
  getDataSize() const
  {
    return sizeof(T);
  }

  virtual void*
  getBasePointer() const
  {
    return d_value.get();
  }

  virtual bool
  copyOut(void* dst) const
  {
    void* src       = (void*)(&d_value);
    size_t numBytes = getDataSize();
    void* retVal    = std::memcpy(dst, src, numBytes);
    return (retVal == dst);
  }

  virtual void
  emitNormal(std::ostream& out,
             const IntVector& /*l*/,
             const IntVector& /*h*/,
             ProblemSpecP /*varnode*/,
             bool /*outputDoubleAsFloat*/)
  {
    ssize_t linesize = (ssize_t)(sizeof(T));
    out.write((char*)(d_value.get()), linesize);
  }

  virtual void
  readNormal(std::istream& in, bool swapBytes)
  {
    ssize_t linesize = (ssize_t)(sizeof(T));
    T val;
    in.read((char*)&val, linesize);

    if (swapBytes) {
      Uintah::swapbytes(val);
    }

    d_value = std::make_shared<T>(val);
  }

  virtual void
  print(std::ostream& out) const
  {
    out << *(d_value.get());
  }

  virtual void
  reduce(const ReductionVariableBase&);

  virtual void
  getMPIInfo(int& count, MPI_Datatype& datatype, MPI_Op& op);
  virtual void
  getMPIData(std::vector<char>& buf, int& index);
  virtual void
  putMPIData(std::vector<char>& buf, int& index);

  //! Sets the value to a harmless value that will have no impact
  //! on a reduction.
  virtual void
  setBenignValue()
  {
    Op op;
    d_value = std::make_shared<T>(op.getBenignValue());
  }

  // check if the value is benign value
  virtual bool
  isBenignValue() const
  {
    Op op;
    return (*(d_value.get()) == op.getBenignValue());
  }

  // Static variable whose entire purpose is to cause the
  // (instantiated) type of this class to be registered with the
  // Core/Disclosure/TypeDescription class when this class' object
  // code is originally loaded from the shared library.  The
  // 'registerMe' variable is not used for anything else in the
  // program.
  static TypeDescription::Register s_registerMe;

private:
  static TypeDescription* s_td;
  static Variable*
  maker()
  {
    return scinew ReductionVariable<T, Op>();
  }

private:
  std::shared_ptr<T> d_value;
};

template<class T, class Op>
TypeDescription* ReductionVariable<T, Op>::s_td = nullptr;

// The following line is the initialization (creation) of the
// 'registerMe' static variable (for each version of CCVariable
// (double, int, etc)).  Note, the 'registerMe' variable is created
// when the object code is initially loaded (usually during intial
// program load by the operating system).
template<class T, class Op>
TypeDescription::Register ReductionVariable<T, Op>::s_registerMe(
  getTypeDescription());

template<class T, class Op>
const TypeDescription*
ReductionVariable<T, Op>::getTypeDescription()
{
  if (!s_td) {
    // this is a hack to get a non-null ReductionVariable var for some
    // functions the ReductionVariables are used in (i.e., task->computes).
    // Since they're not fully-qualified variables, maker would fail
    // anyway.  And since most instances use Handle, it would be
    // difficult.
    T* tmp = nullptr;
    s_td   = scinew TypeDescription(TypeDescription::Type::ReductionVariable,
                                  "ReductionVariable",
                                  &maker,
                                  fun_getTypeDescription(tmp));
  }
  return s_td;
}

template<class T, class Op>
void
ReductionVariable<T, Op>::copyPointer(Variable& copy)
{
  const ReductionVariable<T, Op>* c =
    dynamic_cast<const ReductionVariable<T, Op>*>(&copy);
  if (!c) {
    std::ostringstream out;
    out << "Type mismatch in reduction variable";
    throw TypeMismatchException(out.str(), __FILE__, __LINE__);
  }
  *this = *c;
}

template<class T, class Op>
void
ReductionVariable<T, Op>::reduce(const ReductionVariableBase& other)
{
  const ReductionVariable<T, Op>* c =
    dynamic_cast<const ReductionVariable<T, Op>*>(&other);
  if (!c) {
    std::ostringstream out;
    out << "Type mismatch in reduction variable";
    throw TypeMismatchException(out.str(), __FILE__, __LINE__);
  }
  Op op;
  T val   = op(*(d_value.get()), *(c->d_value.get()));
  d_value = std::make_shared<T>(val);
}

} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_REDUCTION_VARIABLE_H__
