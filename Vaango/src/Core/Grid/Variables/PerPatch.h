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

#ifndef __CORE_GRID_VARIABLES_PERPATCH_H__
#define __CORE_GRID_VARIABLES_PERPATCH_H__

#include <Core/Grid/Variables/PerPatchBase.h>

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>
#include <Core/Exceptions/TypeMismatchException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/Endian.h>

#include <cstring>
#include <iosfwd>
#include <iostream>
#include <memory>

namespace Uintah {

template<class T>
class PerPatch : public PerPatchBase
{
public:
  // Static variable whose entire purpose is to cause the
  // (instantiated) type of this class to be registered with the
  // Core/Disclosure/TypeDescription class when this class' object
  // code is originally loaded from the shared library.  The
  // 'registerMe' variable is not used for anything else in the
  // program.
  static TypeDescription::Register registerMe;

public:
  inline PerPatch()
    : PerPatchBase()
    , value(std::make_shared<T>())
  {
  }

  inline PerPatch(T value)
    : PerPatchBase()
    , value(std::make_shared<T>(value))
  {
  }

  inline PerPatch(const PerPatch<T>& copy)
    : PerPatchBase()
    , value(copy.value)
  {
  }

  virtual void
  copyPointer(Variable&) override;

  virtual ~PerPatch() = default;

  const TypeDescription*
  virtualGetTypeDescription() const override
  {
    return getTypeDescription();
  }

  static const TypeDescription*
  getTypeDescription();

  inline operator T() const { return *value; }

  inline T&
  get()
  {
    return *value;
  }

  inline const T&
  get() const
  {
    return *value;
  }

  void
  setData(const T& val)
  {
    value = std::make_shared<T>(val);
  }

  virtual PerPatchBase*
  clone() const override
  {
    return scinew PerPatch<T>(*this);
  }

  PerPatch<T>&
  operator=(const PerPatch<T>& copy)
  {
    value = copy.value;
    return *this;
  }

  virtual void
  getSizeInfo(std::string& elems, unsigned long& totsize, void*& ptr) const override
  {
    elems   = "1";
    totsize = getDataSize();
    ptr     = getBasePointer();
  }

  virtual size_t
  getDataSize() const override
  {
    return sizeof(T);
  }

  virtual void*
  getBasePointer() const override
  {
    return value.get();
  }

  virtual bool
  copyOut(void* dst) const override
  {
    void* src       = (void*)(&value);
    size_t numBytes = getDataSize();
    void* retVal    = std::memcpy(dst, src, numBytes);
    return (retVal == dst) ? true : false;
  }

  virtual void
  emitNormal(std::ostream& out,
             [[maybe_unused]] const IntVector& l,
             [[maybe_unused]] const IntVector& h,
             [[maybe_unused]] ProblemSpecP varnode,
             [[maybe_unused]] bool outputDoubleAsFloat) override
  {
    ssize_t linesize = (ssize_t)(sizeof(T));

    out.write((char*)(value.get()), linesize);
  }

  virtual void
  readNormal(std::istream& in, bool swapBytes) override
  {
    ssize_t linesize = (ssize_t)(sizeof(T));

    T val;

    in.read((char*)&val, linesize);

    if (swapBytes) {
      Uintah::swapbytes(val);
    }

    value = std::make_shared<T>(val);
  }

  void
  print(std::ostream& out) const
  {
    out << *(value.get());
  }

private:
  inline static TypeDescription* td{ nullptr };

  static Variable*
  maker()
  {
    return scinew PerPatch<T>();
  }

private:
  std::shared_ptr<T> value;
}; // end class PerPatch

// The following line is the initialization (creation) of the
// 'registerMe' static variable (for each version of CCVariable
// (double, int, etc)).  Note, the 'registerMe' variable is created
// when the object code is initially loaded (usually during intial
// program load by the operating system).
template<class T>
TypeDescription::Register PerPatch<T>::registerMe(getTypeDescription());

template<class T>
const TypeDescription*
PerPatch<T>::getTypeDescription()
{
  if (!td) {
    // this is a hack to get a non-null perpatch
    // var for some functions the perpatches are used in (i.e., task->computes).
    // Since they're not fully-qualified variables, maker
    // would fail anyway.  And since most instances use Handle, it would be
    // difficult.
    td = scinew TypeDescription(TypeDescription::Type::PerPatch,
                                "PerPatch",
                                &maker,
                                fun_getTypeDescription((int*)nullptr));
  }
  return td;
}

// Manually list the int and double basic data type.  If others are
// needed, list them here too.  For everything else, a hacky
// solution is used which defaults their internal type as if it were
// int
template<>
inline const TypeDescription*
PerPatch<int>::getTypeDescription()
{
  if (!td) {
    TypeDescription* sub_td;
    sub_td = scinew TypeDescription(
      TypeDescription::Type::int_type, "int", true, MPI_INT);

    td = scinew TypeDescription(
      TypeDescription::Type::PerPatch, "PerPatch", &maker, sub_td);
  }
  return td;
}

template<>
inline const TypeDescription*
PerPatch<double>::getTypeDescription()
{
  if (!td) {
    TypeDescription* sub_td;
    sub_td = scinew TypeDescription(
      TypeDescription::Type::double_type, "double", true, MPI_DOUBLE);

    td = scinew TypeDescription(
      TypeDescription::Type::PerPatch, "PerPatch", &maker, sub_td);
  }
  return td;
}

template<class T>
void
PerPatch<T>::copyPointer(Variable& copy)
{
  const PerPatch<T>* c = dynamic_cast<const PerPatch<T>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in PerPatch variable", __FILE__, __LINE__));
  }
  *this = *c;
}

template<>
inline void
PerPatch<double>::copyPointer(Variable& copy)
{
  const PerPatch<double>* c = dynamic_cast<const PerPatch<double>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in PerPatch variable", __FILE__, __LINE__));
  }
  *this = *c;
}

template<>
inline void
PerPatch<double*>::copyPointer(Variable& copy)
{
  const PerPatch<double*>* c = dynamic_cast<const PerPatch<double*>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in PerPatch variable", __FILE__, __LINE__));
  }
  *this = *c;
}

} // End namespace Uintah

#endif // __CORE_GRID_VARIABLES_PERPATCH_H__
