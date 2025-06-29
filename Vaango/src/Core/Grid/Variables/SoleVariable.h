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

#ifndef UINTAH_HOMEBREW_SoleVARIABLE_H
#define UINTAH_HOMEBREW_SoleVARIABLE_H

#include <Core/Disclosure/TypeDescription.h>
#include <Core/Disclosure/TypeUtils.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/TypeMismatchException.h>
#include <Core/Grid/Variables/DataItem.h>
#include <Core/Grid/Variables/SoleVariableBase.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Util/Endian.h>

#include <cstring>
#include <iosfwd>
#include <iostream>

namespace Uintah {

/**************************************

CLASS
   SoleVariable

   Short description...

GENERAL INFORMATION

   SoleVariable.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   Sole_Variable

DESCRIPTION
   Long description...

WARNING

****************************************/

template<class T>
class SoleVariable : public SoleVariableBase
{
public:
  inline SoleVariable()
    : value(std::make_shared<T>())
  {
  }

  inline SoleVariable(T value)
    : value(std::make_shared<T>(value))
  {
  }

  inline SoleVariable(const SoleVariable<T>& copy)
    : value(copy.value)
  {
  }

  virtual ~SoleVariable() = default;

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

  virtual SoleVariableBase*
  clone() const override
  {
    return scinew SoleVariable<T>(*this);
  }

  SoleVariable<T>&
  operator=(const SoleVariable<T>& copy)
  {
    value = copy.value;
    return *this;
  };

  virtual void
  copyPointer(Variable&) override;

  virtual void
  getSizeInfo(std::string& elems,
              unsigned long& totsize,
              void*& ptr) const override
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
  getBasePointer() const
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
             ProblemSpecP /*varnode*/,
             [[maybe_unused]] bool outputDoubleAsFloat)
  {
    ssize_t linesize = (ssize_t)(sizeof(T));

    out.write((char*)(value.get()), linesize);
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

    value = std::make_shared<T>(val);
  }

  void
  print(std::ostream& out) const
  {
    out << *(value.get());
  }

  // Static variable whose entire purpose is to cause the
  // (instantiated) type of this class to be registered with the
  // Core/Disclosure/TypeDescription class when this class' object
  // code is originally loaded from the shared library.  The
  // 'registerMe' variable is not used for anything else in the
  // program.
  static TypeDescription::Register registerMe;

private:
  inline static TypeDescription* td{ nullptr };
  static Variable*
  maker()
  {
    return scinew SoleVariable<T>();
  };

  std::shared_ptr<T> value;
};

// The following line is the initialization (creation) of the
// 'registerMe' static variable (for each version of CCVariable
// (double, int, etc)).  Note, the 'registerMe' variable is created
// when the object code is initially loaded (usually during intial
// program load by the operating system).
template<class T>
TypeDescription::Register SoleVariable<T>::registerMe(getTypeDescription());

template<class T>
const TypeDescription*
SoleVariable<T>::getTypeDescription()
{
  static TypeDescription* td;
  if (!td) {
    td = scinew TypeDescription(TypeDescription::Type::SoleVariable,
                                "SoleVariable",
                                &maker,
                                fun_getTypeDescription((int*)0));
  }
  return td;
}

template<class T>
void
SoleVariable<T>::copyPointer(Variable& copy)
{
  SoleVariable<T>* c = dynamic_cast<SoleVariable<T>*>(&copy);
  if (!c) {
    SCI_THROW(TypeMismatchException(
      "Type mismatch in sole variable", __FILE__, __LINE__));
  }
  *this = *c;
}
} // End namespace Uintah

#endif
