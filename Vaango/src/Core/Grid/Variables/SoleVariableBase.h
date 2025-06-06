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

#ifndef UINTAH_HOMEBREW_SoleVariableBase_H
#define UINTAH_HOMEBREW_SoleVariableBase_H

#include <Core/Grid/Variables/Variable.h>
#include <iosfwd>
#include <sci_defs/mpi_defs.h> // For MPIPP_H on SGI

namespace Uintah {

/**************************************

CLASS
   SoleVariableBase

   Short description...

GENERAL INFORMATION

   SoleVariableBase.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   SoleVariableBase

DESCRIPTION
   Long description...

WARNING

****************************************/

class SoleVariableBase : public Variable
{

public:
  virtual ~SoleVariableBase();

  virtual void
  copyPointer(Variable&) override = 0;

  virtual SoleVariableBase*
  clone() const = 0;

  virtual const TypeDescription*
  virtualGetTypeDescription() const override = 0;

  virtual RefCounted*
  getRefCounted() override;

  virtual void
  getSizeInfo(std::string& elems,
              unsigned long& totsize,
              void*& ptr) const override = 0;

  virtual size_t
  getDataSize() const override = 0;

  virtual bool
  copyOut(void* dst) const override = 0;

  virtual void*
  getBasePointer() const = 0;

  virtual void
  emitNormal(std::ostream& out,
             const IntVector& l,
             const IntVector& h,
             ProblemSpecP varnode,
             bool outputDoubleAsFloat) override = 0;

  virtual void
  readNormal(std::istream& in, bool swapbytes) override = 0;

  virtual void
  allocate(const Patch* patch, const IntVector& boundary) override;

  virtual void
  print(std::ostream&) const = 0;

protected:
  SoleVariableBase(const SoleVariableBase&) = default;
  SoleVariableBase();

private:
  SoleVariableBase&
  operator=(const SoleVariableBase&) = delete;
};
} // End namespace Uintah

#endif
