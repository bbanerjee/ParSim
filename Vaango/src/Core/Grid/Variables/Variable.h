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

#ifndef __CORE_GRID_VARIABLES_Variable_H__
#define __CORE_GRID_VARIABLES_Variable_H__

#include <Core/ProblemSpec/ProblemSpec.h>

#include <sci_defs/pidx_defs.h>

#include <iosfwd>
#include <string>

namespace Uintah {

class TypeDescription;
class InputContext;
class OutputContext;
class PIDXOutputContext;
class Patch;
class RefCounted;
class VarLabel;

class Variable
{

public:
  virtual ~Variable();

  // eliminate copy, assignment and move
  Variable(const Variable&) = delete;
  Variable(Variable&&)      = delete;

  Variable&
  operator=(const Variable&) = delete;
  Variable&
  operator=(Variable&&) = delete;

  virtual const TypeDescription*
  virtualGetTypeDescription() const = 0;

  void
  setForeign();

  bool
  isForeign() const
  {
    return d_foreign;
  }

  // marks a variable as invalid (for example, it is in the process of receiving
  // mpi)
  void
  setValid()
  {
    d_valid = true;
  }

  void
  setInvalid()
  {
    d_valid = false;
  }

  // returns if a variable is marked valid or invalid
  bool
  isValid() const
  {
    return d_valid;
  }

  size_t
  emit(OutputContext&,
       const IntVector& l,
       const IntVector& h,
       const std::string& compressionModeHint);

  void
  read(InputContext&,
       long end,
       bool swapbytes,
       int nByteMode,
       const std::string& compressionMode);

#if HAVE_PIDX
  virtual void
  emitPIDX(PIDXOutputContext& oc,
           unsigned char* buffer,
           const IntVector& l,
           const IntVector& h,
           const size_t pidx_bufferSize // buffer size used for bullet proofing.
  );

  void
  readPIDX(const unsigned char* pidx_buffer,
           const size_t& pidx_bufferSize,
           const bool swapBytes);
#endif

  virtual void
  emitNormal(std::ostream& out,
             const IntVector& l,
             const IntVector& h,
             ProblemSpecP varnode,
             bool outputDoubleAsFloat) = 0;

  virtual void
  readNormal(std::istream& in, bool swapbytes) = 0;

  virtual void
  allocate(const Patch* patch, const IntVector& boundary) = 0;

  virtual void
  getSizeInfo(std::string& elems, unsigned long& totsize, void*& ptr) const = 0;

  // used to get size info of the underlying data; this is for host-->device
  // variable copy
  virtual size_t
  getDataSize() const = 0;

  // used to copy Variables to contiguous buffer prior to bulk host-->device
  // copy
  virtual bool
  copyOut(void* dst) const = 0;

  virtual void
  copyPointer(Variable&) = 0;

  // Only affects grid variables
  virtual void
  offsetGrid(const IntVector& /*offset*/);

  virtual RefCounted*
  getRefCounted() = 0;

protected:
  Variable();

private:
  // Compresses the string pointed to by pUncompressed and but the
  // resulting compressed data into the string pointed to by pBuffer.
  // Returns the pointer to whichever one is shortest and erases the
  // other one.
  std::string*
  gzipCompress(std::string* pUncompressed, std::string* pBuffer);

  // states that the variable is from another node - these variables (ghost
  // cells, slabs, corners) are communicated via MPI
  bool d_foreign{ false };

  // signals of the variable is valid, an mpi variable is not valid until mpi
  // has been recieved
  bool d_valid{ true };
};

} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_Variable_H__
