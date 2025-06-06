#ifndef __CORE_GRID_VARIABLES_REDUCTION_VARIABLE_BASE_H__
#define __CORE_GRID_VARIABLES_REDUCTION_VARIABLE_BASE_H__

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

#include <Core/Grid/Variables/Variable.h>
#include <Core/Parallel/UintahMPI.h>

namespace Uintah {

class ReductionVariableBase : public Variable
{
public:
  virtual ~ReductionVariableBase();

  virtual const TypeDescription* virtualGetTypeDescription() const override = 0;

  virtual void copyPointer(Variable&) override = 0;
  virtual ReductionVariableBase* clone() const = 0;
  virtual RefCounted* getRefCounted() override;

  virtual void getSizeInfo(std::string& elems, unsigned long& totsize,
                           void*& ptr) const override = 0;
  virtual size_t getDataSize() const override = 0;
  virtual bool copyOut(void* dst) const override = 0;
  virtual void* getBasePointer() const = 0;

  virtual void emitNormal(std::ostream& out, const IntVector& l,
                          const IntVector& h, ProblemSpecP varnode,
                          bool outputDoubleAsFloat) override = 0;
  virtual void readNormal(std::istream& in, bool swapbytes) override = 0;
  virtual void allocate(const Patch* patch, const IntVector& boundary) override;

  virtual void print(std::ostream&) const = 0;

  virtual void reduce(const ReductionVariableBase&) = 0;
  virtual void getMPIInfo(int& count, MPI_Datatype& datatype, MPI_Op& op) = 0;
  virtual void getMPIData(std::vector<char>& buf, int& index) = 0;
  virtual void putMPIData(std::vector<char>& buf, int& index) = 0;
  virtual void setBenignValue() = 0;

protected:
  ReductionVariableBase(const ReductionVariableBase&);
  ReductionVariableBase();

  ReductionVariableBase& operator=(const ReductionVariableBase&) = delete;
};
} // End namespace Uintah

#endif //__CORE_GRID_VARIABLES_REDUCTION_VARIABLE_BASE_H__
