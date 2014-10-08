/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#ifndef __VAANGO_FAMILY_COMPUTER_H__
#define __VAANGO_FAMILY_COMPUTER_H__

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Geometry/Point.h>
#include <Core/Geometry/IntVector.h>

namespace Vaango {

  class PeridynamicsMaterial;
  class PeridynamicsFlags;
  class PeridynamicsLabel;

  class FamilyComputer {

  public:
    
    FamilyComputer(PeridynamicsFlags* flags, PeridynamicsLabel* labels);
    virtual ~FamilyComputer();

    /*! Initial computes and requires for the family computer */
    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet* patches) const;

    void createNeighborList(PeridynamicsMaterial* matl,
                            const Uintah::Patch* patch,
                            Uintah::DataWarehouse* new_dw);

  protected:

    void findCellsInHorizon(const Uintah::Patch* patch,
                            const SCIRun::Point& pos,
                            const double& horizon,
                            SCIRun::IntVector& cellLow,
                            SCIRun::IntVector& cellHigh);

  private:

    PeridynamicsLabel* d_label;
    PeridynamicsFlags* d_flags;
    
  };

} // End of namespace Vaango

#endif // __VAANGO_FAMILY_COMPUTER_H__
