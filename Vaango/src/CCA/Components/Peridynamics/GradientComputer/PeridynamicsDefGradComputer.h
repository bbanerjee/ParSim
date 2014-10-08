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

#ifndef __VAANGO_PERIDYNAMICS_DEFORMATION_GRADIENT_COMPUTER_H__
#define __VAANGO_PERIDYNAMICS_DEFORMATION_GRADIENT_COMPUTER_H__

#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Task.h>

#include <vector>

namespace Uintah {
  class Patch;
  class DataWarehouse;
  class VarLabel;
}

namespace Vaango {

  class PeridynamicsFlags;
  class PeridynamicsMaterial;
  class PeridynamicsLabel;

  class PeridynamicsDefGradComputer {

  public:

    PeridynamicsDefGradComputer(PeridynamicsFlags* flags, 
                                PeridynamicsLabel* labels);
    virtual ~PeridynamicsDefGradComputer();

    void addInitialComputesAndRequires(Uintah::Task* task,
                                       const PeridynamicsMaterial* matl,
                                       const Uintah::PatchSet*);

    void addComputesAndRequires(Uintah::Task* task,
                                const PeridynamicsMaterial* matl,
                                const Uintah::PatchSet*);

    void initialize(const Uintah::Patch* patch,    
                    PeridynamicsMaterial* matl, 
                    Uintah::DataWarehouse* new_dw);

    void computeDeformationGradient(const Uintah::Patch* patch,
                                    const PeridynamicsMaterial* matl,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);

   void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                         std::vector<const Uintah::VarLabel*>& to);

  private:

    PeridynamicsLabel* d_labels;
    PeridynamicsFlags* d_flags;
  };

} // end namespace Vaango

#endif
