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

#ifndef __VAANGO_CONTACT_MODEL_BASE_H__
#define __VAANGO_CONTACT_MODEL_BASE_H__

#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <CCA/Ports/SchedulerP.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <cmath>

namespace Uintah {
  class DataWarehouse;
  class ProcessorGroup;
  class Patch;
  class VarLabel;
  class Task;
}

namespace Vaango {

  class PeridynamicsLabel;
  class PeridynamicsFlags;

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   ContactModelBase
    \brief   Base class for contact models for Peridynamics-MPM adapted from MPM.
    \author  Biswajit Banerjee 
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class ContactModelBase {

  public:

    // Constructor
    ContactModelBase(const Uintah::ProcessorGroup* myworld, 
                     const Uintah::MaterialManagerP& mat_manager,
                     PeridynamicsLabel* labels, 
                     PeridynamicsFlags* flags,
                     Uintah::ProblemSpecP ps);

    virtual ~ContactModelBase();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;

    // Basic contact methods
    virtual void exchangeMomentumInterpolated(const Uintah::ProcessorGroup*,
                                              const Uintah::PatchSubset* patches,
                                              const Uintah::MaterialSubset* matls,
                                              Uintah::DataWarehouse* old_dw,
                                              Uintah::DataWarehouse* new_dw) = 0;
         
    virtual void exchangeMomentumIntegrated(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* matls,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw) = 0;
         
    virtual void addComputesAndRequiresInterpolated(Uintah::SchedulerP & sched,
                                                    const Uintah::PatchSet* patches,
                                                    const Uintah::MaterialSet* matls) = 0;
         
    virtual void addComputesAndRequiresIntegrated(Uintah::SchedulerP & sched,
                                                  const Uintah::PatchSet* patches,
                                                  const Uintah::MaterialSet* matls) = 0;
         
  protected:

    const Uintah::MaterialManagerP& d_mat_manager;
    PeridynamicsLabel* d_labels;
    PeridynamicsFlags* d_flags;
         
    Uintah::ContactMaterialSpec d_bodiesThatCanInteract;

  };
      
  inline bool compare(double num1, double num2) {
    double EPSILON=1.e-14;
    return (std::abs(num1-num2) <= EPSILON);
  }

} // End namespace Vaango

#endif // __VAANGO_CONTACT_MODEL_BASE_H__
