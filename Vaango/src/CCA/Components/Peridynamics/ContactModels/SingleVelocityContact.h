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

#ifndef __VAANGO_SINGLE_VELOCITY_CONTACT_H__
#define __VAANGO_SINGLE_VELOCITY_CONTACT_H__

#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   SingleVelocityContact
    \brief   Recaptures single velocity field behavior from objects belonging to 
             multiple velocity fields.  Ensures that one can get the same answer 
             using prescribed contact as can be gotten using "automatic" 
             MPM-type contact on the background grid.
    \author  Steven G. Parker
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class SingleVelocityContact : public ContactModelBase {

  public:

    SingleVelocityContact(const Uintah::ProcessorGroup* myworld,
                          Uintah::ProblemSpecP& ps,
                          Uintah::MaterialManagerP& mat_manager,
                          PeridynamicsLabel* labels,
                          PeridynamicsFlags* flags);
         
    virtual ~SingleVelocityContact();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);

    virtual void exchangeMomentumInterpolated(const Uintah::ProcessorGroup*,
                                              const Uintah::PatchSubset* patches,
                                              const Uintah::MaterialSubset* matls,
                                              Uintah::DataWarehouse* old_dw,
                                              Uintah::DataWarehouse* new_dw);
         
    virtual void exchangeMomentumIntegrated(const Uintah::ProcessorGroup*,
                                            const Uintah::PatchSubset* patches,
                                            const Uintah::MaterialSubset* matls,
                                            Uintah::DataWarehouse* old_dw,
                                            Uintah::DataWarehouse* new_dw);

    virtual void addComputesAndRequiresInterpolated(Uintah::SchedulerP & sched,
                                                    const Uintah::PatchSet* patches,
                                                    const Uintah::MaterialSet* matls);

    virtual void addComputesAndRequiresIntegrated(Uintah::SchedulerP & sched,
                                                  const Uintah::PatchSet* patches,
                                                  const Uintah::MaterialSet* matls);
  private:
         
    // Prevent copying of this class
    SingleVelocityContact(const SingleVelocityContact &con) = delete;
    SingleVelocityContact& operator=(const SingleVelocityContact &con) = delete;

  };
} // End namespace Vaango
      

#endif /* __VAANGO_SINGLE_VELOCITY_CONTACT_H__ */

