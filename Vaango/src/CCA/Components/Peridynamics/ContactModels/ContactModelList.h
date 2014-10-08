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

#ifndef __VAANGO_CONTACT_MODEL_LIST_H__
#define __VAANGO_CONTACT_MODEL_LIST_H__

#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <list>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class   ContactModelList
    \brief   Derived from ContactModelBase but acts just as a container for 
             a list of contact models
    \author  Biswajit Banerjee 
    \original Andrew Brydon
    \warning 
  */
  /////////////////////////////////////////////////////////////////////////////

  class ContactModelList : public ContactModelBase {

  public:

    // Constructor
    ContactModelList(const Uintah::ProcessorGroup* myworld, 
                     PeridynamicsLabel* labels, 
                     PeridynamicsFlags* flags);

    virtual ~ContactModelList();

    void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    void add(ContactModelBase* model);
         
    size_t size() const { return d_modelList.size(); }
         
    void exchangeMomentumInterpolated(const Uintah::ProcessorGroup*,
                                      const Uintah::PatchSubset* patches,
                                      const Uintah::MaterialSubset* matls,
                                      Uintah::DataWarehouse* old_dw,
                                      Uintah::DataWarehouse* new_dw);
         
    void exchangeMomentumIntegrated(const Uintah::ProcessorGroup*,
                                    const Uintah::PatchSubset* patches,
                                    const Uintah::MaterialSubset* matls,
                                    Uintah::DataWarehouse* old_dw,
                                    Uintah::DataWarehouse* new_dw);
         
    void addComputesAndRequiresInterpolated(Uintah::SchedulerP & sched,
                                            const Uintah::PatchSet* patches,
                                            const Uintah::MaterialSet* matls);
         
    void addComputesAndRequiresIntegrated(Uintah::SchedulerP & sched,
                                          const Uintah::PatchSet* patches,
                                          const Uintah::MaterialSet* matls);

  private: 

    std::list<ContactModelBase*> d_modelList;

    // Allow no copies of the list
    ContactModelList(const ContactModelList &);
    ContactModelList& operator=(const ContactModelList &);

  };
      
} // End namespace Vaango

#endif // __VAANGO_CONTACT_MODEL_LIST_H__
