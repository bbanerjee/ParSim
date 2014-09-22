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
