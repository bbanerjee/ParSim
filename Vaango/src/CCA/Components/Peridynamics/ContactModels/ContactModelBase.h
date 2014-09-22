#ifndef __VAANGO_CONTACT_MODEL_BASE_H__
#define __VAANGO_CONTACT_MODEL_BASE_H__

#include <Core/Parallel/UintahParallelComponent.h>
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

  class ContactModelBase : public Uintah::UintahParallelComponent {

  public:

    // Constructor
    ContactModelBase(const Uintah::ProcessorGroup* myworld, 
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
