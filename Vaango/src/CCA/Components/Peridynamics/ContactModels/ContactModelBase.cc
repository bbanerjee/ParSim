#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>

using namespace Vaango;

ContactModelBase::ContactModelBase(const Uintah::ProcessorGroup* myworld, 
                                   PeridynamicsLabel* labels, 
                                   PeridynamicsFlags* flags, 
                                   Uintah::ProblemSpecP ps)
  : Uintah::UintahParallelComponent(myworld), d_labels(labels), d_flags(flags), d_bodiesThatCanInteract(ps)
{
}

ContactModelBase::~ContactModelBase()
{
}
