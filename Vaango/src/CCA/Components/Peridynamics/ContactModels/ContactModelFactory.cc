#include <CCA/Components/Peridynamics/ContactModels/ContactModelFactory.h>
#include <CCA/Components/Peridynamics/ContactModels/ContactModelBase.h>
#include <CCA/Components/Peridynamics/ContactModels/NullContact.h>
#include <CCA/Components/Peridynamics/ContactModels/SingleVelocityContact.h>
#include <CCA/Components/Peridynamics/ContactModels/ContactModelList.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <string>

using namespace Vaango;

ContactModelBase* 
ContactModelFactory::create(const Uintah::ProcessorGroup* myworld,
                            const Uintah::ProblemSpecP& ps, 
                            Uintah::SimulationStateP &ss,
                            PeridynamicsLabel* labels, 
                            PeridynamicsFlags* flags)
{
  // Get the peridynamics material probem spec
  Uintah::ProblemSpecP peri_mat_ps = 
    ps->findBlockWithOutAttribute("MaterialProperties")->findBlock("Peridynamics");
  if (!peri_mat_ps) {
    std::ostringstream out;
    out << "Cannot create contact model.  The input file does not contain a" 
        << " <MaterialProperties> block that contains a <Peridynamics> block.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  ContactModelList* contact_list = scinew ContactModelList(myworld, labels, flags);
   
  for (Uintah::ProblemSpecP contact_ps = peri_mat_ps->findBlock("ContactModel"); contact_ps != 0;
       contact_ps = contact_ps->findNextBlock("ContactModel")) {
     
    std::string contact_type;
    contact_ps->getWithDefault("type", contact_type, "null");
     
    if (contact_type == "null")
      contact_list->add(scinew NullContact(myworld, ss, labels, flags));
      
    else if (contact_type == "single_velocity")
      contact_list->add(scinew SingleVelocityContact(myworld, contact_ps, ss, labels, flags));
      
    else {
      std::ostringstream out;
      out << "**ERROR** Unknown Contact Type (" << contact_type << ")" ;
      throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
      }
   }
   
   if (contact_list->size()==0) {
     std::cout << "Peridynamics: No contact model provided in input file - using null" << std::endl;
     contact_list->add(scinew NullContact(myworld, ss, labels, flags));
   }
   
   return contact_list;
}
