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
                            Uintah::MaterialManagerP &ss,
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
