/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/Contact/ContactFactory.h>

#include <CCA/Components/MPM/Contact/ApproachContact.h>
#include <CCA/Components/MPM/Contact/CompositeContact.h>
#include <CCA/Components/MPM/Contact/FluidContact.h>
#include <CCA/Components/MPM/Contact/FrictionContact.h>
#include <CCA/Components/MPM/Contact/FrictionContactBard.h>
#include <CCA/Components/MPM/Contact/FrictionContactLR.h>
#include <CCA/Components/MPM/Contact/FrictionContactLRGuilkey.h>
#include <CCA/Components/MPM/Contact/NodalSVFContact.h>
#include <CCA/Components/MPM/Contact/NullContact.h>
#include <CCA/Components/MPM/Contact/SingleVelContact.h>
#include <CCA/Components/MPM/Contact/SpecifiedBodyContact.h>
#include <CCA/Components/MPM/Contact/SpecifiedBodyFrictionContact.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/HydroMPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>

using namespace Uintah;

std::unique_ptr<Contact>
ContactFactory::create(const ProcessorGroup* myworld,
                       const ProblemSpecP& ps,
                       const MaterialManagerP& mat_manager,
                       const MPMLabel* labels,
                       const MPMFlags* flags)
{

  ProblemSpecP mpm_ps =
    ps->findBlockWithOutAttribute("MaterialProperties")->findBlock("MPM");

  if (!mpm_ps) {
    string warn = "ERROR: Missing either <MaterialProperties> or <MPM> block "
                  "from input file";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  ProblemSpecP child = mpm_ps->findBlock("contact");

  std::unique_ptr<CompositeContact> contact_list =
    std::make_unique<CompositeContact>(myworld, labels, flags, child);

  for (; child != 0; child = child->findNextBlock("contact")) {

    std::string con_type;
    child->getWithDefault("type", con_type, "null");

    if (con_type == "null") {
      contact_list->add(std::make_unique<NullContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "single_velocity") {
      contact_list->add(std::make_unique<SingleVelContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "nodal_svf") {
      contact_list->add(std::make_unique<NodalSVFContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "fluid") {
      contact_list->add(std::make_unique<FluidContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "friction") {
      contact_list->add(std::make_unique<FrictionContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "friction_bard") {
      contact_list->add(std::make_unique<FrictionContactBard>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "friction_LR") {
      contact_list->add(std::make_unique<FrictionContactLR>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "friction_LR_guilkey") {
      contact_list->add(std::make_unique<FrictionContactLRGuilkey>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "approach") {
      contact_list->add(std::make_unique<ApproachContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "specified_velocity" || con_type == "specified" ||
             con_type == "rigid") {
      contact_list->add(std::make_unique<SpecifiedBodyContact>(
        myworld, mat_manager, labels, flags, child));
    }

    else if (con_type == "specified_friction") {
      contact_list->add(std::make_unique<SpecifiedBodyFrictionContact>(
        myworld, mat_manager, labels, flags, child));

    }

    else {
      std::cerr << "Unknown Contact Type R (" << con_type << ")\n";
      throw ProblemSetupException(
        " E R R O R----->MPM:Unknown Contact type", __FILE__, __LINE__);
    }
  }

  if (contact_list->size() == 0) {
    std::cout << "no contact - using null\n";
    contact_list->add(std::make_unique<NullContact>(
      myworld, mat_manager, labels, flags, child));
  }

  return contact_list;
}
