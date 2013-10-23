/*
 * MPMContactFactory.cc
 *
 *  Created on: 23/10/2013
 *      Author: banerjee
 */

#include <Contact/MPMContactFactory.h>
#include <Exception.h>

#include <Contact/MPMContact.h>
#include <Contact/MPMFreeContact.h>
#include <Contact/MPMFrictionContact.h>
#include <Contact/MPMFrictionlessContact.h>
#include <Contact/MPMVelocityContact.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>
#include <string>

using namespace BrMPM;
namespace BrMPM
{

MPMContactFactory::MPMContactFactory()
{
}

MPMContactFactory::~MPMContactFactory()
{
}

MPMContactP
MPMContactFactory::create(const Uintah::ProblemSpecP& ps,
                          std::vector<int>& dwis,
                          MPMPatchP& patch)
{
  Uintah::ProblemSpecP contact_ps = ps->findBlock("Contact");
  if (!contact_ps) {
    throw Exception("**ERROR** <Contact> tag not found", __FILE__, __LINE__);
  }

  std::string contact_algo;
  if (!contact_ps->getAttribute("algorithm", contact_algo)) {
    throw Exception("**ERROR** <Contact algorithm=?> Algorithm attribute not found",
                    __FILE__, __LINE__);
  }

  // Create contact function
  if (contact_algo == "free") {
    return std::make_shared<MPMFreeContact>(dwis, patch);
  } else if (contact_algo == "frictionless") {
    return std::make_shared<MPMFrictionlessContact>(dwis, patch);
  } else if (contact_algo == "friction") {
    return std::make_shared<MPMFrictionContact>(dwis, patch);
  } else if (contact_algo == "velocity") {
    return std::make_shared<MPMVelocityContact>(dwis, patch);
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown contact algorithm" << contact_algo << " for MPM. " << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}

} /* namespace BrMPM */
