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
    double mu;
    contact_ps->get("friction_coeff", mu);
    return std::make_shared<MPMFrictionContact>(dwis, patch, mu);
  } else if (contact_algo == "velocity") {
    return std::make_shared<MPMVelocityContact>(dwis, patch);
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown contact algorithm" << contact_algo << " for MPM. " << std::endl;
    throw Exception(out.str(), __FILE__, __LINE__);
  }
}

} /* namespace BrMPM */
