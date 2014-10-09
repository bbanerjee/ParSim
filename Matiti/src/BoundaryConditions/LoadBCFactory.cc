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

#include <BoundaryConditions/LoadBCFactory.h>
#include <BoundaryConditions/LoadBC.h>
#include <BoundaryConditions/ForceBC.h>
#include <BoundaryConditions/TractionBC.h>
#include <Pointers/LoadBCSP.h>

#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

LoadBCSP 
LoadBCFactory::create(Uintah::ProblemSpecP& load_ps)
{
  // Get the attribute of the type
  std::string load_type;
  if (!load_ps->getAttribute("type", load_type)) {
    std::ostringstream out;
    out << "**ERROR** Load BC does not have type attribute to determine if it is a force/traction/pressure etc." 
        << load_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (load_type == "force") {
    return (std::make_shared<ForceBC>());
  } else if (load_type == "traction") {
    return (std::make_shared<TractionBC>());
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown load boundary condition type" << load_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
