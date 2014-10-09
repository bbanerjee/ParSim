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

#include <MaterialModels/DamageModelFactory.h>
#include <MaterialModels/DamageModelSimple.h>
#include <Pointers/DamageModelSP.h>

#include <Core/Exception.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <iostream>

using namespace Matiti;

DamageModelSP 
DamageModelFactory::create(Uintah::ProblemSpecP& damage_ps)
{
  // Get the attribute of the type
  std::string damage_type;
  if (!damage_ps->getAttribute("type", damage_type)) {
    std::ostringstream out;
    out << "**ERROR** Damage model does not have type attribute.";
    out << " Cannot determine which damage model to apply." ;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }

  // Create geometry
  if (damage_type == "isotropic") {
    return std::make_shared<DamageModelSimple>();
  } else {
    std::ostringstream out;
    out << "**ERROR** Unknown damage model type" << damage_type;
    throw Exception(out.str(), __FILE__, __LINE__); 
  }
}
