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

#include "ViscoPlasticityModelFactory.h"
#include "SuvicI.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>
using std::cerr;
using std::ifstream;
using std::ofstream;

using namespace Uintah;

ViscoPlasticityModel*
ViscoPlasticityModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("viscoplastic_flow_model");
  if (!child)
    throw ProblemSetupException("Cannot find viscoplastic_flow_model tag",
                                __FILE__, __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for viscoplastic_flow_model", __FILE__,
                                __LINE__);

  if (mat_type == "suvic_i")
    return (scinew SuvicI(child));
  else
    throw ProblemSetupException(
      "Unknown ViscoPlasticity Model (" + mat_type + ")", __FILE__, __LINE__);
}

ViscoPlasticityModel*
ViscoPlasticityModelFactory::createCopy(const ViscoPlasticityModel* pm)
{

  if (dynamic_cast<const SuvicI*>(pm))
    return (scinew SuvicI(dynamic_cast<const SuvicI*>(pm)));

  else
    throw ProblemSetupException(
      "Cannot create copy of unknown Viscoplasticity model", __FILE__,
      __LINE__);
}
