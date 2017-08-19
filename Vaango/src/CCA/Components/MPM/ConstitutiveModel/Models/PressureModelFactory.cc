/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/PressureModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Air.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Borja.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Granite.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Hyperelastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Hypoelastic.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_MieGruneisen.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Water.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>
using std::cerr;
using std::ifstream;
using std::ofstream;

using namespace Vaango;
using namespace Uintah;

PressureModel*
PressureModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("pressure_model");
  if (!child) {
    throw ProblemSetupException("Cannot find pressure_model tag.", __FILE__,
                                __LINE__);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for pressure_model", __FILE__,
                                __LINE__);

  if (mat_type == "mie_gruneisen")
    return (scinew Pressure_MieGruneisen(child));
  else if (mat_type == "default_hypo")
    return (scinew Pressure_Hypoelastic(child));
  else if (mat_type == "default_hyper")
    return (scinew Pressure_Hyperelastic(child));
  else if (mat_type == "borja_pressure")
    return (scinew Pressure_Borja(child));
  else if (mat_type == "air")
    return (scinew Pressure_Air(child));
  else if (mat_type == "water")
    return (scinew Pressure_Water(child));
  else if (mat_type == "granite")
    return (scinew Pressure_Granite(child));
  else {
    throw ProblemSetupException("Cannot create pressure_model.", __FILE__,
                                __LINE__);
  }
}

PressureModel*
PressureModelFactory::createCopy(const PressureModel* eos)
{
  if (dynamic_cast<const Pressure_Borja*>(eos))
    return (scinew Pressure_Borja(dynamic_cast<const Pressure_Borja*>(eos)));

  else if (dynamic_cast<const Pressure_MieGruneisen*>(eos))
    return (scinew Pressure_MieGruneisen(
      dynamic_cast<const Pressure_MieGruneisen*>(eos)));

  else if (dynamic_cast<const Pressure_Hypoelastic*>(eos))
    return (scinew Pressure_Hypoelastic(
      dynamic_cast<const Pressure_Hypoelastic*>(eos)));

  else if (dynamic_cast<const Pressure_Air*>(eos))
    return (scinew Pressure_Air(dynamic_cast<const Pressure_Air*>(eos)));

  else if (dynamic_cast<const Pressure_Water*>(eos))
    return (scinew Pressure_Water(dynamic_cast<const Pressure_Water*>(eos)));

  else if (dynamic_cast<const Pressure_Granite*>(eos))
    return (
      scinew Pressure_Granite(dynamic_cast<const Pressure_Granite*>(eos)));

  else {
    throw ProblemSetupException("Cannot create copy of pressure_model.",
                                __FILE__, __LINE__);
  }
}
