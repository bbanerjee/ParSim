/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#include "MPMEquationOfStateFactory.h"
#include "AirEOS.h"
#include "BorjaEOS.h"
#include "DefaultHypoElasticEOS.h"
#include "GraniteEOS.h"
#include "HyperElasticEOS.h"
#include "MieGruneisenEOS.h"
#include "MieGruneisenEOSEnergy.h"
#include "WaterEOS.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace Uintah;
using namespace Vaango;

MPMEquationOfState*
MPMEquationOfStateFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("equation_of_state");
  if (!child) {
    proc0cout
      << "**WARNING** Creating default hyperelastic equation of state\n";
    return scinew HyperElasticEOS(ps);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "No type for equation_of_state", __FILE__, __LINE__);

  if (mat_type == "mie_gruneisen")
    return scinew MieGruneisenEOS(child);
  else if (mat_type == "mie_gruneisen_energy")
    return scinew MieGruneisenEOSEnergy(child);
  else if (mat_type == "default_hypo")
    return scinew DefaultHypoElasticEOS(child);
  else if (mat_type == "default_hyper")
    return scinew HyperElasticEOS(child);
  else if (mat_type == "borja_pressure")
    return scinew BorjaEOS(child);
  else if (mat_type == "air")
    return scinew AirEOS(child);
  else if (mat_type == "water")
    return scinew WaterEOS(child);
  else if (mat_type == "granite")
    return scinew GraniteEOS(child);
  else {
    proc0cout
      << "**WARNING** Creating default hyperelastic equation of state\n";
    return scinew HyperElasticEOS(ps);
  }

  return nullptr;
}

MPMEquationOfState*
MPMEquationOfStateFactory::createCopy(const MPMEquationOfState* eos)
{
  if (dynamic_cast<const MieGruneisenEOS*>(eos))
    return scinew MieGruneisenEOS(dynamic_cast<const MieGruneisenEOS*>(eos));

  else if (dynamic_cast<const MieGruneisenEOSEnergy*>(eos))
    return scinew MieGruneisenEOSEnergy(
      dynamic_cast<const MieGruneisenEOSEnergy*>(eos));

  else if (dynamic_cast<const DefaultHypoElasticEOS*>(eos))
    return scinew DefaultHypoElasticEOS(
      dynamic_cast<const DefaultHypoElasticEOS*>(eos));

  else if (dynamic_cast<const HyperElasticEOS*>(eos))
    return scinew HyperElasticEOS(dynamic_cast<const HyperElasticEOS*>(eos));

  else if (dynamic_cast<const BorjaEOS*>(eos))
    return scinew BorjaEOS(dynamic_cast<const BorjaEOS*>(eos));

  else if (dynamic_cast<const AirEOS*>(eos))
    return scinew AirEOS(dynamic_cast<const AirEOS*>(eos));

  else if (dynamic_cast<const WaterEOS*>(eos))
    return scinew WaterEOS(dynamic_cast<const WaterEOS*>(eos));

  else if (dynamic_cast<const GraniteEOS*>(eos))
    return scinew GraniteEOS(dynamic_cast<const GraniteEOS*>(eos));

  else {
    proc0cout << "**WARNING** Creating a copy of the default hyperelastic "
                 "equation of state\n";
    return scinew HyperElasticEOS(dynamic_cast<const HyperElasticEOS*>(eos));
  }

  return nullptr;
}
