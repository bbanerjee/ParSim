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

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfStateFactory.h>

#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/AirEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/BorjaEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/DefaultHypoElasticEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/GraniteEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/HyperElasticEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MieGruneisenEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MieGruneisenEOSEnergy.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/WaterEOS.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Uintah;
using namespace Vaango;

std::unique_ptr<MPMEquationOfState>
MPMEquationOfStateFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("equation_of_state");
  if (!child) {
    proc0cout
      << "**WARNING** Creating default hyperelastic equation of state\n";
    return std::make_unique<HyperElasticEOS>(ps);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw ProblemSetupException(
      "No type for equation_of_state", __FILE__, __LINE__);
  }

  if (mat_type == "mie_gruneisen") {
    return std::make_unique<MieGruneisenEOS>(child);
  } else if (mat_type == "mie_gruneisen_energy") {
    return std::make_unique<MieGruneisenEOSEnergy>(child);
  } else if (mat_type == "default_hypo") {
    return std::make_unique<DefaultHypoElasticEOS>(child);
  } else if (mat_type == "default_hyper") {
    return std::make_unique<HyperElasticEOS>(child);
  } else if (mat_type == "borja_pressure") {
    return std::make_unique<BorjaEOS>(child);
  } else if (mat_type == "air") {
    return std::make_unique<AirEOS>(child);
  } else if (mat_type == "water") {
    return std::make_unique<WaterEOS>(child);
  } else if (mat_type == "granite") {
    return std::make_unique<GraniteEOS>(child);
  } else {
    proc0cout
      << "**WARNING** Creating default hyperelastic equation of state\n";
    return std::make_unique<HyperElasticEOS>(ps);
  }

  return nullptr;
}

std::unique_ptr<MPMEquationOfState>
MPMEquationOfStateFactory::createCopy(const MPMEquationOfState* eos)
{
  if (dynamic_cast<const MieGruneisenEOS*>(eos)) {
    return std::make_unique<MieGruneisenEOS>(
      dynamic_cast<const MieGruneisenEOS*>(eos));
  }

  else if (dynamic_cast<const MieGruneisenEOSEnergy*>(eos)) {
    return std::make_unique<MieGruneisenEOSEnergy>(
      dynamic_cast<const MieGruneisenEOSEnergy*>(eos));
  }

  else if (dynamic_cast<const DefaultHypoElasticEOS*>(eos)) {
    return std::make_unique<DefaultHypoElasticEOS>(
      dynamic_cast<const DefaultHypoElasticEOS*>(eos));
  }

  else if (dynamic_cast<const HyperElasticEOS*>(eos)) {
    return std::make_unique<HyperElasticEOS>(
      dynamic_cast<const HyperElasticEOS*>(eos));
  }

  else if (dynamic_cast<const BorjaEOS*>(eos)) {
    return std::make_unique<BorjaEOS>(dynamic_cast<const BorjaEOS*>(eos));
  }

  else if (dynamic_cast<const AirEOS*>(eos)) {
    return std::make_unique<AirEOS>(dynamic_cast<const AirEOS*>(eos));
  }

  else if (dynamic_cast<const WaterEOS*>(eos)) {
    //return std::make_unique<WaterEOS>(dynamic_cast<const WaterEOS*>(eos));
    return (dynamic_cast<const WaterEOS*>(eos))->clone();
  }

  else if (dynamic_cast<const GraniteEOS*>(eos)) {
    return std::make_unique<GraniteEOS>(dynamic_cast<const GraniteEOS*>(eos));
  }

  else {
    proc0cout << "**WARNING** Creating a copy of the default hyperelastic "
                 "equation of state\n";
    return std::make_unique<HyperElasticEOS>(
      dynamic_cast<const HyperElasticEOS*>(eos));
  }

  return nullptr;
}
