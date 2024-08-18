/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#include "DamageModelFactory.h"
#include "HancockMacKenzieDamage.h"
#include "JohnsonCookDamage.h"
#include "NullDamage.h"
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Parallel/Parallel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>

using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;

using namespace Uintah;

std::unique_ptr<DamageModel>
DamageModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("damage_model");
  if (!child) {
    proc0cout << "**WARNING** Creating default null damage model" << std::endl;
    return (std::make_unique<NullDamage>());
    // throw ProblemSetupException("Cannot find damage_model tag", __FILE__,
    // __LINE__);
  }
  string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw ProblemSetupException("No type for damage_model", __FILE__, __LINE__);
  }

  if (mat_type == "johnson_cook") {
    return (std::make_unique<JohnsonCookDamage>(child));
  } else if (mat_type == "hancock_mackenzie") {
    return (std::make_unique<HancockMacKenzieDamage>(child));
  } else {
    proc0cout << "**WARNING** Creating default null damage model" << std::endl;
    return (std::make_unique<NullDamage>(child));
    // throw ProblemSetupException("Unknown Damage Model ("+mat_type+")",
    // __FILE__, __LINE__);
  }

  // return 0;
}

std::unique_ptr<DamageModel>
DamageModelFactory::createCopy(const DamageModel* dm)
{
  if (dynamic_cast<const JohnsonCookDamage*>(dm)) {
    return (std::make_unique<JohnsonCookDamage>(
      dynamic_cast<const JohnsonCookDamage*>(dm)));
  }

  else if (dynamic_cast<const HancockMacKenzieDamage*>(dm)) {
    return (std::make_unique<HancockMacKenzieDamage>(
      dynamic_cast<const HancockMacKenzieDamage*>(dm)));
  }

  else {
    proc0cout << "**WARNING** Creating copy of default null damage model"
              << std::endl;
    return (std::make_unique<NullDamage>(dynamic_cast<const NullDamage*>(dm)));
    // throw ProblemSetupException("Cannot create copy of unknown damage model",
    // __FILE__, __LINE__);
  }

  // return 0;
}
