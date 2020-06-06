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

#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardeningModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_Armstrong.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_None.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardening_Prager.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>
#include <string>
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;

using namespace Uintah;
using namespace Vaango;

KinematicHardeningModel*
KinematicHardeningModelFactory::create(ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("kinematic_hardening_model");
  if (!child) {
    cerr << "**WARNING** Creating default (no kinematic hardening) model"
         << endl;
    return (scinew KinematicHardening_None());
  }

  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for kinematic hardening model",
                                __FILE__, __LINE__);

  if (mat_type == "none")
    return (scinew KinematicHardening_None(child));
  else if (mat_type == "prager")
    return (scinew KinematicHardening_Prager(child));
  else if (mat_type == "armstrong_frederick")
    return (scinew KinematicHardening_Armstrong(child));
  else {
    cerr << "**WARNING** Creating default (no kinematic hardening) model"
         << endl;
    return (scinew KinematicHardening_None(child));
  }
}

KinematicHardeningModel*
KinematicHardeningModelFactory::create(ProblemSpecP& ps,
                                       InternalVariableModel* intvar)
{
  ProblemSpecP child = ps->findBlock("kinematic_hardening_model");
  if (!child) {
    cerr << "**WARNING** Creating default (no kinematic hardening) model"
         << endl;
    return (scinew KinematicHardening_None());
  }

  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException("No type for kinematic hardening model",
                                __FILE__, __LINE__);

  if (mat_type == "arena")
    return (scinew KinematicHardening_Arena(child, intvar));
  else {
    cerr << "**WARNING** Creating default (no kinematic hardening) model"
         << endl;
    return (scinew KinematicHardening_None(child));
  }
}

KinematicHardeningModel*
KinematicHardeningModelFactory::createCopy(const KinematicHardeningModel* pm)
{
  if (dynamic_cast<const KinematicHardening_None*>(pm)) {
    return (scinew KinematicHardening_None(
      dynamic_cast<const KinematicHardening_None*>(pm)));

  } else if (dynamic_cast<const KinematicHardening_Prager*>(pm)) {
    return (scinew KinematicHardening_Prager(
      dynamic_cast<const KinematicHardening_Prager*>(pm)));

  } else if (dynamic_cast<const KinematicHardening_Armstrong*>(pm)) {
    return (scinew KinematicHardening_Armstrong(
      dynamic_cast<const KinematicHardening_Armstrong*>(pm)));

  } else if (dynamic_cast<const KinematicHardening_Arena*>(pm)) {
    return (scinew KinematicHardening_Arena(
      dynamic_cast<const KinematicHardening_Arena*>(pm)));

  } else {
    cerr
      << "**WARNING** Creating copy of default (no kinematic hardening) model"
      << endl;
    return (scinew KinematicHardening_None(
      dynamic_cast<const KinematicHardening_None*>(pm)));
    // throw ProblemSetupException("Cannot create copy of unknown
    // kinematic_hardening model", __FILE__, __LINE__);
  }
}
