/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_ArenaMixture.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Arenisca3.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_CamClay.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Gurson.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Rousselier.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_vonMises.h>
#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldConditionFactory.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>

using namespace Uintah;
using namespace Vaango;

/// Create an instance of a Yield Condition.
std::unique_ptr<YieldCondition>
YieldConditionFactory::create(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP child = ps->findBlock("yield_condition");
  if (!child)
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Cannot find yield condition.",
      __FILE__,
      __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:No type for yield condition.",
      __FILE__,
      __LINE__);

  if (mat_type == "arena")
    return (std::make_unique<YieldCond_Arena>(child));
  else if (mat_type == "arena_mixture")
    return (std::make_unique<YieldCond_ArenaMixture>(child));
  else if (mat_type == "arenisca3")
    return (std::make_unique<YieldCond_Arenisca3>(child));
  else
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Unknown Yield Condition (" + mat_type + ")",
      __FILE__,
      __LINE__);
}

// Create an instance of a Yield Condition with an associated internal
// variable model of type IntVar_Metal.
std::unique_ptr<YieldCondition>
YieldConditionFactory::create(Uintah::ProblemSpecP& ps,
                              IntVar_Metal* intvar,
                              const Uintah::FlowStressModel* flow)
{
  ProblemSpecP child = ps->findBlock("yield_condition");
  if (!child)
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Cannot find yield condition.",
      __FILE__, __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:No type for yield condition.",
      __FILE__, __LINE__);

  if (mat_type == "gurson")
    return (std::make_unique<YieldCond_Gurson>(child, intvar, flow));
  else if (mat_type == "rousselier")
    return (std::make_unique<YieldCond_Rousselier>(child, intvar, flow));
  else if (mat_type == "von_mises")
    return (std::make_unique<YieldCond_vonMises>(child, intvar, flow));
  else
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Unknown Yield Condition (" + mat_type + ")",
      __FILE__, __LINE__);
}

// Create an instance of a Yield Condition with an associated internal
// variable model of type IntVar_BorjaPressure.
std::unique_ptr<YieldCondition>
YieldConditionFactory::create(Uintah::ProblemSpecP& ps,
                              IntVar_BorjaPressure* intvar)
{
  ProblemSpecP child = ps->findBlock("yield_condition");
  if (!child)
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Cannot find yield condition.",
      __FILE__, __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:No type for yield condition.",
      __FILE__, __LINE__);

  else if (mat_type == "camclay")
    return (std::make_unique<YieldCond_CamClay>(child, intvar));
  else
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Unknown Yield Condition (" + mat_type + ")",
      __FILE__, __LINE__);
}

// Create an instance of a Yield Condition with an associated internal
// variable model of type IntVar_TabularCap.
std::unique_ptr<YieldCondition>
YieldConditionFactory::create(Uintah::ProblemSpecP& ps,
                              IntVar_TabularCap* intvar)
{
  ProblemSpecP child = ps->findBlock("yield_condition");
  if (!child)
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Cannot find yield condition.",
      __FILE__, __LINE__);
  string mat_type;
  if (!child->getAttribute("type", mat_type))
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:No type for yield condition.",
      __FILE__, __LINE__);

  else if (mat_type == "tabular")
    return (std::make_unique<YieldCond_Tabular>(child, nullptr));
  else if (mat_type == "tabular_cap")
    return (std::make_unique<YieldCond_TabularCap>(child, intvar));
  else
    throw ProblemSetupException(
      "MPM::ConstitutiveModel:Unknown Yield Condition (" + mat_type + ")",
      __FILE__, __LINE__);
}

std::unique_ptr<YieldCondition>
YieldConditionFactory::createCopy(const YieldCondition* yc)
{
  if (dynamic_cast<const YieldCond_Arena*>(yc))
    return (std::make_unique<YieldCond_Arena>(dynamic_cast<const YieldCond_Arena*>(yc)));

  else if (dynamic_cast<const YieldCond_ArenaMixture*>(yc))
    return (std::make_unique<YieldCond_ArenaMixture>(
      dynamic_cast<const YieldCond_ArenaMixture*>(yc)));

  else if (dynamic_cast<const YieldCond_Arenisca3*>(yc))
    return (
      std::make_unique<YieldCond_Arenisca3>(dynamic_cast<const YieldCond_Arenisca3*>(yc)));

  else if (dynamic_cast<const YieldCond_CamClay*>(yc))
    return (
      std::make_unique<YieldCond_CamClay>(dynamic_cast<const YieldCond_CamClay*>(yc)));

  else if (dynamic_cast<const YieldCond_Gurson*>(yc))
    return (std::make_unique<YieldCond_Gurson>(dynamic_cast<const YieldCond_Gurson*>(yc)));

  else if (dynamic_cast<const YieldCond_Rousselier*>(yc))
    return (std::make_unique<YieldCond_Rousselier>(dynamic_cast<const YieldCond_Rousselier*>(yc)));

  else if (dynamic_cast<const YieldCond_Tabular*>(yc))
    return (
      std::make_unique<YieldCond_Tabular>(dynamic_cast<const YieldCond_Tabular*>(yc)));

  else if (dynamic_cast<const YieldCond_TabularCap*>(yc))
    return (std::make_unique<YieldCond_TabularCap>(
      dynamic_cast<const YieldCond_TabularCap*>(yc)));

  else if (dynamic_cast<const YieldCond_vonMises*>(yc))
    return (
      std::make_unique<YieldCond_vonMises>(dynamic_cast<const YieldCond_vonMises*>(yc)));

  else
    throw ProblemSetupException(
      "Cannot create copy of unknown yield condition", __FILE__, __LINE__);
}
