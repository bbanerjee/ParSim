/*
 * The MIT License
 *
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

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModelFactory.h>

#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/MaterialModels/IsotropicElasticNeoHookeanStateModel.h>
#include <CCA/Components/Peridynamics/MaterialModels/LinearElasticBondModel.h>
#include <CCA/Components/Peridynamics/MaterialModels/PolarOrthotropicLinearElasticStateModel.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace Vaango;

PeridynamicsMaterialModel*
PeridynamicsMaterialModelFactory::create(Uintah::ProblemSpecP& ps,
                                         PeridynamicsFlags* flags)
{
  Uintah::ProblemSpecP child = ps->findBlock("material_model");
  if (!child) {
    throw Uintah::ProblemSetupException(
      "Cannot find material_model tag", __FILE__, __LINE__);
  }
  std::string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    throw Uintah::ProblemSetupException(
      "No type for material_model", __FILE__, __LINE__);
  }

  if (mat_type == "linear_elastic_bond") {
    return (scinew LinearElasticBondModel(child, flags));
  } else if (mat_type == "elastic_neo_hookean_state") {
    return (scinew IsotropicElasticNeoHookeanStateModel(child, flags));
  } else if (mat_type == "polar_orthotropic_linear_elastic_state") {
    return (scinew PolarOrthotropicLinearElasticStateModel(child, flags));
  } else {
    throw Uintah::ProblemSetupException("Unknown peridynamic material type (" +
                                          mat_type + ")",
                                        __FILE__,
                                        __LINE__);
  }

  return 0;
}
