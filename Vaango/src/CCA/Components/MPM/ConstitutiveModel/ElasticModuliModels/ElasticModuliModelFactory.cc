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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_ArenaMixture.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arenisca.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Constant.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_Bulk.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_BulkPressure.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_NeuralNet.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_NeuralNet_Bulk.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_SupportVector.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <string>

using namespace Vaango;

ElasticModuliModel*
ElasticModuliModelFactory::create(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP child = ps->findBlock("elastic_moduli_model");
  if (!child) {
    std::ostringstream out;
    out << "**Error** No Elastic modulus model provided."
        << " Default (constant elasticity) model needs at least two input "
           "parameters."
        << std::endl;
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::string mat_type;
  if (!child->getAttribute("type", mat_type)) {
    std::ostringstream out;
    out << "MPM::ConstitutiveModel:No type provided for elasticity model.";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  if (mat_type == "constant")
    return (scinew ElasticModuli_Constant(child));
  else if (mat_type == "arenisca")
    return (scinew ElasticModuli_Arenisca(child));
  else if (mat_type == "arena")
    return (scinew ElasticModuli_Arena(child));
  else if (mat_type == "arena_mixture")
    return (scinew ElasticModuli_ArenaMixture(child));
  else if (mat_type == "tabular")
    return (scinew ElasticModuli_Tabular(child));
  else if (mat_type == "tabular_bulk")
    return (scinew ElasticModuli_Tabular_Bulk(child));
  else if (mat_type == "tabular_bulk_pressure")
    return (scinew ElasticModuli_Tabular_BulkPressure(child));
  else if (mat_type == "neural_net")
    return (scinew ElasticModuli_NeuralNet(child));
  else if (mat_type == "neural_net_bulk")
    return (scinew ElasticModuli_NeuralNet_Bulk(child));
  else if (mat_type == "support_vector")
    return (scinew ElasticModuli_SupportVector(child));
  else {
    std::cerr << "**WARNING** No elasticity model provided. "
              << "Creating default (constant elasticity) model" << std::endl;
    return (scinew ElasticModuli_Constant(child));
  }
}

ElasticModuliModel*
ElasticModuliModelFactory::createCopy(const ElasticModuliModel* model)
{
  if (dynamic_cast<const ElasticModuli_Constant*>(model))
    return (scinew ElasticModuli_Constant(
      dynamic_cast<const ElasticModuli_Constant*>(model)));
  else if (dynamic_cast<const ElasticModuli_Arenisca*>(model))
    return (scinew ElasticModuli_Arenisca(
      dynamic_cast<const ElasticModuli_Arenisca*>(model)));
  else if (dynamic_cast<const ElasticModuli_Arena*>(model))
    return (scinew ElasticModuli_Arena(
      dynamic_cast<const ElasticModuli_Arena*>(model)));
  else if (dynamic_cast<const ElasticModuli_ArenaMixture*>(model))
    return (scinew ElasticModuli_ArenaMixture(
      dynamic_cast<const ElasticModuli_ArenaMixture*>(model)));
  else if (dynamic_cast<const ElasticModuli_Tabular*>(model))
    return (scinew ElasticModuli_Tabular(
      dynamic_cast<const ElasticModuli_Tabular*>(model)));
  else if (dynamic_cast<const ElasticModuli_Tabular_Bulk*>(model))
    return (scinew ElasticModuli_Tabular_Bulk(
      dynamic_cast<const ElasticModuli_Tabular_Bulk*>(model)));
  else if (dynamic_cast<const ElasticModuli_Tabular_BulkPressure*>(model))
    return (scinew ElasticModuli_Tabular_BulkPressure(
      dynamic_cast<const ElasticModuli_Tabular_BulkPressure*>(model)));
  else if (dynamic_cast<const ElasticModuli_NeuralNet*>(model))
    return (scinew ElasticModuli_NeuralNet(
      dynamic_cast<const ElasticModuli_NeuralNet*>(model)));
  else if (dynamic_cast<const ElasticModuli_NeuralNet_Bulk*>(model))
    return (scinew ElasticModuli_NeuralNet_Bulk(
      dynamic_cast<const ElasticModuli_NeuralNet_Bulk*>(model)));
  else if (dynamic_cast<const ElasticModuli_SupportVector*>(model))
    return (scinew ElasticModuli_SupportVector(
      dynamic_cast<const ElasticModuli_SupportVector*>(model)));
  else {
    std::cerr << "**WARNING** No elasticity model provided. "
              << "Creating default (constant elasticity) model" << std::endl;
    return (scinew ElasticModuli_Constant(
      dynamic_cast<const ElasticModuli_Constant*>(model)));
  }
}
