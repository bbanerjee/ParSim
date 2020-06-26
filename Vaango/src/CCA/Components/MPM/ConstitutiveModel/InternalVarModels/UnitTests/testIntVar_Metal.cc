#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_Metal.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Math/Matrix3.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Malloc/Allocator.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;
using Uintah::Matrix3;

TEST(IntVarMetalTest, constructors)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "internal_variable_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_plasticity");
  xmlDocSetRootElement(doc, rootNode);

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors
  IntVar_Metal model(ps);
  IntVar_Metal model_copy(model);

  // Get labels
  auto labels = model.getLabels();
  for (auto label : labels) {
    std::cout << "Label = " << label->getName() << "\n";
  }

  // Create an empty problem spec
  xmlDocPtr out_doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr out_rootNode = xmlNewNode(nullptr, BAD_CAST "OutputPS");
  xmlDocSetRootElement(out_doc, out_rootNode);
  ProblemSpecP out_ps = scinew ProblemSpec(xmlDocGetRootElement(out_doc), false);
  model.outputProblemSpec(out_ps);
  xmlSaveFormatFileEnc("-", out_ps->getNode()->doc, "ISO-8859-1", 1);
  

  /*
  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model_copy.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;

  ASSERT_NEAR(moduli.bulkModulus, 1.0e6, 1.0e-6);
  ASSERT_NEAR(moduli.shearModulus, 0.7e6, 1.0e-6);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_NEAR(moduli.bulkModulus, 1.0e6, 1.0e-6);
  ASSERT_NEAR(moduli.shearModulus, 0.7e6, 1.0e-6);

  moduli = model.getElasticModuliLowerBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_NEAR(moduli.bulkModulus, 1.0e6, 1.0e-6);
  ASSERT_NEAR(moduli.shearModulus, 0.7e6, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1000.0;
  state.meltingTemp = 500.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };

  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        state.pressure = p;
        state.density = rho;
        state.temperature = T;
        moduli = model.getCurrentElasticModuli(&state);
        //std::cout << std::setprecision(16)
        //          << " p = " << state.pressure
        //          << " rho = " << state.density 
        //          << " T = " << state.temperature << " K = " << moduli.bulkModulus
        //          << " G = " << moduli.shearModulus << "\n";
        ASSERT_NEAR(moduli.bulkModulus, 1.0e6, 1.0e-6);
        ASSERT_NEAR(moduli.shearModulus, 0.7e6, 1.0e-6);
      }
    }
  }
  */
}

