#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_NeuralNet.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;
using Uintah::ProblemSetupException;
using Uintah::InvalidValue;
using nlohmann::json;

TEST(ElasticModuliNeuralNetTest, constructorTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string json_file = std::string(currPath) + "/" + "mlp_regression_keras_total_scaled.h5";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "neural_net");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "mlp_regression_keras_total_scaled.h5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "min_strain", 
              BAD_CAST "0.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "max_strain", 
              BAD_CAST "0.452");
  xmlNewChild(rootNode, nullptr, BAD_CAST "min_pressure", 
              BAD_CAST "0.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "max_pressure", 
              BAD_CAST "1.0e6");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", 
              BAD_CAST "1.0e4");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu", 
              BAD_CAST "0.2");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }


  // Create a model
  ElasticModuli_NeuralNet model(ps);
  try {
    ElasticModuli moduli = model.getInitialElasticModuli();
    ASSERT_NEAR(moduli.bulkModulus, 4.44376499e+08, 1);

    ModelState_Tabular state_init;
    state_init.elasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    state_init.plasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state_init);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_NEAR(KG.bulkModulus, 4.44376499e+08, 1);
    EXPECT_NEAR(KG.shearModulus, 3.33282374e8, 1);
    EXPECT_NEAR(dKdG.bulkModulus, -5398166832, 1);
    ASSERT_NEAR(dKdG.shearModulus, -4048625123, 1);

  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  // Copy
  ElasticModuli_NeuralNet modelCopy(&model);
  try {
    ElasticModuli moduli = modelCopy.getInitialElasticModuli();
    ASSERT_NEAR(moduli.bulkModulus, 4.44376499e+08, 1);
    ASSERT_NEAR(moduli.shearModulus, 3.33282374e8, 1);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  // Modelstate test
  ModelState_Tabular state;
  state.elasticStrainTensor = Uintah::Matrix3(-0.03, 0, 0, 0, -0.03, 0, 0, 0, -0.03);
  state.plasticStrainTensor = Uintah::Matrix3(-0.02, 0, 0, 0, -0.02, 0, 0, 0, -0.02);
  try {
    ElasticModuli moduli = model.getCurrentElasticModuli(&state);
    ASSERT_NEAR(moduli.bulkModulus,  5.71826722e+09, 2.0);
    ASSERT_NEAR(moduli.shearModulus, 4.288700413e+09, 1.0);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;

    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_NEAR(KG.bulkModulus,  5.71826722e+09, 2.0);
    EXPECT_NEAR(KG.shearModulus, 4.288700413e+09, 1.0);
    EXPECT_NEAR(dKdG.bulkModulus, -72637200355, 1.0);
    ASSERT_NEAR(dKdG.shearModulus, -54477900266, 1.0);
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  state.elasticStrainTensor = Uintah::Matrix3(0.03, 0, 0, 0, 0.03, 0, 0, 0, 0.03);
  state.plasticStrainTensor = Uintah::Matrix3(0.02, 0, 0, 0, 0.02, 0, 0, 0, 0.02);
  try {
    ElasticModuli moduli = model.getCurrentElasticModuli(&state);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
    //ASSERT_NEAR(moduli.bulkModulus,  93178465.39, 1.0);
    //ASSERT_NEAR(moduli.shearModulus, 69883849.037, 1.0);
    ASSERT_NEAR(moduli.bulkModulus,  35461413.9, 1.0);
    ASSERT_NEAR(moduli.shearModulus, 26596060.4, 1.0);

    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_NEAR(KG.bulkModulus,  35461413.9, 1.0);
    EXPECT_NEAR(KG.shearModulus, 26596060.4, 1.0);
    EXPECT_NEAR(dKdG.bulkModulus, -548735260, 1.0);
    ASSERT_NEAR(dKdG.shearModulus, -411551445, 1.0);
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }
}
