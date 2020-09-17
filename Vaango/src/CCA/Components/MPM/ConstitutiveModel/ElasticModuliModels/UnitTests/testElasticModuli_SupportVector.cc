#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_SupportVector.h>
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

TEST(ElasticModuliSVRTest, constructorTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string json_file = std::string(currPath) + "/" + "ARL_Sand_SVR_fit_10_001.json";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "support_vector");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "svr_filename", 
              BAD_CAST "ARL_Sand_SVR_fit_10_001.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", 
              BAD_CAST "1.0e4");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu", 
              BAD_CAST "0.2");

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a model
  try {
    ElasticModuli_SupportVector model_test(ps);
  } catch (Uintah::ProblemSetupException& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::InvalidValue& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::Exception& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (...) {
    std::cout << "**ERROR** Unknown exception\n";
    throw;
  }

  ElasticModuli_SupportVector model(ps);
  try {
    ElasticModuli moduli = model.getInitialElasticModuli();
    std::cout << "Initial moduli = " << moduli.bulkModulus << ", " << moduli.shearModulus << "\n";
    EXPECT_DOUBLE_EQ(moduli.bulkModulus, 336515295.5149824);

    ModelState_Tabular state_init;
    state_init.elasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    state_init.plasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state_init);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_DOUBLE_EQ(KG.bulkModulus, 336515295.59960222);
    EXPECT_DOUBLE_EQ(KG.shearModulus, 252386471.69970167);
    EXPECT_DOUBLE_EQ(dKdG.bulkModulus, 25346817562.290344);
    ASSERT_DOUBLE_EQ(dKdG.shearModulus, 19010113171.717758);
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  // Copy
  ElasticModuli_SupportVector modelCopy(&model);
  try {
    ElasticModuli moduli = modelCopy.getInitialElasticModuli();
    EXPECT_DOUBLE_EQ(moduli.bulkModulus, 336515295.5149824);
    EXPECT_DOUBLE_EQ(moduli.shearModulus, 252386471.63623676);
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
    EXPECT_NEAR(moduli.bulkModulus, 336474416.11733419, 1.0e-3);
    EXPECT_NEAR(moduli.shearModulus, 252355812.08800063, 1.0e-3);

    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_NEAR(KG.bulkModulus, 11440, 1.0e-7);
    EXPECT_NEAR(KG.shearModulus, 8580, 1.0e-7);
    EXPECT_NEAR(dKdG.bulkModulus, -24000, 1.0);
    ASSERT_NEAR(dKdG.shearModulus, -18000, 1.0);

    // Compute tangent modulus
    auto tangent = model.computeElasticTangentModulus(&state);
    double K = KG.bulkModulus;
    double G = KG.shearModulus;
    double K43G = K + 4.0 * G / 3.0;
    double K23G = K - 2.0 * G / 3.0;
    //std::cout << "K = " << K << " G = " << G
    //          << "K43G = " << K43G << " K23G = " << K23G
    //          << "\nTangent = \n" << tangent << "\n";
    ASSERT_DOUBLE_EQ(tangent(1,1), K43G);
    ASSERT_DOUBLE_EQ(tangent(1,2), K23G);
    ASSERT_DOUBLE_EQ(tangent(4,4), G);

  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }
  
}
