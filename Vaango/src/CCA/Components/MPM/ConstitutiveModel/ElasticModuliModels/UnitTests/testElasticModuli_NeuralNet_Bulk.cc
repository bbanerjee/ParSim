#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_NeuralNet_Bulk.h>
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

TEST(ElasticModuliNeuralNetBulkTest, constructorTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string json_file = std::string(currPath) + "/" + "mlp_regression_keras_bulk_elastic_ten_com_scaled_relu_32_400.h5";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "neural_net_bulk");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "mlp_regression_keras_bulk_elastic_ten_com_scaled_relu_32_400.h5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "mean_elastic_strain", 
              BAD_CAST "0.03547612");
  xmlNewChild(rootNode, nullptr, BAD_CAST "std_dev_elastic_strain", 
              BAD_CAST "0.05561166");
  xmlNewChild(rootNode, nullptr, BAD_CAST "mean_plastic_strain", 
              BAD_CAST "0.21299634");
  xmlNewChild(rootNode, nullptr, BAD_CAST "std_dev_plastic_strain", 
              BAD_CAST "0.14141187");
  xmlNewChild(rootNode, nullptr, BAD_CAST "mean_bulk_modulus", 
              BAD_CAST "1.93588029e+10");
  xmlNewChild(rootNode, nullptr, BAD_CAST "std_dev_bulk_modulus", 
              BAD_CAST "1.89843194e+10");
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

  // Test model creation
  try {
    ElasticModuli_NeuralNet_Bulk model(ps);
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
    std::cout << "**ERROR** Unknown exception during model creation\n";
    throw;
  }

  // Create a model
  ElasticModuli_NeuralNet_Bulk model(ps);
  try {
    ElasticModuli moduli = model.getInitialElasticModuli();
    ASSERT_DOUBLE_EQ(moduli.bulkModulus, 531873380.8029747);

    ModelState_Tabular state_init;
    state_init.elasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    state_init.plasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state_init);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_DOUBLE_EQ(KG.bulkModulus, 531873380.8029747);
    EXPECT_DOUBLE_EQ(KG.shearModulus, 398905035.60223103);
    EXPECT_DOUBLE_EQ(dKdG.bulkModulus, -2000964164.7338867);
    ASSERT_DOUBLE_EQ(dKdG.shearModulus, -1500723123.550415);

  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
    throw;
  }

  // Copy
  ElasticModuli_NeuralNet_Bulk modelCopy(&model);
  try {
    ElasticModuli moduli = modelCopy.getInitialElasticModuli();
    EXPECT_DOUBLE_EQ(moduli.bulkModulus, 531873380.8029747);
    EXPECT_DOUBLE_EQ(moduli.shearModulus, 398905035.60223103);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
    throw;
  }

  // Modelstate test
  ModelState_Tabular state;
  state.elasticStrainTensor = Uintah::Matrix3(-0.03, 0, 0, 0, -0.03, 0, 0, 0, -0.03);
  state.plasticStrainTensor = Uintah::Matrix3(-0.02, 0, 0, 0, -0.02, 0, 0, 0, -0.02);
  try {
    ElasticModuli moduli = model.getCurrentElasticModuli(&state);
    ASSERT_DOUBLE_EQ(moduli.bulkModulus,  2979061231.3637466);
    ASSERT_DOUBLE_EQ(moduli.shearModulus, 2234295923.52281);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;

    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_DOUBLE_EQ(KG.bulkModulus,  2979061231.3637466);
    EXPECT_DOUBLE_EQ(KG.shearModulus, 2234295923.52281);
    EXPECT_DOUBLE_EQ(dKdG.bulkModulus, 1492763710.0219727);
    ASSERT_DOUBLE_EQ(dKdG.shearModulus, 1119572782.5164795);
  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
    throw;
  }

  state.elasticStrainTensor = Uintah::Matrix3(0.03, 0, 0, 0, 0.03, 0, 0, 0, 0.03);
  state.plasticStrainTensor = Uintah::Matrix3(0.02, 0, 0, 0, 0.02, 0, 0, 0, 0.02);
  try {
    ElasticModuli moduli = model.getCurrentElasticModuli(&state);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
    ASSERT_DOUBLE_EQ(moduli.bulkModulus,  531873380.8029747);
    ASSERT_DOUBLE_EQ(moduli.shearModulus, 398905035.60223103);

    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_DOUBLE_EQ(KG.bulkModulus,  531873380.8029747);
    EXPECT_DOUBLE_EQ(KG.shearModulus, 398905035.60223103);
    EXPECT_DOUBLE_EQ(dKdG.bulkModulus, 0.0);
    ASSERT_DOUBLE_EQ(dKdG.shearModulus, 0.0);

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
  } catch (const Uintah::InvalidValue& e) {
    std::cout << e.message() << std::endl;
    throw;
  }
}
