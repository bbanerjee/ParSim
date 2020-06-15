#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_MetalIso.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Math/Matrix3.h>

#include <Core/Exceptions/ProblemSetupException.h>
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

TEST(ElasticModuliMetalIsoTest, constantBulkConstantShear)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "default_hypo");
  xmlNewChild(eosNode, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e6");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "0.7e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors
  ElasticModuli_MetalIso model(ps);
  ElasticModuli_MetalIso model_copy(model);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["bulk_modulus"] =  1.0e6;
  test_params["shear_modulus"] =  0.7e6;

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

  //std::map<double, std::pair<double, double>> I1_KG;
  //I1_KG[-0.01] = std::make_pair(116348059.2816119, 173990384.2029271);
  //I1_KG[-1] = std::make_pair(116407266.5941095, 174078654.5535643);
  //I1_KG[-1e+06] = std::make_pair(201335456.1999204, 300413920.3418134);
  //I1_KG[-1e+08] = std::make_pair(1046521964.819405, 1527736301.547844);
  //I1_KG[-1e+11] = std::make_pair(74552697327.11647, 75257479257.76637);
  //ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
  //ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);

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
}

TEST(ElasticModuliMetalIsoTest, constantBulkShearBorja)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "default_hypo");
  xmlNewChild(eosNode, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e6");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "0.7e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors
  ElasticModuli_MetalIso model(ps);
  ElasticModuli_MetalIso model_copy(model);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["bulk_modulus"] =  1.0e6;
  test_params["shear_modulus"] =  0.7e6;

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

  //std::map<double, std::pair<double, double>> I1_KG;
  //I1_KG[-0.01] = std::make_pair(116348059.2816119, 173990384.2029271);
  //I1_KG[-1] = std::make_pair(116407266.5941095, 174078654.5535643);
  //I1_KG[-1e+06] = std::make_pair(201335456.1999204, 300413920.3418134);
  //I1_KG[-1e+08] = std::make_pair(1046521964.819405, 1527736301.547844);
  //I1_KG[-1e+11] = std::make_pair(74552697327.11647, 75257479257.76637);
  //ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
  //ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);

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
}
