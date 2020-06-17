#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_MetalIso.h>
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
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "borja_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu0", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors (check wrong init)
  EXPECT_THROW({
    try {
      ElasticModuli_MetalIso model(ps);
    } catch (Uintah::ProblemSetupException e) {
      // std::cout << e.message() << std::endl;
      throw;
    }
  }, Uintah::ProblemSetupException);
}

TEST(ElasticModuliMetalIsoTest, constantBulkShearMTS)
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
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "mts_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "28.0e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "D", BAD_CAST "4.5e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "T_0", BAD_CAST "294.0");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  28.0e9;
  test_params["D"] =  4.5e9;
  test_params["T_0"] =  294;
  test_params["bulk_modulus"] =  1.0e6;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;

  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 28000000000);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 28000000000);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1000.0;
  state.meltingTemp = 500.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27748825708.72905);
  p_KG[keys[ii++]] = std::make_pair(1000000, 22377699014.68184);
  p_KG[keys[ii++]] = std::make_pair(1000000, 14833786051.01098);

  //p_KG[-1e+11] = std::make_pair(74552697327.11647, 75257479257.76637);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, constantBulkShearNadal)
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
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "np_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "26.5e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "zeta", BAD_CAST "0.04");
  xmlNewChild(shearNode, nullptr, BAD_CAST "slope_mu_p_over_mu0", BAD_CAST "65.0e-12");
  xmlNewChild(shearNode, nullptr, BAD_CAST "C", BAD_CAST "0.047");
  xmlNewChild(shearNode, nullptr, BAD_CAST "m", BAD_CAST "26.98");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  26.5e9;
  test_params["zeta"] =  0.04;
  test_params["slope_mu_p_over_mu0"] =  65.0e-12;
  test_params["C"] =  0.047;
  test_params["m"] =  26.98;
  test_params["bulk_modulus"] =  1.0e6;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 26500000000);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 26506185736.60703);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 339800594.3758016);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1000.0;
  state.meltingTemp = 500.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
  p_KG[keys[ii++]] = std::make_pair(1000000, 21265571354.59132);
  p_KG[keys[ii++]] = std::make_pair(1000000, 163920965.387606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 21855685239.16529);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1639209653.87606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27756839254.62606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 16392096538.7606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 21265568385.78031);
  p_KG[keys[ii++]] = std::make_pair(1000000, 163920965.387606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 21855683861.16529);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1639209653.87606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27756838615.01512);
  p_KG[keys[ii++]] = std::make_pair(1000000, 16392096538.7606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 21265565416.96931);
  p_KG[keys[ii++]] = std::make_pair(1000000, 163920965.387606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 21855682483.16529);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1639209653.87606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);
  p_KG[keys[ii++]] = std::make_pair(1000000, 27756837975.40418);
  p_KG[keys[ii++]] = std::make_pair(1000000, 16392096538.7606);
  p_KG[keys[ii++]] = std::make_pair(1000000, 1e-08);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, constantBulkShearPTW)
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
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "ptw_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "895e8");
  xmlNewChild(shearNode, nullptr, BAD_CAST "alpha", BAD_CAST "0.23");
  xmlNewChild(shearNode, nullptr, BAD_CAST "alphap", BAD_CAST "0.66");
  xmlNewChild(shearNode, nullptr, BAD_CAST "slope_mu_p_over_mu0", BAD_CAST "26.0e-12");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  895e8;
  test_params["alpha"] =  0.23;
  test_params["alphap"] =  0.66;
  test_params["slope_mu_p_over_mu0"] =  26.0e-12;
  test_params["bulk_modulus"] =  1.0e6;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 89500000000);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 89499498662.45692);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 30435906999.96326);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 2300.0;
  state.meltingTemp = 1700.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025300478.40388);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126475921.30559);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752945224.93272);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025297070.04886);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126473063.62726);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752943055.60025);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025295488.03056);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126471737.21045);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752942048.68531);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025294117.64706);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126470588.23529);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752941176.47059);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025294117.64706);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126470588.23529);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752941176.47059);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025294117.64706);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126470588.23529);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752941176.47059);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025287756.89023);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126465255.16501);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752937128.00847);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025291165.24527);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126468112.84334);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752939297.34093);
   p_KG[keys[ii++]] = std::make_pair(1000000, 86025292747.26353);
   p_KG[keys[ii++]] = std::make_pair(1000000, 72126469439.26013);
   p_KG[keys[ii++]] = std::make_pair(1000000, 54752940304.25587);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, constantBulkShearSCG)
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
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "scg_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "27.6e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "A", BAD_CAST "65.0e-12");
  xmlNewChild(shearNode, nullptr, BAD_CAST "B", BAD_CAST "0.62e-3");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  27.6e9;
  test_params["A"] =  65.0e-12;
  test_params["B"] =  0.62e-3;
  test_params["bulk_modulus"] =  1.0e6;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 27600000000);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 32733042374.41661);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1500.0;
  state.meltingTemp = 700.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022404424.38447);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177604424.38446);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621604424.38446);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022402053.61735);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177602053.61735);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621602053.61735);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022400953.20474);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177600953.20474);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621600953.20474);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022395575.61554);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177595575.61554);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621595575.61554);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022397946.38265);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177597946.38265);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621597946.38265);
   p_KG[keys[ii++]] = std::make_pair(1000000, 31022399046.79527);
   p_KG[keys[ii++]] = std::make_pair(1000000, 24177599046.79526);
   p_KG[keys[ii++]] = std::make_pair(1000000, 15621599046.79527);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, bulkAirShearConstant)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For Air EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "air");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors (check wrong init)
  EXPECT_THROW({
    try {
      ElasticModuli_MetalIso model(ps);
    } catch (Uintah::ProblemSetupException e) {
      //std::cout << e.message() << std::endl;
      throw;
    }
  }, Uintah::ProblemSetupException);
}

TEST(ElasticModuliMetalIsoTest, bulkWaterShearConstant)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For Air EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "water");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors (check wrong init)
  EXPECT_THROW({
    try {
      ElasticModuli_MetalIso model(ps);
    } catch (Uintah::ProblemSetupException e) {
      //std::cout << e.message() << std::endl;
      throw;
    }
  }, Uintah::ProblemSetupException);
}

TEST(ElasticModuliMetalIsoTest, bulkBorjaShearConstant)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For Air EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "borja_pressure");
  xmlNewChild(eosNode, nullptr, BAD_CAST "p0", BAD_CAST "-9.0");
  xmlNewChild(eosNode, nullptr, BAD_CAST "alpha", BAD_CAST "60.0");
  xmlNewChild(eosNode, nullptr, BAD_CAST "kappatilde", BAD_CAST "0.018");
  xmlNewChild(eosNode, nullptr, BAD_CAST "epse_v0", BAD_CAST "0.0");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors (check wrong init)
  EXPECT_THROW({
    try {
      ElasticModuli_MetalIso model(ps);
    } catch (Uintah::ProblemSetupException e) {
      //std::cout << e.message() << std::endl;
      throw;
    }
  }, Uintah::ProblemSetupException);
}

TEST(ElasticModuliMetalIsoTest, graniteBulkConstantShear)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "granite");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "28.0e9");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  EXPECT_THROW({
    try {
      ElasticModuli_MetalIso model(ps);
    } catch (Uintah::ProblemSetupException e) {
      //std::cout << e.message() << "\n";
      throw;
    }
  }, Uintah::ProblemSetupException);
}

TEST(ElasticModuliMetalIsoTest, hyperBulkShearSCG)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "default_hyper");
  xmlNewChild(eosNode, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e6");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "scg_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "27.6e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "A", BAD_CAST "65.0e-12");
  xmlNewChild(shearNode, nullptr, BAD_CAST "B", BAD_CAST "0.62e-3");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  27.6e9;
  test_params["A"] =  65.0e-12;
  test_params["B"] =  0.62e-3;
  test_params["bulk_modulus"] =  1.0e6;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1000000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 27600000000);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 5000500000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 32733042374.41661);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 500050);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1500.0;
  state.meltingTemp = 700.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 31022404424.38447);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 24177604424.38446);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 15621604424.38446);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 31022402053.61735);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 24177602053.61735);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 15621602053.61735);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 31022400953.20474);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 24177600953.20474);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 15621600953.20474);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 31022395575.61554);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 24177595575.61554);
   p_KG[keys[ii++]] = std::make_pair(502222.2222222222, 15621595575.61554);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 31022397946.38265);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 24177597946.38265);
   p_KG[keys[ii++]] = std::make_pair(722222.2222222222, 15621597946.38265);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 31022399046.79527);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 24177599046.79526);
   p_KG[keys[ii++]] = std::make_pair(22722222.22222222, 15621599046.79527);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, mieBulkShearSCG)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "mie_gruneisen");
  xmlNewChild(eosNode, nullptr, BAD_CAST "C_0", BAD_CAST "5386");
  xmlNewChild(eosNode, nullptr, BAD_CAST "Gamma_0", BAD_CAST "1.99");
  xmlNewChild(eosNode, nullptr, BAD_CAST "S_alpha", BAD_CAST "1.339");
  xmlNewChild(eosNode, nullptr, BAD_CAST "rho_0", BAD_CAST "1500");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "scg_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "27.6e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "A", BAD_CAST "65.0e-12");
  xmlNewChild(shearNode, nullptr, BAD_CAST "B", BAD_CAST "0.62e-3");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  27.6e9;
  test_params["A"] =  65.0e-12;
  test_params["B"] =  0.62e-3;
  test_params["C_0"] =  5386;
  test_params["Gamma_0"] =  1.99;
  test_params["S_alpha"] =  1.339;
  test_params["rho_0"] =  1500;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 43513494000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 27600000000);

  // Get the moduli upper bound at zero pressure
  try {
    moduli = model.getElasticModuliUpperBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 1.490529443465273e+19);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 32733042374.41661);

  try {
  moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 796886.9730089129);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 1500.0;
  state.meltingTemp = 700.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 31022404424.38447);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 24177604424.38446);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 15621604424.38446);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 31022402053.61735);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 24177602053.61735);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 15621602053.61735);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 31022400953.20474);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 24177600953.20474);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 15621600953.20474);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 31022395575.61554);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 24177595575.61554);
   p_KG[keys[ii++]] = std::make_pair(57162275.46470965, 15621595575.61554);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 31022397946.38265);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 24177597946.38265);
   p_KG[keys[ii++]] = std::make_pair(12394939722.54329, 15621597946.38265);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 31022399046.79527);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 24177599046.79526);
   p_KG[keys[ii++]] = std::make_pair(2.235794165197909e+19, 15621599046.79527);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

TEST(ElasticModuliMetalIsoTest, mieEnergyBulkShearSCG)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "metal_iso");
  xmlDocSetRootElement(doc, rootNode);

  // For EOS model
  auto eosNode = xmlNewChild(rootNode, nullptr, BAD_CAST "equation_of_state", BAD_CAST "");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "mie_gruneisen_energy");
  xmlNewChild(eosNode, nullptr, BAD_CAST "C_0", BAD_CAST "1680");
  xmlNewChild(eosNode, nullptr, BAD_CAST "Gamma_0", BAD_CAST "0.59");
  xmlNewChild(eosNode, nullptr, BAD_CAST "S_alpha", BAD_CAST "1.123");
  xmlNewChild(eosNode, nullptr, BAD_CAST "S_2", BAD_CAST "3.98");
  xmlNewChild(eosNode, nullptr, BAD_CAST "S_3", BAD_CAST "-5.8");
  xmlNewChild(eosNode, nullptr, BAD_CAST "rho_0", BAD_CAST "2200");

  // For shear model
  auto shearNode = xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus_model", BAD_CAST "");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "scg_shear");
  xmlNewChild(shearNode, nullptr, BAD_CAST "mu_0", BAD_CAST "27.6e9");
  xmlNewChild(shearNode, nullptr, BAD_CAST "A", BAD_CAST "65.0e-12");
  xmlNewChild(shearNode, nullptr, BAD_CAST "B", BAD_CAST "0.62e-3");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Check whether the combination fails
  try {
    ElasticModuli_MetalIso model(ps);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }

  ElasticModuli_MetalIso model(ps);

  // Get initial parameters
  std::map<std::string, double> test_params;
  test_params["mu_0"] =  27.6e9;
  test_params["A"] =  65.0e-12;
  test_params["B"] =  0.62e-3;
  test_params["C_0"] =  1680;
  test_params["Gamma_0"] =  0.59;
  test_params["S_alpha"] =  1.123;
  test_params["S_2"] =  3.98;
  test_params["S_3"] =  -5.8;
  test_params["rho_0"] =  2200;

  std::map<std::string, double> params = model.getParameters();
  for (const auto& param : params) {
    //std::cout << "params = " << param.first  << " : " << param.second << "\n";
    ASSERT_DOUBLE_EQ(test_params[param.first], param.second);
  }

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 6209280000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 27600000000);

  // Get the moduli upper bound at zero pressure
  try {
    moduli = model.getElasticModuliUpperBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 94299588488.5387);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 32733042374.41661);

  try {
    moduli = model.getElasticModuliLowerBound();
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << "\n";
  }
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_DOUBLE_EQ(moduli.bulkModulus, 2822400000);
  ASSERT_DOUBLE_EQ(moduli.shearModulus, 1.0e-6);

  // Create a model state
  ModelStateBase state;
  state.initialDensity = 2200.0;
  state.meltingTemp = 700.0;
  state.energy = 0.0;

  // Set up pressures, densities, temperatures
  std::vector<double> pressures = { -1000, 0, 1000 };
  std::vector<double> densities = { 100, 1000, 10000 }; // not consistent with pressures
  std::vector<double> temperatures = { 100, 500, 1000 };
  std::vector<std::tuple<double, double, double>> keys;
  for (double p : pressures) {
    for (double rho : densities) {
      for (double T : temperatures) {
        keys.push_back(std::make_tuple(p, rho, T));
      }
    }
  }
  std::map<std::tuple<double, double, double>, std::pair<double, double>> p_KG;
  auto ii = 0u;
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022405026.85857);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177605026.85856);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621605026.85856);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022402333.26105);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177602333.26105);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621602333.26106);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 31022401083.00385);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 24177601083.00385);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 15621601083.00385);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 31022400000);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 24177600000);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 15621600000);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022394973.14144);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177594973.14144);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621594973.14144);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 31022397666.73894);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 24177597666.73895);
   p_KG[keys[ii++]] = std::make_pair(6209280000, 15621597666.73895);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 31022398916.99615);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 24177598916.99615);
   p_KG[keys[ii++]] = std::make_pair(207459094674.7852, 15621598916.99615);

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
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].first, moduli.bulkModulus);
        ASSERT_DOUBLE_EQ(p_KG[std::make_tuple(p, rho, T)].second, moduli.shearModulus);
      }
    }
  }
}

