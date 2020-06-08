#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>

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

TEST(ElasticModuliArenaTest, computeTests)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "arena");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "b0", BAD_CAST "0.00290858181781614");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b1", BAD_CAST "0.47312420784766");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b2", BAD_CAST "1.50567779375549");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b3", BAD_CAST "2.57284042409447");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b4", BAD_CAST "2.07992105987609");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", BAD_CAST "1.0e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu1", BAD_CAST "0.35");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu2", BAD_CAST "-0.35");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  /*
  // Run thru the ProblemSpec machinery
  xmlNode* d_node;
  d_node = const_cast<xmlNode*>(xmlDocGetRootElement(doc));
  if (d_node == 0) {
    std::cout << "**Error** Could not create node." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }
  const xmlNode* d_child = d_node->children;
  while (d_child != 0) {
    std::string child_name((const char *)(d_child->name));
    std::cout << "child_name = " << child_name << std::endl;
    d_child = d_child->next;
  }
  */

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors
  ElasticModuli_Arena model(ps);
  ElasticModuli_Arena model_copy(model);

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  //std::cout << std::setprecision(16) 
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_NEAR(moduli.bulkModulus, 116343272.7126456, 1.0e-6);
  ASSERT_NEAR(moduli.shearModulus, 173983253.4171182, 1.0e-6);

  // Get the moduli upper bound at zero pressure
  moduli = model.getElasticModuliUpperBound();
  //std::cout << std::setprecision(16)
  //          << " K = " << moduli.bulkModulus << " G = " << moduli.shearModulus
  //          << std::endl;
  ASSERT_NEAR(moduli.bulkModulus, std::numeric_limits<double>::max(), 1.0e-6);
  ASSERT_NEAR(moduli.shearModulus, std::numeric_limits<double>::max(), 1.0e-6);

  // Set up list of pressures in log10 scale
  std::vector<double> pressures = { -2, 0, 6, 8, 11 };

  // Set up volumetric plastic strains
  std::vector<double> volPlasticStrains = { -1.0, -0.1, 0.0, 0.1, 0.3, 1.0 };
  std::vector<double> porosities = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
  std::vector<double> saturations = { 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 };

  // Set up plastic strain matrices
  std::vector<Matrix3> plasticStrains;
  for (double ep_v : volPlasticStrains) {
    double ep = -ep_v / 3.0;
    plasticStrains.emplace_back(
      Matrix3(ep, 0.0, 0.0, 0.0, ep, 0.0, 0.0, 0.0, ep));
  }

  // Drained (dry) sand
  ModelState_Arena state;
  state.porosity = porosities[4];
  state.saturation = saturations[0];

  std::map<double, std::pair<double, double>> I1_KG;
  I1_KG[-0.01] = std::make_pair(116348059.2816119, 173990384.2029271);
  I1_KG[-1] = std::make_pair(116407266.5941095, 174078654.5535643);
  I1_KG[-1e+06] = std::make_pair(201335456.1999204, 300413920.3418134);
  I1_KG[-1e+08] = std::make_pair(1046521964.819405, 1527736301.547844);
  I1_KG[-1e+11] = std::make_pair(74552697327.11647, 75257479257.76637);

  for (double pp : pressures) {
    // ** Compression tests **
    state.plasticStrainTensor = plasticStrains[0];
    state.I1_eff = -std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  I1_KG.clear();
  I1_KG[0.01] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+06] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+08] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+11] = std::make_pair(116343272.7126456, 173983253.4171182);
  for (double pp : pressures) {
    // ** Tension tests **
    state.plasticStrainTensor = plasticStrains[5];
    state.I1_eff = std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  // Partially saturated sand
  state.porosity = porosities[4];
  state.saturation = saturations[3];

  I1_KG.clear();
  I1_KG[-0.01] = std::make_pair(117053161.3779722, 173990384.2029271);
  I1_KG[-1] = std::make_pair(117112366.5970196, 174078654.5535643);
  I1_KG[-1e+06] = std::make_pair(202037556.7507421, 300413920.3418134);
  I1_KG[-1e+08] = std::make_pair(1047194681.441464, 1527736301.547844);
  I1_KG[-1e+11] = std::make_pair(74552927664.48474, 75257479257.76637);
  for (double pp : pressures) {
    // ** Compression tests **
    state.plasticStrainTensor = plasticStrains[1];
    state.I1_eff = -std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  I1_KG.clear();
  I1_KG[0.01] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+06] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+08] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+11] = std::make_pair(116343272.7126456, 173983253.4171182);
  for (double pp : pressures) {
    // ** Tension tests **
    state.plasticStrainTensor = plasticStrains[4];
    state.I1_eff = std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  // Fully saturated sand
  state.porosity = porosities[4];
  state.saturation = saturations[6];

  I1_KG.clear();
  I1_KG[-0.01] = std::make_pair(5190747007.780967, 173990384.2029271);
  I1_KG[-1] = std::make_pair(5190792107.501109, 174078654.5535643);
  I1_KG[-1e+06] = std::make_pair(5255515989.859488, 300413920.3418134);
  I1_KG[-1e+08] = std::make_pair(5903390471.241636, 1527736301.547844);
  I1_KG[-1e+11] = std::make_pair(76337395480.34607, 75257479257.76637);
  for (double pp : pressures) {
    // ** Compression tests **
    state.plasticStrainTensor = plasticStrains[1];
    state.I1_eff = -std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  I1_KG.clear();
  I1_KG[0.01] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+06] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+08] = std::make_pair(116343272.7126456, 173983253.4171182);
  I1_KG[1e+11] = std::make_pair(116343272.7126456, 173983253.4171182);
  for (double pp : pressures) {
    // ** Tension tests **
    state.plasticStrainTensor = plasticStrains[4];
    state.I1_eff = std::pow(10, pp);
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(I1_KG[state.I1_eff].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(I1_KG[state.I1_eff].second, moduli.shearModulus, 1.0e-3);
  }

  // Vary saturation
  //std::cout << "Varying saturation" << std::endl;
  state.porosity = porosities[4];
  state.plasticStrainTensor = plasticStrains[1];
  state.I1_eff = -std::pow(10, pressures[4]);

  std::map<double, std::pair<double, double>> Sw_KG;
  Sw_KG[0] = std::make_pair(74552697327.11647, 75257479257.76637);
  Sw_KG[0.1] = std::make_pair(74552825299.66171, 75257479257.76637);
  Sw_KG[0.3] = std::make_pair(74552861859.87523, 75257479257.76637);
  Sw_KG[0.5] = std::make_pair(74552927664.48474, 75257479257.76637);
  Sw_KG[0.7] = std::make_pair(74553081189.70212, 75257479257.76637);
  Sw_KG[0.9] = std::make_pair(74553848419.70689, 75257479257.76637);
  Sw_KG[1.0] = std::make_pair(76337395480.34607, 75257479257.76637);
  for (double sw : saturations) {
    state.saturation = sw;
    moduli = model.getCurrentElasticModuli(&state);
    //std::cout << std::setprecision(16)
    //          << "ep_v = " << state.plasticStrainTensor.Trace()
    //          << " phi = " << state.porosity << " S_w = " << state.saturation
    //          << " I1_eff = " << state.I1_eff << " K = " << moduli.bulkModulus
    //          << " G = " << moduli.shearModulus << std::endl;
    ASSERT_NEAR(Sw_KG[sw].first, moduli.bulkModulus, 1.0e-3);
    ASSERT_NEAR(Sw_KG[sw].second, moduli.shearModulus, 1.0e-3);
  }

  /*
  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  std::cout << "params[Ks] = " << params["Ks"]  << " Pa" << std::endl;

  // Set up densities
  double rho_orig  = 1.0;
  std::vector<double> rho_cur;
  for (int ii = -2; ii < 10; ii++) {
    double rho = rho_orig*(1.0 + (double) ii/10.0);
    rho_cur.emplace_back(rho);
    std::cout << "rho = " << rho << std::endl;
  }

  // Compute the pressure
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);
    std::cout << "p = " << pp << std::endl;
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (int ii = -2; ii < 12; ii++) {
    pressures.emplace_back((double) ii);
  }

  // Convert to Pa and Compute bulk modulus
  std::cout << "Compression:" << std::endl;
  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  ModelState_Arenisca3 state;
  for (double pp : pressures) {
    state.I1_eff = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1_eff = " << state.I1_eff << " K = " << K << std::endl;
  }

  // Test tension states
  // Convert to Pa and Compute bulk modulus
  std::cout << "Tension:" << std::endl;
  for (double pp : pressures) {
    pp = -std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  for (double pp : pressures) {
    state.I1_eff = std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1_eff = " << state.I1_eff << " K = " << K << std::endl;
  }
  */

  // Free the document (*WARNING** `doc` seems to be getting freed elsewhere)
  // xmlFreeDoc(doc);
}
