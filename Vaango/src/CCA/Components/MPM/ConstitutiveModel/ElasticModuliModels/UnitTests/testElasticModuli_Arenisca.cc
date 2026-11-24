#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arenisca.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arenisca3.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <gtest/gtest.h>
#include <cmath>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;

TEST(ElasticModuliAreniscaTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "arenisca");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "B0", BAD_CAST "1.0e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "B1", BAD_CAST "2.0e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "B2", BAD_CAST "1.0e6");
  xmlNewChild(rootNode, nullptr, BAD_CAST "B3", BAD_CAST "0.5e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "B4", BAD_CAST "0.1");
  xmlNewChild(rootNode, nullptr, BAD_CAST "B01", BAD_CAST "0.01");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", BAD_CAST "0.6e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G1", BAD_CAST "0.2");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G2", BAD_CAST "0.1");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G3", BAD_CAST "0.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G4", BAD_CAST "0.0");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  ElasticModuli_Arenisca model(ps);
  
  // Initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, 1.0e8);
  EXPECT_DOUBLE_EQ(moduli.shearModulus, 0.6e8);

  // Bounds
  moduli = model.getElasticModuliLowerBound();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, 1.0e8);
  
  moduli = model.getElasticModuliUpperBound();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, 1.0e8 + 2.0e8);

  // Current moduli
  ModelState_Arenisca3 state;
  state.I1 = -1.0e7; // Compressive
  state.plasticStrainTensor = Uintah::Matrix3(0,0,0,0,0,0,0,0,0);
  
  moduli = model.getCurrentElasticModuli(&state);
  
  double I1_bar = 1.0e7;
  double exp_B2_by_I1 = std::exp(-1.0e6 / I1_bar);
  double expected_K = 1.0e8 + 0.01 * I1_bar + 2.0e8 * exp_B2_by_I1;
  double expected_nu = 0.2 + 0.1 * exp_B2_by_I1;
  double expected_G = 1.5 * expected_K * (1.0 - 2.0 * expected_nu) / (1.0 + expected_nu);
  
  EXPECT_NEAR(moduli.bulkModulus, expected_K, 1.0e-5);
  EXPECT_NEAR(moduli.shearModulus, expected_G, 1.0e-5);

  // With plastic strain
  state.plasticStrainTensor = Uintah::Matrix3(-0.01, 0, 0, 0, -0.01, 0, 0, 0, -0.01);
  double ev_p_bar = 0.03;
  // expected_G does not change because it is calculated before KK is updated with ev_p in the model
  expected_K -= 0.5e8 * std::exp(-0.1 / ev_p_bar);
  
  moduli = model.getCurrentElasticModuli(&state);
  EXPECT_NEAR(moduli.bulkModulus, expected_K, 1.0e-5);
  EXPECT_NEAR(moduli.shearModulus, expected_G, 1.0e-5);

  xmlFreeDoc(doc);
}
