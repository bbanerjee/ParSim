#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/HyperElasticEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
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

TEST(HyperElasticEOSTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "equation_of_state");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "default_hyper");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  double bulk = 1.0e9;
  xmlNewChild(rootNode, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e9");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  HyperElasticEOS model(ps);
  HyperElasticEOS model_copy(&model);

  // Parameters
  auto params = model.getParameters();
  EXPECT_DOUBLE_EQ(params["bulk_modulus"], bulk);

  // Pressure (rho_orig, rho_cur)
  double rho_orig = 1000.0;
  double rho_cur = 1100.0;
  double p = model.computePressure(rho_orig, rho_cur);
  
  double J = rho_orig / rho_cur;
  double expected_p = 0.5 * bulk * (J - 1.0 / J);
  EXPECT_DOUBLE_EQ(p, expected_p);

  // Bulk Modulus
  double K = model.computeBulkModulus(rho_orig, rho_cur);
  double expected_K = 0.5 * bulk * (1.0 + 1.0 / (J * J));
  EXPECT_DOUBLE_EQ(K, expected_K);

  // Strain Energy
  double U = model.computeStrainEnergy(rho_orig, rho_cur);
  double expected_U = 0.5 * bulk * (0.5 * (J * J - 1.0) - std::log(J));
  EXPECT_NEAR(U, expected_U, 1.0e-7);

  // Density
  double rho_calc = model.computeDensity(rho_orig, expected_p);
  EXPECT_NEAR(rho_calc, rho_cur, 1.0e-7);

  xmlFreeDoc(doc);
}
