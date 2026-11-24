#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MieGruneisenEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
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

TEST(MieGruneisenEOSTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "equation_of_state");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "mie_gruneisen");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  double C_0 = 3500.0;
  double Gamma_0 = 1.5;
  double S_alpha = 1.2;
  double rho_0 = 8900.0;
  xmlNewChild(rootNode, nullptr, BAD_CAST "C_0", BAD_CAST "3500.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "Gamma_0", BAD_CAST "1.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "S_alpha", BAD_CAST "1.2");
  xmlNewChild(rootNode, nullptr, BAD_CAST "rho_0", BAD_CAST "8900.0");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  MieGruneisenEOS model(ps);
  MieGruneisenEOS model_copy(&model);

  // Parameters
  auto params = model.getParameters();
  EXPECT_DOUBLE_EQ(params["C_0"], C_0);
  EXPECT_DOUBLE_EQ(params["Gamma_0"], Gamma_0);
  EXPECT_DOUBLE_EQ(params["S_alpha"], S_alpha);
  EXPECT_DOUBLE_EQ(params["rho_0"], rho_0);

  // Pressure (rho_orig, rho_cur)
  double rho_cur = 9000.0;
  double p = model.computePressure(rho_0, rho_cur);
  
  double J = rho_0 / rho_cur;
  double J_one = J - 1.0;
  double expected_p = rho_0 * C_0 * C_0 * J_one * (2.0 + Gamma_0 * J_one) / (2.0 * std::pow(1.0 + J_one * S_alpha, 2));
  EXPECT_NEAR(p, expected_p, 1.0e-7);

  // Bulk Modulus
  double K = model.computeBulkModulus(rho_0, rho_cur);
  double denom = 1.0 + J_one * S_alpha;
  double dp_dJ = (C_0 * C_0 * rho_0 * (1.0 + J_one * (Gamma_0 - S_alpha))) / (denom * denom * denom);
  double expected_K = J * dp_dJ;
  EXPECT_NEAR(K, expected_K, 1.0e-7);

  // Strain Energy
  double U = model.computeStrainEnergy(rho_0, rho_cur);
  double expected_U = rho_0 * C_0 * C_0 / (2.0 * std::pow(S_alpha, 3)) * 
                      ((J_one * S_alpha * (-2.0 * S_alpha + Gamma_0 * (2.0 + J_one * S_alpha))) / (1.0 + J_one * S_alpha) +
                       2.0 * (S_alpha - Gamma_0) * std::log(1.0 + J_one * S_alpha));
  EXPECT_NEAR(U, expected_U, 1.0e-7);

  // Density
  double rho_calc = model.computeDensity(rho_0, expected_p);
  EXPECT_NEAR(rho_calc, rho_cur, 1.0e-2);

  xmlFreeDoc(doc);
}
