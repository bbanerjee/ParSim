#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/BorjaEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
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

TEST(BorjaEOSTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "equation_of_state");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "borja_pressure");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  double p0 = -1.0e4;
  double alpha = 0.5;
  double kappatilde = 0.02;
  double epse_v0 = 0.0;
  xmlNewChild(rootNode, nullptr, BAD_CAST "p0", BAD_CAST "-1.0e4");
  xmlNewChild(rootNode, nullptr, BAD_CAST "alpha", BAD_CAST "0.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "kappatilde", BAD_CAST "0.02");
  xmlNewChild(rootNode, nullptr, BAD_CAST "epse_v0", BAD_CAST "0.0");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  BorjaEOS model(ps);
  
  double kappahat = kappatilde / (1.0 - kappatilde);

  // Pressure (rho_orig, rho_cur)
  double rho_orig = 1000.0;
  double rho_cur = 1100.0;
  double p = model.computePressure(rho_orig, rho_cur);
  
  double epse_v = rho_orig / rho_cur - 1.0;
  double expected_p = p0 * std::exp(-epse_v / kappahat);
  EXPECT_NEAR(p, expected_p, 1.0e-7);

  // Bulk Modulus
  double K = model.computeBulkModulus(rho_orig, rho_cur);
  double expected_K = -expected_p / kappahat;
  EXPECT_NEAR(K, expected_K, 1.0e-7);

  // Strain Energy
  double U = model.computeStrainEnergy(rho_orig, rho_cur);
  double expected_U = -p0 * kappahat * std::exp(-epse_v / kappahat);
  EXPECT_NEAR(U, expected_U, 1.0e-7);

  // Density
  double rho_calc = model.computeDensity(rho_orig, expected_p);
  EXPECT_NEAR(rho_calc, rho_cur, 1.0e-7);

  xmlFreeDoc(doc);
}
