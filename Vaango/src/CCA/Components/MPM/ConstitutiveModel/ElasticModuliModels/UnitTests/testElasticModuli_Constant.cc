#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Constant.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;

TEST(ElasticModuliConstantTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "constant");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  double bulk = 1.0e9;
  double shear = 0.6e9;
  xmlNewChild(rootNode, nullptr, BAD_CAST "bulk_modulus", BAD_CAST "1.0e9");
  xmlNewChild(rootNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "0.6e9");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  ElasticModuli_Constant model(ps);
  ElasticModuli_Constant model_copy(&model);

  // Get the initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, bulk);
  EXPECT_DOUBLE_EQ(moduli.shearModulus, shear);

  // Get current moduli
  moduli = model.getCurrentElasticModuli(nullptr);
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, bulk);
  EXPECT_DOUBLE_EQ(moduli.shearModulus, shear);

  // Bounds
  moduli = model.getElasticModuliLowerBound();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, bulk);
  EXPECT_DOUBLE_EQ(moduli.shearModulus, shear);

  moduli = model.getElasticModuliUpperBound();
  EXPECT_DOUBLE_EQ(moduli.bulkModulus, bulk);
  EXPECT_DOUBLE_EQ(moduli.shearModulus, shear);

  // Parameters
  auto params = model.getParameters();
  EXPECT_DOUBLE_EQ(params["K"], bulk);
  EXPECT_DOUBLE_EQ(params["G"], shear);

  // Derivatives
  auto derivs = model.computeDModuliDIntVar(nullptr);
  EXPECT_TRUE(derivs.empty());

  auto moduli_derivs = model.computeModuliAndDModuliDIntVar(nullptr);
  EXPECT_DOUBLE_EQ(moduli_derivs.first.bulkModulus, bulk);
  EXPECT_TRUE(moduli_derivs.second.empty());

  xmlFreeDoc(doc);
}
