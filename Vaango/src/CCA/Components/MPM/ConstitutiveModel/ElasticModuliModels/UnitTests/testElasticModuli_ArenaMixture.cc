#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_ArenaMixture.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
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

TEST(ElasticModuliArenaMixtureTest, basicTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "arena_mixture");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "vol_frac.phase1", BAD_CAST "0.5");

  xmlNewChild(rootNode, nullptr, BAD_CAST "b0.phase1", BAD_CAST "0.003");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b1.phase1", BAD_CAST "0.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b2.phase1", BAD_CAST "1.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b3.phase1", BAD_CAST "2.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b4.phase1", BAD_CAST "2.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0.phase1", BAD_CAST "1.0e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu1.phase1", BAD_CAST "0.35");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu2.phase1", BAD_CAST "-0.35");

  xmlNewChild(rootNode, nullptr, BAD_CAST "b0.phase2", BAD_CAST "0.003");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b1.phase2", BAD_CAST "0.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b2.phase2", BAD_CAST "1.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b3.phase2", BAD_CAST "2.5");
  xmlNewChild(rootNode, nullptr, BAD_CAST "b4.phase2", BAD_CAST "2.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0.phase2", BAD_CAST "1.0e8");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu1.phase2", BAD_CAST "0.35");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu2.phase2", BAD_CAST "-0.35");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  ElasticModuli_ArenaMixture model(ps);
  
  // Initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  // Standard values for granite bulk modulus are around 30 GPa
  EXPECT_GT(moduli.bulkModulus, 0.0);
  EXPECT_GT(moduli.shearModulus, 0.0);

  // Current moduli
  ModelState_Arena state;
  state.I1_eff = -1.0e7; // Compressive
  state.plasticStrainTensor = Uintah::Matrix3(0,0,0,0,0,0,0,0,0);
  state.porosity = 0.4;
  state.saturation = 0.0;
  
  moduli = model.getCurrentElasticModuli(&state);
  EXPECT_GT(moduli.bulkModulus, 0.0);
  EXPECT_GT(moduli.shearModulus, 0.0);

  // With saturation
  state.saturation = 0.5;
  state.pbar_w = 1.0e6;
  moduli = model.getCurrentElasticModuli(&state);
  EXPECT_GT(moduli.bulkModulus, 0.0);
  EXPECT_GT(moduli.shearModulus, 0.0);

  xmlFreeDoc(doc);
}
