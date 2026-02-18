#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_Bulk.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <gtest/gtest.h>
#include <fstream>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;

TEST(ElasticModuliTabularBulkTest, basicTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string table_file = std::string(currPath) + "/table_bulk.json";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_bulk");
  xmlDocSetRootElement(doc, rootNode);

  // Create child nodes for TabularData
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST table_file.c_str());
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", BAD_CAST "PlasticStrainVol, ElasticStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", BAD_CAST "BulkModulus");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation", nullptr);
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Create child nodes for ElasticModuli_Tabular_Bulk
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", BAD_CAST "0.6e9");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu", BAD_CAST "0.2");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  ElasticModuli_Tabular_Bulk model(ps);
  
  // Initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  // K at ev_e=1e-6, ev_p=0 should be around 1.0e9 (interpolated)
  EXPECT_NEAR(moduli.bulkModulus, 1.0e9, 1.0e5);
  EXPECT_GT(moduli.shearModulus, 0.0);

  // Current moduli
  ModelState_Tabular state;
  state.elasticStrainTensor = Uintah::Matrix3(-0.05, 0, 0, 0, -0.05, 0, 0, 0, -0.05); // ev_e = 0.15
  state.plasticStrainTensor = Uintah::Matrix3(-0.1, 0, 0, 0, -0.1, 0, 0, 0, -0.1);   // ev_p = 0.3
  
  // Hand interpolation:
  // ev_p = 0.3 is between 0 and 0.5. weight = (0.3 - 0) / 0.5 = 0.6
  // At ev_p = 0: K = 1.0e9 + (1.1e9 - 1.0e9) * (0.15 - 0) / 0.1 = 1.15e9
  // At ev_p = 0.5: K = 1.2e9 + (1.3e9 - 1.2e9) * (0.15 - 0) / 0.1 = 1.35e9
  // K_interp = 1.15e9 + (1.35e9 - 1.15e9) * 0.6 = 1.15e9 + 0.2e9 * 0.6 = 1.27e9
  
  // Wait, my table only goes to ev_e = 0.1. ev_e = 0.15 is out of bounds.
  // TabularData::interpolate usually throws or clamps.
  
  state.elasticStrainTensor = Uintah::Matrix3(-0.02, 0, 0, 0, -0.02, 0, 0, 0, -0.02); // ev_e = 0.06
  // At ev_p = 0: K = 1.0e9 + (1.1e9 - 1.0e9) * (0.06 - 0) / 0.1 = 1.06e9
  // At ev_p = 0.5: K = 1.2e9 + (1.3e9 - 1.2e9) * (0.06 - 0) / 0.1 = 1.26e9
  // K_interp = 1.06e9 + (1.26e9 - 1.06e9) * 0.6 = 1.06e9 + 0.2e9 * 0.6 = 1.18e9

  moduli = model.getCurrentElasticModuli(&state);
  EXPECT_NEAR(moduli.bulkModulus, 1.18e9, 1.0e-7);

  xmlFreeDoc(doc);
}
