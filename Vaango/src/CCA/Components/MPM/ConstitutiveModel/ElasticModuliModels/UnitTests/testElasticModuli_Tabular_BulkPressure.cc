#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_BulkPressure.h>
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

TEST(ElasticModuliTabularBulkPressureTest, basicTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string table_file = std::string(currPath) + "/table_bulk_pressure.json";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_bulk_pressure");
  xmlDocSetRootElement(doc, rootNode);

  // Create child nodes for TabularData
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST table_file.c_str());
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", BAD_CAST "PlasticStrainVol, Pressure");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", BAD_CAST "BulkModulus");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation", nullptr);
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Create child nodes for ElasticModuli_Tabular_BulkPressure
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", BAD_CAST "0.6e9");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu", BAD_CAST "0.2");

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    FAIL() << "**Error** Could not create ProblemSpec.";
  }

  // Constructors
  ElasticModuli_Tabular_BulkPressure model(ps);
  
  // Initial moduli
  ElasticModuli moduli = model.getInitialElasticModuli();
  EXPECT_NEAR(moduli.bulkModulus, 1.0e9, 1.0e5);
  EXPECT_GT(moduli.shearModulus, 0.0);

  // Current moduli
  ModelState_Tabular state;
  state.stressTensor = Uintah::Matrix3(-2.0e5, 0, 0, 0, -2.0e5, 0, 0, 0, -2.0e5); // p_bar = 2.0e5
  state.plasticStrainTensor = Uintah::Matrix3(-0.1, 0, 0, 0, -0.1, 0, 0, 0, -0.1);   // ev_p = 0.3
  
  // ev_p = 0.3 is between 0 and 0.5. weight = 0.6
  // p_bar = 2.0e5 is between 0 and 1.0e6. weight = 0.2
  // At ev_p = 0: K = 1.0e9 + (1.1e9 - 1.0e9) * 0.2 = 1.02e9
  // At ev_p = 0.5: K = 1.2e9 + (1.3e9 - 1.2e9) * 0.2 = 1.22e9
  // K_interp = 1.02e9 + (1.22e9 - 1.02e9) * 0.6 = 1.02e9 + 0.2e9 * 0.6 = 1.14e9

  moduli = model.getCurrentElasticModuli(&state);
  EXPECT_NEAR(moduli.bulkModulus, 1.14e9, 1.0e-7);

  xmlFreeDoc(doc);
}
