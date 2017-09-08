#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;
using Uintah::ProblemSetupException;
using Uintah::InvalidValue;
using nlohmann::json;

TEST(ElasticModuliTabularTest, constructorTest)
{
  char currPath[2000];
  if (!getcwd(currPath, sizeof(currPath))) {
    std::cout << "Current path not found\n";
  }
  std::string json_file = std::string(currPath) + "/" + "table_elastic.json";

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "table_elastic.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "PlasticStrainVol, TotalStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");
  xmlNewChild(rootNode, nullptr, BAD_CAST "G0", 
              BAD_CAST "1.0e4");
  xmlNewChild(rootNode, nullptr, BAD_CAST "nu", 
              BAD_CAST "0.2");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a model
  ElasticModuli_Tabular model(ps);
  try {
    ElasticModuli moduli = model.getInitialElasticModuli();
    EXPECT_DOUBLE_EQ(moduli.bulkModulus, 1.0e4);
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  // Copy
  ElasticModuli_Tabular modelCopy(&model);
  try {
    ElasticModuli moduli = modelCopy.getInitialElasticModuli();
    EXPECT_DOUBLE_EQ(moduli.bulkModulus, 1.0e4);
    EXPECT_DOUBLE_EQ(moduli.shearModulus, 7500);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }

  // Modelstate test
  ModelState_Tabular state;
  state.elasticStrainTensor = Uintah::Matrix3(-0.03, 0, 0, 0, -0.03, 0, 0, 0, -0.03);
  state.plasticStrainTensor = Uintah::Matrix3(-0.02, 0, 0, 0, -0.02, 0, 0, 0, -0.02);
  try {
    ElasticModuli moduli = model.getCurrentElasticModuli(&state);
    EXPECT_NEAR(moduli.bulkModulus, 9040, 1.0e-7);
    EXPECT_NEAR(moduli.shearModulus, 6780, 1.0e-7);
    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }
}
