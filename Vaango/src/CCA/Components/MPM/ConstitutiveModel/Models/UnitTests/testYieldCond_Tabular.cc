#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCond_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xmlstring.h>

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

class YieldCondTabularTest : public ::testing::Test {

protected:

  static void SetUpTestCase() {
    char currPath[2000];
    if (!getcwd(currPath, sizeof(currPath))) {
      std::cout << "Current path not found\n";
    }
    std::string json_file = std::string(currPath) + "/" + "table_yield.json";

    // Create a new document
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

    // Create root node
    xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "plastic_yield_condition");
    xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular");
    xmlDocSetRootElement(doc, rootNode);

    // Create a child node
    xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
                BAD_CAST "table_yield.json");
    xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "Pressure");
    xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "SqrtJ2");
    auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                              BAD_CAST "");
    xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

    // Print the document to stdout
    xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    if (!ps) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }

    // Update the json file and create a circular yield function
    auto doc1 = xmlCopyDoc(doc, 1);
    auto node = xmlDocGetRootElement(doc1);
    auto xpathCtx = xmlXPathNewContext(doc1);
    const xmlChar* xpathExpr = xmlStrncatNew(BAD_CAST ".//", BAD_CAST "filename", -1);
    auto xpathObj = xmlXPathNodeEval(node, xpathExpr, xpathCtx);
    xmlXPathFreeContext(xpathCtx);
    auto jsonNode = xpathObj->nodesetval->nodeTab[0];

    xmlNodeSetContent(jsonNode, BAD_CAST "table_yield_circle.json");
    xmlSaveFormatFileEnc("-", doc1, "ISO-8859-1", 1);

    ps_circle = scinew ProblemSpec(xmlDocGetRootElement(doc1), false);
    if (!ps_circle) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }
  }

  static void TearDownTestCase() {}

  static ProblemSpecP ps;
  static ProblemSpecP ps_circle;
};

ProblemSpecP YieldCondTabularTest::ps = nullptr;
ProblemSpecP YieldCondTabularTest::ps_circle = nullptr;

TEST_F(YieldCondTabularTest, constructorTest)
{
  // Create a model
  YieldCond_Tabular model(ps);
  //std::cout << model;

  // Copy
  YieldCond_Tabular modelCopy(&model);
  //std::cout << modelCopy;

  YieldCond_Tabular model_circle(ps_circle);
  //std::cout << model_circle;

  auto params = model.getParameters();
  ASSERT_DOUBLE_EQ(params.at("I1_min"), -19200);
  ASSERT_DOUBLE_EQ(params.at("I1_max"), 30);
  ASSERT_DOUBLE_EQ(params.at("sqrtJ2_max"), 900);
}

TEST_F(YieldCondTabularTest, evalYieldCondition)
{
  YieldCond_Tabular model(ps);
  ModelState_Tabular state;

  state.I1 = 300*3; // Tension
  state.sqrt_J2 = 1000;
  EXPECT_EQ(model.evalYieldCondition(&state), 1);

  state.I1 = 2*3;  // Tension
  state.sqrt_J2 = 39;
  EXPECT_EQ(model.evalYieldCondition(&state), -1);

  state.I1 = -1000*3;  // Compression
  state.sqrt_J2 = 625;
  EXPECT_EQ(model.evalYieldCondition(&state), -1);

  state.I1 = -1000*3;  // Compression
  state.sqrt_J2 = 605;
  EXPECT_EQ(model.evalYieldCondition(&state), -1);

  state.I1 = -1000*3;  // Compression
  state.sqrt_J2 = 635;
  EXPECT_EQ(model.evalYieldCondition(&state), 1);

  state.I1 = -7000*3;  // Compression
  state.sqrt_J2 = 1000;
  EXPECT_THROW(model.evalYieldCondition(&state), Uintah::InvalidValue);

  EXPECT_EQ(model.evalYieldConditionMax(&state), 900);
}

TEST_F(YieldCondTabularTest, eval_df_dsigma)
{
  YieldCond_Tabular model(ps);
  ModelState_Tabular state;
  Matrix3 zero(0.0);
  Matrix3 df_dsigma(0.0);
  //std::cout << model;

  // Zero everything (elastic)
  model.eval_df_dsigma(zero, &state, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Tension (zero stress)
  state.I1 = 300*3; 
  state.sqrt_J2 = 1000;
  model.eval_df_dsigma(zero, &state, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Circular yield function
  YieldCond_Tabular model_circle(ps_circle);
  ModelState_Tabular state_circle;

  state_circle.stressTensor = Matrix3(2, 4, 0, 4, 2, 0, 0, 0, 2);
  state_circle.updateStressInvariants();
  model_circle.eval_df_dsigma(zero, &state_circle, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), 0.22177, 1.0e-5); // Exact value = 0.20569
  EXPECT_NEAR(df_dsigma(0,1), 0.65286, 1.0e-5); // Exact value = 0.66071

  state_circle.stressTensor = Matrix3(-2, 4, 0, 4, -2, 0, 0, 0, -2);
  state_circle.updateStressInvariants();
  model_circle.eval_df_dsigma(zero, &state_circle, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), -0.22177, 1.0e-5); // Exact value = -0.20569
  EXPECT_NEAR(df_dsigma(0,1), 0.65286, 1.0e-5); // Exact value = 0.66071

  state_circle.stressTensor = Matrix3(3, 0, 0, 0, 3, 0, 0, 0, 3);
  state_circle.updateStressInvariants();
  model_circle.eval_df_dsigma(zero, &state_circle, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  state_circle.stressTensor = Matrix3(3, 1, 0, 1, 3, 0, 0, 0, 3);
  state_circle.updateStressInvariants();
  model_circle.eval_df_dsigma(zero, &state_circle, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), 0.53852, 1.0e-5); // Exact value = 0.53644
  EXPECT_NEAR(df_dsigma(0,1), 0.25494, 1.0e-5); // Exact value = 0.26145

  state_circle.stressTensor = Matrix3(-3, 1, 0, 1, -3, 0, 0, 0, -3);
  state_circle.updateStressInvariants();
  model_circle.eval_df_dsigma(zero, &state_circle, df_dsigma);
  EXPECT_NEAR(df_dsigma(0,0), -0.53852, 1.0e-5); // Exact value = -0.53644
  EXPECT_NEAR(df_dsigma(0,1), 0.25494, 1.0e-5); // Exact value = 0.26145

  /*
  try {
    state_circle.stressTensor = Matrix3(3, 0, 0, 0, 3, 0, 0, 0, 3);
    std::cout << "df_dsigma = " << df_dsigma << std::endl;
  } catch (Uintah::InvalidValue e)  {
    std::cout <<  e.message() << std::endl;
  }
  */
}

TEST_F(YieldCondTabularTest, getClosestPoint)
{
  YieldCond_Tabular model(ps);
  ModelState_Tabular state;
  Matrix3 zero(0.0);
  Matrix3 df_dsigma(0.0);
  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;
  double sqrtKG = std::sqrt(1.5 * state.bulkModulus / state.shearModulus);
  //std::cout << model;

  state.stressTensor = Matrix3(2000, 4000, 0, 4000, 2000, 0, 0, 0, 2000);
  state.updateStressInvariants();
  double z = state.zz;
  double rprime = state.rr*sqrtKG;
  double z_close = 0.0, rprime_close = 0.0;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -664.8953, 1.0e-4);
  EXPECT_NEAR(rprime_close, 781.4511, 1.0e-4);

  state.stressTensor = Matrix3(-2000, 4000, 0, 4000, -2000, 0, 0, 0, -2000);
  state.updateStressInvariants();
  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -3839.0782, 1.0e-4);
  EXPECT_NEAR(rprime_close, 1278.5554, 1.0e-4);

  state.stressTensor = Matrix3(3000, 0, 0, 0, 3000, 0, 0, 0, 3000);
  state.updateStressInvariants();
  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, 17.2339, 1.0e-4);
  EXPECT_NEAR(rprime_close, 0, 1.0e-10);

  state.stressTensor = Matrix3(3000, 1000, 0, 1000, 3000, 0, 0, 0, 3000);
  state.updateStressInvariants();
  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -1.67258, 1.0e-4);
  EXPECT_NEAR(rprime_close, 93.47397, 1.0e-4);

  state.stressTensor = Matrix3(-3000, 1000, 0, 1000, -3000, 0, 0, 0, -3000);
  state.updateStressInvariants();
  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -5214.01272, 1.0e-4);
  EXPECT_NEAR(rprime_close, 1355.79387, 1.0e-4);
}
