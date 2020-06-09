#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>

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

class YieldCondTabularCapTest : public ::testing::Test {

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
    xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_cap");
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
    xmlNewChild(rootNode, nullptr, BAD_CAST "cap_ellipticity_ratio", 
                BAD_CAST "0.7");

    // Print the document to stdout
    xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    if (!ps) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }
  }

  static void TearDownTestCase() {}

  static ProblemSpecP ps;
};

ProblemSpecP YieldCondTabularCapTest::ps = nullptr;

TEST_F(YieldCondTabularCapTest, constructorTest)
{
  // Create a model
  YieldCond_TabularCap model(ps);
  //std::cout << model;

  // Copy
  YieldCond_TabularCap modelCopy(&model);
  //std::cout << modelCopy;

  auto params = model.getParameters();
  ASSERT_DOUBLE_EQ(params.at("I1_min"), -19200);
  ASSERT_DOUBLE_EQ(params.at("I1_max"), 30);
  ASSERT_DOUBLE_EQ(params.at("sqrtJ2_max"), 900);
  ASSERT_DOUBLE_EQ(params.at("R"), 0.7);
}

TEST_F(YieldCondTabularCapTest, computeCapPoints)
{
  YieldCond_TabularCap model(ps);
  ModelState_TabularCap state;

  Polyline p_q_2000 = {
                       Uintah::Point(10,-100, 0),
                       Uintah::Point(-9.800000000000001,-1, 0),
                       Uintah::Point(-10,0,0),
                       Uintah::Point(-9.800000000000001,1, 0),
                       Uintah::Point(1.00000000000000e+01, 1.00000000000000e+02, 0),
                       Uintah::Point(1.52579301505928e+02, 2.46235181031721e+02, 0), 
                       Uintah::Point(2.51297661006088e+02, 3.47484780519064e+02, 0), 
                       Uintah::Point(4.00000000000000e+02, 5.00000000000000e+02, 0),
                       Uintah::Point(8.00000000000000e+02, 6.00000000000000e+02, 0),
                       Uintah::Point(1.39700000000000e+03, 6.74625000000000e+02, 0),
                       Uintah::Point(1.44955491287684e+03, 6.78602213895966e+02, 0),
                       Uintah::Point(1.50170985113316e+03, 6.77265814527447e+02, 0),
                       Uintah::Point(1.55306788419682e+03, 6.70481460560263e+02, 0),
                       Uintah::Point(1.60323814642538e+03, 6.57975013443945e+02, 0),
                       Uintah::Point(1.65183881182964e+03, 6.37351820852651e+02, 0),
                       Uintah::Point(1.69850000000000e+03, 6.11549251541155e+02, 0),
                       Uintah::Point(1.74286659111968e+03, 5.80720772263312e+02, 0),
                       Uintah::Point(1.78460092864098e+03, 5.45069392407036e+02, 0),
                       Uintah::Point(1.82338538905549e+03, 5.04847079544279e+02, 0),
                       Uintah::Point(1.85892479920074e+03, 4.60353430078503e+02, 0),
                       Uintah::Point(1.89094868270626e+03, 4.11933587232463e+02, 0),
                       Uintah::Point(1.91921331848202e+03, 3.59975416202563e+02, 0),
                       Uintah::Point(1.94350359558310e+03, 3.04905963997678e+02, 0),
                       Uintah::Point(1.96363465033390e+03, 2.47187248779554e+02, 0),
                       Uintah::Point(1.97945327325231e+03, 1.87311439937026e+02, 0),
                       Uintah::Point(1.99083907506636e+03, 1.25795505188495e+02, 0),
                       Uintah::Point(1.99770540294932e+03, 6.31754142851610e+01, 0),
                       Uintah::Point(2.00000000000000e+03, 0.00000000000000e+00, 0),
                       Uintah::Point(1.99770540294932e+03, -6.31754142851610e+01, 0),
                       Uintah::Point(1.99083907506636e+03, -1.25795505188495e+02, 0)
                      };
  Polyline p_q_2000_all;
  model.computeCapPoints(3.0*2000, p_q_2000_all);
  state.yield_f_pts = p_q_2000_all;

  int index = 0;
  for (const auto& p_q : p_q_2000_all) {
    ASSERT_NEAR(p_q.x(), p_q_2000[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_2000[index].y(), 1.0e-8);
    ++index;
  }
  ASSERT_NEAR(model.evalYieldConditionMax(&state), 6.78602213895966e+02, 1.0e-8);

  Polyline p_q_6400 = {
                       Uintah::Point(10,-100, 0),
                       Uintah::Point(-9.800000000000001,-1, 0),
                       Uintah::Point(-10,0,0),
                       Uintah::Point(-9.800000000000001,1, 0),
                       Uintah::Point(1.00000000000000e+01, 1.00000000000000e+02, 0),
                       Uintah::Point(1.52579301505928e+02, 2.46235181031721e+02, 0),
                       Uintah::Point(2.51297661006088e+02, 3.47484780519064e+02, 0),
                       Uintah::Point(4.00000000000000e+02, 5.00000000000000e+02, 0),
                       Uintah::Point(8.00000000000000e+02, 6.00000000000000e+02, 0),
                       Uintah::Point(1.60000000000000e+03, 7.00000000000000e+02, 0),
                       Uintah::Point(3.20000000000000e+03, 8.00000000000000e+02, 0),
                       Uintah::Point(4.47700000000000e+03, 8.39906250000000e+02, 0),
                       Uintah::Point(4.64460049330375e+03, 8.41927738232456e+02, 0),
                       Uintah::Point(4.81092544565351e+03, 8.37422823297398e+02, 0),
                       Uintah::Point(4.97470902373215e+03, 8.26310576036603e+02, 0),
                       Uintah::Point(5.13470473561526e+03, 8.08567464236765e+02, 0),
                       Uintah::Point(5.28969491732737e+03, 7.84230816382592e+02, 0),
                       Uintah::Point(5.43850000000000e+03, 7.53401506351659e+02, 0),
                       Uintah::Point(5.57998748710306e+03, 7.16245748413778e+02, 0),
                       Uintah::Point(5.71308057342722e+03, 6.72995911007188e+02, 0),
                       Uintah::Point(5.83676634022173e+03, 6.23950279935964e+02, 0),
                       Uintah::Point(5.95010346411779e+03, 5.69471726252073e+02, 0),
                       Uintah::Point(6.05222938116773e+03, 5.09985260459147e+02, 0),
                       Uintah::Point(6.14236685147748e+03, 4.45974482054336e+02, 0),
                       Uintah::Point(6.21982987447148e+03, 3.77976961026990e+02, 0),
                       Uintah::Point(6.28402890977130e+03, 3.06578614964925e+02, 0),
                       Uintah::Point(6.33447536395388e+03, 2.32407171100639e+02, 0),
                       Uintah::Point(6.37078530904248e+03, 1.56124826217554e+02, 0),
                       Uintah::Point(6.39268240443043e+03, 7.84202381454863e+01, 0),
                       Uintah::Point(6.40000000000000e+03, 0.00000000000000e+00, 0),
                       Uintah::Point(6.39268240443043e+03, -7.84202381454863e+01, 0),
                       Uintah::Point(6.37078530904248e+03, -1.56124826217554e+02, 0)
                       };

  Polyline p_q_6400_all;
  model.computeCapPoints(3.0*6400, p_q_6400_all);
  state.yield_f_pts = p_q_6400_all;

  index = 0;
  for (const auto& p_q : p_q_6400_all) {
    ASSERT_NEAR(p_q.x(), p_q_6400[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_6400[index].y(), 1.0e-8);
    ++index;
  }
  ASSERT_NEAR(model.evalYieldConditionMax(&state), 8.41927738232456e+02, 1.0e-8);

  Polyline p_q_10000 = {
                       Uintah::Point(10,-100, 0),
                       Uintah::Point(-9.800000000000001,-1, 0),
                       Uintah::Point(-10,0,0),
                       Uintah::Point(-9.800000000000001,1, 0),
                       Uintah::Point(1.00000000000000e+01, 1.00000000000000e+02, 0),
                       Uintah::Point(1.52579301505928e+02, 2.46235181031721e+02, 0),
                       Uintah::Point(2.51297661006088e+02, 3.47484780519064e+02, 0),
                       Uintah::Point(4.00000000000000e+02, 5.00000000000000e+02, 0),
                       Uintah::Point(8.00000000000000e+02, 6.00000000000000e+02, 0),
                       Uintah::Point(1.60000000000000e+03, 7.00000000000000e+02, 0),
                       Uintah::Point(3.20000000000000e+03, 8.00000000000000e+02, 0),
                       Uintah::Point(6.40000000000000e+03, 9.00000000000000e+02, 0),
                       Uintah::Point(6.72000000000000e+03, 9.10000000000000e+02, 0),
                       Uintah::Point(6.75200000000000e+03, 9.11000000000000e+02, 0),
                       Uintah::Point(6.78400000000000e+03, 9.12000000000000e+02, 0),
                       Uintah::Point(6.81600000000000e+03, 9.13000000000000e+02, 0),
                       Uintah::Point(6.84800000000000e+03, 9.14000000000000e+02, 0),
                       Uintah::Point(6.88000000000000e+03, 9.15000000000000e+02, 0),
                       Uintah::Point(6.91200000000000e+03, 9.16000000000000e+02, 0),
                       Uintah::Point(6.94400000000000e+03, 9.17000000000000e+02, 0),
                       Uintah::Point(6.97600000000000e+03, 9.18000000000000e+02, 0),
                       Uintah::Point(6.99700000000000e+03, 9.18656250000000e+02, 0),
                       Uintah::Point(7.25872869547122e+03, 9.23308383705311e+02, 0),
                       Uintah::Point(7.51846547753379e+03, 9.20748023765733e+02, 0),
                       Uintah::Point(7.77423359244287e+03, 9.10814734856869e+02, 0),
                       Uintah::Point(8.02408649040698e+03, 8.93415299037112e+02, 0),
                       Uintah::Point(8.26612264000732e+03, 8.68529554589357e+02, 0),
                       Uintah::Point(8.49850000000000e+03, 8.36215185588546e+02, 0),
                       Uintah::Point(8.71945003836219e+03, 7.96611284877295e+02, 0),
                       Uintah::Point(8.92729119188868e+03, 7.49940541734890e+02, 0),
                       Uintah::Point(9.12044166390320e+03, 6.96509938954404e+02, 0),
                       Uintah::Point(9.29743146268629e+03, 6.36709881346974e+02, 0),
                       Uintah::Point(9.45691358899984e+03, 5.71011717797553e+02, 0),
                       Uintah::Point(9.59767428756467e+03, 4.99963660743197e+02, 0),
                       Uintah::Point(9.71864228447106e+03, 4.24185149116703e+02, 0),
                       Uintah::Point(9.81889694022008e+03, 3.44359742165283e+02, 0),
                       Uintah::Point(9.89767525634607e+03, 2.61226670902461e+02, 0),
                       Uintah::Point(9.95437768229566e+03, 1.75571210127445e+02, 0),
                       Uintah::Point(9.98857267836951e+03, 8.82140658849937e+01, 0),
                       Uintah::Point(1.00000000000000e+04, 0.00000000000000e+00, 0),
                       Uintah::Point(9.98857267836951e+03, -8.82140658849937e+01, 0),
                       Uintah::Point(9.95437768229566e+03, -1.75571210127445e+02, 0)
                       };

  Polyline p_q_10000_all;
  model.computeCapPoints(3.0*10000, p_q_10000_all);
  state.yield_f_pts = p_q_10000_all;

  index = 0;
  for (const auto& p_q : p_q_10000_all) {
    //std::cout << std::setprecision(16) << "(" << p_q.x() << "," << p_q.y() << "),";
    ASSERT_NEAR(p_q.x(), p_q_10000[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_10000[index].y(), 1.0e-8);
    ++index;
  }
  //std::cout << std::endl;

  ASSERT_NEAR(model.evalYieldConditionMax(&state), 9.23308383705311e+02, 1.0e-8);
}

TEST_F(YieldCondTabularCapTest, evalYieldCondition)
{
  Matrix3 one(0.0); one.Identity();

  YieldCond_TabularCap model(ps);
  ModelState_TabularCap state;

  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;

  state.capX = -2000*3;
  Polyline p_q_2000_all;
  model.computeCapPoints(-state.capX, p_q_2000_all);
  state.yield_f_pts = p_q_2000_all;

  double p = 300; // Tension
  double sqrt_J2 = 1000;
  Matrix3 s(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  Matrix3 sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  p = 2;  // Tension
  sqrt_J2 = 39;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, -1);

  p = -1000;  // Compression
  sqrt_J2 = 625;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  p = -1000;  // Compression
  sqrt_J2 = 605;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, -1);

  p = -1000;  // Compression
  sqrt_J2 = 635;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  p = -7000;  // Compression
  sqrt_J2 = 1000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);
  //EXPECT_THROW(model.evalYieldCondition(&state), Uintah::InvalidValue);

  p = -1700;  // Compression
  sqrt_J2 = 6.10612759097964e+02;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  p = -1700.1;  // Compression
  sqrt_J2 = 6.10612759097964e+02;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  state.capX = -10000*3;
  Polyline p_q_10000_all;
  model.computeCapPoints(-state.capX, p_q_10000_all);
  state.yield_f_pts = p_q_10000_all;

  p = -9700;  // Compression
  sqrt_J2 = 4.37045077325265e+02;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

  p = -9700.1;  // Compression
  sqrt_J2 = 4.37045077325265e+02;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  EXPECT_EQ(model.evalYieldCondition(&state).first, 1);

}

TEST_F(YieldCondTabularCapTest, df_dsigma)
{
  YieldCond_TabularCap model(ps);
  ModelState_TabularCap state;
  Matrix3 zero(0.0);
  Matrix3 one(0.0); one.Identity();
  Matrix3 df_dsigma(0.0);
  //std::cout << model;

  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;
  state.stressTensor = zero;
  state.updateStressInvariants();

  state.capX = -2000*3;
  Polyline p_q_2000_all;
  model.computeCapPoints(-state.capX, p_q_2000_all);
  state.yield_f_pts = p_q_2000_all;

  // Zero everything (elastic)
  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  //ASSERT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,0), -0.57735, 1.0e-5);

  // Tension (p = 2000, J2 = 0)
  double p = 2000; 
  double sqrt_J2 = 0;
  Matrix3 s(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  Matrix3 sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);
  
  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  ASSERT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Tension (p = 300, J2 = 1000)
  p = 300; 
  sqrt_J2 = 1000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  ASSERT_NEAR(df_dsigma(0,0), 0.37068, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,1), 0.54212, 1.0e-5);

  // Tension (p = 2000, J2 = 4000)
  p = 2000; 
  sqrt_J2 = 4000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  //ASSERT_NEAR(df_dsigma(0,0), 0.266676, 1.0e-5);
  //ASSERT_NEAR(df_dsigma(0,1), 0.627157, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,0), 0.276741, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,1), 0.620581, 1.0e-5);

  // Compression (p = -2000, J2 = 4000)
  p = -2000; 
  sqrt_J2 = 4000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  //ASSERT_NEAR(df_dsigma(0,0), -0.0874635, 1.0e-5);
  //ASSERT_NEAR(df_dsigma(0,1), 0.698946, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,0), -0.065979, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,1),  0.702474, 1.0e-5);

  // Compression (p = -3000, J2 = 0)
  p = -3000; 
  sqrt_J2 = 0;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  ASSERT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Compression (p = -3000, J2 = 1000)
  p = -3000; 
  sqrt_J2 = 1000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  model.df_dsigma(zero, &state, df_dsigma);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  //ASSERT_NEAR(df_dsigma(0,0), -0.471837, 1.0e-5);
  //ASSERT_NEAR(df_dsigma(0,1), 0.407499, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,0), -0.477547, 1.0e-5);
  ASSERT_NEAR(df_dsigma(0,1),  0.397394, 1.0e-5);
}

TEST_F(YieldCondTabularCapTest, getClosestPoint)
{
  YieldCond_TabularCap model(ps);
  ModelState_TabularCap state;
  Matrix3 zero(0.0);
  Matrix3 df_dsigma(0.0);
  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;
  state.capX = -2000*3;
  double sqrtKG = std::sqrt(1.5 * state.bulkModulus / state.shearModulus);
  //std::cout << model;

  Polyline p_q_2000_all;
  model.computeCapPoints(-state.capX, p_q_2000_all);
  state.yield_f_pts = p_q_2000_all;

  state.stressTensor = Matrix3(2000, 4000, 0, 4000, 2000, 0, 0, 0, 2000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*2000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 4000, 1.0e-8);

  double z = state.zz;
  double rprime = state.rr*sqrtKG;
  double z_close = 0.0, rprime_close = 0.0;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -638.55547416470529, 1.0e-8);
  ASSERT_NEAR(rprime_close, 794.83369595956685, 1.0e-8);

  state.stressTensor = Matrix3(-2000, 4000, 0, 4000, -2000, 0, 0, 0, -2000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, -3*2000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 4000, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -2652.3249941689264, 1.0e-8);
  ASSERT_NEAR(rprime_close, 1166.2539083254812, 1.0e-8);

  state.stressTensor = Matrix3(3000, 0, 0, 0, 3000, 0, 0, 0, 3000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*3000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 0, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, 17.23390553531033, 1.0e-8);
  ASSERT_NEAR(rprime_close, 0, 1.0e-10);

  state.stressTensor = Matrix3(3000, 1000, 0, 1000, 3000, 0, 0, 0, 3000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*3000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 1000, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -4.8512584455542367, 1.0e-8);
  ASSERT_NEAR(rprime_close, 105.62075244805428, 1.0e-8);

  state.stressTensor = Matrix3(-3000, 1000, 0, 1000, -3000, 0, 0, 0, -3000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, -3*3000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 1000, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -3294.8725022226445, 1.0e-8);
  ASSERT_NEAR(rprime_close, 677.2749660674408, 1.0e-8);

  state.stressTensor = Matrix3(-3000, 0, 0, 0, -3000, 0, 0, 0, -3000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, -3*3000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 0, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -3463.108025469086, 1.0e-8);
  ASSERT_NEAR(rprime_close, 0, 1.0e-8);
}
