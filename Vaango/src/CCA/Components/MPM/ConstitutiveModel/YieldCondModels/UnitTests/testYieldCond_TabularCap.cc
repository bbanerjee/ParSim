#include <CCA/Components/MPM/ConstitutiveModel/YieldCondModels/YieldCond_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_TabularCap.h>
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
    xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "yield_condition");
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
    //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    if (!ps) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }

    // For internal var
    xmlDocPtr doc_intvar = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr rootNode_intvar = xmlNewNode(nullptr, BAD_CAST "internal_variable_model");
    xmlNewProp(rootNode_intvar, BAD_CAST "type", BAD_CAST "tabular_cap");
    xmlDocSetRootElement(doc_intvar, rootNode_intvar);
    xmlNewChild(rootNode_intvar, nullptr, BAD_CAST "filename", 
                BAD_CAST "DrySand_HydrostaticCapData.json");
    xmlNewChild(rootNode_intvar, nullptr, BAD_CAST "independent_variables", 
                BAD_CAST "PlasticStrainVol");
    xmlNewChild(rootNode_intvar, nullptr, BAD_CAST "dependent_variables", 
                BAD_CAST "Pressure");
    auto interp_intvar = xmlNewChild(rootNode_intvar, nullptr, BAD_CAST "interpolation",
                              BAD_CAST "");
    xmlNewProp(interp_intvar, BAD_CAST "type", BAD_CAST "linear");

    // Print the document to stdout
    //xmlSaveFormatFileEnc("-", doc_intvar, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ps_intvar = scinew ProblemSpec(xmlDocGetRootElement(doc_intvar), false);
    if (!ps_intvar) {
      std::cout << "**Error** Could not create ProblemSpec for internal var." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }
  }

  static void TearDownTestCase() {}

  static ProblemSpecP ps;
  static ProblemSpecP ps_intvar;
};

ProblemSpecP YieldCondTabularCapTest::ps = nullptr;
ProblemSpecP YieldCondTabularCapTest::ps_intvar = nullptr;

TEST_F(YieldCondTabularCapTest, constructorTest)
{
  try {
    IntVar_TabularCap intvar(ps_intvar);
  } catch (Uintah::ProblemSetupException e) {
    std::cout << e.message() << "\n";
  }
  // Create a model
  IntVar_TabularCap intvar(ps_intvar);
  YieldCond_TabularCap model(ps, &intvar);
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
  IntVar_TabularCap intvar(ps_intvar);
  YieldCond_TabularCap model(ps, &intvar);
  ModelState_TabularCap state;

  Polyline p_q_2000 = {
                      Uintah::Point(10, -100, 0),
                      Uintah::Point(-9.800000000000001, -1, 0),
                      Uintah::Point(-10, 0, 0),
                      Uintah::Point(-9.800000000000001, 1, 0),
                      Uintah::Point(10, 100, 0),
                      Uintah::Point(152.5793015059276, 246.2351810317206, 0),
                      Uintah::Point(251.2976610060878, 347.4847805190644, 0),
                      Uintah::Point(400, 500, 0),
                      Uintah::Point(800, 600, 0),
                      Uintah::Point(1397, 674.625, 0),
                      Uintah::Point(1402.262100905519, 675.2570499147716, 0),
                      Uintah::Point(1407.523801081682, 675.837526126856, 0),
                      Uintah::Point(1412.784699829647, 676.3662342125414, 0),
                      Uintah::Point(1418.044396511608, 676.8429837819892, 0),
                      Uintah::Point(1423.302490581298, 677.2675885394633, 0),
                      Uintah::Point(1428.558581614495, 677.6398663432042, 0),
                      Uintah::Point(1433.812269339519, 677.9596392649335, 0),
                      Uintah::Point(1439.063153667707, 678.2267336489664, 0),
                      Uintah::Point(1444.310834723891, 678.4409801709153, 0),
                      Uintah::Point(1449.554912876838, 678.6022138959662, 0),
                      Uintah::Point(1454.794988769695, 678.7102743367114, 0),
                      Uintah::Point(1460.030663350395, 678.7650055105157, 0),
                      Uintah::Point(1465.261537902048, 678.7662559964023, 0),
                      Uintah::Point(1470.487214073304, 678.7138789914394, 0),
                      Uintah::Point(1475.707293908691, 678.607732366608, 0),
                      Uintah::Point(1480.921379878919, 678.4476787221349, 0),
                      Uintah::Point(1486.129074911155, 678.2335854422712, 0),
                      Uintah::Point(1491.329982419259, 677.9653247495007, 0),
                      Uintah::Point(1496.523706333989, 677.6427737581589, 0),
                      Uintah::Point(1501.709851133159, 677.2658145274471, 0),
                      Uintah::Point(1506.888021871765, 676.8343341138213, 0),
                      Uintah::Point(1512.057824212057, 676.3482246227426, 0),
                      Uintah::Point(1517.21886445357, 675.8073832597667, 0),
                      Uintah::Point(1522.370749563109, 675.2117123809629, 0),
                      Uintah::Point(1527.513087204676, 674.5611195426377, 0),
                      Uintah::Point(1532.645485769351, 673.8555175503527, 0),
                      Uintah::Point(1537.767554405111, 673.0948245072168, 0),
                      Uintah::Point(1542.8789030466, 672.2789638614373, 0),
                      Uintah::Point(1547.979142444828, 671.4078644531147, 0),
                      Uintah::Point(1553.06788419682, 670.4814605602627, 0),
                      Uintah::Point(1558.144740775189, 669.49969194404, 0),
                      Uintah::Point(1563.209325557651, 668.4625038931773, 0),
                      Uintah::Point(1568.261252856465, 667.3698472675836, 0),
                      Uintah::Point(1573.30013794781, 666.2216785411166, 0),
                      Uintah::Point(1578.325597101077, 665.0179598435027, 0),
                      Uintah::Point(1583.337247608093, 663.7586590013905, 0),
                      Uintah::Point(1588.334707812271, 662.4437495785251, 0),
                      Uintah::Point(1593.317597137665, 661.0732109150248, 0),
                      Uintah::Point(1598.285536117964, 659.6470281657517, 0),
                      Uintah::Point(1603.238146425378, 657.9750134439453, 0),
                      Uintah::Point(1608.175050899459, 656.1491164003291, 0),
                      Uintah::Point(1613.095873575816, 654.2704267504326, 0),
                      Uintah::Point(1618.000239714751, 652.3390312910731, 0),
                      Uintah::Point(1622.887775829795, 650.3550217102718, 0),
                      Uintah::Point(1627.758109716149, 648.3184945971477, 0),
                      Uintah::Point(1632.610870479032, 646.2295514511618, 0),
                      Uintah::Point(1637.445688561924, 644.0882986907118, 0),
                      Uintah::Point(1642.262195774708, 641.8948476610698, 0),
                      Uintah::Point(1647.060025321712, 639.649314641657, 0),
                      Uintah::Point(1651.838811829642, 637.3518208526507, 0),
                      Uintah::Point(1656.598191375402, 635.0024924609183, 0),
                      Uintah::Point(1661.337801513814, 632.6014605852729, 0),
                      Uintah::Point(1666.057281305215, 630.1488613010456, 0),
                      Uintah::Point(1670.756271342947, 627.6448356439705, 0),
                      Uintah::Point(1675.434413780725, 625.0895296133772, 0),
                      Uintah::Point(1680.091352359892, 622.4830941746852, 0),
                      Uintah::Point(1684.726732436544, 619.8256852611987, 0),
                      Uintah::Point(1689.340201008541, 617.1174637751935, 0),
                      Uintah::Point(1693.931406742391, 614.3585955882978, 0),
                      Uintah::Point(1698.5, 611.5492515411549, 0),
                      Uintah::Point(1703.045632865305, 608.6896074423742, 0),
                      Uintah::Point(1707.567959170763, 605.779844066756, 0),
                      Uintah::Point(1712.066634523717, 602.8201471527966, 0),
                      Uintah::Point(1716.541316332623, 599.8107073994635, 0),
                      Uintah::Point(1720.991663833135, 596.7517204622409, 0),
                      Uintah::Point(1725.417338114061, 593.643386948442, 0),
                      Uintah::Point(1729.818002143171, 590.485912411785, 0),
                      Uintah::Point(1734.19332079286, 587.2795073462299, 0),
                      Uintah::Point(1738.542960865674, 584.0243871790758, 0),
                      Uintah::Point(1742.866591119681, 580.7207722633121, 0),
                      Uintah::Point(1747.163882293697, 577.3688878692273, 0),
                      Uintah::Point(1751.434507132361, 573.968964175268, 0),
                      Uintah::Point(1755.678140411059, 570.5212362581497, 0),
                      Uintah::Point(1759.894458960685, 567.0259440822176, 0),
                      Uintah::Point(1764.083141692259, 563.4833324880542, 0),
                      Uintah::Point(1768.243869621372, 559.8936511803334, 0),
                      Uintah::Point(1772.376325892485, 556.2571547149215, 0),
                      Uintah::Point(1776.480195803052, 552.5741024852219, 0),
                      Uintah::Point(1780.555166827492, 548.8447587077642, 0),
                      Uintah::Point(1784.600928640983, 545.0693924070364, 0),
                      Uintah::Point(1788.617173143101, 541.2482773995608, 0),
                      Uintah::Point(1792.603594481276, 537.3816922772129, 0),
                      Uintah::Point(1796.55988907409, 533.4699203897821, 0),
                      Uintah::Point(1800.485755634392, 529.5132498267774, 0),
                      Uintah::Point(1804.380895192243, 525.5119733984748, 0),
                      Uintah::Point(1808.245011117687, 521.4663886162097, 0),
                      Uintah::Point(1812.077809143334, 517.3767976719145, 0),
                      Uintah::Point(1815.878997386775, 513.2435074169002, 0),
                      Uintah::Point(1819.64828637281, 509.0668293398867, 0),
                      Uintah::Point(1823.385389055488, 504.8470795442789, 0),
                      Uintah::Point(1827.090020839971, 500.5845787246945, 0),
                      Uintah::Point(1830.761899604207, 496.2796521427406, 0),
                      Uintah::Point(1834.400745720409, 491.9326296020448, 0),
                      Uintah::Point(1838.00628207636, 487.5438454225396, 0),
                      Uintah::Point(1841.578234096505, 483.1136384140048, 0),
                      Uintah::Point(1845.116329762869, 478.6423518488687, 0),
                      Uintah::Point(1848.620299635768, 474.1303334342702, 0),
                      Uintah::Point(1852.089876874332, 469.577935283387, 0),
                      Uintah::Point(1855.524797256819, 464.985513886029, 0),
                      Uintah::Point(1858.924799200744, 460.3534300785033, 0),
                      Uintah::Point(1862.289623782795, 455.6820490127525, 0),
                      Uintah::Point(1865.619014758553, 450.971740124769, 0),
                      Uintah::Point(1868.912718582006, 446.2228771022906, 0),
                      Uintah::Point(1872.170484424853, 441.4358378517799, 0),
                      Uintah::Point(1875.392064195615, 436.611004464691, 0),
                      Uintah::Point(1878.577212558518, 431.74876318303, 0),
                      Uintah::Point(1881.725686952182, 426.8495043642097, 0),
                      Uintah::Point(1884.837247608093, 421.9136224452072, 0),
                      Uintah::Point(1887.91165756886, 416.9415159060249, 0),
                      Uintah::Point(1890.948682706262, 411.9335872324628, 0),
                      Uintah::Point(1893.948091739076, 406.8902428782058, 0),
                      Uintah::Point(1896.90965625069, 401.8118932262307, 0),
                      Uintah::Point(1899.833150706502, 396.6989525495384, 0),
                      Uintah::Point(1902.718352471091, 391.5518389712186, 0),
                      Uintah::Point(1905.56504182517, 386.3709744238483, 0),
                      Uintah::Point(1908.373001982325, 381.1567846082346, 0),
                      Uintah::Point(1911.142019105518, 375.9096989515052, 0),
                      Uintah::Point(1913.871882323374, 370.630150564552, 0),
                      Uintah::Point(1916.56238374624, 365.3185761988369, 0),
                      Uintah::Point(1919.213318482017, 359.9754162025629, 0),
                      Uintah::Point(1921.82448465176, 354.6011144762197, 0),
                      Uintah::Point(1924.395683405055, 349.1961184275082, 0),
                      Uintah::Point(1926.926718935165, 343.7608789256537, 0),
                      Uintah::Point(1929.417398493933, 338.2958502551116, 0),
                      Uintah::Point(1931.867532406468, 332.8014900686749, 0),
                      Uintah::Point(1934.276934085586, 327.2782593399901, 0),
                      Uintah::Point(1936.645420046021, 321.7266223154891, 0),
                      Uintah::Point(1938.972809918398, 316.1470464657444, 0),
                      Uintah::Point(1941.258926462966, 310.5400024362561, 0),
                      Uintah::Point(1943.5035955831, 304.9059639976777, 0),
                      Uintah::Point(1945.706646338555, 299.2454079954891, 0),
                      Uintah::Point(1947.867910958488, 293.5588142991256, 0),
                      Uintah::Point(1949.98722485423, 287.8466657505695, 0),
                      Uintah::Point(1952.064426631822, 282.1094481124154, 0),
                      Uintah::Point(1954.099358104306, 276.3476500154148, 0),
                      Uintah::Point(1956.091864303773, 270.5617629055108, 0),
                      Uintah::Point(1958.041793493161, 264.7522809903713, 0),
                      Uintah::Point(1959.948997177813, 258.9197011854289, 0),
                      Uintah::Point(1961.813330116784, 253.0645230594376, 0),
                      Uintah::Point(1963.634650333903, 247.1872487795544, 0),
                      Uintah::Point(1965.412819128584, 241.2883830559565, 0),
                      Uintah::Point(1967.147701086388, 235.3684330860024, 0),
                      Uintah::Point(1968.839164089338, 229.4279084979463, 0),
                      Uintah::Point(1970.487079325977, 223.4673212942173, 0),
                      Uintah::Point(1972.091321301181, 217.4871857942707, 0),
                      Uintah::Point(1973.65176784571, 211.4880185770229, 0),
                      Uintah::Point(1975.168300125521, 205.4703384228793, 0),
                      Uintah::Point(1976.640802650806, 199.4346662553657, 0),
                      Uintah::Point(1978.0691632848, 193.3815250823732, 0),
                      Uintah::Point(1979.453273252308, 187.3114399370263, 0),
                      Uintah::Point(1980.793027147999, 181.2249378181852, 0),
                      Uintah::Point(1982.088322944426, 175.1225476305933, 0),
                      Uintah::Point(1983.339061999799, 169.0048001246791, 0),
                      Uintah::Point(1984.545149065497, 162.8722278360244, 0),
                      Uintah::Point(1985.70649229332, 156.7253650245086, 0),
                      Uintah::Point(1986.823003242485, 150.5647476131408, 0),
                      Uintah::Point(1987.89459688636, 144.3909131265896, 0),
                      Uintah::Point(1988.921191618942, 138.2044006294225, 0),
                      Uintah::Point(1989.902709261065, 132.0057506640659, 0),
                      Uintah::Point(1990.839075066362, 125.7955051884951, 0),
                      Uintah::Point(1991.730217726951, 119.5742075136691, 0),
                      Uintah::Point(1992.576069378868, 113.3424022407166, 0),
                      Uintah::Point(1993.376565607236, 107.1006351978896, 0),
                      Uintah::Point(1994.131645451167, 100.8494533772916, 0),
                      Uintah::Point(1994.841251408408, 94.58940487139482, 0),
                      Uintah::Point(1995.505329439717, 88.32103880935666, 0),
                      Uintah::Point(1996.123828972982, 82.0449052931466, 0),
                      Uintah::Point(1996.696702907069, 75.76155533349623, 0),
                      Uintah::Point(1997.223907615409, 69.47154078568322, 0),
                      Uintah::Point(1997.705402949322, 63.17541428516103, 0),
                      Uintah::Point(1998.141152241076, 56.87372918304628, 0),
                      Uintah::Point(1998.531122306674, 50.56703948147535, 0),
                      Uintah::Point(1998.875283448386, 44.25589976884206, 0),
                      Uintah::Point(1999.173609457008, 37.94086515492818, 0),
                      Uintah::Point(1999.42607761386, 31.62249120593868, 0),
                      Uintah::Point(1999.632668692515, 25.30133387945356, 0),
                      Uintah::Point(1999.793366960261, 18.97794945930804, 0),
                      Uintah::Point(1999.908160179304, 12.65289449041318, 0),
                      Uintah::Point(1999.977039607695, 6.326725713528697, 0),
                      Uintah::Point(2000, 0, 0),
                      Uintah::Point(1999.977039607695, -6.326725713528697, 0),
                      Uintah::Point(1999.908160179304, -12.65289449041318, 0)
                      };
  Polyline p_q_2000_all;
  model.computeCapPoints(3.0*2000, p_q_2000_all);
  //state.yield_f_pts = p_q_2000_all;
  state.updateYieldSurface(p_q_2000_all);

  int index = 0;
  for (const auto& p_q : p_q_2000_all) {
    ASSERT_NEAR(p_q.x(), p_q_2000[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_2000[index].y(), 1.0e-8);
    ++index;
  }
  ASSERT_NEAR(model.evalYieldConditionMax(&state), 6.7876625599640226e+02, 1.0e-8);

  /*
  for (const auto& p_q : p_q_2000_all) {
    std::cout << std::setprecision(16) 
               << "Uintah::Point(" << p_q.x() << ", " << p_q.y() << ", 0)," << "\n";
  }
  std::cout << "p_q_2000_x = ";
  for (const auto& p_q : p_q_2000_all) {
    std::cout << p_q.x() << ", ";
  }
  std::cout << "\n";
  std::cout << "p_q_2000_y = ";
  for (const auto& p_q : p_q_2000_all) {
    std::cout << p_q.y() << ", ";
  }
  std::cout << "\n";
  */


  Polyline p_q_6400 = {
                       Uintah::Point(10, -100, 0),
                       Uintah::Point(-9.800000000000001, -1, 0),
                       Uintah::Point(-10, 0, 0),
                       Uintah::Point(-9.800000000000001, 1, 0),
                       Uintah::Point(10, 100, 0),
                       Uintah::Point(152.5793015059276, 246.2351810317206, 0),
                       Uintah::Point(251.2976610060878, 347.4847805190644, 0),
                       Uintah::Point(400, 500, 0),
                       Uintah::Point(800, 600, 0),
                       Uintah::Point(1600, 700, 0),
                       Uintah::Point(3200, 800, 0),
                       Uintah::Point(4477, 839.90625, 0),
                       Uintah::Point(4493.781127763373, 840.3986592182869, 0),
                       Uintah::Point(4510.560977578896, 840.8269490249307, 0),
                       Uintah::Point(4527.33827159604, 841.1909670404889, 0),
                       Uintah::Point(4544.11173215891, 841.490565853038, 0),
                       Uintah::Point(4560.880081903541, 841.7256030658582, 0),
                       Uintah::Point(4577.642043855181, 841.8959413447041, 0),
                       Uintah::Point(4594.39634152553, 842.0014484646426, 0),
                       Uintah::Point(4611.141699009952, 842.04199735645, 0),
                       Uintah::Point(4627.876841084646, 842.0174661525457, 0),
                       Uintah::Point(4644.600493303747, 841.9277382324562, 0),
                       Uintah::Point(4661.311382096391, 841.7727022677878, 0),
                       Uintah::Point(4678.008234863698, 841.5522522666985, 0),
                       Uintah::Point(4694.689780075685, 841.2662876178521, 0),
                       Uintah::Point(4711.354747368098, 840.9147131338405, 0),
                       Uintah::Point(4728.001867639159, 840.4974390940617, 0),
                       Uintah::Point(4744.629873146206, 840.0143812870385, 0),
                       Uintah::Point(4761.237497602241, 839.4654610521609, 0),
                       Uintah::Point(4777.823476272364, 838.8506053208447, 0),
                       Uintah::Point(4794.386546070083, 838.1697466570852, 0),
                       Uintah::Point(4810.925445653507, 837.4228232973983, 0),
                       Uintah::Point(4827.438915521399, 836.6097791901318, 0),
                       Uintah::Point(4843.925698109096, 835.7305640341332, 0),
                       Uintah::Point(4860.38453788427, 834.7851333167647, 0),
                       Uintah::Point(4876.814181442552, 833.7734483512462, 0),
                       Uintah::Point(4893.213377602972, 832.6954763133167, 0),
                       Uintah::Point(4909.580877503252, 831.5511902772021, 0),
                       Uintah::Point(4925.915434694906, 830.3405692508727, 0),
                       Uintah::Point(4942.215805238161, 829.0635982105809, 0),
                       Uintah::Point(4958.480747796691, 827.7202681346641, 0),
                       Uintah::Point(4974.709023732147, 826.3105760366027, 0),
                       Uintah::Point(4990.899397198488, 824.834524997318, 0),
                       Uintah::Point(5007.05063523609, 823.2921241967008, 0),
                       Uintah::Point(5023.161507865643, 821.6833889443554, 0),
                       Uintah::Point(5039.230788181823, 820.0083407095508, 0),
                       Uintah::Point(5055.257252446717, 818.2670071503634, 0),
                       Uintah::Point(5071.239680183024, 816.4594221420017, 0),
                       Uintah::Point(5087.176854266992, 814.5856258043035, 0),
                       Uintah::Point(5103.067561021113, 812.6456645283893, 0),
                       Uintah::Point(5118.910590306541, 810.6395910024656, 0),
                       Uintah::Point(5134.704735615262, 808.5674642367645, 0),
                       Uintah::Point(5150.448794161955, 806.429349587611, 0),
                       Uintah::Point(5166.141566975612, 804.2253187806042, 0),
                       Uintah::Point(5181.781858990824, 801.9554499329053, 0),
                       Uintah::Point(5197.368479138799, 799.619827574621, 0),
                       Uintah::Point(5212.900240438068, 797.2185426692726, 0),
                       Uintah::Point(5228.375960084873, 794.7516926333391, 0),
                       Uintah::Point(5243.794459543248, 792.2193813548691, 0),
                       Uintah::Point(5259.154564634764, 789.6217192111469, 0),
                       Uintah::Point(5274.455105627948, 786.9588230854076, 0),
                       Uintah::Point(5289.694917327365, 784.2308163825918, 0),
                       Uintah::Point(5304.872839162352, 781.4378290441275, 0),
                       Uintah::Point(5319.987715275396, 778.5799975617368, 0),
                       Uintah::Point(5335.038394610162, 775.6574649902532, 0),
                       Uintah::Point(5350.023730999148, 772.670380959446, 0),
                       Uintah::Point(5364.94258325097, 769.6189016848416, 0),
                       Uintah::Point(5379.793815237268, 766.5031899775328, 0),
                       Uintah::Point(5394.576295979227, 763.323415252974, 0),
                       Uintah::Point(5409.288899733706, 760.0797535387483, 0),
                       Uintah::Point(5423.930506078967, 756.7723874813051, 0),
                       Uintah::Point(5438.5, 753.4015063516592, 0),
                       Uintah::Point(5452.996271973434, 749.9673060500442, 0),
                       Uintah::Point(5467.418218052034, 746.4699891095165, 0),
                       Uintah::Point(5481.76473994877, 742.9097646985017, 0),
                       Uintah::Point(5496.034745120453, 739.2868486222784, 0),
                       Uintah::Point(5510.227146850943, 735.6014633233959, 0),
                       Uintah::Point(5524.340864333897, 731.8538378810163, 0),
                       Uintah::Point(5538.374822755088, 728.0442080091814, 0),
                       Uintah::Point(5552.327953374246, 724.1728160539944, 0),
                       Uintah::Point(5566.199193606453, 720.2399109897179, 0),
                       Uintah::Point(5579.987487103062, 716.2457484137775, 0),
                       Uintah::Point(5593.691783832137, 712.190590540673, 0),
                       Uintah::Point(5607.311040158425, 708.0747061947892, 0),
                       Uintah::Point(5620.84421892283, 703.8983708021065, 0),
                       Uintah::Point(5634.290289521389, 699.6618663808057, 0),
                       Uintah::Point(5647.64822798377, 695.3654815307646, 0),
                       Uintah::Point(5660.917017051241, 691.0095114219449, 0),
                       Uintah::Point(5674.095646254143, 686.5942577816667, 0),
                       Uintah::Point(5687.183111988837, 682.1200288807676, 0),
                       Uintah::Point(5700.17841759414, 677.5871395186462, 0),
                       Uintah::Point(5713.080573427215, 672.9959110071878, 0),
                       Uintah::Point(5725.888596938943, 668.3466711535705, 0),
                       Uintah::Point(5738.601512748745, 663.639754241952, 0),
                       Uintah::Point(5751.218352718863, 658.8755010140355, 0),
                       Uintah::Point(5763.738156028085, 654.0542586485162, 0),
                       Uintah::Point(5776.159969244914, 649.1763807394037, 0),
                       Uintah::Point(5788.482846400185, 644.2422272732271, 0),
                       Uintah::Point(5800.705849059089, 639.2521646051181, 0),
                       Uintah::Point(5812.828046392651, 634.2065654337745, 0),
                       Uintah::Point(5824.848515248614, 629.1058087753065, 0),
                       Uintah::Point(5836.766340221731, 623.9502799359636, 0),
                       Uintah::Point(5848.580613723491, 618.7403704837488, 0),
                       Uintah::Point(5860.290436051226, 613.4764782189166, 0),
                       Uintah::Point(5871.894915456629, 608.1590071433608, 0),
                       Uintah::Point(5883.393168213664, 602.7883674288935, 0),
                       Uintah::Point(5894.784318685868, 597.364975384416, 0),
                       Uintah::Point(5906.067499393029, 591.8892534219884, 0),
                       Uintah::Point(5917.241851077251, 586.3616300217972, 0),
                       Uintah::Point(5928.306522768391, 580.7825396960276, 0),
                       Uintah::Point(5939.26067184886, 575.1524229516428, 0),
                       Uintah::Point(5950.103464117795, 569.4717262520735, 0),
                       Uintah::Point(5960.834073854586, 563.7409019778248, 0),
                       Uintah::Point(5971.451683881755, 557.9604083860011, 0),
                       Uintah::Point(5981.955485627192, 552.1307095687595, 0),
                       Uintah::Point(5992.344679185726, 546.2522754106903, 0),
                       Uintah::Point(6002.618473380046, 540.3255815451352, 0),
                       Uintah::Point(6012.776085820944, 534.3511093094467, 0),
                       Uintah::Point(6022.816742966908, 528.3293456991955, 0),
                       Uintah::Point(6032.739680183024, 522.2607833213309, 0),
                       Uintah::Point(6042.544141799202, 516.1459203463022, 0),
                       Uintah::Point(6052.229381167731, 509.9852604591474, 0),
                       Uintah::Point(6061.794660720136, 503.779312809556, 0),
                       Uintah::Point(6071.239252023345, 497.5285919609134, 0),
                       Uintah::Point(6080.562435835164, 491.2336178383346, 0),
                       Uintah::Point(6089.76350215905, 484.8949156756951, 0),
                       Uintah::Point(6098.841750298179, 478.5130159616664, 0),
                       Uintah::Point(6107.796488908807, 472.0884543847656, 0),
                       Uintah::Point(6116.62703605292, 465.6217717774272, 0),
                       Uintah::Point(6125.332719250162, 459.1135140591057, 0),
                       Uintah::Point(6133.912875529054, 452.5642321784179, 0),
                       Uintah::Point(6142.366851477475, 445.9744820543355, 0),
                       Uintah::Point(6150.694003292427, 439.3448245164365, 0),
                       Uintah::Point(6158.893696829058, 432.6758252442252, 0),
                       Uintah::Point(6166.965307648959, 425.9680547055319, 0),
                       Uintah::Point(6174.908221067717, 419.2220880940018, 0),
                       Uintah::Point(6182.721832201721, 412.4385052656835, 0),
                       Uintah::Point(6190.405546014232, 405.6178906747284, 0),
                       Uintah::Point(6197.958777360695, 398.7608333082126, 0),
                       Uintah::Point(6205.380951033299, 391.8679266200905, 0),
                       Uintah::Point(6212.671501804782, 384.9397684642946, 0),
                       Uintah::Point(6219.829874471478, 377.9769610269899, 0),
                       Uintah::Point(6226.855523895592, 370.9801107579979, 0),
                       Uintah::Point(6233.747915046722, 363.9498283013999, 0),
                       Uintah::Point(6240.506523042593, 356.8867284253336, 0),
                       Uintah::Point(6247.130833189043, 349.7914299509945, 0),
                       Uintah::Point(6253.620341019205, 342.6645556808558, 0),
                       Uintah::Point(6259.974552331933, 335.5067323261185, 0),
                       Uintah::Point(6266.192983229434, 328.3185904334072, 0),
                       Uintah::Point(6272.275160154119, 321.1007643107212, 0),
                       Uintah::Point(6278.220619924668, 313.8538919526587, 0),
                       Uintah::Point(6284.028909771302, 306.5786149649252, 0),
                       Uintah::Point(6289.699587370259, 299.27557848814, 0),
                       Uintah::Point(6295.232220877486, 291.9454311209569, 0),
                       Uintah::Point(6300.626388961521, 284.5888248425111, 0),
                       Uintah::Point(6305.88168083558, 277.2064149342086, 0),
                       Uintah::Point(6310.99769628884, 269.7988599008712, 0),
                       Uintah::Point(6315.974045716917, 262.3668213912537, 0),
                       Uintah::Point(6320.810350151535, 254.9109641179462, 0),
                       Uintah::Point(6325.506241289388, 247.4319557766786, 0),
                       Uintah::Point(6330.061361520182, 239.9304669650414, 0),
                       Uintah::Point(6334.475363953879, 232.407171100639, 0),
                       Uintah::Point(6338.747912447101, 224.8627443386899, 0),
                       Uintah::Point(6342.878681628741, 217.2978654890907, 0),
                       Uintah::Point(6346.867356924732, 209.7132159329607, 0),
                       Uintah::Point(6350.713634582007, 202.1094795386799, 0),
                       Uintah::Point(6354.417221691632, 194.4873425774398, 0),
                       Uintah::Point(6357.977836211108, 186.8474936383219, 0),
                       Uintah::Point(6361.395206985855, 179.1906235429202, 0),
                       Uintah::Point(6364.669073769857, 171.5174252595248, 0),
                       Uintah::Point(6367.799187245485, 163.8285938168829, 0),
                       Uintah::Point(6370.785309042476, 156.1248262175536, 0),
                       Uintah::Point(6373.627211756097, 148.4068213508747, 0),
                       Uintah::Point(6376.32467896445, 140.6752799055562, 0),
                       Uintah::Point(6378.877505244966, 132.9309042819205, 0),
                       Uintah::Point(6381.285496190039, 125.1743985038038, 0),
                       Uintah::Point(6383.548468421837, 117.4064681301375, 0),
                       Uintah::Point(6385.666249606262, 109.6278201662259, 0),
                       Uintah::Point(6387.638678466077, 101.8391629747392, 0),
                       Uintah::Point(6389.465604793189, 94.04120618643744, 0),
                       Uintah::Point(6391.146889460085, 86.23466061064399, 0),
                       Uintah::Point(6392.682404430427, 78.42023814548631, 0),
                       Uintah::Point(6394.072032768805, 70.59865168792096, 0),
                       Uintah::Point(6395.315668649642, 62.77061504356121, 0),
                       Uintah::Point(6396.413217365251, 54.93684283632492, 0),
                       Uintah::Point(6397.364595333045, 47.09805041792032, 0),
                       Uintah::Point(6398.169730101912, 39.25495377718778, 0),
                       Uintah::Point(6398.828560357721, 31.40826944931547, 0),
                       Uintah::Point(6399.341035927997, 23.55871442494679, 0),
                       Uintah::Point(6399.707117785741, 15.70700605919762, 0),
                       Uintah::Point(6399.926778052401, 7.853861980601386, 0),
                       Uintah::Point(6400, 0, 0),
                       Uintah::Point(6399.926778052401, -7.853861980601386, 0),
                       Uintah::Point(6399.707117785741, -15.70700605919762, 0)
                       };

  Polyline p_q_6400_all;
  model.computeCapPoints(3.0*6400, p_q_6400_all);
  //state.yield_f_pts = p_q_6400_all;
  state.updateYieldSurface(p_q_6400_all);

  index = 0;
  for (const auto& p_q : p_q_6400_all) {
    ASSERT_NEAR(p_q.x(), p_q_6400[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_6400[index].y(), 1.0e-8);
    ++index;
  }
  ASSERT_NEAR(model.evalYieldConditionMax(&state), 842.04199735645, 1.0e-8);

  /*
  for (const auto& p_q : p_q_6400_all) {
    std::cout << std::setprecision(16) 
              << "Uintah::Point(" << p_q.x() << ", " << p_q.y() << ", 0)," << "\n";
  }
  std::cout << "p_q_6400_x = ";
  for (const auto& p_q : p_q_6400_all) {
    std::cout << p_q.x() << ", ";
  }
  std::cout << "\n";
  std::cout << "p_q_6400_y = ";
  for (const auto& p_q : p_q_6400_all) {
    std::cout << p_q.y() << ", ";
  }
  std::cout << "\n";
  */


  Polyline p_q_10000 = {
                       Uintah::Point(10, -100, 0),
                       Uintah::Point(-9.800000000000001, -1, 0),
                       Uintah::Point(-10, 0, 0),
                       Uintah::Point(-9.800000000000001, 1, 0),
                       Uintah::Point(10, 100, 0),
                       Uintah::Point(152.5793015059276, 246.2351810317206, 0),
                       Uintah::Point(251.2976610060878, 347.4847805190644, 0),
                       Uintah::Point(400, 500, 0),
                       Uintah::Point(800, 600, 0),
                       Uintah::Point(1600, 700, 0),
                       Uintah::Point(3200, 800, 0),
                       Uintah::Point(6400, 900, 0),
                       Uintah::Point(6720.000000000001, 910, 0),
                       Uintah::Point(6751.999999999999, 910.9999999999998, 0),
                       Uintah::Point(6783.999999999997, 911.9999999999995, 0),
                       Uintah::Point(6815.999999999995, 912.9999999999991, 0),
                       Uintah::Point(6847.999999999994, 913.9999999999991, 0),
                       Uintah::Point(6879.999999999992, 914.9999999999991, 0),
                       Uintah::Point(6911.99999999999, 915.9999999999982, 0),
                       Uintah::Point(6943.999999999988, 916.9999999999982, 0),
                       Uintah::Point(6975.999999999986, 917.9999999999982, 0),
                       Uintah::Point(6997, 918.6562499999982, 0),
                       Uintah::Point(7023.205786101616, 919.4401700182183, 0),
                       Uintah::Point(7049.409576531162, 920.153884025349, 0),
                       Uintah::Point(7075.609375768543, 920.797150643911, 0),
                       Uintah::Point(7101.803188597611, 921.3697339752219, 0),
                       Uintah::Point(7127.989020258104, 921.8714036742936, 0),
                       Uintah::Point(7154.164876597561, 922.3019350242655, 0),
                       Uintah::Point(7180.328764223175, 922.6611090103238, 0),
                       Uintah::Point(7206.478690653608, 922.9487123931093, 0),
                       Uintah::Point(7232.612664470718, 923.1645377815805, 0),
                       Uintah::Point(7258.728695471217, 923.308383705308, 0),
                       Uintah::Point(7284.824794818232, 923.3800546861789, 0),
                       Uintah::Point(7310.898975192763, 923.3793613094956, 0),
                       Uintah::Point(7336.949250945025, 923.3061202944319, 0),
                       Uintah::Point(7362.973638245658, 923.1601545638335, 0),
                       Uintah::Point(7388.970155236815, 922.9412933133532, 0),
                       Uintah::Point(7414.936822183076, 922.6493720798455, 0),
                       Uintah::Point(7440.871661622221, 922.2842328091036, 0),
                       Uintah::Point(7466.772698515813, 921.8457239227848, 0),
                       Uintah::Point(7492.637960399615, 921.3337003845995, 0),
                       Uintah::Point(7518.465477533792, 920.7480237657256, 0),
                       Uintah::Point(7544.253283052919, 920.0885623093635, 0),
                       Uintah::Point(7569.999413115765, 919.3551909945239, 0),
                       Uintah::Point(7595.701907054843, 918.5477915989053, 0),
                       Uintah::Point(7621.358807525732, 917.6662527609333, 0),
                       Uintah::Point(7646.968160656123, 916.7104700408834, 0),
                       Uintah::Point(7672.528016194627, 915.6803459811002, 0),
                       Uintah::Point(7698.036427659285, 914.5757901652889, 0),
                       Uintah::Point(7723.491452485802, 913.3967192768201, 0),
                       Uintah::Point(7748.891152175488, 912.1430571560937, 0),
                       Uintah::Point(7774.23359244287, 910.8147348568618, 0),
                       Uintah::Point(7799.516843363006, 909.4116907015989, 0),
                       Uintah::Point(7824.738979518448, 907.9338703357706, 0),
                       Uintah::Point(7849.89808014588, 906.381226781099, 0),
                       Uintah::Point(7874.992229282379, 904.7537204877021, 0),
                       Uintah::Point(7900.019515911332, 903.0513193852072, 0),
                       Uintah::Point(7924.978034107968, 901.273998932673, 0),
                       Uintah::Point(7949.865883184491, 899.4217421674788, 0),
                       Uintah::Point(7974.681167834842, 897.4945397529509, 0),
                       Uintah::Point(7999.421998279014, 895.4923900249315, 0),
                       Uintah::Point(8024.086490406984, 893.4152990371056, 0),
                       Uintah::Point(8048.672765922181, 891.2632806051261, 0),
                       Uintah::Point(8073.178952484537, 889.0363563495536, 0),
                       Uintah::Point(8097.603183853065, 886.7345557375349, 0),
                       Uintah::Point(8121.943600027984, 884.3579161232464, 0),
                       Uintah::Point(8146.198347392365, 881.9064827870495, 0),
                       Uintah::Point(8170.365578853289, 879.380308973427, 0),
                       Uintah::Point(8194.443453982514, 876.7794559275133, 0),
                       Uintah::Point(8218.430139156628, 874.103992930423, 0),
                       Uintah::Point(8242.323807696686, 871.3539973331917, 0),
                       Uintah::Point(8266.122640007321, 868.5295545893566, 0),
                       Uintah::Point(8289.824823715311, 865.6307582862573, 0),
                       Uintah::Point(8313.4285538076, 862.6577101749016, 0),
                       Uintah::Point(8336.932032768756, 859.6105201984502, 0),
                       Uintah::Point(8360.333470717858, 856.4893065193482, 0),
                       Uintah::Point(8383.631085544806, 853.2941955449974, 0),
                       Uintah::Point(8406.82310304603, 850.0253219520272, 0),
                       Uintah::Point(8429.907757059604, 846.6828287091719, 0),
                       Uintah::Point(8452.883289599751, 843.2668670986037, 0),
                       Uintah::Point(8475.747950990712, 839.7775967359133, 0),
                       Uintah::Point(8498.5, 836.2151855885398, 0),
                       Uintah::Point(8521.137703970995, 832.5798099927795, 0),
                       Uintah::Point(8543.659338954893, 828.8716546692946, 0),
                       Uintah::Point(8566.063189841994, 825.0909127370601, 0),
                       Uintah::Point(8588.347550492314, 821.237785725888, 0),
                       Uintah::Point(8610.510723865513, 817.3124835873939, 0),
                       Uintah::Point(8632.551022150126, 813.3152247044242, 0),
                       Uintah::Point(8654.46676689211, 809.2462358989718, 0),
                       Uintah::Point(8676.256289122652, 805.1057524385058, 0),
                       Uintah::Point(8697.917929485273, 800.8940180408343, 0),
                       Uintah::Point(8719.450038362193, 796.6112848772887, 0),
                       Uintah::Point(8740.850975999952, 792.2578135744628, 0),
                       Uintah::Point(8762.119112634296, 787.833873214292, 0),
                       Uintah::Point(8783.252828614279, 783.3397413325838, 0),
                       Uintah::Point(8804.250514525602, 778.7757039159796, 0),
                       Uintah::Point(8825.110571313187, 774.1420553973215, 0),
                       Uintah::Point(8845.831410402952, 769.4390986493747, 0),
                       Uintah::Point(8866.411453822771, 764.6671449770358, 0),
                       Uintah::Point(8886.849134322661, 759.8265141078828, 0),
                       Uintah::Point(8907.142895494126, 754.9175341810972, 0),
                       Uintah::Point(8927.291191888678, 749.9405417348787, 0),
                       Uintah::Point(8947.292489135541, 744.8958816921071, 0),
                       Uintah::Point(8967.145264058494, 739.783907344501, 0),
                       Uintah::Point(8986.84800479186, 734.6049803351066, 0),
                       Uintah::Point(9006.399210895652, 729.3594706391951, 0),
                       Uintah::Point(9025.797393469828, 724.0477565434834, 0),
                       Uintah::Point(9045.041075267683, 718.6702246238586, 0),
                       Uintah::Point(9064.128790808343, 713.2272697213093, 0),
                       Uintah::Point(9083.059086488369, 707.7192949163843, 0),
                       Uintah::Point(9101.830520692452, 702.1467115019586, 0),
                       Uintah::Point(9120.441663903202, 696.5099389544039, 0),
                       Uintah::Point(9138.891098810007, 690.8094049031135, 0),
                       Uintah::Point(9157.177420416971, 685.0455450984971, 0),
                       Uintah::Point(9175.2992361499, 679.2188033782194, 0),
                       Uintah::Point(9193.255165962368, 673.3296316319427, 0),
                       Uintah::Point(9211.043842440802, 667.3784897644307, 0),
                       Uintah::Point(9228.663910908615, 661.3658456570882, 0),
                       Uintah::Point(9246.114029529374, 655.2921751277269, 0),
                       Uintah::Point(9263.392869408985, 649.1579618890385, 0),
                       Uintah::Point(9280.499114696893, 642.9636975051778, 0),
                       Uintah::Point(9297.431462686291, 636.7098813469552, 0),
                       Uintah::Point(9314.188623913324, 630.3970205453642, 0),
                       Uintah::Point(9330.769322255284, 624.0256299435472, 0),
                       Uintah::Point(9347.172295027798, 617.5962320472154, 0),
                       Uintah::Point(9363.396293080987, 611.1093569734822, 0),
                       Uintah::Point(9379.440080894579, 604.5655423981918, 0),
                       Uintah::Point(9395.30243667202, 597.9653335016234, 0),
                       Uintah::Point(9410.982152433504, 591.309282912733, 0),
                       Uintah::Point(9426.478034107968, 584.5979506518349, 0),
                       Uintah::Point(9441.788901624026, 577.8319040717731, 0),
                       Uintah::Point(9456.913588999843, 571.0117177975449, 0),
                       Uintah::Point(9471.850944431913, 564.1379736645196, 0),
                       Uintah::Point(9486.599830382791, 557.2112606550443, 0),
                       Uintah::Point(9501.159123667707, 550.2321748336578, 0),
                       Uintah::Point(9515.527715540109, 543.2013192808381, 0),
                       Uintah::Point(9529.704511776095, 536.1193040252145, 0),
                       Uintah::Point(9543.688432757746, 528.9867459744211, 0),
                       Uintah::Point(9557.478413555338, 521.8042688444788, 0),
                       Uintah::Point(9571.073404008443, 514.5725030877607, 0),
                       Uintah::Point(9584.472368805902, 507.292085819547, 0),
                       Uintah::Point(9597.674287564671, 499.9636607431893, 0),
                       Uintah::Point(9610.678154907519, 492.5878780739276, 0),
                       Uintah::Point(9623.482980539606, 485.1653944612561, 0),
                       Uintah::Point(9636.087789323883, 477.6968729100446, 0),
                       Uintah::Point(9648.491621355359, 470.1829827002528, 0),
                       Uintah::Point(9660.6935320342, 462.6243993053137, 0),
                       Uintah::Point(9672.692592137668, 455.0218043092925, 0),
                       Uintah::Point(9684.487887890882, 447.3758853226534, 0),
                       Uintah::Point(9696.078521036399, 439.6873358968394, 0),
                       Uintah::Point(9707.463608902632, 431.9568554375313, 0),
                       Uintah::Point(9718.64228447106, 424.1851491166968, 0),
                       Uintah::Point(9729.613696442259, 416.3729277834264, 0),
                       Uintah::Point(9740.37700930073, 408.5209078735415, 0),
                       Uintah::Point(9750.931403378527, 400.6298113180139, 0),
                       Uintah::Point(9761.276074917678, 392.700365450232, 0),
                       Uintah::Point(9771.410236131394, 384.733302912124, 0),
                       Uintah::Point(9781.333115264062, 376.7293615591091, 0),
                       Uintah::Point(9791.043956650021, 368.6892843639942, 0),
                       Uintah::Point(9800.542020771096, 360.6138193197099, 0),
                       Uintah::Point(9809.826584312937, 352.5037193410529, 0),
                       Uintah::Point(9818.896940220082, 344.3597421652776, 0),
                       Uintah::Point(9827.752397749811, 336.1826502517581, 0),
                       Uintah::Point(9836.392282524748, 327.9732106805771, 0),
                       Uintah::Point(9844.815936584217, 319.7321950500964, 0),
                       Uintah::Point(9853.022718434346, 311.4603793736654, 0),
                       Uintah::Point(9861.012003096926, 303.1585439752506, 0),
                       Uintah::Point(9868.783182156996, 294.8274733842358, 0),
                       Uintah::Point(9876.335663809183, 286.4679562292568, 0),
                       Uintah::Point(9883.668872902772, 278.0807851311974, 0),
                       Uintah::Point(9890.782250985494, 269.6667565953075, 0),
                       Uintah::Point(9897.675256346072, 261.2266709024573, 0),
                       Uintah::Point(9904.347364055458, 252.7613319996304, 0),
                       Uintah::Point(9910.798066006817, 244.2715473895741, 0),
                       Uintah::Point(9917.026870954223, 235.7581280197154, 0),
                       Uintah::Point(9923.033304550061, 227.2218881703228, 0),
                       Uintah::Point(9928.816909381159, 218.6636453419336, 0),
                       Uintah::Point(9934.377245003619, 210.0842201421229, 0),
                       Uintah::Point(9939.713887976352, 201.484436171526, 0),
                       Uintah::Point(9944.826431893334, 192.8651199093154, 0),
                       Uintah::Point(9949.714487414556, 184.2271005979626, 0),
                       Uintah::Point(9954.377682295661, 175.5712101274424, 0),
                       Uintah::Point(9958.815661416305, 166.8982829188632, 0),
                       Uintah::Point(9963.028086807199, 158.2091558075489, 0),
                       Uintah::Point(9967.014637675837, 149.5046679255703, 0),
                       Uintah::Point(9970.775010430936, 140.7856605838179, 0),
                       Uintah::Point(9974.308918705552, 132.0529771535702, 0),
                       Uintah::Point(9977.616093378891, 123.3074629476225, 0),
                       Uintah::Point(9980.696282596793, 114.549965101013, 0),
                       Uintah::Point(9983.549251790924, 105.7813324513139, 0),
                       Uintah::Point(9986.174783696639, 97.00241541859018, 0),
                       Uintah::Point(9988.572678369512, 88.21406588499246, 0),
                       Uintah::Point(9990.742753200582, 79.41713707404178, 0),
                       Uintah::Point(9992.684842930252, 70.612483429611, 0),
                       Uintah::Point(9994.398799660867, 61.80096049465629, 0),
                       Uintah::Point(9995.884492867985, 52.98342478969321, 0),
                       Uintah::Point(9997.141809410319, 44.16073369107422, 0),
                       Uintah::Point(9998.170653538346, 35.33374530906848, 0),
                       Uintah::Point(9998.970946901598, 26.503318365791, 0),
                       Uintah::Point(9999.542628554642, 17.67031207298813, 0),
                       Uintah::Point(9999.885654961707, 8.835586009727233, 0),
                       Uintah::Point(10000, 0, 0),
                       Uintah::Point(9999.885654961707, -8.835586009727233, 0),
                       Uintah::Point(9999.542628554642, -17.67031207298813, 0)
                       };

  Polyline p_q_10000_all;
  model.computeCapPoints(3.0*10000, p_q_10000_all);
  //state.yield_f_pts = p_q_10000_all;
  state.updateYieldSurface(p_q_10000_all);

  index = 0;
  for (const auto& p_q : p_q_10000_all) {
    ASSERT_NEAR(p_q.x(), p_q_10000[index].x(), 1.0e-8);
    ASSERT_NEAR(p_q.y(), p_q_10000[index].y(), 1.0e-8);
    ++index;
  }
  ASSERT_NEAR(model.evalYieldConditionMax(&state), 923.38005468617894, 1.0e-8);
 
  /*
  for (const auto& p_q : p_q_10000_all) {
    std::cout << std::setprecision(16) 
              << "Uintah::Point(" << p_q.x() << ", " << p_q.y() << ", 0)," << "\n";
  }
  std::cout << "p_q_10000_x = ";
  for (const auto& p_q : p_q_10000_all) {
    std::cout << p_q.x() << ", ";
  }
  std::cout << "\n";
  std::cout << "p_q_10000_y = ";
  for (const auto& p_q : p_q_10000_all) {
    std::cout << p_q.y() << ", ";
  }
  std::cout << "\n";
  */
}

TEST_F(YieldCondTabularCapTest, evalYieldCondition)
{
  Matrix3 one(0.0); one.Identity();

  IntVar_TabularCap intvar(ps_intvar);
  YieldCond_TabularCap model(ps, &intvar);
  ModelState_TabularCap state;

  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;

  state.capX = -2000*3;
  Polyline p_q_2000_all;
  model.computeCapPoints(-state.capX, p_q_2000_all);
  //state.yield_f_pts = p_q_2000_all;
  state.updateYieldSurface(p_q_2000_all);

  /*
  std::cout << "p_q_2000_x = (";
  for (auto pt : state.yield_f_pts) {
    std::cout << pt.x() << ",";
  }
  std::cout << ")\n";
  std::cout << "p_q_2000_y = (";
  for (auto pt : state.yield_f_pts) {
    std::cout << pt.y() << ",";
  }
  std::cout << ")\n";
  */

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
  //state.yield_f_pts = p_q_10000_all;
  state.updateYieldSurface(p_q_10000_all);

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

void
updateClosestAndTangent(YieldCond_TabularCap& model,
                        ModelState_TabularCap& state)
{
  double sqrtKG = std::sqrt(1.5 * state.bulkModulus / state.shearModulus);

  double z = 0.0, rprime = 0.0;
  Vaango::Util::convertToZRprime(sqrtKG, -state.I1/3.0, state.sqrt_J2, z, rprime);

  double closest_z = 0.0, closest_rprime = 0.0;
  double tangent_z = 0.0, tangent_rprime = 0.0;
  model.getClosestPointAndTangent(
    &state, z, rprime, closest_z, closest_rprime, tangent_z, tangent_rprime);

  double closest_p_bar = 0.0, closest_sqrt_J2 = 0.0;
  double tangent_p_bar = 0.0, tangent_sqrt_J2 = 0.0;
  Vaango::Util::revertFromZRprime(
    sqrtKG, closest_z, closest_rprime, closest_p_bar, closest_sqrt_J2);
  Vaango::Util::revertFromZRprime(
    sqrtKG, tangent_z, tangent_rprime, tangent_p_bar, tangent_sqrt_J2);

  state.closest = Uintah::Point(closest_p_bar, closest_sqrt_J2, 0.0);
  state.tangent = Uintah::Point(tangent_p_bar, tangent_sqrt_J2, 0.0);
}

TEST_F(YieldCondTabularCapTest, df_dsigma)
{
  IntVar_TabularCap intvar(ps_intvar);
  YieldCond_TabularCap model(ps, &intvar);
  ModelState_TabularCap state;
  Matrix3 zero(0.0);
  Matrix3 one(0.0); one.Identity();
  //std::cout << model;

  state.bulkModulus = 1.0e5;
  state.shearModulus = 1.0e5;
  state.stressTensor = zero;
  state.updateStressInvariants();

  state.capX = -2000*3;

  Polyline p_q_2000_all;
  model.computeCapPoints(-state.capX, p_q_2000_all);
  //state.yield_f_pts = p_q_2000_all;
  state.updateYieldSurface(p_q_2000_all);

  // Zero everything (elastic)
  updateClosestAndTangent(model, state);
  Matrix3 df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Tension (p = 2000, J2 = 0)
  double p = 2000; 
  double sqrt_J2 = 0;
  Matrix3 s(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  Matrix3 sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);
  
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), 0.57735, 1.0e-5);

  // Tension (p = 300, J2 = 1000)
  p = 300; 
  sqrt_J2 = 1000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  updateClosestAndTangent(model, state);
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), 0.37068, 1.0e-5);
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

  updateClosestAndTangent(model, state);
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), 0.276741, 1.0e-5);
  EXPECT_NEAR(df_dsigma(0,1), 0.620581, 1.0e-5);

  // Compression (p = -2000, J2 = 4000)
  p = -2000; 
  sqrt_J2 = 4000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  updateClosestAndTangent(model, state);
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), -0.064404186516701101, 1.0e-5);
  EXPECT_NEAR(df_dsigma(0,1),  0.70269349729358022, 1.0e-5);

  // Compression (p = -3000, J2 = 0)
  p = -3000; 
  sqrt_J2 = 0;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  updateClosestAndTangent(model, state);
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), -0.57735, 1.0e-5);

  // Compression (p = -3000, J2 = 1000)
  p = -3000; 
  sqrt_J2 = 1000;
  s = Matrix3(0, sqrt_J2, 0, sqrt_J2, 0, 0, 0, 0, 0);
  sigma = s + one * p;
  state.stressTensor = sigma;
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, 3*p, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, sqrt_J2, 1.0e-8);

  updateClosestAndTangent(model, state);
  df_dsigma = model.df_dsigma(zero, &state);
  df_dsigma /= df_dsigma.Norm();
  //std::cout << "df_dsigma = " << df_dsigma << "\n";
  EXPECT_NEAR(df_dsigma(0,0), -0.47597677375342551, 1.0e-5);
  EXPECT_NEAR(df_dsigma(0,1),  0.40021140197515703, 1.0e-5);
}

TEST_F(YieldCondTabularCapTest, getClosestPoint)
{
  IntVar_TabularCap intvar(ps_intvar);
  YieldCond_TabularCap model(ps, &intvar);
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
  //state.yield_f_pts = p_q_2000_all;
  state.updateYieldSurface(p_q_2000_all);

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
  EXPECT_NEAR(z_close, -2654.6198513435979, 1.0e-8);
  ASSERT_NEAR(rprime_close, 1167.1519933578322, 1.0e-8);

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
  EXPECT_NEAR(z_close, -3295.6048589736765, 1.0e-8);
  ASSERT_NEAR(rprime_close, 678.18767889517812, 1.0e-8);

  state.stressTensor = Matrix3(-3000, 0, 0, 0, -3000, 0, 0, 0, -3000);
  state.updateStressInvariants();
  ASSERT_NEAR(state.I1, -3*3000, 1.0e-8);
  ASSERT_NEAR(state.sqrt_J2, 0, 1.0e-8);

  z = state.zz;
  rprime = state.rr*sqrtKG;
  model.getClosestPoint(&state, z, rprime, z_close, rprime_close);
  EXPECT_NEAR(z_close, -3464.1016151377544, 1.0e-8);
  ASSERT_NEAR(rprime_close, 0, 1.0e-8);
}
