#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>

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

TEST(TabularDataTest, parseVariableNames)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "table_eos");
  xmlNewProp(rootNode, BAD_CAST "interpolation", BAD_CAST "linear");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "table_eos.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", BAD_CAST "temperature, density, other 1, other 2, other 3");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", BAD_CAST "pressure, volume 1, other 2");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TabularData eos(ps);

}

// Create a 1D test JSON document
// {"Vaango_tabular_data": {
//   "Meta" : {
//     "title" : "1D Test data"
//   },
//   "Data" : {
//     "Volume" : [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8],
//     "Pressure" : [10 20 30 40 50 60 70 80],
//     "Density" : [1.1 2.1 3.1 4.1 5.1 6.1 7.1 8.1]
//   }
// }
TEST(TabularDataTest, readJSONTableFromStream1D)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "constitutive_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_eos");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "table_eos.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "Volume");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure, Density");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TabularData eos(ps);

  json docJSON1D = {
    {"Vaango_tabular_data", {
      {"Meta" , {
        {"title" , "Test data"}
      }},
      {"Data" , {
        {"Volume" , {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}},
        {"Pressure" , {10, 20, 30, 40, 50, 60, 70, 80}},
        {"Density" , {1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1}}
      }}
   }}
  };
  //std::cout << docJSON1D << std::endl;


  try {
    eos.readJSONTable<1>(docJSON1D, "test_dummy");
  } catch (ProblemSetupException e) {
    std::cout << e.message() << std::endl;
  }

  EXPECT_EQ(eos.getNumIndependents(), 1u);
  EXPECT_EQ(eos.getNumDependents(), 2u);

  auto indepData = eos.getIndependentVarData("Volume", TableContainers::IndexKey(0, 0, 0, 0));
  auto depData = eos.getDependentVarData("Pressure", TableContainers::IndexKey(0, 0, 0, 0));

  auto val = eos.interpolateLinearSpline<1>({{0.1}}, 
             eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 10);
  EXPECT_DOUBLE_EQ(val[1], 1.1);
  val = eos.interpolateLinearSpline<1>({{0.8}}, 
        eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 80);
  val = eos.interpolateLinearSpline<1>({{0.625}}, 
        eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 62.5);
  EXPECT_DOUBLE_EQ(val[1], 6.35);

  val = eos.interpolateLinearSpline<1>({{-0.1}}, 
        eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], -10);
  EXPECT_NEAR(val[1], -0.8999, 1.0e-4);

  //EXPECT_THROW(eos.interpolateLinearSpline<1>({{-0.1}}, 
  //      eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  val = eos.interpolateLinearSpline<1>({{0.9}}, 
        eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_NEAR(val[0], 89.9999, 1.0e-4);
  EXPECT_NEAR(val[1], 9.0999, 1.0e-4);

  //EXPECT_THROW(eos.interpolateLinearSpline<1>({{0.9}}, 
  //      eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  //---------------------------------------------------
  // Translation test
  //---------------------------------------------------
  std::vector<Uintah::Vector> normals;
  for (auto val : indepData) {
    normals.push_back(Uintah::Vector(val,0,0));
  }
  eos.translateAlongNormals<1>(normals, -0.1);
  auto data = eos.getIndependentVarData("Volume", 
                                        TableContainers::IndexKey(0, 0, 0, 0));
  EXPECT_DOUBLE_EQ(data[0], 0.09);
  
  //std::cout << "Volume:" ;
  //std::copy(data.begin(), data.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;
  
}

// Create a 2D test JSON document
// {"Vaango_tabular_data": {
//    "Meta" : {
//      "title" : "Test data"
//    },
//    "Data" : {
//      "Temperature" : [100, 200, 300],
//      "Data" : [{
//        "Volume" : [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8],
//        "Pressure" : [10 20 30 40 50 60 70 80]
//      }, {
//        "Volume" : [0.15 0.25 0.35],
//        "Pressure" : [100 200 300]
//      }, {
//        "Volume" : [0.05 0.45 0.75],
//        "Pressure" : [1000 2000 3000]
//      }]
//    }
//  }}
TEST(TabularDataTest, readJSONTableFromStream2D)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "table_eos");
  xmlNewProp(rootNode, BAD_CAST "interpolation", BAD_CAST "linear");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "table_eos.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "Temperature, Volume");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TabularData eos(ps);

  json docJSON2D = {
    {"Vaango_tabular_data", {
      {"Meta" , {
        {"title" , "Test data"}
      }},
      {"Data" , {
        {"Temperature" , {100, 200, 300}},
        {"Data" , {{
          {"Volume" , {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}},
          {"Pressure" , {10, 20, 30, 40, 50, 60, 70, 80}}
        }, {
          {"Volume" , {0.15, 0.25, 0.35}},
          {"Pressure" , {100, 200, 300}}
        }, {
          {"Volume" , {0.05, 0.45, 0.75}},
          {"Pressure" , {1000, 2000, 3000}}
        }}}
      }}
    }}
  };
  //std::cout << docJSON2D;

  try {
    eos.readJSONTable<2>(docJSON2D, "test_dummy");
  } catch (ProblemSetupException e) {
    std::cout << e.message() << std::endl;
  }

  EXPECT_EQ(eos.getNumIndependents(), 2u);
  EXPECT_EQ(eos.getNumDependents(), 1u);

  auto tempData = eos.getIndependentVarData("Temperature", 
                                             TableContainers::IndexKey(0, 0, 0, 0));
  std::vector<std::vector<double>> indepData;
  std::vector<std::vector<double>> depData;
  for (auto ii = 0u; ii < tempData.size(); ii++) {
    indepData.push_back(eos.getIndependentVarData("Volume", 
                                             TableContainers::IndexKey(ii, 0, 0, 0)));
    depData.push_back(eos.getDependentVarData("Pressure", 
                                         TableContainers::IndexKey(ii, 0, 0, 0)));
    /*
    std::cout << "Volume[" << ii << "]";
    std::copy(indepData[ii].begin(), indepData[ii].end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Pressure[" << ii << "]";
    std::copy(depData[ii].begin(), depData[ii].end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
  }
                                          
  //---------------------------------------------------
  // Interpolation test
  //---------------------------------------------------
  std::array<double, 2> indepVals = {{150, 0.25}};
  auto val = eos.interpolateLinearSpline<2>(indepVals, 
             eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 112.5);
  indepVals = {{20, 0.25}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], -115);
  //ASSERT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{500, 0.25}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 4100);
  //EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.01}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 148);
  //EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.4}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 655);
  //EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.349}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 588.7);

  //---------------------------------------------------
  // Translation test
  //---------------------------------------------------
  eos.translateIndepVar1ByIndepVar0<2>();
  for (auto ii = 0u; ii < tempData.size(); ii++) {
    auto data = eos.getIndependentVarData("Volume", 
                                          TableContainers::IndexKey(ii, 0, 0, 0));
    if (ii == 0) {
      EXPECT_DOUBLE_EQ(data[0], -99.9);
    } else if (ii == 1) {
      EXPECT_DOUBLE_EQ(data[1], -199.75);
    } else if (ii == 2) {
      EXPECT_DOUBLE_EQ(data[2], -299.25);
    }
    /*
    std::cout << "Volume[" << ii << "]";
    std::copy(data.begin(), data.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
  }

}

// Create a 3D test JSON document
// {"Vaango_tabular_data": {
//    "Meta" : {
//      "title" : "Test data"
//    },
//    "Data" : {
//      "Salinity" : [0.1, 0.2],
//      "Data" : [{
//        "Temperature" : [100, 200, 300],
//        "Data" : [{
//          "Volume" : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
//          "Pressure" : [10, 20, 30, 40, 50, 60, 70, 80]
//        }, {
//          "Volume" : [0.15, 0.25, 0.35],
//          "Pressure" : [100, 200, 300]
//        }, {
//          "Volume" : [0.05, 0.45, 0.75],
//          "Pressure" : [1000, 2000, 3000]
//        }]
//      },{
//        "Temperature" : [0, 400],
//        "Data" : [{
//          "Volume" : [0.1, 0.2, 0.3, 0.4],
//          "Pressure" : [15, 25, 35, 45]
//        }, {
//          "Volume" : [0.1, 0.45, 0.65],
//          "Pressure" : [150, 250, 350]
//        }]
//      }]
//    }
//  }
TEST(TabularDataTest, readJSONTableFromStream3D)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "table_eos");
  xmlNewProp(rootNode, BAD_CAST "interpolation", BAD_CAST "linear");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "table_eos.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "Salinity, Temperature, Volume");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure, Density");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TabularData eos(ps);

  json docJSON3D = {
    {"Vaango_tabular_data", {
      {"Meta" , {
        {"title" , "Test data"}
      }},
      {"Data" , {
        {"Salinity" , {0.1, 0.2}},
        {"Data" , {{
          {"Temperature" , {100, 200, 300}},
          {"Data" , {{
            {"Volume" , {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}},
            {"Pressure" , {10, 20, 30, 40, 50, 60, 70, 80}},
            {"Density" , {1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80}}
          }, {
            {"Volume" , {0.15, 0.25, 0.35, 0.45, 0.55}},
            {"Pressure" , {100, 200, 300, 400, 500}},
            {"Density" , {2.100, 2.200, 2.300, 2.400, 2.500}}
          }, {
            {"Volume" , {0.05, 0.75}},
            {"Pressure" , {1000, 2000}},
            {"Density" , {3.1000, 3.2000}}
          }}}
        }, {
          {"Temperature" , {0, 400}},
          {"Data" , {{
            {"Volume" , {0.1, 0.2, 0.3, 0.8}},
            {"Pressure" , {15, 25, 35, 45}},
            {"Density" , {4.15, 4.25, 4.35, 4.45}}
          }, {
            {"Volume" , {0.1, 0.45, 0.65}},
            {"Pressure" , {150, 250, 350}},
            {"Density" , {5.150, 5.250, 5.350}}
          }}}
        }}}
      }}
   } }
  };
  //std::cout << docJSON3D;

  try {
    eos.readJSONTable<3>(docJSON3D, "test_dummy");
  } catch (ProblemSetupException e) {
    std::cout << e.message() << std::endl;
  }

  EXPECT_EQ(eos.getNumIndependents(), 3u);
  EXPECT_EQ(eos.getNumDependents(), 2u);

  auto salData = eos.getIndependentVarData("Salinity", 
                                           TableContainers::IndexKey(0, 0, 0, 0));
  /*
  std::cout << "Salinity ";
  std::copy(salData.begin(), salData.end(),
            std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  */

  std::vector<std::vector<double>> indepData;
  std::vector<std::vector<double>> depData;
  auto numSal = salData.size();
  for (auto ii = 0u; ii < numSal; ii++) {
    auto tempData = eos.getIndependentVarData("Temperature", 
                                              TableContainers::IndexKey(ii, 0, 0, 0));
    /*
    std::cout << "Temperature[" << ii << "]";
    std::copy(tempData.begin(), tempData.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */

    auto numTemp = tempData.size();
    for (auto jj = 0u; jj < numTemp; jj++) {
      indepData.push_back(eos.getIndependentVarData("Volume", 
                                             TableContainers::IndexKey(ii, jj, 0, 0)));
      depData.push_back(eos.getDependentVarData("Pressure", 
                                         TableContainers::IndexKey(ii, jj, 0, 0)));
      /*
      std::cout << "Volume[" << ii << ", " << jj << "]";
      std::copy(indepData[jj].begin(), indepData[jj].end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
      std::cout << "Pressure[" << ii << ", " << jj << "]";
      std::copy(depData[jj].begin(), depData[jj].end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
      */
    }
  }

  std::array<double, 3> indepVals = {{0.12, 150, 0.25}};
  auto val = eos.interpolateLinearSpline<3>(indepVals, 
                                            eos.getIndependentVars(),
                                            eos.getDependentVars());
  EXPECT_EQ(val.size(), eos.getNumDependents());
  EXPECT_NEAR(val[0], 1.08214285714286e+02, 1.0e-10);
  ASSERT_NEAR(val[1], 2.30696428571429e+00, 1.0e-10);

  //std::cout << "vals = " << val[0] << " " << val[1] << std::endl;

  indepVals = {{0.12, 230, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], 6.81225714285714e+02, 1.0e-10);
  ASSERT_NEAR(val[1], 3.11120357142857e+00, 1.0e-10);
  //std::cout << "vals = " << val[0] << " " << val[1] << std::endl;

  indepVals = {{0.05, 230, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], 1124.4357142857143, 1.0e-10);
  ASSERT_NEAR(val[1], 1.5469910714285704, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  indepVals = {{0.3, 230, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], -458.45714285714263, 1.0e-10);
  ASSERT_NEAR(val[1], 7.133464285714286, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  indepVals = {{0.12, 301, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], 1367.1465714285714, 1.0e-10);
  ASSERT_NEAR(val[1], 3.5483353571428569, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  indepVals = {{0.12, 0, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], -272.19999999999993, 1.0e-10);
  ASSERT_NEAR(val[1], 1.3179999999999998, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  indepVals = {{0.12, 230, 0.6}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], 777.43142857142846, 1.0e-10);
  ASSERT_NEAR(val[1], 3.1780821428571424, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);

  indepVals = {{0.12, 230, 0.05}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
  EXPECT_NEAR(val[0], 256.45714285714286, 1.0e-10);
  ASSERT_NEAR(val[1], 2.8031071428571428, 1.0e-10);
  //EXPECT_THROW(eos.interpolateLinearSpline<3>(indepVals, 
  //             eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
}

TEST(TabularDataTest, findIntersectionTable)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "constitutive_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_hydrostat");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "DrySand_HydrostatData.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "TotalStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table hydrostat
  TabularData hydrostat(ps);
  hydrostat.setup();

  auto x_coords = hydrostat.getIndependentVarData("TotalStrainVol", IndexKey(0, 0, 0, 0));
  auto y_coords = hydrostat.getDependentVarData("Pressure", IndexKey(0, 0, 0, 0));

  bool status;
  std::size_t index;
  double t_p, t_q, x_q, y_q;
  std::tie(status, index, t_p, t_q, x_q, y_q) = 
    Vaango::Util::findIntersectionTableLinearSearch(x_coords, y_coords,
                                                    -1, 100, 1, 100);
  //std::cout << "status = " << status << " index = " << index
  //          << " t_p = " << t_p << " t_q = " << t_q 
  //          << " x_q = " << x_q << " y_q = " << y_q << "\n";
  ASSERT_NEAR(t_q, 0.00241924, 1.0e-8);

  std::tie(status, index, t_p, t_q, x_q, y_q) = 
    Vaango::Util::findIntersectionTableLinearSearch(x_coords, y_coords,
                                                    -1, -1.0e6, 1, -1.0e6);
  ASSERT_NEAR(t_q, -24.1924, 1.0e-4);
  //std::cout << "status = " << status << " index = " << index
  //          << " t_p = " << t_p << " t_q = " << t_q 
  //          << " x_q = " << x_q << " y_q = " << y_q << "\n";

  std::tie(status, index, t_p, t_q, x_q, y_q) = 
    Vaango::Util::findIntersectionTableLinearSearch(x_coords, y_coords,
                                                    -1, 1.0e8, 1, 1.0e8);
  ASSERT_NEAR(t_q, 0.504185, 1.0e-6);
  //std::cout << "status = " << status << " index = " << index
  //          << " t_p = " << t_p << " t_q = " << t_q 
  //          << " x_q = " << x_q << " y_q = " << y_q << "\n";

  std::tie(status, index, t_p, t_q, x_q, y_q) = 
    Vaango::Util::findIntersectionTableLinearSearch(x_coords, y_coords,
                                                    -1, 6.0e9, 1, 6.0e9);
  ASSERT_NEAR(t_q, 15.0588, 1.0e-4);
  //std::cout << "status = " << status << " index = " << index
  //          << " t_p = " << t_p << " t_q = " << t_q 
  //          << " x_q = " << x_q << " y_q = " << y_q << "\n";
}

TEST(TabularDataTest, intersect1D)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "constitutive_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_hydrostat");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "DrySand_HydrostatData.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "TotalStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table hydrostat
  TabularData hydrostat(ps);
  hydrostat.setup();

  bool status;
  double x_p, y_p;

  std::tie(status, x_p, y_p) = hydrostat.intersect1D(1, -1, 100, 100);
  //std::cout << "status = " << status 
  //          << " x_p = " << x_p << " y_p = " << y_p << "\n";
  ASSERT_NEAR(x_p, 4.83848e-07, 1.0e-10);

  std::tie(status, x_p, y_p) = hydrostat.intersect1D(1, -1, -1.0e6, -1.0e6);
  //std::cout << "status = " << status 
  //          << " x_p = " << x_p << " y_p = " << y_p << "\n";
  ASSERT_NEAR(x_p, -4.83848e-03, 1.0e-8);

  std::tie(status, x_p, y_p) = hydrostat.intersect1D(1, -1, 1.0e8, 1.0e8);
  //std::cout << "status = " << status 
  //          << " x_p = " << x_p << " y_p = " << y_p << "\n";
  ASSERT_NEAR(x_p, 0.210452, 1.0e-6);

  std::tie(status, x_p, y_p) = hydrostat.intersect1D(1, -1, 6.0e9, 6.0e9);
  //std::cout << "status = " << status 
  //          << " x_p = " << x_p << " y_p = " << y_p << "\n";
  ASSERT_NEAR(x_p, 0.466059, 1.0e-6);
}

TEST(TabularDataTest, computeDerivative)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "constitutive_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_elastic");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", BAD_CAST "DrySand_ElasticData.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "PlasticStrainVol, TotalStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure");
  auto interp = xmlNewChild(rootNode, nullptr, BAD_CAST "interpolation",
                            BAD_CAST "");
  xmlNewProp(interp, BAD_CAST "type", BAD_CAST "linear");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table hydrostat
  TabularData elastic(ps);
  elastic.setup();

  double deriv = elastic.computeDerivative(0, 0); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(-0.5, 0); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.1, 1.0e6); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.1, -1.0e6); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.2, 1.0e8); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.3, 5.0e8); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.5, 5.0e8); 
  //std::cout << "deriv = " << deriv << "\n";

  deriv = elastic.computeDerivative(0.4, 1.5e9); 
  //std::cout << "deriv = " << deriv << "\n";

}
