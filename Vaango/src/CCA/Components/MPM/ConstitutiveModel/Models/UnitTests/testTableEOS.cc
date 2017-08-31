#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>

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

TEST(TableEOSTest, parseVariableNames)
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

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TableEOS eos(ps);

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
TEST(TableEOSTest, readJSONTableFromStream1D)
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
              BAD_CAST "Volume");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "Pressure, Density");

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TableEOS eos(ps);

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

  auto indepData = eos.getIndependentVarData("Volume", TableEOS::IndexKey(0, 0, 0, 0));

  auto depData = eos.getDependentVarData("Pressure", TableEOS::IndexKey(0, 0, 0, 0));
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline1D(0.1, indepData, depData), 10);
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline1D(0.8, indepData, depData), 80);
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline1D(0.625, indepData, depData), 62.5);
  EXPECT_THROW(eos.interpolateLinearSpline1D(-0.1, indepData, depData), InvalidValue);
  EXPECT_THROW(eos.interpolateLinearSpline1D(0.9, indepData, depData), InvalidValue);

  depData = eos.getDependentVarData("Density", TableEOS::IndexKey(0, 0, 0, 0));
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline1D(0.625, indepData, depData), 6.35);

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
  EXPECT_THROW(eos.interpolateLinearSpline<1>({{-0.1}}, 
        eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  EXPECT_THROW(eos.interpolateLinearSpline<1>({{0.9}}, 
        eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
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
TEST(TableEOSTest, readJSONTableFromStream2D)
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

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TableEOS eos(ps);

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
                                             TableEOS::IndexKey(0, 0, 0, 0));
  std::vector<std::vector<double>> indepData;
  std::vector<std::vector<double>> depData;
  for (auto ii = 0u; ii < tempData.size(); ii++) {
    indepData.push_back(eos.getIndependentVarData("Volume", 
                                             TableEOS::IndexKey(ii, 0, 0, 0)));
    depData.push_back(eos.getDependentVarData("Pressure", 
                                         TableEOS::IndexKey(ii, 0, 0, 0)));
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
                                          
  std::array<double, 2> indepVals = {{150, 0.25}};
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline2D(indepVals, tempData, 
                   indepData, depData), 112.5);

  indepVals = {{20, 0.25}};
  EXPECT_THROW(eos.interpolateLinearSpline2D(indepVals, tempData, 
               indepData, depData), InvalidValue);
  /*
  try {
  val = eos.interpolateLinearSpline2D(indepVals, tempData, indepData, depData);
  } catch (InvalidValue e) {
    std::cout << e.message() << std::endl;
  }
  std::cout << "val = " << val << " input: " << indepVals[0] << "," << indepVals[1] << std::endl;
  */

  indepVals = {{500, 0.25}};
  EXPECT_THROW(eos.interpolateLinearSpline2D(indepVals, tempData, 
               indepData, depData), InvalidValue);

  indepVals = {{220, 0.01}};
  EXPECT_THROW(eos.interpolateLinearSpline2D(indepVals, tempData, 
               indepData, depData), InvalidValue);

  indepVals = {{220, 0.4}};
  EXPECT_THROW(eos.interpolateLinearSpline2D(indepVals, tempData, 
               indepData, depData), InvalidValue);

  indepVals = {{220, 0.349}};
  EXPECT_DOUBLE_EQ(eos.interpolateLinearSpline2D(indepVals, tempData, 
                   indepData, depData), 588.7);

  indepVals = {{150, 0.25}};
  auto val = eos.interpolateLinearSpline<2>(indepVals, 
             eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 112.5);
  indepVals = {{20, 0.25}};
  EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{500, 0.25}};
  EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.01}};
  EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.4}};
  EXPECT_THROW(eos.interpolateLinearSpline<2>(indepVals, 
               eos.getIndependentVars(), eos.getDependentVars()), InvalidValue);
  indepVals = {{220, 0.349}};
  val = eos.interpolateLinearSpline<2>(indepVals, 
             eos.getIndependentVars(), eos.getDependentVars());
  EXPECT_DOUBLE_EQ(val[0], 588.7);
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
TEST(TableEOSTest, readJSONTableFromStream4D)
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

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a table eos
  TableEOS eos(ps);

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
                                           TableEOS::IndexKey(0, 0, 0, 0));
  std::cout << "Salinity ";
  std::copy(salData.begin(), salData.end(),
            std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;

  std::vector<std::vector<double>> indepData;
  std::vector<std::vector<double>> depData;
  auto numSal = salData.size();
  for (auto ii = 0u; ii < numSal; ii++) {
    auto tempData = eos.getIndependentVarData("Temperature", 
                                              TableEOS::IndexKey(ii, 0, 0, 0));
    std::cout << "Temperature[" << ii << "]";
    std::copy(tempData.begin(), tempData.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    auto numTemp = tempData.size();
    for (auto jj = 0u; jj < numTemp; jj++) {
      indepData.push_back(eos.getIndependentVarData("Volume", 
                                             TableEOS::IndexKey(ii, jj, 0, 0)));
      depData.push_back(eos.getDependentVarData("Pressure", 
                                         TableEOS::IndexKey(ii, jj, 0, 0)));
      std::cout << "Volume[" << ii << ", " << jj << "]";
      std::copy(indepData[jj].begin(), indepData[jj].end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
      std::cout << "Pressure[" << ii << ", " << jj << "]";
      std::copy(depData[jj].begin(), depData[jj].end(),
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << std::endl;
    }
  }

  std::array<double, 3> indepVals = {{0.12, 150, 0.25}};
  auto val = eos.interpolateLinearSpline<3>(indepVals, 
                                            eos.getIndependentVars(),
                                            eos.getDependentVars());

  indepVals = {{0.12, 230, 0.5}};
  val = eos.interpolateLinearSpline<3>(indepVals, 
                                       eos.getIndependentVars(),
                                       eos.getDependentVars());
}