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
  std::cout << docJSON1D << std::endl;


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
}

// Create a 4D test JSON document
// {"Vaango_tabular_data": {
//    "Meta" : {
//      "title" : "Test data"
//    },
//    "Data" : {
//      "Salinity" : [0.1, 0.2],
//      "Data" : {
//        "Temperature" : [100, 200, 300],
//        "Data" : {
//          "Volume" : [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8],
//          "Pressure" : [10 20 30 40 50 60 70 80]
//        }, {
//          "Volume" : [0.15 0.25 0.35],
//          "Pressure" : [100 200 300]
//        }, {
//          "Volume" : [0.05 0.45 0.75],
//          "Pressure" : [1000 2000 3000]
//        }
//      },{
//        "Temperature" : [0, 400],
//        "Data" : {
//          "Volume" : [0.1 0.2 0.3 0.4],
//          "Pressure" : [15 25 35 45]
//        }, {
//          "Volume" : [0.1 0.45 0.65],
//          "Pressure" : [150 250 350]
//        }
//      }
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

  json docJSON4D = {
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
            {"Pressure" , {10, 20, 30, 40, 50, 60, 70, 80}}
          }, {
            {"Volume" , {0.15, 0.25, 0.35}},
            {"Pressure" , {100, 200, 300}}
          }, {
            {"Volume" , {0.05, 0.45, 0.75}},
            {"Pressure" , {1000, 2000, 3000}}
          }}}
        }, {
          {"Temperature" , {0, 400}},
          {"Data" , {{
            {"Volume" , {0.1, 0.2, 0.3, 0.4}},
            {"Pressure" , {15, 25, 35, 45}}
          }, {
            {"Volume" , {0.1, 0.45, 0.65}},
            {"Pressure" , {150, 250, 350}}
          }}}
        }}}
      }}
   } }
  };
  //std::cout << docJSON4D;

  try {
    eos.readJSONTable<4>(docJSON4D, "test_dummy");
  } catch (ProblemSetupException e) {
    std::cout << e.message() << std::endl;
  }

}