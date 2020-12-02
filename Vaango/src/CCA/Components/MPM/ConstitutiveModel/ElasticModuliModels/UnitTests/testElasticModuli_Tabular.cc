#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_Bulk.h>
#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Tabular_BulkPressure.h>
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/YieldCondUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <fstream>
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

/*
Data:
{"Vaango_tabular_data": {
  "Meta" : {
    "title" : "Test elastic data"
  },
  "Data" : {
    "PlasticStrainVol" : [0, 0.5, 1.0],
    "Data" : [{
      "TotalStrainVol" : [0, 0.1],
      "Pressure" : [0, 1000]
    }, {
      "TotalStrainVol" : [0.5, 0.55, 0.75, 0.85 ],
      "Pressure" : [0, 1100, 1500, 2000]
    }, {
      "TotalStrainVol" : [1.0, 1.45, 1.75],
      "Pressure" : [0, 2000, 3000]
    }]
  }
}}
*/

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

    ModelState_Tabular state_init;
    state_init.elasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    state_init.plasticStrainTensor = Uintah::Matrix3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state_init);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_DOUBLE_EQ(KG.bulkModulus, 1.0e4);
    EXPECT_DOUBLE_EQ(KG.shearModulus, 7500);
    EXPECT_NEAR(dKdG.bulkModulus, -24000, 1.0e-3);
    ASSERT_NEAR(dKdG.shearModulus, -18000, 1.0e-3);
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
    EXPECT_NEAR(moduli.bulkModulus, 11439.99999999994, 1.0e-3);
    EXPECT_NEAR(moduli.shearModulus, 8579.999999999955, 1.0e-3);

    //std::cout << "K,G = " << moduli.bulkModulus << "," 
    //            << moduli.shearModulus << std::endl;
    auto moduli_derivs = model.getElasticModuliAndDerivatives(&state);
    auto KG = moduli_derivs.first;
    auto dKdG = moduli_derivs.second;
    EXPECT_NEAR(KG.bulkModulus, 11440, 1.0e-7);
    EXPECT_NEAR(KG.shearModulus, 8580, 1.0e-7);
    EXPECT_NEAR(dKdG.bulkModulus, -24000, 1.0);
    ASSERT_NEAR(dKdG.shearModulus, -18000, 1.0);

    // Compute tangent modulus
    auto tangent = model.computeElasticTangentModulus(&state);
    double K = KG.bulkModulus;
    double G = KG.shearModulus;
    double K43G = K + 4.0 * G / 3.0;
    double K23G = K - 2.0 * G / 3.0;
    //std::cout << "K = " << K << " G = " << G
    //          << "K43G = " << K43G << " K23G = " << K23G
    //          << "\nTangent = \n" << tangent << "\n";
    ASSERT_DOUBLE_EQ(tangent(1,1), K43G);
    ASSERT_DOUBLE_EQ(tangent(1,2), K23G);
    ASSERT_DOUBLE_EQ(tangent(4,4), G);

  } catch (Uintah::InvalidValue e) {
    std::cout << e.message() << std::endl;
  }
  
}

TEST(ElasticModuliTabularTest, predictorTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "DrySand_ElasticData_Ten_Com.json");
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
    throw;
  }

  try {

    // Create a model
    ElasticModuli_Tabular model(ps);

    // Create an array of total strains
    auto eps_v = Vaango::Util::linspace(0, 0.5, 100);

    /*
    std::cout << "eps_v:\n";
    std::copy(eps_v.begin(), eps_v.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create an array of plastic strains
    std::vector<double> eps_v_p  = {0.0, 0.05, 0.0628, 0.10, 0.1451, 0.20, 0.25, 
                                    0.3019, 0.32, 0.35, 0.225};

    /*
    std::cout << "eps_v_p:\n";
    std::copy(eps_v_p.begin(), eps_v_p.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create a key-value map to store elastic strains and bulk moduli
    std::map<double, std::vector<double>> elasticStrains, bulkModuli;
    std::map<double, std::vector<double>> totalStrains, pressures;

    // Loop thru plastic strains
    for (auto eps_p : eps_v_p) {

      // Compute elastic strains
      std::vector<double> eps_v_e;
      std::transform(eps_v.begin(), eps_v.end(), std::back_inserter(eps_v_e), 
                     [eps_p](double strain) -> double 
                     {
                       return (strain - eps_p);
                     });
      /*
      std::cout << "eps_v_e:\n";
      std::copy(eps_v_e.begin(), eps_v_e.end(), 
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << "\n";
      */

      // Save the elastic strains
      elasticStrains[eps_p] = eps_v_e; 

      // Save the total strains
      totalStrains[eps_p] = eps_v; 

      // Create vectors to store bulk moduli and pressures
      std::vector<double> K_vals, p_vals;

      // Loop thru elastic strains
      for (auto eps_e : eps_v_e) {

        // Create a state
        ModelState_Tabular state;
        state.elasticStrainTensor = Uintah::Matrix3(-eps_e, 0, 0, 0, 0, 0, 0, 0, 0);
        state.plasticStrainTensor = Uintah::Matrix3(-eps_p, 0, 0, 0, 0, 0, 0, 0, 0);

        // Compute moduli
        auto moduli = model.getCurrentElasticModuli(&state);

        // Save bulk modulus
        K_vals.push_back(moduli.bulkModulus);

        // Compute and save pressures
        auto pressure = model.getPressure(eps_e, eps_p);
        p_vals.push_back(pressure);
      }

      // Add element to the maps
      bulkModuli[eps_p] = K_vals;
      pressures[eps_p] = p_vals;

    }

    // Print the moduli
    /*
    for (auto K = bulkModuli.begin(); K != bulkModuli.end(); ++K) {
      std::cout << "eps_v_p = " << K->first << "\n"; 
      std::cout << "K:\n";
      std::copy(K->second.begin(), K->second.end(), 
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << "\n";
    }
    */

    // Print the moduli to a file
    std::ofstream outfile;
    outfile.open("Fox_DrySand_BulkModulus_Pred_tabular.csv");
    for (auto& eps_p : eps_v_p) {
      outfile << "Elastic Strain, Bulk modulus\n";
      outfile << eps_p << ", " << eps_p << "\n";
      auto eps_e = elasticStrains[eps_p];
      auto K = bulkModuli[eps_p];

      auto eps_iter = eps_e.begin();
      auto K_iter = K.begin(); 
      for (; K_iter != K.end(); ++eps_iter, ++K_iter) {
        outfile << *eps_iter << ", " << *K_iter << "\n";
      }
    }
    outfile.close();

    // Print the pressures to a file
    outfile.open("Fox_DrySand_Pressure_Pred_tabular.csv");
    for (auto& eps_p : eps_v_p) {
      outfile << "Total Strain, Pressure\n";
      outfile << eps_p << ", " << eps_p << "\n";
      auto eps = totalStrains[eps_p];
      auto p = pressures[eps_p];

      auto eps_iter = eps.begin();
      auto p_iter = p.begin(); 
      for (; p_iter != p.end(); ++eps_iter, ++p_iter) {
        outfile << *eps_iter << ", " << *p_iter << "\n";
      }
    }
    outfile.close();

  } catch (...) {
    throw;
  }
}

TEST(ElasticModuliTabularTest, bulkPredictorTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_bulk");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "DrySand_BulkModulusData_Ten_Com.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "PlasticStrainVol, ElasticStrainVol");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "BulkModulus");
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
    throw;
  }

  // Test model creation
  try {
    ElasticModuli_Tabular_Bulk model(ps);
  } catch (Uintah::ProblemSetupException& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::InvalidValue& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::Exception& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (...) {
    std::cout << "**ERROR** Unknown exception during model creation\n";
    throw;
  }

  try {

    // Create a model
    ElasticModuli_Tabular_Bulk model(ps);

    // Create an array of total strains
    auto eps_v = Vaango::Util::linspace(0, 0.5, 100);

    /*
    std::cout << "eps_v:\n";
    std::copy(eps_v.begin(), eps_v.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create an array of plastic strains
    std::vector<double> eps_v_p  = {0.0, 0.05, 0.0628, 0.10, 0.1451, 0.20, 0.25, 
                                    0.3019, 0.32, 0.35, 0.225};

    /*
    std::cout << "eps_v_p:\n";
    std::copy(eps_v_p.begin(), eps_v_p.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create a key-value map to store elastic strains and bulk moduli
    std::map<double, std::vector<double>> elasticStrains, bulkModuli;

    // Loop thru plastic strains
    for (auto eps_p : eps_v_p) {

      // Compute elastic strains
      std::vector<double> eps_v_e;
      std::transform(eps_v.begin(), eps_v.end(), std::back_inserter(eps_v_e), 
                     [eps_p](double strain) -> double 
                     {
                       return (strain - eps_p);
                     });
      /*
      std::cout << "eps_v_e:\n";
      std::copy(eps_v_e.begin(), eps_v_e.end(), 
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << "\n";
      */

      // Save the elastic strains
      elasticStrains[eps_p] = eps_v_e; 

      // Create vectors to store bulk moduli and pressures
      std::vector<double> K_vals;

      // Loop thru elastic strains
      for (auto eps_e : eps_v_e) {

        // Create a state
        ModelState_Tabular state;
        state.elasticStrainTensor = Uintah::Matrix3(-eps_e, 0, 0, 0, 0, 0, 0, 0, 0);
        state.plasticStrainTensor = Uintah::Matrix3(-eps_p, 0, 0, 0, 0, 0, 0, 0, 0);

        // Compute moduli
        auto moduli = model.getCurrentElasticModuli(&state);

        // Save bulk modulus
        K_vals.push_back(moduli.bulkModulus);
      }

      // Add element to the maps
      bulkModuli[eps_p] = K_vals;
    }

    // Print the moduli
    /*
    for (auto K = bulkModuli.begin(); K != bulkModuli.end(); ++K) {
      std::cout << "eps_v_p = " << K->first << "\n"; 
      std::cout << "K:\n";
      std::copy(K->second.begin(), K->second.end(), 
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << "\n";
    }
    */

    // Print the moduli to a file
    std::ofstream outfile;
    outfile.open("Fox_DrySand_BulkModulus_Pred_tabular_bulk.csv");
    for (auto& eps_p : eps_v_p) {
      outfile << "Elastic Strain, Bulk modulus\n";
      outfile << eps_p << ", " << eps_p << "\n";
      auto eps_e = elasticStrains[eps_p];
      auto K = bulkModuli[eps_p];

      auto eps_iter = eps_e.begin();
      auto K_iter = K.begin(); 
      for (; K_iter != K.end(); ++eps_iter, ++K_iter) {
        outfile << *eps_iter << ", " << *K_iter << "\n";
      }
    }
    outfile.close();

  } catch (...) {
    throw;
  }
}

TEST(ElasticModuliTabularTest, bulkPressurePredictorTest)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "tabular_bulk_pressure");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, nullptr, BAD_CAST "filename", 
              BAD_CAST "DrySand_BulkModulusPressureData_Ten_Com.json");
  xmlNewChild(rootNode, nullptr, BAD_CAST "independent_variables", 
              BAD_CAST "PlasticStrainVol, Pressure");
  xmlNewChild(rootNode, nullptr, BAD_CAST "dependent_variables", 
              BAD_CAST "BulkModulus");
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
    throw;
  }

  // Test model creation
  try {
    ElasticModuli_Tabular_BulkPressure model(ps);
  } catch (Uintah::ProblemSetupException& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::InvalidValue& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (Uintah::Exception& e) {
    std::cout << e.message() << "\n";
    throw;
  } catch (...) {
    std::cout << "**ERROR** Unknown exception during model creation\n";
    throw;
  }

  try {

    // Create a model
    ElasticModuli_Tabular_BulkPressure model(ps);

    // Create an array of pressures
    auto p_bar_array = Vaango::Util::linspace(0, 10, 100);
    std::transform(p_bar_array.begin(), p_bar_array.end(), p_bar_array.begin(), 
                     [](double p_bar) -> double 
                     {
                       return std::pow(10, p_bar);
                     });

    /*
    std::cout << "p_bar_array:\n";
    std::copy(p_bar_array.begin(), p_bar_array.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create an array of plastic strains
    std::vector<double> eps_v_p  = {0.0, 0.05, 0.0628, 0.10, 0.1451, 0.20, 0.25, 
                                    0.3019, 0.32, 0.35, 0.225};

    /*
    std::cout << "eps_v_p:\n";
    std::copy(eps_v_p.begin(), eps_v_p.end(), 
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";
    */

    // Create a key-value map to store elastic strains and bulk moduli
    std::map<double, std::vector<double>> pressures, bulkModuli;

    // Loop thru plastic strains
    for (auto eps_p : eps_v_p) {

      // Save the elastic strains
      pressures[eps_p] = p_bar_array; 

      // Create vectors to store bulk moduli and pressures
      std::vector<double> K_vals;

      // Loop thru pressures
      for (auto p_bar : p_bar_array) {

        // Create a state
        ModelState_Tabular state;
        state.stressTensor = Uintah::Matrix3(-p_bar, 0, 0, 0, -p_bar, 0, 0, 0, -p_bar);
        state.plasticStrainTensor = Uintah::Matrix3(-eps_p, 0, 0, 0, 0, 0, 0, 0, 0);

        // Compute moduli
        auto moduli = model.getCurrentElasticModuli(&state);

        // Save bulk modulus
        K_vals.push_back(moduli.bulkModulus);
      }

      // Add element to the maps
      bulkModuli[eps_p] = K_vals;
    }

    // Print the moduli
    /*
    for (auto K = bulkModuli.begin(); K != bulkModuli.end(); ++K) {
      std::cout << "eps_v_p = " << K->first << "\n"; 
      std::cout << "K:\n";
      std::copy(K->second.begin(), K->second.end(), 
                std::ostream_iterator<double>(std::cout, " "));
      std::cout << "\n";
    }
    */


    // Print the moduli to a file
    std::ofstream outfile;
    outfile.open("Fox_DrySand_BulkModulus_Pred_tabular_bulk_pressure.csv");
    for (auto& eps_p : eps_v_p) {
      outfile << "Pressure, Bulk modulus\n";
      outfile << eps_p << ", " << eps_p << "\n";
      auto p_bar = pressures[eps_p];
      auto K = bulkModuli[eps_p];

      auto p_iter = p_bar.begin();
      auto K_iter = K.begin(); 
      for (; K_iter != K.end(); ++p_iter, ++K_iter) {
        outfile << *p_iter << ", " << *K_iter << "\n";
      }
    }
    outfile.close();

  } catch (...) {
    throw;
  }
}
