#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulus_BorjaT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ShearModulusModels/ShearModulus_ConstantT.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EOS_BorjaT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_BorjaT.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_MetalT.h>

#include <iostream>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <gtest/gtest.h>

TEST(TestShearModulus_ConstantT, Constructors)
{
  // Create a new document
  xmlDocPtr shear_doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr shearNode = xmlNewNode(nullptr, BAD_CAST "shear_modulus_model");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "constant_shear");
  xmlDocSetRootElement(shear_doc, shearNode);

  xmlNewChild(shearNode, nullptr, BAD_CAST "shear_modulus", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", shear_doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  Uintah::ProblemSpecP shear_ps = scinew Uintah::ProblemSpec(xmlDocGetRootElement(shear_doc), false);
  if (!shear_ps) {
    std::cout << "**Error** Could not create Shear Modulus ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create an shear modulus model
  Vaango::ShearModulus_ConstantT model(shear_ps);
}

TEST(TestShearModulus_BorjaT, Constructors)
{
  // Create a new document
  xmlDocPtr eos_doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr eosNode = xmlNewNode(nullptr, BAD_CAST "equation_of_state");
  xmlNewProp(eosNode, BAD_CAST "type", BAD_CAST "borja_pressure");
  xmlDocSetRootElement(eos_doc, eosNode);

  xmlNewChild(eosNode, nullptr, BAD_CAST "p0", BAD_CAST "-9.0");
  xmlNewChild(eosNode, nullptr, BAD_CAST "alpha", BAD_CAST "60.0");
  xmlNewChild(eosNode, nullptr, BAD_CAST "kappatilde", BAD_CAST "0.018");
  xmlNewChild(eosNode, nullptr, BAD_CAST "epse_v0", BAD_CAST "0.0");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", eos_doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  Uintah::ProblemSpecP eos_ps = scinew Uintah::ProblemSpec(xmlDocGetRootElement(eos_doc), false);
  if (!eos_ps) {
    std::cout << "**Error** Could not create EOS ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a new document
  xmlDocPtr shear_doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr shearNode = xmlNewNode(nullptr, BAD_CAST "shear_modulus_model");
  xmlNewProp(shearNode, BAD_CAST "type", BAD_CAST "borja_shear");
  xmlDocSetRootElement(shear_doc, shearNode);

  xmlNewChild(shearNode, nullptr, BAD_CAST "mu0", BAD_CAST "5.4e6");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", shear_doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  Uintah::ProblemSpecP shear_ps = scinew Uintah::ProblemSpec(xmlDocGetRootElement(shear_doc), false);
  if (!shear_ps) {
    std::cout << "**Error** Could not create Shear Modulus ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create an EOS first
  Vaango::EOS_BorjaT eos(eos_ps);

  // Create an shear modulus model
  Vaango::ShearModulus_BorjaT model(shear_ps, &eos);

  // Check parameters
  auto params = model.getParameters();
  //for (auto param : params) {
  //  std::cout << "param: " << param.first << " value:" << param.second << "\n";
  //}
  ASSERT_EQ(params["kappatilde"], 0.018);
  ASSERT_EQ(params["mu0"], 5.4e6);
  //ASSERT_DOUBLE_EQ(state.pressure, state_epcopy.pressure);

  // Create a copy
  Vaango::ShearModulus_BorjaT model_copy(&model);
  params = model_copy.getParameters();
  ASSERT_EQ(params["kappatilde"], 0.018);
  ASSERT_EQ(params["mu0"], 5.4e6);

  // Compute shear modulus
  double mu = model.computeInitialShearModulus();
  //std::cout << std::setprecision(16) << "mu = " << mu << "\n";
  ASSERT_DOUBLE_EQ(mu, 5400540);

  // Create high model state
  Vaango::ModelState_BorjaT state_ten, state_com;
  auto strain_ten = Uintah::Matrix3(1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0);
  auto strain_com = Uintah::Matrix3(-1.0, -2.0, -3.0, -2.0, -4.0, -5.0, -3.0, -5.0, -6.0);
  auto dev_strains_ten = state_ten.updateStrainScalars(strain_ten, strain_ten);
  auto dev_strains_com = state_com.updateStrainScalars(strain_com, strain_com);

  mu = model.computeShearModulus(&state_ten);
  //std::cout << std::setprecision(16) << "mu = " << mu << "\n";
  ASSERT_DOUBLE_EQ(mu, 5400000);
  mu = model.computeShearModulus(&state_com);
  //std::cout << std::setprecision(16) << "mu = " << mu << "\n";
  ASSERT_DOUBLE_EQ(mu, 1.363255144153479e+268);

  double se = model.computeStrainEnergy(&state_ten);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 478800000);
  se = model.computeStrainEnergy(&state_com);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 1.208752894482752e+270);

  double q = model.computeQ(&state_ten);
  //std::cout << std::setprecision(16) << "q = " << q << "\n";
  ASSERT_DOUBLE_EQ(q, 124551676.0224446);
  q = model.computeQ(&state_com);
  //std::cout << std::setprecision(16) << "q = " << q << "\n";
  ASSERT_DOUBLE_EQ(q, 3.144365056491392e+269);

  double deriv = model.computeDqDepse_v(&state_ten);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -2.740903121985701e-260);
  deriv = model.computeDqDepse_v(&state_com);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -1.746869475828551e+271);
  
  deriv = model.computeDqDepse_s(&state_ten);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 16200000);
  deriv = model.computeDqDepse_s(&state_com);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 4.089765432460438e+268);

  strain_ten = Uintah::Matrix3(1.0e-2, 2.0e-2, 3.0e-2, 2.0e-2, 4.0e-2, 5.0e-2, 3.0e-2, 5.0e-2, 6.0e-2);
  strain_com = Uintah::Matrix3(-1.0e-2, -2.0e-2, -3.0e-2, -2.0e-2, -4.0e-2, -5.0e-2, -3.0e-2, -5.0e-2, -6.0e-2);
  dev_strains_ten = state_ten.updateStrainScalars(strain_ten, strain_ten);
  dev_strains_com = state_com.updateStrainScalars(strain_com, strain_com);

  mu = model.computeShearModulus(&state_ten);
  //std::cout << std::setprecision(16) << "mu = " << mu << "\n";
  ASSERT_DOUBLE_EQ(mu, 5400001.197765849);
  mu = model.computeShearModulus(&state_com);
  //std::cout << std::setprecision(16) << "mu = " << mu << "\n";
  ASSERT_DOUBLE_EQ(mu, 5643453.259588134);

  se = model.computeStrainEnergy(&state_ten);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 47880.01062019054);
  se = model.computeStrainEnergy(&state_com);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 50038.61890168147);

  q = model.computeQ(&state_ten);
  //std::cout << std::setprecision(16) << "q = " << q << "\n";
  ASSERT_DOUBLE_EQ(q, 1245517.036490639);
  q = model.computeQ(&state_com);
  //std::cout << std::setprecision(16) << "q = " << q << "\n";
  ASSERT_DOUBLE_EQ(q, 1301669.559325982);

  deriv = model.computeDqDepse_v(&state_ten);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -15.34812179962064);
  deriv = model.computeDqDepse_v(&state_com);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -3119599.950085324);
  
  deriv = model.computeDqDepse_s(&state_ten);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 16200003.59329754);
  deriv = model.computeDqDepse_s(&state_com);
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 16930359.7787644);
}

