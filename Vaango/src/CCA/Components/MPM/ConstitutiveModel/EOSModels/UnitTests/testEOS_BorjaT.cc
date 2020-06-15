#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Borja.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/EOS_BorjaT.h>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;

TEST(TestEOS_Borja, computePressure)
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "equation_of_state");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "borja_pressure");
  xmlDocSetRootElement(doc, rootNode);

  xmlNewChild(rootNode, nullptr, BAD_CAST "p0", BAD_CAST "-9.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "alpha", BAD_CAST "60.0");
  xmlNewChild(rootNode, nullptr, BAD_CAST "kappatilde", BAD_CAST "0.018");
  xmlNewChild(rootNode, nullptr, BAD_CAST "epse_v0", BAD_CAST "0.0");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  Uintah::ProblemSpecP ps = scinew Uintah::ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Compute pressure
  Uintah::MPMMaterial matl;
  Uintah::Matrix3 Zero(0.0);
  EOS_BorjaT eos(ps);
  EOS_BorjaT eos_copy(&eos);

  ModelState_Borja state_ten, state_com;
  auto strain_ten = Uintah::Matrix3(1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0);
  auto strain_com = Uintah::Matrix3(-1.0, -2.0, -3.0, -2.0, -4.0, -5.0, -3.0, -5.0, -6.0);
  auto dev_strains_ten = state_ten.updateStrainScalars(strain_ten, strain_ten);
  auto dev_strains_com = state_com.updateStrainScalars(strain_com, strain_com);

  double p = eos.computePressure(&matl, &state_ten, Zero, Zero, 0.0);
  //std::cout << std::setprecision(16) << "p = " << p << "\n";
  ASSERT_DOUBLE_EQ(p, -1.053658125670924e-259);
  p = eos_copy.computePressure(&matl, &state_ten, Zero, Zero, 0.0);
  //std::cout << "p = " << p << "\n";
  ASSERT_DOUBLE_EQ(p, -1.053658125670924e-259);
  p = eos.computePressure(&matl, &state_com, Zero, Zero, 0.0);
  //std::cout << "p = " << p << "\n";
  ASSERT_DOUBLE_EQ(p, -6.715316579156578e+271);

  // Compute dp/depse_v
  double deriv = eos.computeDpDepse_v(&state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 5.853656253727355e-258);
  deriv = eos.computeDpDepse_v(&state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 3.730731432864766e+273);
  
  // Compute dp/depse_s
  deriv = eos.computeDpDepse_s(&state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -2.740903121985701e-260);
  deriv = eos.computeDpDepse_s(&state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -1.746869475828551e+271);

  // Compute dp_dJ
  deriv = eos.eval_dp_dJ(&matl, 0.0, &state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 5.853656253727355e-258);
  deriv = eos.eval_dp_dJ(&matl, 0.0, &state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 3.730731432864766e+273);

  double K = eos.computeInitialBulkModulus();
  //std::cout << std::setprecision(16) << "K0 = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 500.0000000000001);
  K = eos.computeBulkModulus(&state_ten);
  //std::cout << std::setprecision(16) << "K = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 1.980553685477988e-263);
  K = eos.computeBulkModulus(&state_com);
  //std::cout << std::setprecision(16) << "K = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 1.262273281623592e+268);

  double se = eos.computeStrainEnergy(&state_ten);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 6.416993940948682e-267);
  se = eos.computeStrainEnergy(&state_com);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 4.089765432460437e+264);
  
  strain_ten = Uintah::Matrix3(1.0e-2, 2.0e-2, 3.0e-2, 
                               2.0e-2, 4.0e-2, 5.0e-2, 
                               3.0e-2, 5.0e-2, 6.0e-2);
  strain_com = Uintah::Matrix3(-1.0e-2, -2.0e-2, -3.0e-2, 
                               -2.0e-2, -4.0e-2, -5.0e-2, 
                               -3.0e-2, -5.0e-2, -6.0e-2);
  dev_strains_ten = state_ten.updateStrainScalars(strain_ten, strain_ten);
  dev_strains_com = state_com.updateStrainScalars(strain_com, strain_com);

  p = eos.computePressure(&matl, &state_ten, Zero, Zero, 0.0);
  //std::cout << std::setprecision(16) << "p = " << p << "\n";
  ASSERT_DOUBLE_EQ(p, -0.6099733486880704);
  p = eos.computePressure(&matl, &state_com, Zero, Zero, 0.0);
  //std::cout << "p = " << p << "\n";
  ASSERT_DOUBLE_EQ(p, -123980.8266421052);

  // Compute dp/depse_v
  deriv = eos.computeDpDepse_v(&state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 33.88740826044836);
  deriv = eos.computeDpDepse_v(&state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 6887823.702339181);
  
  // Compute dp/depse_s
  deriv = eos.computeDpDepse_s(&state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -15.34812179962064);
  deriv = eos.computeDpDepse_s(&state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, -3119599.950085324);

  // Compute dp_dJ
  deriv = eos.eval_dp_dJ(&matl, 0.0, &state_ten); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 33.88740826044836);
  deriv = eos.eval_dp_dJ(&matl, 0.0, &state_com); 
  //std::cout << std::setprecision(16) << "deriv = " << deriv << "\n";
  ASSERT_DOUBLE_EQ(deriv, 6887823.702339181);

  K = eos.computeInitialBulkModulus();
  //std::cout << std::setprecision(16) << "K0 = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 500.0000000000001);
  K = eos.computeBulkModulus(&state_ten);
  //std::cout << std::setprecision(16) << "K = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 1.109042452160128);
  K = eos.computeBulkModulus(&state_com);
  //std::cout << std::setprecision(16) << "K = " << K << "\n";
  ASSERT_DOUBLE_EQ(K, 225419.6848038276);

  se = eos.computeStrainEnergy(&state_ten);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se, 0.0003593297544998812);
  se = eos.computeStrainEnergy(&state_com);
  //std::cout << std::setprecision(16) << "se = " << se << "\n";
  ASSERT_DOUBLE_EQ(se,  73.03597787644013);
}
