#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Arenisca3.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include <libxml/parser.h>
#include <libxml/tree.h>

using namespace Vaango;
using Uintah::ProblemSpec;
using Uintah::ProblemSpecP;



int main()
{
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(NULL, BAD_CAST "elastic_moduli_model");
  xmlNewProp(rootNode, BAD_CAST "type", BAD_CAST "mason_sand");
  xmlDocSetRootElement(doc, rootNode);

  // Create a child node
  xmlNewChild(rootNode, NULL, BAD_CAST "b0",
              BAD_CAST "0.001");
  xmlNewChild(rootNode, NULL, BAD_CAST "b1",
              BAD_CAST "0.051");
  xmlNewChild(rootNode, NULL, BAD_CAST "b2",
              BAD_CAST "2.094");
  xmlNewChild(rootNode, NULL, BAD_CAST "alpha0",
              BAD_CAST "-19.705e6");
  xmlNewChild(rootNode, NULL, BAD_CAST "alpha1",
              BAD_CAST "149.814e6");
  xmlNewChild(rootNode, NULL, BAD_CAST "alpha2",
              BAD_CAST "0.006");
  xmlNewChild(rootNode, NULL, BAD_CAST "alpha3",
              BAD_CAST "7.570");
  xmlNewChild(rootNode, NULL, BAD_CAST "alpha4",
              BAD_CAST "12.734");
  xmlNewChild(rootNode, NULL, BAD_CAST "G0",
              BAD_CAST "1.0e8");
  xmlNewChild(rootNode, NULL, BAD_CAST "G1",
              BAD_CAST "0.35");
  xmlNewChild(rootNode, NULL, BAD_CAST "G2",
              BAD_CAST "-0.35");
  xmlNewChild(rootNode, NULL, BAD_CAST "G3",
              BAD_CAST "0.0");
  xmlNewChild(rootNode, NULL, BAD_CAST "G4",
              BAD_CAST "0.0");

  // Print the document to stdout
  xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Run thru the ProblemSpec machinery
  xmlNode* d_node;
  d_node = const_cast<xmlNode*>(xmlDocGetRootElement(doc));
  if (d_node == 0) {
    std::cout << "**Error** Could not create node." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }
  const xmlNode* d_child = d_node->children;
  while (d_child != 0) {
    std::string child_name((const char *)(d_child->name));
    std::cout << "child_name = " << child_name << std::endl;
    d_child = d_child->next;
  } 

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Constructors
  ElasticModuli_MasonSand model(ps);

  // Free the document (*WARNING** Doc seems to be getting freed elasewhere)
  //xmlFreeDoc(doc);
  
  
  
  
  /*
  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  std::cout << "params[Ks] = " << params["Ks"]  << " Pa" << std::endl;

  // Set up densities
  double rho_orig  = 1.0;
  std::vector<double> rho_cur;
  for (int ii = -2; ii < 10; ii++) {
    double rho = rho_orig*(1.0 + (double) ii/10.0);
    rho_cur.emplace_back(rho);
    std::cout << "rho = " << rho << std::endl;
  }

  // Compute the pressure
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);  
    std::cout << "p = " << pp << std::endl;
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (int ii = -2; ii < 12; ii++) {
    pressures.emplace_back((double) ii);
  }

  // Convert to Pa and Compute bulk modulus
  std::cout << "Compression:" << std::endl;
  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  ModelState_Arenisca3 state; 
  for (double pp : pressures) {
    state.I1 = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1 = " << state.I1 << " K = " << K << std::endl;
  }
  
  // Test tension states
  // Convert to Pa and Compute bulk modulus
  std::cout << "Tension:" << std::endl;
  for (double pp : pressures) {
    pp = -std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  for (double pp : pressures) {
    state.I1 = std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1 = " << state.I1 << " K = " << K << std::endl;
  }
  */
  
}
