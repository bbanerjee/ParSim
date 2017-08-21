#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>

#include <Core/Malloc/Allocator.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

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