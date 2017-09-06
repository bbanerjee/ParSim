#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>

#include <Core/Parallel/Parallel.h>
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

class VaangoEnv : public ::testing::Environment {
public:

  int d_argc;
  char** d_argv;

  explicit VaangoEnv(int argc, char**argv) {
    d_argc = argc;
    d_argv = argv;
  }

  virtual ~VaangoEnv() {}

  virtual void SetUp() {
    Uintah::Parallel::determineIfRunningUnderMPI(d_argc, d_argv);
    Uintah::Parallel::initializeManager(d_argc, d_argv);



  }

  virtual void TearDown() {
    Uintah::Parallel::finalizeManager();

  }

  ProblemSpecP createInput() {

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
    //xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

    // Create a ProblemSpec
    ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
    if (!ps) {
      std::cout << "**Error** Could not create ProblemSpec." << std::endl;
      std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      exit(-1);
    }

    return ps;
  }
};

int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new VaangoEnv(argc, argv));
  return RUN_ALL_TESTS();
}

TEST(TabularDataTest, parseVariableNames)
{
}
