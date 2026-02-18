#include <CCA/Components/DataArchiver/DataArchiver.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Exceptions/Exception.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <memory>

using namespace Uintah;

class DataArchiverTest : public ::testing::Test {
protected:
  virtual void SetUp() override {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testDataArchiver", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
    world = Uintah::Parallel::getRootProcessorGroup();
    mat_manager = std::make_shared<MaterialManager>();
    mat_manager->finalizeMaterials();
    data_archiver = std::make_unique<DataArchiver>(world, -1);
  }

  const ProcessorGroup* world;
  MaterialManagerP mat_manager;
  std::unique_ptr<DataArchiver> data_archiver;
};

TEST_F(DataArchiverTest, ConstructorTest) {
  EXPECT_EQ(data_archiver->getOutputInterval(), 0.0);
  EXPECT_EQ(data_archiver->getOutputTimestepInterval(), 0);
}

TEST_F(DataArchiverTest, ProblemSetupTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr daNode = xmlNewChild(root, nullptr, BAD_CAST "DataArchiver", nullptr);
  xmlNewChild(daNode, nullptr, BAD_CAST "filebase", BAD_CAST "test.uda");
  xmlNewChild(daNode, nullptr, BAD_CAST "outputInterval", BAD_CAST "0.1");
  
  xmlNodePtr checkpointNode = xmlNewChild(daNode, nullptr, BAD_CAST "checkpoint", nullptr);
  xmlNewProp(checkpointNode, BAD_CAST "interval", BAD_CAST "1.0");
  
  xmlNodePtr saveNode = xmlNewChild(daNode, nullptr, BAD_CAST "save", nullptr);
  xmlNewProp(saveNode, BAD_CAST "label", BAD_CAST "p.x");

  ProblemSpecP ps = scinew ProblemSpec(root, false);
  
  try {
    data_archiver->problemSetup(ps, nullptr, mat_manager);
  } catch (Exception& e) {
    FAIL() << "DataArchiver::problemSetup failed: " << e.message();
  } catch (...) {
    FAIL() << "DataArchiver::problemSetup failed with unknown exception";
  }
  
  EXPECT_DOUBLE_EQ(data_archiver->getOutputInterval(), 0.1);
  EXPECT_DOUBLE_EQ(data_archiver->getCheckpointInterval(), 1.0);
  
  xmlFreeDoc(doc);
}
