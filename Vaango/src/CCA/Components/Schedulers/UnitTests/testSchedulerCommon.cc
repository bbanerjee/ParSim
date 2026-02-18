#include <CCA/Components/Schedulers/SchedulerCommon.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <memory>

using namespace Uintah;

class MockScheduler : public SchedulerCommon {
public:
  MockScheduler(const ProcessorGroup* myworld) : SchedulerCommon(myworld) {}
  virtual ~MockScheduler() {}

  virtual void execute([[maybe_unused]] int tgnum = 0, [[maybe_unused]] int iteration = 0) override {}
  virtual SchedulerP createSubScheduler() override { return nullptr; }
  virtual void verifyChecksum() override {}
};

class SchedulerCommonTest : public ::testing::Test {
protected:
  virtual void SetUp() override {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testSchedulerCommon", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
    world = Uintah::Parallel::getRootProcessorGroup();
    scheduler = std::make_unique<MockScheduler>(world);
  }

  const ProcessorGroup* world;
  std::unique_ptr<MockScheduler> scheduler;
};

TEST_F(SchedulerCommonTest, ConstructorTest) {
  EXPECT_EQ(scheduler->getNumTaskGraphs(), 0);
  EXPECT_FALSE(scheduler->isRestartInitTimestep());
}

TEST_F(SchedulerCommonTest, ProblemSetupTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr schNode = xmlNewChild(root, nullptr, BAD_CAST "Scheduler", nullptr);
  xmlNewChild(schNode, nullptr, BAD_CAST "small_messages", BAD_CAST "false");

  ProblemSpecP ps = scinew ProblemSpec(root, false);
  MaterialManagerP mat_manager = std::make_shared<MaterialManager>();
  
  scheduler->problemSetup(ps, mat_manager);
  
  EXPECT_FALSE(scheduler->useSmallMessages());
  
  xmlFreeDoc(doc);
}

TEST_F(SchedulerCommonTest, OverrideVariableBehaviorTest) {
  scheduler->overrideVariableBehavior("test_var", true, true, true, true, true);
  
  EXPECT_TRUE(scheduler->getNoScrubVars().count("test_var"));
  EXPECT_TRUE(scheduler->getCopyDataVars().count("test_var"));
  EXPECT_TRUE(scheduler->getNotCheckPointVars().count("test_var"));
}
