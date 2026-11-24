#include <CCA/Components/Regridder/RegridderCommon.h>
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

class MockRegridder : public RegridderCommon {
public:
  MockRegridder(const ProcessorGroup* pg) : RegridderCommon(pg) {}
  virtual ~MockRegridder() {}

  virtual std::string getName() override { return "MockRegridder"; }
  virtual Grid* regrid([[maybe_unused]] Grid* oldGrid, [[maybe_unused]] int time_step) override { return nullptr; }
  virtual std::vector<IntVector> getMinPatchSize() override { return std::vector<IntVector>(); }
};

class RegridderCommonTest : public ::testing::Test {
protected:
  virtual void SetUp() override {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testRegridderCommon", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
    world = Uintah::Parallel::getRootProcessorGroup();
    regridder = std::make_unique<MockRegridder>(world);
  }

  const ProcessorGroup* world;
  std::unique_ptr<MockRegridder> regridder;
};

TEST_F(RegridderCommonTest, ConstructorTest) {
  EXPECT_FALSE(regridder->isAdaptive());
  EXPECT_EQ(regridder->maxLevels(), 1);
}

TEST_F(RegridderCommonTest, AdaptivityTest) {
  regridder->setAdaptivity(true);
  EXPECT_TRUE(regridder->isAdaptive());
  regridder->setAdaptivity(false);
  EXPECT_FALSE(regridder->isAdaptive());
}
