#include <CCA/Components/LoadBalancers/LoadBalancerCommon.h>
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

class MockLoadBalancer : public LoadBalancerCommon {
public:
  MockLoadBalancer(const ProcessorGroup* myworld) : LoadBalancerCommon(myworld) {}
  virtual ~MockLoadBalancer() {}

  virtual bool needRecompile(const GridP&) override { return false; }
};

class LoadBalancerCommonTest : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testLoadBalancerCommon", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
    world = Uintah::Parallel::getRootProcessorGroup();
    load_balancer = std::make_unique<MockLoadBalancer>(world);
  }

  static void TearDownTestSuite() {
    load_balancer.reset();
  }

  static const ProcessorGroup* world;
  static std::unique_ptr<MockLoadBalancer> load_balancer;
};

const ProcessorGroup* LoadBalancerCommonTest::world = nullptr;
std::unique_ptr<MockLoadBalancer> LoadBalancerCommonTest::load_balancer = nullptr;

TEST_F(LoadBalancerCommonTest, ConstructorTest) {
  // Since we don't call problemSetup, we might need to set it to 1 if it's uninitialized
  load_balancer->setNthRank(1);
  EXPECT_EQ(load_balancer->getNthRank(), 1);
}

TEST_F(LoadBalancerCommonTest, NthRankTest) {
  load_balancer->setNthRank(2);
  EXPECT_EQ(load_balancer->getNthRank(), 2);
}
