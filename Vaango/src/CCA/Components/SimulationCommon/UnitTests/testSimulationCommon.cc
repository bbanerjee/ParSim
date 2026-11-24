#include <CCA/Components/SimulationCommon/SimulationCommon.h>
#include <CCA/Components/SimulationCommon/SimulationReductionVariable.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>
#include <memory>

using namespace Uintah;

class MockSimulationCommon : public SimulationCommon {
public:
  MockSimulationCommon(const ProcessorGroup* myworld,
                       const MaterialManagerP materialManager)
    : SimulationCommon(myworld, materialManager) {}

  // SimulationCommon pure virtuals
  virtual void problemSetup(const ProblemSpecP& params,
                           const ProblemSpecP& restart_prob_spec,
                           GridP& grid,
                           const std::string& input_ups_dir = "") override {}

  virtual void outputProblemSpec(ProblemSpecP& ps) override {}

  virtual void scheduleInitialize(const LevelP& level, SchedulerP& scheduler) override {}

  virtual void scheduleRestartInitialize(const LevelP& level, SchedulerP& scheduler) override {}

  virtual void scheduleComputeStableTimestep(const LevelP& level, SchedulerP& scheduler) override {}

  // SimulationInterface pure virtuals (not already covered by SimulationCommon)
  virtual void restartInitialize() override {}
  
  // These are already implemented in SimulationCommon, but we need to satisfy SimulationInterface
  // if SimulationCommon didn't override them as public.
  // Actually SimulationCommon DOES implement them.

  // Expose private/protected methods for testing
  using SimulationCommon::problemSetup;
  using SimulationCommon::problemSetupDeltaT;
  using SimulationCommon::setDelT;
  using SimulationCommon::validateNextDelT;
  using SimulationCommon::isLastTimestep;
};

class SimulationCommonTest : public ::testing::Test {
protected:
  virtual void SetUp() override {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testSimulationCommon", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
    world = Uintah::Parallel::getRootProcessorGroup();
    mat_manager = std::make_shared<MaterialManager>();
    sim_common = std::make_unique<MockSimulationCommon>(world, mat_manager);
  }

  const ProcessorGroup* world;
  MaterialManagerP mat_manager;
  std::unique_ptr<MockSimulationCommon> sim_common;
};

TEST_F(SimulationCommonTest, ReductionVariableTest) {
  sim_common->addReductionVariable("test_bool", bool_or_vartype::getTypeDescription(), true);
  // Default variables are: outputTimestep, checkpointTimestep, outputInterval, checkpointInterval, 
  // recomputeTimestep, abortTimestep, endSimulation. (Total 7)
  EXPECT_EQ(sim_common->numReductionVariable(), 8u); 
  EXPECT_TRUE(sim_common->activeReductionVariable("test_bool"));
  EXPECT_EQ(sim_common->getReductionVariable("test_bool"), 0.0);
}

TEST_F(SimulationCommonTest, TimeStepTest) {
  sim_common->setTimestepsMax(10);
  EXPECT_EQ(sim_common->getTimestepsMax(), 10);
  
  // We can't easily test isLastTimestep because it might call Bcast if d_wallTimeMax > 0
}

TEST_F(SimulationCommonTest, ProblemSetupTest) {
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
  xmlNodePtr root = xmlNewNode(nullptr, BAD_CAST "Uintah_specification");
  xmlDocSetRootElement(doc, root);
  
  xmlNodePtr timeNode = xmlNewChild(root, nullptr, BAD_CAST "Time", nullptr);
  xmlNewChild(timeNode, nullptr, BAD_CAST "initTime", BAD_CAST "0.0");
  xmlNewChild(timeNode, nullptr, BAD_CAST "maxTime", BAD_CAST "1.0");
  xmlNewChild(timeNode, nullptr, BAD_CAST "timestep_multiplier", BAD_CAST "0.5");
  xmlNewChild(timeNode, nullptr, BAD_CAST "delt_min", BAD_CAST "1e-6");
  xmlNewChild(timeNode, nullptr, BAD_CAST "delt_max", BAD_CAST "0.1");
  xmlNewChild(timeNode, nullptr, BAD_CAST "max_timesteps", BAD_CAST "10");

  ProblemSpecP ps = scinew ProblemSpec(root, false);
  
  sim_common->problemSetup(ps);
  
  EXPECT_EQ(sim_common->getSimTimeMax(), 1.0);
  EXPECT_EQ(sim_common->getDelTMultiplier(), 0.5);
  EXPECT_EQ(sim_common->getDelTMin(), 1e-6);
  EXPECT_EQ(sim_common->getDelTMax(), 0.1);
  EXPECT_EQ(sim_common->getTimestepsMax(), 10);
  
  xmlFreeDoc(doc);
}

TEST_F(SimulationCommonTest, ValidateNextDelTTest) {
  sim_common->setDelTMax(0.1);
  sim_common->setDelTMin(1e-6);
  sim_common->setDelTMaxIncrease(0.1); 
  sim_common->setDelT(0.05);
  
  double nextDelT = 0.06; 
  sim_common->validateNextDelT(nextDelT, 0);
  EXPECT_NEAR(nextDelT, 0.055, 1e-9);
  
  nextDelT = 1e-7; 
  sim_common->validateNextDelT(nextDelT, 0);
  EXPECT_NEAR(nextDelT, 1e-6, 1e-9);
  
  nextDelT = 0.2; 
  sim_common->setDelTMaxIncrease(0.0); // Disable max increase check to test maxDelT clamping
  sim_common->validateNextDelT(nextDelT, 0);
  EXPECT_NEAR(nextDelT, 0.1, 1e-9);
}
