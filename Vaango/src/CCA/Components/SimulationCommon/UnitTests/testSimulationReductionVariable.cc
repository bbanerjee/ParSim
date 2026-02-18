#include <CCA/Components/SimulationCommon/SimulationReductionVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Malloc/Allocator.h>

#include <gtest/gtest.h>
#include <limits>

using namespace Uintah;

class SimulationReductionVariableTest : public ::testing::Test {
protected:
  virtual void SetUp() override {
    if (!Uintah::Parallel::isInitialized()) {
       static char* argv[] = {(char*)"testSimulationReductionVariable", nullptr};
       static int argc = 1;
       char** argv_ptr = argv;
       Uintah::Parallel::initializeManager(argc, argv_ptr);
    }
  }
};

TEST_F(SimulationReductionVariableTest, BasicTest) {
  SimulationReductionVariable var("test_min", min_vartype::getTypeDescription(), true);
  
  EXPECT_EQ(var.getName(), "test_min");
  EXPECT_TRUE(var.getActive());
  EXPECT_TRUE(var.isBenignValue());
  EXPECT_FALSE(var.overridden());
  
  var.setActive(false);
  EXPECT_FALSE(var.getActive());
}

TEST_F(SimulationReductionVariableTest, ValueTest) {
  SimulationReductionVariable var("test_max", max_vartype::getTypeDescription(), true);
  
  // Before reduction, getValue returns the benign value of d_max_var
  // For max reduction, the benign value is the smallest possible double.
  EXPECT_DOUBLE_EQ(var.getValue(), -std::numeric_limits<double>::max());
  
  // setValue sets internal d_double_value and d_overrideValue = true
  var.setValue(nullptr, 100.0);
  
  // getValue still returns the internal d_max_var because reduce() hasn't been called.
  EXPECT_DOUBLE_EQ(var.getValue(), -std::numeric_limits<double>::max());
}
