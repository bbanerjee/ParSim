#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Exceptions/InternalError.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(ParallelTest, Initialization) {
  EXPECT_TRUE(Parallel::usingMPI());
}

TEST(ParallelTest, Accessors) {
  Parallel::setNumThreads(4);
  EXPECT_EQ(Parallel::getNumThreads(), 4);
  
  Parallel::setUsingCPU(true);
  EXPECT_TRUE(Parallel::usingCPU());
}

TEST(ProcessorGroupTest, Uninitialized) {
  // Parallel is not initialized in this test environment, 
  // so it should throw InternalError.
  EXPECT_THROW(Parallel::getRootProcessorGroup(), InternalError);
}
