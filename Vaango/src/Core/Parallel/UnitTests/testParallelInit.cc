#include <Core/Parallel/Parallel.h>
#include <gtest/gtest.h>
#include <cstring>

using namespace Uintah;

TEST(ParallelInitTest, ParallelInitialize) {
  int argc = 1;
  char* argv_0 = strdup("testParallelInit");
  char* argv[] = {argv_0, nullptr};
  char** argv_ptr = argv;
  
  if (!Parallel::isInitialized()) {
    Parallel::initializeManager(argc, argv_ptr);
  }
  
  EXPECT_TRUE(Parallel::isInitialized());
  EXPECT_TRUE(Parallel::usingMPI());
  
  free(argv_0);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
