#include <Core/Util/Utility.h>
#include <gtest/gtest.h>

using namespace util;

TEST(LinspaceTest, intLinspace1) {

  // Test case where there are two items
  int low = 0;
  int high = 1;
  int spacing = 1;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, spacing);
  EXPECT_EQ(vec[0], 0);
  EXPECT_EQ(vec[1], 1);

  int numSpaces = 1;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], 0);
  EXPECT_EQ(vec[1], 1);
}

TEST(LinspaceTest, intLinspace2) {
  int low = 1;
  int high = 0;
  int spacing = 1;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, spacing);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 0);

  int numSpaces = -1;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 0);

  numSpaces = 0;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 0);

  numSpaces = 5;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], 1);
  EXPECT_EQ(vec[1], 0);
}

TEST(LinspaceTest, intLinspace3) {
  int low = 0;
  int high = 1;
  int spacing = 2;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, spacing);
  EXPECT_EQ(vec[0], 0);
  EXPECT_EQ(vec[1], 1);
}

TEST(LinspaceTest, intLinspace4) {
  int low = -1;
  int high = 5;
  double rspacing = 0.5;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, rspacing);
  EXPECT_EQ(vec[0], -1);
  EXPECT_EQ(vec[1], 5);
}

TEST(LinspaceTest, intLinspace5) {
  int low = -1;
  int high = 5;
  int rspacing = 1;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, rspacing);
  EXPECT_EQ(vec[0], -1);
  EXPECT_EQ(vec[1], 0);
  EXPECT_EQ(vec[2], 1);
  EXPECT_EQ(vec[3], 2);
  EXPECT_EQ(vec[4], 3);
  EXPECT_EQ(vec[5], 4);
  EXPECT_EQ(vec[6], 5);

  int numSpaces = 6;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], -1);
  EXPECT_EQ(vec[1], 0);
  EXPECT_EQ(vec[2], 1);
  EXPECT_EQ(vec[3], 2);
  EXPECT_EQ(vec[4], 3);
  EXPECT_EQ(vec[5], 4);
  EXPECT_EQ(vec[6], 5);
}

TEST(LinspaceTest, intLinspace6) {
  int low = -1;
  int high = 5;
  double rspacing = 2.6;
  std::vector<int> vec = util::linspaceApprox<int>(low, high, rspacing);
  EXPECT_EQ(vec[0], -1);
  EXPECT_EQ(vec[1], 1);
  EXPECT_EQ(vec[2], 3);
  EXPECT_EQ(vec[3], 5);

  int numSpaces = 3;
  vec = util::linspace<int>(low, high, numSpaces);
  EXPECT_EQ(vec[0], -1);
  EXPECT_EQ(vec[1], 1);
  EXPECT_EQ(vec[2], 3);
  EXPECT_EQ(vec[3], 5);
}

TEST(LinspaceTest, REALLinspace1) {

  // Test case where there are two items
  REAL low = 0;
  REAL high = 1;
  REAL spacing = 1;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, spacing);
  EXPECT_DOUBLE_EQ(vec[0], 0);
  EXPECT_DOUBLE_EQ(vec[1], 1);

  int numSpaces = 1;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], 0);
  EXPECT_DOUBLE_EQ(vec[1], 1);
}

TEST(LinspaceTest, REALLinspace2) {
  REAL low = 1;
  REAL high = 0;
  REAL spacing = 1;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, spacing);
  EXPECT_DOUBLE_EQ(vec[0], 1);
  EXPECT_DOUBLE_EQ(vec[1], 0);

  int numSpaces = -1;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], 1);
  EXPECT_DOUBLE_EQ(vec[1], 0);

  numSpaces = 0;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], 1);
  EXPECT_DOUBLE_EQ(vec[1], 0);

  numSpaces = 5;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], 1);
  EXPECT_DOUBLE_EQ(vec[1], 0.8);
  EXPECT_DOUBLE_EQ(vec[2], 0.6);
  EXPECT_DOUBLE_EQ(vec[3], 0.4);
  EXPECT_DOUBLE_EQ(vec[4], 0.2);
  EXPECT_DOUBLE_EQ(vec[5], 0);
}

TEST(LinspaceTest, REALLinspace3) {
  REAL low = 0;
  REAL high = 1;
  REAL spacing = 2;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, spacing);
  EXPECT_DOUBLE_EQ(vec[0], 0);
  EXPECT_DOUBLE_EQ(vec[1], 1);
}

TEST(LinspaceTest, REALLinspace4) {
  REAL low = -1;
  REAL high = 5;
  double rspacing = 0.5;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, rspacing);
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], -0.5);
  EXPECT_DOUBLE_EQ(vec[3], 0.5);
  EXPECT_DOUBLE_EQ(vec[5], 1.5);
  EXPECT_DOUBLE_EQ(vec[12], 5);

  int numSpaces = 12;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], -0.5);
  EXPECT_DOUBLE_EQ(vec[3], 0.5);
  EXPECT_DOUBLE_EQ(vec[5], 1.5);
  EXPECT_DOUBLE_EQ(vec[12], 5);
}

TEST(LinspaceTest, REALLinspace5) {
  REAL low = -1;
  REAL high = 5;
  REAL rspacing = 1;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, rspacing);
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], 0);
  EXPECT_DOUBLE_EQ(vec[2], 1);
  EXPECT_DOUBLE_EQ(vec[3], 2);
  EXPECT_DOUBLE_EQ(vec[4], 3);
  EXPECT_DOUBLE_EQ(vec[5], 4);
  EXPECT_DOUBLE_EQ(vec[6], 5);

  int numSpaces = 6;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], 0);
  EXPECT_DOUBLE_EQ(vec[2], 1);
  EXPECT_DOUBLE_EQ(vec[3], 2);
  EXPECT_DOUBLE_EQ(vec[4], 3);
  EXPECT_DOUBLE_EQ(vec[5], 4);
  EXPECT_DOUBLE_EQ(vec[6], 5);
}

TEST(LinspaceTest, REALLinspace6) {
  REAL low = -1;
  REAL high = 5;
  double rspacing = 2.6;
  std::vector<REAL> vec = util::linspaceApprox<REAL>(low, high, rspacing);
  //for (auto& val : vec) std::cout << val << " "; std::cout << "\n";
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], 1.6);
  EXPECT_DOUBLE_EQ(vec[2], 4.2);
  EXPECT_DOUBLE_EQ(vec[3], 5);

  int numSpaces = 3;
  vec = util::linspace<REAL>(low, high, numSpaces);
  EXPECT_DOUBLE_EQ(vec[0], -1);
  EXPECT_DOUBLE_EQ(vec[1], 1);
  EXPECT_DOUBLE_EQ(vec[2], 3);
  EXPECT_DOUBLE_EQ(vec[3], 5);
}
