#include <Core/Math/Matrix3.h>
#include <Core/Math/MiscMath.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(Matrix3Test, Construction) {
  Matrix3 m(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(m(2, 2), 9.0);
}

TEST(Matrix3Test, Identity) {
  Matrix3 m;
  m.Identity();
  EXPECT_DOUBLE_EQ(m(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(m(1, 1), 1.0);
  EXPECT_DOUBLE_EQ(m(2, 2), 1.0);
  EXPECT_DOUBLE_EQ(m(0, 1), 0.0);
}

TEST(Matrix3Test, Determinant) {
  Matrix3 m(1, 2, 3, 0, 1, 4, 5, 6, 0);
  // Det = 1*(0-24) - 2*(0-20) + 3*(0-5) = -24 + 40 - 15 = 1
  EXPECT_DOUBLE_EQ(m.Determinant(), 1.0);
}

TEST(Matrix3Test, Inverse) {
  Matrix3 m(1, 2, 3, 0, 1, 4, 5, 6, 0);
  Matrix3 inv = m.Inverse();
  Matrix3 res = m * inv;
  Matrix3 id; id.Identity();
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      EXPECT_NEAR(res(i, j), id(i, j), 1e-12);
}

TEST(MiscMathTest, MinMax) {
  EXPECT_EQ(Min(10, 20), 10);
  EXPECT_EQ(Max(10, 20), 20);
  EXPECT_DOUBLE_EQ(Clamp(15.0, 0.0, 10.0), 10.0);
}
