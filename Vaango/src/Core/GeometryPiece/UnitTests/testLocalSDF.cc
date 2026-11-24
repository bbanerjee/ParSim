#include <Core/GeometryPiece/LocalSDF.h>
#include <Core/Grid/Box.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Math/Matrix3.h>
#include <gtest/gtest.h>
#include <cmath>

using namespace Uintah;

TEST(LocalSDFTest, BasicTransformation) {
  Box box(Point(0, 0, 0), Point(1, 1, 1));
  IntVector res(11, 11, 11);
  LocalSDF sdf(box, res);

  Point center(5, 5, 5);
  Matrix3 orientation;
  orientation.Identity();

  Point worldP(6, 6, 6);
  Point localP = sdf.worldToLocal(worldP, center, orientation);
  EXPECT_NEAR(localP.x(), 1.0, 1e-12);
  EXPECT_NEAR(localP.y(), 1.0, 1e-12);
  EXPECT_NEAR(localP.z(), 1.0, 1e-12);

  Vector localV(1, 0, 0);
  Vector worldV = sdf.localToWorld(localV, orientation);
  EXPECT_NEAR(worldV.x(), 1.0, 1e-12);
  EXPECT_NEAR(worldV.y(), 0.0, 1e-12);
  EXPECT_NEAR(worldV.z(), 0.0, 1e-12);
}

TEST(LocalSDFTest, Interpolation) {
  Box box(Point(0, 0, 0), Point(1, 1, 1));
  IntVector res(2, 2, 2);
  LocalSDF sdf(box, res);

  // Set values at corners
  sdf.setDistance(IntVector(0, 0, 0), 0.0);
  sdf.setDistance(IntVector(1, 0, 0), 1.0);
  sdf.setDistance(IntVector(0, 1, 0), 0.0);
  sdf.setDistance(IntVector(1, 1, 0), 1.0);
  sdf.setDistance(IntVector(0, 0, 1), 0.0);
  sdf.setDistance(IntVector(1, 0, 1), 1.0);
  sdf.setDistance(IntVector(0, 1, 1), 0.0);
  sdf.setDistance(IntVector(1, 1, 1), 1.0);

  // Center should be 0.5
  double dist = sdf.getDistance(Point(0.5, 0.5, 0.5));
  EXPECT_NEAR(dist, 0.5, 1e-12);

  // Gradient should be (1, 0, 0)
  Vector grad = sdf.getGradient(Point(0.5, 0.5, 0.5));
  EXPECT_NEAR(grad.x(), 1.0, 1e-7);
  EXPECT_NEAR(grad.y(), 0.0, 1e-7);
  EXPECT_NEAR(grad.z(), 0.0, 1e-7);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
