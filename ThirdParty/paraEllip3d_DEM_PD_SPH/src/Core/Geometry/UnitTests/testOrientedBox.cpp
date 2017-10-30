#include <Core/Geometry/OrientedBox.h>
#include <Core/Geometry/Box.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(OrientedBoxTest, intersect) {

  Box box1(Vec(0, 0, 0), Vec(1, 1, 1));
  Box box2(Vec(1, 1, 1), Vec(2, 2, 2));

  OrientedBox obox1(box1);
  OrientedBox obox2(box2);

  bool intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, true);

  obox2.translate(Vec(0.001, 0.001, 0.001));
  intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, false);

  obox1.rotate(30*3.1412147/180, Vec(0, 0, 1));
  intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, false);

  obox2.translate(Vec(-0.1, -0.1, -0.1));
  obox2.rotate(30*3.1412147/180, Vec(0, 0, 1));
  intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, false);

  std::vector<Vec> vert1 = obox1.vertices();
  std::vector<Vec> vert2 = obox2.vertices();

  obox2.rotate(30*3.1412147/180, Vec(0, 0, 1));
  intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, false);

  obox1.rotate(-18*3.1412147/180, Vec(0, 0, 1));
  obox2.rotate(30*3.1412147/180, Vec(0, 0, 1));
  intersects = obox1.intersects(obox2);
  EXPECT_EQ(intersects, true);

  vert1 = obox1.vertices();
  vert2 = obox2.vertices();
  EXPECT_NEAR(vert1[0].x(), -0.0930199, 1.0e-7);
  EXPECT_NEAR(vert1[0].y(), 0.1148671, 1.0e-7);
  EXPECT_NEAR(vert2[0].x(), 0.9009055, 1.0e-7);
  EXPECT_NEAR(vert2[0].y(), 1.9009055, 1.0e-7);
  //std::copy(vert1.begin(), vert1.end(),
  //          std::ostream_iterator<Vec>(std::cout, "\n"));
  //std::copy(vert2.begin(), vert2.end(),
  //          std::ostream_iterator<Vec>(std::cout, "\n"));
}
