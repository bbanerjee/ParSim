#include <Core/GeometryPiece/BoxGeometryPiece.h>
#include <Core/GeometryPiece/SphereGeometryPiece.h>
#include <Core/Geometry/Point.h>
#include <Core/Grid/Box.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(BoxGeometryPieceTest, Construction) {
  Point p1(0.0, 0.0, 0.0);
  Point p2(1.0, 1.0, 1.0);
  BoxGeometryPiece box(p1, p2);
  
  Box bb = box.getBoundingBox();
  EXPECT_EQ(bb.lower(), p1);
  EXPECT_EQ(bb.upper(), p2);
  EXPECT_DOUBLE_EQ(box.volume(), 1.0);
  EXPECT_TRUE(box.inside(Point(0.5, 0.5, 0.5)));
  EXPECT_FALSE(box.inside(Point(1.5, 1.5, 1.5)));
}

TEST(SphereGeometryPieceTest, Construction) {
  Point center(0.0, 0.0, 0.0);
  double radius = 1.0;
  SphereGeometryPiece sphere(center, radius);
  
  EXPECT_EQ(sphere.getCenter(), center);
  EXPECT_DOUBLE_EQ(sphere.radius(), radius);
  EXPECT_DOUBLE_EQ(sphere.volume(), 4.0/3.0 * M_PI);
  EXPECT_TRUE(sphere.inside(Point(0.5, 0.5, 0.5)));
  EXPECT_FALSE(sphere.inside(Point(1.5, 1.5, 1.5)));
}
