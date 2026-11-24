#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/BBox.h>
#include <gtest/gtest.h>

using namespace Uintah;

TEST(PointTest, Construction) {
  Point p1(1.0, 2.0, 3.0);
  EXPECT_EQ(p1.x(), 1.0);
  EXPECT_EQ(p1.y(), 2.0);
  EXPECT_EQ(p1.z(), 3.0);
  
  Point p2 = p1;
  EXPECT_EQ(p2, p1);
}

TEST(VectorTest, Construction) {
  Vector v1(1.0, 0.0, 0.0);
  EXPECT_EQ(v1.x(), 1.0);
  EXPECT_EQ(v1.length(), 1.0);
  
  Vector v2(0.0, 1.0, 0.0);
  EXPECT_EQ(Dot(v1, v2), 0.0);
}

TEST(PointVectorTest, Arithmetic) {
  Point p(1.0, 1.0, 1.0);
  Vector v(1.0, 2.0, 3.0);
  
  Point p2 = p + v;
  EXPECT_EQ(p2.x(), 2.0);
  EXPECT_EQ(p2.y(), 3.0);
  EXPECT_EQ(p2.z(), 4.0);
  
  Vector v2 = p2 - p;
  EXPECT_EQ(v2, v);
}

TEST(BBoxTest, ConstructionAndExtend) {
  BBox box;
  EXPECT_FALSE(box.valid());
  
  Point p1(0.0, 0.0, 0.0);
  Point p2(1.0, 1.0, 1.0);
  
  box.extend(p1);
  EXPECT_TRUE(box.valid());
  EXPECT_EQ(box.min(), p1);
  EXPECT_EQ(box.max(), p1);
  
  box.extend(p2);
  EXPECT_EQ(box.min(), p1);
  EXPECT_EQ(box.max(), p2);
  
  EXPECT_TRUE(box.inside(Point(0.5, 0.5, 0.5)));
  EXPECT_FALSE(box.inside(Point(1.5, 1.5, 1.5)));
}
