#include <Core/Containers/RangeTree.h>
#include <gtest/gtest.h>
#include <list>
#include <vector>

using namespace Uintah;

// Minimal point class for testing RangeTree
struct TestPoint {
  int coords[3];
  int operator[](int i) const { return coords[i]; }
};

TEST(RangeTreeTest, ConstructionAndBasicQuery) {
  TestPoint p1 = {{1, 2, 3}};
  TestPoint p2 = {{4, 5, 6}};
  TestPoint p3 = {{7, 8, 9}};
  
  std::list<TestPoint*> points = {&p1, &p2, &p3};
  RangeTree<TestPoint, int> tree(points, 3);
  
  std::list<TestPoint*> found;
  TestPoint low = {{0, 0, 0}};
  TestPoint high = {{5, 5, 5}};
  
  tree.query(low, high, found);
  
  EXPECT_EQ(found.size(), 1u);
  EXPECT_EQ((*found.front())[0], 1);
}

TEST(RangeTreeTest, QueryAll) {
  TestPoint p1 = {{1, 2, 3}};
  TestPoint p2 = {{4, 5, 6}};
  
  std::list<TestPoint*> points = {&p1, &p2};
  RangeTree<TestPoint, int> tree(points, 3);
  
  std::list<TestPoint*> found;
  TestPoint low = {{0, 0, 0}};
  TestPoint high = {{10, 10, 10}};
  
  tree.query(low, high, found);
  EXPECT_EQ(found.size(), 2u);
}

TEST(RangeTreeTest, SphereQuery) {
  TestPoint p1 = {{1, 1, 1}};
  TestPoint p2 = {{10, 10, 10}};
  
  std::list<TestPoint*> points = {&p1, &p2};
  RangeTree<TestPoint, int> tree(points, 3);
  
  std::list<TestPoint*> found;
  TestPoint center = {{0, 0, 0}};
  int radius = 2; // Should include p1 but not p2
  
  tree.querySphere(center, radius, found);
  EXPECT_EQ(found.size(), 1u);
  EXPECT_EQ((*found.front())[0], 1);
}
