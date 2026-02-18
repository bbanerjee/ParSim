#include <Core/Containers/ConsecutiveRangeSet.h>
#include <gtest/gtest.h>
#include <vector>
#include <string>

using namespace Uintah;

TEST(ConsecutiveRangeSetTest, DefaultConstructor) {
  ConsecutiveRangeSet set;
  EXPECT_EQ(set.size(), 0u);
  EXPECT_TRUE(set.begin() == set.end());
}

TEST(ConsecutiveRangeSetTest, RangeConstructor) {
  ConsecutiveRangeSet set(1, 5);
  EXPECT_EQ(set.size(), 5u);
  auto it = set.begin();
  EXPECT_EQ(*it, 1);
  EXPECT_EQ(*(++it), 2);
  EXPECT_EQ(*(++it), 3);
  EXPECT_EQ(*(++it), 4);
  EXPECT_EQ(*(it++), 4);
  EXPECT_EQ(*it, 5);
  EXPECT_TRUE(++it == set.end());
}

TEST(ConsecutiveRangeSetTest, StringConstructor) {
  ConsecutiveRangeSet set("1, 3-5, 10");
  EXPECT_EQ(set.size(), 5u);
  
  std::vector<int> expected = {1, 3, 4, 5, 10};
  std::vector<int> actual;
  for (int val : set) {
    actual.push_back(val);
  }
  EXPECT_EQ(actual, expected);
}

TEST(ConsecutiveRangeSetTest, AddInOrder) {
  ConsecutiveRangeSet set;
  set.addInOrder(1);
  set.addInOrder(2);
  set.addInOrder(5);
  set.addInOrder(6);
  set.addInOrder(7);
  
  EXPECT_EQ(set.size(), 5u);
  EXPECT_EQ(set.toString(), "1 - 2, 5 - 7");
}

TEST(ConsecutiveRangeSetTest, UnionIntersection) {
  ConsecutiveRangeSet set1("1-5");
  ConsecutiveRangeSet set2("3-7");
  
  ConsecutiveRangeSet unionSet = set1.unioned(set2);
  EXPECT_EQ(unionSet.toString(), "1 - 7");
  
  ConsecutiveRangeSet intersectSet = set1.intersected(set2);
  EXPECT_EQ(intersectSet.toString(), "3 - 5");
}
