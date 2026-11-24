#include <Core/Lockfree/Lockfree_Pool.hpp>
#include <gtest/gtest.h>
#include <string>

using namespace Lockfree;

TEST(LockfreePoolTest, Construction) {
  Pool<int> pool(10);
  EXPECT_EQ(pool.size(), 0u);
  EXPECT_TRUE(pool.empty());
  EXPECT_EQ(pool.num_levels(), 10u);
}

TEST(LockfreePoolTest, InsertAndFind) {
  Pool<int> pool(1);
  auto it1 = pool.insert(10);
  auto it2 = pool.insert(20);
  
  EXPECT_EQ(pool.size(), 2u);
  EXPECT_FALSE(pool.empty());
  
  auto found = pool.find_any([](int v){ return v == 20; });
  ASSERT_TRUE(static_cast<bool>(found));
  EXPECT_EQ(*found, 20);
}

TEST(LockfreePoolTest, Erase) {
  Pool<int> pool(1);
  auto it = pool.insert(42);
  EXPECT_EQ(pool.size(), 1u);
  
  pool.erase(it);
  EXPECT_EQ(pool.size(), 0u);
  EXPECT_TRUE(pool.empty());
}

TEST(LockfreePoolTest, CopyAndMove) {
  Pool<int> pool1(1);
  pool1.insert(100);
  
  // Shallow copy
  Pool<int> pool2 = pool1;
  EXPECT_EQ(pool2.size(), 1u);
  EXPECT_EQ(pool1.ref_count(), 2u);
  
  // Move
  Pool<int> pool3 = std::move(pool1);
  EXPECT_EQ(pool3.size(), 1u);
  EXPECT_EQ(pool3.ref_count(), 2u);
  EXPECT_EQ(pool1.num_levels(), 0u); // Invalidated
}
