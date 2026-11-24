#include <Core/Malloc/Allocator.h>
#include <gtest/gtest.h>
#include <vector>

#include <sci_defs/malloc_defs.h>

#ifndef DISABLE_SCI_MALLOC

using namespace Uintah;

TEST(AllocatorTest, DefaultAllocator) {
  if (!default_allocator) {
    MakeDefaultAllocator();
  }
  ASSERT_NE(default_allocator, nullptr);
  
  void* p = malloc(100);
  ASSERT_NE(p, nullptr);
  free(p);
}

TEST(AllocatorTest, SciNew) {
  if (!default_allocator) {
    MakeDefaultAllocator();
  }
  
  int* p = scinew int;
  ASSERT_NE(p, nullptr);
  *p = 42;
  EXPECT_EQ(*p, 42);
  delete p;
}

TEST(AllocatorTest, ArraySciNew) {
  if (!default_allocator) {
    MakeDefaultAllocator();
  }
  
  int* p = scinew int[10];
  ASSERT_NE(p, nullptr);
  for(int i=0; i<10; ++i) p[i] = i;
  EXPECT_EQ(p[5], 5);
  delete[] p;
}

#else

TEST(AllocatorTest, Disabled) {
  // Test that scinew still works (it should be defined to new)
  int* p = scinew int;
  ASSERT_NE(p, nullptr);
  *p = 42;
  EXPECT_EQ(*p, 42);
  delete p;
}

#endif
