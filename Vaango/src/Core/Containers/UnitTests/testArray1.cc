#include <Core/Containers/Array1.h>
#include <gtest/gtest.h>
#include <string>

using namespace Uintah;

TEST(Array1Test, DefaultConstructor) {
  Array1<int> arr;
  EXPECT_EQ(arr.size(), 0);
}

TEST(Array1Test, SizeConstructor) {
  Array1<int> arr(5);
  EXPECT_EQ(arr.size(), 5);
}

TEST(Array1Test, AddAndAccess) {
  Array1<int> arr;
  arr.add(10);
  arr.add(20);
  arr.add(30);
  
  EXPECT_EQ(arr.size(), 3);
  EXPECT_EQ(arr[0], 10);
  EXPECT_EQ(arr[1], 20);
  EXPECT_EQ(arr[2], 30);
}

TEST(Array1Test, Resize) {
  Array1<int> arr;
  arr.resize(3);
  EXPECT_EQ(arr.size(), 3);
  arr[0] = 1;
  arr[1] = 2;
  arr[2] = 3;
  
  arr.resize(5);
  EXPECT_EQ(arr.size(), 5);
  EXPECT_EQ(arr[0], 1);
  EXPECT_EQ(arr[1], 2);
  EXPECT_EQ(arr[2], 3);
  
  arr.resize(2);
  EXPECT_EQ(arr.size(), 2);
  EXPECT_EQ(arr[0], 1);
  EXPECT_EQ(arr[1], 2);
}

TEST(Array1Test, CopyConstructor) {
  Array1<int> arr1;
  arr1.add(1);
  arr1.add(2);
  
  Array1<int> arr2(arr1);
  EXPECT_EQ(arr2.size(), 2);
  EXPECT_EQ(arr2[0], 1);
  EXPECT_EQ(arr2[1], 2);
  
  arr1[0] = 10;
  EXPECT_EQ(arr2[0], 1); // Deep copy check
}

TEST(Array1Test, AssignmentOperator) {
  Array1<int> arr1;
  arr1.add(1);
  Array1<int> arr2;
  arr2 = arr1;
  
  EXPECT_EQ(arr2.size(), 1);
  EXPECT_EQ(arr2[0], 1);
}

TEST(Array1Test, InsertAndRemove) {
  Array1<int> arr;
  arr.add(1);
  arr.add(3);
  arr.insert(1, 2);
  
  EXPECT_EQ(arr.size(), 3);
  EXPECT_EQ(arr[0], 1);
  EXPECT_EQ(arr[1], 2);
  EXPECT_EQ(arr[2], 3);
  
  arr.remove(1);
  EXPECT_EQ(arr.size(), 2);
  EXPECT_EQ(arr[0], 1);
  EXPECT_EQ(arr[1], 3);
}

TEST(Array1Test, Initialize) {
  Array1<int> arr(5);
  arr.initialize(42);
  for(int i=0; i<5; ++i) {
    EXPECT_EQ(arr[i], 42);
  }
}
