#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>

#include <Core/Geometry/Point.h>

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::Point;

TEST(YieldCondUtilsTest, reverse)
{
  std::vector<double> array = {{1,2,3,4,5,6,7}};
  std::size_t index = 0;
  for (auto i : Vaango::Util::reverse(array)) {
    if (index == 0) EXPECT_EQ(i, 7);
    if (index == array.size()-1) EXPECT_EQ(i, 1);
    index++;
  }
}

TEST(YieldCondUtilsTest, convexHull2D)
{
  std::vector<Point> points = {{Point(-10,0,0), Point(10,100,0),
                                Point(200,200,0), Point(300,300,0),
                                Point(400,500,0), Point(800,600,0),
                                Point(1600,700,0), Point(3200,800,0),
                                Point(6400,900,0)}};
  //std::copy(points.begin(), points.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;

  std::vector<Point> hull = Vaango::Util::convexHull2D(points);
  //std::copy(hull.begin(), hull.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;

  EXPECT_EQ(hull.size(), 7);
  EXPECT_EQ(hull[0], Point(-10, 0, 0));
  EXPECT_EQ(hull[1], Point(6400, 900, 0));
  EXPECT_EQ(hull[2], Point(3200, 800, 0));
  EXPECT_EQ(hull[3], Point(1600, 700, 0));
  EXPECT_EQ(hull[4], Point(800, 600, 0));
  EXPECT_EQ(hull[5], Point(400, 500, 0));
  EXPECT_EQ(hull[6], Point(10, 100, 0));
}

TEST(YieldCondUtilsTest, linspace)
{
  std::vector<double> array;
  Vaango::Util::linspace(-1, 10, -1, array);
  EXPECT_EQ(array.size(), 0);
  Vaango::Util::linspace(-1, 10, 0, array);
  EXPECT_EQ(array.size(), 0);
  Vaango::Util::linspace(-1, 10, 1, array);
  EXPECT_EQ(array.size(), 2);
  Vaango::Util::linspace(-1, 10, 15, array);
  EXPECT_EQ(array.size(), 16);
  EXPECT_EQ(array.front(), -1);
  EXPECT_EQ(array.back(), 10);
  EXPECT_NEAR(array[4], 1.9333333, 1.0e-6);
  Vaango::Util::linspace(10, -1, 6, array);
  EXPECT_EQ(array.size(), 7);
  EXPECT_EQ(array.front(), 10);
  EXPECT_EQ(array.back(), -1);
  EXPECT_NEAR(array[4], 2.6666667, 1.0e-6);

  std::vector<double> array1 = Vaango::Util::linspace(-1, 10, -1);
  EXPECT_EQ(array1.size(), 0);
  array1 = Vaango::Util::linspace(-1, 10, 0);
  EXPECT_EQ(array1.size(), 0);
  array1 = Vaango::Util::linspace(-1, 10, 1);
  EXPECT_EQ(array1.size(), 2);
  array1 = Vaango::Util::linspace(-1, 10, 15);
  EXPECT_EQ(array1.size(), 16);
  EXPECT_EQ(array1.front(), -1);
  EXPECT_EQ(array1.back(), 10);
  EXPECT_NEAR(array1[4], 1.9333333, 1.0e-6);
  array1 = Vaango::Util::linspace(10, -1, 6);
  EXPECT_EQ(array1.size(), 7);
  EXPECT_EQ(array1.front(), 10);
  EXPECT_EQ(array1.back(), -1);
  EXPECT_NEAR(array1[4], 2.6666667, 1.0e-6);
  //std::copy(array1.begin(), array1.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;
}
