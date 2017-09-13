#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldCondUtils.h>

#include <Core/Geometry/Point.h>

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;
using Uintah::Point;
using Uintah::Vector;

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

TEST(YieldCondUtilsTest, computeNormals)
{
  std::vector<Point> points = {{Point(-10,0,0), Point(10,100,0),
                                Point(400,500,0), Point(800,600,0),
                                Point(1600,700,0), Point(3200,800,0),
                                Point(6400,900,0)}};
  std::vector<Point> polyline;
  polyline.push_back(Point(points[2].x(), -points[2].y(), 0));
  polyline.push_back(Point(points[1].x(), -points[1].y(), 0));
  polyline.insert(polyline.end(), points.begin(), points.end());
  Point last = points[points.size()-1];
  Point secondlast = points[points.size()-2];
  double t = 1.1;
  Vector extra1 = secondlast*(1 - t) + last*t;
  Vector extra2 = secondlast*(1 - t)*t + last*(t*t + 1 - t);
  polyline.push_back(Point(extra1));
  polyline.push_back(Point(extra2));
  //std::copy(polyline.begin(), polyline.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;

  std::vector<Vector> normals = Vaango::Util::computeNormals(polyline);

  //std::cout << "Normals:";
  //std::copy(normals.begin(), normals.end(),
  //          std::ostream_iterator<Uintah::Vector>(std::cout, " "));
  //std::cout << std::endl;
  EXPECT_EQ(points.size(), normals.size());
  EXPECT_NEAR(normals[0].x(), 
              Vector(-1.00000000000000e+00, -0.00000000000000e+00, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[0].y(), 
              Vector(-1.00000000000000e+00, -0.00000000000000e+00, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[1].x(), 
              Vector(-8.76013011703533e-01, 4.82287469592675e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[2].x(), 
              Vector(-5.00835396432814e-01, 8.65542549895720e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[3].x(), 
              Vector(-3.18081036180333e-01, 9.48063528684890e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[4].x(), 
              Vector(-1.03211498859446e-01, 9.94659432420558e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[5].x(), 
              Vector(-5.71593617013789e-02, 9.98365067182286e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[6].x(), 
              Vector(-3.64331413805806e-02, 9.99336092718132e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(normals[1].y(), 
              Vector(-8.76013011703533e-01, 4.82287469592675e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[2].y(), 
              Vector(-5.00835396432814e-01, 8.65542549895720e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[3].y(), 
              Vector(-3.18081036180333e-01, 9.48063528684890e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[4].y(), 
              Vector(-1.03211498859446e-01, 9.94659432420558e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[5].y(), 
              Vector(-5.71593617013789e-02, 9.98365067182286e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(normals[6].y(), 
              Vector(-3.64331413805806e-02, 9.99336092718132e-01, 0).y(), 1.0e-10);
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

TEST(YieldCondUtilsTest, computeBSpline)
{
  std::vector<Point> polyline = 
    {{Point(1, -1 ,0), Point(2, 1, 0), Point(3, -1, 0), Point(4, 1, 0),
      Point(5, -1, 0), Point(6, 1, 0)}};
  std::vector<Point> spline;
  Vaango::Util::computeOpenUniformQuadraticBSpline(polyline, 5, spline);
  //std::copy(spline.begin(), spline.end(),
  //          std::ostream_iterator<Point>(std::cout, " "));
  //std::cout << std::endl;
  EXPECT_EQ(spline.size(), 21);
  EXPECT_NEAR(spline[0].x(), Point(1.000e+00, -1.000e+00, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[1].x(), Point(1.380e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[2].x(), Point(1.720e+00,  1.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[3].x(), Point(2.020e+00,  3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[4].x(), Point(2.280e+00,  2.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[5].x(), Point(2.500e+00,  0.000e+00, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[6].x(), Point(2.700e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[7].x(), Point(2.900e+00, -4.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[8].x(), Point(3.100e+00, -4.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[9].x(), Point(3.300e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[10].y(), Point(3.500e+00,  0.000e+00, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[11].y(), Point(3.700e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[12].y(), Point(3.900e+00,  4.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[13].y(), Point(4.100e+00,  4.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[14].y(), Point(4.300e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[15].y(), Point(4.500e+00,  0.000e+00, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[16].y(), Point(4.720e+00, -2.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[17].y(), Point(4.980e+00, -3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[18].y(), Point(5.280e+00, -1.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[19].y(), Point(5.620e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[20].y(), Point(6.000e+00,  1.000e+00, 0).y(), 1.0e-10);
}