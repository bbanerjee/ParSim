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

TEST(YieldCondUtilsTest, findClosestPoint)
{
  std::vector<Point> poly = {{
    Point(-692.82, -866.025, 0), Point(-17.3205, -173.205, 0), 
    Point(17.3205, 0, 0),        Point(-17.3205, 173.205, 0), 
    Point(-692.82, 866.025, 0),  Point(-1385.64, 1039.23, 0), 
    Point(-2771.28, 1212.44, 0), Point(-5542.56, 1385.64, 0), 
    Point(-11085.1, 1558.85, 0), Point(-11639.4, 1576.17, 0), 
    Point(-11694.8, 1577.9, 0) 
  }};
  {
    Point min_pt(0, 0, 0);
    Point pt(3464.1, 6928.2, 0);
    Vaango::Util::findClosestPoint(pt, poly, min_pt);
    EXPECT_NEAR(min_pt.x(), -6.9282e+02, 1.0e-10);
    EXPECT_NEAR(min_pt.y(), 8.66025e+02, 1.0e-10);
  }
  {
    Point min_pt(0, 0, 0);
    Point pt(-3464.1, 6928.2, 0);
    Vaango::Util::findClosestPoint(pt, poly, min_pt);
    EXPECT_NEAR(min_pt.x(), -3.8172391454861686e+03, 1.0e-10);
    EXPECT_NEAR(min_pt.y(), 1.2778105594520239e+03, 1.0e-10);
  }
  {
    Point min_pt(0, 0, 0);
    Point pt(5196.15, 0, 0);
    Vaango::Util::findClosestPoint(pt, poly, min_pt);
    EXPECT_NEAR(min_pt.x(), 1.73205e+01, 1.0e-10);
    EXPECT_NEAR(min_pt.y(), 0, 1.0e-10);
  }
  {
    Point min_pt(0, 0, 0);
    Point pt(5196.15, 1732.05, 0);
    Vaango::Util::findClosestPoint(pt, poly, min_pt);
    EXPECT_NEAR(min_pt.x(), -1.73205e+01, 1.0e-10);
    EXPECT_NEAR(min_pt.y(), 1.73205e+02, 1.0e-10);
  }
  {
    Point min_pt(0, 0, 0);
    Point pt(-5196.15, 1732.05, 0);
    Vaango::Util::findClosestPoint(pt, poly, min_pt);
    EXPECT_NEAR(min_pt.x(), -5.2190635849151222e+03, 1.0e-10);
    EXPECT_NEAR(min_pt.y(), 1.3654220577160372e+03, 1.0e-10);
  }
}

TEST(YieldCondUtilsTest, findClosestSegments)
{
  std::vector<Point> poly = {{
    Point(17.3205, 0, 0),        Point(-17.3205, 173.205, 0), 
    Point(-692.82, 866.025, 0),  Point(-1385.64, 1039.23, 0), 
    Point(-2771.28, 1212.44, 0), Point(-5542.56, 1385.64, 0), 
    Point(-11085.1, 1558.85, 0)
  }};
  {
    std::vector<Point> min_seg;
    Point pt(3464.1, 6928.2, 0);
    auto index = Vaango::Util::getClosestSegments(pt, poly, min_seg);
    EXPECT_EQ(index, 2);
  }
  {
    std::vector<Point> min_seg;
    Point pt(-3464.1, 6928.2, 0);
    auto index = Vaango::Util::getClosestSegments(pt, poly, min_seg);
    EXPECT_EQ(index, 4);
  }
  {
    std::vector<Point> min_seg;
    Point pt(5196.15, 0, 0);
    auto index = Vaango::Util::getClosestSegments(pt, poly, min_seg);
    EXPECT_EQ(index, 1);
  }
  {
    std::vector<Point> min_seg;
    Point pt(5196.15, 1732.05, 0);
    auto index = Vaango::Util::getClosestSegments(pt, poly, min_seg);
    EXPECT_EQ(index, 1);
  }
  {
    std::vector<Point> min_seg;
    Point pt(-5196.15, 1732.05, 0);
    auto index = Vaango::Util::getClosestSegments(pt, poly, min_seg);
    EXPECT_EQ(index, 5);
    //std::cout << "index = " << index;
    //std::copy(min_seg.begin(), min_seg.end(),
    //          std::ostream_iterator<Point>(std::cout, " "));
    //std::cout << std::endl;
  }
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

  auto array2 = Vaango::Util::linspace(0, 0.5, 1);
  //std::copy(array2.begin(), array2.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  array2 = Vaango::Util::linspace(0, 0.5, 5);
  //std::copy(array2.begin(), array2.end(),
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
  EXPECT_EQ(spline.size(), 24);
  EXPECT_NEAR(spline[0].x(), Point(1.000e+00, -1.000e+00, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[1].x(), Point(1.380e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[2].x(), Point(1.720e+00,  1.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[3].x(), Point(2.020e+00,  3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[4].x(), Point(2.280e+00,  2.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[5].x(), Point(2.500e+00,  0.000e+00, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[7].x(), Point(2.700e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[8].x(), Point(2.900e+00, -4.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[9].x(), Point(3.100e+00, -4.800e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[10].x(), Point(3.300e+00, -3.200e-01, 0).x(), 1.0e-10);
  EXPECT_NEAR(spline[11].y(), Point(3.500e+00,  0.000e+00, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[13].y(), Point(3.700e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[14].y(), Point(3.900e+00,  4.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[15].y(), Point(4.100e+00,  4.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[16].y(), Point(4.300e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[17].y(), Point(4.500e+00,  0.000e+00, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[19].y(), Point(4.720e+00, -2.800e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[20].y(), Point(4.980e+00, -3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[21].y(), Point(5.280e+00, -1.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[22].y(), Point(5.620e+00,  3.200e-01, 0).y(), 1.0e-10);
  EXPECT_NEAR(spline[23].y(), Point(6.000e+00,  1.000e+00, 0).y(), 1.0e-10);
}

TEST(YieldCondUtilsTest, intersectSegments)
{
  bool status;
  double t1, t2;
  Point intersect;

  // Parallel segments
  Point seg_1_start(0, 0, 0);
  Point seg_1_end(1, 1, 0);
  Point seg_2_start(1, 0, 0);
  Point seg_2_end(2, 1, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, false);
  ASSERT_EQ(t1, 0);
  ASSERT_EQ(intersect.x(), seg_2_start.x());
  ASSERT_EQ(intersect.y(), seg_2_start.y());

  // Collinear segments
  seg_2_start = Point(2, 2, 0);
  seg_2_end = Point(3, 3, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_1_start, seg_1_end, seg_2_start, seg_2_end);
  ASSERT_EQ(status, false);
  ASSERT_EQ(t1, 0);
  ASSERT_EQ(intersect.x(), seg_1_start.x());
  ASSERT_EQ(intersect.y(), seg_1_start.y());

  // Segments that don't intersect within start and end
  seg_2_start = Point(1, 2, 0);
  seg_2_end = Point(3, 0, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, false);
  ASSERT_EQ(t1, 0.25);
  ASSERT_EQ(t2, 1.5);
  ASSERT_EQ(intersect.x(), 1.5);
  ASSERT_EQ(intersect.y(), 1.5);

  // Segments that intersect within start and end
  seg_2_start = Point(0, 1.5, 0);
  seg_2_end = Point(2.5, 0, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(t1, 0.375);
  ASSERT_EQ(t2, 0.9375);
  ASSERT_EQ(intersect.x(), 0.9375);
  ASSERT_EQ(intersect.y(), 0.9375);

  // Segments that intersect at start pt
  seg_2_start = Point(0, 0, 0);
  seg_2_end = Point(2.5, 0, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(t2, 0);
  ASSERT_EQ(status, true);

  // Segments that intersect at end pt
  seg_2_start = Point(1, 1, 0);
  seg_2_end = Point(2.5, 0, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(t2, 1);

  // Segments that graze at end pt
  seg_1_start = Point(1, 0, 0);
  seg_1_end = Point(2, 0, 0);
  seg_2_start = Point(2, -1, 0);
  seg_2_end = Point(2, 1, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(t2, 1);
  //std::cout << status << " " << t1 << " " << t2 << " " << intersect << "\n";

  // Segments that graze at end pt
  seg_1_start = Point(0, 0, 0);
  seg_1_end = Point(1, 2, 0);
  seg_2_start = Point(1, 1, 0);
  seg_2_end = Point(1, 3, 0);
  std::tie(status, t1, t2, intersect) = 
    Vaango::Util::intersectionPoint(seg_2_start, seg_2_end, seg_1_start, seg_1_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(t2, 1);
  //std::cout << status << " " << t1 << " " << t2 << " " << intersect << "\n";
}

TEST(YieldCondUtilsTest, intersectPolylineSegment)
{
  std::vector<Point> polyline = 
    {{Point(0, 0 ,0), Point(1, 2, 0), Point(2, 3, 0), Point(3, 2, 0),
      Point(4, 3, 0), Point(5, 0, 0)}};

  bool status;
  std::size_t index;
  double t;
  Point intersect;

  // Outside segment
  Point seg_start(-1, 1, 0);
  Point seg_end(-1, 2, 0);
  std::tie(status, index, t, intersect) = 
    Vaango::Util::intersectionPoint(polyline, seg_start, seg_end);
  ASSERT_EQ(status, false);
  ASSERT_EQ(index, 5);

  // Parallel segment
  seg_start = Point(1, 1, 0);
  seg_end = Point(2, 2, 0);
  std::tie(status, index, t, intersect) = 
    Vaango::Util::intersectionPoint(polyline, seg_start, seg_end);
  ASSERT_EQ(status, false);
  ASSERT_EQ(index, 5);

  // Vertical end pt intersection segment
  seg_start = Point(1, 1, 0);
  seg_end = Point(1, 3, 0);
  std::tie(status, index, t, intersect) = 
    Vaango::Util::intersectionPoint(polyline, seg_start, seg_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(index, 0);
  ASSERT_EQ(t, 0.5);

  // Interior intersection
  seg_start = Point(2.5, 2, 0);
  seg_end = Point(2.5, 3, 0);
  std::tie(status, index, t, intersect) = 
    Vaango::Util::intersectionPoint(polyline, seg_start, seg_end);
  ASSERT_EQ(status, true);
  ASSERT_EQ(index, 2);
  ASSERT_EQ(t, 0.5);
  //std::cout << status << " " << index << " " << t << " " << intersect << "\n";
}

TEST(YieldCondUtilsTest, evalFunctinJacobianInverse)
{
  std::vector<Point> polyline = 
    {{Point(0, 0 ,0), Point(1, 2, 0), Point(2, 3, 0), Point(3, 2, 0),
      Point(4, 3, 0), Point(5, 0, 0), Point(6, 0, 0), Point(7, 0, 0)}};

  // Parallel segment
  Point bezier_p0 = polyline[0];
  Point bezier_p1 = polyline[1];
  Point bezier_p2 = polyline[2];
  Point seg_start(1, 1, 0);
  Point seg_end(2, 2, 0);
  Vector t(0.6, 0.5, 0);
  auto F_Jinv = 
    Vaango::Util::evalFunctionJacobianInverse(bezier_p0, bezier_p1, bezier_p2, 
                                               seg_start, seg_end, t);
  ASSERT_NEAR(F_Jinv.first.x(), -0.4, 1.0e-6);
  ASSERT_NEAR(F_Jinv.first.y(), 0.52, 1.0e-6);
  ASSERT_NEAR(F_Jinv.second(0,0), -2.5, 1.0e-6);
  ASSERT_NEAR(F_Jinv.second(0,1), 2.5, 1.0e-6);
  ASSERT_NEAR(F_Jinv.second(1,0), -3.5, 1.0e-6);

  // Perpendicular segment
  seg_start = Point(2, 1, 0);
  seg_end = Point(0, 3, 0);
  t = Vector(0.6, 0.5, 0);
  F_Jinv = 
    Vaango::Util::evalFunctionJacobianInverse(bezier_p0, bezier_p1, bezier_p2, 
                                               seg_start, seg_end, t);
  ASSERT_NEAR(F_Jinv.first.x(), 0.1, 1.0e-6);
  ASSERT_NEAR(F_Jinv.first.y(), 0.02, 1.0e-6);
  ASSERT_NEAR(F_Jinv.second(0,0), 0.41667, 1.0e-5);
  ASSERT_NEAR(F_Jinv.second(0,1), 0.41667, 1.0e-5);
  ASSERT_NEAR(F_Jinv.second(1,0), 0.291667, 1.0e-5);

  // Linear segment
  bezier_p0 = polyline[5];
  bezier_p1 = polyline[6];
  bezier_p2 = polyline[7];
  F_Jinv = 
    Vaango::Util::evalFunctionJacobianInverse(bezier_p0, bezier_p1, bezier_p2, 
                                               seg_start, seg_end, t);
  ASSERT_NEAR(F_Jinv.first.x(), 5.1, 1.0e-6);
  ASSERT_NEAR(F_Jinv.first.y(), -2, 1.0e-6);
  ASSERT_NEAR(F_Jinv.second(0,0), 1, 1.0e-5);
  ASSERT_NEAR(F_Jinv.second(0,1), 1, 1.0e-5);
  ASSERT_NEAR(F_Jinv.second(1,1), -0.5, 1.0e-5);
  //std::cout << "F = " << F_Jinv.first << " Jinv = " << F_Jinv.second << "\n";
}

TEST(YieldCondUtilsTest, intersectionPointBSpline)
{
  bool status;
  Vector t;
  Point intersection;
  std::vector<Point> polyline = 
    {{Point(0, 0 ,0), Point(1, 2, 0), Point(2, 3, 0), Point(3, 2, 0),
      Point(4, 3, 0), Point(5, 0, 0), Point(6, 0, 0), Point(7, 0, 0)}};

  // Parallel segment
  Point bezier_p0 = polyline[0];
  Point bezier_p1 = polyline[1];
  Point bezier_p2 = polyline[2];
  Point seg_start(1, 1, 0);
  Point seg_end(2, 2, 0);
  std::tie(status, t, intersection) = 
    Vaango::Util::intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, 
                                           seg_start, seg_end);
  ASSERT_EQ(status, false);
  ASSERT_NEAR(t.x(), -0.414214, 1.0e-6);
  ASSERT_NEAR(t.y(), -0.914214, 1.0e-6);
  ASSERT_NEAR(intersection.x(), 0.0857864, 1.0e-6);
  ASSERT_NEAR(intersection.y(), 0.0857864, 1.0e-6);

  // Perpendicular segment
  seg_start = Point(2, 1, 0);
  seg_end = Point(0, 3, 0);
  std::tie(status, t, intersection) = 
    Vaango::Util::intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, 
                                           seg_start, seg_end);
  ASSERT_EQ(status, true);
  ASSERT_NEAR(t.x(), 0.55051, 1.0e-5);
  ASSERT_NEAR(t.y(), 0.474745, 1.0e-6);
  ASSERT_NEAR(intersection.x(), 1.05051, 1.0e-6);
  ASSERT_NEAR(intersection.y(), 1.94949, 1.0e-6);

  // Linear segment
  bezier_p0 = polyline[5];
  bezier_p1 = polyline[6];
  bezier_p2 = polyline[7];
  std::tie(status, t, intersection) = 
    Vaango::Util::intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, 
                                           seg_start, seg_end);
  ASSERT_EQ(status, false);

  // Short perpendicular segment
  seg_start = Point(2, 0, 0);
  seg_end = Point(1, 1, 0);
  std::tie(status, t, intersection) = 
    Vaango::Util::intersectionPointBSpline(bezier_p0, bezier_p1, bezier_p2, 
                                           seg_start, seg_end);
  ASSERT_EQ(status, false);
  //std::cout << "Status = " << status << " t = " << t << " intersection = " << intersection << "\n";
}
