#include <Core/Geometry/Ellipsoid.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(EllipsoidTest, intersect) {

  Vec axis(0, 0, 1);
  REAL angle = -0.5*M_PI;
  REAL cx = 2; REAL cy = 1.5; REAL cz = 1.0;
  REAL lx = 1; REAL ly = 2; REAL lz = 3;
  Ellipsoid ell(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  //std::cout << "Ell: " << ell << "\n";

  REAL mx = 0; REAL my = 0; REAL mz = 0;
  REAL ex = 1; REAL ey = 3; REAL ez = 2;
  OrientedBox box(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  //std::cout << "Box: " << box << "\n";

  ell.rotate(angle, axis);
  //std::cout << "Ell: " << ell << "\n";
  Vec u0(0.00000000000000e+00,   1.00000000000000e+00,   0.00000000000000e+00);
  Vec u1(-1.00000000000000e+00,   1.11022302462516e-16,   0.00000000000000e+00);
  Vec u2(0.00000000000000e+00,   0.00000000000000e+00,   1.00000000000000e+00);
  EXPECT_NEAR(u0.x(), ell.axis(0).x(), 1.0e-10);
  EXPECT_NEAR(u1.x(), ell.axis(1).x(), 1.0e-10);
  EXPECT_NEAR(u2.x(), ell.axis(2).x(), 1.0e-10);

  box.rotate(0.5*angle, axis);
  //std::cout << "Box: " << box << "\n";
  Vec v0(7.07106781186548e-01,   7.07106781186547e-01,   0.00000000000000e+00);
  Vec v1(-7.07106781186547e-01,   7.07106781186548e-01,   0.00000000000000e+00);
  Vec v2(0.00000000000000e+00,   0.00000000000000e+00,   1.00000000000000e+00);
  EXPECT_NEAR(v0.x(), box.axis(0).x(), 1.0e-10);
  EXPECT_NEAR(v1.x(), box.axis(1).x(), 1.0e-10);
  EXPECT_NEAR(v2.x(), box.axis(2).x(), 1.0e-10);

  OrientedBox ell_box = ell.getOrientedBoundingBox();
  bool intersects = box.intersects(ell_box);
  EXPECT_EQ(intersects, true);

  //-----------------------------------------------
  angle = 0;
  cx = -2.0000001; cy = -1.5;
  Ellipsoid ell1(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  OrientedBox box1(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  ell1.rotate(angle, axis);
  box1.rotate(angle, axis);
  ell_box = ell1.getOrientedBoundingBox();
  intersects = box1.intersects(ell_box);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";
  EXPECT_EQ(intersects, false);

  //-----------------------------------------------
  angle = 0.5*M_PI;
  cx = -0.5; cy = -2.5;
  Ellipsoid ell2(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  OrientedBox box2(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  ell2.rotate(angle, axis);
  box2.rotate(angle, axis);
  ell_box = ell2.getOrientedBoundingBox();
  intersects = box2.intersects(ell_box);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";
  EXPECT_EQ(intersects, false);

  //-----------------------------------------------
  angle = 40*M_PI/180;
  cx = 2; cy = -2;
  Ellipsoid ell3(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  OrientedBox box3(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  ell3.rotate(angle, axis);
  ell_box = ell3.getOrientedBoundingBox();
  intersects = box3.intersects(ell_box);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";
  EXPECT_EQ(intersects, true);

  //std::vector<Vec> vert1 = box3.vertices();
  //std::copy(vert1.begin(), vert1.end(),
  //          std::ostream_iterator<Vec>(std::cout, "\n"));

  intersects = ell3.intersects(box3);
  EXPECT_EQ(intersects, true);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";

  //-----------------------------------------------
  angle = 40*M_PI/180;
  cx = 2.6; cy = -2;
  Ellipsoid ell4(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  OrientedBox box4(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  ell4.rotate(angle, axis);
  ell_box = ell4.getOrientedBoundingBox();
  intersects = box4.intersects(ell_box);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";
  EXPECT_EQ(intersects, true);

  intersects = ell4.intersects(box3);
  EXPECT_EQ(intersects, false);
  //std::cout << "Intersects = " << std::boolalpha << intersects << "\n";

  //-----------------------------------------------
  angle = 40*M_PI/180;
  cx = 2; cy = -2;
  Ellipsoid ell5(Vec(cx, cy, cz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                 lx, ly, lz);
  OrientedBox box5(Vec(mx, my, mz), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1),
                  ex, ey, ez);
  ell5.rotate(angle, axis);
  std::vector<Vec> vertices = box5.vertices();

  constexpr std::array<std::array<int, 4>, 6> faces = {{
    {{0, 4, 7, 3}}, // x-
    {{1, 2, 6, 5}}, // x+
    {{0, 1, 5, 4}}, // y-
    {{2, 3, 7, 6}}, // y+
    {{0, 3, 2, 1}}, // z-
    {{4, 5, 6, 7}}  // z+
  }};

  int faceID = 0;
  std::array<bool, 6> faceIntersects;
  for (const auto& face : faces) {
    int v0 = face[0];
    int v1 = face[1];
    int v2 = face[2];
    int v3 = face[3];
    Face ff(vertices[v0], vertices[v1], vertices[v2], vertices[v3]);
    faceIntersects[faceID] = (ell5.intersects(ff)).first;
    //std::cout << "Intersects face: " << faceID 
    //          << " = " << std::boolalpha << faceIntersects[faceID] << "\n";
    ++faceID;
  }

  EXPECT_EQ(faceIntersects[0], false);
  EXPECT_EQ(faceIntersects[1], true);
  EXPECT_EQ(faceIntersects[2], true);
  EXPECT_EQ(faceIntersects[3], false);
  EXPECT_EQ(faceIntersects[4], false);
  EXPECT_EQ(faceIntersects[5], true);
}
