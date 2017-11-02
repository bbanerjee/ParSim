#include <Core/Geometry/FaceEdge.h>
#include <gtest/gtest.h>

using namespace dem;

TEST(FaceEdgeTest, edgeDistance) {

  Vertex v1(0, 0, 0);
  Vertex v2(1, 0, 0);
  Point pt(0.5, 0.5, 0);

  Edge edge(v1, v2);
  //std::cout << edge << "\n";
  auto td = edge.distance(pt);
  //std::cout << "(t,d) = (" << td.first << ": " << td.second << ")\n";
  EXPECT_EQ(td.first, 0.5);
  EXPECT_EQ(td.second, 0.5);

  pt = Point(-0.5, 0.5, 0);
  td = edge.distance(pt);
  EXPECT_EQ(td.first, -0.5);

  REAL x1 = -8.57106781186548e+00;
  REAL x2 = -7.15685424949238e+00;
  REAL x3 =  5.57106781186548e+00;
  REAL x4 =  4.15685424949238e+00;
  REAL y1 = -1.82842712474619e+00;
  REAL y2 = -2.53553390593274e+00;
  REAL y3 =  3.82842712474619e+00;
  REAL y4 =  4.53553390593274e+00;
  v1 = Vertex(x1, y1, 0);
  v2 = Vertex(x2, y2, 0);
  Vertex v3(x3, y3, 0);
  Vertex v4(x4, y4, 0);

  pt = Point(0, 0, 0);
  Edge e1(v1, v2);
  Edge e2(v2, v3);
  Edge e3(v3, v4);
  Edge e4(v4, v1);

  REAL t1 = 4.33137084989848e+00;
  REAL d1 = 5.46849299055311e+00;
  REAL t2 = 5.29520602772138e-01;
  REAL d2 = 9.32792052216177e-01;
  REAL t3 = 2.06862915010152e+00;
  REAL d3 = 5.91570658605306e+00;
  REAL t4 = 4.03812730561196e-01;
  REAL d4 = 2.19770311628353e+00;

  td = e1.distance(pt);
  EXPECT_NEAR(td.first, t1, 1.0e-10);
  EXPECT_NEAR(td.second, d1, 1.0e-10);

  td = e2.distance(pt);
  EXPECT_NEAR(td.first, t2, 1.0e-10);
  EXPECT_NEAR(td.second, d2, 1.0e-10);

  td = e3.distance(pt);
  EXPECT_NEAR(td.first, t3, 1.0e-10);
  EXPECT_NEAR(td.second, d3, 1.0e-10);

  td = e4.distance(pt);
  EXPECT_NEAR(td.first, t4, 1.0e-10);
  EXPECT_NEAR(td.second, d4, 1.0e-10);
}

TEST(FaceEdgeTest, faceDistance) {
  Vertex v1(-8.57106781186548e+00, -1.82842712474619e+00, -1.66666666666667e+00);
  Vertex v2(-7.15685424949238e+00, -2.53553390593274e+00, -1.66666666666667e+00);
  Vertex v3(5.57106781186548e+00, 3.82842712474619e+00, -1.66666666666667e+00);
  Vertex v4(4.15685424949238e+00, 4.53553390593274e+00, -1.66666666666667e+00);
  Vertex v5(-8.57106781186548e+00, -1.82842712474619e+00, 1.00000000000000e+00);
  Vertex v6(-7.15685424949238e+00, -2.53553390593274e+00, 1.00000000000000e+00);
  Vertex v7(5.57106781186548e+00, 3.82842712474619e+00, 1.00000000000000e+00);
  Vertex v8(4.15685424949238e+00, 4.53553390593274e+00, 1.00000000000000e+00);

  Face f1(v1, v4, v3, v2);
  Face f2(v5, v6, v7, v8);
  Face f3(v1, v5, v8, v4);
  Face f4(v2, v3, v7, v6);
  Face f5(v1, v2, v6, v5);
  Face f6(v3, v4, v8, v7);

  Vec n1( 0.00000000000000e+00, 0.00000000000000e+00, -1.00000000000000e+00);
  Vec n2( -0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e+00);
  Vec n3( -4.47213595499958e-01, 8.94427190999916e-01, 0.00000000000000e+00);
  Vec n4( 4.47213595499958e-01, -8.94427190999916e-01, 0.00000000000000e+00);
  Vec n5( -4.47213595499958e-01, -8.94427190999916e-01, 0.00000000000000e+00);
  Vec n6( 4.47213595499958e-01, 8.94427190999916e-01, 0.00000000000000e+00);

  Vec nn = f1.normal();
  EXPECT_NEAR(n1.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n1.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n1.z(), nn.z(), 1.0e-10);

  nn = f2.normal();
  EXPECT_NEAR(n2.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n2.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n2.z(), nn.z(), 1.0e-10);

  nn = f3.normal();
  EXPECT_NEAR(n3.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n3.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n3.z(), nn.z(), 1.0e-10);

  nn = f4.normal();
  EXPECT_NEAR(n4.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n4.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n4.z(), nn.z(), 1.0e-10);

  nn = f5.normal();
  EXPECT_NEAR(n5.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n5.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n5.z(), nn.z(), 1.0e-10);

  nn = f6.normal();
  EXPECT_NEAR(n6.x(), nn.x(), 1.0e-10);
  EXPECT_NEAR(n6.y(), nn.y(), 1.0e-10);
  EXPECT_NEAR(n6.z(), nn.z(), 1.0e-10);

  Point pt(0, 0, 0);
  auto td = f1.distance(pt);
  Point facePt = pt - f1.normal()*td.second;

  Point point(0.00000000000000e+00, 0.00000000000000e+00, -1.66666666666667e+00);
  REAL dist = -1.66666666666667e+00;
  std::array<REAL, 4> tt = {{
     5.96187269438804e-01,  -1.06862915010152e+00,   
     4.70479397227862e-01,  -3.33137084989848e+00}};
  std::array<REAL, 4> dd = {{
     2.19770311628353e+00,   5.91570658605306e+00,   
     9.32792052216177e-01,   5.46849299055311e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

  td = f2.distance(pt);
  facePt = pt - f2.normal()*td.second;

  point = Point( 0.00000000000000e+00, 0.00000000000000e+00, 1.00000000000000e+00);
  dist = -1.00000000000000e+00;
  tt = {{ 4.33137084989848e+00,   5.29520602772138e-01,   
          2.06862915010152e+00,   4.03812730561196e-01}};
  dd = {{ 5.46849299055311e+00,   9.32792052216177e-01,   
          5.91570658605306e+00,   2.19770311628353e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

  td = f3.distance(pt);
  facePt = pt - f3.normal()*td.second;

  point = Point( -9.82842712474619e-01, 1.96568542494924e+00, 0.00000000000000e+00);
  dist = -2.19770311628353e+00;
  tt = {{ 6.25000000000000e-01,   5.96187269438804e-01,   
          3.75000000000000e-01,   4.03812730561196e-01}};
  dd = {{ 8.48389357540403e+00,   1.00000000000000e+00,   
          5.74635589535368e+00,   1.66666666666667e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

  td = f4.distance(pt);
  facePt = pt - f4.normal()*td.second;

  point = Point( -4.17157287525381e-01, 8.34314575050762e-01, 0.00000000000000e+00);
  dist = 9.32792052216177e-01;
  tt = {{ 5.29520602772138e-01,   6.25000000000000e-01,   
          4.70479397227862e-01,   3.75000000000000e-01}};
  dd = {{ 1.66666666666667e+00,   6.69503919340420e+00,   
          1.00000000000000e+00,   7.53521027735352e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

  td = f5.distance(pt);
  facePt = pt - f5.normal()*td.second;

  point = Point( -2.44558441227157e+00, -4.89116882454314e+00, 0.00000000000000e+00);
  dist = -5.46849299055311e+00;
  tt = {{ 4.33137084989848e+00,   6.25000000000000e-01,  
          -3.33137084989848e+00,   3.75000000000000e-01}};
  dd = {{ 1.66666666666667e+00,   5.26735980818505e+00,   
          1.00000000000000e+00,   6.84849863826924e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

  td = f6.distance(pt);
  facePt = pt - f6.normal()*td.second;

  point = Point( 2.64558441227157e+00, 5.29116882454314e+00, 0.00000000000000e+00);
  dist = -5.91570658605306e+00;
  tt = {{ 2.06862915010152e+00,   6.25000000000000e-01,  
          -1.06862915010152e+00,   3.75000000000000e-01}};
  dd = {{ 1.66666666666667e+00,   1.68965104418539e+00,   
          1.00000000000000e+00,   3.27078987426957e+00}};
  EXPECT_NEAR(facePt.x(), point.x(), 1.0e-10);
  EXPECT_NEAR(facePt.y(), point.y(), 1.0e-10);
  EXPECT_NEAR(facePt.z(), point.z(), 1.0e-10);
  EXPECT_NEAR(td.second, dist, 1.0e-10);
  EXPECT_NEAR(td.first[0].first, tt[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].first, tt[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].first, tt[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].first, tt[3], 1.0e-10);
  EXPECT_NEAR(td.first[0].second, dd[0], 1.0e-10);
  EXPECT_NEAR(td.first[1].second, dd[1], 1.0e-10);
  EXPECT_NEAR(td.first[2].second, dd[2], 1.0e-10);
  EXPECT_NEAR(td.first[3].second, dd[3], 1.0e-10);

}