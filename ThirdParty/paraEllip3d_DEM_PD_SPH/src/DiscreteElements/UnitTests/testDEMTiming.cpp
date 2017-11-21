#include <Core/Util/Utility.h>
#include <Core/Math/Vec.h>
#include <Core/Math/Matrix3.h>
#include <chrono>
#include <gtest/gtest.h>

using namespace dem;

class timer {
private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
public:
  void start() {
    start_time = std::chrono::high_resolution_clock::now();
  }
  double stop() {
    auto stop_time = std::chrono::high_resolution_clock::now();
    return double(std::chrono::duration_cast<std::chrono::nanoseconds>(stop_time - start_time).count());
  }
};

TEST(DEMParticleTimingTest, getRadius) {

  REAL ra = 3.1, rb = 2.1, rc = 1.1;
  REAL x = 2.9, y = 1.8, z = 0.9;
  REAL B = 0, C = 0;

  timer t;
  t.start();
  std::size_t iter = 0;
  const std::size_t N = 10000000;
  while (iter < N) {
    auto raSq = ra * ra;
    auto rbSq = rb * rb;
    auto rcSq = rc * rc;
    auto rcSq_raSq = rcSq/raSq;
    auto rcSq_rbSq = rcSq/rbSq;
    auto xSq = x * x;
    auto ySq = y * y;
    auto zSq = z * z;

    auto p = -(rcSq_raSq * x) / z;
    auto q = -(rcSq_rbSq * y) / z;
    auto r = -rcSq_raSq * (1 + rcSq_raSq * xSq / zSq) / z;
    auto t = -rcSq_rbSq * (1 + rcSq_rbSq * ySq / zSq) / z;
    auto s = -(rcSq * rcSq * x * y) / (raSq * rbSq * zSq * z);

    auto pSq = p * p; 
    auto qSq = q * q; 

    auto n = std::sqrt( 1 + pSq + qSq);

    B = n * ( 2 * p * q * s - (1 + pSq) * t - (1 + qSq) * r);
    C = n * n * n * n;
    ++iter;
  }
  double time = t.stop();
  std::cout << "New Time = " << time << " ns. (B, C) = " << B << ", " << C << std::endl;


  t.start();
  iter = 0;
  while (iter < N) {
    REAL p = -rc * rc / ra / ra * x / z;
    REAL q = -rc * rc / rb / rb * y / z;
    REAL r = -rc * rc / ra / ra * (1 / z + rc * rc / ra / ra * x * x / pow(z, 3));
    REAL t = -rc * rc / rb / rb * (1 / z + rc * rc / rb / rb * y * y / pow(z, 3));
    REAL s = -pow(rc, 4) / ra / ra / rb / rb * x * y / pow(z, 3);
    REAL n = sqrt(1 + p * p + q * q);
    B = n * (2 * p * q * s - (1 + p * p) * t - (1 + q * q) * r);
    C = n * n * n * n;
    ++iter;
  }
  time = t.stop();
  std::cout << "Old Time = " << time << " ns. (B, C) = " << B << ", " << C << std::endl;

  /*
  t.start();
  iter = 0;
  while (iter < N) {

    auto p = -(rc * rc * x) / (z * ra * ra);
    auto q = -(rc * rc * y) / (z * rb * rb);
    auto r = -rc * rc * (1 + rc * rc * x * x / (z * z * ra * ra)) / (z * ra * ra);
    auto t = -rc * rc * (1 + rc * rc * y * y / (z * z * rb * rb)) / (z * rb * rb);
    auto s = -(rc * rc * rc * rc * x * y) / (ra * ra * rb * rb * z * z * z);

    auto n = std::sqrt( 1 + p * p + q * q);

    B = n * ( 2 * p * q * s - (1 + p * p) * t - (1 + q * q) * r);
    C = n * n * n * n;
    ++iter;
  }
  time = t.stop();
  std::cout << "New Time Alt= " << time << " ns. (B, C) = " << B << ", " << C << std::endl;
  */
}

TEST(DEMParticleTimingTest, computeAndSetGlobalCoefOld)
{
  REAL d_a = 3.1, d_b = 2.1, d_c = 1.1;
  Vec d_currPos(5, 5, 5);
  REAL d_coef[10];
  Vec d_currDirecA(2.677339e+00,   1.121432e+00,   1.462010e+00);
  Vec d_currDirecB(1.122425e+00,   4.497600e-01,   1.603758e+00);
  Vec d_currDirecC(1.683123e+00,   1.553197e+00,   3.027884e+00);
  Vec v1 = vcos(d_currDirecA);
  Vec v2 = vcos(d_currDirecB);
  Vec v3 = vcos(d_currDirecC);

  timer t;
  t.start();
  std::size_t iter = 0;
  const std::size_t N = 10000000;
  while (iter < N) {

    REAL X0 = d_currPos.x();
    REAL Y0 = d_currPos.y();
    REAL Z0 = d_currPos.z();
    REAL l1 = v1.x();
    REAL m1 = v1.y();
    REAL n1 = v1.z();
    REAL l2 = v2.x();
    REAL m2 = v2.y();
    REAL n2 = v2.z();
    REAL l3 = v3.x();
    REAL m3 = v3.y();
    REAL n3 = v3.z();
    d_coef[0] = l1 * l1 / d_a / d_a + l2 * l2 / d_b / d_b + l3 * l3 / d_c / d_c;
    d_coef[1] = m1 * m1 / d_a / d_a + m2 * m2 / d_b / d_b + m3 * m3 / d_c / d_c;
    d_coef[2] = n1 * n1 / d_a / d_a + n2 * n2 / d_b / d_b + n3 * n3 / d_c / d_c;
    d_coef[3] = (2 * l1 * m1) / d_a / d_a + (2 * l2 * m2) / d_b / d_b +
                (2 * l3 * m3) / d_c / d_c;
    d_coef[4] = (2 * m1 * n1) / d_a / d_a + (2 * m2 * n2) / d_b / d_b +
                (2 * m3 * n3) / d_c / d_c;
    d_coef[5] = (2 * l1 * n1) / d_a / d_a + (2 * l2 * n2) / d_b / d_b +
                (2 * l3 * n3) / d_c / d_c;
    d_coef[6] =
      -2 * l1 * m1 * Y0 * pow(d_a, -2) - 2 * l1 * n1 * Z0 * pow(d_a, -2) -
      2 * l2 * m2 * Y0 * pow(d_b, -2) - 2 * l2 * n2 * Z0 * pow(d_b, -2) -
      2 * l3 * m3 * Y0 * pow(d_c, -2) - 2 * l3 * n3 * Z0 * pow(d_c, -2) -
      2 * X0 * pow(d_a, -2) * pow(l1, 2) - 2 * X0 * pow(d_b, -2) * pow(l2, 2) -
      2 * X0 * pow(d_c, -2) * pow(l3, 2);
    d_coef[7] = (-2 * l1 * m1 * X0) / d_a / d_a - (2 * l2 * m2 * X0) / d_b / d_b -
                (2 * l3 * m3 * X0) / d_c / d_c - (2 * m1 * m1 * Y0) / d_a / d_a -
                (2 * m2 * m2 * Y0) / d_b / d_b - (2 * m3 * m3 * Y0) / d_c / d_c -
                (2 * m1 * n1 * Z0) / d_a / d_a - (2 * m2 * n2 * Z0) / d_b / d_b -
                (2 * m3 * n3 * Z0) / d_c / d_c;
    d_coef[8] = (-2 * l1 * n1 * X0) / d_a / d_a - (2 * l2 * n2 * X0) / d_b / d_b -
                (2 * l3 * n3 * X0) / d_c / d_c - (2 * m1 * n1 * Y0) / d_a / d_a -
                (2 * m2 * n2 * Y0) / d_b / d_b - (2 * m3 * n3 * Y0) / d_c / d_c -
                (2 * n1 * n1 * Z0) / d_a / d_a - (2 * n2 * n2 * Z0) / d_b / d_b -
                (2 * n3 * n3 * Z0) / d_c / d_c;
    d_coef[9] = -1 + 2 * l1 * m1 * X0 * Y0 * pow(d_a, -2) +
                2 * l1 * n1 * X0 * Z0 * pow(d_a, -2) +
                2 * m1 * n1 * Y0 * Z0 * pow(d_a, -2) +
                2 * l2 * m2 * X0 * Y0 * pow(d_b, -2) +
                2 * l2 * n2 * X0 * Z0 * pow(d_b, -2) +
                2 * m2 * n2 * Y0 * Z0 * pow(d_b, -2) +
                2 * l3 * m3 * X0 * Y0 * pow(d_c, -2) +
                2 * l3 * n3 * X0 * Z0 * pow(d_c, -2) +
                2 * m3 * n3 * Y0 * Z0 * pow(d_c, -2) +
                pow(d_a, -2) * pow(l1, 2) * pow(X0, 2) +
                pow(d_b, -2) * pow(l2, 2) * pow(X0, 2) +
                pow(d_c, -2) * pow(l3, 2) * pow(X0, 2) +
                pow(d_a, -2) * pow(m1, 2) * pow(Y0, 2) +
                pow(d_b, -2) * pow(m2, 2) * pow(Y0, 2) +
                pow(d_c, -2) * pow(m3, 2) * pow(Y0, 2) +
                pow(d_a, -2) * pow(n1, 2) * pow(Z0, 2) +
                pow(d_b, -2) * pow(n2, 2) * pow(Z0, 2) +
                pow(d_c, -2) * pow(n3, 2) * pow(Z0, 2);

    REAL divd = d_coef[0];
    for (auto& coef : d_coef) {
      coef /= divd;
    }
    ++iter;
  }
  auto time = t.stop();
  std::cout << "Old Time = " << time << " ns. Coefs = ";
  for (double& i : d_coef) std::cout << i << ", ";
  std::cout << std::endl;
}

TEST(DEMParticleTimingTest, computeAndSetGlobalCoefOldNoPow)
{
  REAL d_a = 3.1, d_b = 2.1, d_c = 1.1;
  Vec d_currPos(5, 5, 5);
  REAL d_coef[10];
  Vec d_currDirecA(2.677339e+00,   1.121432e+00,   1.462010e+00);
  Vec d_currDirecB(1.122425e+00,   4.497600e-01,   1.603758e+00);
  Vec d_currDirecC(1.683123e+00,   1.553197e+00,   3.027884e+00);
  Vec v1 = vcos(d_currDirecA);
  Vec v2 = vcos(d_currDirecB);
  Vec v3 = vcos(d_currDirecC);

  timer t;
  t.start();
  std::size_t iter = 0;
  const std::size_t N = 10000000;
  while (iter < N) {


    REAL X0 = d_currPos.x();
    REAL Y0 = d_currPos.y();
    REAL Z0 = d_currPos.z();

    REAL l1 = v1.x() / d_a;
    REAL m1 = v1.y() / d_a;
    REAL n1 = v1.z() / d_a;
    REAL l2 = v2.x() / d_b;
    REAL m2 = v2.y() / d_b;
    REAL n2 = v2.z() / d_b;
    REAL l3 = v3.x() / d_c;
    REAL m3 = v3.y() / d_c;
    REAL n3 = v3.z() / d_c;

    auto xsq = l1 * l1 + l2 * l2 + l3 * l3 ;
    auto ysq = m1 * m1 + m2 * m2 + m3 * m3 ;
    auto zsq = n1 * n1 + n2 * n2 + n3 * n3 ;
    auto xy = 2 * (l1 * m1 +  l2 * m2 + l3 * m3) ;
    auto yz = 2 * (m1 * n1 +  m2 * n2 + m3 * n3) ;
    auto zx = 2 * (l1 * n1 +  l2 * n2 + l3 * n3) ;
    auto x = 2 * ( -l1 * m1 * Y0 - l1 * n1 * Z0 -
                    l2 * m2 * Y0 - l2 * n2 * Z0 - 
                    l3 * m3 * Y0 - l3 * n3 * Z0 -
                    X0 * (l1 * l1) - X0 * (l2 * l2) - X0 * (l3 * l3) );
    auto y = 2 * ( (-l1 * m1 * X0) - ( l2 * m2 * X0) -
                   ( l3 * m3 * X0) - ( m1 * m1 * Y0) -
                   ( m2 * m2 * Y0) - ( m3 * m3 * Y0) -
                   ( m1 * n1 * Z0) - ( m2 * n2 * Z0) -
                   ( m3 * n3 * Z0) );
    auto z = 2 * ( (-l1 * n1 * X0) - ( l2 * n2 * X0) -
                   ( l3 * n3 * X0) - ( m1 * n1 * Y0) -
                   ( m2 * n2 * Y0) - ( m3 * n3 * Y0) -
                   ( n1 * n1 * Z0) - ( n2 * n2 * Z0) -
                   ( n3 * n3 * Z0) );
    auto c = -1 + 2 * (l1 * m1 * X0 * Y0 +
                       l1 * n1 * X0 * Z0 +
                       m1 * n1 * Y0 * Z0 +
                       l2 * m2 * X0 * Y0 +
                       l2 * n2 * X0 * Z0 +
                       m2 * n2 * Y0 * Z0 +
                       l3 * m3 * X0 * Y0 +
                       l3 * n3 * X0 * Z0 +
                       m3 * n3 * Y0 * Z0) +
                      (l1 * l1) * (X0 * X0) +
                      (l2 * l2) * (X0 * X0) +
                      (l3 * l3) * (X0 * X0) +
                      (m1 * m1) * (Y0 * Y0) +
                      (m2 * m2) * (Y0 * Y0) +
                      (m3 * m3) * (Y0 * Y0) +
                      (n1 * n1) * (Z0 * Z0) +
                      (n2 * n2) * (Z0 * Z0) +
                      (n3 * n3) * (Z0 * Z0);

    REAL fac = 1.0/xsq;
    d_coef[0] = 1;
    d_coef[1] = ysq * fac;
    d_coef[2] = zsq * fac;
    d_coef[3] = xy * fac;
    d_coef[4] = yz * fac;
    d_coef[5] = zx * fac;
    d_coef[6] = x * fac;
    d_coef[7] = y * fac;
    d_coef[8] = z * fac;
    d_coef[9] = c * fac;
    ++iter;
  }

  auto time = t.stop();

  std::cout << "New Time = " << time << " ns. Coefs = ";
  for (double& i : d_coef) std::cout << i << ", ";
  std::cout << std::endl;
}

/*
TEST(DEMParticleTimingTest, computeAndSetGlobalCoefNew)
{
  REAL d_a = 3.1, d_b = 2.1, d_c = 1.1;
  Vec d_currPos(5, 5, 5);
  REAL d_coef[10];
  for (double& coef : d_coef) coef = 0;
  Vec d_currDirecA(2.677339e+00,   1.121432e+00,   1.462010e+00);
  Vec d_currDirecB(1.122425e+00,   4.497600e-01,   1.603758e+00);
  Vec d_currDirecC(1.683123e+00,   1.553197e+00,   3.027884e+00);
  Vec v1 = vcos(d_currDirecA);
  Vec v2 = vcos(d_currDirecB);
  Vec v3 = vcos(d_currDirecC);

  timer t;
  t.start();
  std::size_t iter = 0;
  const std::size_t N = 10000000;

  while (iter < N) {
    Vec v11 = v1 * d_a; Vec v22 = v2 * d_b; Vec v33 = v3 * d_c;

    //Matrix3 L(v11.x(), v11.y(), v11.z(), 
    //          v22.x(), v22.y(), v22.z(),
    //          v33.x(), v33.y(), v33.z());
    //auto f123 = L.Determinant();
    //auto f23_x = L.Minor(0, 0);
    //auto f23_y = -L.Minor(0, 1);
    //auto f23_z = L.Minor(0, 2);
    //auto f13_x = L.Minor(1, 0);
    //auto f13_y = -L.Minor(1, 1);
    //auto f13_z = L.Minor(1, 2);
    //auto f12_x = L.Minor(2, 0);
    //auto f12_y = -L.Minor(2, 1);
    //auto f12_z = L.Minor(2, 2);

    //REAL X0 = d_currPos.x();
    //REAL Y0 = d_currPos.y();
    //REAL Z0 = d_currPos.z();

    REAL v1x = v11.x();
    REAL v1y = v11.y();
    REAL v1z = v11.z();
    REAL v2x = v22.x();
    REAL v2y = v22.y();
    REAL v2z = v22.z();
    REAL v3x = v33.x();
    REAL v3y = v33.y();
    REAL v3z = v33.z();

    REAL f23_x = (v2y * v3z) - (v2z * v3y);
    REAL f23_y = -(v2x * v3z) + (v2z * v3x);
    REAL f23_z = (v2x * v3y) - (v2y * v3x);
    REAL f13_x = (v1y * v3z) - (v1z * v3y);
    REAL f13_y = -(v1x * v3z) + (v1z * v3x);
    REAL f13_z = (v1x * v3y) - (v1y * v3x);
    REAL f12_x = (v1y * v2z) - (v1z * v2y);
    REAL f12_y = -(v1x * v2z) + (v1z * v2x);
    REAL f12_z = (v1x * v2y) - (v1y * v2x);

    REAL xsq = f23_x * f23_x + f12_x * f12_x +  f13_x * f13_x;
    REAL ysq = f23_y * f23_y + f12_y * f12_y +  f13_y * f13_y;
    REAL zsq = f23_z * f23_z + f12_z * f12_z +  f13_z * f13_z;

    REAL xy = 2 * (f23_x * f23_y + f12_x * f12_y + f13_x * f13_y);
    REAL yz = 2 * (f23_y * f23_z + f12_y * f12_z + f13_y * f13_z);
    REAL zx = 2 * (f23_z * f23_x + f12_z * f12_x + f13_z * f13_x);

    //std::cout << "f23_x = " << f23_x << "\n";
    //d_coef[0] = (f23_x * f23_x) + (f12_x * f12_x) +  (f13_x * f13_x);
    //d_coef[1] = (f23_y * f23_y) + (f12_y * f12_y) +  (f13_y * f13_y);
    //d_coef[2] = (f23_z * f23_z) + (f12_z * f12_z) +  (f13_z * f13_z);

    //d_coef[3] = 2 * ((f23_x * f23_y) + (f12_x * f12_y) + (f13_x * f13_y));
    //d_coef[4] = 2 * ((f23_y * f23_z) + (f12_y * f12_z) + (f13_y * f13_z));
    //d_coef[5] = 2 * ((f23_z * f23_x) + (f12_z * f12_x) + (f13_z * f13_x));

    REAL fac = 1.0/xsq;

    d_coef[0] = 1;
    d_coef[1] = ysq * fac;
    d_coef[2] = zsq * fac;
    d_coef[3] = xy * fac;
    d_coef[4] = yz * fac;
    d_coef[5] = zx * fac;

    ++iter;
  }
  auto time = t.stop();

  //REAL divd = d_coef[0];
  //for (double& i : d_coef) i /= divd;

  std::cout << "New Time = " << time << " ns. Coefs = ";
  for (double& i : d_coef) std::cout << i << ", ";
  std::cout << std::endl;
}

TEST(DEMParticleTimingTest, computeAndSetGlobalCoefNewExtra)
{
  REAL d_a = 3.1, d_b = 2.1, d_c = 1.1;
  Vec d_currPos(5, 5, 5);
  REAL d_coef[10];
  Vec d_currDirecA(2.677339e+00,   1.121432e+00,   1.462010e+00);
  Vec d_currDirecB(1.122425e+00,   4.497600e-01,   1.603758e+00);
  Vec d_currDirecC(1.683123e+00,   1.553197e+00,   3.027884e+00);
  Vec v1 = vcos(d_currDirecA);
  Vec v2 = vcos(d_currDirecB);
  Vec v3 = vcos(d_currDirecC);
  for (double& i : d_coef) i = 0;

  timer t;
  t.start();
  std::size_t iter = 0;
  const std::size_t N = 10000000;
  while (iter < N) {

    Vec v11 = v1 * d_a; Vec v22 = v2 * d_b; Vec v33 = v3 * d_c;

    REAL v1x = v11.x();
    REAL v1y = v11.y();
    REAL v1z = v11.z();
    REAL v2x = v22.x();
    REAL v2y = v22.y();
    REAL v2z = v22.z();
    REAL v3x = v33.x();
    REAL v3y = v33.y();
    REAL v3z = v33.z();

    d_coef[0] = (v2y * v3z - v2z * v3y) * (v2y * v3z - v2z * v3y) + 
                (v1y * v2z - v1z * v2y) * (v1y * v2z - v1z * v2y) +  
                (v1y * v3z - v1z * v3y) * (v1y * v3z - v1z * v3y);
    d_coef[1] = (v2x * v3z - v2z * v3x) * (v2x * v3z - v2z * v3x) + 
                (v1x * v2z - v1z * v2x)* (v1x * v2z - v1z * v2x) +  
                (v1x * v3z - v1z * v3x) * (v1x * v3z - v1z * v3x);
    d_coef[2] = (v2x * v3y - v2y * v3x) * (v2x * v3y - v2y * v3x) + 
                (v1x * v2y - v1y * v2x) * (v1x * v2y - v1y * v2x) +  
                (v1x * v3y - v1y * v3x) * (v1x * v3y - v1y * v3x);

    d_coef[3] = 2 * (-(v2y * v3z - v2z * v3y) * (v2x * v3z - v2z * v3x) - 
                (v1y * v2z - v1z * v2y) * (v1x * v2z - v1z * v2x) 
                - (v1y * v3z - v1z * v3y) * (v1x * v3z - v1z * v3x));
    d_coef[4] = 2 * (-(v2x * v3z - v2z * v3x) * (v2x * v3y - v2y * v3x) - 
                (v1x * v2z - v1z * v2x) * (v1x * v2y - v1y * v2x) - 
                (v1x * v3z - v1z * v3x) * (v1x * v3y - v1y * v3x));
    d_coef[5] = 2 * ((v2x * v3y - v2y * v3x) * (v2y * v3z - v2z * v3y) + 
                (v1x * v2y - v1y * v2x) * (v1y * v2z - v1z * v2y) + 
                (v1x * v3y - v1y * v3x) * (v1y * v3z - v1z * v3y));

    ++iter;
  }

  auto time = t.stop();

    REAL divd = 1.0/d_coef[0];
    for (double& i : d_coef) i *= divd;

  std::cout << "New Time = " << time << " ns. Coefs = ";
  for (double& i : d_coef) std::cout << i << ", ";
  std::cout << std::endl;
}
*/