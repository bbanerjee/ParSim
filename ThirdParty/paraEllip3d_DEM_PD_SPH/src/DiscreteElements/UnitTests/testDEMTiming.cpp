#include <Core/Util/Utility.h>
#include <chrono>
#include <gtest/gtest.h>


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
  std::cout << "Orig Time = " << time << " ns. (B, C) = " << B << ", " << C << std::endl;

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
}

