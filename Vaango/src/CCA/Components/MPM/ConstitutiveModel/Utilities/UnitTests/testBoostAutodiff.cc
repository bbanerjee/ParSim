#include "autodiff.hpp"
#include <Eigen/Dense>

#include <iostream>

#include <gtest/gtest.h>

//using namespace Vaango;

using namespace autodiff_v1;

template <typename T>
T fourth_power(T const& x) {
  T x4 = x * x;  // retval in operator*() uses x4's memory via NRVO.
  x4 *= x4;      // No copies of x4 are made within operator*=() even when squaring.
  return x4;     // x4 uses y's memory in main() via NRVO.
}

template <typename W, typename X, typename Y, typename Z>
promote<W, X, Y, Z> 
f(const W& w, const X& x, const Y& y, const Z& z) {
  
  return exp(w * sin(x * log(y) / z) + sqrt(w * z / (x * y))) + w * w / tan(z);
}

TEST(BoostAutoDiffTest, scalars)
{
  constexpr unsigned Order = 5;                  // Highest order derivative to be calculated.
  auto const x = make_fvar<double, Order>(2.0);  // Find derivatives at x=2.
  auto const y = fourth_power(x);
  for (unsigned i = 0; i <= Order; ++i) {
    std::cout << "y.derivative(" << i << ") = " << y.derivative(i) << std::endl;
  }

  ASSERT_NEAR(y.derivative(0), 16., 1.0e-10);
}

TEST(BoostAutoDiffTest, multivariable)
{
  constexpr unsigned Nw = 3;  // Max order of derivative to calculate for w
  constexpr unsigned Nx = 2;  // Max order of derivative to calculate for x
  constexpr unsigned Ny = 4;  // Max order of derivative to calculate for y
  constexpr unsigned Nz = 3;  // Max order of derivative to calculate for z
  // Declare 4 independent variables together into a std::tuple.
  auto const variables = make_ftuple<double, Nw, Nx, Ny, Nz>(11.0, 12.0, 13.0, 14.0);
  auto const& w = std::get<0>(variables);  // Up to Nw derivatives at w=11
  auto const& x = std::get<1>(variables);  // Up to Nx derivatives at x=12
  auto const& y = std::get<2>(variables);  // Up to Ny derivatives at y=13
  auto const& z = std::get<3>(variables);  // Up to Nz derivatives at z=14
  auto const v = f(w, x, y, z);
  // Calculated from Mathematica symbolic differentiation.
  double const answer(1976.319600747797727779881875290418720908121189218755);
  ASSERT_NEAR(v.derivative(Nw, Nx, Ny, Nz), answer, 1.0e-11);
  std::cout << std::setprecision(std::numeric_limits<double>::digits)
            << "mathematica   : " << answer << '\n'
            << "autodiff      : " << v.derivative(Nw, Nx, Ny, Nz) << '\n'
            << std::setprecision(3)
            << "relative error: " << (v.derivative(Nw, Nx, Ny, Nz) / answer - 1) << '\n';
}

