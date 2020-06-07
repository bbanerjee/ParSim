#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/PressureModels/Pressure_Air.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace Vaango;

TEST(PressureAirTest, computePressure)
{
  Pressure_Air model;

  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  //std::cout << "params[Ka] = " << params["Ka"] << " Pa" << std::endl;
  ASSERT_EQ(params["Ka"], model.computeInitialBulkModulus());

  // Set up densities
  double rho_orig = 1.0;
  std::vector<double> rho_cur;
  for (unsigned int ii = 0; ii < 10; ii++) {
    double rho = rho_orig * (1.0 + (double)ii / 10.0);
    rho_cur.emplace_back(rho);
    //std::cout << "rho = " << rho << std::endl;
  }

  // Compute the pressure
  std::map<double, double> rho_p;
  rho_p[1] = 0;
  rho_p[1.1] = 14463.8;
  rho_p[1.2] = 29463.7;
  rho_p[1.3] = 44972.6;
  rho_p[1.4] = 60966.5;
  rho_p[1.5] = 77424.3;
  rho_p[1.6] = 94327.1;
  rho_p[1.7] = 111658;
  rho_p[1.8] = 129402;
  rho_p[1.9] = 147544;
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);
    //std::cout << "rho = " << rho << " p = " << pp << std::endl;
    ASSERT_NEAR(pp, rho_p[rho], 1.0);
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (unsigned int ii = 0; ii < 12; ii++) {
    pressures.emplace_back((double)ii);
  }

  // Convert to Pa and Compute bulk modulus
  std::map<double, double> p_K;
  p_K[1] = 141856;
  p_K[10] = 141869;
  p_K[100] = 141995;
  p_K[1000] = 143255;
  p_K[10000] = 155855;
  p_K[100000] = 281855;
  p_K[1e+06] = 1.54186e+06;
  p_K[1e+07] = 1.41419e+07;
  p_K[1e+08] = 1.40141855e+08;
  p_K[1e+09] = 1.400141855e+09;
  p_K[1e+10] = 1.4000141855e+10;
  p_K[1e+11] = 1.40000141855e+11;
  for (double pp : pressures) {
    pp = std::pow(10, pp);
    params = model.getParameters();
    //std::cout << "Before: params[Ka] = " << params["Ka"] << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    //std::cout << "After: params[Ka] = " << params["Ka"] << " Pa" << std::endl;
    //std::cout << "p = " << pp << " K = " << K << std::endl;
    ASSERT_NEAR(K, p_K[pp], 1.0e2);
  }

  // Test model state input
  ModelState_Arena state;
  for (double pp : pressures) {
    state.pbar_w = std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    //std::cout << "After: params[Ka] = " << params["Ka"] << " Pa" << std::endl;
    //std::cout << "pbar_w = " << state.pbar_w << " K = " << K << std::endl;
    EXPECT_NEAR(K, p_K[state.pbar_w], 1.0e2);
  }
}
