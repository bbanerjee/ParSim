#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/GraniteEOS.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <gtest/gtest.h>

using namespace Vaango;

TEST(PressureGraniteTest, computePressure)
{

  GraniteEOS model;

  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  //std::cout << "params[Ks] = " << params["Ks"] << " Pa" << std::endl;
  ASSERT_EQ(params["Ks"], model.computeInitialBulkModulus());

  // Set up densities
  double rho_orig = 1.0;
  std::vector<double> rho_cur;
  for (int ii = -2; ii < 10; ii++) {
    double rho = rho_orig * (1.0 + (double)ii / 10.0);
    rho_cur.emplace_back(rho);
    //std::cout << "rho = " << rho << std::endl;
  }

  // Compute the pressure
  std::map<double, double> rho_p;
  rho_p[0.8] = -5.9039e+09;
  rho_p[0.9] = -3.4389e+09;
  rho_p[1] = 101325;
  rho_p[1.1] = 4.6411e+09;
  rho_p[1.2] = 1.07361e+10;
  rho_p[1.3] = 1.85611e+10;
  rho_p[1.4] = 2.84161e+10;
  rho_p[1.5] = 4.06251e+10;
  rho_p[1.6] = 5.55361e+10;
  rho_p[1.7] = 7.35211e+10;
  rho_p[1.8] = 9.49761e+10;
  rho_p[1.9] = 1.20321e+11;
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);
    //std::cout << "rho = " << rho << " p = " << pp << std::endl;
    ASSERT_NEAR(pp, rho_p[rho], 1.0e6);
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (int ii = -2; ii < 12; ii++) {
    pressures.emplace_back((double)ii);
  }

  // Convert to Pa and Compute bulk modulus
  //std::cout << "Compression:" << std::endl;
  std::map<double, double> p_K;
  p_K[0.01] = 3.99996e+10;
  p_K[0.1] = 3.99996e+10;
  p_K[1] = 3.99996e+10;
  p_K[10] = 3.99996e+10;
  p_K[100] = 3.99996e+10;
  p_K[1000] = 3.99996e+10;
  p_K[10000] = 3.99996e+10;
  p_K[100000] = 4e+10;
  p_K[1e+06] = 4.00036e+10;
  p_K[1e+07] = 4.00396e+10;
  p_K[1e+08] = 4.03996e+10;
  p_K[1e+09] = 4.39996e+10;
  p_K[1e+10] = 7.99996e+10;
  p_K[1e+11] = 4.4e+11;
  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    //std::cout << "Before: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    //std::cout << "After: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    //std::cout << "p = " << pp << " K = " << K << std::endl;
    ASSERT_NEAR(K, p_K[pp], 1.0e6);
  }

  // Test model state input
  std::map<double, double> I1_K;
  I1_K[-0.01] = 3.99996e+10;
  I1_K[-0.1] = 3.99996e+10;
  I1_K[-1] = 3.99996e+10;
  I1_K[-10] = 3.99996e+10;
  I1_K[-100] = 3.99996e+10;
  I1_K[-1000] = 3.99996e+10;
  I1_K[-10000] = 3.99996e+10;
  I1_K[-100000] = 3.99997e+10;
  I1_K[-1e+06] = 4.00009e+10;
  I1_K[-1e+07] = 4.00129e+10;
  I1_K[-1e+08] = 4.01329e+10;
  I1_K[-1e+09] = 4.13329e+10;
  I1_K[-1e+10] = 5.33329e+10;
  I1_K[-1e+11] = 1.73333e+11;

  ModelState_Arena state;
  for (double pp : pressures) {
    state.I1_eff = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    //std::cout << "After: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    //std::cout << "I1_eff = " << state.I1_eff << " K = " << K << std::endl;
    ASSERT_NEAR(K, I1_K[state.I1_eff], 1.0e6);
  }

  // Test tension states
  // Convert to Pa and Compute bulk modulus
  //std::cout << "Tension:" << std::endl;
  for (double pp : pressures) {
    pp = -std::pow(10, pp);

    params = model.getParameters();
    //std::cout << "Before: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    //std::cout << "After: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    //std::cout << "p = " << pp << " K = " << K << std::endl;
    ASSERT_NEAR(K, 4.0e10, 1.0e6);
  }

  // Test model state input
  for (double pp : pressures) {
    state.I1_eff = std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    //std::cout << "After: params[Ks] = " << params["Ks"] << " Pa" << std::endl;
    //std::cout << "I1_eff = " << state.I1_eff << " K = " << K << std::endl;
    EXPECT_NEAR(K, 4.0e10, 1.0e6);
  }
}
