#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/WaterEOS.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <gtest/gtest.h>

using namespace Vaango;

TEST(PressureWaterTest, computePressure)
{

  WaterEOS model;

  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  //std::cout << "params[Kw] = " << params["Kw"] << " Pa" << std::endl;
  ASSERT_EQ(params["Kw"], model.computeInitialBulkModulus());

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
  rho_p[1] = 0.0001;
  rho_p[1.1] = 2.84622e+08;
  rho_p[1.2] = 7.33788e+08;
  rho_p[1.3] = 1.41627e+09;
  rho_p[1.4] = 2.42054e+09;
  rho_p[1.5] = 3.85819e+09;
  rho_p[1.6] = 5.86772e+09;
  rho_p[1.7] = 8.61855e+09;
  rho_p[1.8] = 1.23154e+10;
  rho_p[1.9] = 1.72027e+10;
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);
    //std::cout << "rho = " << rho << " p = " << pp << std::endl;
    ASSERT_NEAR(pp, rho_p[rho], 1.0e5);
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (unsigned int ii = 0; ii < 12; ii++) {
    pressures.emplace_back((double)ii);
  }

  // Convert to Pa and Compute bulk modulus
  std::map<double, double> p_K;
  p_K[1] = 2.21e+09;
  p_K[10] = 2.21e+09;
  p_K[100] = 2.21e+09;
  p_K[1000] = 2.21001e+09;
  p_K[10000] = 2.21006e+09;
  p_K[100000] = 2.2106e+09;
  p_K[1e+06] = 2.21603e+09;
  p_K[1e+07] = 2.27029e+09;
  p_K[1e+08] = 2.8129e+09;
  p_K[1e+09] = 8.239e+09;
  p_K[1e+10] = 6.25e+10;
  p_K[1e+11] = 6.0511e+11;

  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    //std::cout << "Before: params[Kw] = " << params["Kw"] << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    //std::cout << "After: params[Kw] = " << params["Kw"] << " Pa" << std::endl;
    //std::cout << "p = " << pp << " K = " << K << std::endl;
    ASSERT_NEAR(K, p_K[pp], 1.0e6);
  }

  // Test model state input
  ModelState_Arena state;
  for (double pp : pressures) {
    state.pbar_w = std::pow(10, pp);
    state.updateStressInvariants(); 
    //state.I1_eff = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    //std::cout << "After: params[Kw] = " << params["Kw"] << " Pa" << std::endl;
    //std::cout << "pbar_w = " << state.pbar_w << " K = " << K << std::endl;
    ASSERT_NEAR(K, p_K[state.pbar_w], 1.0e6);
  }
}
