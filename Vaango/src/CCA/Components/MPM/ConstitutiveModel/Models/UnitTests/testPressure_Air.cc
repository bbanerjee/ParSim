#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Air.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Arena.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>

using namespace Vaango;



int main()
{
  
  Pressure_Air model;

  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  std::cout << "params[Ka] = " << params["Ka"]  << " Pa" << std::endl;

  // Set up densities
  double rho_orig  = 1.0;
  std::vector<double> rho_cur;
  for (unsigned int ii = 0; ii < 10; ii++) {
    double rho = rho_orig*(1.0 + (double) ii/10.0);
    rho_cur.emplace_back(rho);
    std::cout << "rho = " << rho << std::endl;
  }

  // Compute the pressure
  for (double rho : rho_cur) {
    double pp = model.computePressure(rho_orig, rho);  
    std::cout << "p = " << pp << std::endl;
  }

  // Set up list of pressures in log10 scale
  std::vector<double> pressures;
  for (unsigned int ii = 0; ii < 12; ii++) {
    pressures.emplace_back((double) ii);
  }

  // Convert to Pa and Compute bulk modulus
  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ka] = " << params["Ka"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ka] = " << params["Ka"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  ModelState_Arena state; 
  for (double pp : pressures) {
    state.I1_eff = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ka] = " << params["Ka"]  << " Pa" << std::endl;
    std::cout << "I1_eff = " << state.I1_eff << " K = " << K << std::endl;
  }
  
}
