#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Granite.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>

using namespace Vaango;



int main()
{
  
  Pressure_Granite model;

  // Get initial parameters
  std::map<std::string, double> params = model.getParameters();
  std::cout << "params[Ks] = " << params["Ks"]  << " Pa" << std::endl;

  // Set up densities
  double rho_orig  = 1.0;
  std::vector<double> rho_cur;
  for (int ii = -2; ii < 10; ii++) {
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
  for (int ii = -2; ii < 12; ii++) {
    pressures.emplace_back((double) ii);
  }

  // Convert to Pa and Compute bulk modulus
  std::cout << "Compression:" << std::endl;
  for (double pp : pressures) {
    pp = std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  ModelState_MasonSand state; 
  for (double pp : pressures) {
    state.I1 = -std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1 = " << state.I1 << " K = " << K << std::endl;
  }
  
  // Test tension states
  // Convert to Pa and Compute bulk modulus
  std::cout << "Tension:" << std::endl;
  for (double pp : pressures) {
    pp = -std::pow(10, pp);

    params = model.getParameters();
    std::cout << "Before: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    double K = model.computeBulkModulus(pp);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "p = " << pp << " K = " << K << std::endl;
  }

  // Test model state input
  for (double pp : pressures) {
    state.I1 = std::pow(10, pp);
    double K = model.computeBulkModulus(&state);
    params = model.getParameters();
    std::cout << "After: params[Ks] = " << params["Ks"]  << " Pa" << std::endl;
    std::cout << "I1 = " << state.I1 << " K = " << K << std::endl;
  }
  
}
