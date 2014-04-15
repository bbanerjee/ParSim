#include <MaterialModels/Density.h>
#include <Pointers/DensitySP.h>
#include <Core/Node.h>

#include <Core/ProblemSpec/ProblemSpec.h>
#include <InputOutput/ProblemSpecUtil.h>

#include <iostream>
//#include <fstream>
//#include <string>
#include <cmath>
#include <vector>
#include <memory>

using namespace Matiti;

Density::Density()
  :d_ring_width(0.0)
{
  d_poly_coeffs.reserve(6);
}

Density::Density(const Density& den)
  :d_ring_width(den.d_ring_width)
{
  d_poly_coeffs = den.d_poly_coeffs;
} 

Density::~Density()
{
}

void
Density::clone(const DensitySP& den)
{
  d_ring_width = den->d_ring_width;
  d_poly_coeffs = den->d_poly_coeffs;
}

void
Density::initialize(Uintah::ProblemSpecP& ps)
{
  ps->require("ring_width", d_ring_width);
  
  //  Get the coefficients of density polynomial fitting function
  Uintah::ProblemSpecP poly_coeffs_ps = ps->findBlock("polynomial_coefficients");
  if (poly_coeffs_ps) {
    Matiti_ProblemSpecUtil::readVector(poly_coeffs_ps, d_poly_coeffs);
  } 
}


double
Density::remind (double lengthPeriod, double nodePos)
{
  std::cout << "nodePos_y= " << nodePos;
  int num_of_periods = nodePos / lengthPeriod;
  std::cout << "  num of periods= " << num_of_periods;
  double remind = nodePos - num_of_periods*lengthPeriod;
  remind *= 18.88;
  std::cout << " remind=  " << remind;
  return remind;
}


double
Density::density (const std::vector<double>& polyCoeff, double reminder)
{
  int power = 0;
  double density = 0.0;
  for (auto coeff_iter = polyCoeff.begin(); coeff_iter != polyCoeff.end() ; ++coeff_iter)
  {
    double cur_coeff = *coeff_iter;
    density += cur_coeff*pow(reminder, power);
    ++power;
  }
  return density;
} 


void
Density::nodeDensity (const NodeP& node, double& node_density)
{
  Point3D node_pos = node->position();
  double node_pos_y = node_pos.y()-2 ;
  //   std::cout << "y of node= " << node_pos_y << "  ";
  //   std::cout << d_ring_width << "   ";
  double reminder = remind(d_ring_width, node_pos_y);
  node_density = density(d_poly_coeffs, reminder);
}


 



