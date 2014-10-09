/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

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
Density::clone(const DensitySP den)
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
//  std::cout << "nodePos_y= " << nodePos;
  int num_of_periods = nodePos / lengthPeriod;
//  std::cout << "  num of periods= " << num_of_periods;
  double remind = nodePos - num_of_periods*lengthPeriod;
  remind *= 18.88;
//  std::cout << " remind=  " << remind << std::endl;
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
Density::nodeDensity (const NodeP node, double& node_density)
{
  Point3D node_pos = node->position();
  double node_pos_y = node_pos.y()-2 ;
  //   std::cout << "y of node= " << node_pos_y << "  ";
  //   std::cout << d_ring_width << "   ";
  double reminder = remind(d_ring_width, node_pos_y);
  node_density = density(d_poly_coeffs, reminder);
}


 



