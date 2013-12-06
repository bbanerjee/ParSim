/*
 * FastMarchingMethod.cc
 *
 *  Created on: 5/12/2013
 *      Author: banerjee
 */

#include <Contact/FastMarching/FastMarchingMethod.h>

namespace BrMPM
{

FastMarchingMethod::FastMarchingMethod()
{
}

FastMarchingMethod::~FastMarchingMethod()
{
}

void FastMarchingMethod::extension_velocities(std::vector<double>& phi,
    std::vector<double>& speed, double* dx, bool self_test, int order,
    std::vector<int>& ext_mask, std::vector<double>& dd,
    std::vector<double>& f_ext)
{

  // Check input and add masks if needed
  preprocess(phi, dx, ext_mask);
}

} /* namespace BrMPM */
