/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2023 Biswajit Banerjee
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

/*
 *  Weibull.cc: support choosing a random value from a 1D Weibull
 *               distribution (rand), as well as evaluate the probability
 *               of a particular value occuring.
 *
 *  Written by:
 *   Scot Swan
 *   Department of Mechanical Engineering
 *   University of Utah
 *   March 2009
 *
 */

#include <Core/Math/Weibull.h>

namespace Uintah {

Weibull::Weibull(double mean,
                 double weibullModulus,
                 double referenceVolume,
                 int seed,
                 double sizeEffectExponent) 
  : d_mean(mean),
    d_weibullModulus(weibullModulus),
    d_referenceVolume(referenceVolume),
    d_sizeEffectExponent(sizeEffectExponent),
    d_uniformMusilRNG(new MusilRNG(seed)),
    d_uniformMersenneRNG(seed)
{
  // Shape parameter
  double a = d_weibullModulus;
     
  // Scaling parameter
  double b = d_mean/std::tgamma(1./d_weibullModulus + 1.0);

  // Create the distribution
  d_weibullDist.param(std::weibull_distribution<double>::param_type(a, b));
}

Weibull::~Weibull() {
  delete d_uniformMusilRNG;
}

} // End namespace Uintah
