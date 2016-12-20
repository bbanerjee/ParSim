/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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
 *  Weibull.h: support for Weibull distributions
 *
 *  Written by:
 *   Scot Swan
 *   Department of Mechanical Engineering
 *   University of Utah
 *   April 2009
 *  Modified by:
 *   Jim Guilkey
 *   Schlumberger
 *   September 2011
 *
 */


#ifndef SCI_WEIBULL_H__
#define SCI_WEIBULL_H__

#include <cmath>
#include <random>

#include <Core/Math/MusilRNG.h>
#include <Core/Math/share.h>

namespace Uintah {

class SCISHARE Weibull {

  public:
    double d_mean;               // Mean of the Weibull distributed variable
    double d_weibullModulus;
    double d_referenceVolume;
    double d_sizeEffectExponent;
    MusilRNG *d_uniformMusilRNG;
    std::mt19937 d_uniformMersenneRNG;                // The Mersenne twister RNG
    std::weibull_distribution<double> d_weibullDist;

    Weibull(double mean = 0,
            double weibullModulus = 1,
            double referenceVolume = 1,
            int seed = 0,
            double sizeEffectExponent = 1);

    ~Weibull();

    inline double rand(double particleVolume) {
   
      // Compute the volume scaling factor
      double C = std::pow(d_referenceVolume/particleVolume, 1./d_sizeEffectExponent);

      // Get the weibull distributed random variable value
      double weibull_val = d_weibullDist(d_uniformMersenneRNG);
     
      // Return scaled value
      return C*weibull_val;
    }

    // Equation taken from
    // "Theory and Results for Incorporating Aleatory Uncertainty and Scale
    // Effects in Damage Models for Failure and Fragmentation"
    // by R.M. Brannon and O.E. Strack, Sandia National Laboratories
    inline double randOld(double particleVolume) {

      // Get the uniformly distributed random #
      double y = (*d_uniformMusilRNG)();

      // Include a volume scaling factor
      double C = std::pow(d_referenceVolume/particleVolume, 1./d_sizeEffectExponent);

      double eta = d_mean/std::tgamma(1./d_weibullModulus + 1.0);
      //double eta = WeibMed_/pow(log(2.0),1./WeibMod_);

      // New version, easy to read and comprehend!
      return C*eta*std::pow(-log(1.0 - y), 1./d_weibullModulus);

      // Old way, hard to read
      //return WeibMed_*pow(log((*mr_)())/((PartVol/WeibRefVol_)*log(.5)),1/WeibMod_);
    }

    // Probability that x was picked from this Weibull distribution
    // found on http://www.weibull.com/LifeDataWeb/weibull_probability_density_function.htm
    double prob(double x, double particleVolume) {

      // The following is new and hopefully correct
      double C = pow(d_referenceVolume/particleVolume, 1./d_sizeEffectExponent);

      double eta = d_mean/std::tgamma(1./d_weibullModulus + 1.0);
      //double eta = WeibMed_/pow(log(2.0),1./WeibMod_);

      return d_weibullModulus/(C*eta)*std::pow(x/(C*eta), d_weibullModulus-1.)
                     *std::exp(-std::pow(x/(C*eta), d_weibullModulus));

      // Old and evidently wrong
      //   return (WeibMod_/(PartVol/WeibRefVol_))*
      //           pow((x-WeibMed_)/(PartVol/WeibRefVol_),WeibMod_-1)*
      //           exp(-pow((x-WeibMed_)/(PartVol/WeibRefVol_),WeibMod_));

    }
};

} // End namespace Uintah

#endif //SCI_WEIBULL_H__
