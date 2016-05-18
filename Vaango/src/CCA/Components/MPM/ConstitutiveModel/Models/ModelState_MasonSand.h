/*
 * The MIT License
 *
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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

#ifndef __MODEL_STATE_ARENISCA3_PARTIALLY_SATURATED_H__
#define __MODEL_STATE_ARENISCA3_PARTIALLY_SATURATED_H__

#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Default.h>

namespace Vaango {

  /////////////////////////////////////////////////////////////////////////////
  /*!
    \class ModelState_MasonSand
    \brief A structure that stores the state data that is specialized for
           the PartallySaturated model.
           ** Derived from PlasticityState:ModelState
    \author Biswajit Banerjee \n
  */
  /////////////////////////////////////////////////////////////////////////////

  class ModelState_MasonSand: public ModelState_Default {

  public:

    static const Uintah::Matrix3 Identity;
    static const double sqrtTwo;
    static const double sqrtThree;
 
    double capX;      // The cap hydrostatic compressive strength X 
    double kappa;     // The cap kappa parameter (branch point)
    double pbar_w;    // The back stress parameter (-1/3*trace of isotropic backstress)

    Uintah::Matrix3 stressTensor;           // The tensor form of the total stress
    Uintah::Matrix3 deviatoricStressTensor; // The deviatoric part of the total stress
    double I1_eff;    // I1_eff = Tr(sigma_eff) = Tr(sigma) + 3*pbar_w
    double J2;
    double sqrt_J2;   // sqrt(J2) 
    double rr;        // Lode coordinate 'r'
    double zz_eff;    // Lode coordinate 'z'

    Uintah::Matrix3 plasticStrainTensor;  // The tensor form of plastic strain
    double ep_v;      // ep_v = Tr(ep) : Volumetric part of the plastic strain
    double dep_v;     // Increment of the volumetric plastic strain

    double ev_0;      // Volumetric strain at zero pressure.  This is
                      // non-zero if the initial fluid pressure is non-zero

    double phi0;        // Initial porosity
    double Sw0;         // Initial saturation
    double saturation;  // Water saturation

    double p3;        // P3 used by disaggregation algorithm

    std::vector<double> yieldParams;  // The yield parameters for a single particle
                                      // (variability)

    // Defined in base class
    // double bulkModulus;   // Bulk and shear moduli
    // double shearModulus;
    // double porosity;    // Porosity
    // Matrix3 backStress; // Back stress

    ModelState_MasonSand();

    ModelState_MasonSand(const ModelState_MasonSand& state);
    ModelState_MasonSand(const ModelState_MasonSand* state);

    ~ModelState_MasonSand();

    ModelState_MasonSand& operator=(const ModelState_MasonSand& state);
    ModelState_MasonSand* operator=(const ModelState_MasonSand* state);

    void updateStressInvariants();
    void updateVolumetricPlasticStrain();

    friend std::ostream& operator<<(std::ostream& os, 
                                    const ModelState_MasonSand& state) {
      os << "\t I1_eff = " << state.I1_eff << ", sqrt_J2 = " << state.sqrt_J2
         << ", r = " << state.rr << ", z_eff = " << state.zz_eff
         << ", evp = " << state.ep_v << ", p3 = " << state.p3 << "\n"
         << "\t K = " << state.bulkModulus << ", G = " << state.shearModulus << "\n"
         << "\t X = " << state.capX << ", kappa = " << state.kappa
         << ", pbar_w = " << state.pbar_w  << "\n"
         << "\t phi = " << state.porosity << ", Sw = " << state.saturation << std::endl;
      os << "\t Yield parameters: ";
      for (double val : state.yieldParams) {
         os << val << ", ";
      }
      os << std::endl;
      return os;
    }
    
  };

} // End namespace Uintah

#endif  // __MODEL_STATE_ARENISCA3_PARTIALLY_SATURATED_H__ 
