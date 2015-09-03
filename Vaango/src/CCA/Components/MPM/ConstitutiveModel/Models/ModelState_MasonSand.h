/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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
 
    double capX;      // The cap hydrostatic compressive strength X 
    double kappa;     // The cap kappa parameter (branch point)
    double zeta;      // The back stress parameter (trace of isotropic backstress)

    Uintah::Matrix3* stressTensor;           // The tensor form of the unrotated stress
    Uintah::Matrix3  deviatoricStressTensor; // The deviatoric part of the stress
    double I1;        // I1 = Tr(sigma)
    double J2;
    double sqrt_J2;   // sqrt(J2) 

    Uintah::Matrix3* plasticStrainTensor;  // The tensor form of plastic strain
    double ev_p;      // ev_p = Tr(ep) : Volumetric part of the plastic strain

    double ev_0;      // Volumetric strain at zero pressure.  This is
                      // non-zero if the initial fluid pressure is non-zero

    double porosity;    // Porosity
    double saturation;  // Water saturation

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
      os << "I1 = " << state.I1 << "sqrt_J2 = " << state.sqrt_J2
         << ", evp = " << state.ev_p
         << ", X = " << state.capX << ", kappa = " << state.kappa
         << ", zeta = " << state.zeta 
         << ", phi = " << state.porosity << ", Sw = " << state.saturation << std::endl;
      return os;
    }
    
  };

} // End namespace Uintah

#endif  // __MODEL_STATE_ARENISCA3_PARTIALLY_SATURATED_H__ 