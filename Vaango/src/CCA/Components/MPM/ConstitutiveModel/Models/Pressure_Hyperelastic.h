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

#ifndef __HYPERELASTIC_EOS_MODEL_H__
#define __HYPERELASTIC_EOS_MODEL_H__



#include <CCA/Components/MPM/ConstitutiveModel/Models/PressureModel.h> 
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelStateBase.h> 
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  ////////////////////////////////////////////////////////////////////////////
  /*! 
    \class Pressure_Hyperelastic
    \brief Hyperelastic relation for pressure from Simo and Hughes, 1998.

    The equation of state is given by
    \f[
    p = Tr(D) K \Delta T
    \f]
    where \n
    \f$p\f$ = pressure\n
    \f$D\f$ = rate of deformation tensor\n
    \f$K\f$ = bulk modulus\n
    \f$\Delta T\f$ = time increment
  */
  ////////////////////////////////////////////////////////////////////////////

  class Pressure_Hyperelastic : public PressureModel {

  private:

    // Prevent copying of this class
    // copy constructor
    //Pressure_Hyperelastic(const Pressure_Hyperelastic &cm);
    Pressure_Hyperelastic& operator=(const Pressure_Hyperelastic &cm);

  public:
    // constructors
    Pressure_Hyperelastic(); // This constructor is used when there is
                             // no equation_of_state tag in the input
                             // file  ** WARNING **
    Pressure_Hyperelastic(Uintah::ProblemSpecP& ps); 
    Pressure_Hyperelastic(const Pressure_Hyperelastic* cm);
         
    // destructor 
    virtual ~Pressure_Hyperelastic();

    virtual void outputProblemSpec(Uintah::ProblemSpecP& ps);
         
    /*! Get parameters */
    std::map<std::string, double> getParameters() const {
      std::map<std::string, double> params;
      params["K"] = d_bulk;
      return params;
    }

    //////////
    // Calculate the pressure using a equation of state
    double computePressure(const Uintah::MPMMaterial* matl,
                           const ModelStateBase* state,
                           const Uintah::Matrix3& deformGrad,
                           const Uintah::Matrix3& rateOfDeformation,
                           const double& delT);
  
    double eval_dp_dJ(const Uintah::MPMMaterial* matl,
                      const double& detF,
                      const ModelStateBase* state);

    // Compute bulk modulus
    double computeBulkModulus(const ModelStateBase* state) {return 0.0;};

    // Compute strain energy
    double computeStrainEnergy(const ModelStateBase* state) {return 0.0;};

    // Compute pressure (option 1)
    double computePressure(const double& rho_orig,
                           const double& rho_cur);

    // Compute pressure (option 2)
    void computePressure(const double& rho_orig,
                         const double& rho_cur,
                         double& pressure,
                         double& dp_drho,
                         double& csquared);

    // Compute bulk modulus
    double computeInitialBulkModulus();
    double computeBulkModulus(const double& rho_orig,
                              const double& rho_cur);

    // Compute strain energy
    double computeStrainEnergy(const double& rho_orig,
                               const double& rho_cur);

    // Compute density given pressure
    double computeDensity(const double& rho_orig,
                          const double& pressure);

    ////////////////////////////////////////////////////////////////////////
    /*! Calculate the derivative of p with respect to epse_v
        where epse_v = tr(epse)
              epse = total elastic strain */
    ////////////////////////////////////////////////////////////////////////
    double computeDpDepse_v(const ModelStateBase* state) const
    {
      return 0.0;
    };

    ////////////////////////////////////////////////////////////////////////
    /*! Calculate the derivative of p with respect to epse_s
        where epse_s = sqrt{2}{3} ||ee||
              ee = epse - 1/3 tr(epse) I
              epse = total elastic strain */
    ////////////////////////////////////////////////////////////////////////
    double computeDpDepse_s(const ModelStateBase* state) const
    {
      return 0.0;
    };
  };

} // End namespace Uintah

#endif  // __HYPERELASTIC_EOS_MODEL_H__ 
