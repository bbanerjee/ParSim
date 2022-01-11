/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#ifndef __ZERILLI_ARMSTRONG_POLYMER_MODEL_H__
#define __ZERILLI_ARMSTRONG_POLYMER_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/FlowStressModel.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

////////////////////////////////////////////////////////////////////////////////
/*!
  \class ZAPolymerFlow
  \brief Zerilli-Armstrong
  \author Todd Harman/Scott Bardenhagen
  Department of Mechanical Engineering,
  University of Utah

  <<<< Scott: Please add a reference and any documention>>>>>>

*/
/////////////////////////////////////////////////////////////////////////////

class ZAPolymerFlow : public FlowStressModel
{

public:
  // Create datatype for storing model parameters
  struct CMData
  {
    double sigma_g;
    double B_pa;
    double B_pb;
    double B_pn;
    double beta_0;
    double beta_1;
    double T_0;
    double B_0pa;
    double B_0pb;
    double B_0pn;
    double omega_a;
    double omega_b;
    double omega_p;
    double alpha_0;
    double alpha_1;
  };

private:
  CMData d_CM;

  // Prevent copying of this class
  // copy constructor
  ZAPolymerFlow& operator=(const ZAPolymerFlow& cm);

public:
  // constructors
  ZAPolymerFlow(ProblemSpecP& ps);
  ZAPolymerFlow(const ZAPolymerFlow* cm);

  // destructor
  ~ZAPolymerFlow() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  ///////////////////////////////////////////////////////////////////////////
  /*! \brief  compute the flow stress */
  ///////////////////////////////////////////////////////////////////////////
  double computeFlowStress(const ModelStateBase* state, const double& delT,
                           const double& tolerance, const MPMMaterial* matl,
                           const particleIndex idx) override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to plastic strain
  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTPlasticStrain(const ModelStateBase* state,
                                        const particleIndex idx) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to strain rate.
  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTStrainRate(const ModelStateBase* state,
                                     const particleIndex idx) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the shear modulus.
  */
  ///////////////////////////////////////////////////////////////////////////
  double computeShearModulus(const ModelStateBase* state) override;

  //______________________________________________________________________
  //  Empty functions
  double computeEpdot(const ModelStateBase*, const double&, const double&,
                      const MPMMaterial*, const particleIndex) override
  {
    throw InternalError("ZAPolymerFlow::ComputeEpdot has not been implemented",
                        __FILE__, __LINE__);
  };

  double computeMeltingTemp(const ModelStateBase*) override
  {
    throw InternalError(
      "ZAPolymerFlow::computeMeltingTemp has not been implemented", __FILE__,
      __LINE__);
  };

  void evalDerivativeWRTScalarVars(const ModelStateBase*, const particleIndex,
                                   Vector&) const override
  {
    throw InternalError(
      "ZAPolymerFlow::evalDerivativeWRTScalarVars has not been implemented",
      __FILE__, __LINE__);
  }

  double evalDerivativeWRTTemperature(const ModelStateBase*,
                                      const particleIndex) const
  {
    throw InternalError(
      "ZAPolymerFlow::evalDerivativeWRTTemperature has not been implemented",
      __FILE__, __LINE__);
  }
};

} // End namespace Uintah

#endif // __ZERILLI_ARMSTRONG_POLYMER_MODEL_H__
