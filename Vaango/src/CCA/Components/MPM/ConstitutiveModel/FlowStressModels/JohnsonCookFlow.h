/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#ifndef __JOHNSONCOOK_FLOW_MODEL_H__
#define __JOHNSONCOOK_FLOW_MODEL_H__

#include "FlowStressModel.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class JohnsonCookFlow
  \brief Johnson-Cook Strain rate dependent plasticity model
  \author Biswajit Banerjee,
  Department of Mechanical Engineering,
  University of Utah

  Johnson-Cook Plasticity Model \n
  (Johnson and Cook, 1983, Proc. 7th Intl. Symp. Ballistics, The Hague) \n
  The flow rule is given by
  \f[
  f(\sigma) = [A + B (\epsilon_p)^n][1 + C \ln(\dot{\epsilon_p^*})]
  [1 - (T^*)^m]
  \f]

  where \f$ f(\sigma)\f$  = equivalent stress \n
  \f$ \epsilon_p\f$  = plastic strain \n
  \f$ \dot{\epsilon_p^{*}} = \dot{\epsilon_p}/\dot{\epsilon_{p0}}\f$
  where \f$ \dot{\epsilon_{p0}}\f$  = a user defined plastic
  strain rate,  \n
  A, B, C, n, m are material constants \n
  (for HY-100 steel tubes :
  A = 316 MPa, B = 1067 MPa, C = 0.0277, n = 0.107, m = 0.7) \n
  A is interpreted as the initial yield stress - \f$ \sigma_0 \f$ \n
  \f$ T^* = (T-T_{room})/(T_{melt}-T_{room}) \f$ \n
*/
/////////////////////////////////////////////////////////////////////////////

class JohnsonCookFlow : public FlowStressModel
{

public:
  // Create datatype for storing model parameters
  struct CMData
  {
    double A;
    double B;
    double C;
    double n;
    double m;
    double TRoom;
    double TMelt;
    double epdot_0;
  };

private:
  CMData d_CM;

  // Prevent copying of this class
  // copy constructor
  // JohnsonCookFlow(const JohnsonCookFlow &cm);
  JohnsonCookFlow& operator=(const JohnsonCookFlow& cm);

public:
  // constructors
  JohnsonCookFlow(ProblemSpecP& ps);
  JohnsonCookFlow(const JohnsonCookFlow* cm);

  // destructor
  ~JohnsonCookFlow() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  // Computes and requires for internal evolution variables
  void addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patches) override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, bool recurse,
                              bool SchedParent) override;

  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch, MPMLabel* lb) override;

  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  void initializeInternalVars(ParticleSubset* pset,
                              DataWarehouse* new_dw) override;

  void getInternalVars(ParticleSubset* pset, DataWarehouse* old_dw) override;

  void allocateAndPutInternalVars(ParticleSubset* pset,
                                  DataWarehouse* new_dw) override;

  void allocateAndPutRigid(ParticleSubset* pset,
                           DataWarehouse* new_dw) override;

  void updateElastic(const particleIndex idx) override;

  void updatePlastic(const particleIndex idx, const double& delGamma) override;

  ///////////////////////////////////////////////////////////////////////////
  /*! \brief  compute the flow stress */
  ///////////////////////////////////////////////////////////////////////////
  double computeFlowStress(const ModelStateBase* state, const double& delT,
                           const double& tolerance, const MPMMaterial* matl,
                           const particleIndex idx) override;

  //////////
  /*! \brief Calculate the plastic strain rate [epdot(tau,ep,T)] */
  //////////
  double computeEpdot(const ModelStateBase* state, const double& delT,
                      const double& tolerance, const MPMMaterial* matl,
                      const particleIndex idx) override;

  ///////////////////////////////////////////////////////////////////////////
  /*! Compute the elastic-plastic tangent modulus
  **WARNING** Assumes vonMises yield condition and the
  associated flow rule */
  ///////////////////////////////////////////////////////////////////////////
  void computeTangentModulus(const Matrix3& stress,
                             const ModelStateBase* state, const double& delT,
                             const MPMMaterial* matl, const particleIndex idx,
                             TangentModulusTensor& Ce,
                             TangentModulusTensor& Cep) override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to scalar and
    internal variables.

    \return Three derivatives in Vector deriv
    (deriv[0] = \f$d\sigma_Y/d\dot\epsilon\f$,
    deriv[1] = \f$d\sigma_Y/dT\f$,
    deriv[2] = \f$d\sigma_Y/d\epsilon\f$)
  */
  ///////////////////////////////////////////////////////////////////////////
  void evalDerivativeWRTScalarVars(const ModelStateBase* state,
                                   const particleIndex idx,
                                   Vector& derivs) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to plastic strain

    The Johnson-Cook yield stress is given by :
    \f[
    \sigma_Y(\epsilon_p) := D\left[A+B\epsilon_p^n\right]
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{d\epsilon_p} := nDB\epsilon_p^{n-1}
    \f]

    \return Derivative \f$ d\sigma_Y / d\epsilon_p\f$.

    \warning Expect error when \f$ \epsilon_p = 0 \f$.
  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTPlasticStrain(const ModelStateBase* state,
                                        const particleIndex idx) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to strain rate.

    The Johnson-Cook yield stress is given by :
    \f[
    \sigma_Y(\dot\epsilon_p) := E\left[1+C\ln\left(\frac{\dot\epsilon_p}
    {\dot\epsilon_{p0}}\right)\right]
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{d\dot\epsilon_p} := \frac{EC}{\dot\epsilon_p}
    \f]

    \return Derivative \f$ d\sigma_Y / d\dot\epsilon_p \f$.
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

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the melting temperature
  */
  ///////////////////////////////////////////////////////////////////////////
  double computeMeltingTemp(const ModelStateBase* state) override;

protected:
  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to temperature.

    The Johnson-Cook yield stress is given by :
    \f[
    \sigma_Y(T) := F\left[1-\left(\frac{T-T_r}{T_m-T_r}\right)^m\right]
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{dT} := -\frac{mF\left(\frac{T-T_r}{T_m-T_r}\right)^m}
    {T-T_r}
    \f]

    \return Derivative \f$ d\sigma_Y / dT \f$.

    \warning Expect error when \f$ T < T_{room} \f$.
  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTTemperature(const ModelStateBase* state,
                                      const particleIndex idx) const;
};

} // End namespace Uintah

#endif // __JOHNSONCOOK_FLOW_MODEL_H__
