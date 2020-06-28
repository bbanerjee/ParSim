/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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

#ifndef __ISOHARDENING_FLOW_MODEL_H__
#define __ISOHARDENING_FLOW_MODEL_H__

#include "FlowStressModel.h"
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

/////////////////////////////////////////////////////////////////////////////
/*!
  \class IsoHardeningFlow
  \brief Isotropic Hardening flow model.
  (Simo and Hughes, 1998, Computational Inelasticity, p. 319)
  \author Biswajit Banerjee,
  \author Department of Mechanical Engineering,
  \author University of Utah,

  The flow rule is given by
  \f[
  f(\sigma) = K \alpha + \sigma_0
  \f]

  where \f$f(\sigma)\f$ = flow stress \n
  \f$K\f$ = hardening modulus \n
  \f$\alpha\f$ = evolution parameter for hardening law \n
  \f$\sigma_0\f$ = initial yield stress \n

  Evolution of alpha is given by
  \f[
  d\alpha/dt = \sqrt{2/3}*\gamma
  \f]
  where \f$\gamma\f$ = consistency parameter
*/
/////////////////////////////////////////////////////////////////////////////

class IsoHardeningFlow : public FlowStressModel
{

public:
  struct CMData
  {
    double K;
    double sigma_0;
  };

  constParticleVariable<double> pAlpha;
  ParticleVariable<double> pAlpha_new;

  const VarLabel* pAlphaLabel;          // For Isotropic Hardening Plasticity
  const VarLabel* pAlphaLabel_preReloc; // For Isotropic Hardening Plasticity

  // constructors
  IsoHardeningFlow(ProblemSpecP& ps);
  IsoHardeningFlow(const IsoHardeningFlow* cm);
  IsoHardeningFlow& operator=(const IsoHardeningFlow& cm) = delete;

  // destructor
  ~IsoHardeningFlow() override;

  void outputProblemSpec(ProblemSpecP& ps) override;

  // Computes and requires for internal evolution variables
  // Only one internal variable for Isotropic-Hardening :: plastic strain
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
  /*! compute the flow stress */
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
    deriv[2] = \f$d\sigma_Y/d\alpha\f$)
  */
  ///////////////////////////////////////////////////////////////////////////
  void evalDerivativeWRTScalarVars(const ModelStateBase* state,
                                   const particleIndex idx,
                                   Vector& derivs) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to plastic
    strain.

    The Isotropic-Hardening yield stress is given by :
    \f[
    \sigma_Y(\alpha) := \sigma_0 + K\alpha
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{d\epsilon_p} := K
    \f]

    \return Derivative \f$ d\sigma_Y / d\alpha\f$.

    \warning Not implemented yet.

  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTPlasticStrain(const ModelStateBase* state,
                                        const particleIndex idx) const override;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to strain rate.

    The Isotropic-Hardening yield stress is given by :
    \f[
    \sigma_Y(\dot\epsilon_p) := C
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{d\dot\epsilon_p} := 0
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

    The Isotropic-Hardening yield stress is given by :
    \f[
    \sigma_Y(T) := C
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{dT} := 0
    \f]

    \return Derivative \f$ d\sigma_Y / dT \f$.
  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTTemperature(const ModelStateBase* state,
                                      const particleIndex idx) const;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to alpha

    The Isotropic-Hardening yield stress is given by :
    \f[
    \sigma_Y(\alpha) := \sigma_0 + K\alpha
    \f]

    The derivative is given by
    \f[
    \frac{d\sigma_Y}{d\alpha} := K
    \f]

    \return Derivative \f$ d\sigma_Y / d\alpha\f$.

  */
  ///////////////////////////////////////////////////////////////////////////
  double evalDerivativeWRTAlpha(const ModelStateBase* state,
                                const particleIndex idx) const;

private:

  CMData d_CM;

};

} // End namespace Uintah

#endif // __ISOHARDENING_FLOW_MODEL_H__
