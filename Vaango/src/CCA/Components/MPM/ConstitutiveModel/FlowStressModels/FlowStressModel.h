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

#ifndef __PLASTIC_FLOW_STRESS_MODEL_H__
#define __PLASTIC_FLOW_STRESS_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/TangentModulusTensor.h>
#include <vector>

namespace Uintah {

///////////////////////////////////////////////////////////////////////////
/*!
  \class  FlowStressModel
  \brief  Abstract Base class for flow models (calculate yield stress)
  \author Biswajit Banerjee, \n
          C-SAFE and Department of Mechanical Engineering, \n
          University of Utah,\n
  \warn   Assumes vonMises yield condition and the associated flow rule for
          all cases other than Gurson plasticity.
*/
///////////////////////////////////////////////////////////////////////////

using Vaango::ModelStateBase;

class FlowStressModel
{

private:
public:
  FlowStressModel();
  virtual ~FlowStressModel();

  virtual void outputProblemSpec(ProblemSpecP& ps) = 0;

  // Computes and requires for internal evolution variables
  virtual void
  addInitialComputesAndRequires([[maybe_unused]] Task* task,
                                [[maybe_unused]] const MPMMaterial* matl,
                                [[maybe_unused]] const PatchSet* patches) {};

  virtual void
  addComputesAndRequires([[maybe_unused]] Task* task,
                         [[maybe_unused]] const MPMMaterial* matl,
                         [[maybe_unused]] const PatchSet* patches) {};

  virtual void
  addComputesAndRequires([[maybe_unused]] Task* task,
                         [[maybe_unused]] const MPMMaterial* matl,
                         [[maybe_unused]] const PatchSet* patches,
                         [[maybe_unused]] bool recurse,
                         [[maybe_unused]] bool SchedParent) {};

  virtual void
  allocateCMDataAddRequires([[maybe_unused]] Task* task,
                            [[maybe_unused]] const MPMMaterial* matl,
                            [[maybe_unused]] const PatchSet* patch,
                            [[maybe_unused]] MPMLabel* lb) {};

  virtual void
  allocateCMDataAdd([[maybe_unused]] DataWarehouse* new_dw,
                    [[maybe_unused]] ParticleSubset* addset,
                    [[maybe_unused]] ParticleLabelVariableMap* newState,
                    [[maybe_unused]] ParticleSubset* delset,
                    [[maybe_unused]] DataWarehouse* old_dw) {};

  virtual void
  addParticleState([[maybe_unused]] std::vector<const VarLabel*>& from,
                   [[maybe_unused]] std::vector<const VarLabel*>& to) {};

  virtual void
  initializeInternalVars([[maybe_unused]] ParticleSubset* pset,
                         [[maybe_unused]] DataWarehouse* new_dw) {};

  virtual void
  getInternalVars([[maybe_unused]] ParticleSubset* pset,
                  [[maybe_unused]] DataWarehouse* old_dw) {};

  virtual void
  allocateAndPutInternalVars([[maybe_unused]] ParticleSubset* pset,
                             [[maybe_unused]] DataWarehouse* new_dw) {};

  virtual void
  allocateAndPutRigid([[maybe_unused]] ParticleSubset* pset,
                      [[maybe_unused]] DataWarehouse* new_dw) {};

  virtual void
  updateElastic([[maybe_unused]] const particleIndex idx) {};

  virtual void
  updatePlastic([[maybe_unused]] const particleIndex idx,
                [[maybe_unused]] const double& delGamma) {};

  //////////
  /*! \brief Calculate the flow stress */
  //////////
  virtual double computeFlowStress(const ModelStateBase* state,
                                   const double& delT, const double& tolerance,
                                   const MPMMaterial* matl,
                                   const particleIndex idx) = 0;

  //////////
  /*! \brief Calculate the plastic strain rate [epdot(tau,ep,T)] */
  //////////
  virtual double computeEpdot(const ModelStateBase* state, const double& delT,
                              const double& tolerance, const MPMMaterial* matl,
                              const particleIndex idx) = 0;

  /*! Compute the elastic-plastic tangent modulus
    This is given by
    \f[
    C_{ep} = C_e - \frac{(C_e:r) (x) (f_s:C_e)}
    {-f_q.h + f_s:C_e:r}
    \f]
    where \n
    \f$ C_{ep} \f$ is the continuum elasto-plastic tangent modulus \n
    \f$ C_{e} \f$ is the continuum elastic tangent modulus \n
    \f$ r \f$ is the plastic flow direction \f$ d\phi/d\sigma = r \f$\n
    \f$ h \f$ gives the evolution of \f$ q \f$ \n
    \f$ f_s = \partial f /\partial \sigma \f$ \n
    \f$ f_q = \partial f /\partial q \f$
  */
  virtual void
  computeTangentModulus([[maybe_unused]] const Matrix3& stress,
                        [[maybe_unused]] const ModelStateBase* state,
                        [[maybe_unused]] const double& delT,
                        [[maybe_unused]] const MPMMaterial* matl,
                        [[maybe_unused]] const particleIndex idx,
                        [[maybe_unused]] TangentModulusTensor& Ce,
                        [[maybe_unused]] TangentModulusTensor& Cep) {};

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to scalar and
      internal variables.

    \return Three derivatives in Vector derivs
      (derivs[0] = \f$d\sigma_Y/d\dot\epsilon\f$,
       derivs[1] = \f$d\sigma_Y/dT\f$,
       derivs[2] = \f$d\sigma_Y/d(int. var.)\f$)
  */
  ///////////////////////////////////////////////////////////////////////////
  virtual void evalDerivativeWRTScalarVars([[maybe_unused]] const ModelStateBase* state,
                                           [[maybe_unused]] const particleIndex idx,
                                           [[maybe_unused]] Vector& derivs) const {};

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to plastic
      strain.

    \return \f$d\sigma_Y/d\epsilon_p\f$
  */
  ///////////////////////////////////////////////////////////////////////////
  virtual double evalDerivativeWRTPlasticStrain(const ModelStateBase* state,
                                                const particleIndex idx) const = 0;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Evaluate derivative of flow stress with respect to plastic
      strain rate.

    \return \f$d\sigma_Y/d\dot{\epsilon_p}\f$
  */
  ///////////////////////////////////////////////////////////////////////////
  virtual double evalDerivativeWRTStrainRate(const ModelStateBase* state,
                                             const particleIndex idx) const  = 0;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the shear modulus.
  */
  ///////////////////////////////////////////////////////////////////////////
  virtual double computeShearModulus(const ModelStateBase* state) = 0;

  ///////////////////////////////////////////////////////////////////////////
  /*!
    \brief Compute the melting temperature
  */
  ///////////////////////////////////////////////////////////////////////////
  virtual double computeMeltingTemp(const ModelStateBase* state) = 0;
};
} // End namespace Uintah

#endif // __PLASTIC_FLOW_STRESS_MODEL_H__
