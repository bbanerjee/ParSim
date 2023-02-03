/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef __INTERNAL_VARIABLE_MODEL_H__
#define __INTERNAL_VARIABLE_MODEL_H__

#include <CCA/Components/MPM/Core/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Math/Matrix3.h>
#include <map>
#include <vector>

namespace Vaango {

using ParameterDict = std::map<std::string, double>;

class ElasticModuliModel;
class ShearModulusModel;
class MPMEquationOfState;

///////////////////////////////////////////////////////////////////////////
/*!
  \class  InternalVariableModel
  \brief  Abstract base class for the evolution of internal variables
          for plasticity models
*/
///////////////////////////////////////////////////////////////////////////

using constParticleDouble     = Uintah::constParticleVariable<double>;
using constParticleMatrix3    = Uintah::constParticleVariable<Uintah::Matrix3>;
using ParticleDouble          = Uintah::ParticleVariable<double>;
using ParticleMatrix3         = Uintah::ParticleVariable<Uintah::Matrix3>;
using constParticleDoubleVec  = std::vector<constParticleDouble>;
using constParticleMatrix3Vec = std::vector<constParticleMatrix3>;
using ParticleDoublePVec      = std::vector<ParticleDouble*>;
using ParticleMatrix3PVec     = std::vector<ParticleMatrix3*>;
using ParamMap                = std::map<std::string, double>;

class InternalVariableModel
{

public:
  InternalVariableModel();
  virtual ~InternalVariableModel();

  virtual void
  outputProblemSpec(Uintah::ProblemSpecP& ps) = 0;

  virtual void
  addParticleState(std::vector<const Uintah::VarLabel*>& from,
                   std::vector<const Uintah::VarLabel*>& to) = 0;

  /* Get the internal variable labels */
  virtual std::vector<const Uintah::VarLabel*>
  getLabels() const = 0;

  /* Get the model parameters */
  virtual ParamMap
  getParameters() const = 0;

  /* Initialialize
   *   Sometimes extra parameters from the YieldCondition/Modulus models may
   *   need to be passed
   */
  virtual void
  addInitialComputesAndRequires(Uintah::Task* task,
                                const Uintah::MPMMaterial* matl,
                                const Uintah::PatchSet* patches) = 0;
  virtual void
  initializeInternalVariable(Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw) = 0;
  virtual void
  initializeInternalVariable(const Uintah::Patch* patch,
                             const Uintah::MPMMaterial* matl,
                             Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw,
                             Uintah::MPMLabel* lb,
                             ParamMap& params) = 0;

  /* Compute */
  virtual void
  addComputesAndRequires(Uintah::Task* task,
                         const Uintah::MPMMaterial* matl,
                         const Uintah::PatchSet* patches) = 0;

  /*! Update the internal variable */
  template <typename T>
  void
  evolveInternalVariable(Uintah::particleIndex pidx,
                         const ModelStateBase* state,
                         Uintah::constParticleVariable<T>& var_old,
                         Uintah::ParticleVariable<T>& var_new);

  /* Compute the internal variable */
  virtual double
  computeInternalVariable(const std::string& label,
                          const ModelStateBase* state_old,
                          const ModelStateBase* state_cur) const = 0;

  /* Compute derivative of internal variable with respect to volumetric
     elastic strain */
  virtual double
  computeVolStrainDerivOfInternalVariable(
    const std::string& label,
    const ModelStateBase* state) const = 0;

  /* For material conversion */
  virtual void
  allocateCMDataAddRequires(Uintah::Task* task,
                            const Uintah::MPMMaterial* matl,
                            const Uintah::PatchSet* patch,
                            Uintah::MPMLabel* lb) = 0;
  virtual void
  allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                    Uintah::ParticleSubset* addset,
                    Uintah::ParticleLabelVariableMap* newstate,
                    Uintah::ParticleSubset* delset,
                    Uintah::DataWarehouse* old_dw) = 0;

  /* For RigidMPM */
  virtual void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleVariableBase& intvar) = 0;
  virtual void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleLabelVariableMap& intvars) = 0;

protected:
  ElasticModuliModel* d_elastic;
  ShearModulusModel* d_shear;
  MPMEquationOfState* d_eos;
};
} // End namespace Uintah

#endif // __INTERNAL_VARIABLE_MODEL_H__
