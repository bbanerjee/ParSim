/*
 * The MIT License
 *
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

#ifndef __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__
#define __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class IntVar_MetalPlastic
  \brief The evolution of equivalent plastic strain and porosity
         for the hypoelasticplastic model
*/
////////////////////////////////////////////////////////////////////////////

class IntVar_MetalPlastic : public InternalVariableModel
{
public:
  const Uintah::VarLabel* pEqPlasticStrainLabel;
  const Uintah::VarLabel* pEqPlasticStrainLabel_preReloc;

  const Uintah::VarLabel* pPlasticPorosityLabel;
  const Uintah::VarLabel* pPlasticPorosityLabel_preReloc;

  /* constructors/destructor */
  IntVar_MetalPlastic(Uintah::ProblemSpecP& ps, ElasticModuliModel* elastic);
  IntVar_MetalPlastic(const IntVar_MetalPlastic* cm);
  ~IntVar_MetalPlastic() override;
  IntVar_MetalPlastic&
  operator=(const IntVar_MetalPlastic& cm) = delete;

  /* for restart */
  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  ParameterDict
  getParameters() const override
  {
    ParameterDict params;
    return params;
  }

  /* for relocation */
  void
  addParticleState(std::vector<const Uintah::VarLabel*>& from,
                   std::vector<const Uintah::VarLabel*>& to) override;

  /* Return the internal variable labels */
  std::vector<const Uintah::VarLabel*>
  getLabels() const override;

  /* initialize */
  void
  addInitialComputesAndRequires(Uintah::Task* task,
                                const Uintah::MPMMaterial* matl,
                                const Uintah::PatchSet* patches) override;
  void
  initializeInternalVariable(Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw) override;
  void
  initializeInternalVariable(const Uintah::Patch* patch,
                             const Uintah::MPMMaterial* matl,
                             Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw,
                             Uintah::MPMLabel* lb,
                             ParameterDict& params) override;

  /*! Compute the internal variable */
  void
  addComputesAndRequires(Uintah::Task* task,
                         const Uintah::MPMMaterial* matl,
                         const Uintah::PatchSet* patches) override;

  // Get the internal variables
  std::vector<Uintah::constParticleVariable<double>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw,
                       const double& dummy) override;

  std::vector<Uintah::constParticleVariable<Uintah::Matrix3>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw,
                       const Uintah::Matrix3& dummy);

  // Allocate and put the local particle internal variables
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 vectorParticleDoubleP& pVars) override;

  // Allocate and put the local <Matrix3> particle variables
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 vectorParticleMatrix3P& pVars) override;

  // Actually compute
  double
  computeInternalVariable(const Uintah::VarLabel* label,
                          const ModelStateBase* state) const override;
  double
  computeVolStrainDerivOfInternalVariable(const Uintah::VarLabel* label,
                                          const ModelStateBase*) const override;

  /* For material conversion */
  void
  allocateCMDataAddRequires(Uintah::Task* task,
                            const Uintah::MPMMaterial* matl,
                            const Uintah::PatchSet* patch,
                            Uintah::MPMLabel* lb) override;
  void
  allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                    Uintah::ParticleSubset* addset,
                    Uintah::ParticleLabelVariableMap* newState,
                    Uintah::ParticleSubset* delset,
                    Uintah::DataWarehouse* old_dw) override;

  /* For RigidMPM */
  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleLabelVariableMap& intvar) override;

private:
  /* Initialize local VarLabels */
  void
  initializeLocalMPMLabels();

  /**
   * Function: computeEqPlasticStrain
   */
  double
  computeEqPlasticStrain(const ModelStateBase* state) const;

  /**
   * Function: computePorosity
   */
  double
  computePorosity(const ModelStateBase* state) const;
};

} // End namespace Vaango

#endif // __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__
