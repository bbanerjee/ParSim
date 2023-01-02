/*
 * The MIT License
 *
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

#ifndef __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__
#define __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <Core/Util/Endian.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class IntVar_Metal
  \brief The evolution of equivalent plastic strain and porosity
         for the hypoelasticplastic model
*/
////////////////////////////////////////////////////////////////////////////
class IntVar_Metal : public InternalVariableModel
{
public:
  const Uintah::VarLabel* pIntVarLabel;
  const Uintah::VarLabel* pIntVarLabel_preReloc;

  /* constructors/destructor */
  IntVar_Metal(Uintah::ProblemSpecP& ps);
  IntVar_Metal(const IntVar_Metal* cm);
  IntVar_Metal(const IntVar_Metal& cm) = delete;
  ~IntVar_Metal() override;

  IntVar_Metal&
  operator=(const IntVar_Metal& cm) = delete;

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

  /* Get one (possibly composite) internal variable */
  template <typename T>
  void
  getInternalVariable(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::constParticleVariable<T>& intvar);

  /* Get multiple local <int/double/Vector/Matrix3> internal variables */
  template <typename T>
  std::vector<Uintah::constParticleVariable<T>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw);

  /* Allocate one (possibly composite) internal variable */
  template <typename T>
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 Uintah::ParticleVariable<T>& var);

  /* Allocate multiple local <int/double/Vector/Matrix3> internal variables */
  template <typename T>
  void
  allocateAndPutInternalVariable(
    Uintah::ParticleSubset* pset,
    Uintah::DataWarehouse* new_dw,
    std::vector<Uintah::ParticleVariable<T>>& pVars);

  // Actually compute
  template <typename T>
  void
  evolveInternalVariable(Uintah::particleIndex pidx,
                         const ModelStateBase* state,
                         Uintah::constParticleVariable<T>& var_old,
                         Uintah::ParticleVariable<T>& var);
  double
  computeInternalVariable(const std::string& label,
                          const ModelStateBase* state_old,
                          const ModelStateBase* state_cur) const override;

  // Hardening moduli
  template <typename T>
  void
  computeHardeningModulus(const ModelStateBase* state,
                          T& hardeningModulus) const;

  double
  computeVolStrainDerivOfInternalVariable(const std::string& label,
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
                      Uintah::constParticleVariableBase& intvar) override;
  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleLabelVariableMap& intvar) override;

private:
  // Create datatype for storing porosity parameters
  struct PoreNucleationModelParams
  {
    double phi_n;      /*< Volume fraction of void nucleating particles */
    double eps_mean_n; /*< Mean strain for nucleation */
    double eps_std_n;  /*< Standard deviation of strain for nucleation */
  };

  PoreNucleationModelParams d_poreNucleation;

  /* Initialize local VarLabels */
  void
  initializeLocalMPMLabels();

  void
  getPorosityModelParams(Uintah::ProblemSpecP& ps);
  void
  copyPorosityModelParams(const IntVar_Metal* cm);

  double
  computeEqPlasticStrain(double eqPlasticStrain_old,
                         const ModelStateBase* state) const;

  double
  computePlasticPorosity(double plasticPorosity_old,
                         const ModelStateBase* state) const;

  double
  eqPlasticStrainHardeningModulus(const ModelStateBase* state) const;

  double
  plasticPorosityHardeningModulus(const ModelStateBase* state) const;

  /* Calculate hardening modulus due to void growth */
  double
  poreGrowthHardeningModulus(const ModelStateBase* state) const;

  /* Calculate hardening modulus due to void nucleation */
  double
  poreNucleationHardeningModulus(const ModelStateBase* state) const;

  /* Calculate the void nucleation factor */
  double
  voidNucleationFactor(double eqPlasticStrain) const;
};

} // End namespace Vaango

#endif // __ELASTIC_PLASTIC_HP_INT_VAR_MODEL_H__
