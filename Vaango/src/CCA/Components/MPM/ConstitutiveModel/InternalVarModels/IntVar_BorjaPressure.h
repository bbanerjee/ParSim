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

#ifndef __BORJA_PRESSURE_INT_VAR_MODEL_H__
#define __BORJA_PRESSURE_INT_VAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_CamClay.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class IntVar_BorjaPressure
  \brief The evolution of the consolidation pressure internal variable in the
         Borja model

  Reference:Borja, R.I. and Tamagnini, C.(1998) Cam-Clay plasticity Part III:
  Extension of the infinitesimal model to include finite strains,
  Computer Methods in Applied Mechanics and Engineering, 155 (1-2),
  pp. 73-95.

  The consolidation presure (p_c) is defined by the rate equation

         1/p_c dp_c/dt = 1/(lambdatilde - kappatilde) depsp_v/dt

         where lambdatilde = material constant
               kappatilde = material constant
               epsp_v = volumetric plastic strain

  The incremental update of the consolidation pressure is given by

     (p_c)_{n+1} = (p_c)_n exp[((epse_v)_trial - (epse_v)_{n+1})/(lambdabar -
  kappabar)]
*/
////////////////////////////////////////////////////////////////////////////

class IntVar_BorjaPressure : public InternalVariableModel
{

public:
  // Internal variables
  const Uintah::VarLabel* pPcLabel;
  const Uintah::VarLabel* pPcLabel_preReloc;

  // Return the internal variable labels
  std::vector<const Uintah::VarLabel*>
  getLabels() const override
  {
    std::vector<const Uintah::VarLabel*> labels;
    labels.push_back(pPcLabel); // Preconsolidation pressure
    labels.push_back(pPcLabel_preReloc);

    return labels;
  }

private:
  // Model parameters
  double d_pc0;
  double d_lambdatilde;
  double d_kappatilde;

  // Prevent copying of this class
  // copy constructor
  // IntVar_BorjaPressure(const IntVar_BorjaPressure &cm);
  IntVar_BorjaPressure&
  operator=(const IntVar_BorjaPressure& cm);

public:
  // constructors
  IntVar_BorjaPressure(Uintah::ProblemSpecP& ps, ShearModulusModel* shear);
  IntVar_BorjaPressure(const IntVar_BorjaPressure* cm);

  // destructor
  ~IntVar_BorjaPressure() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double>
  getParameters() const override
  {
    std::map<std::string, double> params;
    params["pc0"]         = d_pc0;
    params["lambdatilde"] = d_lambdatilde;
    params["kappatilde"]  = d_kappatilde;
    return params;
  }

  void
  addParticleState(std::vector<const Uintah::VarLabel*>& from,
                   std::vector<const Uintah::VarLabel*>& to) override;

  // Computes and requires for internal evolution variables
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
                             ParamMap& params) override {}

  void
  addComputesAndRequires(Uintah::Task* task,
                         const Uintah::MPMMaterial* matl,
                         const Uintah::PatchSet* patches) override;


  /* Get one (possibly composite) internal variable */
  template<typename T>
  void
  getInternalVariable(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* old_dw,
                      Uintah::constParticleVariable<T>& intvar);

  /* Get multiple local <int/double/Vector/Matrix3> internal variables */
  template<typename T>
  std::vector<Uintah::constParticleVariable<T>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw);

  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 Uintah::ParticleVariableBase& intvar) override;

  /* Allocate and put the local <double> particle variables */
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 ParticleDoublePVec& pVars) override {}

  /* Allocate and put the local <Matrix3> particle variables */
  void
  allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                 Uintah::DataWarehouse* new_dw,
                                 ParticleMatrix3PVec& pVars) override {}

  /*! \brief Compute the internal variable */
  template <typename T>
  void
  evolveInternalVariable(Uintah::particleIndex pidx,
                         const ModelStateBase* state,
                         Uintah::constParticleVariable<T>& var_old,
                         Uintah::ParticleVariable<T>& var);
  double
  computeInternalVariable(const std::string& label,
                          const ModelStateBase* state) const override;

  // Compute derivative of internal variable with respect to volumetric
  // elastic strain
  double
  computeVolStrainDerivOfInternalVariable(
    const std::string& label,
    const ModelStateBase* state) const override;

  void
  allocateCMDataAddRequires(Uintah::Task* task,
                            const Uintah::MPMMaterial* matl,
                            const Uintah::PatchSet* patch,
                            Uintah::MPMLabel* lb) override;
  void
  allocateCMDataAdd(
    Uintah::DataWarehouse* new_dw,
    Uintah::ParticleSubset* addset,
    std::map<const Uintah::VarLabel*, Uintah::ParticleVariableBase*>* newState,
    Uintah::ParticleSubset* delset,
    Uintah::DataWarehouse* old_dw) override;

  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleVariableBase& intvar) override;

  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleLabelVariableMap& intvars) override {}
};

} // End namespace Uintah

#endif // __BORJA_PRESSURE_INT_VAR_MODEL_H__
