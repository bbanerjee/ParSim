/*
 * The MIT License
 *
 * Copyright (c) 2015-2018 Parresia Research Limited, New Zealand
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

#ifndef __INTERNAL_VAR_MODEL_TABULAR_PLASTICITY_CAP_H__
#define __INTERNAL_VAR_MODEL_TABULAR_PLASTICITY_CAP_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_TabularCap.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

  using Polyline = std::vector<Uintah::Point>;

////////////////////////////////////////////////////////////////////////////
/*!
  \class InternalVar_TabularCap
  \brief The evolution of the hydrostatic strength (X) internal variable
         in the tabular plasticity model with cap
*/
////////////////////////////////////////////////////////////////////////////

class InternalVar_TabularCap : public InternalVariableModel
{

public:

  const Uintah::VarLabel* pCapXLabel; // Hydrostatic strength
  const Uintah::VarLabel* pCapXLabel_preReloc;

  // constructors
  explicit InternalVar_TabularCap(Uintah::ProblemSpecP& ps);
  InternalVar_TabularCap(const InternalVar_TabularCap* cm);

  // destructor
  ~InternalVar_TabularCap() override;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  ParameterDict getParameters() const override
  {
    ParameterDict params;
    return params;
  }

  // Computes and requires for internal evolution variables
  void addInitialComputesAndRequires(Uintah::Task* task,
                                     const Uintah::MPMMaterial* matl,
                                     const Uintah::PatchSet* patches) override;

  void initializeInternalVariable(Uintah::ParticleSubset* pset,
                                  Uintah::DataWarehouse* new_dw) override;

  void initializeInternalVariable(const Uintah::Patch* patch,
                                  const Uintah::MPMMaterial* matl,
                                  Uintah::ParticleSubset* pset,
                                  Uintah::DataWarehouse* new_dw,
                                  Uintah::MPMLabel* lb,
                                  ParameterDict& params) override
  {
  }

  void addComputesAndRequires(Uintah::Task* task,
                              const Uintah::MPMMaterial* matl,
                              const Uintah::PatchSet* patches) override;

  void allocateCMDataAddRequires(Uintah::Task* task,
                                 const Uintah::MPMMaterial* matl,
                                 const Uintah::PatchSet* patch,
                                 Uintah::MPMLabel* lb) override;

  void allocateCMDataAdd(Uintah::DataWarehouse* new_dw,
                         Uintah::ParticleSubset* addset,
                         Uintah::ParticleLabelVariableMap* newState,
                         Uintah::ParticleSubset* delset,
                         Uintah::DataWarehouse* old_dw) override;

  void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                        std::vector<const Uintah::VarLabel*>& to) override;

  void allocateAndPutRigid(
    Uintah::ParticleSubset* pset, Uintah::DataWarehouse* new_dw,
    Uintah::constParticleLabelVariableMap& intvar) override;

  ///////////////////////////////////////////////////////////////////////////
  /*! \brief Compute the internal variable */
  double computeInternalVariable(const ModelStateBase* state) const override;

  ///////////////////////////////////////////////////////////////////////////
  // Compute derivative of internal variable with respect to volumetric
  // plastic strain
  double computeVolStrainDerivOfInternalVariable(const ModelStateBase*) const 
    override;

  void getInternalVariable(Uintah::ParticleSubset* pset,
                           Uintah::DataWarehouse* old_dw,
                           Uintah::constParticleVariableBase& pCapX) override
  {
    old_dw->get(pCapX, pCapXLabel, pset);
  }

  void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset, 
                                      Uintah::DataWarehouse* new_dw,
                                      Uintah::ParticleVariableBase& pCapX_new) override
  {
    new_dw->allocateAndPut(pCapX_new, pCapXLabel_preReloc, pset);
  }

  // Get the internal variables
  std::vector<Uintah::constParticleVariable<double>> getInternalVariables(
    Uintah::ParticleSubset* pset, Uintah::DataWarehouse* old_dw,
    const double& dummy) override
  {
    std::vector<Uintah::constParticleVariable<double>> pIntVars;
    return pIntVars;
  }

  std::vector<Uintah::constParticleVariable<Uintah::Matrix3>>
  getInternalVariables(Uintah::ParticleSubset* pset,
                       Uintah::DataWarehouse* old_dw,
                       const Uintah::Matrix3& dummy) override
  {
    std::vector<Uintah::constParticleVariable<Uintah::Matrix3>> pIntVars;
    return pIntVars;
  }

  // Allocate and put the local particle internal variables
  void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                      Uintah::DataWarehouse* new_dw,
                                      vectorParticleDoubleP& pVars) override
  {
  }

  // Allocate and put the local <Matrix3> particle variables
  void allocateAndPutInternalVariable(Uintah::ParticleSubset* pset,
                                      Uintah::DataWarehouse* new_dw,
                                      vectorParticleMatrix3P& pVars) override
  {
  }

  // Return the internal variable labels
  std::vector<const Uintah::VarLabel*> getLabels() const override
  {
    std::vector<const Uintah::VarLabel*> labels;

    labels.push_back(pCapXLabel); // Hydrostatic strength
    labels.push_back(pCapXLabel_preReloc);

    return labels;
  }

private:

  /**
   *  These are the parameters that are read from the input file
   */
  struct HydrostaticStrengthFunction
  {
    TabularData table;
    HydrostaticStrengthFunction() = default;
    HydrostaticStrengthFunction(Uintah::ProblemSpecP& ps)
      : table(ps)
    {
      table.setup();
    }
    HydrostaticStrengthFunction(const HydrostaticStrengthFunction& hsf)
    {
      table = hsf.table;
    }
    HydrostaticStrengthFunction& operator=(const HydrostaticStrengthFunction& hsf)
    {
      if (this != &hsf) {
        table = hsf.table;
      }
      return *this;
    }
  };

  HydrostaticStrengthFunction d_capX_fn;

  // Prevent copying of this class
  // copy constructor
  InternalVar_TabularCap(const InternalVar_TabularCap &cm) = delete;
  InternalVar_TabularCap& operator=(const InternalVar_TabularCap& cm);

  // Initialize local VarLabels
  void initializeLocalMPMLabels()
  {
    pCapXLabel = Uintah::VarLabel::create(
      "p.capX", Uintah::ParticleVariable<double>::getTypeDescription());
    pCapXLabel_preReloc = Uintah::VarLabel::create(
      "p.capX+", Uintah::ParticleVariable<double>::getTypeDescription());
  }

  /**
   * Function: computeDrainedHydrostaticStrength
   *
   * Purpose:
   *   Compute the drained hydrostatic strength (X) using the crush curve
   *   and the initial porosity
   *
   * Inputs:
   *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
   *
   * Returns:
   *   Xbar     = hydrostatic compressive strength
   */
  double computeDrainedHydrostaticStrength(const double& ep_v_bar) const;

};

} // End namespace Vaango

#endif // __INTERNAL_VAR_MODEL_TABULAR_PLASTICITY_CAP_H__