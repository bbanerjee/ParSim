/*
 * The MIT License
 *
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

#ifndef __ARENA_INT_VAR_MODEL_H__
#define __ARENA_INT_VAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVariableModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Vaango {

////////////////////////////////////////////////////////////////////////////
/*!
  \class IntVar_Arena
  \brief The evolution of the kappa, X, porosity, and saturation internal
  variables
         in the partially saturated Arenisca model
*/
////////////////////////////////////////////////////////////////////////////

class IntVar_Arena : public InternalVariableModel
{
public:
  // Internal variables
  const Uintah::VarLabel* pKappaLabel; // Branch point
  const Uintah::VarLabel* pKappaLabel_preReloc;

  const Uintah::VarLabel* pCapXLabel; // Hydrostatic strength
  const Uintah::VarLabel* pCapXLabel_preReloc;

  const Uintah::VarLabel* pPlasticStrainLabel; // Plastic Strain
  const Uintah::VarLabel* pPlasticStrainLabel_preReloc;

  const Uintah::VarLabel* pPlasticVolStrainLabel; // Plastic Volumetric Strain
  const Uintah::VarLabel* pPlasticVolStrainLabel_preReloc;

  const Uintah::VarLabel* pP3Label; // Evolution of parameter P3
  const Uintah::VarLabel* pP3Label_preReloc;

  // constructors
  IntVar_Arena(Uintah::ProblemSpecP& ps, ElasticModuliModel* elastic);
  IntVar_Arena(const IntVar_Arena* cm);

  IntVar_Arena&
  operator=(const IntVar_Arena& cm) = delete;

  // destructor
  ~IntVar_Arena() override;

  void
  outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  void
  addParticleState(std::vector<const Uintah::VarLabel*>& from,
                   std::vector<const Uintah::VarLabel*>& to) override;

  /*! Get parameters */
  ParameterDict
  getParameters() const override
  {
    ParameterDict params;
    params["p0"] = d_crushParam.p0;
    params["p1"] = d_crushParam.p1;
    params["p2"] = d_crushParam.p2;
    params["p3"] = d_crushParam.p3;
    return params;
  }

  // Return the internal variable labels
  std::vector<const Uintah::VarLabel*>
  getLabels() const override
  {
    std::vector<const Uintah::VarLabel*> labels;
    labels.push_back(pKappaLabel); // Branch point
    labels.push_back(pKappaLabel_preReloc);

    labels.push_back(pCapXLabel); // Hydrostatic strength
    labels.push_back(pCapXLabel_preReloc);

    labels.push_back(pPlasticStrainLabel); // Plastic Strain
    labels.push_back(pPlasticStrainLabel_preReloc);

    labels.push_back(pPlasticVolStrainLabel); // Plastic Volumetric Strain
    labels.push_back(pPlasticVolStrainLabel_preReloc);

    labels.push_back(pP3Label); // Evolution of parameter P3
    labels.push_back(pP3Label_preReloc);

    return labels;
  }

  // Computes and requires for internal evolution variables
  void
  addInitialComputesAndRequires(Uintah::Task* task,
                                const Uintah::MPMMaterial* matl,
                                const Uintah::PatchSet* patches) override;
  void
  initializeInternalVariable(Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw) override
  {
  }
  void
  initializeInternalVariable(const Uintah::Patch* patch,
                             const Uintah::MPMMaterial* matl,
                             Uintah::ParticleSubset* pset,
                             Uintah::DataWarehouse* new_dw,
                             Uintah::MPMLabel* lb,
                             ParameterDict& params) override;

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
                                 Uintah::ParticleVariable<T>& intvar);

  /* Allocate multiple local <int/double/Vector/Matrix3> internal variables */
  template <typename T>
  void
  allocateAndPutInternalVariable(
    Uintah::ParticleSubset* pset,
    Uintah::DataWarehouse* new_dw,
    std::vector<Uintah::ParticleVariable<T>>& pVars);

  /*! \brief Compute the internal variable */
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

  // Compute derivative of internal variable with respect to volumetric
  // elastic strain
  double
  computeVolStrainDerivOfInternalVariable(const std::string& label,
                                          const ModelStateBase*) const override
  {
    return 0.0;
  }

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
  virtual void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleVariableBase& intvar) override
  {
  }
  void
  allocateAndPutRigid(Uintah::ParticleSubset* pset,
                      Uintah::DataWarehouse* new_dw,
                      Uintah::constParticleLabelVariableMap& intvar) override;

public:
  /**
   * Function: computeDrainedHydrostaticStrength
   *
   * Purpose:
   *   Compute the drained hydrostatic strength (X) using the crush curve
   *   and the initial porosity
   *
   * Inputs:
   *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
   *   phi0     = initial porosity
   *
   * Returns:
   *   Xbar     = hydrostatic compressive strength
   */
  double
  computeDrainedHydrostaticStrength(const double& ep_v_bar,
                                    const double& phi0) const;

  /**
   * Function: computeP3
   *
   * Purpose:
   *   Compute the crush curve parameter P3 from the initial porosity
   *
   * Inputs:
   *   phi0 = initial porosity
   *
   * Returns:
   *   p3 = crush curve parameter
   */
  inline double
  computeP3(const double& phi0) const
  {
    double p3 = -std::log(1.0 - phi0);
    return p3;
  };

  /**
   * Function: computePorosity
   *
   * Purpose:
   *   Compute the porosity from crush curve parameter P3
   *
   * Inputs:
   *   ep_v = volumetric plastic strain
   *   p3 = crush
   *
   * Returns:
   *   porosity = porosity from crush curve parameter
   */
  inline double
  computePorosity(const double& ep_v, const double& p3) const
  {
    double porosity = 1.0 - std::exp(-p3 + ep_v);
    return porosity;
  }

  /**
   * Function: computeElasticVolStrainAtYield
   *
   * Purpose:
   *   Compute the elastic volumetric strain at yield
   *
   * Inputs:
   *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
   *   phi0     = initial porosity
   *
   * Returns:
   *   ev_e_yield = elastic volumetric strain at yield
   */
  double
  computeElasticVolStrainAtYield(const double& ep_v_bar,
                                 const double& phi0) const;

  /**
   * Function: computePartSatHydrostaticStrength
   *
   * Purpose:
   *   Compute the partially saturated hydrostatic strength (X_sat)
   *
   * Inputs:
   *   I1_eff_bar   = -tr(sigma) = isotropic part of effective stress tensor
   *   pw_bar   = pore pressure
   *   ep_v_bar = -tr(ep) = volumetric part of plastic strain tensor
   *   phi      = porosity
   *   Sw       = saturation
   *   phi0     = initial porosity
   *
   * Returns:
   *   Xbar_sat = hydrostatic compressive strength
   */
  double
  computePartSatHydrostaticStrength(const double& I1_bar,
                                    const double& pw_bar,
                                    const double& ep_v_bar,
                                    const double& phi,
                                    const double& Sw,
                                    const double& phi0) const;

private:
  // Crush Curve Model parameters
  struct CrushParameters
  {
    double p0;
    double p1;
    double p2;
    double p3;
  };

  CrushParameters d_crushParam;
  bool d_use_disaggregation_algorithm;

  // Initialize local VarLabels
  void
  initializeLocalMPMLabels();
};

} // End namespace Uintah

#endif // __ARENA_INT_VAR_MODEL_H__
