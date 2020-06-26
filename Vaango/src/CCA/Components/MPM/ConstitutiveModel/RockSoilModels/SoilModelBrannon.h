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

#ifndef __SOIL_MODEL_BRANNON_H__
#define __SOIL_MODEL_BRANNON_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_SoilModelBrannonKappa.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cmath>

namespace Uintah {

class MPMLabel;
class MPMFlags;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class SoilModelBrannon
  \brief Nonlinear Drucker-Prager model with cap
*/
/////////////////////////////////////////////////////////////////////////////

class SoilModelBrannon : public ConstitutiveModel
{

public:
  const VarLabel* pPlasticStrainLabel;
  const VarLabel* pPlasticStrainVolLabel;
  const VarLabel* pElasticStrainVolLabel;
  const VarLabel* pBackStressLabel;
  const VarLabel* pBackStressIsoLabel;
  const VarLabel* pKappaStateLabel;
  const VarLabel* pLocalizedLabel;

  const VarLabel* pPlasticStrainLabel_preReloc;
  const VarLabel* pPlasticStrainVolLabel_preReloc;
  const VarLabel* pElasticStrainVolLabel_preReloc;
  const VarLabel* pBackStressLabel_preReloc;
  const VarLabel* pBackStressIsoLabel_preReloc;
  const VarLabel* pKappaStateLabel_preReloc;
  const VarLabel* pLocalizedLabel_preReloc;

  // Create datatype for storing model parameters
  struct CMData
  {
    double fSlope;
    double fSlope_p;
    double hardening_modulus;
    double cap_ratio;
    double p0_crush_curve;
    double p1_crush_curve;
    double p3_crush_curve;
    double p4_fluid_effect;
    double kinematic_hardening_constant;
    double fluid_B0;
    double fluid_pressure_initial;
    double subcycling_characteristic_number;
    double peakI1;
    double B0;
    double G0;
  };

  // constructor
  SoilModelBrannon(ProblemSpecP& ps, MPMFlags* flag);
  SoilModelBrannon(const SoilModelBrannon* cm);
  SoilModelBrannon& operator=(const SoilModelBrannon& cm) = delete;

  // destructor
  ~SoilModelBrannon() override;

  ModelType modelType() const override
  {
    return ModelType::RATE_FORM;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  // clone

  SoilModelBrannon* clone() override;

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  // compute stress at each particle in the patch
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  void computeInvariants(Matrix3& stress, Matrix3& S, double& I1, double& J2);

  void computeInvariants(const Matrix3& stress, Matrix3& S, double& I1,
                         double& J2);

  double YieldFunction(const Matrix3& stress, const double& fSlope,
                       const double& kappa, const double& cap_radius,
                       const double& peakI1);

  double YieldFunction(Matrix3& stress, const double& fSlope,
                       const double& kappa, const double& cap_radius,
                       const double& peakI1);

  ////////////////////////////////////////////////////////////////////////
  /* Make the value for pLocalized computed locally available outside of the
   * model. */
  ////////////////////////////////////////////////////////////////////////
  void addRequiresDamageParameter(Task* task, const MPMMaterial* matl,
                                  const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /* Make the value for pLocalized computed locally available outside of the
   * model */
  ////////////////////////////////////////////////////////////////////////
  void getDamageParameter(const Patch* patch, ParticleVariable<int>& damage,
                          int dwi, DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override;

  // carry forward CM data for RigidMPM
  void carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                    DataWarehouse* old_dw, DataWarehouse* new_dw) override;

  // initialize  each particle's constitutive model data
  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;

  void addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patches) const override;

  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, const bool recursion,
                              const bool dummy) const override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  double computeRhoMicroCM(double pressure, const double p_ref,
                           const MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const MPMMaterial* matl, double temperature) override;

  double getCompressibility() override;

private:

  CMData d_cm;
  Vaango::IntVar_SoilModelBrannonKappa* d_intvar;

  void initializeLocalMPMLabels();

  void computeEffectiveModuli(const double& eps_v, double& bulk_modulus,
                              double& lame_modulus) const;
};
} // End namespace Uintah

#endif // __SOIL_MODEL_BRANNON_H__
