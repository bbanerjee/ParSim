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

#ifndef __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__
#define __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__

#include <Core/Math/Short27.h>
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <CCA/Components/MPM/ConstitutiveModel/EOSModels/MPMEquationOfState.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cmath>
#include <vector>
#include <Core/Math/Matrix3.h>

namespace Uintah {

class UCNH : public ConstitutiveModel, public ImplicitCM
{

public:
  // Create datatype for storing model parameters
  struct CMData
  {
    double Bulk;
    double tauDev;
    // For Plasticity
    double FlowStress;
    double K;
    double Alpha;
  };

  struct YieldDistribution
  {
    std::string dist;
    double range;
    int seed;
  };

  const VarLabel* bElBarLabel;
  const VarLabel* bElBarLabel_preReloc;
  const VarLabel* pDeformRateLabel;
  const VarLabel* pDeformRateLabel_preReloc;

  // Plasticity Requirements //
  /////////////////////////////
  const VarLabel* pPlasticStrain_label;
  const VarLabel* pPlasticStrain_label_preReloc;
  const VarLabel* pYieldStress_label;
  const VarLabel* pYieldStress_label_preReloc;

protected:
  // Flags indicating if plasticity should be used
  bool d_usePlasticity;

  // Basic Requirements //
  ////////////////////////
  CMData d_initialData;
  bool d_useModifiedEOS;
  int d_8or27;

  // Damage Requirments //
  ////////////////////////
  YieldDistribution d_yield;

  // Initial stress state
  bool d_useInitialStress;
  double d_init_pressure; // Initial pressure

  // Model factories
  // bool d_useEOSFactory;
  Vaango::MPMEquationOfState* d_eos;

public:
  UCNH(ProblemSpecP& ps, MPMFlags* flag);
  UCNH(ProblemSpecP& ps, MPMFlags* flag, bool plas, bool dam);
  UCNH(const UCNH* cm);
  UCNH& operator=(const UCNH& cm) = delete;
  ~UCNH() override;

  UCNH* clone() override;

  ModelType modelType() const override { return ModelType::TOTAL_FORM; }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  void addInitialComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet* patches) const override;
  void initializeCMData(const Patch* patch,
                        const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  void addComputesAndRequires(Task* task,
                              const MPMMaterial* matl,
                              const PatchSet* patches) const override;
  void computeStressTensor(const PatchSubset* patches,
                           const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  void addComputesAndRequires(Task* task,
                              const MPMMaterial* matl,
                              const PatchSet* patches,
                              const bool recursion,
                              const bool schedPar = true) const override;
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   Solver* solver,
                                   const bool) override;

  ////////////////////////////////////////////////////////////////////////
  /*! \\brief Get the flag that marks a failed particle. */
  ////////////////////////////////////////////////////////////////////////
  void addRequiresDamageParameter(Task* task,
                                  const MPMMaterial* matl,
                                  const PatchSet* patches) const override;
  void getDamageParameter(const Patch* patch,
                          ParticleVariable<int>& damage,
                          int dwi,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override;

  /* For particle material change after failure */
  void allocateCMDataAddRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;
  void allocateCMDataAdd(DataWarehouse* new_dw,
                         ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  // carry forward CM data for RigidMPM
  void carryForward(const PatchSubset* patches,
                    const MPMMaterial* matl,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw) override;

  void computePressEOSCM(double rho_m,
                         double& press_eos,
                         double p_ref,
                         double& dp_drho,
                         double& ss_new,
                         const MPMMaterial* matl,
                         double temperature) override;
  double computeRhoMicroCM(double pressure,
                           const double p_ref,
                           const MPMMaterial* matl,
                           double temperature,
                           double rho_guess) override;
  double getCompressibility() override;

private:
  // Damage requirements //
  /////////////////////////
  void getYieldStressDistribution(ProblemSpecP& ps);

  void setYieldStressDistribution(const UCNH* cm);

protected:
  // compute stress at each particle in the patch
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw) override;

  /*! Compute tangent stiffness matrix */
  void computeTangentStiffnessMatrix(const Matrix3& sigDev,
                                     const double& mubar,
                                     const double& J,
                                     const double& bulk,
                                     double D[6][6]);
  /*! Compute BT*Sig*B (KGeo) */
  void BnlTSigBnl(const Matrix3& sig,
                  const double Bnl[3][24],
                  double BnTsigBn[24][24]) const;

  /*! Compute K matrix */
  void computeStiffnessMatrix(const double B[6][24],
                              const double Bnl[3][24],
                              const double D[6][6],
                              const Matrix3& sig,
                              const double& vol_old,
                              const double& vol_new,
                              double Kmatrix[24][24]);
};
} // End namespace Uintah

#endif // __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__
