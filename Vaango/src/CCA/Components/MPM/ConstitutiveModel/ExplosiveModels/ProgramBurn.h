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

#ifndef __PROGRAM_BURN_CONSTITUTIVE_MODEL_H__
#define __PROGRAM_BURN_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {
class ProgramBurn : public ConstitutiveModel
{

public:

  struct Params
  {
    // These two parameters are used for the unburned Murnahan EOS
    double d_K;
    double d_n;

    // These parameters are used for the product JWL EOS
    double d_A;
    double d_B;
    double d_C;
    double d_R1;
    double d_R2;
    double d_om;
    double d_rho0;

    // These parameters are needed for the reaction model
    Point d_start_place; // Starting point of the detonation
    Vector d_direction;  // Direction if starting from a plane (point-normal)
    double d_D;          // Detonation velocity
    double d_T0;         // Detonation initiation time
  };

  const VarLabel* pProgressFLabel;
  const VarLabel* pProgressFLabel_preReloc;
  const VarLabel* pLocalizedLabel;
  const VarLabel* pLocalizedLabel_preReloc;

  ProgramBurn(ProblemSpecP& ps, MPMFlags* flag);
  ProgramBurn(const ProgramBurn* cm);
  ProgramBurn& operator=(const ProgramBurn& cm) = delete;
  ~ProgramBurn() override;
  
  std::unique_ptr<ConstitutiveModel> clone() override;

  ModelType modelType() const override { return ModelType::TOTAL_FORM; }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  /* initialize */
  void addInitialComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet* patches) const override;
  void initializeCMData(const Patch* patch,
                        const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  /* compute stress */
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
                              const bool SchedParent) const override;

  /* for particle deletion */
  void addRequiresDamageParameter(Task* task,
                                  const MPMMaterial* matl,
                                  const PatchSet* patches) const override;
  void getDamageParameter(const Patch* patch,
                          ParticleVariable<int>& damage,
                          int dwi,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override;

  /* for particle conversion */
  void allocateCMDataAddRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;
  void allocateCMDataAdd(DataWarehouse* new_dw,
                         ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  /* for RigidMPM */
  void carryForward(const PatchSubset* patches,
                    const MPMMaterial* matl,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw) override;

  /* For MPMICE */
  double computeRhoMicroCM(double pressure,
                           const double p_ref,
                           const MPMMaterial* matl,
                           double temperature,
                           double rho_guess) override;
  void computePressEOSCM(double rho_m,
                         double& press_eos,
                         double p_ref,
                         double& dp_drho,
                         double& ss_new,
                         const MPMMaterial* matl,
                         double temperature) override;
  double getCompressibility() override;

protected:

  Params d_cm;
  bool d_useModifiedEOS;
  int d_8or27;

};
} // End namespace Uintah

#endif // __PROGRAM_BURN_CONSTITUTIVE_MODEL_H__
