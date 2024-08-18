/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

#ifndef __TransIsoHyper_CONSTITUTIVE_MODEL_H__
#define __TransIsoHyper_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {
class TransIsoHyperState
{
public:
  double J;
  double I1_bar;
  double I2_bar;
  double I4_bar;
  double I4;
  double lambda_bar;
  double dWdI4_bar;
  double d2WdI4_bar2;
  double shearModulus;
  double U;
  double W;
  Vector fiberDir_new;
  Matrix3 C;
  Matrix3 B;
  Matrix3 C_bar;
  Matrix3 B_bar;
};

class TransIsoHyper : public ConstitutiveModel
{
public:
  struct CMData
  {
    double bulkModulus;
    Vector a0;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double c6;
    double lambdaStar;
    double failure;
    double critShear;
    double critStretch;
  };

  const VarLabel* pStretchLabel;
  const VarLabel* pStretchLabel_preReloc;

  const VarLabel* pFailureLabel;
  const VarLabel* pFailureLabel_preReloc;

  TransIsoHyper(ProblemSpecP& ps, MPMFlags* flag);
  TransIsoHyper(const TransIsoHyper* cm);
  TransIsoHyper& operator=(const TransIsoHyper& cm) = delete;
  virtual ~TransIsoHyper() override;
  
  std::unique_ptr<ConstitutiveModel> clone() override;

  ModelType modelType() const override { return ModelType::TOTAL_FORM; }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  void addInitialComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet*) const override;
  void initializeCMData(const Patch* patch,
                        const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;
  void computeStableTimestep(const Patch* patch,
                             const MPMMaterial* matl,
                             DataWarehouse* new_dw);

  void addComputesAndRequires(Task* task,
                              const MPMMaterial* matl,
                              const PatchSet* patches) const override;
  void computeStressTensor(const PatchSubset* patches,
                           const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  virtual void addComputesAndRequires(
    Task* task,
    const MPMMaterial* matl,
    const PatchSet* patches,
    const bool recursion,
    const bool schedParent = true) const override;

  // carry forward CM data for RigidMPM
  void carryForward(const PatchSubset* patches,
                    const MPMMaterial* matl,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw) override;

  void allocateCMDataAddRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;
  void allocateCMDataAdd(DataWarehouse* new_dw,
                         ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

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

  Vector getInitialFiberDir() override;

protected:
  bool d_useModifiedEOS;
  CMData d_param;

  TransIsoHyperState computeDeformationState(const Matrix3& defGrad,
                                             const Vector& fiberDir) const;
  void computeStressWithFailure(const TransIsoHyperState& ss,
                                double pStretch,
                                double pFailure_old,
                                double& pFailure,
                                Matrix3& hydrostatic_stress,
                                Matrix3& deviatoric_stress,
                                Matrix3& fiber_stress) const;
  Vector computeEigenvalues(const Matrix3& C) const;
  Matrix3 computeHydrostaticStress(const TransIsoHyperState& ss) const;
  Matrix3 computeDeviatoricStress(const TransIsoHyperState& ss) const;
  Matrix3 computeFiberStress(const TransIsoHyperState& ss) const;
};
} // End namespace Uintah

#endif // __TransIsoHyper_CONSTITUTIVE_MODEL_H__
