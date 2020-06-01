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

#ifndef __TransIsoHyper_CONSTITUTIVE_MODEL_H__
#define __TransIsoHyper_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {
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
  TransIsoHyper* clone() override;

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

  void addComputesAndRequires(Task* task,
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

private:
  bool d_useModifiedEOS;
  CMData d_modelParam;

  Vector computeEigenvalues(const Matrix3& C) const;

};
} // End namespace Uintah

#endif // __TransIsoHyper_CONSTITUTIVE_MODEL_H__
