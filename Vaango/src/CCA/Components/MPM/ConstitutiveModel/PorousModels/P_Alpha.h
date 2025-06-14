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

#ifndef __P_ALPHA_CONSTITUTIVE_MODEL_H__
#define __P_ALPHA_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {

class P_Alpha : public ConstitutiveModel
{
public:
  struct Params
  {
    // For P-alpha response
    double Ps;     // Press. at which material reaches full density (alpha=1)
    double Pe;     // Press. at which initial elastic response starts to yield
    double rhoS;   // Solid material density (corresponds to Ps)
    double alpha0; // Initial value of alpha for virgin material
    double K0;     // Initial bulk modulus in elastic region
    double Ks;     // Bulk modulus of fully densified material
    double Ku;     // Bulk modulus in unloading for alpha > alpha0, or p <= 0
                   // Ku defaults to .1*K0
    // For Mie-Gruneisen response
    double T_0;
    double C_0;
    double Gamma_0;
    double S_alpha;
  };

  const VarLabel* pAlphaLabel;
  const VarLabel* pAlphaMinLabel;
  const VarLabel* pAlphaMinLabel_preReloc;
  const VarLabel* pTempAlpha1Label;
  const VarLabel* pTempAlpha1Label_preReloc;

public:

  P_Alpha(ProblemSpecP& ps, MPMFlags* flag);
  P_Alpha(const P_Alpha* cm);
  P_Alpha& operator=(const P_Alpha& cm) = delete;
  virtual ~P_Alpha() override;

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

  void addComputesAndRequires([[maybe_unused]] Task* task,
                              [[maybe_unused]] const MPMMaterial* matl,
                              [[maybe_unused]] const PatchSet* patches,
                              [[maybe_unused]] const bool recursion,
                              [[maybe_unused]] const bool schedParent = true) const override {};

  /* For material conversion */
  void allocateCMDataAddRequires(Task* task,
                                 const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;
  void allocateCMDataAdd(DataWarehouse* new_dw,
                         ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

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

private:

  Params d_modelParam;

};
} // End namespace Uintah

#endif // __P_ALPHA_CONSTITUTIVE_MODEL_H__
