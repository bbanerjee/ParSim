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

//  JWLppMPM.h
//     Author: Joseph R. Peterson

#ifndef __JWL_PLUSPLUS_CONSTITUTIVE_MODEL_H__
#define __JWL_PLUSPLUS_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {

class JWLppMPM : public ConstitutiveModel
{

public:

  struct Params
  {
    // Igniition pressure
    double ignition_pressure;

    // Shear viscosity
    double mu;

    // These two parameters are used for the unburned Murnahan EOS
    double K;
    double n;

    // These parameters are used for the product JWL EOS
    double A;
    double B;
    double C;
    double R1;
    double R2;
    double omega;
    // double rho0;

    // These parameters are needed for the reaction model
    double G; // rate coefficient, JWL++
    double b; // pressure exponenet, JWL++
    double
      max_burn_timestep;    // Maximum time increment for burn model subcycling
    double max_burned_frac; // Limit on the fraction that remains unburned
  };

  const VarLabel* pProgressFLabel;
  const VarLabel* pProgressFLabel_preReloc;
  const VarLabel* pProgressdelFLabel;
  const VarLabel* pProgressdelFLabel_preReloc;
  const VarLabel* pLocalizedLabel;
  const VarLabel* pLocalizedLabel_preReloc;

public:
  
  JWLppMPM(ProblemSpecP& ps, MPMFlags* flag);
  JWLppMPM(const JWLppMPM* cm);
  JWLppMPM& operator=(const JWLppMPM& cm) = delete;
  JWLppMPM* clone() override;
  virtual ~JWLppMPM() override;

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
                              const bool schedPar = true) const override;

  /* for material conversion */
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

  /* for MPMICE */
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

  // Initial stress state
  bool d_useInitialStress;
  double d_init_pressure; // Initial pressure

  // Stress update algorithm flag
  bool d_fastCompute;     // true = two stage backward Euler
                          // false = backward Euler with newton iterations
  double d_newtonIterTol; // Tolerance on the norm of [frac pressure] vector
                          // to stop iterations
  int d_newtonIterMax;    // Maximum number of Newton iterations

private:

  //------------------------------------------------------------------
  // Do Newton iterations or two step Backward Euler
  //------------------------------------------------------------------
  void computeUpdatedFractionAndPressure(const double& J_old,
                                         const double& J,
                                         const double& f_old,
                                         const double& p_old,
                                         const double& delT,
                                         double& f_new,
                                         double& p_new) const;

  //------------------------------------------------------------------
  // Two step Backward Euler
  //------------------------------------------------------------------
  void computeWithTwoStageBackwardEuler(const double& J,
                                        const double& f_old,
                                        const double& p_old,
                                        const double& delT,
                                        const double& pM,
                                        const double& pJWL,
                                        double& f_new,
                                        double& p_new) const;

  //------------------------------------------------------------------
  // Newton iterations
  //------------------------------------------------------------------
  void computeWithNewtonIterations(const double& J,
                                   const double& f_old,
                                   const double& p_old,
                                   const double& delT,
                                   const double& pM,
                                   const double& pJWL,
                                   double& f_new,
                                   double& p_new) const;

  //------------------------------------------------------------------
  // Compute G
  //  G = [F_n+1 P_n+1]^T
  //   F_n+1 = 0 = f_n+1 - f_n - G*(1 - f_n+1)*(p_n+1)^b*Delta t
  //   P_n+1 = 0 = p_n+1 - (1 - f_n+1) p_m - f_n+1 p_jwl
  //------------------------------------------------------------------
  void computeG(const double& J,
                const double& f_old,
                const double& f_new,
                const double& p_new,
                const double& pM,
                const double& pJWL,
                const double& delT,
                vector<double>& G) const;

  //------------------------------------------------------------------
  // Compute the Jacobian of G
  //  J_G = [[dF_n+1/df_n+1 dF_n+1/dp_n+1];[dP_n+1/df_n+1 dP_n+1/dp_n+1]]
  //   F_n+1 = 0 = f_n+1 - f_n - G*(1 - f_n+1)*(p_n+1)^b*Delta t
  //   P_n+1 = 0 = p_n+1 - (1 - f_n+1) p_m - f_n+1 p_jwl
  //   dF_n+1/df_n+1 = 1 + G*(p_n+1)^b*Delta t
  //   dF_n+1/dp_n+1 =  b*G*(1 - f_n+1)*(p_n+1)^(b-1)*Delta t
  //   dP_n+1/df_n+1 =  p_m - p_jwl
  //   dP_n+1/dp_n+1 = 1
  //------------------------------------------------------------------
  void computeJacobianG(const double& J,
                        const double& f_new,
                        const double& p_new,
                        const double& pM,
                        const double& pJWL,
                        const double& delT,
                        FastMatrix& JacobianG) const;

  //------------------------------------------------------------------
  //  df/dt = G (1-f) p^b
  //------------------------------------------------------------------
  double computeBurnRate(const double& f, const double& p) const;

  //------------------------------------------------------------------
  //  p_m = (1/nK) [J^(-n) - 1]
  //------------------------------------------------------------------
  double computePressureMurnaghan(const double& J) const;

  //------------------------------------------------------------------
  // p_jwl = A exp(-R1 J) + B exp(-R2 J) + C J^[-(1+omega)]
  //------------------------------------------------------------------
  double computePressureJWL(const double& J) const;

};
} // End namespace Uintah

#endif // __JWL_PLUSPLUS_CONSITUTIVE_MODEL_H__
