/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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

#ifndef __Arenisca3_H__
#define __Arenisca3_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <cmath>

namespace Uintah {
class MPMLabel;
class MPMFlags;

/**************************************

****************************************/

class Arenisca3 : public ConstitutiveModel
{

public:
  static const double one_third;
  static const double two_third;
  static const double four_third;
  static const double sqrt_two;
  static const double one_sqrt_two;
  static const double sqrt_three;
  static const double one_sqrt_three;
  static const double one_sixth;
  static const double one_ninth;
  static const double pi;
  static const double pi_fourth;
  static const double pi_half;
  static const Matrix3 Identity;

  // For usage instructions, see the 'WeibullParser' function
  // header in Kayenta.cc
  struct WeibParameters
  {
    bool Perturb; // 'True' for perturbed parameter
    double
      WeibMed; // Medain distrib. value OR const value depending on bool Perturb
    int WeibSeed;         // seed for random number generator
    double WeibMod;       // Weibull modulus
    double WeibRefVol;    // Reference Volume
    std::string WeibDist; // String for Distribution
  };

  // Create datatype for storing model parameters
  struct CMData
  {
    double PEAKI1;
    double FSLOPE;
    double STREN;
    double YSLOPE;
    double BETA_nonassociativity;
    double B0;
    double B01;
    double B1;
    double B2;
    double B3;
    double B4;
    double G0;
    double G1;
    double G2;
    double G3;
    double G4;
    double p0_crush_curve;
    double p1_crush_curve;
    double p2_crush_curve;
    double p3_crush_curve;
    double CR;
    double fluid_B0;
    double fluid_pressure_initial;
    double T1_rate_dependence;
    double T2_rate_dependence;
    double subcycling_characteristic_number;
    bool Use_Disaggregation_Algorithm;
    double K0_Murnaghan_EOS;
    double n_Murnaghan_EOS;
  };

  const VarLabel* pLocalizedLabel;
  const VarLabel* pLocalizedLabel_preReloc;
  const VarLabel* pAreniscaFlagLabel; // 0: ok, 1: pevp<-p3
  const VarLabel* pAreniscaFlagLabel_preReloc;
  const VarLabel* pScratchDouble1Label;
  const VarLabel* pScratchDouble1Label_preReloc;
  const VarLabel* pScratchDouble2Label;
  const VarLabel* pScratchDouble2Label_preReloc;
  const VarLabel* pPorePressureLabel;
  const VarLabel* pPorePressureLabel_preReloc;
  const VarLabel* pepLabel; // Plastic Strain
  const VarLabel* pepLabel_preReloc;
  const VarLabel* pevpLabel; // Plastic Volumetric Strain
  const VarLabel* pevpLabel_preReloc;
  const VarLabel* peveLabel; // Elastic Volumetric Strain
  const VarLabel* peveLabel_preReloc;
  const VarLabel* pCapXLabel;
  const VarLabel* pCapXLabel_preReloc;
  const VarLabel* pKappaLabel;
  const VarLabel* pKappaLabel_preReloc;
  const VarLabel* pZetaLabel;
  const VarLabel* pZetaLabel_preReloc;
  const VarLabel* pP3Label;
  const VarLabel* pP3Label_preReloc;
  const VarLabel* pStressQSLabel;
  const VarLabel* pStressQSLabel_preReloc;
  const VarLabel* pScratchMatrixLabel;
  const VarLabel* pScratchMatrixLabel_preReloc;

  // weibull parameter set
  WeibParameters wdist;
  const VarLabel* peakI1IDistLabel;
  const VarLabel* peakI1IDistLabel_preReloc;

private:
  double small_number;
  double big_number;

  double d_Kf, d_Km, d_phi_i, d_ev0, d_C1;

  CMData d_cm;

  // For initialization with body force
  bool d_initializeWithBodyForce;
  Point d_surfaceRefPoint;

  void initializeLocalMPMLabels();

public:
  // constructor
  Arenisca3(ProblemSpecP& ps, MPMFlags* flag);
  Arenisca3(const Arenisca3* cm);
  Arenisca3& operator=(const Arenisca3& cm) = delete;

  // destructor
  ~Arenisca3() override;

  ModelType modelType() const override
  {
    return ModelType::RATE_FORM;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  // clone
  Arenisca3* clone() override;

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  // compute stress at each particle in the patch
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

private: // non-uintah mpm constitutive specific functions
  struct AreniscaState
  {
    double capX;   // Hydrostatic compressive strength
    double kappa;  // Branch point
    double zeta;   // Trace of isotropic backstress
    Matrix3 sigma; // Unrotated stress
    Matrix3 ep;    // Plastic strain

    AreniscaState()
    {
      capX = 0.0;
      kappa = 0.0;
      zeta = 0.0;
      sigma = Matrix3(0.0);
      ep = Matrix3(0.0);
    }

    AreniscaState(const double& capX_in, const double& kappa_in,
                  const double& zeta_in, const Matrix3& sigma_in,
                  const Matrix3& ep_in)
    {
      capX = capX_in;
      kappa = kappa_in;
      zeta = zeta_in;
      sigma = sigma_in;
      ep = ep_in;
    }

    AreniscaState(const AreniscaState& state)
    {
      capX = state.capX;
      kappa = state.kappa;
      zeta = state.zeta;
      sigma = state.sigma;
      ep = state.ep;
    }

    AreniscaState& operator=(const AreniscaState& state) = default;

    friend std::ostream& operator<<(std::ostream& os,
                                    const AreniscaState& state)
    {
      os << "I1 = " << state.sigma.Trace() << ", evp = " << state.ep.Trace()
         << ", X = " << state.capX << ", kappa = " << state.kappa
         << ", zeta = " << state.zeta << std::endl;
      return os;
    }
  };

  struct Invariants
  {
    Matrix3 S;  // Deviatoric part
    double I1;  // Hydrostatic part
    double J2;  // 2nd Deviatoric Invariant
    double rJ2; // sqrt(J2)

    Invariants() {}

    Invariants(const Matrix3& stress)
    {

      // Compute the first invariants
      I1 = stress.Trace(); // Pa

      // Compute the deviatoric part of the tensor
      S = stress - one_third * Identity * I1; // Pa

      // Compute the second invariant
      J2 = 0.5 * S.Contract(S); // Pa^2
      J2 = (J2 < 1e-16 * (I1 * I1 + J2)) ? 0.0 : J2;
      rJ2 = sqrt(J2);
    }

    friend std::ostream& operator<<(std::ostream& os, const Invariants& inv)
    {
      os << "I1 = " << inv.I1 << ", sqrt(J2) = " << inv.rJ2 << std::endl;
      return os;
    }
  };

  bool computeStep(particleIndex idx, long64 particleID, const Matrix3& D,
                   const double& dt, const AreniscaState& state_n,
                   const double& coher, const double& P3,
                   AreniscaState& state_p, long64 ParticleID);

  void computeElasticProperties(double& bulk, double& shear);

  void computeElasticProperties(const AreniscaState& state, const double& P3,
                                double& bulk, double& shear);

  void computeElasticProperties(const Matrix3& sigma, const Matrix3& ep,
                                const double& P3, double& bulk, double& shear);

  Matrix3 computeTrialStress(const Matrix3& sigma_old, // old stress
                             const Matrix3& d_e,       // Strain increment
                             const double& bulk,       // bulk modulus
                             const double& shear);     // shear modulus

  int computeStepDivisions(particleIndex idx, long64 particleID,
                           const AreniscaState& state, const double& P3,
                           const Matrix3& sigma_trial);

  bool computeSubstep(
    particleIndex idx, long64 particleID,
    const Matrix3& d_e,             // total strain increment for substep
    const AreniscaState& state_old, // state at start of substep
    const double& coher,            // scalar valued coher
    const double& P3,               // initial disaggregation strain
    AreniscaState& state_new        // state at end of substep
    );

  double computeX(const double& evp, const double& P3);

  double computedZetadevp(double Zeta, double evp);

  double computePorePressure(const double ev);

  int nonHardeningReturn(const Invariants& invar_trial,
                         const Invariants& invar_old, const Matrix3& d_e,
                         const AreniscaState& state, const double& coher,
                         const double& bulk, const double& shear,
                         Invariants& invar_new, Matrix3& d_ep_new,
                         double& kappa);

  void transformedBisection(double& z_0, double& r_0, const double& z_trial,
                            const double& r_trial, const AreniscaState& state,
                            const double& coher,
                            const double limitParameters[4], // XXX
                            const double& r_to_rJ2, double& kappa);

  int transformedYieldFunction(const double& z, const double& r,
                               const AreniscaState& state, const double& coher,
                               const double limitParameters[4], // XXX
                               const double& r_to_rJ2, double& kappa);

  int computeYieldFunction(const Invariants& invariants,
                           const AreniscaState& state, const double& coher,
                           const double limitParameters[4],
                           double& kappa // XXX
                           );

  int computeYieldFunction(const double& I1, const double& rJ2,
                           const AreniscaState& state, const double& coher,
                           const double limitParameters[4],
                           double& kappa // XXX
                           );

  void computeLimitParameters(double* limitParameters,
                              const double& coher // XXX
                              );
  void checkInputParameters();

public: // Uintah MPM constitutive model specific functions
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

  ///////////////////////////////////////////////////////////////////////
  /*!
   * Actually initialize the stress and deformation gradient assuming linear
   * elastic behavior after computing the body force acceleration
   *
   * **WARNING** 1) Assumes zero shear stresses and that body forces are aligned
   *                with coordinate directions
   *             2) Needs the model to have a "initializeWithBodyForce" flag
   *                set as true.  A more general implementation is not worth
   *                the significant extra effort.
   */
  ///////////////////////////////////////////////////////////////////////
  void initializeStressAndDefGradFromBodyForce(
    const Patch* patch, const MPMMaterial* matl,
    DataWarehouse* new_dw) const override;

  void addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patches) const override;

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

  // Weibull input parser that accepts a structure of input
  // parameters defined as:
  //
  // bool Perturb        'True' for perturbed parameter
  // double WeibMed       Medain distrib. value OR const value
  //                         depending on bool Perturb
  // double WeibMod       Weibull modulus
  // double WeibScale     Scale parameter
  // std::string WeibDist  String for Distribution
  virtual void WeibullParser(WeibParameters& iP);

  /*! This is for adding/deleting particles when a particle is switched
      from one material to another */
  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;
};
} // End namespace Uintah

#endif // __Arenisca3_H__
