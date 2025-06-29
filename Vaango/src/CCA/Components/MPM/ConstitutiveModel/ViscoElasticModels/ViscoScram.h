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

#ifndef __VISCOSCRAM_CONSTITUTIVE_MODEL_H__
#define __VISCOSCRAM_CONSTITUTIVE_MODEL_H__

#include <Core/Math/Short27.h>
#include <Core/Math/Matrix3.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Util/Endian.h>

#include <Core/Disclosure/TypeDescription.h>
#include <Core/ProblemSpec/ProblemSpec.h>


#include <cmath>
#include <vector>

namespace Uintah {

class MPMLabel;
class MPMFlags;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class ViscoScram
  \brief Light version of ViscoSCRAM
  \author Scott Bardenhagen \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n
*/
/////////////////////////////////////////////////////////////////////////////


class ViscoScram : public ConstitutiveModel
{

public:
  struct CMData
  {
    double PR;
    double CoefThermExp;
    double CrackParameterA;
    double CrackPowerValue;
    double CrackMaxGrowthRate;
    double StressIntensityF;
    double CrackFriction;
    double InitialCrackRadius;
    double CrackGrowthRate;
    double G[5];
    double RTau[5];
    double Beta, Gamma;
    double DCp_DTemperature;
    int LoadCurveNumber, NumberOfPoints;
  };

  struct TimeTemperatureData
  {
    double T0_WLF;
    double C1_WLF;
    double C2_WLF;
  };

  // Murnaghan Equation of State Variables
  struct MurnaghanEOS
  {
    double P0;
    double gamma;
    double bulkPrime;
  };
  // JWL Equation of State Variables
  struct JWLEOS
  {
    double A;  // Pa
    double B;  // Pa
    double C;  // Pa
    double Cv; // Pa/K
    double R1;
    double R2;
    double om;
  };

  typedef ViscoScramStateData StateData;

  const VarLabel* pVolChangeHeatRateLabel;
  const VarLabel* pViscousHeatRateLabel;
  const VarLabel* pCrackHeatRateLabel;
  const VarLabel* pCrackRadiusLabel;
  const VarLabel* pStatedataLabel;
  const VarLabel* pRandLabel;
  const VarLabel* pStrainRateLabel;

  const VarLabel* pVolChangeHeatRateLabel_preReloc;
  const VarLabel* pViscousHeatRateLabel_preReloc;
  const VarLabel* pCrackHeatRateLabel_preReloc;
  const VarLabel* pCrackRadiusLabel_preReloc;
  const VarLabel* pStatedataLabel_preReloc;
  const VarLabel* pRandLabel_preReloc;
  const VarLabel* pStrainRateLabel_preReloc;

protected:
  friend const Uintah::TypeDescription* fun_getTypeDescription(
    ViscoScramStateData*);

  // Create datatype for storing model parameters
  bool d_useJWLEOS;
  bool d_useJWLCEOS;
  bool d_useModifiedEOS;
  bool d_useMurnahanEOS;
  bool d_useBirchMurnaghanEOS;
  bool d_random;
  bool d_doTimeTemperature;
  bool d_useObjectiveRate;
  double d_bulk;

  CMData d_initialData;
  TimeTemperatureData d_tt;
  MurnaghanEOS d_murnahanEOSData;
  JWLEOS d_JWLEOSData;


public:
  // constructors
  ViscoScram(ProblemSpecP& ps, MPMFlags* flag);
  ViscoScram(const ViscoScram* cm);
  ViscoScram& operator=(const ViscoScram& cm) = delete;

  // destructor
  ~ViscoScram() override;

  ModelType modelType() const override
  {
    return ModelType::RATE_FORM;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  std::unique_ptr<ConstitutiveModel> clone() override;

  /*! Computes and requires for initialization of history variables */
  void addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patches) const override;

  /*! initialize  each particle's constitutive model data */
  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  /*! compute stable timestep for this patch */
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  /*! Set up data required by and computed in computeStressTensor */
  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;

  /*! Set up data required by and computed in computeStressTensor
      for implicit methods */
  void addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                              [[maybe_unused]] const bool recursion,
                              [[maybe_unused]] const bool schedParent = true) const override
  {
  }

  /*! compute stress at each particle in the patch */
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  /*! carry forward CM data (computed in computeStressTensor) for RigidMPM */
  void carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                    DataWarehouse* old_dw, DataWarehouse* new_dw) override;

  /*! Set up data required in the particle conversion process */
  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;

  /*! Copy data from the delset to the addset in the particle
      conversion process */
  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  /*! Add the particle data that have to be saved at the end of each
      timestep */
  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  /*! Used by MPMICE for pressure equilibriation */
  double computeRhoMicroCM(double pressure, const double p_ref,
                           const MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  /*! Used by MPMICE for pressure equilibriation */
  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const MPMMaterial* matl, double temperature) override;

  /*! Used by MPMICE for pressure equilibriation */
  double getCompressibility() override;

private:
  // Functions and variables for solving the BirchMurnaghan equation of state
  double computePBirchMurnaghan(double v);
  double computedPdrhoBirchMurnaghan(double v, double rho0);

  // Functions and variables for solving JWL temperature dependend form of
  // equation of state
  typedef struct
  {
    double Pressure;
    double Temperature;
    double SpecificHeat;
    double IL, IR;
  } IterationVariables;

  void setInterval(double f, double rhoM, IterationVariables*);
  double computePJWL(double rhoM, double rho0, IterationVariables*);
  double computedPdrhoJWL(double rhoM, double rho0, IterationVariables*);

  void computeRhoRef(const double rho_orig, const double p_ref,
                     const double temperature, const double pressure,
                     double& rho_refrr, double& K0);
};

} // End namespace Uintah

#endif // __VISCOSCRAM_CONSTITUTIVE_MODEL_H__
