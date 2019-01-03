/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

//  UCNH.h
//  class ConstitutiveModel ConstitutiveModel data type -- 3D -
//  holds ConstitutiveModel
//  information for the FLIP technique:
//    This is for Compressible NeoHookean materials
//    Features:
//      Usage:

#ifndef __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__
#define __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__

#include <Core/Math/Short27.h>
#include <Core/Math/Matrix3.h>

namespace Uintah {
// Structures for Plasticitity

struct UCNHStateData
{
  double Alpha;
};
class TypeDescription;
const TypeDescription* fun_getTypeDescription(UCNHStateData*);
}

#include <Core/Util/Endian.h>

namespace Uintah {
using namespace Uintah;
inline void
swapbytes(Uintah::UCNHStateData& d)
{
  swapbytes(d.Alpha);
}
} // namespace Uintah

#include "ConstitutiveModel.h"
#include "ImplicitCM.h"
#include "PlasticityModels/MPMEquationOfState.h"
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cmath>
#include <vector>

namespace Uintah {
// Classes needed by UCNH
class TypeDescription;

class UCNH : public ConstitutiveModel, public ImplicitCM
{

  ///////////////
  // Variables //
  ///////////////
public:
  // Basic Requirements //
  ////////////////////////
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
  MPMEquationOfState* d_eos;

public:
  // constructors
  UCNH(ProblemSpecP& ps, MPMFlags* flag);
  UCNH(ProblemSpecP& ps, MPMFlags* flag, bool plas, bool dam);
  UCNH(const UCNH* cm);
  UCNH& operator=(const UCNH& cm) = delete;

  ModelType modelType() const override
  {
    return ModelType::TOTAL_FORM;
  }

  // specify what to output from the constitutive model to an .xml file
  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  // clone
  UCNH* clone() override;

  // destructor
  ~UCNH() override;

  // Initialization Functions //
  //////////////////////////////
  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;

  // carry forward CM data for RigidMPM
  void carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                    DataWarehouse* old_dw, DataWarehouse* new_dw) override;

  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  // Scheduling Functions //
  //////////////////////////
  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, const bool recursion,
                              const bool schedPar = true) const override;

  void addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                     const PatchSet* patches) const override;

  ////////////////////////////////////////////////////////////////////////
  /*! \\brief Add the requires for failure simulation. */
  ////////////////////////////////////////////////////////////////////////
  void addRequiresDamageParameter(Task* task, const MPMMaterial* matl,
                                  const PatchSet* patches) const override;

  // Compute Functions //
  ///////////////////////
  // main computation of pressure from constitutive model's equation of state
  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const MPMMaterial* matl, double temperature) override;

  // main computation of density from constitutive model's equation of state
  double computeRhoMicroCM(double pressure, const double p_ref,
                           const MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  // compute stress at each particle in the patch
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  // Damage specific CST for solver
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw, DataWarehouse* new_dw,
                                   Solver* solver, const bool) override;

  // Helper Functions //
  //////////////////////
  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  // Returns the compressibility of the material
  double getCompressibility() override;

  ////////////////////////////////////////////////////////////////////////
  /*! \\brief Get the flag that marks a failed particle. */
  ////////////////////////////////////////////////////////////////////////
  void getDamageParameter(const Patch* patch, ParticleVariable<int>& damage,
                          int dwi, DataWarehouse* old_dw,
                          DataWarehouse* new_dw) override;

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
  void computeTangentStiffnessMatrix(const Matrix3& sigDev, const double& mubar,
                                     const double& J, const double& bulk,
                                     double D[6][6]);
  /*! Compute BT*Sig*B (KGeo) */
  void BnlTSigBnl(const Matrix3& sig, const double Bnl[3][24],
                  double BnTsigBn[24][24]) const;

  /*! Compute K matrix */
  void computeStiffnessMatrix(const double B[6][24], const double Bnl[3][24],
                              const double D[6][6], const Matrix3& sig,
                              const double& vol_old, const double& vol_new,
                              double Kmatrix[24][24]);
};
} // End namespace Uintah

#endif // __UNIFIED_NEOHOOK_CONSTITUTIVE_MODEL_H__
