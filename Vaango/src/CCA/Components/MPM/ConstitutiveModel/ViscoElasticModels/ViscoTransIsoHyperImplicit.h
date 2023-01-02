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

//  ViscoTransIsoHyperImplicit.h
//  class ConstitutiveModel ConstitutiveModel data type -- 3D -
//  holds ConstitutiveModel
//  information for the FLIP technique:
//    This is for Compressible Transversely isotropic viscoelastic materials
//    Features:
//      Usage:

#ifndef __Visco_Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__
#define __Visco_Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <Core/Disclosure/TypeDescription.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {
class ViscoTransIsoHyperImplicit : public ConstitutiveModel, public ImplicitCM
{
private:
  // Create datatype for storing model parameters
  bool d_useModifiedEOS;
  double d_active;
  string d_StrainEnergy;

public:
  struct CMData
  { //______________________________modified here
    double Bulk;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double lambda_star;
    Vector a0;
    double failure;
    double crit_shear;
    double crit_stretch;
    double y1; // visco properties
    double y2;
    double y3;
    double y4;
    double y5;
    double y6;
    double t1;
    double t2;
    double t3;
    double t4;
    double t5;
    double t6;
  };

  const VarLabel* pStretchLabel;          // For diagnostic
  const VarLabel* pStretchLabel_preReloc; // For diagnostic

  const VarLabel* pFailureLabel; // fail_labels
  const VarLabel* pFailureLabel_preReloc;

  const VarLabel* pElasticStressLabel;
  const VarLabel* pElasticStressLabel_preReloc; // visco stress

  const VarLabel* pHistory1Label;
  const VarLabel* pHistory1Label_preReloc;

  const VarLabel* pHistory2Label;
  const VarLabel* pHistory2Label_preReloc;

  const VarLabel* pHistory3Label;
  const VarLabel* pHistory3Label_preReloc;

  const VarLabel* pHistory4Label;
  const VarLabel* pHistory4Label_preReloc;

  const VarLabel* pHistory5Label;
  const VarLabel* pHistory5Label_preReloc;

  const VarLabel* pHistory6Label;
  const VarLabel* pHistory6Label_preReloc;

private:
  CMData d_initialData;


public:
  // constructors
  ViscoTransIsoHyperImplicit(ProblemSpecP& ps, MPMFlags* flag);
  ViscoTransIsoHyperImplicit(const ViscoTransIsoHyperImplicit* cm);
  ViscoTransIsoHyperImplicit& operator=(const ViscoTransIsoHyperImplicit& cm) = delete;

  // destructor
  ~ViscoTransIsoHyperImplicit() override;

  ModelType modelType() const override
  {
    return ModelType::TOTAL_FORM;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  // clone
  ViscoTransIsoHyperImplicit* clone() override;

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Patch* patch,
                                     const MPMMaterial* matl,
                                     DataWarehouse* new_dw);

  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw, DataWarehouse* new_dw,
                                   Solver* solver,
                                   const bool recursion) override;

  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw) override;

  // initialize  each particle's constitutive model data
  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  void allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                 const PatchSet* patch,
                                 MPMLabel* lb) const override;

  void allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* subset,
                         ParticleLabelVariableMap* newState,
                         ParticleSubset* delset,
                         DataWarehouse* old_dw) override;

  /*virtual void addInitialComputesAndRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*) const;*/

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, const bool recursion,
                              const bool SchedParent) const override;

  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;

  double computeRhoMicroCM(double pressure, const double p_ref,
                           const MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const MPMMaterial* matl, double temperature) override;

  double getCompressibility() override;

  Vector getInitialFiberDir() override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  // const VarLabel* bElBarLabel;
  // const VarLabel* bElBarLabel_preReloc;
};
} // End namespace Uintah

#endif // __Visco_Trans_Iso_Hyper_Implicit_CONSTITUTIVE_MODEL_H__
