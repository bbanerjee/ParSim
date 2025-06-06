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

#ifndef __RIGID_CONSTITUTIVE_MODEL_H__
#define __RIGID_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <vector>

namespace Uintah {

class MPMLabel;
class MPMFlags;

/////////////////////////////////////////////////////////////////////////////
/*!
  \class RigidMaterial
  \brief Rigid material - no stresses or deformation.
  \author Biswajit Banerjee \n
  C-SAFE and Department of Mechanical Engineering \n
  University of Utah \n

  The material does not deform and does not develop any stresses or
  internal heating.

  Shear and bulk moduli are used to compute the wave speed and contact
  and interaction with ICE.
*/
/////////////////////////////////////////////////////////////////////////////

class RigidMaterial : public ConstitutiveModel, public ImplicitCM
{

public:
  struct CMData
  {
    double G;
    double K;
  };

private:
  CMData d_initialData;

  // Prevent assignment of this class
  RigidMaterial& operator=(const RigidMaterial& cm);

public:
  // constructors
  RigidMaterial(ProblemSpecP& ps, MPMFlags* flag);
  RigidMaterial(const RigidMaterial* cm);

  // destructor
  ~RigidMaterial() override;

  ModelType modelType() const override
  {
    return ModelType::RATE_FORM;
  }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  std::unique_ptr<ConstitutiveModel> clone() override;

  /*! initialize  each particle's constitutive model data */
  void initializeCMData(const Patch* patch, const MPMMaterial* matl,
                        DataWarehouse* new_dw) override;

  /*! Computes and requires for compute stress tensor added to
    the taskgraph */
  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches) const override;

  /*! compute stress at each particle in the patch */
  void computeStressTensor(const PatchSubset* patches, const MPMMaterial* matl,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw) override;

  /* Add computes and requires for the implicit code */
  void addComputesAndRequires(Task* task, const MPMMaterial* matl,
                              const PatchSet* patches, const bool recursion,
                              const bool SchedParent) const override;

  /* Computes stress tensor for the implicit code */
  void computeStressTensorImplicit(const PatchSubset*, const MPMMaterial*,
                                   DataWarehouse*, DataWarehouse*, Solver*,
                                   const bool) override;

  /*! carry forward CM data for RigidMPM */
  void carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                    DataWarehouse* old_dw, DataWarehouse* new_dw) override;

  /*! Add requires to task graph for particle conversion */
  void allocateCMDataAddRequires(Task*, const MPMMaterial*, const PatchSet*,
                                 MPMLabel*) const override
  {
  }

  /*! Add requires to task graph for particle conversion */
  void allocateCMDataAdd(DataWarehouse*, ParticleSubset*,
                         ParticleLabelVariableMap*,
                         ParticleSubset*, DataWarehouse*) override
  {
  }

  /*! Add particle state for relocation */
  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  /*! Function that interacts with ice */
  double computeRhoMicroCM(double pressure, const double p_ref,
                           const MPMMaterial* matl, double temperature,
                           double rho_guess) override;

  /*! Function that interacts with ice */
  void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                         double& dp_drho, double& ss_new,
                         const MPMMaterial* matl, double temperature) override;

  /*! Function that interacts with ice */
  double getCompressibility() override;

protected:
  /*! compute stress at each particle in the patch (replacement for
      standard compute stress tensor without the recursion flag) */
  void computeStressTensorImplicit(const PatchSubset* patches,
                                   const MPMMaterial* matl,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw) override;
};

} // End namespace Uintah

#endif // __RIGID_CONSTITUTIVE_MODEL_H__
