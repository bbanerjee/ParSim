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

#ifndef __MURNAGHAN_CONSTITUTIVE_MODEL_H__
#define __MURNAGHAN_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <Core/Math/Matrix3.h>
#include <cmath>
#include <vector>

namespace Uintah {

/**
 *  Equation of state outlined in:
 *  F.D. Murnaghan, 'The Compressibility of Media under Extreme Pressures',
 *  in Proceedings of the National Academy of Sciences, vol. 30, pp. 244-247,
 *  1944.
 * 
 *  P = P0 + (K/K')*((V/V0)^(-K')-1)
 *
 *  K          (Pa)         Bulk modulus of the material
 *  K' = dK/dp (unitless)   Slope of bulk modulus vs. pressure
 *  V0         (m^3)        Reference volume of material in normal state
 *  V          (m^3)        Current volume of material in normal state
 *  P0         (Pa)         Reference Pressure
 *  rho0       (kg/m^3)     Reference density
 *  viscosity
 */ 

class MurnaghanMPM : public ConstitutiveModel
{

public:
  struct Params
  {             
    double d_K;
    double d_Kprime;
    double d_P0;
    double d_rho0;
    double d_viscosity;
  };

  MurnaghanMPM(ProblemSpecP& ps, MPMFlags* flag);
  MurnaghanMPM(const MurnaghanMPM* cm);
  MurnaghanMPM& operator=(const MurnaghanMPM& cm) = delete;
  virtual ~MurnaghanMPM() override = default;
  MurnaghanMPM* clone() override;

  ModelType modelType() const override { return ModelType::TOTAL_FORM; }

  void outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag = true) override;

  void addParticleState(std::vector<const VarLabel*>& from,
                        std::vector<const VarLabel*>& to) override;

  /* initialize */
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
                              const bool schedParent = true) const override;

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
  Params d_modelParam;
  bool d_useModifiedEOS;
  int d_8or27;

};
} // End namespace Uintah

#endif // __MURNAGHAN_CONSTITUTIVE_MODEL_H__
