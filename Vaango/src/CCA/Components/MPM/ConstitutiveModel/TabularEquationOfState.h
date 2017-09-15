/*
 * The MIT License
 *
 * Copyright (c) 2017- Parresia Research Limited, New Zealand
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

#ifndef __TABULAR_EOS_CONSTITUTIVE_MODEL_H__
#define __TABULAR_EOS_CONSTITUTIVE_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Math/Matrix3.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <cmath>

namespace Uintah {

class MPMLabel;
class MPMFlags;
}

namespace Vaango {

class TabularEquationOfState : public Uintah::ConstitutiveModel
{

private:
  Vaango::TabularData d_table;

  // Prevent copying of this class
  // copy constructor
  // TabularEquationOfState(const TabularEquationOfState &cm);
  TabularEquationOfState& operator=(const TabularEquationOfState& cm) = delete;

public:
  // constructor
  TabularEquationOfState(Uintah::ProblemSpecP& ps, Uintah::MPMFlags* flag);
  TabularEquationOfState(const TabularEquationOfState* cm);

  // destructor
  virtual ~TabularEquationOfState();

  virtual void outputProblemSpec(Uintah::ProblemSpecP& ps, bool output_cm_tag = true);

  // clone
  TabularEquationOfState* clone();

  // compute stable timestep for this patch
  virtual void computeStableTimestep(const Uintah::Patch* patch,
                                     const Uintah::MPMMaterial* matl,
                                     Uintah::DataWarehouse* new_dw);

  // compute stress at each particle in the patch
  virtual void computeStressTensor(const Uintah::PatchSubset* patches,
                                   const Uintah::MPMMaterial* matl,
                                   Uintah::DataWarehouse* old_dw,
                                   Uintah::DataWarehouse* new_dw);

  // carry forward CM data for RigidMPM
  virtual void carryForward(const Uintah::PatchSubset* patches, const Uintah::MPMMaterial* matl,
                            Uintah::DataWarehouse* old_dw, Uintah::DataWarehouse* new_dw);

  virtual double computeRhoMicroCM(double pressure, const double p_ref,
                                   const Uintah::MPMMaterial* matl, double temperature,
                                   double rho_guess);

  virtual void computePressEOSCM(double rho_m, double& press_eos, double p_ref,
                                 double& dp_drho, double& ss_new,
                                 const Uintah::MPMMaterial* matl, double temperature);

  virtual double getCompressibility();

  double computeBulkModulus(const double& rho0,
                            const double& rho) const;

  // initialize  each particle's constitutive model data
  virtual void initializeCMData(const Uintah::Patch* patch, const Uintah::MPMMaterial* matl,
                                Uintah::DataWarehouse* new_dw);

  virtual void allocateCMDataAddRequires(Uintah::Task* task, const Uintah::MPMMaterial* matl,
                                         const Uintah::PatchSet* patch,
                                         Uintah::MPMLabel* lb) const;

  virtual void allocateCMDataAdd(Uintah::DataWarehouse* new_dw, Uintah::ParticleSubset* addset,
                                 Uintah::ParticleLabelVariableMap* newState,
                                 Uintah::ParticleSubset* delset, Uintah::DataWarehouse* old_dw);

  virtual void addComputesAndRequires(Uintah::Task* task, const Uintah::MPMMaterial* matl,
                                      const Uintah::PatchSet* patches) const;

  virtual void addComputesAndRequires(Uintah::Task* task, const Uintah::MPMMaterial* matl,
                                      const Uintah::PatchSet* patches,
                                      const bool recursion,
                                      const bool schedParent = true) const;

  virtual void addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                std::vector<const Uintah::VarLabel*>& to);

};
} // End namespace Vaango

#endif // __TABULAR_EOS_CONSTITUTIVE_MODEL_H__
