/*
 * The MIT License
 *
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

#ifndef __DEFORMATION_GRADIENT_COMPUTER_H__
#define __DEFORMATION_GRADIENT_COMPUTER_H__

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Variables/ComputeSet.h>
#include <Core/Grid/Variables/ParticleSubset.h>
#include <Core/Grid/Variables/ParticleVariableBase.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <vector>

namespace Uintah {

class Task;
class Patch;

//////////////////////////////////////////////////////////////////////////
/*!
  \class DeformationGradientComputer
  \brief Class for computing deformation gradients
*/
//////////////////////////////////////////////////////////////////////////

class DeformationGradientComputer
{

public:
  DeformationGradientComputer(MaterialManagerP& ss,
                              const MPMLabel* mpm_labels,
                              const MPMFlags* mpm_flags);
  DeformationGradientComputer(const DeformationGradientComputer* gc);
  virtual ~DeformationGradientComputer();

  // Make a clone of the gradient computer
  DeformationGradientComputer*
  clone();

  // Computes and requires
  void
  addInitialComputesAndRequires(Task* task,
                                const MPMMaterial* mpm_matl,
                                const PatchSet*);

  void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* mpm_matl,
                         const PatchSet*);

  void
  addComputesOnly(Task* task, const MPMMaterial* mpm_matl, const PatchSet*);

  void
  addComputesAndRequires(Task* task,
                         const MPMMaterial* matl,
                         const PatchSet* patches,
                         const bool /*recurse*/,
                         const bool SchedParent) const;

  void
  addRequiresForConvert(Task* task, const MPMMaterial* mpm_matl);

  void
  copyAndDeleteForConvert(
    DataWarehouse* new_dw,
    ParticleSubset* addset,
    std::map<const VarLabel*, ParticleVariableBase*>* newState,
    ParticleSubset* delset,
    DataWarehouse* old_dw);

  void
  initializeGradient(const Patch* patch,
                     const MPMMaterial* mpm_matl,
                     DataWarehouse* new_dw);

  void
  computeDeformationGradient(const PatchSubset* patches,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw);

  void
  computeDeformationGradient(const PatchSubset* patches,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             const bool recurse);

protected:
  void
  addComputesAndRequiresExplicit(Task* task, const MPMMaterial* mpm_matl);

  void
  addComputesAndRequiresImplicit(Task* task, const MPMMaterial* mpm_matl);

  void
  initializeGradientExplicit(const Patch* patch,
                             const MPMMaterial* mpm_matl,
                             DataWarehouse* new_dw);

  void
  initializeGradientImplicit(const Patch* patch,
                             const MPMMaterial* mpm_matl,
                             DataWarehouse* new_dw);

  void
  computeDeformationGradientExplicit(const Patch* patch,
                                     const MPMMaterial* mpm_matl,
                                     const double& delT,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);

  void
  computeDeformationGradientImplicit(const Patch* patch,
                                     const MPMMaterial* mpm_matl,
                                     const double& delT,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw);

  void
  computeDeformationGradientImplicit(const Patch* patch,
                                     const MPMMaterial* mpm_matl,
                                     const double& delT,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* parent_old_dw,
                                     DataWarehouse* new_dw);

  void
  computeDeformationGradientFromVelocity(const Matrix3& velGrad_old,
                                         const Matrix3& velGrad_new,
                                         const Matrix3& defGrad_old,
                                         const double& delT,
                                         Matrix3& defGrad_new,
                                         Matrix3& defGrad_inc);

  void
  computeDeformationGradientFromTotalDisplacement(const Matrix3& dispGrad_new,
                                                  const Matrix3& defGrad_old,
                                                  Matrix3& defGrad_new,
                                                  Matrix3& defGrad_inc);

  void
  seriesUpdateConstantVelGrad(const Matrix3& velGrad_new,
                              const Matrix3& defGrad_old,
                              const double& delT,
                              Matrix3& defGrad_new,
                              Matrix3& defGrad_inc);

  void
  cayleyUpdateConstantVelGrad(const Matrix3& velGrad_new,
                              const Matrix3& defGrad_old,
                              const double& delT,
                              Matrix3& defGrad_new,
                              Matrix3& defGrad_inc);

  void
  subcycleUpdateConstantVelGrad(const Matrix3& velGrad_new,
                                const Matrix3& defGrad_old,
                                const double& delT,
                                Matrix3& defGrad_new,
                                Matrix3& defGrad_inc);

  void
  computeDeformationGradientFromIncrementalDisplacement(
    const Matrix3& dispGrad_new,
    const Matrix3& defGrad_old,
    Matrix3& defGrad_new,
    Matrix3& defGrad_inc);

protected:
  const MPMLabel* lb;
  const MPMFlags* flag;
  int NGP;
  int NGN;
  MaterialManagerP d_mat_manager;

  static const Matrix3 Identity;
  static const Matrix3 Zero;
};

} // End namespace Uintah

#endif // __DEFORMATION_GRADIENT_COMPUTER_H__
