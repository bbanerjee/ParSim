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

#ifndef __DEVSTRESSMODEL_H__
#define __DEVSTRESSMODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ModelState/DeformationState.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Math/Matrix3.h>

namespace Uintah {

///////////////////////////////////////////////////////////////////////////
/*!
  \class  DevStressModel
  \brief  Abstract Base class for Deviatoric Stress models
  \author Todd Harman
*/
///////////////////////////////////////////////////////////////////////////

using Vaango::ModelStateBase;

class DevStressModel
{

private:
public:
  DevStressModel();
  virtual ~DevStressModel();

  virtual void outputProblemSpec(ProblemSpecP& ps) = 0;

  // Computes and requires for internal evolution variables
  virtual void
  addInitialComputesAndRequires([[maybe_unused]] Task* task,
                                [[maybe_unused]] const MPMMaterial* matl) {};

  virtual void
  addComputesAndRequires([[maybe_unused]] Task* task,
                         [[maybe_unused]] const MPMMaterial* matl) {};

  virtual void
  addComputesAndRequires([[maybe_unused]] Task* task,
                         [[maybe_unused]] const MPMMaterial* matl,
                         [[maybe_unused]] bool SchedParent) {};

  virtual void
  addParticleState([[maybe_unused]] std::vector<const VarLabel*>& from,
                   [[maybe_unused]] std::vector<const VarLabel*>& to) {};

  virtual void
  initializeInternalVars([[maybe_unused]] ParticleSubset* pset,
                         [[maybe_unused]] DataWarehouse* new_dw) {};

  virtual void
  getInternalVars([[maybe_unused]] ParticleSubset* pset,
                  [[maybe_unused]] DataWarehouse* old_dw) {};

  virtual void
  allocateAndPutInternalVars([[maybe_unused]] ParticleSubset* pset,
                             [[maybe_unused]] DataWarehouse* new_dw) {};

  virtual void
  allocateAndPutRigid([[maybe_unused]] ParticleSubset* pset,
                      [[maybe_unused]] DataWarehouse* new_dw) {};

  //__________________________________
  //  where the work is done
  virtual void
  computeDeviatoricStressInc([[maybe_unused]] const particleIndex idx,
                             [[maybe_unused]] const ModelStateBase* plaState,
                             [[maybe_unused]] DeformationState* defState,
                             [[maybe_unused]] const double delT) {};

  virtual void
  updateInternalStresses([[maybe_unused]] const particleIndex idx,
                         [[maybe_unused]] const Matrix3&,
                         [[maybe_unused]] DeformationState* defState,
                         [[maybe_unused]] const double delT) {};

  virtual void
  rotateInternalStresses([[maybe_unused]] const particleIndex idx,
                         [[maybe_unused]] const Matrix3&) {};
};
} // End namespace Uintah

#endif // __DEVSTRESSMODEL_H__
