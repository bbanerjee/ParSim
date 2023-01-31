/*
 * The MIT License
 *
 * Copyright (c) 1997-2022 The University of Utah
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

// FrictionContactLR.h

#ifndef __FRICTIONLR_H__
#define __FRICTIONLR_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {
/**************************************

CLASS
   FrictionContactLR

   This version of contact is based on John Nairn and Chad
   Hammerquist's 2019 manuscript that describes the use of logistic
   regression to find a common normal between objects, and uses
   particle geometry to find the most prominent portion of a particle
   at each node, and applies contact when contacting materials'
   prominences overlap.

GENERAL INFORMATION

   FrictionContactLR.h

   Steven G. Parker
   Department of Computer Science
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
   Contact_Model_Friction

DESCRIPTION
  One of the derived Contact classes.  This particular
  version is used to apply Coulombic frictional contact.

WARNING

****************************************/

class FrictionContactLR : public Contact
{

public:
  FrictionContactLR(const ProcessorGroup* myworld,
                    const MaterialManagerP& mat_manager,
                    const MPMLabel* labels,
                    const MPMFlags* flags,
                    ProblemSpecP& ps);

  FrictionContactLR(const FrictionContactLR& con) = delete;
  FrictionContactLR&
  operator=(const FrictionContactLR& con) = delete;

  virtual ~FrictionContactLR() = default;

  virtual void
  setContactMaterialAttributes() override;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) override;

  void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;

  void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

private:
  // Coefficient of friction
  double d_mu{ 0.0 };

  // Nodal volume fraction that must occur before contact is applied
  double d_vol_const{ 0.0 };
  double d_sep_fac{ 1.0e200 };
  bool d_rigid_material{ false };

  void
  exMomInterpolated(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* matls,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw);

  void
  exMomIntegrated(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* matls,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw);
};
} // End namespace Uintah

#endif /* __FRICTIONLR_H__ */
