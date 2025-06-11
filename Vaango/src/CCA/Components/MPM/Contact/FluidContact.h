/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

// Fluid.h

#ifndef __CCA_COMPONENTS_MPM_CONTACT_FLUIDCONTACT_H__
#define __CCA_COMPONENTS_MPM_CONTACT_FLUIDCONTACT_H__

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>

#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

namespace Uintah {

class DataWarehouse;
class MPMLabel;
class HydroMPMLabel;
class MPMFlags;
class ProcessorGroup;
class Patch;
class VarLabel;
class Task;

/**************************************

CLASS
   FluidContact

   Short description...

GENERAL INFORMATION

   FluidContact.h

   Hilde Aas Nøst
   Department of Civil Engineering
   Norwegian University of Science and Technology


KEYWORDS
   Contact_Model_Fluid

DESCRIPTION
  One of the derived Contact classes.  This particular
  version is used to apply contact with pore fluid.

WARNING

****************************************/

class FluidContact : public Contact
{
public:
  // Constructor
  FluidContact(const ProcessorGroup* myworld,
               const MaterialManagerP& mat_manager,
               const MPMLabel* labels,
               const MPMFlags* flags,
               ProblemSpecP& ps);

  // Destructor
  virtual ~FluidContact() = default;

  // Prevent copying/move of this class
  FluidContact(const FluidContact& con) = delete;
  FluidContact(FluidContact&& con)      = delete;
  FluidContact&
  operator=(const FluidContact& con) = delete;
  FluidContact&
  operator=(FluidContact&& con) = delete;

  virtual void
  setContactMaterialAttributes() override;

  virtual void
  outputProblemSpec(ProblemSpecP& ps) override;

  // Basic contact methods
  virtual void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

  virtual void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;

private:
  int d_rigid_material{ 0 };
  double d_vol_const{ 0.0 };
  double d_sep_fac{ 1.0e100 };
  bool d_comp_collinear_norms{ true };

  std::unique_ptr<HydroMPMLabel> d_hydro_mpm_labels;
};
} // End namespace Uintah

#endif /* __CCA_COMPONENTS_MPM_CONTACT_FLUIDCONTACT_H__ */
