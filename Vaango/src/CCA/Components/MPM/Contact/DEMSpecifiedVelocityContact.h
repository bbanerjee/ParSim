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

#ifndef __DEM_SPECIFIED_VELOCITY_CONTACT_H_
#define __DEM_SPECIFIED_VELOCITY_CONTACT_H_

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactMaterialSpec.h>
#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/MPM/Core/MPMUtils.h>
#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>
#include <vector>

namespace Uintah {

class VarLabel;

class DEMSpecifiedVelocityContact : public Contact
{
public:
  DEMSpecifiedVelocityContact(const ProcessorGroup* myworld,
                              const MaterialManagerP& d_mat_manager,
                              const MPMLabel* lb,
                              const MPMFlags* flag,
                              ProblemSpecP& ps);

  DEMSpecifiedVelocityContact(const DEMSpecifiedVelocityContact& con) = delete;
  DEMSpecifiedVelocityContact(DEMSpecifiedVelocityContact&& con)      = delete;
  DEMSpecifiedVelocityContact&
  operator=(const DEMSpecifiedVelocityContact& con) = delete;
  DEMSpecifiedVelocityContact&
  operator=(DEMSpecifiedVelocityContact&& con) = delete;

  virtual ~DEMSpecifiedVelocityContact() = default;

  virtual void
  setContactMaterialAttributes() override;

  void
  outputProblemSpec(ProblemSpecP& ps) override;

  void
  exchangeMomentum(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* matls,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw,
                   const VarLabel* label) override;

  void
  addComputesAndRequires(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* matls,
                         const VarLabel* label) override;

protected:
  int d_material{ 0 };
  std::string d_filename{ "" };
  std::vector<std::pair<double, Vector>> d_vel_profile;

  void
  readSpecifiedVelocityFile();

  Vector
  findVelFromProfile(double t) const;
};

} // end namespace Uintah

#endif /* __DEM_SPECIFIED_VELOCITY_CONTACT_H_ */
