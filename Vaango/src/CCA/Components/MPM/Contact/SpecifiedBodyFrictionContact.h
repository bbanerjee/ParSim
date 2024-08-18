/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2015-2024 Biswajit Banerjee
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

// SpecifiedBodyFrictionContact.h

#ifndef __SPECIFIED_BODY_FRICTION_H_
#define __SPECIFIED_BODY_FRICTION_H_

#include <CCA/Components/MPM/Contact/SpecifiedBodyContact.h>

#include <CCA/Ports/DataWarehouseP.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/LevelP.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>

namespace Uintah {

class VarLabel;
class Output;

/**************************************

CLASS
 SpecifiedBodyFrictionContact

 Short description...

GENERAL INFORMATION

 SpecifiedBodyFrictionContact.h

 Jim Guilkey
 Department of Mechanical Engineering
 University of Utah

 Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
 Contact_Model_specified_velocity

DESCRIPTION
One of the derived Contact classes.  Allow motion of a body
to be specified by an input file, while allowing for frictional
sliding with other materials.

the format of the input is
<contact>
  <type>specified_friction</type>
  <filename>fname.txt</filename>
  <master_material>0</master_material>
  <mu>0.1</mu>
</contact>

where filename points to an existing test file (which much exist from all
mpi processors), storing a list of text columns
   {simtime}  {uvel} {vvel} {wvel}

the times must be given in ascending order.
linear interpolation is performed on values, and the end values are used if out
of the time range.

the direction can be used to limit the directions which the rigid region affects
the velocities of the other material. This can be used to impose a normal
velocity with slip in the other directions. The default of [1,1,1] applies
sticky contact.

the material is the rigid material to use, and is optional. default is 0.

There are two alternate formats (which exist for compatability with
RigidBodyContact)

when t>stop_time, impose velocity_after_stop.
This can be combined with either rigid velocity (as shown) or a velocity profile
to produce wild and wacky results.

the velocity_after_stop is optional, with default (0,0,0).

****************************************/

class SpecifiedBodyFrictionContact : public SpecifiedBodyContact
{
private:
  double d_mu;

public:
  // Constructor
  SpecifiedBodyFrictionContact(const ProcessorGroup* myworld,
                               const MaterialManagerP& mat_manager,
                               const MPMLabel* labels,
                               const MPMFlags* flags,
                               ProblemSpecP& ps);

  // Destructor
  virtual ~SpecifiedBodyFrictionContact() = default;

  // Prevent copying/move of this class
  SpecifiedBodyFrictionContact(const SpecifiedBodyFrictionContact& con) =
    delete;
  SpecifiedBodyFrictionContact(SpecifiedBodyFrictionContact&& con) = delete;
  SpecifiedBodyFrictionContact&
  operator=(const SpecifiedBodyFrictionContact& con) = delete;
  SpecifiedBodyFrictionContact&
  operator=(SpecifiedBodyFrictionContact&& con) = delete;

  virtual void
  outputProblemSpec(ProblemSpecP& ps);

  // Currently, setting if any materials are rigid
  virtual void
  setContactMaterialAttributes();

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

};

} // end namespace Uintah

#endif /* __SPECIFIED_BODY_FRICTION_H_ */
