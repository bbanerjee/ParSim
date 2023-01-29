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

// SpecifiedBody.h

#ifndef __SPECIFIED_BODY_H_
#define __SPECIFIED_BODY_H_

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
#include <Core/Parallel/UintahParallelComponent.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/ProblemSpec/ProblemSpecP.h>
#include <string>

namespace Uintah {

class VarLabel;

/**************************************

CLASS
 SpecifiedBodyContact

 Short description...

GENERAL INFORMATION

 SpecifiedBodyContact.h

 Andrew Brydon
 andrew@lanl.gov

 based on RigidBodyContact.
   Jim Guilkey
   Department of Mechanical Engineering
   University of Utah

   Center for the Simulation of Accidental Fires and Explosions (C-SAFE)


KEYWORDS
 Contact_Model_specified_velocity

DESCRIPTION
One of the derived Contact classes.  Allow motion of a body
to be specified by an input file.

the format of the input is
<contact>
  <type>specified</type>
  <filename>fname.txt</filename>
  <direction>[1,1,1]</direction>
  <material>0</material>
</contact>

where filename points to an existing test file (which much exist from all
mpi processors), storing a list of text columns
   {simtime}  {uvel} {vvel} {wvel}

the times must be given in ascending order.
linear interpolation is performed on values, and the end values are used if out
of the time range.

the direction can be used to limit the directions which the rigid region affects
the velocities of the other material. This can be used to impose a normal
velocity
with slip in the other directions. The default of [1,1,1] applies sticky
contact.

the material is the rigid material to use, and is optional. default is 0.

There are two alternate formats (which exist for compatability with
RigidBodyContact)

<contact>
  <direction>[0,0,1]</direction>
</contact>

Apply center-of-mass velocity to objects in contact with material (default of 0
assumed).


<contact>
  <direction>[1,1,1]</direction>
  <stop_time>2.0</stop_time>
  <velocity_after_stop>[0,1,0]</velocity_after_stop>
<contact>

when t>stop_time, impose velocity_after_stop.
This can be combined with either rigid velocity (as shown) or a velocity profile
to
produce wild and wacky results.

the velocity_after_stop is optional, with default (0,0,0).


Notes: (Probably out of date - Jim ??)

   Contact conditions are specified in two stages, and deal with 4 velocity
fields

   v^k        velocity at start of exMomInterpolated
   v*^k       velocity coming out of exMomInterpolated (with rigid cells set)
   v^k+1      velocity coming in to exMomIntegrated
   v*^k+1     velocity coming out of exMomIntegrated (with rigid cells set)

****************************************/

class SpecifiedBodyContact : public Contact
{
public:
  SpecifiedBodyContact(const ProcessorGroup* myworld,
                       const MaterialManagerP& d_mat_manager,
                       const MPMLabel* lb,
                       const MPMFlags* flag,
                       ProblemSpecP& ps);

  SpecifiedBodyContact(const SpecifiedBodyContact& con) = delete;
  SpecifiedBodyContact(SpecifiedBodyContact&& con)      = delete;
  SpecifiedBodyContact&
  operator=(const SpecifiedBodyContact& con) = delete;
  SpecifiedBodyContact&
  operator=(SpecifiedBodyContact&& con) = delete;

  virtual ~SpecifiedBodyContact() = default;

  // Currently, setting if any materials are rigid
  virtual void
  setContactMaterialAttributes();

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

private:
  double d_stop_time{ std::numeric_limits<double>::max() };
  double d_vol_const{ 0.0 };
  Vector d_vel_after_stop{ 0.0, 0.0, 0.0 };
  int d_material{ 0 };
  bool d_normal_only{ false };
  bool d_include_rotation{ false };
  int d_rotation_axis{ -99 };
  bool d_rigid_velocity{ true };
  int d_exclude_material{ -999 };
  std::string d_filename{ "none" };
  IntVector d_direction{ 0, 0, 1 };
  std::vector<std::pair<double, Vector>> d_vel_profile;
  std::vector<std::pair<double, Vector>> d_rot_profile;
  std::vector<std::pair<double, Vector>> d_ori_profile;

  struct ImposedData
  {
    Vector velocity{ 0.0, 0.0, 0.0 };
    Vector omega{ 0.0, 0.0, 0.0 };
    Vector origin{ 0.0, 0.0, 0.0 };
  };

  struct ReactionData
  {
    std::map<int, Vector> force;
    std::map<int, Vector> torque;

    ReactionData();
    ~ReactionData() = default;
  };

  struct TransmittedData
  {
    std::map<int, Vector> force;
    Vector all_material_force{ 0.0, 0.0, 0.0 };

    TransmittedData();
    ~TransmittedData() = default;
  };

  void
  readSpecifiedVelocityFile();

  void
  writeSpecifiedVelocityFile();

  Vector
  findVelFromProfile(double t) const;

  Vector
  findValueFromProfile(
    double t,
    const std::vector<std::pair<double, Vector>>& profile) const;

  void
  computeNormalBasedExchange(const Patch* patch,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             constNCdoubleArray& gMass,
                             constNCdoubleArray& gVolume,
                             constNCVectorArray& gInternalForce,
                             const ImposedData& imposed,
                             double delT,
                             NCVectorArray& gVelocity_star,
                             ReactionData& reaction,
                             TransmittedData& transmitted);

  std::pair<Vector, Vector>
  getRotationComponent(const Patch* patch,
                       const IntVector& node,
                       const ImposedData& imposed,
                       double delT);

  void
  computeDirectionBasedExchange(const Patch* patch,
                                DataWarehouse* old_dw,
                                constNCdoubleArray& gMass,
                                constNCdoubleArray& gVolume,
                                constNCVectorArray& gInternalForce,
                                const ImposedData& imposed,
                                double delT,
                                NCVectorArray& gVelocity_star,
                                ReactionData& reaction,
                                TransmittedData& transmitted);
};

} // end namespace Uintah

#endif /* __SPECIFIED_BODY_H_ */
