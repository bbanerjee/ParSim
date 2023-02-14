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

// SpecifiedBodyContact.cc
#include <CCA/Components/MPM/Contact/SpecifiedBodyContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

using std::cerr;

using namespace Uintah;

SpecifiedBodyContact::SpecifiedBodyContact(const ProcessorGroup* myworld,
                                           const MaterialManagerP& mat_manager,
                                           const MPMLabel* labels,
                                           const MPMFlags* flags,
                                           ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  d_need_normals = true;

  IntVector defaultDir(0, 0, 1);
  ps->getWithDefault("direction", d_direction, defaultDir);

  ps->getWithDefault("master_material", d_material, 0);
  d_matls.add(d_material); // always need specified material

  ps->get("volume_constraint", d_vol_const);

  ps->getWithDefault("normal_only", d_normal_only, false);

  ps->getWithDefault("include_rotation", d_include_rotation, false);

  // read a list of values from a file
  ps->get("filename", d_filename);
  readSpecifiedVelocityFile();

  // disable all changes after this time
  ps->getWithDefault(
    "stop_time", d_stop_time, std::numeric_limits<double>::max());
  ps->getWithDefault("velocity_after_stop", d_vel_after_stop, Vector(0, 0, 0));
}

void
SpecifiedBodyContact::readSpecifiedVelocityFile()
{
  if (d_filename != "") {
    std::ifstream is(d_filename.c_str());
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Could not open MPM specified contact motion file "
          << d_filename << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
    double t0(-1.e9);
    if (d_include_rotation) {
      while (is) {
        double t1;
        double vx, vy, vz, ox, oy, oz, wx, wy, wz;
        is >> t1 >> vx >> vy >> vz >> ox >> oy >> oz >> wx >> wy >> wz;
        if (is) {
          if (t1 <= t0) {
            std::ostringstream err;
            err << "**ERROR** Time in specified contact profile file "
                << "is not monotomically increasing";
            throw ProblemSetupException(err.str(), __FILE__, __LINE__);
          }
          d_vel_profile.push_back(
            std::pair<double, Vector>(t1, Vector(vx, vy, vz)));
          d_rot_profile.push_back(
            std::pair<double, Vector>(t1, Vector(wx, wy, wz)));
          d_ori_profile.push_back(
            std::pair<double, Vector>(t1, Vector(ox, oy, oz)));
        }
        t0 = t1;
      }
    } else {
      while (is) {
        double t1;
        double vx, vy, vz;
        is >> t1 >> vx >> vy >> vz;
        if (is) {
          if (t1 <= t0) {
            std::ostringstream err;
            err << "**ERROR** Time in specified contact profile file "
                << "is not monotomically increasing";
            throw ProblemSetupException(err.str(), __FILE__, __LINE__);
          }
          d_vel_profile.push_back(
            std::pair<double, Vector>(t1, Vector(vx, vy, vz)));
        }
        t0 = t1;
      }
    }
    if (d_vel_profile.size() < 2) {
      std::ostringstream err;
      err << "**ERROR** Specified contact: failed to generate valid "
          << "velocity profile.";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }
  return;
}

void
SpecifiedBodyContact::setContactMaterialAttributes()
{
  MPMMaterial* mpm_matl =
    static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", d_material));
  mpm_matl->setIsRigid(true);
}

void
SpecifiedBodyContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "specified");
  contact_ps->appendElement("filename", d_filename);
  contact_ps->appendElement("direction", d_direction);
  contact_ps->appendElement("master_material", d_material);
  contact_ps->appendElement("stop_time", d_stop_time);
  contact_ps->appendElement("velocity_after_stop", d_vel_after_stop);
  contact_ps->appendElement("volume_constraint", d_vol_const);
  contact_ps->appendElement("include_rotation", d_include_rotation);
  contact_ps->appendElement("normal_only", d_normal_only);
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  contact_ps->appendElement("exclude_material", d_exclude_material);

  d_matls.outputProblemSpec(contact_ps);

  writeSpecifiedVelocityFile();
}

void
SpecifiedBodyContact::writeSpecifiedVelocityFile()
{
  if (d_filename != "") {
    std::string udaDir = d_mpm_flags->d_output->getOutputLocation();

    //  Bulletproofing
    Uintah::Dir dir(udaDir);
    if (!dir.exists()) {
      std::ostringstream warn;
      warn
        << "ERROR:SpecifiedBodyContact The main uda directory does not exist.";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    std::ostringstream fname;
    fname << udaDir << "/" << d_filename;
    std::string filename = fname.str();

    std::ofstream fp(filename.c_str());

    int smax = (int)(d_vel_profile.size());

    if (d_include_rotation) {
      for (int i = 0; i < smax; i++) {
        fp << d_vel_profile[i].first << " " << d_vel_profile[i].second.x()
           << " " << d_vel_profile[i].second.y() << " "
           << d_vel_profile[i].second.z() << " " << d_ori_profile[i].second.x()
           << " " << d_ori_profile[i].second.y() << " "
           << d_ori_profile[i].second.z() << " " << d_rot_profile[i].second.x()
           << " " << d_rot_profile[i].second.y() << " "
           << d_rot_profile[i].second.z() << endl;
      }
    } else {
      for (int i = 0; i < smax; i++) {
        fp << d_vel_profile[i].first << " " << d_vel_profile[i].second.x()
           << " " << d_vel_profile[i].second.y() << " "
           << d_vel_profile[i].second.z() << endl;
      }
    }
  }
}

void
SpecifiedBodyContact::addComputesAndRequires(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls,
                                             const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  Task* t = scinew Task("SpecifiedBodyContact::exchangeMomentum",
                        this,
                        &SpecifiedBodyContact::exchangeMomentum,
                        gVelocity_label);

  MaterialSubset* zero_matl = scinew MaterialSubset();
  zero_matl->add(0);
  zero_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gSurfNormLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gInternalForceLabel, Ghost::None);
  t->requires(
    Task::OldDW, d_mpm_labels->NC_CCweightLabel, zero_matl, Ghost::None);

  t->modifies(gVelocity_label, mss);

  //  Create reductionMatlSubSet that includes all mss matls
  //  and the global matlsubset
  const MaterialSubset* global_mss = t->getGlobalMatlSubset();

  MaterialSubset* reduction_mss = scinew MaterialSubset();
  reduction_mss->addReference();
  reduction_mss->add(global_mss->get(0));

  size_t numMatls = mss->size();
  if (numMatls > 1) { // ignore for single matl problems
    for (size_t m = 0; m < numMatls; m++) {
      reduction_mss->add(mss->get(m));
    }
  }

  t->computes(d_mpm_labels->RigidReactionForceLabel, reduction_mss);
  t->computes(d_mpm_labels->RigidReactionTorqueLabel, reduction_mss);

  if (d_mpm_flags->d_reductionVars->sumTransmittedForce) {
    t->computes(
      d_mpm_labels->SumTransmittedForceLabel, reduction_mss, Task::OutOfDomain);
  }

  sched->addTask(t, patches, matls);

  if (zero_matl && zero_matl->removeReference()) {
    delete zero_matl;
  }
  if (reduction_mss && reduction_mss->removeReference()) {
    delete reduction_mss;
  };
}

void
SpecifiedBodyContact::exchangeMomentum(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* matls,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw,
                                       const VarLabel* gVelocity_label)
{
  Ghost::GhostType gnone = Ghost::None;

  int numMatls = d_mat_manager->getNumMaterials("MPM");

  simTime_vartype simTime;
  old_dw->get(simTime, d_mpm_labels->simulationTimeLabel);

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  // rigid_velocity just means that the master_material's initial velocity
  // remains constant through the simulation, until d_stop_time is reached.
  // If the velocity comes from a profile specified in a file, or after
  // d_stop_time, rigid_velocity is false
  ImposedData imposed;
  const double tcurr = simTime;
  if (tcurr > d_stop_time) {
    d_rigid_velocity = false;
    imposed.velocity = d_vel_after_stop;
  } else if (d_vel_profile.size() > 0) {
    d_rigid_velocity = false;
    imposed.velocity = findVelFromProfile(tcurr);
    imposed.omega    = findValueFromProfile(tcurr, d_rot_profile);
    imposed.origin   = findValueFromProfile(tcurr, d_ori_profile);
  }

  // If rotation axis is aligned with a ordinal direction,
  // use the exact treatment, otherwise default to the approximate
  if (d_include_rotation) {
    double ROL = imposed.omega.length();
    if (std::abs(Dot(imposed.omega / ROL, Vector(1., 0., 0.))) > 0.99) {
      d_rotation_axis = 0;
    } else if (std::abs(Dot(imposed.omega / ROL, Vector(0., 1., 0.))) > 0.99) {
      d_rotation_axis = 1;
    } else if (std::abs(Dot(imposed.omega / ROL, Vector(0., 0., 1.))) > 0.99) {
      d_rotation_axis = 2;
    }
  }

  // Initialize reaction forces and torques for each material
  ReactionData reaction;

  // get material id of master material
  int matID_master =
    d_mat_manager->getMaterial("MPM", d_material)->getDWIndex();

  // Set up loop
  constNCdoubleArray gMass(numMatls), gVolume(numMatls);
  constNCVectorArray gVelocity(numMatls), gInternalForce(numMatls);
  NCVectorArray gVelocity_star(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    TransmittedData transmitted;

    for (int m = 0; m < matls->size(); m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->get(gInternalForce[m],
                  d_mpm_labels->gInternalForceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, dwi, patch);
    }

    if (d_normal_only) {
      computeNormalBasedExchange(patch,
                                 old_dw,
                                 new_dw,
                                 gMass,
                                 gVolume,
                                 gInternalForce,
                                 imposed,
                                 delT,
                                 gVelocity_star,
                                 reaction,
                                 transmitted);
    } else {
      computeDirectionBasedExchange(patch,
                                    old_dw,
                                    gMass,
                                    gVolume,
                                    gInternalForce,
                                    imposed,
                                    delT,
                                    gVelocity_star,
                                    reaction,
                                    transmitted);
    }

    // Put the sumTransmittedForce contribution into the reduction variables
    if (d_mpm_flags->d_reductionVars->sumTransmittedForce) {
      new_dw->put(sumvec_vartype(transmitted.all_material_force),
                  d_mpm_labels->SumTransmittedForceLabel,
                  nullptr,
                  -1);
      new_dw->put_sum_vartype(
        transmitted.force, d_mpm_labels->SumTransmittedForceLabel, matls);
    }

  } // end patch for

  // Collect reactions from all patches and materials
  for (int n = 0; n < numMatls; n++) {
    if (n != d_material) {
      int matID = d_mat_manager->getMaterial("MPM", n)->getDWIndex();
      reaction.force[matID_master] += reaction.force[matID];
      reaction.torque[matID_master] += reaction.torque[matID];
    }
  }

  for (int n = 0; n < numMatls; n++) {
    if (numMatls > 1) { // ignore for single matl problems
      int matID = d_mat_manager->getMaterial("MPM", n)->getDWIndex();
      new_dw->put(sumvec_vartype(reaction.force[matID]),
                  d_mpm_labels->RigidReactionForceLabel,
                  nullptr,
                  matID);
      new_dw->put(sumvec_vartype(reaction.torque[matID]),
                  d_mpm_labels->RigidReactionTorqueLabel,
                  nullptr,
                  matID);
    }
  }

  new_dw->put(sumvec_vartype(reaction.force[matID_master]),
              d_mpm_labels->RigidReactionForceLabel,
              nullptr,
              -1);
  new_dw->put(sumvec_vartype(reaction.torque[matID_master]),
              d_mpm_labels->RigidReactionTorqueLabel,
              nullptr,
              -1);
}

void
SpecifiedBodyContact::computeNormalBasedExchange(
  const Patch* patch,
  DataWarehouse* old_dw,
  DataWarehouse* new_dw,
  constNCdoubleArray& gMass,
  constNCdoubleArray& gVolume,
  constNCVectorArray& gInternalForce,
  const ImposedData& imposed,
  double delT,
  NCVectorArray& gVelocity_star,
  ReactionData& reaction,
  TransmittedData& transmitted)
{
  Ghost::GhostType gnone = Ghost::None;
  int numMatls           = d_mat_manager->getNumMaterials("MPM");

  Vector dx       = patch->dCell();
  double cell_vol = dx.x() * dx.y() * dx.z();

  constNCdouble NC_CCweight;
  old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

  constNCVector gSurfNorm;
  new_dw->get(
    gSurfNorm, d_mpm_labels->gSurfNormLabel, d_material, patch, gnone, 0);

  int matID_master =
    d_mat_manager->getMaterial("MPM", d_material)->getDWIndex();

  for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
    IntVector node = *iter;

    // Determine nodal volume
    double totalNodalVol = 0.0;
    for (int mat = 0; mat < numMatls; mat++) {
      totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
    }

    // Compute the master velocity
    auto rotation_part = getRotationComponent(patch, node, imposed, delT);
    Vector master_vel  = rotation_part.second + imposed.velocity;
    if (d_rigid_velocity) {
      master_vel = gVelocity_star[d_material][node];
    }

    double exclude_mass = 0.0;
    if (d_exclude_material >= 0) {
      exclude_mass = gMass[d_exclude_material][node];
    }

    for (int mat = 0; mat < numMatls; mat++) {

      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();

      if (!d_matls.requested(mat) || exclude_mass >= 1.e-99) {
        continue;
      }

      Vector new_vel        = gVelocity_star[mat][node];
      Vector normal         = gSurfNorm[node];
      Vector deltaV         = new_vel - master_vel;
      double normalDeltaVel = Dot(normal, deltaV);
      if (normalDeltaVel < 0.0) {
        Vector normal_normaldV = normal * normalDeltaVel;
        new_vel                = gVelocity_star[mat][node] - normal_normaldV;
      }

      if (!compare(gMass[d_material][node], 0.) &&
          (totalNodalVol / cell_vol) > d_vol_const) {
        Vector old_vel          = gVelocity_star[mat][node];
        Vector transmittedForce = gMass[mat][node] * (new_vel - old_vel) / delT;

        reaction.force[matID] -= gInternalForce[mat][node];
        reaction.torque[matID] -= Cross(rotation_part.first, transmittedForce);
        transmitted.force[matID_master] -= transmittedForce;
        transmitted.all_material_force -= transmittedForce;

        gVelocity_star[mat][node] = new_vel;
      }

      // std::cout << "After rigid contact: Node = " << c << " material = " << n
      //           << " gVel = " << gVelocity_star[n][c] << "\n";
    } // end for matls
  }   // end for NodeIterator
}

std::pair<Vector, Vector>
SpecifiedBodyContact::getRotationComponent(const Patch* patch,
                                           const IntVector& node,
                                           const ImposedData& imposed,
                                           double delT)
{
  Vector rotation_part(0.0, 0.0, 0.0);
  Vector r_vec(0.0, 0.0, 0.0);

  if (d_include_rotation) {
    Point nodePos     = patch->getNodePosition(node);
    Point nodePos_new = nodePos;

    // vector from node to a point on the axis of rotation
    r_vec = nodePos - imposed.origin.asPoint();
    if (d_rotation_axis == 0) { // rotation about x-axis
      double posz      = nodePos.z() - imposed.origin.z();
      double posy      = nodePos.y() - imposed.origin.y();
      double theta     = std::atan2(posz, posy);
      double thetaPlus = theta + imposed.omega[0] * delT;
      double R         = std::sqrt(posy * posy + posz * posz);
      nodePos_new      = Point(nodePos.x(),
                          R * std::cos(thetaPlus) + imposed.origin.y(),
                          R * std::sin(thetaPlus) + imposed.origin.z());
    } else if (d_rotation_axis == 1) { // rotation about y-axis
      double posx      = nodePos.x() - imposed.origin.x();
      double posz      = nodePos.z() - imposed.origin.z();
      double theta     = std::atan2(posx, posz);
      double thetaPlus = theta + imposed.omega[1] * delT;
      double R         = std::sqrt(posz * posz + posx * posx);
      nodePos_new      = Point(R * std::sin(thetaPlus) + imposed.origin.x(),
                          nodePos.y(),
                          R * std::cos(thetaPlus) + imposed.origin.z());
    } else if (d_rotation_axis == 2) { // rotation about z-axis
      double posx      = nodePos.x() - imposed.origin.x();
      double posy      = nodePos.y() - imposed.origin.y();
      double theta     = std::atan2(posy, posx);
      double thetaPlus = theta + imposed.omega[2] * delT;
      double R         = std::sqrt(posx * posx + posy * posy);
      nodePos_new      = Point(R * std::cos(thetaPlus) + imposed.origin.x(),
                          R * std::sin(thetaPlus) + imposed.origin.y(),
                          nodePos.z());
    }
    rotation_part = (nodePos_new - nodePos) / delT;
    if (d_rotation_axis == -99) {
      // normal vector from the axis of rotation to the node
      // Vector axis_norm=requested_omega/(requested_omega.length()+1.e-100);
      // Vector rad = r - Dot(r,axis_norm)*axis_norm;
      rotation_part = Cross(imposed.omega, r_vec);
    }
  }

  return std::make_pair(r_vec, rotation_part);
}

void
SpecifiedBodyContact::computeDirectionBasedExchange(
  const Patch* patch,
  DataWarehouse* old_dw,
  constNCdoubleArray& gMass,
  constNCdoubleArray& gVolume,
  constNCVectorArray& gInternalForce,
  const ImposedData& imposed,
  double delT,
  NCVectorArray& gVelocity_star,
  ReactionData& reaction,
  TransmittedData& transmitted)
{
  Ghost::GhostType gnone = Ghost::None;
  int numMatls           = d_mat_manager->getNumMaterials("MPM");

  Vector dx       = patch->dCell();
  double cell_vol = dx.x() * dx.y() * dx.z();

  constNCdouble NC_CCweight;
  old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

  int matID_master =
    d_mat_manager->getMaterial("MPM", d_material)->getDWIndex();

  for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
    IntVector node = *iter;

    // Determine nodal volume
    double totalNodalVol = 0.0;
    for (int mat = 0; mat < numMatls; mat++) {
      totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
    }

    // Compute vector from node to a point on the axis of rotation
    Vector rVec{ 0.0, 0.0, 0.0 };
    if (d_include_rotation) {
      rVec = patch->getNodePosition(node) - imposed.origin.asPoint();
    }

    // also updates material d_material to new velocity.
    for (int mat = 0; mat < numMatls; mat++) {

      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();

      Vector master_vel = imposed.velocity;

      Vector new_vel = gVelocity_star[mat][node];
      if (mat == d_material || d_direction[0]) {
        new_vel.x(master_vel.x());
      }
      if (mat == d_material || d_direction[1]) {
        new_vel.y(master_vel.y());
      }
      if (mat == d_material || d_direction[2]) {
        new_vel.z(master_vel.z());
      }

      if (!compare(gMass[d_material][node], 0.) &&
          (totalNodalVol / cell_vol) > d_vol_const) {
        Vector old_vel          = gVelocity_star[mat][node];
        Vector transmittedForce = gMass[mat][node] * (new_vel - old_vel) / delT;

        reaction.force[matID] -= gInternalForce[mat][node];
        reaction.torque[matID] -= Cross(rVec, transmittedForce);
        transmitted.force[matID_master] -= transmittedForce;
        transmitted.all_material_force -= transmittedForce;

        gVelocity_star[mat][node] = new_vel;
      }

      // std::cout << "After rigid contact: Node = " << c << " material = " << n
      //           << " gVel = " << gVelocity_star[n][c] << "\n";
    } // end for matls
  }   // end for NodeIterator
}

// find velocity from table of values
Vector
SpecifiedBodyContact::findVelFromProfile(double t) const
{
  return findValueFromProfile(t, d_vel_profile);
}

// find value from table of values
Vector
SpecifiedBodyContact::findValueFromProfile(
  double t,
  const std::vector<std::pair<double, Vector>>& profile) const
{
  auto iter = std::find_if(
    profile.begin(), profile.end(), [t](const std::pair<double, Vector>& data) {
      return data.first > t;
    });
  if (iter == profile.begin()) {
    return iter->second;
  } else if (iter == profile.end()) {
    return (iter - 1)->second;
  } else {
    double t_val = (iter->first - t) / (iter->first - (iter - 1)->first);
    Vector vel   = (1.0 - t_val) * (iter - 1)->second + t_val * iter->second;
    return vel;
  }
}

SpecifiedBodyContact::ReactionData::ReactionData()
{
  force  = MPMCommon::initializeMap(Vector(0.));
  torque = MPMCommon::initializeMap(Vector(0.));
}

SpecifiedBodyContact::TransmittedData::TransmittedData()
{
  force = MPMCommon::initializeMap(Vector(0.));
}
