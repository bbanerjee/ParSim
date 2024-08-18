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

// SpecifiedBodyFrictionContact.cc
#include <CCA/Components/MPM/Contact/SpecifiedBodyFrictionContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMCommon.h>
#include <CCA/Components/MPM/Core/MPMFlags.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Output.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <dirent.h>
#include <fstream>
#include <limits>
#include <vector>

using namespace Uintah;
using std::cerr;

SpecifiedBodyFrictionContact::SpecifiedBodyFrictionContact(
  const ProcessorGroup* myworld,
  const MaterialManagerP& mat_manager,
  const MPMLabel* labels,
  const MPMFlags* flags,
  ProblemSpecP& ps)
  : SpecifiedBodyContact(myworld, mat_manager, labels, flags, ps)
{
  // Constructor
  // read a list of values from a file
  ps->get("filename", d_filename);
  ps->require("mu", d_mu);

  ps->getWithDefault("master_material", d_material, 0);
  d_matls.add(d_material); // always need specified material

  ps->getWithDefault("include_rotation", d_include_rotation, false);
  ps->getWithDefault("exclude_material", d_exclude_material, -999);

  readSpecifiedVelocityFile();

  // disable all changes after this time
  ps->getWithDefault(
    "stop_time", d_stop_time, std::numeric_limits<double>::max());
  ps->getWithDefault("velocity_after_stop", d_vel_after_stop, Vector(0, 0, 0));
}

void
SpecifiedBodyFrictionContact::setContactMaterialAttributes()
{
  SpecifiedBodyContact::setContactMaterialAttributes();
}

void
SpecifiedBodyFrictionContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "specified_friction");
  contact_ps->appendElement("filename", d_filename);
  contact_ps->appendElement("master_material", d_material);
  contact_ps->appendElement("stop_time", d_stop_time);
  contact_ps->appendElement("velocity_after_stop", d_vel_after_stop);
  contact_ps->appendElement("include_rotation", d_include_rotation);
  contact_ps->appendElement("mu", d_mu);
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  contact_ps->appendElement("exclude_material", d_exclude_material);

  d_matls.outputProblemSpec(contact_ps);

  writeSpecifiedVelocityFile();
}

void
SpecifiedBodyFrictionContact::addComputesAndRequires(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSet* matls,
  const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  Task* t = scinew Task("SpecifiedBodyFrictionContact::exchangeMomentum",
                        this,
                        &SpecifiedBodyFrictionContact::exchangeMomentum,
                        gVelocity_label);

  auto* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMatlProminenceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gAlphaMaterialLabel, Ghost::None);
  t->requires(Task::NewDW,
              d_mpm_labels->gNormAlphaToBetaLabel,
              z_matl,
              Ghost::None);
  t->requires(
    Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);

  t->modifies(gVelocity_label, mss);

  //  Create reductionMatlSubSet that includes all mss matls
  //  and the global matlsubset
  const MaterialSubset* global_mss = t->getGlobalMatlSubset();

  auto* reduction_mss = scinew MaterialSubset();
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
    t->computes(d_mpm_labels->SumTransmittedForceLabel,
                reduction_mss,
                Task::OutOfDomain);
  }

  sched->addTask(t, patches, matls);

  if (z_matl && z_matl->removeReference()) {
    delete z_matl;
  }
  if (reduction_mss && reduction_mss->removeReference()) {
    delete reduction_mss;
  }
}

// apply boundary conditions to the interpolated velocity v^k+1
void
SpecifiedBodyFrictionContact::exchangeMomentum(const ProcessorGroup*,
                                               const PatchSubset* patches,
                                               const MaterialSubset* matls,
                                               DataWarehouse* old_dw,
                                               DataWarehouse* new_dw,
                                               const VarLabel* gVelocity_label)
{
  Ghost::GhostType gnone = Ghost::None;

  simTime_vartype simTime;
  old_dw->get(simTime, d_mpm_labels->simulationTimeLabel);

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  int numMatls = d_mat_manager->getNumMaterials("MPM");

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

  // Retrieve necessary data from DataWarehouse
  constNCdoubleArray gMass(numMatls);
  NCVectorArray gVelocity_star(numMatls);
  constNCVectorArray gInternalForce(numMatls);
  constNCdoubleArray gMatlProminence(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    TransmittedData transmitted;

    constNCVariable<double> NC_CCweight;
    constNCVariable<int> alphaMaterial;
    constNCVariable<Vector> normAlphaToBeta;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->get(
      alphaMaterial, d_mpm_labels->gAlphaMaterialLabel, 0, patch, gnone, 0);
    new_dw->get(
      normAlphaToBeta, d_mpm_labels->gNormAlphaToBetaLabel, 0, patch, gnone, 0);

    for (int m = 0; m < matls->size(); m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gInternalForce[m],
                  d_mpm_labels->gInternalForceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(gMatlProminence[m],
                  d_mpm_labels->gMatlProminenceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, dwi, patch);
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;

      auto rotation_part = getRotationComponent(patch, node, imposed, delT);
      Vector rigid_vel   = rotation_part.second + imposed.velocity;
      if (d_rigid_velocity) {
        rigid_vel = gVelocity_star[d_material][node];
      }

      double exclude_mass = 0.;
      if (d_exclude_material >= 0) {
        exclude_mass = gMass[d_exclude_material][node];
      }

      int alpha = alphaMaterial[node];
      if (alpha >= 0) { // Only work on nodes where alpha!=-99
        for (int n = 0; n < numMatls; n++) {
          int dwi = d_mat_manager->getMaterial("MPM", n)->getDWIndex();
          if (!d_matls.requested(n) || exclude_mass >= 1.e-99) {
            continue;
          }
          Vector new_vel = rigid_vel;

          if (n == d_material) {
            gVelocity_star[n][node] = new_vel;
          } else if (!compare(gMass[d_material][node], 0.) &&
                     !compare(gMass[n][node], 0)) {
            double separation =
              gMatlProminence[n][node] - gMatlProminence[alpha][node];
            if (separation <= 0.0) {
              Vector old_vel = gVelocity_star[n][node];

              Vector deltaVelocity  = gVelocity_star[n][node] - new_vel;
              Vector normal         = -1.0 * normAlphaToBeta[node];
              double normalDeltaVel = Dot(deltaVelocity, normal);
              Vector Dv(0., 0., 0.);
              if (normalDeltaVel > 0.0) {

                Vector normal_normaldV = normal * normalDeltaVel;
                Vector dV_normalDV     = deltaVelocity - normal_normaldV;
                Vector surfaceTangent =
                  dV_normalDV / (dV_normalDV.length() + 1.e-100);
                double tangentDeltaVelocity =
                  Dot(deltaVelocity, surfaceTangent);
                double frictionCoefficient =
                  Min(d_mu,
                      tangentDeltaVelocity / (fabs(normalDeltaVel) + 1.e-100));
                // Calculate velocity change needed to enforce contact
                Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                          fabs(normalDeltaVel);

                gVelocity_star[n][node] += Dv;

                reaction.force[dwi] -= gInternalForce[n][node];
                reaction.torque[dwi] +=
                  Cross(rotation_part.first,
                        gMass[n][node] * (new_vel - old_vel) / delT);
                transmitted.force[matID_master] -= gMass[n][node] * (Dv / delT);
                transmitted.all_material_force -= gMass[n][node] * (Dv / delT);
              } // if normalDeltaVel > 0
            }   // if separation
          }     // if mass of both matls>0
        }       // for matls
      } else {  // alpha>=0
        for (int n = 0; n < numMatls; n++) {
          if (n == d_material && gMass[n][node] > 1.e-99) {
            gVelocity_star[n][node] = rigid_vel;
          }
        }
      }
    } // for Node Iterator

    // Put the sumTransmittedForce contribution into the reduction variables
    if (d_mpm_flags->d_reductionVars->sumTransmittedForce) {
      new_dw->put(sumvec_vartype(transmitted.all_material_force),
                  d_mpm_labels->SumTransmittedForceLabel,
                  nullptr,
                  -1);
      new_dw->put_sum_vartype(
        transmitted.force, d_mpm_labels->SumTransmittedForceLabel, matls);
    }

  } // loop over patches

  //__________________________________
  //  reduction Vars
  for (int n = 0; n < numMatls; n++) {
    if (n != d_material) {
      int dwi = d_mat_manager->getMaterial("MPM", n)->getDWIndex();
      reaction.force[matID_master] += reaction.force[dwi];
      reaction.torque[matID_master] += reaction.torque[dwi];
    }
  }

  for (int n = 0; n < numMatls; n++) {
    int dwi = d_mat_manager->getMaterial("MPM", n)->getDWIndex();

    if (numMatls > 1) { // ignore for single matl problems
      new_dw->put(sumvec_vartype(reaction.force[dwi]),
                  d_mpm_labels->RigidReactionForceLabel,
                  nullptr,
                  dwi);
      new_dw->put(sumvec_vartype(reaction.torque[dwi]),
                  d_mpm_labels->RigidReactionTorqueLabel,
                  nullptr,
                  dwi);
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
