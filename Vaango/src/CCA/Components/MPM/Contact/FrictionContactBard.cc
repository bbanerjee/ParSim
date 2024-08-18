/*
 * The MIT License
 *
 * Copyright (c) 1997-2022 The University of Utah
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

#include <CCA/Components/MPM/Contact/FrictionContactBard.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>
#include <iostream>
#include <vector>

using namespace Uintah;

FrictionContactBard::FrictionContactBard(const ProcessorGroup* myworld,
                                         const MaterialManagerP& mat_manager,
                                         const MPMLabel* labels,
                                         const MPMFlags* flags,
                                         ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  d_need_normals = true;

  ps->require("mu", d_mu);
  ps->get("volume_constraint", d_vol_const);
  ps->get("separation_factor", d_sep_fac);
  ps->get("one_or_two_step", d_one_or_two_step);
}

void
FrictionContactBard::setContactMaterialAttributes()
{
}

void
FrictionContactBard::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "friction_bard");
  contact_ps->appendElement("mu", d_mu);
  contact_ps->appendElement("volume_constraint", d_vol_const);
  contact_ps->appendElement("separation_factor", d_sep_fac);
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  d_matls.outputProblemSpec(contact_ps);
}

void
FrictionContactBard::addComputesAndRequires(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls,
                                            const VarLabel* label)
{
  Task* t = scinew Task("Friction::exchangeMomentum",
                        this,
                        &FrictionContactBard::exchangeMomentum,
                        label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gSurfNormLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gPositionLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gNormTractionLabel, Ghost::None);

  t->modifies(d_mpm_labels->frictionalWorkLabel, mss);
  if (label == d_mpm_labels->gVelocityLabel) {
    t->modifies(d_mpm_labels->gVelocityLabel, mss);
  } else {
    t->modifies(d_mpm_labels->gVelocityStarLabel, mss);
  }

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}

void
FrictionContactBard::exchangeMomentum(const ProcessorGroup* pg,
                                      const PatchSubset* patches,
                                      const MaterialSubset* matls,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw,
                                      const VarLabel* label)
{
  if (label == d_mpm_labels->gVelocityLabel) {
    if (d_one_or_two_step == 2) {
      exMomInterpolated(pg, patches, matls, old_dw, new_dw);
    }
  } else {
    exMomIntegrated(pg, patches, matls, old_dw, new_dw);
  }
  /*
  if (label->getName() == "g.velocity") {
  } else {
  }
  */
}

void
FrictionContactBard::exMomInterpolated(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* matls,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once
  std::vector<constNCdouble> gMass(numMatls);
  std::vector<constNCdouble> gVolume(numMatls);
  std::vector<constNCPoint> gPosition(numMatls);
  std::vector<constNCVector> gSurfNorm(numMatls);
  std::vector<constNCdouble> gNormTraction(numMatls);
  std::vector<NCVector> gVelocity(numMatls);
  std::vector<NCdouble> gFrictionWork(numMatls);

  Ghost::GhostType gnone = Ghost::None;

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double cell_vol    = dx.x() * dx.y() * dx.z();

    constNCdouble NC_CCweight;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

    // First, calculate the gradient of the mass everywhere
    // normalize it, and save it in surfNorm
    for (int mat = 0; mat < numMatls; mat++) {
      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();
      new_dw->get(gMass[mat], d_mpm_labels->gMassLabel, matID, patch, gnone, 0);
      new_dw->get(
        gVolume[mat], d_mpm_labels->gVolumeLabel, matID, patch, gnone, 0);
      new_dw->get(
        gSurfNorm[mat], d_mpm_labels->gSurfNormLabel, matID, patch, gnone, 0);
      new_dw->get(
        gPosition[mat], d_mpm_labels->gPositionLabel, matID, patch, gnone, 0);
      new_dw->get(gNormTraction[mat],
                  d_mpm_labels->gNormTractionLabel,
                  matID,
                  patch,
                  gnone,
                  0);

      new_dw->getModifiable(
        gVelocity[mat], d_mpm_labels->gVelocityLabel, matID, patch);
      new_dw->getModifiable(
        gFrictionWork[mat], d_mpm_labels->frictionalWorkLabel, matID, patch);
    } // loop over matls

    double sepDis = d_sep_fac * std::cbrt(cell_vol);
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      Vector centerOfMassMom(0., 0., 0.);
      Point centerOfMassPos(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      for (int mat = 0; mat < numMatls; mat++) {
        if (!d_matls.requested(mat)) {
          continue;
        }
        centerOfMassMom += gVelocity[mat][node] * gMass[mat][node];
        centerOfMassPos += gPosition[mat][node].asVector() * gMass[mat][node];
        centerOfMassMass += gMass[mat][node];
        totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
      }
      centerOfMassPos /= centerOfMassMass;

      // Apply Coulomb friction contact
      // For grid points with mass calculate velocity
      if (!compare(centerOfMassMass, 0.0)) {
        Vector centerOfMassVelocity = centerOfMassMom / centerOfMassMass;

        if (d_mpm_flags->d_axisymmetric) {
          // Nodal volume isn't constant for axisymmetry
          // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
          double r = std::min((patch->getNodePosition(node)).x(), .5 * dx.x());
          cell_vol = r * dx.x() * dx.y();
        }

        // Only apply contact if the node is full relative to a constraint
        if ((totalNodalVol / cell_vol) > d_vol_const) {
          double scale_factor = 1.0; // Currently not used, should test again.

          // 2. This option uses only cell volumes.  The idea is that a cell
          //    is full if (totalNodalVol/cell_vol >= 1.0), and the contraint
          //    should only be applied when cells are full.  This logic is used
          //    when d_vol_const=0.
          //    For d_vol_const > 0 the contact forces are ramped up linearly
          //    from 0 for (totalNodalVol/cell_vol <= 1.0-d_vol_const)
          //    to 1.0 for (totalNodalVol/cell_vol = 1).
          //    Ramping the contact influence seems to help remove a "switching"
          //    instability.  A good value seems to be d_vol_const=.05

          //      double scale_factor=0.0;
          //      if(d_vol_const > 0.0){
          //        scale_factor=
          //          (totalNodalVol/cell_vol-1.+d_vol_const)/d_vol_const;
          //        scale_factor=Max(0.0,scale_factor);
          //      }
          //      else if(totalNodalVol/cell_vol > 1.0){
          //        scale_factor=1.0;
          //      }

          //      if(scale_factor > 0.0){
          //        scale_factor=Min(1.0,scale_factor);
          //      }

          // Loop over velocity fields.  Only proceed if velocity field mass
          // is nonzero (not numerical noise) and the difference from
          // the centerOfMassVelocity is nonzero (More than one velocity
          // field is contributing to grid vertex).
          for (int mat = 0; mat < numMatls; mat++) {
            if (!d_matls.requested(mat)) {
              continue;
            }
            double mass          = gMass[mat][node];
            Vector deltaVelocity = gVelocity[mat][node] - centerOfMassVelocity;
            if (!compare(mass / centerOfMassMass, 0.0) &&
                !compare(mass - centerOfMassMass, 0.0)) {

              // Apply frictional contact IF the surface is in compression
              // OR the surface is stress free and approaching.
              // Otherwise apply free surface conditions (do nothing).
              Vector normal = gSurfNorm[mat][node];
              Vector sepVec = (centerOfMassMass / (centerOfMassMass - mass)) *
                              (centerOfMassPos - gPosition[mat][node]);
              double sepScalar = sepVec.length();
              if (sepScalar < sepDis) {
                double normalDeltaVel = Dot(deltaVelocity, normal);
                Vector Dv(0., 0., 0.);
                double Tn = gNormTraction[mat][node];
                if ((Tn < -1.e-12) || (normalDeltaVel > 0.0)) {

                  Vector normal_normaldV = normal * normalDeltaVel;
                  Vector dV_normalDV     = deltaVelocity - normal_normaldV;

                  if (compare(dV_normalDV.length2(), 0.0)) {

                    // Simplify algorithm in case where approach velocity
                    // is in direction of surface normal (no slip).
                    // Calculate velocity change needed to enforce contact
                    Dv = -normal_normaldV;

                  } else if (!compare(std::abs(normalDeltaVel), 0.0)) {

                    // General algorithm, including frictional slip.  The
                    // contact velocity change and frictional work are both
                    // zero if normalDeltaVel is zero.
                    Vector surfaceTangent = dV_normalDV / dV_normalDV.length();
                    double tangentDeltaVelocity =
                      Dot(deltaVelocity, surfaceTangent);
                    double frictionCoefficient = Min(
                      d_mu, tangentDeltaVelocity / std::abs(normalDeltaVel));

                    // Calculate velocity change needed to enforce contact
                    Dv = -normal_normaldV - surfaceTangent *
                                              frictionCoefficient *
                                              std::abs(normalDeltaVel);

                    // Calculate work done by the frictional force (only) if
                    // contact slips.  Because the frictional force opposes
                    // motion it is dissipative and should always be negative
                    // per the conventional definition.  However, here it is
                    // calculated as positive (Work=-force*distance).
                    if (compare(frictionCoefficient, d_mu)) {
                      gFrictionWork[mat][node] =
                        mass * frictionCoefficient *
                        (normalDeltaVel * normalDeltaVel) *
                        (tangentDeltaVelocity / fabs(normalDeltaVel) -
                         frictionCoefficient);
                    }
                  }

                  // Define contact algorithm imposed strain, find maximum
                  Vector epsilon     = (Dv / dx) * delT;
                  double epsilon_max = Max(std::abs(epsilon.x()),
                                           std::abs(epsilon.y()),
                                           std::abs(epsilon.z()));
                  if (!compare(epsilon_max, 0.0)) {
                    epsilon_max *= Max(1.0, mass / (centerOfMassMass - mass));

                    // Scale velocity change if contact algorithm
                    // imposed strain is too large.
                    double ff = Min(epsilon_max, .5) / epsilon_max;
                    Dv        = Dv * ff;
                  }
                  Dv *= scale_factor;
                  gVelocity[mat][node] += Dv;
                } // if traction
              }   // if sepScalar
            }     // if !compare && !compare
          }       // matls
        }         // if (volume constraint)
      }           // if(!compare(centerOfMassMass,0.0))
    }             // NodeIterator
  }               // patches
}

void
FrictionContactBard::exMomIntegrated(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* matls,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;

  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once, so store in
  // vectors of NCVariables
  std::vector<constNCdouble> gMass(numMatls);
  std::vector<constNCdouble> gVolume(numMatls);
  std::vector<constNCPoint> gPosition(numMatls);
  std::vector<NCVector> gVelocity_star(numMatls);
  std::vector<constNCdouble> gNormTraction(numMatls);
  std::vector<NCdouble> gFrictionWork(numMatls);
  std::vector<constNCVector> gSurfNorm(numMatls);

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double cell_vol    = dx.x() * dx.y() * dx.z();

    constNCdouble NC_CCweight;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

    for (int mat = 0; mat < matls->size(); mat++) {
      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();
      new_dw->get(gMass[mat], d_mpm_labels->gMassLabel, matID, patch, gnone, 0);
      new_dw->get(gNormTraction[mat],
                  d_mpm_labels->gNormTractionLabel,
                  matID,
                  patch,
                  gnone,
                  0);
      new_dw->get(
        gSurfNorm[mat], d_mpm_labels->gSurfNormLabel, matID, patch, gnone, 0);
      new_dw->get(
        gPosition[mat], d_mpm_labels->gPositionLabel, matID, patch, gnone, 0);
      new_dw->get(
        gVolume[mat], d_mpm_labels->gVolumeLabel, matID, patch, gnone, 0);

      new_dw->getModifiable(
        gVelocity_star[mat], d_mpm_labels->gVelocityStarLabel, matID, patch);
      new_dw->getModifiable(
        gFrictionWork[mat], d_mpm_labels->frictionalWorkLabel, matID, patch);
    }

    double epsilon_max_max = 0.0;
    double sepDis          = d_sep_fac * cbrt(cell_vol);

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      Vector centerOfMassMom(0., 0., 0.);
      Point centerOfMassPos(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;

      for (int mat = 0; mat < numMatls; mat++) {
        if (!d_matls.requested(mat)) {
          continue;
        }
        double mass = gMass[mat][node];
        centerOfMassMom += gVelocity_star[mat][node] * mass;
        centerOfMassPos += gPosition[mat][node].asVector() * gMass[mat][node];
        centerOfMassMass += mass;
        totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
      }
      centerOfMassPos /= centerOfMassMass;

      // Apply Coulomb friction contact
      // For grid points with mass calculate velocity
      if (!compare(centerOfMassMass, 0.0)) {
        Vector centerOfMassVelocity = centerOfMassMom / centerOfMassMass;

        if (d_mpm_flags->d_axisymmetric) {
          // Nodal volume isn't constant for axisymmetry
          // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
          double r = std::min((patch->getNodePosition(node)).x(), .5 * dx.x());
          cell_vol = r * dx.x() * dx.y();
        }

        // Only apply contact if the node is full relative to a constraint
        if ((totalNodalVol / cell_vol) > d_vol_const) {
          double scale_factor = 1.0;

          // 2. This option uses only cell volumes.  The idea is that a cell
          //    is full if (totalNodalVol/cell_vol >= 1.0), and the contraint
          //    should only be applied when cells are full.  This logic is used
          //    when d_vol_const=0.
          //    For d_vol_const > 0 the contact forces are ramped up linearly
          //    from 0 for (totalNodalVol/cell_vol <= 1.0-d_vol_const)
          //    to 1.0 for (totalNodalVol/cell_vol = 1).
          //    Ramping the contact influence seems to help remove a "switching"
          //    instability.  A good value seems to be d_vol_const=.05

          //      double scale_factor=0.0;
          //      if(d_vol_const > 0.0){
          //        scale_factor=
          //          (totalNodalVol/cell_vol-1.+d_vol_const)/d_vol_const;
          //        scale_factor=Max(0.0,scale_factor);
          //      }
          //      else if(totalNodalVol/cell_vol > 1.0){
          //        scale_factor=1.0;
          //      }

          //      if(scale_factor > 0.0){
          //        scale_factor=Min(1.0,scale_factor);
          //      }

          // Loop over velocity fields.  Only proceed if velocity field mass
          // is nonzero (not numerical noise) and the difference from
          // the centerOfMassVelocity is nonzero (More than one velocity
          // field is contributing to grid vertex).
          for (int mat = 0; mat < numMatls; mat++) {
            if (!d_matls.requested(mat)) {
              continue;
            }
            Vector deltaVelocity =
              gVelocity_star[mat][node] - centerOfMassVelocity;
            double mass = gMass[mat][node];
            if (!compare(mass / centerOfMassMass, 0.0) &&
                !compare(mass - centerOfMassMass, 0.0)) {

              // Apply frictional contact IF the surface is in compression
              // OR the surface is stress free and approaching.
              // Otherwise apply free surface conditions (do nothing).
              Vector normal = gSurfNorm[mat][node];
              Vector sepVec = (centerOfMassMass / (centerOfMassMass - mass)) *
                              (centerOfMassPos - gPosition[mat][node]);
              double sepScalar = sepVec.length();

              if (sepScalar < sepDis) {
                double normalDeltaVel = Dot(deltaVelocity, normal);

                Vector Dv(0., 0., 0.);
                double Tn = gNormTraction[mat][node];
                if ((Tn < -1.e-12) || (normalDeltaVel > 0.0)) {

                  Vector normal_normaldV = normal * normalDeltaVel;
                  Vector dV_normaldV     = deltaVelocity - normal_normaldV;
                  if (compare(dV_normaldV.length2(), 0.0)) {

                    // Simplify algorithm in case where approach velocity
                    // is in direction of surface normal (no slip).
                    // Calculate velocity change needed to enforce contact
                    Dv = -normal_normaldV;

                  } else if (!compare(fabs(normalDeltaVel), 0.0)) {

                    // General algorithm, including frictional slip.  The
                    // contact velocity change and frictional work are both
                    // zero if normalDeltaVel is zero.
                    Vector surfaceTangent = dV_normaldV / dV_normaldV.length();
                    double tangentDeltaVelocity =
                      Dot(deltaVelocity, surfaceTangent);
                    double frictionCoefficient = Min(
                      d_mu, tangentDeltaVelocity / std::abs(normalDeltaVel));

                    // Calculate velocity change needed to enforce contact
                    Dv = -normal_normaldV - surfaceTangent *
                                              frictionCoefficient *
                                              std::abs(normalDeltaVel);

                    // Calculate work done by the frictional force (only) if
                    // contact slips.  Because the frictional force opposes
                    // motion it is dissipative and should always be negative
                    // per the conventional definition.  However, here it is
                    // calculated as positive (Work=-force*distance).
                    if (compare(frictionCoefficient, d_mu)) {
                      gFrictionWork[mat][node] +=
                        mass * frictionCoefficient *
                        (normalDeltaVel * normalDeltaVel) *
                        (tangentDeltaVelocity / std::abs(normalDeltaVel) -
                         frictionCoefficient);
                    }
                  }

                  // Define contact algorithm imposed strain, find maximum
                  Vector epsilon     = (Dv / dx) * delT;
                  double epsilon_max = Max(std::abs(epsilon.x()),
                                           std::abs(epsilon.y()),
                                           std::abs(epsilon.z()));
                  epsilon_max_max    = std::max(epsilon_max, epsilon_max_max);
                  if (!compare(epsilon_max, 0.0)) {
                    epsilon_max *= Max(1.0, mass / (centerOfMassMass - mass));

                    // Scale velocity change if contact algorithm imposed strain
                    // is too large.
                    double ff = Min(epsilon_max, .5) / epsilon_max;
                    Dv *= ff;
                  }
                  Dv *= scale_factor;
                  gVelocity_star[mat][node] += Dv;
                } // traction
              }   // if sepScalar
            }     // if !compare && !compare
          }       // for numMatls
        }         // volume constraint
      }           // if centerofmass > 0
    }             // nodeiterator

    //  print out epsilon_max_max
    //  static int ts=0;
    //  static ofstream tmpout("max_strain.dat");
    //  tmpout << ts << " " << epsilon_max_max << std::endl;
    //  ts++;

    // This converts frictional work into a temperature rate
    for (int mat = 0; mat < matls->size(); mat++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", mat));
      if (!d_matls.requested(mat)) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          gFrictionWork[mat][*iter] = 0;
        }
      } else {
        double c_v = mpm_matl->getSpecificHeat();
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gFrictionWork[mat][node] /= (c_v * gMass[mat][node] * delT);
          if (gFrictionWork[mat][node] < 0.0) {
            std::cout << "dT/dt is negative: " << gFrictionWork[mat][node]
                      << std::endl;
          }
        }
      }
    }
  } // end patch loop
}
