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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/FrictionContactLR.h>

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
using std::string;
using std::vector;

FrictionContactLR::FrictionContactLR(const ProcessorGroup* myworld,
                                     const MaterialManagerP& mat_manager,
                                     const MPMLabel* labels,
                                     const MPMFlags* flags,
                                     ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  d_use_logistic_regression = true;

  ps->require("mu", d_mu);
  ps->get("volume_constraint", d_vol_const);
  ps->get("one_or_two_step", d_one_or_two_step);
}

void
FrictionContactLR::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "friction_LR");
  contact_ps->appendElement("mu", d_mu);
  contact_ps->appendElement("volume_constraint", d_vol_const);
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  d_matls.outputProblemSpec(contact_ps);
}

void
FrictionContactLR::addComputesAndRequires(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls,
                                          const VarLabel* label)
{
  Task* t = scinew Task("Friction::exchangeMomentum",
                        this,
                        &FrictionContactLR::exchangeMomentum,
                        label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMatlProminenceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gAlphaMaterialLabel, Ghost::None);
  t->requires(
    Task::NewDW, d_mpm_labels->gNormAlphaToBetaLabel, z_matl, Ghost::None);

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
FrictionContactLR::exchangeMomentum(const ProcessorGroup* pg,
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
}

void
FrictionContactLR::exMomInterpolated(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* matls,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  std::vector<constNCdouble> gMass(numMatls);
  std::vector<constNCdouble> gVolume(numMatls);
  std::vector<constNCdouble> gMatlProminence(numMatls);
  std::vector<NCVector> gVelocity(numMatls);

  Ghost::GhostType gnone = Ghost::None;

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double cell_vol    = dx.x() * dx.y() * dx.z();

    constNCdouble NC_CCweight;
    constNCint gAlphaMaterial;
    constNCVector gNormAlphaToBeta;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->get(
      gAlphaMaterial, d_mpm_labels->gAlphaMaterialLabel, 0, patch, gnone, 0);
    new_dw->get(gNormAlphaToBeta,
                d_mpm_labels->gNormAlphaToBetaLabel,
                0,
                patch,
                gnone,
                0);

    // First, calculate the gradient of the mass everywhere
    // normalize it, and save it in surfNorm
    for (int mat = 0; mat < numMatls; mat++) {
      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();
      new_dw->get(gMass[mat], d_mpm_labels->gMassLabel, matID, patch, gnone, 0);
      new_dw->get(
        gVolume[mat], d_mpm_labels->gVolumeLabel, matID, patch, gnone, 0);
      new_dw->get(gMatlProminence[mat],
                  d_mpm_labels->gMatlProminenceLabel,
                  matID,
                  patch,
                  gnone,
                  0);

      new_dw->getModifiable(
        gVelocity[mat], d_mpm_labels->gVelocityLabel, matID, patch);
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      Vector centerOfMassVelocity(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;

      int alpha = gAlphaMaterial[node];
      // Need to think whether centerOfMass(Stuff) should
      // only include current material and alpha material
      // Why include materials that may be putting mass on the node
      // but aren't near enough to be in proper contact.
      for (int mat = 0; mat < numMatls; mat++) {
        if (!d_matls.requested(mat)) {
          continue;
        }
        centerOfMassVelocity += gVelocity[mat][node] * gMass[mat][node];
        centerOfMassMass += gMass[mat][node];
        totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
      }

      if (alpha >= 0) { // Only work on nodes where alpha!=-99
        centerOfMassVelocity /= centerOfMassMass;

        if (d_mpm_flags->d_axisymmetric) {
          // Nodal volume isn't constant for axisymmetry
          // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
          double r = std::min((patch->getNodePosition(node)).x(), .5 * dx.x());
          cell_vol = r * dx.x() * dx.y();
        }

        // Only apply contact if the node is full relative to a constraint
        if ((totalNodalVol / cell_vol) > d_vol_const) {

          // Loop over materials.  Only proceed if velocity field mass
          // is nonzero (not numerical noise) and the difference from
          // the centerOfMassVelocity is nonzero (More than one velocity
          // field is contributing to grid vertex).
          for (int mat = 0; mat < numMatls; mat++) {
            if (!d_matls.requested(mat) || mat == alpha) {
              continue;
            }

            double mass_mat   = gMass[mat][node];
            double mass_alpha = gMass[alpha][node];
            if (mass_mat > 1.e-16 && mass_alpha > 1.0e-16) {
              // There is mass of material beta at this node
              // Check relative separation of the material prominence
              double separation =
                gMatlProminence[mat][node] - gMatlProminence[alpha][node];
              // If that separation is negative, the matls have overlapped
              if (separation <= 0.01 * dx.x()) {
                Vector deltaVelocity =
                  gVelocity[mat][node] - centerOfMassVelocity;
                Vector normal         = -1.0 * gNormAlphaToBeta[node];
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
                    Min(d_mu, tangentDeltaVelocity / std::abs(normalDeltaVel));

                  // Calculate velocity change needed to enforce contact
                  Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                            std::abs(normalDeltaVel);

                  double ff =
                    std::max(1.0, (.01 * dx.x() - separation) / .01 * dx.x());
                  Dv *= ff;
                  Vector DvAlpha = -Dv * mass_mat / mass_alpha;
                  gVelocity[mat][node] += Dv;
                  gVelocity[alpha][node] += DvAlpha;
                } // if (relative velocity) * normal < 0
              }   // if separation
            }     // if !compare && !compare
          }       // matls
        }         // if (volume constraint)
      }           // if(alpha > 0)
    }             // NodeIterator
  }               // patches
}

void
FrictionContactLR::exMomIntegrated(const ProcessorGroup*,
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
  std::vector<constNCdouble> gMatlProminence(numMatls);
  std::vector<NCVector> gVelocity_star(numMatls);

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double cell_vol    = dx.x() * dx.y() * dx.z();
    constNCdouble NC_CCweight;
    constNCint gAlphaMaterial;
    constNCVector gNormAlphaToBeta;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->get(
      gAlphaMaterial, d_mpm_labels->gAlphaMaterialLabel, 0, patch, gnone, 0);
    new_dw->get(gNormAlphaToBeta,
                d_mpm_labels->gNormAlphaToBetaLabel,
                0,
                patch,
                gnone,
                0);

    for (int mat = 0; mat < matls->size(); mat++) {
      int matID = d_mat_manager->getMaterial("MPM", mat)->getDWIndex();
      new_dw->get(gMass[mat], d_mpm_labels->gMassLabel, matID, patch, gnone, 0);
      new_dw->get(
        gVolume[mat], d_mpm_labels->gVolumeLabel, matID, patch, gnone, 0);
      new_dw->get(gMatlProminence[mat],
                  d_mpm_labels->gMatlProminenceLabel,
                  matID,
                  patch,
                  gnone,
                  0);

      new_dw->getModifiable(
        gVelocity_star[mat], d_mpm_labels->gVelocityStarLabel, matID, patch);
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      Vector centerOfMassVelocity(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      int alpha               = gAlphaMaterial[node];
      for (int mat = 0; mat < numMatls; mat++) {
        if (!d_matls.requested(mat)) {
          continue;
        }
        centerOfMassVelocity += gVelocity_star[mat][node] * gMass[mat][node];
        centerOfMassMass += gMass[mat][node];
        totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
      }

      if (alpha >= 0) { // Only work on nodes where alpha!=-99
        centerOfMassVelocity /= centerOfMassMass;
        if (d_mpm_flags->d_axisymmetric) {
          // Nodal volume isn't constant for axisymmetry
          // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
          double r = std::min((patch->getNodePosition(node)).x(), .5 * dx.x());
          cell_vol = r * dx.x() * dx.y();
        }

        // Only apply contact if the node is full relative to a constraint
        if ((totalNodalVol / cell_vol) > d_vol_const) {

          // Loop over materials.  Only proceed if velocity field mass
          // is nonzero (not numerical noise) and the difference from
          // the centerOfMassVelocity is nonzero (More than one velocity
          // field is contributing to grid vertex).
          for (int mat = 0; mat < numMatls; mat++) {
            if (!d_matls.requested(mat) || mat == alpha) {
              continue;
            }
            double mass_mat   = gMass[mat][node];
            double mass_alpha = gMass[alpha][node];
            if (mass_mat > 1.e-16 && mass_alpha > 1.0e-16) {
              double separation =
                gMatlProminence[mat][node] - gMatlProminence[alpha][node];
              if (separation <= 0.01 * dx.x()) {
                Vector deltaVelocity =
                  gVelocity_star[mat][node] - centerOfMassVelocity;
                Vector normal         = -1.0 * gNormAlphaToBeta[node];
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
                    Min(d_mu, tangentDeltaVelocity / std::abs(normalDeltaVel));

                  // Calculate velocity change needed to enforce contact
                  Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                            std::abs(normalDeltaVel);

                  double ff =
                    std::max(1.0, (.01 * dx.x() - separation) / .01 * dx.x());
                  Dv *= ff;
                  gVelocity_star[mat][node] += Dv;
                  Vector DvAlpha = -Dv * mass_mat / mass_alpha;
                  gVelocity_star[alpha][node] += DvAlpha;
                } // if (relative velocity) * normal < 0
              }   // if separation
            }     // if mass[beta] > 0
          }       // matls
        }         // if (volume constraint)
      }           // if(alpha > 0)
    }             // nodeiterator
  }               // patches
}
