/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2023 Biswajit Banerjee
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

#include <CCA/Components/MPM/Contact/FrictionContactLRGuilkey.h>

#include <CCA/Components/MPM/Core/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
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

FrictionContactLRGuilkey::FrictionContactLRGuilkey(
  const ProcessorGroup* myworld,
  const MaterialManagerP& mat_manager,
  const MPMLabel* labels,
  const MPMFlags* flags,
  ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  // Constructor
  d_one_or_two_step = 2;
  ps->get("one_or_two_step", d_one_or_two_step);

  ProblemSpecP properties = ps->findBlock("variable_friction");
  if (!properties) {
    throw ProblemSetupException(
      "**ERROR** No variable properties specified.", __FILE__, __LINE__);
  }
  for (ProblemSpecP varProp = properties->findBlock("entry");
       varProp != nullptr;
       varProp = varProp->findNextBlock("entry")) {
    double C = 0.0;
    double M = 0.0;
    varProp->require("color", C);
    varProp->require("mu", M);
    d_color_mu.push_back(std::make_pair(C, M));
  }
  if (d_color_mu.size() < 2) {
    std::cout << "d_color_mu.size() = " << d_color_mu.size() << std::endl;
    throw ProblemSetupException(
      "**ERROR** Need at least two entries in Var model..", __FILE__, __LINE__);
  }

  ps->require("master_material", d_material);
  if (!d_matls.requested(d_material)) {
    throw ProblemSetupException(
      "**ERROR: master_material not one of requested materials",
      __FILE__,
      __LINE__);
  }
}

void
FrictionContactLRGuilkey::setContactMaterialAttributes()
{
}

void
FrictionContactLRGuilkey::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "friction_LRVar");
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);
  contact_ps->appendElement("master_material", d_material);
  d_matls.outputProblemSpec(contact_ps);

  ProblemSpecP lc_ps = contact_ps->appendChild("variable_friction");
  for (auto& color_mu : d_color_mu) {
    ProblemSpecP time_ps = lc_ps->appendChild("entry");
    time_ps->appendElement("color", color_mu.first);
    time_ps->appendElement("mu", color_mu.second);
  }
}

void
FrictionContactLRGuilkey::addComputesAndRequires(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls,
                                                 const VarLabel* label)
{
  Task* t = scinew Task("Friction::exchangeMomentum",
                        this,
                        &FrictionContactLRGuilkey::exchangeMomentum,
                        label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gColorLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gMatlProminenceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->gAlphaMaterialLabel, Ghost::None);
  t->requires(
    Task::NewDW, d_mpm_labels->gNormAlphaToBetaLabel, z_matl, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->modifies(label, mss);

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}

void
FrictionContactLRGuilkey::exchangeMomentum(const ProcessorGroup* pg,
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
FrictionContactLRGuilkey::exMomInterpolated(const ProcessorGroup*,
                                            const PatchSubset* patches,
                                            const MaterialSubset* matls,
                                            DataWarehouse* old_dw,
                                            DataWarehouse* new_dw)
{

  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once
  std::vector<constNCVariable<double>> gMass(numMatls), gColor(numMatls);
  std::vector<constNCVariable<double>> gVolume(numMatls);
  std::vector<constNCVariable<double>> gMatlProminence(numMatls);
  std::vector<NCVariable<Vector>> gVelocity(numMatls);

  Ghost::GhostType gnone = Ghost::None;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    //    Vector dx = patch->dCell();
    constNCVariable<double> NC_CCweight;
    constNCVariable<int> alphaMaterial;
    constNCVariable<Vector> normAlphaToBeta;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->get(
      alphaMaterial, d_mpm_labels->gAlphaMaterialLabel, 0, patch, gnone, 0);
    new_dw->get(
      normAlphaToBeta, d_mpm_labels->gNormAlphaToBetaLabel, 0, patch, gnone, 0);

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    // First, calculate the gradient of the mass everywhere
    // normalize it, and stick it in surfNorm
    for (int m = 0; m < numMatls; m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gColor[m], d_mpm_labels->gColorLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->get(gMatlProminence[m],
                  d_mpm_labels->gMatlProminenceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(
        gVelocity[m], d_mpm_labels->gVelocityLabel, dwi, patch);
    } // loop over matls

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector centerOfMassVelocity(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      int alpha               = alphaMaterial[c];
      // Need to think whether centerOfMass(Stuff) should
      // only include current material and alpha material
      // Why include materials that may be putting mass on the node
      // but aren't near enough to be in proper contact.
      for (int n = 0; n < numMatls; n++) {
        if (!d_matls.requested(n)) {
          continue;
        }
        centerOfMassVelocity += gVelocity[n][c] * gMass[n][c];
        centerOfMassMass += gMass[n][c];
        totalNodalVol += gVolume[n][c] * 8.0 * NC_CCweight[c];
      }

      if (alpha >= 0) { // Only work on nodes where alpha!=-99
        centerOfMassVelocity /= centerOfMassMass;

        double gC = gColor[d_material][c];
        double mu = findMuFromColor(gC);

        // Loop over materials.  Only proceed if velocity field mass
        // is nonzero (not numerical noise) and the difference from
        // the centerOfMassVelocity is nonzero (More than one velocity
        // field is contributing to grid vertex).
        for (int n = 0; n < numMatls; n++) {
          if (!d_matls.requested(n)) {
            continue;
          }
          if (n == alpha) {
            continue;
          }
          double mass = gMass[n][c];
          if (mass > 1.e-16) { // There is mass of material beta at this node
            // Check relative separation of the material prominence
            double separation =
              gMatlProminence[n][c] - gMatlProminence[alpha][c];
            // If that separation is negative, the matls have overlapped
            if (separation <= 0.0) {
              Vector deltaVelocity  = gVelocity[n][c] - centerOfMassVelocity;
              Vector normal         = -1.0 * normAlphaToBeta[c];
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
                  Min(mu, tangentDeltaVelocity / fabs(normalDeltaVel));

                // Calculate velocity change needed to enforce contact
                Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                          fabs(normalDeltaVel);

#if 0
              // Define contact algorithm imposed strain, find maximum
              Vector epsilon=(Dv/dx)*delT;
              double epsilon_max=
                Max(fabs(epsilon.x()),fabs(epsilon.y()),fabs(epsilon.z()));
              if(!compare(epsilon_max,0.0)){
                 epsilon_max *= Max(1.0, mass/(centerOfMassMass-mass));

                 // Scale velocity change if contact algorithm
                 // imposed strain is too large.
                 double ff=Min(epsilon_max,.5)/epsilon_max;
                 Dv=Dv*ff;
              }
#endif
                Vector DvAlpha = -Dv * gMass[n][c] / gMass[alpha][c];
                gVelocity[n][c] += Dv;
                gVelocity[alpha][c] += DvAlpha;
              } // if (relative velocity) * normal < 0
            }   // if separation
          }     // if !compare && !compare
        }       // matls
      }         // if(alpha > 0)
    }           // NodeIterator
  }             // patches
}

void
FrictionContactLRGuilkey::exMomIntegrated(const ProcessorGroup*,
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
  std::vector<constNCVariable<double>> gMass(numMatls), gColor(numMatls);
  std::vector<constNCVariable<double>> gVolume(numMatls);
  std::vector<constNCVariable<double>> gMatlProminence(numMatls);
  std::vector<NCVariable<Vector>> gVelocity_star(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    //    Vector dx = patch->dCell();
    constNCVariable<double> NC_CCweight;
    constNCVariable<int> alphaMaterial;
    constNCVariable<Vector> normAlphaToBeta;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->get(
      alphaMaterial, d_mpm_labels->gAlphaMaterialLabel, 0, patch, gnone, 0);
    new_dw->get(
      normAlphaToBeta, d_mpm_labels->gNormAlphaToBetaLabel, 0, patch, gnone, 0);

    // Retrieve necessary data from DataWarehouse
    for (int m = 0; m < matls->size(); m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gColor[m], d_mpm_labels->gColorLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->get(gMatlProminence[m],
                  d_mpm_labels->gMatlProminenceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(
        gVelocity_star[m], d_mpm_labels->gVelocityStarLabel, dwi, patch);
    }

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector centerOfMassVelocity(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      int alpha               = alphaMaterial[c];
      for (int n = 0; n < numMatls; n++) {
        if (!d_matls.requested(n)) {
          continue;
        }
        centerOfMassVelocity += gVelocity_star[n][c] * gMass[n][c];
        centerOfMassMass += gMass[n][c];
        totalNodalVol += gVolume[n][c] * 8.0 * NC_CCweight[c];
      }

      if (alpha >= 0) { // Only work on nodes where alpha!=-99
        centerOfMassVelocity /= centerOfMassMass;

        double gC = gColor[d_material][c];
        double mu = findMuFromColor(gC);

        // Only apply contact if the node is full relative to a constraint

        // Loop over materials.  Only proceed if velocity field mass
        // is nonzero (not numerical noise) and the difference from
        // the centerOfMassVelocity is nonzero (More than one velocity
        // field is contributing to grid vertex).
        for (int n = 0; n < numMatls; n++) {
          if (!d_matls.requested(n)) {
            continue;
          }
          if (n == alpha) {
            continue;
          }
          double mass = gMass[n][c];
          if (mass > 1.e-16) {
            double separation =
              gMatlProminence[n][c] - gMatlProminence[alpha][c];
            if (separation <= 0.0) {
              Vector deltaVelocity =
                gVelocity_star[n][c] - centerOfMassVelocity;
              Vector normal         = -1.0 * normAlphaToBeta[c];
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
                  Min(mu, tangentDeltaVelocity / fabs(normalDeltaVel));
                // Calculate velocity change needed to enforce contact
                Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                          fabs(normalDeltaVel);

#if 0
              // Define contact algorithm imposed strain, find maximum
              Vector epsilon=(Dv/dx)*delT;
              double epsilon_max=
                Max(fabs(epsilon.x()),fabs(epsilon.y()),fabs(epsilon.z()));
              if(!compare(epsilon_max,0.0)){
                epsilon_max *= Max(1.0, mass/(centerOfMassMass-mass));

                // Scale velocity change if contact algorithm
                // imposed strain is too large.
                double ff=Min(epsilon_max,.5)/epsilon_max;
                Dv=Dv*ff;
              }
#endif
                gVelocity_star[n][c] += Dv;
                Vector DvAlpha = -Dv * gMass[n][c] / gMass[alpha][c];
                gVelocity_star[alpha][c] += DvAlpha;
              } // if (relative velocity) * normal < 0
            }   // if separation
          }     // if mass[beta] > 0
        }       // matls
      }         // if(alpha > 0)
    }           // nodeiterator
  }             // patches
}
