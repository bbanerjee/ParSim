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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/ApproachContact.h>
#include <CCA/Components/MPM/MPMBoundCond.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <vector>

using namespace Uintah;
using std::vector;
using std::string;

ApproachContact::ApproachContact(const ProcessorGroup* myworld,
                                 ProblemSpecP& ps, MaterialManagerP& d_sS,
                                 MPMLabel* Mlb, MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  d_vol_const = 0.;

  ps->require("mu", d_mu);
  ps->get("volume_constraint", d_vol_const);

  d_mat_manager = d_sS;
  d_needNormals = true;

  if (flag->d_8or27 == 8) {
    NGP = 1;
    NGN = 1;
  } else {
    NGP = 2;
    NGN = 2;
  }
}

ApproachContact::~ApproachContact()
{
}

void
ApproachContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "approach");
  contact_ps->appendElement("mu", d_mu);
  contact_ps->appendElement("volume_constraint", d_vol_const);

  d_matls.outputProblemSpec(contact_ps);
}

void
ApproachContact::exchangeMomentum(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw, DataWarehouse* new_dw,
                                  const VarLabel* gVelocity_label)
{
  Ghost::GhostType gnone = Ghost::None;

  int numMatls = d_mat_manager->getNumMPMMatls();
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once, so store in
  // vectors of NCVariables
  std::vector<constNCVariable<double>> gmass(numMatls);
  std::vector<constNCVariable<double>> gvolume(numMatls);
  std::vector<NCVariable<Vector>> gvelocity_star(numMatls);
  std::vector<NCVariable<double>> frictionWork(numMatls);
  std::vector<constNCVariable<Vector>> gsurfnorm(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();
    double cell_vol = dx.x() * dx.y() * dx.z();
    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gnone, 0);

    // Retrieve necessary data from DataWarehouse
    for (int m = 0; m < matls->size(); m++) {
      int dwi = matls->get(m);
      new_dw->get(gmass[m], lb->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gsurfnorm[m], lb->gSurfNormLabel, dwi, patch, gnone, 0);
      new_dw->get(gvolume[m], lb->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->getModifiable(gvelocity_star[m], gVelocity_label, dwi, patch);
      new_dw->getModifiable(frictionWork[m], lb->frictionalWorkLabel, dwi,
                            patch);
    }

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    double epsilon_max_max = 0.0;

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector centerOfMassMom(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol = 0.0;
      for (int n = 0; n < numMatls; n++) {
        if (!d_matls.requested(n))
          continue;
        double mass = gmass[n][c];
        centerOfMassMom += gvelocity_star[n][c] * mass;
        centerOfMassMass += mass;
        totalNodalVol += gvolume[n][c] * 8.0 * NC_CCweight[c];
      }

      // Apply Coulomb friction contact
      // For grid points with mass calculate velocity
      if (!compare(centerOfMassMass, 0.0)) {
        Vector centerOfMassVelocity = centerOfMassMom / centerOfMassMass;

        // Only apply contact if the node is nearly "full".  There are
        // two options:

        if ((totalNodalVol / cell_vol) > d_vol_const) {
          double scale_factor = 1.0;

          // 2. This option uses only cell volumes. The idea is that a cell is
          // full
          // if totalNodalVol/cell_vol .ge. 1.0, and the contraint should only
          // be
          // applied when cells are full.  This logic is used when
          // d_vol_const=0.
          // For d_vol_const > 0 the contact forces are ramped up linearly from
          // 0
          // for totalNodalVol/cell_vol .le. 1.0-d_vol_const to 1.0 for
          // totalNodalVol/cell_vol = 1.  Ramping the contact influence seems to
          // help remove a "switching" instability.  A good value seems to be
          // d_vol_const=.05

          //      double scale_factor=0.0;
          //      if(d_vol_const > 0.0){
          //        scale_factor=
          //                 (totalNodalVol/cell_vol-1.+d_vol_const)/d_vol_const;
          //        scale_factor=Max(0.0,scale_factor);
          //      }
          //      else if(totalNodalVol/cell_vol > 1.0){
          //        scale_factor=1.0;
          //      }

          //      if(scale_factor > 0.0){
          //        scale_factor=Min(1.0,scale_factor);
          //       }

          // Loop over velocity fields.  Only proceed if velocity field mass
          // is nonzero (not numerical noise) and the difference from
          // the centerOfMassVelocity is nonzero (More than one velocity
          // field is contributing to grid vertex).
          for (int n = 0; n < numMatls; n++) {
            if (!d_matls.requested(n))
              continue;
            Vector deltaVelocity = gvelocity_star[n][c] - centerOfMassVelocity;
            double mass = gmass[n][c];
            if (!compare(mass / centerOfMassMass, 0.0) &&
                !compare(mass - centerOfMassMass, 0.0)) {

              // Apply frictional contact IF the surface is in compression
              // OR the surface is stress free and approaching.
              // Otherwise apply free surface conditions (do nothing).
              Vector normal = gsurfnorm[n][c];
              double normalDeltaVel = Dot(deltaVelocity, normal);

              Vector Dv(0., 0., 0.);
              if (normalDeltaVel > 0.0) {

                // Simplify algorithm in case where approach velocity
                // is in direction of surface normal (no slip).
                Vector normal_normaldV = normal * normalDeltaVel;
                Vector dV_normaldV = deltaVelocity - normal_normaldV;
                if (compare(dV_normaldV.length2(), 0.0)) {

                  // Calculate velocity change needed to enforce contact
                  Dv = -normal_normaldV;
                }

                // General algorithm, including frictional slip.  The
                // contact velocity change and frictional work are both
                // zero if normalDeltaVel is zero.
                else if (!compare(fabs(normalDeltaVel), 0.0)) {
                  Vector surfaceTangent = dV_normaldV / dV_normaldV.length();
                  double tangentDeltaVelocity =
                    Dot(deltaVelocity, surfaceTangent);
                  double frictionCoefficient =
                    Min(d_mu, tangentDeltaVelocity / fabs(normalDeltaVel));

                  // Calculate velocity change needed to enforce contact
                  Dv =
                    -normal_normaldV -
                    surfaceTangent * frictionCoefficient * fabs(normalDeltaVel);

                  // Calculate work done by the frictional force (only) if
                  // contact slips.  Because the frictional force opposes motion
                  // it is dissipative and should always be negative per the
                  // conventional definition.  However, here it is calculated
                  // as positive (Work=-force*distance).
                  if (compare(frictionCoefficient, d_mu)) {
                    frictionWork[n][c] +=
                      mass * frictionCoefficient *
                      (normalDeltaVel * normalDeltaVel) *
                      (tangentDeltaVelocity / fabs(normalDeltaVel) -
                       frictionCoefficient);
                  }
                }

                // Define contact algorithm imposed strain, find maximum
                Vector epsilon = (Dv / dx) * delT;
                double epsilon_max =
                  Max(fabs(epsilon.x()), fabs(epsilon.y()), fabs(epsilon.z()));
                epsilon_max_max = std::max(epsilon_max, epsilon_max_max);
                if (!compare(epsilon_max, 0.0)) {
                  epsilon_max *= Max(1.0, mass / (centerOfMassMass - mass));

                  // Scale velocity change if contact algorithm imposed strain
                  // is too large.
                  double ff = Min(epsilon_max, .5) / epsilon_max;
                  Dv = Dv * ff;
                }
                Dv = scale_factor * Dv;
                gvelocity_star[n][c] += Dv;
              } // if approaching
            }   // if !compare && !compare
          }     // for numMatls
        }       // volume constraint
      }         // if centerofmass > 0
    }           // nodeiterator

    //  print out epsilon_max_max
    //  static int ts=0;
    //  static ofstream tmpout("max_strain.dat");

    //  tmpout << ts << " " << epsilon_max_max << "\n";
    //  ts++;

    // This converts frictional work into a temperature rate
    for (int m = 0; m < matls->size(); m++) {
      if (!d_matls.requested(m))
        continue;
      MPMMaterial* mpm_matl = d_mat_manager->getMPMMaterial(m);
      double c_v = mpm_matl->getSpecificHeat();
      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        frictionWork[m][c] /= (c_v * gmass[m][c] * delT);
        if (frictionWork[m][c] < 0.0) {
          std::cout << "dT/dt is negative: " << frictionWork[m][c] << "\n";
        }
      }
    }
  }
}

void
ApproachContact::addComputesAndRequires(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls,
                                        const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("ApproachContact::exchangeMomentum", this,
                        &ApproachContact::exchangeMomentum,
                        gVelocity_label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);
  t->requires(Task::OldDW, lb->NC_CCweightLabel, z_matl, Ghost::None);
  t->requires(Task::NewDW, lb->gSurfNormLabel, Ghost::None);
  t->requires(Task::NewDW, lb->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, lb->gVolumeLabel, Ghost::None);
  t->modifies(gVelocity_label, mss);
  if (gVelocity_label == lb->gVelocityLabel) {
    t->computes(lb->frictionalWorkLabel, mss);
  } else {
    t->modifies(lb->frictionalWorkLabel, mss);
  }

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
}
