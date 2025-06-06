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

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/FrictionContact.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
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
#include <Core/Math/Matrix3.h>
#include <Core/Math/MiscMath.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace Uintah;
using std::string;
using std::vector;

FrictionContact::FrictionContact(const ProcessorGroup* myworld,
                                 const MaterialManagerP& mat_manager,
                                 const MPMLabel* labels,
                                 const MPMFlags* flags,
                                 ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  ps->require("mu", d_mu);
  ps->get("volume_constraint", d_vol_const);
  ps->get("one_or_two_step", d_one_or_two_step);

  // Hardcoded normal for objects that can be represented in
  // special coordinate systems (cylindrical/spherical)
  d_hardcodedNormals = false;
  ps->get("use_hardcoded_normals", d_hardcodedNormals);
  if (d_hardcodedNormals) {

    // Look for hardcoded_normal block (one set of normals per material)
    for (ProblemSpecP material_ps = ps->findBlock("hardcoded_normal");
         material_ps != 0;
         material_ps = material_ps->findNextBlock("hardcoded_normal")) {

      // Get the material index (**TODO** Check upper limit too)
      int matIndex = -1;
      material_ps->get("material_index", matIndex);
      if (matIndex < 0) {
        std::ostringstream out;
        out << "*ERROR** Invalid material index " << matIndex
            << " in hardcoded normals.";
        out << " Choose values between 0 and num_matls-1";
        throw ProblemSetupException(out.str(), __FILE__, __LINE__);
      }
      d_matIndex.push_back(matIndex);

      // Get the coordinate system
      ProblemSpecP normal_ps = material_ps->findBlock("coordinate_system");
      if (normal_ps) {

        std::string type("cartesian");
        Vector axisDir(1.0, 0.0, 0.0);
        Point center(0.0, 0.0, 0.0);

        normal_ps->getAttribute("type", type);
        if (type == "cylindrical") {
          d_coordType.push_back(NormalCoordSystem::CYLINDRICAL);
          normal_ps->require("axis", axisDir);  // axis direction
          normal_ps->require("center", center); // center of axis
        } else if (type == "spherical") {
          d_coordType.push_back(NormalCoordSystem::SPHERICAL);
          normal_ps->require("center", center); // center of sphere
        } else {
          d_coordType.push_back(NormalCoordSystem::CARTESIAN);
          normal_ps->require("axis", axisDir); // axis direction
        }

        d_type.push_back(type);
        d_center.push_back(center);

        if (!(axisDir.length() > 0.0)) {
          std::ostringstream out;
          out << "**ERROR** Invalid axis direction " << axisDir
              << " in hardcoded normals.";
          throw ProblemSetupException(out.str(), __FILE__, __LINE__);
        }
        axisDir.normalize();
        d_axisDir.push_back(axisDir);
      }
    }
  } // End if hardcoded normals
}

void
FrictionContact::setContactMaterialAttributes()
{
}

void
FrictionContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "friction");
  contact_ps->appendElement("mu", d_mu);
  contact_ps->appendElement("volume_constraint", d_vol_const);
  contact_ps->appendElement("one_or_two_step", d_one_or_two_step);

  contact_ps->appendElement("use_hardcoded_normals", d_hardcodedNormals);
  for (unsigned int ii = 0; ii < d_matIndex.size(); ii++) {
    ProblemSpecP hardcoded = contact_ps->appendChild("hardcoded_normal");
    hardcoded->appendElement("material_index", d_matIndex[ii]);
    ProblemSpecP normal_ps = hardcoded->appendChild("coordinate_system");
    normal_ps->setAttribute("type", d_type[ii]);
    normal_ps->appendElement("axis", d_axisDir[ii]);
    normal_ps->appendElement("center", d_center[ii]);
  }

  d_matls.outputProblemSpec(contact_ps);
}

void
FrictionContact::exchangeMomentum(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw,
                                  const VarLabel* gVelocity_label)
{
  // If one_step, only do exchange for gVelocity_star
  if (d_one_or_two_step == 1 &&
      gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  Ghost::GhostType gnone = Ghost::None;

  int numMatls = d_mat_manager->getNumMaterials("MPM");
  ASSERTEQ(numMatls, matls->size());

  // Need access to all velocity fields at once, so store in
  // vectors of NCVariables
  std::vector<constNCVariable<double>> gMass(numMatls);
  std::vector<constNCVariable<double>> gVolume(numMatls);
  std::vector<NCVariable<Vector>> gVelocity_star(numMatls);
  std::vector<constNCVariable<double>> normtraction(numMatls);
  std::vector<NCVariable<double>> frictionWork(numMatls);
  std::vector<constNCVariable<Vector>> gSurfNormal(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();
    double cell_vol    = dx.x() * dx.y() * dx.z();
    constNCVariable<double> NC_CCweight;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

    // Retrieve necessary data from DataWarehouse
    for (int m = 0; m < matls->size(); m++) {
      int dwi = matls->get(m);
      new_dw->get(gMass[m], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(normtraction[m],
                  d_mpm_labels->gNormTractionLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(
        gSurfNormal[m], d_mpm_labels->gSurfNormLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, dwi, patch);
      new_dw->getModifiable(
        frictionWork[m], d_mpm_labels->frictionalWorkLabel, dwi, patch);
    }

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));
    double epsilon_max_max = 0.0;

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector centerOfMassMom(0., 0., 0.);
      double centerOfMassMass = 0.0;
      double totalNodalVol    = 0.0;
      for (int n = 0; n < numMatls; n++) {
        if (!d_matls.requested(n)) {
          continue;
        }
        double mass = gMass[n][c];
        centerOfMassMom += gVelocity_star[n][c] * mass;
        centerOfMassMass += mass;
        totalNodalVol += gVolume[n][c] * 8.0 * NC_CCweight[c];
      }

      // Apply Coulomb friction contact
      // For grid points with mass calculate velocity
      if (!compare(centerOfMassMass, 0.0)) {
        Vector centerOfMassVelocity = centerOfMassMom / centerOfMassMass;

        if (d_mpm_flags->d_axisymmetric) {
          // Nodal volume isn't constant for axisymmetry
          // volume = r*dr*dtheta*dy  (dtheta = 1 radian)
          double r = std::min((patch->getNodePosition(c)).x(), .5 * dx.x());
          cell_vol = r * dx.x() * dx.y();
        }

        // Only apply contact if the node is nearly "full".  There are
        // two options:

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
          for (int n = 0; n < numMatls; n++) {
            if (!d_matls.requested(n)) {
              continue;
            }
            Vector deltaVelocity = gVelocity_star[n][c] - centerOfMassVelocity;
            double mass          = gMass[n][c];
            if (!compare(mass / centerOfMassMass, 0.0) &&
                !compare(mass - centerOfMassMass, 0.0)) {

              // Apply frictional contact IF the surface is in compression
              // OR the surface is stress free and approaching.
              // Otherwise apply free surface conditions (do nothing).
              Vector normal         = gSurfNormal[n][c];
              double normalDeltaVel = Dot(deltaVelocity, normal);

              Vector Dv(0., 0., 0.);
              double Tn = normtraction[n][c];
              if ((Tn < 0.0) ||
                  (compare(fabs(Tn), 0.0) && normalDeltaVel > 0.0)) {

                // Simplify algorithm in case where approach velocity
                // is in direction of surface normal (no slip).
                Vector normal_normaldV = normal * normalDeltaVel;
                Vector dV_normaldV     = deltaVelocity - normal_normaldV;
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
                  Dv = -normal_normaldV - surfaceTangent * frictionCoefficient *
                                            fabs(normalDeltaVel);

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
                  Dv        = Dv * ff;
                }
                Dv = scale_factor * Dv;
                gVelocity_star[n][c] += Dv;
              } // traction
            }   // if !compare && !compare
          }     // for numMatls
        }       // volume constraint
      }         // if centerofmass > 0
    }           // nodeiterator

    //  print out epsilon_max_max
    //  static int ts=0;
    //  static ofstream tmpout("max_strain.dat");

    //  tmpout << ts << " " << epsilon_max_max << std::endl;
    //  ts++;

    // This converts frictional work into a temperature rate
    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", m));

      if (!d_matls.requested(m)) {
        for (NodeIterator iter = patch->getNodeIterator(); !iter.done();
             iter++) {
          frictionWork[m][*iter] = 0;
        }
      } else {
        double c_v = mpm_matl->getSpecificHeat();
        for (NodeIterator iter = patch->getNodeIterator(); !iter.done();
             iter++) {
          IntVector c = *iter;
          frictionWork[m][c] /= (c_v * gMass[m][c] * delT);
          if (frictionWork[m][c] < 0.0) {
            std::cout << "dT/dt is negative: " << frictionWork[m][c]
                      << std::endl;
          }
        }
      }
    }
  }
}

void
FrictionContact::addComputesAndRequires(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls,
                                        const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("Friction::exchangeMomentum",
                        this,
                        &FrictionContact::exchangeMomentum,
                        gVelocity_label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);
  t->needs(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gNormTractionLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gSurfNormLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->modifies(gVelocity_label, mss);
  t->modifies(d_mpm_labels->frictionalWorkLabel, mss);

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}
