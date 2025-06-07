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

// NodalSVF_Contact.cc
// This is a contact model developed by Peter Mackenzie. Details about the
// derivation can be found in "Modeling Strategies for Multiphase Drag
// Interactions Using the Material Point Method" (Mackenzie, et al; 2011).
// Of the interaction models proposed in their paper, this particular
// contact model for MPM in Uintah can simulate the Nodal Bang Bang method
// OR the Nodal Smoothed Volume (SVF) Fraction Method. The Nodal Bang Bang
// method is less expensive than Nodal SVF, but much less accurate. These
// two methods, which are of the Node-based type, register interaction
// proportional to the cell volume (dx*dy*dz) between two or more phases.
// As a result, over-estimates of the interaction occur in cells where the
// interacting materials do not completely occupy the computational cell.
// Interaction in this model is quantified by an interaction parameter,
// mu (N/m^4), and the velocity difference between the phases.  Other
// simple Coulomb friction or viscous fluid interaction models can be
// substituted.

#include <CCA/Components/MPM/Contact/NodalSVFContact.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
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
#include <fstream>
#include <iostream>
#include <vector>

using namespace Uintah;
using std::vector;

NodalSVFContact::NodalSVFContact(const ProcessorGroup* myworld,
                                 const MaterialManagerP& mat_manager,
                                 const MPMLabel* labels,
                                 const MPMFlags* flags,
                                 ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  ps->require("myu", d_myu);
  ps->require("use_svf", d_svf);

  ps->get("materials", d_materials);

  if (d_materials.size() > 2) {
    throw ProblemSetupException(" You may only specify two materials in the "
                                "input file per contact block for Nodal SVF.",
                                __FILE__,
                                __LINE__);
  }
}

void
NodalSVFContact::setContactMaterialAttributes()
{
}

void
NodalSVFContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "nodal_svf");
  contact_ps->appendElement("myu", d_myu);
  contact_ps->appendElement("use_svf", d_svf);
  d_matls.outputProblemSpec(contact_ps);
}

void
NodalSVFContact::exchangeMomentum(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw,
                                  const VarLabel* gVelocity_label)
{
  int numMatls = matls->size();
  int alpha    = 0;
  int beta     = 0;
  int n        = 0;
  for (int m = 0; m < numMatls; m++) {
    if ((d_matls.requested(m)) && (n == 0)) {
      alpha = matls->get(m);
      n++;
    } else {
      beta = matls->get(m);
    }
  }

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch     = patches->get(p);
    Ghost::GhostType gnone = Ghost::None;
    delt_vartype delT;

    double dx      = patch->dCell().x();
    double dy      = patch->dCell().y();
    double dz      = patch->dCell().z();
    double cellVol = dx * dy * dz;
    double coeff   = cellVol * d_myu;
    double factor;

    constNCVariable<double> NC_CCweight;
    std::vector<constNCVariable<double>> gMass(numMatls);
    std::vector<constNCVariable<double>> gVolume(numMatls);
    std::vector<NCVariable<double>> gSVF(numMatls);
    std::vector<NCVariable<Vector>> gVelocity_star(numMatls);
    std::vector<NCVariable<Vector>> gVelocity_old(numMatls);
    std::vector<NCVariable<Vector>> gForce(numMatls);

    //---------- Retrieve necessary data from DataWarehouse
    //------------------------------------------------
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);

    for (int m = 0; m < numMatls; m++) {
      int dwi = d_mat_manager->getMaterial("MPM", m)->getDWIndex();
      new_dw->get(gMass[dwi], d_mpm_labels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(
        gVolume[dwi], d_mpm_labels->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->getModifiable(gVelocity_star[dwi], gVelocity_label, dwi, patch);
      new_dw->allocateTemporary(gSVF[dwi], patch, gnone, 0);
      new_dw->allocateTemporary(gVelocity_old[dwi], patch, gnone, 0);
      new_dw->allocateTemporary(gForce[dwi], patch, gnone, 0);
    } // for m=0:numMatls

    //----------- Calculate Interaction Force
    //-----------------------------------
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {

      IntVector c             = *iter;
      gVelocity_old[beta][c]  = gVelocity_star[beta][c];
      gVelocity_old[alpha][c] = gVelocity_star[alpha][c];
      gSVF[beta][c]  = 8.0 * NC_CCweight[c] * gVolume[beta][c] / cellVol;
      gSVF[alpha][c] = 8.0 * NC_CCweight[c] * gVolume[alpha][c] / cellVol;

      // Calculate the appropriate value of "factor" based on whether using SVF.
      if (d_svf) {
        factor = coeff * gSVF[beta][c] * gSVF[alpha][c];
      } else {
        factor = coeff;
      }

      // "If using the model with svf calculation," or "if mass is present on
      // both nodes," calculate a non-zero interaction force based on velocity
      // difference and the appropriate value of "factor".
      if ((d_svf) ||
          (gMass[beta][c] > 1.0e-100 && gMass[alpha][c] > 1.0e-100)) {
        gForce[beta][c] =
          factor * (gVelocity_old[alpha][c] - gVelocity_old[beta][c]);
        gForce[alpha][c] =
          factor * (gVelocity_old[beta][c] - gVelocity_old[alpha][c]);

      } else {
        gForce[beta][c]  = Vector(0.0, 0.0, 0.0);
        gForce[alpha][c] = Vector(0.0, 0.0, 0.0);
      }

      //-- Calculate Updated Velocity ------------------------------------
      gVelocity_star[beta][c] +=
        (gForce[beta][c] / (8.0 * NC_CCweight[c] * gMass[beta][c])) * delT;
      gVelocity_star[alpha][c] +=
        (gForce[alpha][c] / (8.0 * NC_CCweight[c] * gMass[alpha][c])) * delT;

    } // for nodes
  }   // for patches
} // end exmomentumIntegrated

void
NodalSVFContact::addComputesAndRequires(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls,
                                        const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("NodalSVFContact::exchangeMomentum",
                        this,
                        &NodalSVFContact::exchangeMomentum,
                        gVelocity_label);

  const MaterialSubset* mss = matls->getUnion();
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);
  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->modifies(gVelocity_label, mss);
  sched->addTask(t, patches, matls);
}
