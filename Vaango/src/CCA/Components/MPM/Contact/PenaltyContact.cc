/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
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

#include <CCA/Components/MPM/Contact/PenaltyContact.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <Core/Grid/Labels/MPMLabel.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <CCA/Ports/DataWarehouse.h>
#include <vector>
#include <iostream>

using namespace Uintah;
using std::vector;
using std::string;

PenaltyContact::PenaltyContact(const ProcessorGroup* myworld,
                               ProblemSpecP& ps,
                               SimulationStateP& d_sS,
                               MPMLabel* Mlb,
                               MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  d_vol_const=0.;
  d_oneOrTwoStep = 1;

  ps->require("mu",d_mu);

  d_sharedState = d_sS;

  if (flag->d_8or27 == 8) {
    NGP = 1; NGN = 1;
  } else {
    NGP = 2; NGN = 2;
  }
}

void 
PenaltyContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "penalty");
  contact_ps->appendElement("mu", d_mu);
  d_matls.outputProblemSpec(contact_ps);
}

void 
PenaltyContact::addComputesAndRequires(SchedulerP& sched, 
                                       const PatchSet* patches,
                                       const MaterialSet* matls,
                                       const VarLabel* label)
{
  if (label == lb->gVelocityStarLabel) {
    Task * t = scinew Task("PenaltyContact::exchangeMomentum", 
                           this, 
                           &PenaltyContact::exchangeMomentum,
                           label);

    const MaterialSubset* mss = matls->getUnion();
    t->requires(Task::OldDW, lb->delTLabel);
    t->requires(Task::NewDW, lb->gMassLabel,                  Ghost::None);
    t->requires(Task::NewDW, lb->gVolumeLabel,                Ghost::None);
    t->requires(Task::NewDW, lb->gLSContactForceLabel,        Ghost::None);
    t->modifies(lb->gVelocityStarLabel,  mss);

    sched->addTask(t, patches, matls);
  }
}

void 
PenaltyContact::exchangeMomentum(const ProcessorGroup* pg, 
                                 const PatchSubset* patches,
                                 const MaterialSubset* matls, 
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw, 
                                 const VarLabel* label)
{
  if (label == lb->gVelocityStarLabel) {
    exMomIntegrated(pg, patches, matls, old_dw, new_dw);
  }
}

void 
PenaltyContact::exMomIntegrated(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;

  int numMatls = d_sharedState->getNumMPMMatls();
  ASSERTEQ(numMatls, matls->size());

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  // Need access to all velocity fields at once, so store in
  // vectors of NCVariables
  std::vector<constNCdouble> gMass(numMatls);
  std::vector<constNCdouble> gVolume(numMatls);
  std::vector<constNCVector> gLSContactForce(numMatls);
  std::vector<NCVector>      gVelocity_star(numMatls);

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    for (int mat = 0; mat < matls->size(); mat++) {
      int matID = matls->get(mat);
      new_dw->get(gMass[mat],           lb->gMassLabel,           matID, patch, gnone, 0);
      new_dw->get(gVolume[mat],         lb->gVolumeLabel,         matID, patch, gnone, 0);
      new_dw->get(gLSContactForce[mat], lb->gLSContactForceLabel, matID, patch, gnone, 0);

      new_dw->getModifiable(gVelocity_star[mat], lb->gVelocityStarLabel, matID, patch);
    }

    if (flag->d_axisymmetric) {
      std::ostringstream warn;
      warn << "Penalty contact not implemented for axisymmetry\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      Vector centerOfMassVelocity(0.,0.,0.);
      double centerOfMassMass = 0.0; 
      for (int mat = 0; mat < matls->size(); mat++) {
        if (!d_matls.requested(mat)) continue;
        centerOfMassVelocity += gVelocity_star[mat][node] * gMass[mat][node];
        centerOfMassMass     += gMass[mat][node]; 
      }
      centerOfMassVelocity /= centerOfMassMass;

      // Loop over materials.  Only proceed if velocity field mass
      // is nonzero (not numerical noise) and the difference from
      // the centerOfMassVelocity is nonzero (More than one velocity
      // field is contributing to grid vertex).
      for (int mat = 0; mat < numMatls; mat++) {
        if (!d_matls.requested(mat)) continue;

        double mass = gMass[mat][node];
        if (gLSContactForce[mat][node].length2() > 0.0 &&
            !compare(mass,0.0)) {

          // Dv is the change in velocity due to penalty forces at tracers
          Vector Dv = (gLSContactForce[mat][node]/mass)*delT;
          double normalDv = Dv.length();
          Vector normal = Dv/(normalDv+1.e-100);

          // deltaVel is the difference in velocity of this material
          // relative to the centerOfMassVelocity
          Vector deltaVelocity = gVelocity_star[mat][node] - centerOfMassVelocity;
          double normalDeltaVel = Dot(deltaVelocity, normal);
          
          Vector normal_normaldV = normal*normalDeltaVel;
          Vector dV_normalDV = deltaVelocity - normal_normaldV;
          Vector surfaceTangent = dV_normalDV/(dV_normalDV.length()+1.e-100);

          double tangentDeltaVelocity = Dot(deltaVelocity, surfaceTangent);
          double frictionCoefficient = Min(d_mu, tangentDeltaVelocity/std::abs(normalDv));
          Dv -= surfaceTangent*frictionCoefficient*std::abs(normalDv);
          gVelocity_star[mat][node] +=Dv;
        } // if gLSContactForce>0
      } // matls
    } // nodeiterator
  } // patches
}
