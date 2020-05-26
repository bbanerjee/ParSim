/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/SpecifiedBodyContact.h>
#include <CCA/Components/MPM/MPMBoundCond.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/SimulationStateP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <vector>

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>
using std::cerr;

using namespace Uintah;

SpecifiedBodyContact::SpecifiedBodyContact(const ProcessorGroup* myworld,
                                           ProblemSpecP& ps,
                                           SimulationStateP& d_sS,
                                           MPMLabel* Mlb, MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  d_needNormals = true;

  // Constructor
  // read a list of values from a file
  ps->get("filename", d_filename);

  IntVector defaultDir(0, 0, 1);
  ps->getWithDefault("direction", d_direction, defaultDir);

  ps->getWithDefault("master_material", d_material, 0);
  d_matls.add(d_material); // always need specified material

  d_vol_const = 0.;
  ps->get("volume_constraint", d_vol_const);

  ps->getWithDefault("normal_only", d_NormalOnly, false);

  if (d_filename != "") {
    std::ifstream is(d_filename.c_str());
    if (!is) {
      throw ProblemSetupException("ERROR: opening MPM rigid motion file '" +
                                    d_filename +
                                    "'\nFailed to find profile file",
                                  __FILE__, __LINE__);
    }
    double t0(-1.e9);
    while (is) {
      double t1;
      double vx, vy, vz;
      is >> t1 >> vx >> vy >> vz;
      if (is) {
        if (t1 <= t0) {
          throw ProblemSetupException(
            "ERROR: profile file is not monotomically increasing", __FILE__,
            __LINE__);
        }
        d_vel_profile.push_back(
          std::pair<double, Vector>(t1, Vector(vx, vy, vz)));
      }
      t0 = t1;
    }
    if (d_vel_profile.size() < 2) {
      throw ProblemSetupException(
        "ERROR: Failed to generate valid velocity profile", __FILE__, __LINE__);
    }
  }

  // disable all changes after this time
  ps->getWithDefault("stop_time", d_stop_time,
                     std::numeric_limits<double>::max());
  ps->getWithDefault("velocity_after_stop", d_vel_after_stop, Vector(0, 0, 0));

  d_sharedState = d_sS;
  lb = Mlb;
  flag = MFlag;
  if (flag->d_8or27 == 8) {
    NGP = 1;
    NGN = 1;
  } else {
    NGP = 2;
    NGN = 2;
  }
}

SpecifiedBodyContact::~SpecifiedBodyContact()
{
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

  d_matls.outputProblemSpec(contact_ps);
}

// find velocity from table of values
Vector
SpecifiedBodyContact::findVelFromProfile(double t) const
{
  int smin = 0;
  int smax = (int)(d_vel_profile.size()) - 1;
  double tmin = d_vel_profile[0].first;
  double tmax = d_vel_profile[smax].first;
  if (t <= tmin) {
    return d_vel_profile[0].second;
  } else if (t >= tmax) {
    return d_vel_profile[smax].second;
  } else {
    // bisection search on table
    // could probably speed this up by keeping copy of last successful
    // search, and looking at that point and a couple to the right
    //
    while (smax > smin + 1) {
      int smid = (smin + smax) / 2;
      if (d_vel_profile[smid].first < t) {
        smin = smid;
      } else {
        smax = smid;
      }
    }
    double l = (d_vel_profile[smin + 1].first - d_vel_profile[smin].first);
    double xi = (t - d_vel_profile[smin].first) / l;
    double vx = xi * d_vel_profile[smin + 1].second[0] +
                (1 - xi) * d_vel_profile[smin].second[0];
    double vy = xi * d_vel_profile[smin + 1].second[1] +
                (1 - xi) * d_vel_profile[smin].second[1];
    double vz = xi * d_vel_profile[smin + 1].second[2] +
                (1 - xi) * d_vel_profile[smin].second[2];
    return Vector(vx, vy, vz);
  }
}

void
SpecifiedBodyContact::addComputesAndRequires(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls,
                                             const VarLabel* gVelocity_label)
{
  Task* t = scinew Task("SpecifiedBodyContact::exchangeMomentum", this,
                        &SpecifiedBodyContact::exchangeMomentum,
                        gVelocity_label);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);
  t->requires(Task::NewDW, lb->gMassLabel,               Ghost::None);
  //t->requires(Task::NewDW, lb->gInternalForceLabel,    Ghost::None);
  t->requires(Task::NewDW, lb->gVolumeLabel,             Ghost::None);
  t->requires(Task::OldDW, lb->NC_CCweightLabel, z_matl, Ghost::None);

  Ghost::GhostType gan = Ghost::AroundNodes;

  if (d_NormalOnly) {
    t->requires(Task::OldDW, lb->pXLabel, gan, NGP);
    t->requires(Task::OldDW, lb->pMassLabel, gan, NGP);
    t->requires(Task::OldDW, lb->pVolumeLabel, gan, NGP);
    t->requires(Task::OldDW, lb->pStressLabel, gan, NGP);
    t->requires(Task::OldDW, lb->pSizeLabel, gan, NGP);
    t->requires(Task::OldDW, lb->pDefGradLabel, gan, NGP);
    t->computes(lb->gSurfNormLabel);
  }

  t->modifies(gVelocity_label, mss);
  //t->computes(lb->RigidReactionForceLabel);

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
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
  Ghost::GhostType gan = Ghost::AroundNodes;

  int numMatls = d_sharedState->getNumMPMMatls();

  // Retrieve necessary data from DataWarehouse
  std::vector<constNCVariable<double>> gMass(numMatls);
  std::vector<NCVariable<Vector>>      gVelocity_star(numMatls);
  std::vector<constNCVariable<Vector>> gVelocity(numMatls);
  //std::vector<constNCVariable<Vector>> ginternalForce(numMatls);
  std::vector<constNCVariable<double>> gVolume(numMatls);
  NCVariable<Vector> gSurfNorm;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    Vector dx = patch->dCell();
    Vector oodx = {1.0/dx.x(), 1.0/dx.y(), 1.0/dx.z()};
    double cell_vol = dx.x() * dx.y() * dx.z();

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gnone, 0);

    for (int m = 0; m < matls->size(); m++) {
      int dwi = matls->get(m);
      new_dw->get(gMass[m],          lb->gMassLabel, dwi, patch, gnone, 0);
      //new_dw->get(ginternalForce[m], lb->gInternalForceLabel, dwi, patch, gnone,
      //            0);
      new_dw->get(gVolume[m], lb->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, dwi, patch);
    }

    // Compute the normals for the rigid material
    if (d_NormalOnly) {
      new_dw->allocateAndPut(gSurfNorm, lb->gSurfNormLabel, d_material, patch);
      gSurfNorm.initialize(Vector(0.0, 0.0, 0.0));

      ParticleSubset* pset =
        old_dw->getParticleSubset(d_material, patch, gan, NGP, lb->pXLabel);

      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;

      old_dw->get(pX,       lb->pXLabel,       pset);
      old_dw->get(pMass,    lb->pMassLabel,    pset);
      old_dw->get(pVolume,  lb->pVolumeLabel,  pset);
      old_dw->get(pSize,    lb->pSizeLabel,    pset);
      old_dw->get(pDefGrad, lb->pDefGradLabel, pset);

      auto interpolator = flag->d_interpolator->clone(patch);
      auto numInfluenceNodes = interpolator->size();
      vector<IntVector> ni(numInfluenceNodes);
      vector<double> S(numInfluenceNodes);
      vector<Vector> d_S(numInfluenceNodes);
      string interp_type = flag->d_interpolatorType;

      if (flag->d_axisymmetric) {

        for (auto idx : *pset) {

          interpolator->findCellAndShapeDerivatives(
            pX[idx], ni, d_S, pSize[idx], pDefGrad[idx]);
          double rho = pMass[idx] / pVolume[idx];

          int k = 0;
          for (auto& node : ni) {
            if (patch->containsNode(node)) {
              Vector G(d_S[k].x(), d_S[k].y(), 0.0);
              gSurfNorm[node] += rho * G;
            } 
            ++k;
          } // node for

        } // particle for

      } else {

        for (auto idx : *pset) {

          interpolator->findCellAndShapeDerivatives(
            pX[idx], ni, d_S, pSize[idx], pDefGrad[idx]);

          int k = 0;
          for (auto& node : ni) {
            if (patch->containsNode(node)) {
              Vector grad(d_S[k].x() * oodx[0], d_S[k].y() * oodx[1],
                          d_S[k].z() * oodx[2]);
              gSurfNorm[node] += pMass[idx] * grad;
            } 
          } // node for

        } // particle for

      } 

      MPMBoundCond bc;
      bc.setBoundaryCondition(patch, d_material, "Symmetric", gSurfNorm,
                              interp_type);

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        double length = gSurfNorm[c].length();
        if (length > 1.0e-15) {
          gSurfNorm[c] = gSurfNorm[c] / length;
        }
      }
    } // if(d_NormalOnly)

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    // set velocity to appropriate vel
    const double tcurr = d_sharedState->getElapsedTime(); // FIXME: + dt ?

    bool rigid_velocity = true;
    Vector requested_velocity(0.0, 0.0, 0.0);
    if (tcurr > d_stop_time) {
      rigid_velocity = false;
      requested_velocity = d_vel_after_stop;
    } else if (d_vel_profile.size() > 0) {
      rigid_velocity = false;
      requested_velocity = findVelFromProfile(tcurr);
    }

    //Vector reaction_force(0.0, 0.0, 0.0);

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      // Determine nodal volume
      double totalNodalVol = 0.0;
      for (int n = 0; n < numMatls; n++) {
        totalNodalVol += gVolume[n][c] * 8.0 * NC_CCweight[c];
      }

      // also updates material d_material to new velocity.
      for (int n = 0; n < numMatls; n++) { 
        Vector rigid_vel = requested_velocity;
        if (rigid_velocity) {
          rigid_vel = gVelocity_star[d_material][c];
          if (n == d_material) {
            continue; // compatibility with rigid motion, doesnt affect matl 0
          }
        }

        Vector new_vel(gVelocity_star[n][c]);
        if (d_NormalOnly) {
          Vector normal = gSurfNorm[c];
          double normalDeltaVel =
            Dot(normal, (gVelocity_star[n][c] - rigid_vel));
          if (normalDeltaVel < 0.0) {
            Vector normal_normaldV = normal * normalDeltaVel;
            new_vel = gVelocity_star[n][c] - normal_normaldV;
          }
        } else {
          new_vel = gVelocity_star[n][c];
          if (n == d_material || d_direction[0])
            new_vel.x(rigid_vel.x());
          if (n == d_material || d_direction[1])
            new_vel.y(rigid_vel.y());
          if (n == d_material || d_direction[2])
            new_vel.z(rigid_vel.z());
        }

        if (!compare(gMass[d_material][c], 0.) &&
            (totalNodalVol / cell_vol) > d_vol_const) {
          gVelocity_star[n][c] = new_vel;
          // Vector old_vel = gVelocity_star[n][c];
          // reaction_force += gMass[n][c]*(new_vel-old_vel)/delT;
          // reaction_force -= ginternalForce[n][c];
        } // if
        
        //std::cout << "After rigid contact: Node = " << c << " material = " << n
        //          << " gVel = " << gVelocity_star[n][c] << "\n";
      }   // for matls
    }     // for Node Iterator
    // new_dw->put(sumvec_vartype(reaction_force), lb->RigidReactionForceLabel);
  }
}

