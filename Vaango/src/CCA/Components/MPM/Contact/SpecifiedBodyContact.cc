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
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
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

  IntVector defaultDir(0, 0, 1);
  ps->getWithDefault("direction", d_direction, defaultDir);
  if ((d_direction[0] < 0 || d_direction[0] > 1) ||
      (d_direction[1] < 0 || d_direction[1] > 1) ||
      (d_direction[2] < 0 || d_direction[2] > 1)) {
    std::ostringstream err;
    err << "**ERROR** The direction vector components in specified contact cannot"
        << " have values other than 0 or 1";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  ps->getWithDefault("master_material", d_material, 0);
  d_matls.add(d_material); // always need specified material

  ps->getWithDefault("master_material_is_rigid", d_rigid_master_material, true);

  d_vol_const = 0.;
  ps->get("volume_constraint", d_vol_const);

  ps->getWithDefault("normal_only", d_normalOnly, false);

  // read a list of values from a file
  ps->get("filename", d_filename);
  readSpecifiedVelocityFile();

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
SpecifiedBodyContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "specified");
  contact_ps->appendElement("filename", d_filename);
  contact_ps->appendElement("direction", d_direction);
  contact_ps->appendElement("master_material", d_material);
  contact_ps->appendElement("master_material_is_rigid", d_rigid_master_material);
  contact_ps->appendElement("stop_time", d_stop_time);
  contact_ps->appendElement("velocity_after_stop", d_vel_after_stop);
  contact_ps->appendElement("volume_constraint", d_vol_const);
  contact_ps->appendElement("normal_only", d_normalOnly);

  d_matls.outputProblemSpec(contact_ps);
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
  t->requires(Task::NewDW, lb->gVolumeLabel,             Ghost::None);
  t->requires(Task::NewDW, lb->gSurfNormLabel,           Ghost::None);
  t->requires(Task::OldDW, lb->NC_CCweightLabel, z_matl, Ghost::None);

  t->modifies(gVelocity_label, mss);

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

  int numMatls = d_sharedState->getNumMPMMatls();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  Vector imposed_velocity(0.0, 0.0, 0.0);
  const double tcurr = d_sharedState->getElapsedTime(); // FIXME: + dt ?
  if (tcurr > d_stop_time) {
    imposed_velocity = d_vel_after_stop;
  } else if (d_vel_profile.size() > 0) {
    imposed_velocity = findVelFromProfile(tcurr);
  }

  constNCdoubleArray gMass(numMatls), gVolume(numMatls);
  constNCVectorArray gVelocity(numMatls), gSurfNorm(numMatls);
  NCVectorArray      gVelocity_star(numMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    for (int m = 0; m < matls->size(); m++) {
      int dwi = matls->get(m);
      new_dw->get(gMass[m], lb->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gVolume[m], lb->gVolumeLabel, dwi, patch, gnone, 0);
      new_dw->get(gSurfNorm[m], lb->gSurfNormLabel, dwi, patch, gnone, 0);
      new_dw->getModifiable(gVelocity_star[m], gVelocity_label, dwi, patch);
    }

    if (d_normalOnly) {
      computeNormalBasedExchange(patch, old_dw, gMass, gVolume, gSurfNorm,
                                 imposed_velocity, gVelocity_star);
    } else {
      computeDirectionBasedExchange(patch, old_dw, gMass, gVolume,
                                    imposed_velocity, gVelocity_star);
    }
  } // end patch for
}

void
SpecifiedBodyContact::computeNormalBasedExchange(const Patch* patch,
                                                 DataWarehouse* old_dw,
                                                 constNCdoubleArray& gMass,
                                                 constNCdoubleArray& gVolume,
                                                 constNCVectorArray& gSurfNorm,
                                                 const Vector& imposed_velocity,
                                                 NCVectorArray& gVelocity_star)
{
  Ghost::GhostType gnone = Ghost::None;
  int numMatls = d_sharedState->getNumMPMMatls();

  Vector dx = patch->dCell();
  double cell_vol = dx.x() * dx.y() * dx.z();

  constNCdouble NC_CCweight;
  old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gnone, 0);

  for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
    IntVector node = *iter;

    // Determine nodal volume
    double totalNodalVol = 0.0;
    for (int mat = 0; mat < numMatls; mat++) {
      totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
    }

    // also updates material d_material to new velocity.
    if (d_rigid_master_material) {
      if (d_vel_profile.size() > 0) {
        gVelocity_star[d_material][node] = imposed_velocity;
      }
    }
    for (int mat = 0; mat < numMatls; mat++) { 
      Vector master_vel = imposed_velocity;
      if (d_rigid_master_material) {
        master_vel = gVelocity_star[d_material][node];
        if (mat == d_material) {
          continue; // compatibility with rigid motion, doesnt affect matl 0
        }
      }

      Vector new_vel = gVelocity_star[mat][node];
      Vector normal = gSurfNorm[d_material][node];
      Vector deltaV = new_vel - master_vel;
      double normalDeltaVel = Dot(normal, deltaV);
      if (normalDeltaVel < 0.0) {
        Vector normal_normaldV = normal * normalDeltaVel;
        new_vel = gVelocity_star[mat][node] - normal_normaldV;
      }

      if (!compare(gMass[d_material][node], 0.) &&
          (totalNodalVol / cell_vol) > d_vol_const) {
        gVelocity_star[mat][node] = new_vel;
      } 
      
      //std::cout << "After rigid contact: Node = " << c << " material = " << n
      //          << " gVel = " << gVelocity_star[n][c] << "\n";
    } // end for matls
  } // end for NodeIterator
}

void
SpecifiedBodyContact::computeDirectionBasedExchange(const Patch* patch,
                                                    DataWarehouse* old_dw,
                                                    constNCdoubleArray& gMass,
                                                    constNCdoubleArray& gVolume,
                                                    const Vector& imposed_velocity,
                                                    NCVectorArray& gVelocity_star)
{
  Ghost::GhostType gnone = Ghost::None;
  int numMatls = d_sharedState->getNumMPMMatls();

  Vector dx = patch->dCell();
  double cell_vol = dx.x() * dx.y() * dx.z();

  constNCdouble NC_CCweight;
  old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gnone, 0);

  for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
    IntVector node = *iter;

    // Determine nodal volume
    double totalNodalVol = 0.0;
    for (int mat = 0; mat < numMatls; mat++) {
      totalNodalVol += gVolume[mat][node] * 8.0 * NC_CCweight[node];
    }

    // also updates material d_material to new velocity.
    for (int mat = 0; mat < numMatls; mat++) { 
      Vector master_vel = imposed_velocity;
      if (d_rigid_master_material) {
        master_vel = gVelocity_star[d_material][node];
        if (mat == d_material) {
          continue; // compatibility with rigid motion, doesnt affect matl 0
        }
      }

      Vector new_vel = gVelocity_star[mat][node];
      if (mat == d_material || d_direction[0])
        new_vel.x(master_vel.x());
      if (mat == d_material || d_direction[1])
        new_vel.y(master_vel.y());
      if (mat == d_material || d_direction[2])
        new_vel.z(master_vel.z());

      if (!compare(gMass[d_material][node], 0.) &&
          (totalNodalVol / cell_vol) > d_vol_const) {
        gVelocity_star[mat][node] = new_vel;
      } 
      
      //std::cout << "After rigid contact: Node = " << c << " material = " << n
      //          << " gVel = " << gVelocity_star[n][c] << "\n";
    } // end for matls
  } // end for NodeIterator
}

// find velocity from table of values
Vector
SpecifiedBodyContact::findVelFromProfile(double t) const
{
  auto iter = std::find_if(d_vel_profile.begin(), d_vel_profile.end(),
                           [t](const std::pair<double, Vector>& data) {
                              return data.first > t;    
                           });
  if (iter == d_vel_profile.begin()) {
    return iter->second;
  } else if (iter == d_vel_profile.end()) {
    return (iter-1)->second;
  } else {
    double t_val = (iter->first - t)/(iter->first - (iter-1)->first);
    Vector vel = (1.0 - t_val) * (iter-1)->second + t_val * iter->second;
    return vel;
  }
}

