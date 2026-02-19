/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2026 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/Contact/DEMSpecifiedVelocityContact.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace Uintah;

DEMSpecifiedVelocityContact::DEMSpecifiedVelocityContact(const ProcessorGroup* myworld,
                                                         const MaterialManagerP& mat_manager,
                                                         const MPMLabel* labels,
                                                         const MPMFlags* flags,
                                                         ProblemSpecP& ps)
  : Contact(myworld, mat_manager, labels, flags, ps)
{
  ps->getWithDefault("master_material", d_material, 0);
  d_matls.add(d_material);

  ps->get("filename", d_filename);
  readSpecifiedVelocityFile();
}

void
DEMSpecifiedVelocityContact::readSpecifiedVelocityFile()
{
  if (d_filename != "") {
    std::ifstream is(d_filename.c_str());
    if (!is) {
      std::ostringstream err;
      err << "**ERROR** Could not open DEM specified contact motion file " 
          << d_filename << " ";
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
              << "is not monotonically increasing";
          throw ProblemSetupException(err.str(), __FILE__, __LINE__);
        }
        d_vel_profile.push_back(std::make_pair(t1, Vector(vx, vy, vz)));
      }
      t0 = t1;
    }
    if (d_vel_profile.size() < 2) {
      std::ostringstream err;
      err << "**ERROR** DEM Specified contact: failed to generate valid velocity profile.";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }
  }
}

void
DEMSpecifiedVelocityContact::setContactMaterialAttributes()
{
  auto* mat =
    static_cast<MPMMaterial*>(d_mat_manager->getMaterial("MPM", d_material));
  if (!mat) {
    std::ostringstream err;
    err << "**ERROR** DEMSpecifiedVelocityContact: Material " << d_material
        << " not found.";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  if (!mat->isDiscrete()) {
    std::ostringstream err;
    err << "**ERROR** DEMSpecifiedVelocityContact: Material " << d_material
        << " is not a DEM (discrete) material.";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  mat->setIsRigid(true);
}

void
DEMSpecifiedVelocityContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type", "dem_specified");
  contact_ps->appendElement("filename", d_filename);
  contact_ps->appendElement("master_material", d_material);
  d_matls.outputProblemSpec(contact_ps);
}

void
DEMSpecifiedVelocityContact::addComputesAndRequires(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls,
                                                   const VarLabel* gVelocity_label)
{
  if (gVelocity_label == d_mpm_labels->gVelocityLabel) {
    return;
  }

  Task* t = scinew Task("DEMSpecifiedVelocityContact::exchangeMomentum",
                        this,
                        &DEMSpecifiedVelocityContact::exchangeMomentum,
                        gVelocity_label);

  t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->modifies(gVelocity_label, matls->getUnion());

  sched->addTask(t, patches, matls);
}

void
DEMSpecifiedVelocityContact::exchangeMomentum(const ProcessorGroup*,
                                              const PatchSubset* patches,
                                              [[maybe_unused]] const MaterialSubset* matls,
                                              DataWarehouse* old_dw,
                                              DataWarehouse* new_dw,
                                              const VarLabel* gVelocity_label)
{
  simTime_vartype simTime;
  old_dw->get(simTime, d_mpm_labels->simulationTimeLabel);
  double tcurr = simTime;

  Vector prescribed_vel(0, 0, 0);
  if (d_vel_profile.size() > 0) {
    prescribed_vel = findVelFromProfile(tcurr);
  }

  int matID_master = d_mat_manager->getMaterial("MPM", d_material)->getDWIndex();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    constNCdouble gMass_master;
    new_dw->get(gMass_master, d_mpm_labels->gMassLabel, matID_master, patch, Ghost::None, 0);

    NCVector gVelocity_star_master;
    new_dw->getModifiable(gVelocity_star_master, gVelocity_label, matID_master, patch);

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      if (gMass_master[node] > 1e-100) {
        gVelocity_star_master[node] = prescribed_vel;
      }
    }
  }
}

Vector
DEMSpecifiedVelocityContact::findVelFromProfile(double t) const
{
  auto iter = std::find_if(d_vel_profile.begin(), d_vel_profile.end(), [t](const std::pair<double, Vector>& data) {
    return data.first > t;
  });
  if (iter == d_vel_profile.begin()) {
    return iter->second;
  } else if (iter == d_vel_profile.end()) {
    return (iter - 1)->second;
  } else {
    double t_val = (iter->first - t) / (iter->first - (iter - 1)->first);
    Vector vel = (1.0 - t_val) * (iter - 1)->second + t_val * iter->second;
    return vel;
  }
}
