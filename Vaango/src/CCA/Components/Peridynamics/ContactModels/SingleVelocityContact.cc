/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/ContactModels/SingleVelocityContact.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsDomainBoundCond.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>

#include <vector>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>

#include <Core/Grid/Grid.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace Vaango;

using Uintah::ProcessorGroup;
using Uintah::ProblemSpecP;
using Uintah::SchedulerP;
using Uintah::PatchSubset;
using Uintah::MaterialSubset;
using Uintah::PatchSet;
using Uintah::MaterialSet;
using Uintah::DataWarehouse;
using Uintah::MaterialManagerP;

using Uintah::Task;
using Uintah::Patch;
using Uintah::Ghost;
using Uintah::NodeIterator;

using Uintah::delt_vartype;

using Uintah::NCVariable;
using Uintah::constNCVariable;
using Uintah::Vector;

SingleVelocityContact::SingleVelocityContact(const ProcessorGroup* myworld,
                                             ProblemSpecP& ps, 
                                             MaterialManagerP& ss, 
                                             PeridynamicsLabel* labels,
                                             PeridynamicsFlags* flags)
  : ContactModelBase(myworld, labels, flags, ps)
{
  d_mat_manager = ss;
}

SingleVelocityContact::~SingleVelocityContact()
{
}

void 
SingleVelocityContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("ContactModel");
  contact_ps->appendElement("type", "single_velocity");
  d_bodiesThatCanInteract.outputProblemSpec(contact_ps);
}

void 
SingleVelocityContact::addComputesAndRequiresInterpolated(SchedulerP & sched,
                                                          const PatchSet* patches,
                                                          const MaterialSet* ms)
{
  Task * t = scinew Task("SingleVelocityContact::exchangeMomentumInterpolated", 
                      this, &SingleVelocityContact::exchangeMomentumInterpolated);
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::NewDW, d_labels->gMassLabel, Ghost::None);

  t->modifies(d_labels->gVelocityLabel, mss);

  sched->addTask(t, patches, ms);
}

void 
SingleVelocityContact::exchangeMomentumInterpolated(const ProcessorGroup*,
                                                    const PatchSubset* patches,
                                                    const MaterialSubset* matls,
                                                    DataWarehouse*,
                                                    DataWarehouse* new_dw)
{
  std::string interp_type = "linear";

  // Check that only peridynamics bodies are being considered in this simulation
  // **TODO** Add the possibility of MPM and Peridynamics materials interacting
  int numBodies = d_mat_manager->getNumMaterials("Peridynamics");
  ASSERTEQ(numBodies, matls->size());

  for (int p=0; p<patches->size(); p++) {

    const Patch* patch = patches->get(p);
    Vector centerOfMassVelocity(0.0,0.0,0.0);

    // Retrieve necessary data from DataWarehouse
    std::vector<constNCVariable<double> > gMass(numBodies);
    std::vector<NCVariable<Vector> > gVelocity(numBodies);
    for (int m=0; m<matls->size(); m++) {
      int matlIndex = matls->get(m);
      new_dw->get(gMass[m], d_labels->gMassLabel, matlIndex, patch, Ghost::None, 0);
      new_dw->getModifiable(gVelocity[m], d_labels->gVelocityLabel, matlIndex, patch);
    }

    for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {

      IntVector c = *iter;

      Vector centerOfMassMom(0,0,0);
      double centerOfMassMass=0.0;

      for(int n = 0; n < numBodies; n++){
        if(d_bodiesThatCanInteract.requested(n)) {
          centerOfMassMom += gVelocity[n][c] * gMass[n][c];
          centerOfMassMass += gMass[n][c]; 
        }
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity = centerOfMassMom/centerOfMassMass;
      for(int n = 0; n < numBodies; n++) {
        if(d_bodiesThatCanInteract.requested(n)) {
          gVelocity[n][c] = centerOfMassVelocity;
        }
      }
    }

    for(int m=0;m<matls->size();m++){
      int matlIndex = matls->get(m);
      PeridynamicsDomainBoundCond bc;
      bc.setBoundaryCondition(patch, matlIndex, "Symmetric",gVelocity[m],interp_type);
    }
  }
}

void 
SingleVelocityContact::addComputesAndRequiresIntegrated(SchedulerP & sched,
                                                        const PatchSet* patches,
                                                        const MaterialSet* ms) 
{
  Task * t = scinew Task("SingleVelocityContact::exchangeMomentumIntegrated", 
                      this, &SingleVelocityContact::exchangeMomentumIntegrated);
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::OldDW, d_labels->delTLabel);    
  t->requires(Task::NewDW, d_labels->gMassLabel, Ghost::None);

  t->modifies(d_labels->gVelocityStarLabel, mss);

  sched->addTask(t, patches, ms);
}

void 
SingleVelocityContact::exchangeMomentumIntegrated(const ProcessorGroup*,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset* matls,
                                                  DataWarehouse* old_dw,
                                                  DataWarehouse* new_dw)
{
  int numBodies = d_mat_manager->getNumMaterials("Peridynamics");
  ASSERTEQ(numBodies, matls->size());

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);

    Vector zero(0.0,0.0,0.0);
    Vector centerOfMassVelocity(0.0,0.0,0.0);
    Vector centerOfMassMom(0.0,0.0,0.0);
    double centerOfMassMass;

    // Retrieve necessary data from DataWarehouse
    std::vector<constNCVariable<double> > gMass(numBodies);
    std::vector<NCVariable<Vector> > gVelocity_star(numBodies);

    for(int m=0; m < matls->size();m++){
     int matlIndex = matls->get(m);
     new_dw->get(gMass[m],d_labels->gMassLabel, matlIndex, patch, Ghost::None, 0);
     new_dw->getModifiable(gVelocity_star[m],d_labels->gVelocityStarLabel, matlIndex,patch);
    }

    delt_vartype delT;
    old_dw->get(delT, d_labels->delTLabel, getLevel(patches));
    
    for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++){
      IntVector c = *iter;

      centerOfMassMom=zero;
      centerOfMassMass=0.0; 
      for(int  n = 0; n < numBodies; n++){
        if(d_bodiesThatCanInteract.requested(n)) {
          centerOfMassMom += gVelocity_star[n][c] * gMass[n][c];
          centerOfMassMass += gMass[n][c]; 
        }
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity = centerOfMassMom/centerOfMassMass;
      for(int  n = 0; n < numBodies; n++){
        if(d_bodiesThatCanInteract.requested(n)) {
          //Vector dvdt = (centerOfMassVelocity - gVelocity_star[n][c])/delT;
          gVelocity_star[n][c] = centerOfMassVelocity;
        }
      }
    }
  }
}

