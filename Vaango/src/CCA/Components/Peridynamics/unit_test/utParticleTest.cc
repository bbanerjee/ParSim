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

#include <CCA/Components/Peridynamics/unit_test/utParticleTest.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Malloc/Allocator.h>

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
#include <CCA/Components/Peridynamics/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsDomainBoundCond.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/MPMInterpolators/LinearInterpolator.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/Matrix3.h>
#include <Core/Util/DebugStream.h>


#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/test/unit_test.hpp>

using namespace boost::unit_test;
using namespace Vaango;

using Uintah::DataWarehouse;
using Uintah::Ghost;
using Uintah::MaterialSubset;
using Uintah::ParticleSubset;
using Uintah::Patch;
using Uintah::PatchSet;
using Uintah::PatchSubset;
using Uintah::ProblemSpecP;
using Uintah::ProcessorGroup;
using Uintah::EmptyMaterial;
using Uintah::Task;
using Uintah::UintahParallelComponent;

utParticleTest::utParticleTest(const ProcessorGroup* myworld)
  : UintahParallelComponent(myworld)
{
  d_labels = scinew PeridynamicsLabel();
}

utParticleTest::~utParticleTest()
{
  delete d_labels;
}

void utParticleTest::problemSetup(const ProblemSpecP& params, 
                                  const ProblemSpecP& restart_prob_spec, 
                                  Uintah::GridP& /*grid*/,
                                  Uintah::MaterialManagerP& mat_manager)
{
  BOOST_CHECK_CLOSE(1.0, 1.0, 1.0e-12);
  BOOST_CHECK_CLOSE(1.0, 0.0, 1.0e-12);

  d_mat_manager = sharedState;
  dynamic_cast<Uintah::Scheduler*>(getPort("scheduler"))->setPositionVar(d_labels->pPositionLabel);
  ProblemSpecP pt1 = params->findBlock("ParticleTest1");
  pt1->getWithDefault("doOutput", d_doOutput, 0);
  pt1->getWithDefault("doGhostCells", d_numGhostCells , 0);
  d_mymat = scinew Uintah::EmptyMaterial();
  d_mat_manager->registerEmptyMaterial(d_mymat);
}
 
void utParticleTest::scheduleInitialize(const Uintah::LevelP& level,
			                            Uintah::SchedulerP& sched)
{
  Task* task = scinew Task("initialize", this, &utParticleTest::initialize);
  task->computes(d_labels->pPositionLabel);
  task->computes(d_labels->pMassLabel);
  task->computes(d_labels->pParticleIDLabel);
  sched->addTask(task, level->eachPatch(), d_mat_manager->allMaterials());
}
 
void utParticleTest::scheduleComputeStableTimestep(const Uintah::LevelP& level,
					                               Uintah::SchedulerP& sched)
{
  Task* task = scinew Task("computeStableTimestep",
			   this, &utParticleTest::computeStableTimestep);
  task->computes(d_mat_manager->get_delt_label(),level.get_rep());
  sched->addTask(task, level->eachPatch(), d_mat_manager->allMaterials());
}

void
utParticleTest::scheduleTimeAdvance( const Uintah::LevelP& level, Uintah::SchedulerP& sched)
{
  const Uintah::MaterialSet* matls = d_mat_manager->allMaterials();
  Task* task = scinew Task("timeAdvance", this, &utParticleTest::timeAdvance);

  // set this in problemSetup.  0 is no ghost cells, 1 is all with 1 ghost
  // atound-node, and 2 mixes them
  if (d_numGhostCells == 0) {
    task->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pPositionLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::None, 0);
  }
  
  else if (d_numGhostCells == 1) {
    task->requires(Task::OldDW, d_labels->pPositionLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::AroundNodes, 1);
  }
  else if (d_numGhostCells == 2) {
    task->requires(Task::OldDW, d_labels->pPositionLabel, Ghost::None, 0);
    task->requires(Task::OldDW, d_labels->pMassLabel, Ghost::AroundNodes, 1);
    task->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::None, 0);
  }

  task->computes(d_labels->pPositionLabel_preReloc);
  task->computes(d_labels->pMassLabel_preReloc);
  task->computes(d_labels->pParticleIDLabel_preReloc);
  sched->addTask(task, level->eachPatch(), d_mat_manager->allMaterials());
  
  d_mat_manager->d_particleState.clear();
  d_mat_manager->d_particleState_preReloc.clear();
  for (int m = 0; m < matls->size(); m++) {
    std::vector<const Uintah::VarLabel*> vars;
    std::vector<const Uintah::VarLabel*> vars_preReloc;

    vars.push_back(d_labels->pMassLabel);
    vars.push_back(d_labels->pParticleIDLabel);

    vars_preReloc.push_back(d_labels->pMassLabel_preReloc);
    vars_preReloc.push_back(d_labels->pParticleIDLabel_preReloc);

    d_mat_manager->d_particleState.push_back(vars);
    d_mat_manager->d_particleState_preReloc.push_back(vars_preReloc);
  }
  sched->scheduleParticleRelocation(level, 
                                    d_labels->pPositionLabel_preReloc,
                				    d_mat_manager->d_particleState_preReloc,
				                    d_labels->pPositionLabel, 
                                    d_mat_manager->d_particleState,
                				    d_labels->pParticleIDLabel, 
                                    matls);
}

void utParticleTest::computeStableTimestep(const ProcessorGroup* /*pg*/,
				     const PatchSubset* patches,
				     const MaterialSubset* /*matls*/,
				     DataWarehouse*,
				     DataWarehouse* new_dw)
{
  const Uintah::Level* level = getLevel(patches);
  new_dw->put(Uintah::delt_vartype(1), d_mat_manager->get_delt_label(), level);
}

void utParticleTest::initialize(const ProcessorGroup*,
			                    const PatchSubset* patches,
			                    const MaterialSubset* matls,
			                    DataWarehouse* /*old_dw*/, 
                                DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    Uintah::Point low = patch->cellPosition(patch->getCellLowIndex());
    Uintah::Point high = patch->cellPosition(patch->getCellHighIndex());
    for(int m = 0;m<matls->size();m++){
      srand(1);
      int numParticles = 10;
      int matl = matls->get(m);

      Uintah::ParticleVariable<Uintah::Point> px;
      Uintah::ParticleVariable<double> pmass;
      Uintah::ParticleVariable<Uintah::long64> pids;

      ParticleSubset* subset = new_dw->createParticleSubset(numParticles,matl,patch);
      new_dw->allocateAndPut(px,      d_labels->pPositionLabel,      subset);
      new_dw->allocateAndPut(pmass,   d_labels->pMassLabel,          subset);
      new_dw->allocateAndPut(pids,    d_labels->pParticleIDLabel,    subset);

      for (int i = 0; i < numParticles; i++) {
        Uintah::Point pos( (((float) rand()) / RAND_MAX * ( high.x() - low.x()-1) + low.x()),
          (((float) rand()) / RAND_MAX * ( high.y() - low.y()-1) + low.y()),
          (((float) rand()) / RAND_MAX * ( high.z() - low.z()-1) + low.z()));
        px[i] = pos;
        pids[i] = patch->getID()*numParticles+i;
        pmass[i] = ((float) rand()) / RAND_MAX * 10;
      }
    }
  }
}

void utParticleTest::timeAdvance(const ProcessorGroup*,
			const PatchSubset* patches,
			const MaterialSubset* matls,
			DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    for(int m = 0;m<matls->size();m++){
      int matl = matls->get(m);
      ParticleSubset* pset = old_dw->getParticleSubset(matl, patch);
      ParticleSubset* delset = scinew ParticleSubset(0,matl,patch);

      // Get the arrays of particle values to be changed
      Uintah::constParticleVariable<Uintah::Point> px;
      Uintah::ParticleVariable<Uintah::Point> pxnew;
      Uintah::constParticleVariable<double> pmass;
      Uintah::ParticleVariable<double> pmassnew;
      Uintah::constParticleVariable<Uintah::long64> pids;
      Uintah::ParticleVariable<Uintah::long64> pidsnew;

      old_dw->get(pmass, d_labels->pMassLabel,               pset);
      old_dw->get(px,    d_labels->pPositionLabel,           pset);
      old_dw->get(pids,  d_labels->pParticleIDLabel,         pset);

      new_dw->allocateAndPut(pmassnew, d_labels->pMassLabel_preReloc,       pset);
      new_dw->allocateAndPut(pxnew,    d_labels->pPositionLabel_preReloc,   pset);
      new_dw->allocateAndPut(pidsnew,  d_labels->pParticleIDLabel_preReloc, pset);

      // every timestep, move down the +x axis, and decay the mass a little bit
      for (int i = 0; i < pset->numParticles(); i++) {
        Uintah::Point pos( px[i].x() + .25, px[i].y(), px[i].z());
        pxnew[i] = pos;
        pidsnew[i] = pids[i];
        pmassnew[i] = pmass[i] *.9;
        if (d_doOutput)
          std::cout << " Patch " << patch->getID() << ": ID " 
               << pidsnew[i] << ", pos " << pxnew[i] 
               << ", mass " << pmassnew[i] << std::endl;
      }
      new_dw->deleteParticles(delset);
    }
  }
}
