/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#include <CCA/Components/Models/FluidsBased/TestModel.h>

#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/ProblemSpec/ProblemSpec.h>

using namespace Uintah;
using namespace std;

TestModel::TestModel(const ProcessorGroup* myworld,
                     const MaterialManagerP& materialManager,
                     const ProblemSpecP& params)

  : FluidsBasedModel(myworld, materialManager)
  , d_params(params)
{
  mymatls           = 0;
  Ilb               = scinew ICELabel();
  MIlb              = scinew MPMICELabel();
  totalMassXLabel   = 0;
  totalIntEngXLabel = 0;
}

TestModel::~TestModel()
{
  delete Ilb;
  delete MIlb;
  if (mymatls && mymatls->removeReference()) {
    delete mymatls;
  }

  if (0 != totalMassXLabel) {
    VarLabel::destroy(totalMassXLabel);
  }

  if (0 != totalIntEngXLabel) {
    VarLabel::destroy(totalIntEngXLabel);
  }
}

//______________________________________________________________________
void
TestModel::problemSetup(GridP&, [[maybe_unused]] const bool isRestart)
{
  ProblemSpecP test_ps = d_params->findBlock("Test");
  if (!test_ps) {
    throw ProblemSetupException(
      "TestModel: Couldn't find <Test> tag", __FILE__, __LINE__);
  }

  matl0 = d_materialManager->parseAndLookupMaterial(test_ps, "fromMaterial");
  matl1 = d_materialManager->parseAndLookupMaterial(test_ps, "toMaterial");

  test_ps->require("rate", d_rate);
  test_ps->getWithDefault("startTime", d_startTime, 0.0);

  vector<int> m(2);
  m[0]    = matl0->getDWIndex();
  m[1]    = matl1->getDWIndex();
  mymatls = scinew MaterialSet();
  mymatls->addAll(m);
  mymatls->addReference();

  // What flavor of matl it is.
  Material* matl        = d_materialManager->getMaterial(m[0]);
  ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
  MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
  if (mpm_matl) {
    d_is_mpm_matl = true;
    d_matl        = mpm_matl;
  }
  if (ice_matl) {
    d_is_mpm_matl = false;
    d_matl        = ice_matl;
  }

  totalMassXLabel =
    VarLabel::create("totalMassExchanged", sum_vartype::getTypeDescription());

  totalIntEngXLabel =
    VarLabel::create("totalIntEngExchanged", sum_vartype::getTypeDescription());
}

//______________________________________________________________________
void
TestModel::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP model_ps = ps->appendChild("Model");
  model_ps->setAttribute("type", "Test");

  ProblemSpecP test_ps = model_ps->appendChild("Test");
  test_ps->appendElement("fromMaterial", std::string(matl0->getName()));
  test_ps->appendElement("toMaterial", std::string(matl1->getName()));
  test_ps->appendElement("startTime", d_startTime);
  test_ps->appendElement("rate", d_rate);
}

//______________________________________________________________________
void
TestModel::scheduleInitialize(SchedulerP&, [[maybe_unused]] const LevelP& level)
{
  // None necessary...
}

//______________________________________________________________________
void
TestModel::scheduleComputeStableTimestep(SchedulerP&, const LevelP&)
{
  // None necessary...
}

//__________________________________
void
TestModel::scheduleComputeModelSources(SchedulerP& sched, const LevelP& level)
{
  Task* t = scinew Task(
    "TestModel::computeModelSources", this, &TestModel::computeModelSources);
  t->modifies(Ilb->modelMass_srcLabel);
  t->modifies(Ilb->modelMom_srcLabel);
  t->modifies(Ilb->modelEng_srcLabel);
  t->modifies(Ilb->modelVol_srcLabel);
  Ghost::GhostType gn = Ghost::None;

  Task::WhichDW DW;
  Task::WhichDW NDW = Task::NewDW;
  if (d_is_mpm_matl) { // MPM (pull data from newDW)
    DW = Task::NewDW;
    t->needs(DW, MIlb->cMassLabel, matl0->thisMaterial(), gn);
  } else {
    DW = Task::OldDW; // ICE (pull data from old DW)
    t->needs(DW, Ilb->rho_CCLabel, matl0->thisMaterial(), gn);
    t->needs(NDW, Ilb->specific_heatLabel, matl0->thisMaterial(), gn);
  }
  // All matls
  t->needs(DW, Ilb->velocity_CCLabel, matl0->thisMaterial(), gn);
  t->needs(DW, Ilb->temperature_CCLabel, matl0->thisMaterial(), gn);
  t->needs(NDW, Ilb->specificVolume_CCLabel, matl0->thisMaterial(), gn);

  t->computes(TestModel::totalMassXLabel);
  t->computes(TestModel::totalIntEngXLabel);

  t->needs(Task::OldDW, Ilb->delTLabel, level.get_rep());
  t->needs(Task::OldDW, Ilb->simulationTimeLabel);
  sched->addTask(t, level->eachPatch(), mymatls);
}

//__________________________________
void
TestModel::computeModelSources(const ProcessorGroup*,
                               const PatchSubset* patches,
                               [[maybe_unused]] const MaterialSubset* matls,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, Ilb->simulationTimeLabel);
  double simTime = simTimeVar;

  delt_vartype delT;
  old_dw->get(delT, Ilb->delTLabel, getLevel(patches));
  double dt = delT;

  ASSERT(matls->size() == 2);
  int m0 = matl0->getDWIndex();
  int m1 = matl1->getDWIndex();

  double totalMassX   = 0.0;
  double totalIntEngX = 0.0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    CCVariable<double> mass_src_0, mass_src_1, mass_0, cv;
    CCVariable<Vector> mom_src_0, mom_src_1;
    CCVariable<double> eng_src_0, eng_src_1;
    CCVariable<double> sp_vol_src_0, sp_vol_src_1;

    new_dw->allocateTemporary(cv, patch);
    new_dw->getModifiable(mass_src_0, Ilb->modelMass_srcLabel, m0, patch);
    new_dw->getModifiable(mom_src_0, Ilb->modelMom_srcLabel, m0, patch);
    new_dw->getModifiable(eng_src_0, Ilb->modelEng_srcLabel, m0, patch);
    new_dw->getModifiable(sp_vol_src_0, Ilb->modelVol_srcLabel, m0, patch);

    new_dw->getModifiable(mass_src_1, Ilb->modelMass_srcLabel, m1, patch);
    new_dw->getModifiable(mom_src_1, Ilb->modelMom_srcLabel, m1, patch);
    new_dw->getModifiable(eng_src_1, Ilb->modelEng_srcLabel, m1, patch);
    new_dw->getModifiable(sp_vol_src_1, Ilb->modelVol_srcLabel, m1, patch);

    //__________________________________
    //  Compute the mass and specific heat of matl 0
    new_dw->allocateTemporary(mass_0, patch);
    Vector dx     = patch->dCell();
    double volume = dx.x() * dx.y() * dx.z();
    DataWarehouse* dw;
    Ghost::GhostType gn = Ghost::None;

    if (d_is_mpm_matl) {
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(d_matl);
      dw                    = new_dw; // MPM  (Just grab it)
      constCCVariable<double> cmass;
      dw->get(cmass, MIlb->cMassLabel, m0, patch, gn, 0);
      mass_0.copyData(cmass);

      cv.initialize(mpm_matl->getSpecificHeat());
    } else {
      dw = old_dw; // ICE   (compute it from the density)
      constCCVariable<double> rho_tmp, cv_ice;
      old_dw->get(rho_tmp, Ilb->rho_CCLabel, m0, patch, gn, 0);
      new_dw->get(cv_ice, Ilb->specific_heatLabel, m0, patch, gn, 0);

      cv.copyData(cv_ice);

      for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
           iter++) {
        mass_0[*iter] = rho_tmp[*iter] * volume;
      }
    }

    constCCVariable<Vector> vel_0;  // MPM  pull from new_dw
    constCCVariable<double> temp_0; // ICE  pull from old_dw
    constCCVariable<double> sp_vol_0;
    dw->get(vel_0, Ilb->velocity_CCLabel, m0, patch, gn, 0);
    dw->get(temp_0, Ilb->temperature_CCLabel, m0, patch, gn, 0);
    new_dw->get(sp_vol_0, Ilb->specificVolume_CCLabel, m0, patch, gn, 0);

    double trate = d_rate * dt;
    if (trate > 1) {
      trate = 1;
    }
    //__________________________________
    //  Do some work

    // double simTime  = d_materialManager->getElapsedSimTime();

    if (simTime >= d_startTime) {
      for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
           iter++) {
        IntVector c  = *iter;
        double massx = mass_0[c] * trate;
        mass_src_0[c] -= massx;
        mass_src_1[c] += massx;

        Vector momx = vel_0[c] * massx;
        mom_src_0[c] -= momx;
        mom_src_1[c] += momx;

        double energyx = temp_0[c] * massx * cv[c];
        eng_src_0[c] -= energyx;
        eng_src_1[c] += energyx;

        double vol_sourcex = massx * sp_vol_0[c];
        sp_vol_src_0[c] -= vol_sourcex;
        sp_vol_src_1[c] += vol_sourcex;

        totalMassX += massx;
        totalIntEngX += energyx;
      }
    }

    new_dw->put(sum_vartype(totalMassX), TestModel::totalMassXLabel);
    new_dw->put(sum_vartype(totalIntEngX), TestModel::totalIntEngXLabel);
  }
}
//______________________________________________________________________

void
TestModel::scheduleModifyThermoTransportProperties(SchedulerP&,
                                                   const LevelP&,
                                                   const MaterialSet*)
{
  // do nothing
}
void
TestModel::computeSpecificHeat(CCVariable<double>&,
                               const Patch*,
                               DataWarehouse*,
                               const int)
{
  // do nothing
}
//______________________________________________________________________
//
void
TestModel::scheduleErrorEstimate(const LevelP&, SchedulerP&)
{
  // Not implemented yet
}
//__________________________________
void
TestModel::scheduleTestConservation(SchedulerP&, const PatchSet*)
{
  // Not implemented yet
}
