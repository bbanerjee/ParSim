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

#include <CCA/Components/Models/HEChem/Common.h>
#include <CCA/Components/Models/HEChem/LightTime.h>

#include <CCA/Components/ICE/Core/BoundaryCond.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPMICE/Core/MPMICELabel.h>
#include <Core/Grid/Variables/PerPatchVars.h>

#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <iostream>

using namespace Uintah;
using namespace std;

#define d_SMALL_NUM 1e-100
#define d_TINY_RHO 1e-12

//__________________________________
//  setenv SCI_DEBUG "MODELS_NORMAL_COUT:+,MODELS_DOING_COUT:+"
//  MODELS_DOING_COUT:   dumps when tasks are scheduled and performed
static DebugStream cout_doing("MODELS_DOING_COUT", false);

LightTime::LightTime(const ProcessorGroup* myworld,
                     const MaterialManagerP& materialManager,
                     const ProblemSpecP& params)
  : HEChemModel(myworld, materialManager)
  , d_params(params)
{
  mymatls = 0;
  Ilb     = scinew ICELabel();
  //__________________________________
  //  diagnostic labels
  reactedFractionLabel =
    VarLabel::create("F", CCVariable<double>::getTypeDescription());

  delFLabel =
    VarLabel::create("delF", CCVariable<double>::getTypeDescription());

  mag_grad_Fr_Label =
    VarLabel::create("mag_grad_Fr", CCVariable<double>::getTypeDescription());
}

LightTime::~LightTime()
{
  delete Ilb;

  VarLabel::destroy(reactedFractionLabel);
  VarLabel::destroy(mag_grad_Fr_Label);
  VarLabel::destroy(delFLabel);

  if (mymatls && mymatls->removeReference()) {
    delete mymatls;
  }
}
//__________________________________
void
LightTime::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP model_ps = ps->appendChild("Model");
  model_ps->setAttribute("type", "LightTime");
  ProblemSpecP lt_ps = model_ps->appendChild("LightTime");

  lt_ps->appendElement("fromMaterial", std::string(matl0->getName()));
  lt_ps->appendElement("toMaterial", std::string(matl1->getName()));
  lt_ps->appendElement("starting_location", d_start_place);
  lt_ps->appendElement("direction_if_plane", d_direction);
  lt_ps->appendElement("D", d_D);
  lt_ps->appendElement("E0", d_E0);
  lt_ps->appendElement("react_mixed_cells", d_react_mixed_cells);
}
//__________________________________
void
LightTime::problemSetup(GridP&, [[maybe_unused]] const bool isRestart)
{
  ProblemSpecP lt_ps = d_params->findBlock("LightTime");
  if (!lt_ps) {
    throw ProblemSetupException(
      "LightTime: Couldn't find <LightTime> tag", __FILE__, __LINE__);
  }

  matl0 = d_materialManager->parseAndLookupMaterial(lt_ps, "fromMaterial");
  matl1 = d_materialManager->parseAndLookupMaterial(lt_ps, "toMaterial");

  lt_ps->require("starting_location", d_start_place);
  lt_ps->require("direction_if_plane", d_direction);
  lt_ps->require("D", d_D);
  lt_ps->require("E0", d_E0);
  lt_ps->getWithDefault("react_mixed_cells", d_react_mixed_cells, true);
  lt_ps->getWithDefault("AMR_Refinement_Criteria", d_refineCriteria, 1e100);

  // if point  ignition is desired, direction_if_plane = (0.,0.,0)
  // if planar ignition is desired, direction_if_plane is normal in
  // the direction of burning
  //__________________________________
  //  define the materialSet
  mymatls = scinew MaterialSet();

  vector<int> m;
  m.push_back(0); // needed for the pressure and NC_CCWeight
  m.push_back(matl0->getDWIndex());
  m.push_back(matl1->getDWIndex());

  mymatls->addAll_unique(m); // eliminate duplicate entries
  mymatls->addReference();
}
//______________________________________________________________________
//
void
LightTime::scheduleInitialize(SchedulerP& sched, const LevelP& level)
{
  cout_doing << "LightTime::scheduleInitialize " << endl;
  Task* t = scinew Task("LightTime::initialize", this, &LightTime::initialize);

  const MaterialSubset* react_matl = matl0->thisMaterial();

  t->computes(reactedFractionLabel, react_matl);

  sched->addTask(t, level->eachPatch(), mymatls);
}

//______________________________________________________________________
//
void
LightTime::initialize(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset*,
                      DataWarehouse*,
                      DataWarehouse* new_dw)
{
  cout_doing << "Doing Initialize \t\t\t\t\tLightTime" << endl;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int m0             = matl0->getDWIndex();

    CCVariable<double> Fr;
    new_dw->allocateAndPut(Fr, reactedFractionLabel, m0, patch);
    Fr.initialize(0.);
  }
}
//______________________________________________________________________
//
void
LightTime::scheduleComputeStableTimestep(SchedulerP&, const LevelP&)
{
  // None necessary...
}

//______________________________________________________________________
//
void
LightTime::scheduleComputeModelSources(SchedulerP& sched, const LevelP& level)
{
  Task* t = scinew Task(
    "LightTime::computeModelSources", this, &LightTime::computeModelSources);
  cout_doing << "LightTime::scheduleComputeModelSources " << endl;

  Ghost::GhostType gn              = Ghost::None;
  const MaterialSubset* react_matl = matl0->thisMaterial();
  const MaterialSubset* prod_matl  = matl1->thisMaterial();

  t->needs(Task::OldDW, Ilb->timeStepLabel);
  t->needs(Task::OldDW, Ilb->simulationTimeLabel);
  t->needs(Task::OldDW, Ilb->delTLabel, level.get_rep());
  //__________________________________
  // Products
  t->needs(Task::NewDW, Ilb->rho_CCLabel, prod_matl, gn);
  t->needs(Task::NewDW, Ilb->vol_frac_CCLabel, prod_matl, gn);

  //__________________________________
  // Reactants
  t->needs(Task::NewDW, Ilb->vol_frac_CCLabel, react_matl, gn);
  t->needs(Task::NewDW, Ilb->specificVolume_CCLabel, react_matl, gn);
  t->needs(Task::OldDW, Ilb->velocity_CCLabel, react_matl, gn);
  t->needs(Task::OldDW, Ilb->temperature_CCLabel, react_matl, gn);
  t->needs(Task::NewDW, Ilb->rho_CCLabel, react_matl, gn);
  t->needs(Task::NewDW, Ilb->specific_heatLabel, react_matl, gn);

  t->computes(reactedFractionLabel, react_matl);
  t->computes(delFLabel, react_matl);

  t->modifies(Ilb->modelMass_srcLabel);
  t->modifies(Ilb->modelMom_srcLabel);
  t->modifies(Ilb->modelEng_srcLabel);
  t->modifies(Ilb->modelVol_srcLabel);
  sched->addTask(t, level->eachPatch(), mymatls);
}

//______________________________________________________________________
//
void
LightTime::computeModelSources(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  timeStep_vartype timeStep;
  old_dw->get(timeStep, Ilb->timeStepLabel);

  bool isNotInitialTimestep = (timeStep > 0);

  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, Ilb->simulationTimeLabel);
  double simTime = simTimeVar;

  delt_vartype delT;
  old_dw->get(delT, Ilb->delTLabel, level);

  int m0 = matl0->getDWIndex();
  int m1 = matl1->getDWIndex();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    cout_doing << "Doing computeModelSources on patch " << patch->getID()
               << "\t\t\t\t  LightTime \tL-" << level->getIndex() << endl;
    CCVariable<double> mass_src_0, mass_src_1, mass_0;
    CCVariable<Vector> momentum_src_0, momentum_src_1;
    CCVariable<double> energy_src_0, energy_src_1;
    CCVariable<double> sp_vol_src_0, sp_vol_src_1;

    new_dw->getModifiable(mass_src_0, Ilb->modelMass_srcLabel, m0, patch);
    new_dw->getModifiable(momentum_src_0, Ilb->modelMom_srcLabel, m0, patch);
    new_dw->getModifiable(energy_src_0, Ilb->modelEng_srcLabel, m0, patch);
    new_dw->getModifiable(sp_vol_src_0, Ilb->modelVol_srcLabel, m0, patch);

    new_dw->getModifiable(mass_src_1, Ilb->modelMass_srcLabel, m1, patch);
    new_dw->getModifiable(momentum_src_1, Ilb->modelMom_srcLabel, m1, patch);
    new_dw->getModifiable(energy_src_1, Ilb->modelEng_srcLabel, m1, patch);
    new_dw->getModifiable(sp_vol_src_1, Ilb->modelVol_srcLabel, m1, patch);

    constCCVariable<double> vol_frac_rct, vol_frac_prd;
    constCCVariable<double> cv_reactant;
    constCCVariable<double> rctTemp, rctRho, rctSpvol, prodRho;
    constCCVariable<Vector> rctvel_CC;
    CCVariable<double> Fr;
    CCVariable<double> delF;

    Vector dx           = patch->dCell();
    double cell_vol     = dx.x() * dx.y() * dx.z();
    Ghost::GhostType gn = Ghost::None;

    //__________________________________
    // Reactant data
    old_dw->get(rctTemp, Ilb->temperature_CCLabel, m0, patch, gn, 0);
    old_dw->get(rctvel_CC, Ilb->velocity_CCLabel, m0, patch, gn, 0);
    new_dw->get(rctRho, Ilb->rho_CCLabel, m0, patch, gn, 0);
    new_dw->get(rctSpvol, Ilb->specificVolume_CCLabel, m0, patch, gn, 0);
    new_dw->get(vol_frac_rct, Ilb->vol_frac_CCLabel, m0, patch, gn, 0);
    new_dw->get(cv_reactant, Ilb->specific_heatLabel, m0, patch, gn, 0);
    new_dw->allocateAndPut(Fr, reactedFractionLabel, m0, patch);
    new_dw->allocateAndPut(delF, delFLabel, m0, patch);
    Fr.initialize(0.);
    delF.initialize(0.);

    //__________________________________
    // Product Data,
    new_dw->get(prodRho, Ilb->rho_CCLabel, m1, patch, gn, 0);
    new_dw->get(vol_frac_prd, Ilb->vol_frac_CCLabel, m1, patch, gn, 0);

    const Level* level = patch->getLevel();
    // double simTime = d_materialManager->getElapsedSimTime();
    double delta_L = 1.5 * pow(cell_vol, 1. / 3.) / d_D;
    //    double delta_L = 1.5*dx.x()/d_D;
    double A = d_direction.x();
    double B = d_direction.y();
    double C = d_direction.z();

    double x0 = d_start_place.x();
    double y0 = d_start_place.y();
    double z0 = d_start_place.z();

    double D     = -A * x0 - B * y0 - C * z0;
    double denom = 1.0;
    double plane = 0.;
    if (d_direction.length() > 0.0) {
      plane = 1.0;
      denom = sqrt(A * A + B * B + C * C);
    }

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      Point pos = level->getCellPosition(c);
      double dist_plane =
        Abs(A * pos.x() + B * pos.y() + C * pos.z() + D) / denom;
      double dist_straight = (pos - d_start_place).length();
      double dist          = dist_plane * plane + dist_straight * (1. - plane);

      double t_b    = dist / d_D;
      double VF_SUM = 0.;
      if (!d_react_mixed_cells) {
        VF_SUM = .99;
      }
      if ((vol_frac_rct[c] + vol_frac_prd[c]) > VF_SUM) {
        if (simTime >= t_b && rctRho[c] > d_TINY_RHO) {
          Fr[c] = (simTime - t_b) / delta_L;

          if (Fr[c] > .96) {
            Fr[c] = 1.0;
          }

          double Fr_old = prodRho[c] / (rctRho[c] + prodRho[c]);
          delF[c]       = Fr[c] - Fr_old;

          //__________________________________
          // Insert Burn Model Here
          double rctMass    = rctRho[c] * cell_vol;
          double prdMass    = prodRho[c] * cell_vol;
          double burnedMass = min(delF[c] * (prdMass + rctMass), rctMass);

          //__________________________________
          // conservation of mass, momentum and energy
          mass_src_0[c] -= burnedMass;
          mass_src_1[c] += burnedMass;

          Vector momX = rctvel_CC[c] * burnedMass;
          momentum_src_0[c] -= momX;
          momentum_src_1[c] += momX;

          double energyX      = cv_reactant[c] * rctTemp[c] * burnedMass;
          double releasedHeat = burnedMass * d_E0;
          energy_src_0[c] -= energyX;
          energy_src_1[c] += energyX + releasedHeat;

          double createdVolx = burnedMass * rctSpvol[c];
          sp_vol_src_0[c] -= createdVolx;
          sp_vol_src_1[c] += createdVolx;
        } // if (time to light it)
      }   // if cell only contains rct and prod
      if (rctRho[c] <= d_TINY_RHO) {
        Fr[c]   = 1.0;
        delF[c] = 0.0;
      } // reactant mass is already consumed
    }   // cell iterator

    //__________________________________
    //  set symetric BC
    setBC(mass_src_0,
          "set_if_sym_BC",
          patch,
          d_materialManager,
          m0,
          new_dw,
          isNotInitialTimestep);
    setBC(mass_src_1,
          "set_if_sym_BC",
          patch,
          d_materialManager,
          m1,
          new_dw,
          isNotInitialTimestep);
    setBC(delF,
          "set_if_sym_BC",
          patch,
          d_materialManager,
          m0,
          new_dw,
          isNotInitialTimestep);
    setBC(Fr,
          "set_if_sym_BC",
          patch,
          d_materialManager,
          m0,
          new_dw,
          isNotInitialTimestep);
  }
}
//______________________________________________________________________
//
void
LightTime::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  cout_doing << "LightTime::scheduleErrorEstimate \t\t\tL-"
             << coarseLevel->getIndex() << '\n';

  Task* t =
    scinew Task("LightTime::errorEstimate", this, &LightTime::errorEstimate);

  Ghost::GhostType gac             = Ghost::AroundCells;
  const MaterialSubset* react_matl = matl0->thisMaterial();

  t->needs(Task::NewDW, reactedFractionLabel, react_matl, gac, 1);

  t->computes(mag_grad_Fr_Label, react_matl);
  t->modifies(d_regridder->getRefineFlagLabel(),
              d_regridder->refineFlagMaterials());
  t->modifies(d_regridder->getRefinePatchFlagLabel(),
              d_regridder->refineFlagMaterials());

  sched->addTask(t, coarseLevel->eachPatch(), mymatls);
}
/*_____________________________________________________________________
 Function~  PassiveScalar::errorEstimate--
______________________________________________________________________*/
void
LightTime::errorEstimate(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse*,
                         DataWarehouse* new_dw)
{
  cout_doing << "Doing errorEstimate \t\t\t\t\t LightTime" << endl;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    Ghost::GhostType gac             = Ghost::AroundCells;
    const VarLabel* refineFlagLabel  = d_regridder->getRefineFlagLabel();
    const VarLabel* refinePatchLabel = d_regridder->getRefinePatchFlagLabel();

    CCVariable<int> refineFlag;
    new_dw->getModifiable(refineFlag, refineFlagLabel, 0, patch);

    PerPatch<PatchFlagP> refinePatchFlag;
    new_dw->get(refinePatchFlag, refinePatchLabel, 0, patch);

    constCCVariable<double> Fr;
    CCVariable<double> mag_grad_Fr;
    int m0 = matl0->getDWIndex();

    new_dw->get(Fr, reactedFractionLabel, m0, patch, gac, 1);
    new_dw->allocateAndPut(mag_grad_Fr, mag_grad_Fr_Label, m0, patch);
    mag_grad_Fr.initialize(0.0);

    //__________________________________
    // compute gradient
    Vector dx = patch->dCell();

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      Vector grad_Fr;
      for (int dir = 0; dir < 3; dir++) {
        IntVector r   = c;
        IntVector l   = c;
        double inv_dx = 0.5 / dx[dir];
        r[dir] += 1;
        l[dir] -= 1;
        grad_Fr[dir] = (Fr[r] - Fr[l]) * inv_dx;
      }
      mag_grad_Fr[c] = grad_Fr.length();
    }
    //__________________________________
    // set refinement flag
    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      if (mag_grad_Fr[c] > d_refineCriteria) {
        refineFlag[c] = true;
        refinePatch->set();
      }
    }
  } // patches
}
