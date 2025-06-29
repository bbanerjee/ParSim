/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#include <CCA/Components/ICE/ICE.h>

#include <CCA/Components/ICE/Core/BoundaryCond.h>
#include <CCA/Components/ICE/CustomBCs/C_BC_driver.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/Utils.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/MiscMath.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Util/DebugStream.h>

#include <sci_defs/hypre_defs.h>

#ifdef HAVE_HYPRE
#include <CCA/Components/Solvers/HypreSolver.h>
#endif

#include <cmath>

using namespace Uintah;

static DebugStream cout_norm("ICE_NORMAL_COUT", false);
static DebugStream cout_doing("ICE_DOING_COUT", false);
static DebugStream cout_dbg("IMPICE_DBG", false);

/*___________________________________________________________________
 Function~  ICE::scheduleSetupMatrix--
_____________________________________________________________________*/
void
ICE::scheduleSetupMatrix(SchedulerP& sched,
                         const LevelP& level,
                         const PatchSet* patches,
                         const MaterialSubset* one_matl,
                         const MaterialSet* all_matls)
{
  Task* t;
  Ghost::GhostType gaf          = Ghost::AroundFaces;
  Ghost::GhostType gn           = Ghost::None;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  const MaterialSubset* press_matl = one_matl;

  Task::WhichDW whichDW = Task::OldDW;
  //__________________________________
  //  Form the matrix
  printSchedule(level, cout_doing, "ICE::scheduleSetupMatrix");

  t = scinew Task("ICE::setupMatrix", this, &ICE::setupMatrix);
  t->needs(Task::ParentOldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(whichDW, d_ice_labels->sp_volX_FCLabel, gaf, 1);
  t->needs(whichDW, d_ice_labels->sp_volY_FCLabel, gaf, 1);
  t->needs(whichDW, d_ice_labels->sp_volZ_FCLabel, gaf, 1);
  t->needs(whichDW, d_ice_labels->vol_fracX_FCLabel, gaf, 1);
  t->needs(whichDW, d_ice_labels->vol_fracY_FCLabel, gaf, 1);
  t->needs(whichDW, d_ice_labels->vol_fracZ_FCLabel, gaf, 1);
  t->needs(
    Task::ParentNewDW, d_ice_labels->sumKappaLabel, one_matl, oims, gn, 0);

  t->computes(d_ice_labels->matrixLabel, one_matl, oims);
  t->computes(d_ice_labels->imp_delPLabel, press_matl, oims);
  sched->addTask(t, patches, all_matls);
}

/*___________________________________________________________________
 Function~  ICE::scheduleSetupRHS--
_____________________________________________________________________*/
void
ICE::scheduleSetupRHS(SchedulerP& sched,
                      const PatchSet* patches,
                      const MaterialSubset* one_matl,
                      const MaterialSet* all_matls,
                      bool insideOuterIterLoop,
                      const string& computes_or_modifies)
{
  Task* t;
  Ghost::GhostType gac          = Ghost::AroundCells;
  Ghost::GhostType gn           = Ghost::None;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  printSchedule(patches, cout_doing, "ICE::scheduleSetupRHS");

  t = scinew Task("ICE::setupRHS",
                  this,
                  &ICE::setupRHS,
                  insideOuterIterLoop,
                  computes_or_modifies);

  Task::WhichDW pNewDW;
  Task::WhichDW pOldDW;
  if (insideOuterIterLoop) {
    pNewDW = Task::ParentNewDW;
    pOldDW = Task::ParentOldDW;
  } else {
    pNewDW = Task::NewDW;
    pOldDW = Task::OldDW;
  }

  const MaterialSubset* press_matl = one_matl;
  t->needs(pOldDW, d_ice_labels->delTLabel, getLevel(patches));

  if (d_models.size() > 0) {
    t->needs(pNewDW, d_ice_labels->modelMass_srcLabel, gn, 0);
  }

  t->needs(pNewDW, d_ice_labels->sumKappaLabel, press_matl, oims, gn, 0);
  t->needs(pNewDW, d_ice_labels->specificVolume_CCLabel, gn, 0);
  t->needs(pNewDW, d_ice_labels->speedSound_CCLabel, gn, 0);
  t->needs(pNewDW, d_ice_labels->vol_frac_CCLabel, gac, 2);
  t->needs(Task::NewDW, d_ice_labels->uvel_FCMELabel, gac, 2);
  t->needs(Task::NewDW, d_ice_labels->vvel_FCMELabel, gac, 2);
  t->needs(Task::NewDW, d_ice_labels->wvel_FCMELabel, gac, 2);
  t->needs(
    Task::NewDW, d_ice_labels->sum_imp_delPLabel, press_matl, oims, gn, 0);

  t->computes(d_ice_labels->vol_fracX_FCLabel);
  t->computes(d_ice_labels->vol_fracY_FCLabel);
  t->computes(d_ice_labels->vol_fracZ_FCLabel);
  t->computes(d_ice_labels->term2Label, one_matl, oims);

  if (isAMR()) { // compute refluxing variables if using AMR
    t->computes(d_ice_labels->vol_frac_X_FC_fluxLabel);
    t->computes(d_ice_labels->vol_frac_Y_FC_fluxLabel);
    t->computes(d_ice_labels->vol_frac_Z_FC_fluxLabel);
  }

  if (computes_or_modifies == "computes") {
    t->computes(d_ice_labels->rhsLabel, one_matl, oims);
  }
  if (computes_or_modifies == "modifies") {
    t->modifies(d_ice_labels->rhsLabel, one_matl, oims);
  }

  t->computes(VarLabel::find(abortTimestep_name));
  t->computes(VarLabel::find(recomputeTimestep_name));

  sched->addTask(t, patches, all_matls);
}

/*___________________________________________________________________
 Function~  ICE::scheduleCompute_maxRHS--
_____________________________________________________________________*/
void
ICE::scheduleCompute_maxRHS(SchedulerP& sched,
                            const LevelP& level,
                            const MaterialSubset* one_matl,
                            const MaterialSet* allMatls)
{

  printSchedule(level, cout_doing, "ICE::scheduleCompute_maxRHS");

  Task* t;
  t = scinew Task("ICE::compute_maxRHS", this, &ICE::compute_maxRHS);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  t->needs(
    Task::NewDW, d_ice_labels->rhsLabel, one_matl, oims, Ghost::None, 0);
  t->computes(d_ice_labels->max_RHSLabel);

  sched->addTask(t, level->eachPatch(), allMatls);
}

/*___________________________________________________________________
 Function~  ICE::scheduleUpdatePressure--
_____________________________________________________________________*/
void
ICE::scheduleUpdatePressure(SchedulerP& sched,
                            const LevelP& level,
                            const PatchSet* patches,
                            const MaterialSubset* ice_matls,
                            const MaterialSubset* /*mpm_matls*/,
                            const MaterialSubset* press_matl,
                            const MaterialSet* all_matls)
{
  Task* t;
  Ghost::GhostType gn           = Ghost::None;
  Ghost::GhostType gac          = Ghost::AroundCells;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  //__________________________________
  // update the pressure
  printSchedule(level, cout_doing, "ICE::scheduleUpdatePressure");

  t = scinew Task("ICE::updatePressure", this, &ICE::updatePressure);

  t->needs(Task::ParentOldDW, d_ice_labels->timeStepLabel);
  t->needs(Task::ParentOldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(
    Task::ParentNewDW, d_ice_labels->press_equil_CCLabel, press_matl, oims, gn);
  t->needs(Task::ParentNewDW, d_ice_labels->specificVolume_CCLabel, gn);
  t->needs(
    Task::OldDW, d_ice_labels->sum_imp_delPLabel, press_matl, oims, gn);

  // for setting the boundary conditions
  // you need gac,1 for funky multilevel patch configurations
  t->modifies(d_ice_labels->imp_delPLabel, press_matl, oims);
  if (level->getIndex() > 0) {
    t->needs(Task::NewDW,
                d_ice_labels->imp_delPLabel,
                0,
                Task::CoarseLevel,
                press_matl,
                oims,
                gac,
                1);
  }

  t->computes(d_ice_labels->sum_imp_delPLabel, press_matl, oims);

  computesRequires_CustomBCs(t,
                             "imp_update_press_CC",
                             d_ice_labels.get(),
                             ice_matls,
                             d_BC_globalVars.get());

  t->computes(d_ice_labels->press_CCLabel, press_matl, oims);
  sched->addTask(t, patches, all_matls);
}

/*___________________________________________________________________
 Function~  ICE::--scheduleRecomputeVel_FC
_____________________________________________________________________*/
void
ICE::scheduleRecomputeVel_FC(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSubset* ice_matls,
                             const MaterialSubset* mpm_matls,
                             const MaterialSubset* press_matl,
                             const MaterialSet* all_matls,
                             bool recursion)
{
  Task* t = 0;
  printSchedule(patches, cout_doing, "ICE::scheduleRecomputeVel_FC");

  t = scinew Task(
    "ICE::scheduleUpdateVel_FC", this, &ICE::updateVel_FC, recursion);

  Ghost::GhostType gac          = Ghost::AroundCells;
  Ghost::GhostType gn           = Ghost::None;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  t->needs(Task::ParentOldDW, d_ice_labels->delTLabel, getLevel(patches));

  //__________________________________
  // define parent data warehouse
  Task::WhichDW pNewDW = Task::NewDW;
  if (recursion) {
    pNewDW = Task::ParentNewDW;
  }

  t->needs(pNewDW, d_ice_labels->specificVolume_CCLabel, /*all_matls*/ gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->imp_delPLabel, press_matl, oims, gac, 1);
  t->needs(Task::OldDW, d_ice_labels->uvel_FCLabel, gn, 0);
  t->needs(Task::OldDW, d_ice_labels->vvel_FCLabel, gn, 0);
  t->needs(Task::OldDW, d_ice_labels->wvel_FCLabel, gn, 0);

  t->computes(d_ice_labels->uvel_FCLabel);
  t->computes(d_ice_labels->vvel_FCLabel);
  t->computes(d_ice_labels->wvel_FCLabel);
  t->computes(d_ice_labels->grad_dp_XFCLabel);
  t->computes(d_ice_labels->grad_dp_YFCLabel); // debugging variables
  t->computes(d_ice_labels->grad_dp_ZFCLabel);
  sched->addTask(t, patches, all_matls);

  //__________________________________
  //  add exchange contribution
  d_exchModel->sched_AddExch_VelFC(sched,
                                   patches,
                                   ice_matls,
                                   mpm_matls,
                                   all_matls,
                                   d_BC_globalVars.get(),
                                   recursion);
}
/*___________________________________________________________________
 Function~  ICE::scheduleComputeDel_P--
 Note:      This task is scheduled outside the iteration loop
_____________________________________________________________________*/
void
ICE::scheduleComputeDel_P(SchedulerP& sched,
                          const LevelP& level,
                          const PatchSet* patches,
                          const MaterialSubset* one_matl,
                          const MaterialSubset* press_matl,
                          const MaterialSet* all_matls)
{
  Task* t;
  Ghost::GhostType gn           = Ghost::None;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  //__________________________________
  // update the pressure

  printSchedule(level, cout_doing, "ICE::scheduleComputeDel_P");

  t = scinew Task("ICE::computeDel_P", this, &ICE::computeDel_P);

  t->needs(
    Task::NewDW, d_ice_labels->sum_imp_delPLabel, press_matl, oims, gn);
  t->needs(Task::NewDW, d_ice_labels->sumKappaLabel, one_matl, oims, gn);
  t->needs(Task::NewDW, d_ice_labels->term2Label, one_matl, oims, gn);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, gn);

  t->computes(d_ice_labels->delP_DilatateLabel, press_matl, oims);
  t->computes(d_ice_labels->delP_MassXLabel, press_matl, oims);
  t->computes(d_ice_labels->sum_rho_CCLabel, one_matl, oims);
  // t->computes(d_ice_labels->initialGuessLabel,  one_matl,  oims);

  sched->addTask(t, patches, all_matls);
}
/*___________________________________________________________________
 Function~  ICE::scheduleImplicitPressureSolve--
_____________________________________________________________________*/
void
ICE::scheduleImplicitPressureSolve(SchedulerP& sched,
                                   const LevelP& level,
                                   const PatchSet*,
                                   const MaterialSubset* one_matl,
                                   const MaterialSubset* press_matl,
                                   const MaterialSubset* ice_matls,
                                   const MaterialSubset* mpm_matls,
                                   const MaterialSet* all_matls)
{
  printSchedule(level, cout_doing, "ICE::scheduleImplicitPressureSolve");

  // if we're here, we're compiling the outer taskgraph.  Then we should compile
  // the inner one too.
  d_recompileSubsched = true;
  Task* t             = scinew Task("ICE::implicitPressureSolve",
                        this,
                        &ICE::implicitPressureSolve,
                        level,
                        ice_matls,
                        mpm_matls);

  t->hasSubScheduler();
  Ghost::GhostType gac          = Ghost::AroundCells;
  Ghost::GhostType gaf          = Ghost::AroundFaces;
  Ghost::GhostType gn           = Ghost::None;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  //__________________________________
  // The computes and requires when you're looking up
  // from implicitPressure solve
  // OldDW = ParentOldDW
  // NewDW = ParentNewDW
#ifdef HAVE_HYPRE
  if (d_solver->getName() == "hypre") {
    t->needs(Task::OldDW, hypre_solver_label);
    t->computes(hypre_solver_label);
    //(string var, bool treatAsOld, bool copyData, bool noScrub, bool
    // notCopyData, bool noCheckpoint)
    sched->overrideVariableBehavior(
      hypre_solver_label->getName(), false, true, false, false, true);
  }
#endif

  //__________________________________
  // common Variables
  t->needs(Task::OldDW, d_ice_labels->timeStepLabel);
  t->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, gac, 2);
  t->needs(Task::NewDW, d_ice_labels->specificVolume_CCLabel, gac, 1);
  t->needs(Task::NewDW, d_ice_labels->rhsLabel, one_matl, oims, gn, 0);
  // t->needs( Task::OldDW, d_ice_labels->initialGuessLabel,  one_matl,
  // oims,gac,1);

  //__________________________________
  // SetupRHS
  if (d_models.size() > 0) {
    t->needs(Task::NewDW, d_ice_labels->modelMass_srcLabel, gn, 0);
  }
  t->needs(Task::NewDW, d_ice_labels->speedSound_CCLabel, gn, 0);
  t->needs(Task::NewDW, d_ice_labels->max_RHSLabel);

  //__________________________________
  // setup Matrix
  t->needs(Task::NewDW, d_ice_labels->sp_volX_FCLabel, gaf, 1);
  t->needs(Task::NewDW, d_ice_labels->sp_volY_FCLabel, gaf, 1);
  t->needs(Task::NewDW, d_ice_labels->sp_volZ_FCLabel, gaf, 1);
  t->needs(Task::NewDW, d_ice_labels->vol_fracX_FCLabel, gaf, 1);
  t->needs(Task::NewDW, d_ice_labels->vol_fracY_FCLabel, gaf, 1);
  t->needs(Task::NewDW, d_ice_labels->vol_fracZ_FCLabel, gaf, 1);
  t->needs(
    Task::NewDW, d_ice_labels->sumKappaLabel, press_matl, oims, gn, 0);

  //__________________________________
  // Update Pressure
  t->needs(
    Task::NewDW, d_ice_labels->press_equil_CCLabel, press_matl, oims, gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->sum_imp_delPLabel, press_matl, oims, gac, 1);

  computesRequires_CustomBCs(t,
                             "implicitPressureSolve",
                             d_ice_labels.get(),
                             ice_matls,
                             d_BC_globalVars.get(),
                             true);

  //__________________________________
  //  For Vel_FC exchange
  d_exchModel->addExchangeModelRequires(t, one_matl, ice_matls, mpm_matls);

  //__________________________________
  // ImplicitVel_FC
  t->needs(Task::NewDW, d_ice_labels->uvel_FCLabel, gn, 0);
  t->needs(Task::NewDW, d_ice_labels->vvel_FCLabel, gn, 0);
  t->needs(Task::NewDW, d_ice_labels->wvel_FCLabel, gn, 0);

  //__________________________________
  //  what's produced from this task
  t->computes(d_ice_labels->press_CCLabel, press_matl, oims);
  t->computes(d_ice_labels->matrixLabel, one_matl, oims);
  t->computes(d_ice_labels->grad_dp_XFCLabel, press_matl, oims);
  t->computes(d_ice_labels->grad_dp_YFCLabel, press_matl, oims);
  t->computes(d_ice_labels->grad_dp_ZFCLabel, press_matl, oims);
  t->modifies(d_ice_labels->sum_imp_delPLabel, press_matl, oims);
  t->modifies(d_ice_labels->term2Label, one_matl, oims);
  t->modifies(d_ice_labels->rhsLabel, one_matl, oims);

  t->modifies(d_ice_labels->uvel_FCMELabel);
  t->modifies(d_ice_labels->vvel_FCMELabel);
  t->modifies(d_ice_labels->wvel_FCMELabel);

  t->modifies(d_ice_labels->vol_fracX_FCLabel);
  t->modifies(d_ice_labels->vol_fracY_FCLabel);
  t->modifies(d_ice_labels->vol_fracZ_FCLabel);

  t->computes(VarLabel::find(abortTimestep_name));
  t->computes(VarLabel::find(recomputeTimestep_name));

  const PatchSet* perproc_patches =
    d_loadBalancer->getPerProcessorPatchSet(level);

  sched->addTask(t, perproc_patches, all_matls);
}

/*___________________________________________________________________
 Function~  ICE::setupMatrix--
_____________________________________________________________________*/
void
ICE::setupMatrix(const ProcessorGroup*,
                 const PatchSubset* patches,
                 const MaterialSubset*,
                 DataWarehouse* old_dw,
                 DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ICE::setupMatrix");

    DataWarehouse* whichDW = old_dw;

    DataWarehouse* parent_old_dw =
      new_dw->getOtherDataWarehouse(Task::ParentOldDW);
    DataWarehouse* parent_new_dw =
      new_dw->getOtherDataWarehouse(Task::ParentNewDW);

    delt_vartype delT;
    parent_old_dw->get(delT, d_ice_labels->delTLabel, level);
    Vector dx    = patch->dCell();
    int numMatls = d_materialManager->getNumMaterials();
    CCVariable<Stencil7> A;
    CCVariable<double> imp_delP;
    constCCVariable<double> sumKappa;

    Ghost::GhostType gn  = Ghost::None;
    Ghost::GhostType gaf = Ghost::AroundFaces;
    new_dw->allocateAndPut(A, d_ice_labels->matrixLabel, 0, patch, gn, 0);
    new_dw->allocateAndPut(
      imp_delP, d_ice_labels->imp_delPLabel, 0, patch, gn, 0);
    parent_new_dw->get(sumKappa, d_ice_labels->sumKappaLabel, 0, patch, gn, 0);
    imp_delP.initialize(0.0); // you need to intialize the extra cells

    IntVector right, left, top, bottom, front, back;
    IntVector R_CC, L_CC, T_CC, B_CC, F_CC, BK_CC;

    //__________________________________
    //  Initialize A
    for (CellIterator iter(patch->getExtraCellIterator()); !iter.done();
         iter++) {
      IntVector c     = *iter;
      Stencil7& A_tmp = A[c];
      A_tmp.p         = 0.0;
      A_tmp.n         = 0.0;
      A_tmp.s         = 0.0;
      A_tmp.e         = 0.0;
      A_tmp.w         = 0.0;
      A_tmp.t         = 0.0;
      A_tmp.b         = 0.0;
    }

    for (int m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      constSFCXVariable<double> sp_volX_FC, vol_fracX_FC;
      constSFCYVariable<double> sp_volY_FC, vol_fracY_FC;
      constSFCZVariable<double> sp_volZ_FC, vol_fracZ_FC;

      whichDW->get(
        sp_volX_FC, d_ice_labels->sp_volX_FCLabel, indx, patch, gaf, 1);
      whichDW->get(
        sp_volY_FC, d_ice_labels->sp_volY_FCLabel, indx, patch, gaf, 1);
      whichDW->get(
        sp_volZ_FC, d_ice_labels->sp_volZ_FCLabel, indx, patch, gaf, 1);

      whichDW->get(
        vol_fracX_FC, d_ice_labels->vol_fracX_FCLabel, indx, patch, gaf, 1);
      whichDW->get(
        vol_fracY_FC, d_ice_labels->vol_fracY_FCLabel, indx, patch, gaf, 1);
      whichDW->get(
        vol_fracZ_FC, d_ice_labels->vol_fracZ_FCLabel, indx, patch, gaf, 1);

      //__________________________________
      // Sum (<upwinded volfrac> * sp_vol on faces)
      // +x -x +y -y +z -z
      //  e, w, n, s, t, b

      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c     = *iter;
        Stencil7& A_tmp = A[c];
        right           = c + IntVector(1, 0, 0);
        left            = c;
        top             = c + IntVector(0, 1, 0);
        bottom          = c;
        front           = c + IntVector(0, 0, 1);
        back            = c;

        //  use the upwinded vol_frac
        A_tmp.e += vol_fracX_FC[right] * sp_volX_FC[right];
        A_tmp.w += vol_fracX_FC[left] * sp_volX_FC[left];
        A_tmp.n += vol_fracY_FC[top] * sp_volY_FC[top];
        A_tmp.s += vol_fracY_FC[bottom] * sp_volY_FC[bottom];
        A_tmp.t += vol_fracZ_FC[front] * sp_volZ_FC[front];
        A_tmp.b += vol_fracZ_FC[back] * sp_volZ_FC[back];
      }
    } // matl loop

    //__________________________________
    //  Multiple stencil by delT^2 * area/dx
    double delT_2 = delT * delT;

    double vol     = dx.x() * dx.y() * dx.z();
    double tmp_e_w = dx.y() * dx.z() * delT_2 / dx.x();
    double tmp_n_s = dx.x() * dx.z() * delT_2 / dx.y();
    double tmp_t_b = dx.x() * dx.y() * delT_2 / dx.z();

    for (CellIterator iter(patch->getCellIterator()); !iter.done(); iter++) {
      IntVector c     = *iter;
      Stencil7& A_tmp = A[c];
      A_tmp.e *= -tmp_e_w;
      A_tmp.w *= -tmp_e_w;

      A_tmp.n *= -tmp_n_s;
      A_tmp.s *= -tmp_n_s;

      A_tmp.t *= -tmp_t_b;
      A_tmp.b *= -tmp_t_b;

      A_tmp.p = vol * sumKappa[c] -
                (A_tmp.n + A_tmp.s + A_tmp.e + A_tmp.w + A_tmp.t + A_tmp.b);
    }
    //__________________________________
    //  Boundary conditons on A.e, A.w, A.n, A.s, A.t, A.b
    ImplicitMatrixBC(A, patch);
  }
}

/*___________________________________________________________________
 Function~  ICE::setupRHS--
_____________________________________________________________________*/
void
ICE::setupRHS(const ProcessorGroup*,
              const PatchSubset* patches,
              const MaterialSubset*,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw,
              bool insideOuterIterLoop,
              string computes_or_modifies)
{
  const Level* level = getLevel(patches);
  Vector dx          = level->dCell();
  double vol         = dx.x() * dx.y() * dx.z();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ICE::setupRHS");
    // define parent_new/old_dw

    DataWarehouse* pNewDW;
    DataWarehouse* pOldDW;
    if (insideOuterIterLoop) {
      pNewDW = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      pOldDW = new_dw->getOtherDataWarehouse(Task::ParentOldDW);
    } else {
      pNewDW = new_dw;
      pOldDW = old_dw;
    }

    int numMatls = d_materialManager->getNumMaterials();
    delt_vartype delT;
    pOldDW->get(delT, d_ice_labels->delTLabel, level);

    std::unique_ptr<Advector> advector =
      d_advector->clone(new_dw, patch, isRegridTimestep());

    CCVariable<double> q_advected, rhs;
    CCVariable<double> sumAdvection, massExchTerm;
    constCCVariable<double> press_CC, oldPressure, speedSound, sumKappa;
    constCCVariable<double> sum_imp_delP;
    const IntVector gc(1, 1, 1);
    Ghost::GhostType gn  = Ghost::None;
    Ghost::GhostType gac = Ghost::AroundCells;

    new_dw->get(sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch, gn, 0);
    pNewDW->get(sumKappa, d_ice_labels->sumKappaLabel, 0, patch, gn, 0);

    new_dw->allocateAndPut(massExchTerm, d_ice_labels->term2Label, 0, patch);
    new_dw->allocateTemporary(q_advected, patch);
    new_dw->allocateTemporary(sumAdvection, patch);

    if (computes_or_modifies == "computes") {
      new_dw->allocateAndPut(rhs, d_ice_labels->rhsLabel, 0, patch);
    }
    if (computes_or_modifies == "modifies") {
      new_dw->getModifiable(rhs, d_ice_labels->rhsLabel, 0, patch);
    }

    rhs.initialize(0.0);
    sumAdvection.initialize(0.0);
    massExchTerm.initialize(0.0);

    for (int m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      constSFCXVariable<double> uvel_FC;
      constSFCYVariable<double> vvel_FC;
      constSFCZVariable<double> wvel_FC;
      SFCXVariable<double> vol_fracX_FC;
      SFCYVariable<double> vol_fracY_FC;
      SFCZVariable<double> vol_fracZ_FC;
      constCCVariable<double> vol_frac, burnedMass, sp_vol_CC, speedSound;

      new_dw->allocateAndPut(
        vol_fracX_FC, d_ice_labels->vol_fracX_FCLabel, indx, patch);
      new_dw->allocateAndPut(
        vol_fracY_FC, d_ice_labels->vol_fracY_FCLabel, indx, patch);
      new_dw->allocateAndPut(
        vol_fracZ_FC, d_ice_labels->vol_fracZ_FCLabel, indx, patch);

      // lowIndex is the same for all vel_FC
      IntVector lowIndex(patch->getExtraSFCXLowIndex());
      double nan = getNan();
      vol_fracX_FC.initialize(nan, lowIndex, patch->getExtraSFCXHighIndex());
      vol_fracY_FC.initialize(nan, lowIndex, patch->getExtraSFCYHighIndex());
      vol_fracZ_FC.initialize(nan, lowIndex, patch->getExtraSFCZHighIndex());

      new_dw->get(uvel_FC, d_ice_labels->uvel_FCMELabel, indx, patch, gac, 2);
      new_dw->get(vvel_FC, d_ice_labels->vvel_FCMELabel, indx, patch, gac, 2);
      new_dw->get(wvel_FC, d_ice_labels->wvel_FCMELabel, indx, patch, gac, 2);
      pNewDW->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, gac, 2);
      pNewDW->get(sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, gn, 0);
      pNewDW->get(
        speedSound, d_ice_labels->speedSound_CCLabel, indx, patch, gn, 0);

      //__________________________________
      // Advection preprocessing
      bool bulletProof_test = false;

      //__________________________________
      // common variables that get passed into the advection operators
      std::unique_ptr<advectVarBasket> varBasket         = std::make_unique<advectVarBasket>();
      varBasket->new_dw                  = new_dw;
      varBasket->old_dw                  = old_dw;
      varBasket->indx                    = indx;
      varBasket->patch                   = patch;
      varBasket->level                   = level;
      varBasket->lb                      = d_ice_labels.get();
      varBasket->doRefluxing             = isAMR(); // always reflux with amr
      varBasket->is_Q_massSpecific       = false;
      varBasket->useCompatibleFluxes     = d_useCompatibleFluxes;
      varBasket->AMR_subCycleProgressVar = 0; // for lockstep it's always 0

      advector->inFluxOutFluxVolume(uvel_FC,
                                    vvel_FC,
                                    wvel_FC,
                                    delT,
                                    patch,
                                    indx,
                                    bulletProof_test,
                                    pNewDW,
                                    varBasket.get());

      advector->advectQ(vol_frac,
                        patch,
                        q_advected,
                        varBasket.get(),
                        vol_fracX_FC,
                        vol_fracY_FC,
                        vol_fracZ_FC,
                        new_dw);

      varBasket.reset(nullptr);

      //__________________________________
      //  sum Advecton (<vol_frac> vel_FC )
      //  you need to multiply by vol
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        sumAdvection[c] += q_advected[c] * vol;
      }
      //__________________________________
      //  sum mass Exchange term
      if (d_models.size() > 0) {
        pNewDW->get(
          burnedMass, d_ice_labels->modelMass_srcLabel, indx, patch, gn, 0);
        for (CellIterator iter = patch->getCellIterator(); !iter.done();
             iter++) {
          IntVector c = *iter;
          massExchTerm[c] += burnedMass[c] * sp_vol_CC[c];
        }
      }
    } // matl loop

    //__________________________________
    // form the pressure difference term
    CCVariable<double> term1;
    new_dw->allocateTemporary(term1, patch);

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      term1[c]    = vol * sumKappa[c] * sum_imp_delP[c];
    }

    //__________________________________
    //  Form RHS
    // note:  massExchangeTerm has delT incorporated inside of it
    // We need to include the cell volume in rhs for AMR to be properly scaled
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      rhs[c]      = -term1[c] + massExchTerm[c] + sumAdvection[c];
    }

    // Renormalize massExchangeTerm to be consistent with the rest of ICE
    if (d_models.size() > 0) {
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        massExchTerm[c] /= vol;
      }
    }
    //__________________________________
    // set rhs = 0 under all fine patches
    // For a composite grid we ignore what's happening
    // on the coarse grid
    Level::selectType finePatches;
    if (level->hasFinerLevel()) {
      patch->getFineLevelPatches(finePatches);
    }

    for (size_t i = 0; i < finePatches.size(); i++) {
      const Patch* finePatch = finePatches[i];

      IntVector l, h, fl, fh;
      getFineLevelRange(patch, finePatch, l, h, fl, fh);

      if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
        continue;
      }
      rhs.initialize(0.0, l, h);
    }
  } // patches loop
  //  std::cout << " Level " << level->getIndex() << " rhs "
  //       << rhs_max << " rhs * vol " << rhs_max * vol <<  std::endl;
}

/*___________________________________________________________________
 Function~  ICE::compute_maxRHS--
_____________________________________________________________________*/
void
ICE::compute_maxRHS(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset*,
                    DataWarehouse*,
                    DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  double rhs_max     = 0.0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ICE::compute_maxRHS");

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    constCCVariable<double> rhs;
    new_dw->get(rhs, d_ice_labels->rhsLabel, 0, patch, Ghost::None, 0);

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      rhs_max     = Max(rhs_max, Abs(rhs[c] / vol));
    }
    new_dw->put(max_vartype(rhs_max), d_ice_labels->max_RHSLabel);

    //__________________________________
    // debugging output
    if (cout_dbg.active()) {
      rhs_max = 0.0;
      IntVector maxCell(0, 0, 0);
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        if (Abs(rhs[c] / vol) > rhs_max) {
          maxCell = c;
        }
        rhs_max = Max(rhs_max, Abs(rhs[c] / vol));
      }
      std::cout << " maxRHS: " << maxCell << " " << rhs_max << " \t L-"
                << level->getIndex() << std::endl;
    }
  } // patch loop
}

/*___________________________________________________________________
 Function~  ICE::updatePressure--
  add sum_imp_delP to press_equil
  - set boundary condtions on
_____________________________________________________________________*/
void
ICE::updatePressure(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset*,
                    DataWarehouse* old_dw,
                    DataWarehouse* new_dw)
{
  // define parent_dw
  DataWarehouse* parent_new_dw =
    new_dw->getOtherDataWarehouse(Task::ParentNewDW);
  DataWarehouse* parent_old_dw =
    new_dw->getOtherDataWarehouse(Task::ParentOldDW);

  timeStep_vartype timeStep;
  parent_old_dw->get(timeStep, d_ice_labels->timeStepLabel);

  bool isNotInitialTimestep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ICE::updatePressure");

    int numMatls        = d_materialManager->getNumMaterials();
    Ghost::GhostType gn = Ghost::None;

    CCVariable<double> press_CC;
    CCVariable<double> imp_delP;
    CCVariable<double> sum_imp_delP;
    constCCVariable<double> press_equil;
    constCCVariable<double> sum_imp_delP_old;
    std::vector<CCVariable<double>> placeHolder(0);
    std::vector<constCCVariable<double>> sp_vol_CC(numMatls);

    old_dw->get(
      sum_imp_delP_old, d_ice_labels->sum_imp_delPLabel, 0, patch, gn, 0);
    parent_new_dw->get(
      press_equil, d_ice_labels->press_equil_CCLabel, 0, patch, gn, 0);
    new_dw->getModifiable(imp_delP, d_ice_labels->imp_delPLabel, 0, patch);
    new_dw->allocateAndPut(press_CC, d_ice_labels->press_CCLabel, 0, patch);
    new_dw->allocateAndPut(
      sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch);
    press_CC.initialize(d_EVIL_NUM);

    for (int m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      parent_new_dw->get(
        sp_vol_CC[m], d_ice_labels->specificVolume_CCLabel, indx, patch, gn, 0);
    }
    // set boundary conditions on imp_delP
    set_imp_DelP_BC(imp_delP, patch, d_ice_labels->imp_delPLabel, new_dw);
    //__________________________________
    //  add delP to press_equil
    //  AMR:  hit the extra cells, you need to update the pressure in these
    //  cells
    for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
         iter++) {
      IntVector c     = *iter;
      sum_imp_delP[c] = sum_imp_delP_old[c] + imp_delP[c];
      press_CC[c]     = press_equil[c] + sum_imp_delP[c];
      press_CC[c]     = std::max(1.0e-12, press_CC[c]); // C L A M P
    }
    //__________________________________
    //  set boundary conditions
    std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
      std::make_unique<CustomBCDriver::customBC_localVars>();

    preprocess_CustomBCs("imp_update_press_CC",
                         parent_old_dw,
                         parent_new_dw,
                         d_ice_labels.get(),
                         patch,
                         999,
                         d_BC_globalVars.get(),
                         BC_localVars.get());

    setBC(press_CC,
          placeHolder,
          sp_vol_CC,
          d_surroundingMatl_indx,
          "sp_vol",
          "Pressure",
          patch,
          d_materialManager,
          0,
          new_dw,
          d_BC_globalVars.get(),
          BC_localVars.get(),
          isNotInitialTimestep);

    delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());

    //____ B U L L E T   P R O O F I N G----
    // ignore BP if a recompute time step has already been requested
    IntVector neg_cell;
    bool rts = new_dw->recomputeTimestep();

    if (!areAllValuesPositive(press_CC, neg_cell) && !rts) {
      std::ostringstream warn;
      warn << "ERROR ICE::updatePressure cell " << neg_cell
           << " negative pressure\n ";
      throw InvalidValue(warn.str(), __FILE__, __LINE__);
    }
  } // patch loop
}

/*___________________________________________________________________
 Function~  ICE::computeDel_P--
_____________________________________________________________________*/
void
ICE::computeDel_P(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset*,
                  DataWarehouse*,
                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ICE::computeDel_P");

    int numMatls = d_materialManager->getNumMaterials();

    CCVariable<double> delP_Dilatate;
    CCVariable<double> delP_MassX;
    CCVariable<double> sum_rho_CC;
    // CCVariable<double> initialGuess;
    constCCVariable<double> rho_CC;
    constCCVariable<double> sumKappa;
    constCCVariable<double> massExchTerm;
    constCCVariable<double> press_equil, press_CC;
    constCCVariable<double> sum_imp_delP;

    Ghost::GhostType gn = Ghost::None;
    new_dw->get(sumKappa, d_ice_labels->sumKappaLabel, 0, patch, gn, 0);
    new_dw->get(massExchTerm, d_ice_labels->term2Label, 0, patch, gn, 0);
    new_dw->get(sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch, gn, 0);

    new_dw->allocateAndPut(
      delP_Dilatate, d_ice_labels->delP_DilatateLabel, 0, patch);
    new_dw->allocateAndPut(delP_MassX, d_ice_labels->delP_MassXLabel, 0, patch);
    new_dw->allocateAndPut(sum_rho_CC, d_ice_labels->sum_rho_CCLabel, 0, patch);
    //  new_dw->allocateAndPut(initialGuess, d_ice_labels->initialGuessLabel, 0,
    //  patch);

    sum_rho_CC.initialize(0.0);
    delP_Dilatate.initialize(0.0);
    delP_MassX.initialize(0.0);

    for (int m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, gn, 0);
      //__________________________________
      //  compute sum_rho_CC used by press_FC
      for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        sum_rho_CC[c] += rho_CC[c];
      }
    }
    //__________________________________
    // backout delP_Dilatate and delP_MassX
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c      = *iter;
      delP_MassX[c]    = massExchTerm[c] / sumKappa[c];
      delP_Dilatate[c] = sum_imp_delP[c] - delP_MassX[c];
      // initialGuess[c]  = sum_imp_delP[c];
    }

  } // patch loop
}

/*___________________________________________________________________
 Function~  ICE::implicitPressureSolve--
_____________________________________________________________________*/
void
ICE::implicitPressureSolve([[maybe_unused]] const ProcessorGroup* pg,
                           const PatchSubset* patch_sub,
                           const MaterialSubset*,
                           DataWarehouse* ParentOldDW,
                           DataWarehouse* ParentNewDW,
                           LevelP level,
                           const MaterialSubset* ice_matls,
                           const MaterialSubset* mpm_matls)
{
  int proc = d_myworld->myRank();
  printTask(patch_sub, cout_doing, "Doing implicitPressureSolve");

  //__________________________________
  // define Matl sets and subsets
  const MaterialSet* all_matls        = d_materialManager->allMaterials();
  const MaterialSubset* all_matls_sub = all_matls->getUnion();
  const MaterialSubset* one_matl      = d_press_matl;

  //__________________________________
  //  turn off parentDW scrubbing
  DataWarehouse::ScrubMode ParentOldDW_scrubmode =
    ParentOldDW->setScrubbing(DataWarehouse::ScrubNone);
  DataWarehouse::ScrubMode ParentNewDW_scrubmode =
    ParentNewDW->setScrubbing(DataWarehouse::ScrubNone);

  GridP grid = level->getGrid();
  d_subsched->setParentDWs(ParentOldDW, ParentNewDW);
  d_subsched->advanceDataWarehouse(grid);
  d_subsched->setInitTimestep(true);

  bool recursion            = true;
  bool modifies_X           = true;
  const PatchSet* patch_set = level->eachPatch();

  const VarLabel* initialGuessLabel = nullptr;
  // const VarLabel* initialGuessLabel = d_ice_labels->initialGuessLabel;
  Task::WhichDW initialGuess_DW = Task::OldDW;

  DataWarehouse* subOldDW = d_subsched->get_dw(2);
  DataWarehouse* subNewDW = d_subsched->get_dw(3);

  //__________________________________
  //  Move data from parentOldDW to subSchedNewDW.
  max_vartype max_RHS_old;
  ParentNewDW->get(max_RHS_old, d_ice_labels->max_RHSLabel);
  subNewDW->put(max_RHS_old, d_ice_labels->max_RHSLabel);

#ifdef HAVE_HYPRE
  if (d_solver->getName() == "hypre") {
    if (ParentOldDW->exists(hypre_solver_label)) {
      SoleVariable<hypre_solver_structP> hypre_solverP;
      ParentOldDW->get(hypre_solverP, hypre_solver_label);
      subNewDW->put(hypre_solverP, hypre_solver_label);
    }
    d_solver->getParameters()->setWhichOldDW(Task::ParentOldDW);
  }
#endif

  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->sum_imp_delPLabel, patch_sub, d_press_matl);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->rhsLabel, patch_sub, one_matl);

  const MaterialSubset* all_matls_s = all_matls->getUnion();
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->uvel_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->vvel_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->wvel_FCLabel, patch_sub, all_matls_s);

  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->sp_volX_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->sp_volY_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->sp_volZ_FCLabel, patch_sub, all_matls_s);

  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->vol_fracX_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->vol_fracY_FCLabel, patch_sub, all_matls_s);
  subNewDW->transferFrom(
    ParentNewDW, d_ice_labels->vol_fracZ_FCLabel, patch_sub, all_matls_s);

  // subNewDW->transferFrom( ParentOldDW,d_ice_labels->initialGuessLabel,
  // patch_sub, one_matl );
  //__________________________________
  //   Iteration Loop
  max_vartype max_RHS           = 1 / d_SMALL_NUM;
  double smallest_max_RHS_sofar = max_RHS;
  int counter                   = 0;
  bool recomputeTimestep        = false;
  Vector dx                     = level->dCell();
  double vol                    = dx.x() * dx.y() * dx.z();

  d_solver->getParameters()->setResidualNormalizationFactor(vol);

  // debugging set output file name for hypre objects
  timeStep_vartype timeStepVar;
  ParentOldDW->get(timeStepVar, d_ice_labels->timeStepLabel);
  double timeStep = timeStepVar;

  std::ostringstream fname;
  fname << "." << proc << "." << timeStep << "." << counter;
  d_solver->getParameters()->setOutputFileName(fname.str());

  d_subsched->setInitTimestep(false);

  while (counter < d_max_iter_implicit && max_RHS > d_outer_iter_tolerance &&
         !recomputeTimestep) {
    //__________________________________
    // recompile the subscheduler

    if (counter < 2 && d_recompileSubsched) {

      d_subsched->initialize(3, 1);
      //__________________________________
      // schedule the tasks

      scheduleSetupMatrix(d_subsched, level, patch_set, one_matl, all_matls);

      d_solver->scheduleSolve(level,
                              d_subsched,
                              d_press_matlSet,
                              d_ice_labels->matrixLabel,
                              Task::NewDW,
                              d_ice_labels->imp_delPLabel,
                              modifies_X,
                              d_ice_labels->rhsLabel,
                              Task::OldDW,
                              initialGuessLabel,
                              initialGuess_DW,
                              true);

      scheduleUpdatePressure(d_subsched,
                             level,
                             patch_set,
                             ice_matls,
                             mpm_matls,
                             d_press_matl,
                             all_matls);

      scheduleRecomputeVel_FC(d_subsched,
                              patch_set,
                              ice_matls,
                              mpm_matls,
                              d_press_matl,
                              all_matls,
                              recursion);

      scheduleSetupRHS(
        d_subsched, patch_set, one_matl, all_matls, recursion, "computes");

      scheduleCompute_maxRHS(d_subsched, level, one_matl, all_matls);

      d_subsched->compile();
      // d_recompileSubsched = false;
    }
    //__________________________________
    //  - move subNewDW to subOldDW
    //  - scrub the subScheduler
    //  - execute the tasks
    d_subsched->advanceDataWarehouse(grid);
    subOldDW = d_subsched->get_dw(2);
    subNewDW = d_subsched->get_dw(3);
    subOldDW->setScrubbing(DataWarehouse::ScrubComplete);
    subNewDW->setScrubbing(DataWarehouse::ScrubNone);

    d_subsched->execute();

    counter++;
    //  initialGuessLabel = d_ice_labels->sum_imp_delPLabel;
    //  initialGuess_DW   = Task::OldDW;

    //__________________________________
    // diagnostics
    subNewDW->get(max_RHS, d_ice_labels->max_RHSLabel);
    subOldDW->get(max_RHS_old, d_ice_labels->max_RHSLabel);

    DOUT(proc == 0,
         "  Outer iteration " << counter << " max_rhs before solve "
                              << max_RHS_old << " after solve " << max_RHS);

    //__________________________________
    // recompute timestep
    //  too many outer iterations
    if (counter > d_iters_before_timestep_recompute) {
      recomputeTimestep = true;
      DOUT(proc == 0,
           "\nWARNING: The max iterations occurred before a time step "
           "recompute was reached");
    }
    //  The solver or advection has requested to recompute the time step
    if (subNewDW->recomputeTimestep()) {
      DOUT(proc == 0,
           "\n  WARNING:  impICE:implicitPressureSolve time step recompute.");
      recomputeTimestep = true;
    }

    //  solution is diverging
    if (max_RHS < smallest_max_RHS_sofar) {
      smallest_max_RHS_sofar = max_RHS;
    }

    if (((max_RHS - smallest_max_RHS_sofar) > 100.0 * smallest_max_RHS_sofar)) {
      DOUT(proc == 0,
           "\nWARNING: outer iteration is diverging now "
             << "recomputing the timestep"
             << " Max_RHS " << max_RHS << " smallest_max_RHS_sofar "
             << smallest_max_RHS_sofar);
      recomputeTimestep = true;
    }

    if (recomputeTimestep) {
      ParentNewDW->put(bool_or_vartype(true),
                       VarLabel::find(abortTimestep_name));
      ParentNewDW->put(bool_or_vartype(true),
                       VarLabel::find(recomputeTimestep_name));
    }
  } // outer iteration loop

  //__________________________________
  //  BULLET PROOFING
  if ((counter == d_max_iter_implicit) && (max_RHS > d_outer_iter_tolerance) &&
      counter > 1) {
    std::ostringstream s;
    s << "ERROR ICE::implicitPressureSolve, the maximum number of outer"
      << " iterations was reached. \n "
      << "Try either increasing the max_outer_iterations "
      << " or decrease outer_iteration_tolerance\n ";
    throw ConvergenceFailure(
      s.str(), counter, max_RHS, d_outer_iter_tolerance, __FILE__, __LINE__);
  }

  //__________________________________
  // Move products of iteration (only) from sub_new_dw -> parent_new_dw
  subNewDW     = d_subsched->get_dw(3);
  bool replace = true;

#ifdef HAVE_HYPRE
  if (d_solver->getName() == "hypre") {
    if (subNewDW->exists(hypre_solver_label)) {
      SoleVariable<hypre_solver_structP> hypre_solverP;
      subNewDW->get(hypre_solverP, hypre_solver_label);
      ParentNewDW->put(hypre_solverP, hypre_solver_label);
    }
  }
#endif

  ParentNewDW->transferFrom(subNewDW, // press
                            d_ice_labels->press_CCLabel,
                            patch_sub,
                            d_press_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // press
                            d_ice_labels->matrixLabel,
                            patch_sub,
                            one_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW,
                            d_ice_labels->sum_imp_delPLabel,
                            patch_sub,
                            d_press_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // term2
                            d_ice_labels->term2Label,
                            patch_sub,
                            one_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // rhs
                            d_ice_labels->rhsLabel,
                            patch_sub,
                            one_matl,
                            replace);

  ParentNewDW->transferFrom(subNewDW, // uvel_FC
                            d_ice_labels->uvel_FCMELabel,
                            patch_sub,
                            all_matls_sub,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // vvel_FC
                            d_ice_labels->vvel_FCMELabel,
                            patch_sub,
                            all_matls_sub,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // wvel_FC
                            d_ice_labels->wvel_FCMELabel,
                            patch_sub,
                            all_matls_sub,
                            replace);

  ParentNewDW->transferFrom(subNewDW, // grad_impDelP_XFC
                            d_ice_labels->grad_dp_XFCLabel,
                            patch_sub,
                            d_press_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // grad_impDelP_YFC
                            d_ice_labels->grad_dp_YFCLabel,
                            patch_sub,
                            d_press_matl,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // grad_impDelP_ZFC
                            d_ice_labels->grad_dp_ZFCLabel,
                            patch_sub,
                            d_press_matl,
                            replace);

  ParentNewDW->transferFrom(subNewDW, // vol_fracX_FC
                            d_ice_labels->vol_fracX_FCLabel,
                            patch_sub,
                            all_matls_sub,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // vol_fracY_FC
                            d_ice_labels->vol_fracY_FCLabel,
                            patch_sub,
                            all_matls_sub,
                            replace);
  ParentNewDW->transferFrom(subNewDW, // vol_fracZ_FC
                            d_ice_labels->vol_fracZ_FCLabel,
                            patch_sub,
                            all_matls_sub,
                            replace);

  //__________________________________
  //  Turn scrubbing back on
  ParentOldDW->setScrubbing(ParentOldDW_scrubmode);
  ParentNewDW->setScrubbing(ParentNewDW_scrubmode);
}
