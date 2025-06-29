/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd, NZ
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

#include <CCA/Components/ICE/AMRICE.h>

#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/Models/FluidsBased/FluidsBasedModel.h>
#include <CCA/Components/Models/HEChem/HEChemModel.h>

#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SolverInterface.h>

#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AMR_CoarsenRefine.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/PerPatchVars.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Math/MiscMath.h>
#include <Core/Math/Rand48.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Util/DebugStream.h>

#include <cstdio>
#include <iomanip>

using namespace Uintah;
using namespace std;

static DebugStream cout_doing("AMRICE_DOING_COUT", false);

// setenv SCI_DEBUG AMR:+ you can see the new grid it creates

AMRICE::AMRICE(const ProcessorGroup* myworld,
               const MaterialManagerP materialManager)
  : ICE(myworld, materialManager)
{
}

AMRICE::~AMRICE() {}

//___________________________________________________________________
void
AMRICE::problemSetup(const ProblemSpecP& params,
                     const ProblemSpecP& restart_prob_spec,
                     GridP& grid,
                     [[maybe_unused]] const std::string& input_ups_dir)
{
  cout_doing << d_myworld->myRank() << " Doing problemSetup  \t\t\t AMRICE"
             << '\n';

  ICE::problemSetup(params, restart_prob_spec, grid);
  ProblemSpecP ice_ps;
  ProblemSpecP amr_ps = params->findBlock("AMR");

  ProblemSpecP reg_ps = amr_ps->findBlock("Regridder");
  if (reg_ps) {

    std::string regridder;
    reg_ps->getAttribute("type", regridder);

    if (regridder != "Tiled" && regridder != "SingleLevel") {
      std::ostringstream msg;
      msg << "\n    ERROR:AMRICE With the (" << regridder
          << ") regridder the refine() task will overwrite \n";
      msg << "    all data in any newly created patches on the fine level "
             "patches.  There could be valid data on these patches. \n";
      msg << "    To prevent this you must select the \"Tiled\" regridder\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  }

  if (amr_ps) {
    ice_ps = amr_ps->findBlock("ICE");
  }

  if (!ice_ps) {
    std::string warn;
    warn =
      "\n INPUT FILE ERROR:\n <ICE>  block not found inside of <AMR> block \n";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  ProblemSpecP refine_ps = ice_ps->findBlock("Refinement_Criteria_Thresholds");
  if (!refine_ps) {
    std::string warn;
    warn = "\n INPUT FILE ERROR:\n <Refinement_Criteria_Thresholds> "
           " block not found inside of <ICE> block \n";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }
  ice_ps->require("orderOfInterpolation", d_orderOfInterpolation);
  ice_ps->getWithDefault("do_Refluxing", d_doRefluxing, true);
  ice_ps->getWithDefault("orderOf_CFI_Interpolation",
                         d_orderOf_CFI_Interpolation,
                         d_orderOfInterpolation);

  //__________________________________
  // bulletproofing.  Refluxing and first order advection
  ProblemSpecP cfd_ps    = params->findBlock("CFD");
  ProblemSpecP iceps     = cfd_ps->findBlock("ICE");
  ProblemSpecP advect_ps = iceps->findBlock("advection");
  std::map<std::string, std::string> advect_options;
  advect_ps->getAttributes(advect_options);
  if (advect_options["type"] == "FirstOrder" && d_doRefluxing) {
    throw ProblemSetupException("\n ICE: You cannot use AMR refluxing and the "
                                "first order advection operator together."
                                "  The results are significantly worse.",
                                __FILE__,
                                __LINE__);
  }

  //__________________________________
  // Pull out the refinement threshold criteria
  for (ProblemSpecP var_ps = refine_ps->findBlock("Variable");
       var_ps != nullptr;
       var_ps = var_ps->findNextBlock("Variable")) {
    thresholdVar data;
    std::string name, value, matl;

    std::map<std::string, std::string> input;
    var_ps->getAttributes(input);
    name  = input["name"];
    value = input["value"];
    matl  = input["matl"];

    std::stringstream n_ss(name);
    std::stringstream v_ss(value);
    std::stringstream m_ss(matl);

    n_ss >> data.name;
    v_ss >> data.value;
    m_ss >> data.matl;

    if (!n_ss || !v_ss || (!m_ss && matl != "all")) {
      printf("WARNING: AMRICE.cc: stringstream failed...\n");
    }

    int numMatls = d_materialManager->getNumMaterials();

    //__________________________________
    //  bulletproofing
    VarLabel* label = VarLabel::find(name);

    if (label == nullptr) {
      throw ProblemSetupException("The threshold variable name(" + name +
                                    ") could not be found",
                                  __FILE__,
                                  __LINE__);
    }

    if (data.name != "rho_CC" && data.name != "temp_CC" &&
        data.name != "vol_frac_CC" && data.name != "vel_CC" &&
        data.name != "press_CC") {
      std::ostringstream warn;
      warn << "\n INPUT FILE ERROR:\n The threshold variable name (" << name
           << ") is not valid\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    if (data.value < 0) {
      std::ostringstream warn;
      warn << "\n INPUT FILE ERROR:\n The threshold value (" << value
           << ") cannot be negative\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    if ((data.matl < 0 || data.matl > numMatls) && matl != "all") {
      std::ostringstream warn;
      warn << "\n INPUT FILE ERROR:\n The threshold material (" << matl
           << ") is not valid\n"
           << " select any material < total number of materials or 'all'";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }
    if (data.name == "Pressure") { // ignore what the user input, it's always 0
      data.matl = 0;
    }

    //__________________________________
    // if using "all" matls
    if (matl == "all") {
      for (int m = 0; m < numMatls; m++) {
        data.matl = m;
        d_thresholdVars.push_back(data);
      }

    } else {
      d_thresholdVars.push_back(data);
    }
  }

  //__________________________________
  // manual manipulate the scheduling of copy data

  // overrideVariableBehavior(string var, bool treatAsOld, bool copyData, bool
  // noScrub, bool notCopyData, bool noCheckpoint)
  //                              (1)             (2)             (3) (4) (5)
  //                              (6)
  //   treatAsOld:     var - will be checkpointed, copied, and only scrubbed
  //   from an OldDW copyData:       copy variable between AMR levels noScrub:
  //   set variable not to scrub (normally when needed between a normal
  //   taskgraph
  //                   and the regridding phase)
  //   notCopyData:    ignore copying this variable between AMR levels
  //   noCheckpoint:   do not checkpoint this variable.

  // we need these for AMRICE::refine
  d_scheduler->overrideVariableBehavior(
    "specific_heat", true, true, false, false, false);
  d_scheduler->overrideVariableBehavior(
    "gamma", true, true, false, false, false);
  d_scheduler->overrideVariableBehavior(
    "vol_frac_CC", true, true, false, false, false);
  d_scheduler->overrideVariableBehavior(
    "sp_vol_CC", d_with_mpm, true, false, false, false);
  if (d_with_mpm) {
    d_scheduler->overrideVariableBehavior(
      "temp_CC", true, true, false, false, false);
  }

  // We need these variables from OldDW to use between tasks, but do not
  //  schedule datacopy

  d_scheduler->overrideVariableBehavior(
    "mass_X_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "mass_Y_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "mass_Z_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "mom_X_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "mom_Y_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "mom_Z_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "sp_vol_X_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "sp_vol_Y_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "sp_vol_Z_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "int_eng_X_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "int_eng_Y_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "int_eng_Z_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "vol_frac_X_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "vol_frac_Y_FC_flux", false, false, true, false, false);
  d_scheduler->overrideVariableBehavior(
    "vol_frac_Z_FC_flux", false, false, true, false, false);

  //__________________________________
  // Model with reflux variables.
  if (d_models.size()) {
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());

      if (fb_model) {
        for (auto& rvar : fb_model->d_refluxVars) {

          std::string varLabelX = rvar->var_X_FC_flux->getName();
          std::string varLabelY = rvar->var_Y_FC_flux->getName();
          std::string varLabelZ = rvar->var_Z_FC_flux->getName();

          d_scheduler->overrideVariableBehavior(
            varLabelX, false, false, true, false, false);
          d_scheduler->overrideVariableBehavior(
            varLabelY, false, false, true, false, false);
          d_scheduler->overrideVariableBehavior(
            varLabelZ, false, false, true, false, false);
        }
      }
    }
  }
}
//___________________________________________________________________
void
AMRICE::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  cout_doing << d_myworld->myRank() << " AMRICE::scheduleInitialize \t\tL-"
             << level->getIndex() << '\n';
  ICE::scheduleInitialize(level, sched);
}
//___________________________________________________________________
void
AMRICE::initialize(const ProcessorGroup*,
                   const PatchSubset*,
                   const MaterialSubset*,
                   DataWarehouse*,
                   DataWarehouse*)
{
}
/*___________________________________________________________________
 Function~  AMRICE::scheduleRefineInterface_Variable--
 Purpose:
_____________________________________________________________________*/
void
AMRICE::scheduleRefineInterface_Variable(const LevelP& fineLevel,
                                         SchedulerP& sched,
                                         const VarLabel* variable,
                                         Task::MaterialDomainSpec DS,
                                         const MaterialSet* matls,
                                         bool needCoarseOld,
                                         bool needCoarseNew)
{
  cout_doing << d_myworld->myRank() << " \t scheduleRefineInterface_Variable ("
             << variable->getName() << ") matls: \t" << *matls << endl;

  std::ostringstream taskName;
  taskName << "AMRICE::refineCoarseFineInterface(" << variable->getName()
           << ")";
  Task* t;

  void (AMRICE::*func)(const ProcessorGroup*,
                       const PatchSubset*,
                       const MaterialSubset*,
                       DataWarehouse*,
                       DataWarehouse*,
                       const VarLabel*);

  switch (variable->typeDescription()->getSubType()->getType()) {
    case TypeDescription::Type::double_type:
      func = &AMRICE::refineCoarseFineInterface<double>;
      t    = scinew Task(taskName.str().c_str(), this, func, variable);
      break;
    case TypeDescription::Type::Vector:
      func = &AMRICE::refineCoarseFineInterface<Vector>;
      t    = scinew Task(taskName.str().c_str(), this, func, variable);
      break;
    default:
      throw InternalError(
        "Unknown variable type for AMRICE::scheduleRefineInterface_Variable",
        __FILE__,
        __LINE__);
  }

  Task::SearchTG OldTG =
    Task::SearchTG::OldTG; // possibly search the OldTG for computes

  const MaterialSubset* matls_sub = matls->getUnion();

  if (needCoarseOld) {
    cout_dbg << " requires from CoarseOldDW ";
    t->needs(Task::CoarseOldDW,
                variable,
                0,
                Task::CoarseLevel,
                matls_sub,
                DS,
                d_gac,
                1);
  }
  if (needCoarseNew) {
    cout_dbg << " requires from CoarseNewDW ";
    t->needs(Task::CoarseNewDW,
                variable,
                0,
                Task::CoarseLevel,
                matls_sub,
                DS,
                d_gac,
                1,
                OldTG);
  }

  t->modifies(variable, matls_sub, DS, OldTG);

  sched->addTask(t, fineLevel->eachPatch(), matls);
}
/*___________________________________________________________________
 Function~  AMRICE::scheduleRefineInterface--
 Purpose:
_____________________________________________________________________*/
void
AMRICE::scheduleRefineInterface(const LevelP& fineLevel,
                                SchedulerP& sched,
                                bool needCoarseOld,
                                bool needCoarseNew)
{
  if (fineLevel->getIndex() > 0) {
    cout_doing << d_myworld->myRank()
               << " AMRICE::scheduleRefineInterface \t\t\tL-"
               << fineLevel->getIndex() << " coarseOld: " << needCoarseOld
               << " coarseNew: " << needCoarseNew << endl;

    Task::MaterialDomainSpec ND   = Task::NormalDomain;
    Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice
                                                       // matlSet.
    const MaterialSet* all_matls = d_materialManager->allMaterials();
    const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");

    scheduleRefineInterface_Variable(fineLevel,
                                     sched,
                                     d_ice_labels->press_CCLabel,
                                     oims,
                                     d_press_matlSet,
                                     needCoarseOld,
                                     needCoarseNew);
    scheduleRefineInterface_Variable(fineLevel,
                                     sched,
                                     d_ice_labels->rho_CCLabel,
                                     ND,
                                     ice_matls,
                                     needCoarseOld,
                                     needCoarseNew);
    scheduleRefineInterface_Variable(fineLevel,
                                     sched,
                                     d_ice_labels->specificVolume_CCLabel,
                                     ND,
                                     all_matls,
                                     needCoarseOld,
                                     needCoarseNew);
    scheduleRefineInterface_Variable(fineLevel,
                                     sched,
                                     d_ice_labels->temperature_CCLabel,
                                     ND,
                                     all_matls,
                                     needCoarseOld,
                                     needCoarseNew);
    scheduleRefineInterface_Variable(fineLevel,
                                     sched,
                                     d_ice_labels->velocity_CCLabel,
                                     ND,
                                     ice_matls,
                                     needCoarseOld,
                                     needCoarseNew);

    //__________________________________
    // Model with transported variables.
    if (d_models.size()) {
      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());

        if (fb_model && fb_model->d_transVars.size()) {
          std::vector<TransportedVariable*>::iterator t_iter;

          for (t_iter = fb_model->d_transVars.begin();
               t_iter != fb_model->d_transVars.end();
               t_iter++) {
            TransportedVariable* tvar = *t_iter;

            scheduleRefineInterface_Variable(fineLevel,
                                             sched,
                                             tvar->var,
                                             ND,
                                             tvar->matlSet,
                                             needCoarseOld,
                                             needCoarseNew);
          }
        }
      }
    } // transported Vars
  }   // finer level
}

/*______________________________________________________________________
 Function~  AMRICE::refineCoarseFineInterface
 Purpose~
______________________________________________________________________*/
template<typename T>
void
AMRICE::refineCoarseFineInterface(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* matls,
                                  [[maybe_unused]] DataWarehouse* fine_old_dw,
                                  DataWarehouse* fine_new_dw,
                                  const VarLabel* variable)
{
  double subCycleProgress = getSubCycleProgress(fine_new_dw);
  const Level* fineLevel  = getLevel(patches);
  if (fineLevel->getIndex() > 0) {
    cout_doing << d_myworld->myRank() << " Doing refineCoarseFineInterface("
               << variable->getName() << ")\t\t\t AMRICE L-"
               << fineLevel->getIndex() << " Patches: " << *patches
               << " progressVar " << subCycleProgress << endl;

    for (int p = 0; p < patches->size(); p++) {
      const Patch* finePatch = patches->get(p);

      for (int m = 0; m < matls->size(); m++) {
        int indx = matls->get(m);

        CCVariable<T> Q_CC;
        fine_new_dw->getModifiable(Q_CC, variable, indx, finePatch);

        refineCoarseFineBoundaries(
          finePatch, Q_CC, fine_new_dw, variable, indx, subCycleProgress);
      }
    }
  }
}

/*___________________________________________________________________
 Function~  AMRICE::refineCoarseFineBoundaries--    D O U B L E
_____________________________________________________________________*/
void
AMRICE::refineCoarseFineBoundaries(const Patch* finePatch,
                                   CCVariable<double>& val,
                                   DataWarehouse* fine_new_dw,
                                   const VarLabel* label,
                                   int matl,
                                   double subCycleProgress_var)
{
  const Level* fineLevel   = finePatch->getLevel();
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();

  cout_dbg << "\t refineCoarseFineBoundaries (" << label->getName() << ") \t"
           << " subCycleProgress_var " << subCycleProgress_var << " Level-"
           << fineLevel->getIndex() << '\n';
  DataWarehouse* coarse_old_dw = 0;
  DataWarehouse* coarse_new_dw = 0;

  if (subCycleProgress_var != 1.0) {
    coarse_old_dw = fine_new_dw->getOtherDataWarehouse(Task::CoarseOldDW);
  }
  if (subCycleProgress_var != 0.0) {
    coarse_new_dw = fine_new_dw->getOtherDataWarehouse(Task::CoarseNewDW);
  }

  refine_CF_interfaceOperator<double>(finePatch,
                                      fineLevel,
                                      coarseLevel,
                                      val,
                                      label,
                                      subCycleProgress_var,
                                      matl,
                                      fine_new_dw,
                                      coarse_old_dw,
                                      coarse_new_dw);
}
/*___________________________________________________________________
 Function~  AMRICE::refineCoarseFineBoundaries--    V E C T O R
_____________________________________________________________________*/
void
AMRICE::refineCoarseFineBoundaries(const Patch* finePatch,
                                   CCVariable<Vector>& val,
                                   DataWarehouse* fine_new_dw,
                                   const VarLabel* label,
                                   int matl,
                                   double subCycleProgress_var)
{
  const Level* fineLevel   = finePatch->getLevel();
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();

  cout_dbg << "\t refineCoarseFineBoundaries (" << label->getName() << ") \t"
           << " subCycleProgress_var " << subCycleProgress_var << " Level-"
           << fineLevel->getIndex() << '\n';
  DataWarehouse* coarse_old_dw = 0;
  DataWarehouse* coarse_new_dw = 0;

  if (subCycleProgress_var != 1.0) {
    coarse_old_dw = fine_new_dw->getOtherDataWarehouse(Task::CoarseOldDW);
  }
  if (subCycleProgress_var != 0.0) {
    coarse_new_dw = fine_new_dw->getOtherDataWarehouse(Task::CoarseNewDW);
  }

  refine_CF_interfaceOperator<Vector>(finePatch,
                                      fineLevel,
                                      coarseLevel,
                                      val,
                                      label,
                                      subCycleProgress_var,
                                      matl,
                                      fine_new_dw,
                                      coarse_old_dw,
                                      coarse_new_dw);
}

/*___________________________________________________________________
 Function~  AMRICE::scheduleSetBC_FineLevel--
 Purpose:
_____________________________________________________________________*/
void
AMRICE::scheduleSetBC_FineLevel(const PatchSet* patches, SchedulerP& sched)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx             = fineLevel->getIndex();

  if (L_indx > 0) {
    cout_doing << d_myworld->myRank()
               << " AMRICE::scheduleSetBC_FineLevel \t\t\tL-" << L_indx << " P-"
               << *patches << '\n';

    Task* t;
    t = scinew Task("AMRICE::setBC_FineLevel", this, &AMRICE::setBC_FineLevel);
    Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice
                                                       // matlSet.
    Ghost::GhostType gn         = Ghost::None;
    Task::MaterialDomainSpec ND = Task::NormalDomain;

    // need to interpolate these intermediate values
    t->needs(Task::NewDW,
                d_ice_labels->gammaLabel,
                0,
                Task::CoarseLevel,
                0,
                ND,
                gn,
                0);
    t->needs(Task::NewDW,
                d_ice_labels->specific_heatLabel,
                0,
                Task::CoarseLevel,
                0,
                ND,
                gn,
                0);
    t->needs(Task::NewDW,
                d_ice_labels->vol_frac_CCLabel,
                0,
                Task::CoarseLevel,
                0,
                ND,
                gn,
                0);

    const MaterialSubset* all_matls =
      d_materialManager->allMaterials()->getUnion();

    t->needs(Task::OldDW, d_ice_labels->timeStepLabel);
    t->needs(Task::OldDW, d_ice_labels->simulationTimeLabel);
    t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));

    t->modifies(d_ice_labels->press_CCLabel, d_press_matl, oims);
    t->modifies(d_ice_labels->rho_CCLabel);
    t->modifies(d_ice_labels->specificVolume_CCLabel);
    t->modifies(d_ice_labels->temperature_CCLabel);
    t->modifies(d_ice_labels->velocity_CCLabel);

    // we really only do the ice matls, but we need to tell CopyData to do the
    // right thing
    t->computes(d_ice_labels->gammaLabel, all_matls, oims);
    t->computes(d_ice_labels->specific_heatLabel, all_matls, oims);
    t->computes(d_ice_labels->vol_frac_CCLabel, all_matls, oims);

    //__________________________________
    // Model with transported variables.
    if (d_models.size()) {
      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());

        if (fb_model && fb_model->d_transVars.size()) {
          std::vector<TransportedVariable*>::iterator t_iter;

          for (t_iter = fb_model->d_transVars.begin();
               t_iter != fb_model->d_transVars.end();
               t_iter++) {
            TransportedVariable* tvar = *t_iter;

            t->modifies(tvar->var);
          }
        }
      }
    }
    sched->addTask(t, patches, d_materialManager->allMaterials("ICE"));
  }
}

/*______________________________________________________________________
 Function~  AMRICE::setBC_FineLevel
 Purpose~   set the boundary conditions on the fine level at the edge
 of the computational domain
______________________________________________________________________*/
void
AMRICE::setBC_FineLevel(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset*,
                        DataWarehouse* fine_old_dw,
                        DataWarehouse* fine_new_dw)
{
  timeStep_vartype timeStep;
  fine_old_dw->get(timeStep, d_ice_labels->timeStepLabel);

  bool isNotInitialTimestep = (timeStep > 0);

  const Level* fineLevel   = getLevel(patches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();

  if (fineLevel->getIndex() > 0) {
    cout_doing << d_myworld->myRank() << " Doing setBC_FineLevel"
               << "\t\t\t\t AMRICE L-" << fineLevel->getIndex()
               << " Patches: " << *patches << endl;

    int numICEMatls = d_materialManager->getNumMaterials("ICE");
    //    bool dbg_onOff = cout_dbg.active();      // is cout_dbg switch on or
    //    off (turn off for threaded scheduler)

    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      std::vector<CCVariable<double>> sp_vol_CC(numICEMatls);
      std::vector<constCCVariable<double>> sp_vol_const(numICEMatls);

      for (int m = 0; m < numICEMatls; m++) {
        ICEMaterial* matl =
          (ICEMaterial*)d_materialManager->getMaterial("ICE", m);
        int indx = matl->getDWIndex();
        CCVariable<double> rho_CC, temp_CC, cv, gamma, vol_frac;
        CCVariable<Vector> vel_CC;

        fine_new_dw->getModifiable(
          sp_vol_CC[m], d_ice_labels->specificVolume_CCLabel, indx, patch);
        fine_new_dw->getModifiable(
          rho_CC, d_ice_labels->rho_CCLabel, indx, patch);
        fine_new_dw->getModifiable(
          temp_CC, d_ice_labels->temperature_CCLabel, indx, patch);
        fine_new_dw->getModifiable(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch);
        fine_new_dw->allocateAndPut(
          gamma, d_ice_labels->gammaLabel, indx, patch);
        fine_new_dw->allocateAndPut(
          cv, d_ice_labels->specific_heatLabel, indx, patch);
        fine_new_dw->allocateAndPut(
          vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch);

        //__________________________________
        // interpolate the intermediate variables (cv, gamma,vol_frac)
        // to the finer level along the boundary edge
        // Assumption: cv, gamma, vol_frac on the coarest level are accurate
        // enough
        cv.initialize(d_EVIL_NUM);
        gamma.initialize(d_EVIL_NUM);
        vol_frac.initialize(d_EVIL_NUM);

        IntVector refineRatio = fineLevel->getRefinementRatio();
        DataWarehouse* coarse_new_dw =
          fine_new_dw->getOtherDataWarehouse(Task::CoarseNewDW);

        int orderOfInterpolation = 0;
        //__________________________________
        // Iterate over fine level boundary faces
        std::vector<Patch::FaceType>::const_iterator iter;
        std::vector<Patch::FaceType> bf;
        patch->getBoundaryFaces(bf);

        for (iter = bf.begin(); iter != bf.end(); ++iter) {
          Patch::FaceType face = *iter;
          cout_dbg << " Setting BC on Face " << face << " patch "
                   << patch->getGridIndex() << " Level "
                   << fineLevel->getIndex() << endl;
          //__________________________________
          // fine level hi & lo cell iter limits
          // coarselevel hi and low index
          IntVector cl, ch, fl, fh;
          getCoarseFineFaceRange(patch,
                                 coarseLevel,
                                 face,
                                 Patch::ExtraPlusEdgeCells,
                                 orderOfInterpolation,
                                 cl,
                                 ch,
                                 fl,
                                 fh);

          constCCVariable<double> cv_coarse, gamma_coarse, vol_frac_coarse;
          coarse_new_dw->getRegion(cv_coarse,
                                   d_ice_labels->specific_heatLabel,
                                   indx,
                                   coarseLevel,
                                   cl,
                                   ch);
          coarse_new_dw->getRegion(
            gamma_coarse, d_ice_labels->gammaLabel, indx, coarseLevel, cl, ch);
          coarse_new_dw->getRegion(vol_frac_coarse,
                                   d_ice_labels->vol_frac_CCLabel,
                                   indx,
                                   coarseLevel,
                                   cl,
                                   ch);

          selectInterpolator(cv_coarse,
                             orderOfInterpolation,
                             coarseLevel,
                             fineLevel,
                             refineRatio,
                             fl,
                             fh,
                             cv);

          selectInterpolator(gamma_coarse,
                             orderOfInterpolation,
                             coarseLevel,
                             fineLevel,
                             refineRatio,
                             fl,
                             fh,
                             gamma);

          selectInterpolator(vol_frac_coarse,
                             orderOfInterpolation,
                             coarseLevel,
                             fineLevel,
                             refineRatio,
                             fl,
                             fh,
                             vol_frac);
        } // boundary face loop

        std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
          std::make_unique<CustomBCDriver::customBC_localVars>();
#if 0
        // Worry about this later
        // the problem is that you don't know have delT for the finer level at this point in the cycle
        preprocess_CustomBCs("setBC_FineLevel",fine_old_dw, fine_new_dw, lb,  patch, 999,
                             d_BC_globalVars.get(), BC_localVars.get(), isNotInitialTimestep);
#endif

        constCCVariable<double> placeHolder;

        setBC(rho_CC,
              "Density",
              placeHolder,
              placeHolder,
              patch,
              d_materialManager,
              indx,
              fine_new_dw,
              d_BC_globalVars.get(),
              BC_localVars.get(),
              isNotInitialTimestep);

        setBC(vel_CC,
              "Velocity",
              patch,
              d_materialManager,
              indx,
              fine_new_dw,
              d_BC_globalVars.get(),
              BC_localVars.get(),
              isNotInitialTimestep);

        setBC(temp_CC,
              "Temperature",
              gamma,
              cv,
              patch,
              d_materialManager,
              indx,
              fine_new_dw,
              d_BC_globalVars.get(),
              BC_localVars.get(),
              isNotInitialTimestep);

        setSpecificVolBC(sp_vol_CC[m],
                         "SpecificVol",
                         false,
                         rho_CC,
                         vol_frac,
                         patch,
                         d_materialManager,
                         indx);

        sp_vol_const[m] = sp_vol_CC[m]; // needed by pressure BC

        // worry about this later
        delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());

        //__________________________________
        // Model with transported variables.
        if (d_models.size()) {
          for (auto& model : d_models) {
            FluidsBasedModel* fb_model =
              dynamic_cast<FluidsBasedModel*>(model.get());

            if (fb_model && fb_model->d_transVars.size()) {
              std::vector<TransportedVariable*>::iterator t_iter;

              for (t_iter = fb_model->d_transVars.begin();
                   t_iter != fb_model->d_transVars.end();
                   t_iter++) {
                TransportedVariable* tvar = *t_iter;

                if (tvar->matls->contains(indx)) {

                  std::string q_CC_name = tvar->var->getName();
                  CCVariable<double> q_CC;
                  fine_new_dw->getModifiable(q_CC, tvar->var, indx, patch);

                  setBC(q_CC,
                        q_CC_name,
                        patch,
                        d_materialManager,
                        indx,
                        fine_new_dw,
                        isNotInitialTimestep);
                }
              }
            }
          }
        }
      } // matl loop

      //__________________________________
      //  Pressure boundary condition
      CustomBCDriver::customBC_localVars* notUsed =
        scinew CustomBCDriver::customBC_localVars();

      CCVariable<double> press_CC;
      std::vector<CCVariable<double>> placeHolder(0);

      fine_new_dw->getModifiable(
        press_CC, d_ice_labels->press_CCLabel, 0, patch);

      setBC(press_CC,
            placeHolder,
            sp_vol_const,
            d_surroundingMatl_indx,
            "sp_vol",
            "Pressure",
            patch,
            d_materialManager,
            0,
            fine_new_dw,
            d_BC_globalVars.get(),
            notUsed,
            isNotInitialTimestep);

      delete_CustomBCs(d_BC_globalVars.get(), notUsed);

    } // patches loop
      //  cout_dbg.setActive(dbg_onOff);  // reset on/off switch for cout_dbg,
      //  (turn off for tsanitizer warnings)
  }
}

/*___________________________________________________________________
 Function~  AMRICE::scheduleRefine--
_____________________________________________________________________*/
void
AMRICE::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx             = fineLevel->getIndex();

  if (L_indx > 0) {

    cout_doing << d_myworld->myRank() << " AMRICE::scheduleRefine\t\t\t\tL-"
               << L_indx << " P-" << *patches << '\n';
    Task* task = scinew Task("AMRICE::refine", this, &AMRICE::refine);

    MaterialSubset* subset = scinew MaterialSubset;

    subset->add(0);

    task->needs(Task::NewDW,
                   d_ice_labels->press_CCLabel,
                   0,
                   Task::CoarseLevel,
                   subset,
                   Task::OutOfDomain,
                   d_gac,
                   1);

    task->needs(Task::NewDW,
                   d_ice_labels->rho_CCLabel,
                   0,
                   Task::CoarseLevel,
                   0,
                   Task::NormalDomain,
                   d_gac,
                   1);

    task->needs(Task::NewDW,
                   d_ice_labels->specificVolume_CCLabel,
                   0,
                   Task::CoarseLevel,
                   0,
                   Task::NormalDomain,
                   d_gac,
                   1);

    task->needs(Task::NewDW,
                   d_ice_labels->temperature_CCLabel,
                   0,
                   Task::CoarseLevel,
                   0,
                   Task::NormalDomain,
                   d_gac,
                   1);

    task->needs(Task::NewDW,
                   d_ice_labels->velocity_CCLabel,
                   0,
                   Task::CoarseLevel,
                   0,
                   Task::NormalDomain,
                   d_gac,
                   1);

    //__________________________________
    // Models with transported variables.
    if (d_models.size()) {
      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());

        if (fb_model && fb_model->d_transVars.size()) {
          std::vector<TransportedVariable*>::iterator t_iter;

          for (t_iter = fb_model->d_transVars.begin();
               t_iter != fb_model->d_transVars.end();
               t_iter++) {
            TransportedVariable* tvar = *t_iter;

            task->needs(Task::NewDW,
                           tvar->var,
                           0,
                           Task::CoarseLevel,
                           0,
                           Task::NormalDomain,
                           d_gac,
                           1);
            task->computes(tvar->var);
          }
        }
      }
    }

    //__________________________________
    // Models that need to refine/initialize variables on new patches.
    // This will call both ICE and MPMICE based models
    for (auto& model : d_models) {
      model->scheduleRefine(patches, sched);
    }

    task->computes(d_ice_labels->press_CCLabel, subset, Task::OutOfDomain);
    task->computes(d_ice_labels->rho_CCLabel);
    task->computes(d_ice_labels->specificVolume_CCLabel);
    task->computes(d_ice_labels->temperature_CCLabel);
    task->computes(d_ice_labels->velocity_CCLabel);
    sched->addTask(task, patches, d_materialManager->allMaterials("ICE"));

    //__________________________________
    // Sub Task
    scheduleSetBC_FineLevel(patches, sched);
  }
}

/*___________________________________________________________________
This task initializes the variables on all patches that the regridder
creates.  The BNR and Hierarchical regridders will create patches that
are partially filled with old data.  We don't
want to overwrite these data, thus only use the tiled regridder

_____________________________________________________________________*/
void
AMRICE::refine(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* matls,
               DataWarehouse*,
               DataWarehouse* new_dw)
{

  const Level* fineLevel   = getLevel(patches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();

  cout_doing << d_myworld->myRank() << " Doing refine \t\t\t\t\t AMRICE L-"
             << fineLevel->getIndex();
  IntVector rr(fineLevel->getRefinementRatio());
  double invRefineRatio = 1. / (rr.x() * rr.y() * rr.z());

  for (int p = 0; p < patches->size(); p++) {
    const Patch* finePatch = patches->get(p);
    cout_doing << "  patch " << finePatch->getID() << endl;

    Level::selectType coarsePatches;
    finePatch->getCoarseLevelPatches(coarsePatches);

    // bullet proofing
    iteratorTest(finePatch, fineLevel, coarseLevel, new_dw);

/*`==========TESTING==========*/
#if 0
    Patch::FaceType notUsed;
    testInterpolators<double>(new_dw,d_orderOfInterpolation,coarseLevel,fineLevel,
                                finePatch, notUsed, "wholeDomain");
#endif
    /*===========TESTING==========`*/

    // refine pressure
    CCVariable<double> press_CC;
    new_dw->allocateAndPut(press_CC, d_ice_labels->press_CCLabel, 0, finePatch);
    press_CC.initialize(d_EVIL_NUM);
    CoarseToFineOperator<double>(press_CC,
                                 d_ice_labels->press_CCLabel,
                                 0,
                                 new_dw,
                                 invRefineRatio,
                                 finePatch,
                                 fineLevel,
                                 coarseLevel);

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);
      CCVariable<double> rho_CC, temp, sp_vol_CC;
      CCVariable<Vector> vel_CC;

      new_dw->allocateAndPut(
        rho_CC, d_ice_labels->rho_CCLabel, indx, finePatch);
      new_dw->allocateAndPut(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, finePatch);
      new_dw->allocateAndPut(
        temp, d_ice_labels->temperature_CCLabel, indx, finePatch);
      new_dw->allocateAndPut(
        vel_CC, d_ice_labels->velocity_CCLabel, indx, finePatch);

      rho_CC.initialize(d_EVIL_NUM);
      sp_vol_CC.initialize(d_EVIL_NUM);
      temp.initialize(d_EVIL_NUM);
      vel_CC.initialize(Vector(d_EVIL_NUM));

      // refine
      CoarseToFineOperator<double>(rho_CC,
                                   d_ice_labels->rho_CCLabel,
                                   indx,
                                   new_dw,
                                   invRefineRatio,
                                   finePatch,
                                   fineLevel,
                                   coarseLevel);

      CoarseToFineOperator<double>(sp_vol_CC,
                                   d_ice_labels->specificVolume_CCLabel,
                                   indx,
                                   new_dw,
                                   invRefineRatio,
                                   finePatch,
                                   fineLevel,
                                   coarseLevel);

      CoarseToFineOperator<double>(temp,
                                   d_ice_labels->temperature_CCLabel,
                                   indx,
                                   new_dw,
                                   invRefineRatio,
                                   finePatch,
                                   fineLevel,
                                   coarseLevel);

      CoarseToFineOperator<Vector>(vel_CC,
                                   d_ice_labels->velocity_CCLabel,
                                   indx,
                                   new_dw,
                                   invRefineRatio,
                                   finePatch,
                                   fineLevel,
                                   coarseLevel);

      //__________________________________
      // Model with transported variables.
      if (d_models.size()) {
        for (auto& model : d_models) {
          FluidsBasedModel* fb_model =
            dynamic_cast<FluidsBasedModel*>(model.get());

          if (fb_model && fb_model->d_transVars.size()) {
            std::vector<TransportedVariable*>::iterator t_iter;

            for (t_iter = fb_model->d_transVars.begin();
                 t_iter != fb_model->d_transVars.end();
                 t_iter++) {
              TransportedVariable* tvar = *t_iter;

              if (tvar->matls->contains(indx)) {
                CCVariable<double> q_CC;
                new_dw->allocateAndPut(q_CC, tvar->var, indx, finePatch);

                q_CC.initialize(d_EVIL_NUM);

                CoarseToFineOperator<double>(q_CC,
                                             tvar->var,
                                             indx,
                                             new_dw,
                                             invRefineRatio,
                                             finePatch,
                                             fineLevel,
                                             coarseLevel);
              }
            }
          }
        }
      }
    }
  } // course patch loop
}

/*_____________________________________________________________________
 Function~  AMRICE::iteratorTest--
 Purpose~   Verify that the all of the fine level cells will be accessed
_____________________________________________________________________*/
void
AMRICE::iteratorTest(const Patch* finePatch,
                     const Level* fineLevel,
                     const Level* coarseLevel,
                     DataWarehouse* new_dw)
{
  Level::selectType coarsePatches;
  finePatch->getCoarseLevelPatches(coarsePatches);
  IntVector fl = finePatch->getExtraCellLowIndex();
  IntVector fh = finePatch->getExtraCellHighIndex();

  CCVariable<double> hitCells;
  new_dw->allocateTemporary(hitCells, finePatch);
  hitCells.initialize(d_EVIL_NUM);
  IntVector lo(0, 0, 0);
  IntVector hi(0, 0, 0);

  for (size_t i = 0; i < coarsePatches.size(); i++) {
    const Patch* coarsePatch = coarsePatches[i];
    // iterator should hit the cells over the intersection of the fine and
    // coarse patches

    IntVector cl = coarsePatch->getCellLowIndex();
    IntVector ch = coarsePatch->getCellHighIndex();

    IntVector fl_tmp = coarseLevel->mapCellToFiner(cl);
    IntVector fh_tmp = coarseLevel->mapCellToFiner(ch);

    lo = Max(fl, fl_tmp);
    hi = Min(fh, fh_tmp);

    for (CellIterator iter(lo, hi); !iter.done(); iter++) {
      IntVector c = *iter;
      hitCells[c] = 1.0;
    }
#if 0
    cout << " coarsePatch.size() " << coarsePatches.size()
         << " coarsePatch " << coarsePatch->getID()
         << " finePatch " << finePatch->getID()
         << " fineLevel: fl " << fl << " fh " << fh
         << " coarseLevel: cl " << cl << " ch " << ch
         << " final Iterator: " << lo << " " << hi << endl;
#endif
  }

  //____ B U L L E T   P R O O F I N G_______
  // All cells must be initialized at this point
  IntVector badCell;
  CellIterator iter(lo, hi);
  if (isEqual<double>(d_EVIL_NUM, iter, hitCells, badCell)) {

    IntVector c_badCell = fineLevel->mapCellToCoarser(badCell);
    const Patch* patch  = coarseLevel->selectPatchForCellIndex(c_badCell);

    std::ostringstream warn;
    warn << "ERROR AMRICE::Refine Task:iteratorTest "
         << "detected an fine level cell that won't get initialized " << badCell
         << " Patch " << finePatch->getID() << " Level idx "
         << fineLevel->getIndex() << "\n "
         << "The underlying coarse cell " << c_badCell
         << " belongs to coarse level patch " << patch->getID() << "\n";
    throw InvalidValue(warn.str(), __FILE__, __LINE__);
  }
}

/*___________________________________________________________________
 Function~  AMRICE::scheduleCoarsen--
_____________________________________________________________________*/
void
AMRICE::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  Ghost::GhostType gn    = Ghost::None;
  cout_doing << d_myworld->myRank() << " AMRICE::scheduleCoarsen\t\t\t\tL-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
             << endl;

  Task* task = scinew Task("AMRICE::coarsen", this, &AMRICE::coarsen);

  Task::MaterialDomainSpec ND   = Task::NormalDomain;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  const MaterialSet* all_matls  = d_materialManager->allMaterials();
  const MaterialSet* ice_matls  = d_materialManager->allMaterials("ICE");

  const MaterialSubset* all_matls_sub = all_matls->getUnion();
  const PatchSet* patch_set           = coarseLevel->eachPatch();

  Task::SearchTG OldTG =
    Task::SearchTG::OldTG; // possibly search old TG for computes

  task->needs(Task::NewDW,
                 d_ice_labels->press_CCLabel,
                 0,
                 Task::FineLevel,
                 d_press_matl,
                 oims,
                 gn,
                 0,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mass_advLabel,
                 0,
                 Task::FineLevel,
                 all_matls_sub,
                 ND,
                 gn,
                 0,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->sp_vol_advLabel,
                 0,
                 Task::FineLevel,
                 all_matls_sub,
                 ND,
                 gn,
                 0,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->eng_advLabel,
                 0,
                 Task::FineLevel,
                 all_matls_sub,
                 ND,
                 gn,
                 0,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mom_advLabel,
                 0,
                 Task::FineLevel,
                 all_matls_sub,
                 ND,
                 gn,
                 0,
                 OldTG);

  //__________________________________
  // Model with transported variables.
  if (d_models.size()) {
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());

      if (fb_model && fb_model->d_transVars.size()) {
        std::vector<TransportedVariable*>::iterator t_iter;

        for (t_iter = fb_model->d_transVars.begin();
             t_iter != fb_model->d_transVars.end();
             t_iter++) {
          TransportedVariable* tvar = *t_iter;

          task->needs(Task::NewDW,
                         tvar->var_adv,
                         0,
                         Task::FineLevel,
                         all_matls_sub,
                         ND,
                         gn,
                         0,
                         OldTG);

          task->modifies(tvar->var_adv, OldTG);
        }
      }
    }
  }

  task->modifies(d_ice_labels->press_CCLabel, d_press_matl, oims, OldTG);
  task->modifies(d_ice_labels->mass_advLabel, OldTG);
  task->modifies(d_ice_labels->sp_vol_advLabel, OldTG);
  task->modifies(d_ice_labels->eng_advLabel, OldTG);
  task->modifies(d_ice_labels->mom_advLabel, OldTG);

  sched->addTask(task, patch_set, ice_matls);

  //__________________________________
  // schedule refluxing and bulletproofing
  // the bulletproofing tasks have no computes or requires
  if (d_doRefluxing) {
    Task* t;
    Task* t1;

    scheduleReflux_computeCorrectionFluxes(coarseLevel,
                                           sched); // compute correction

    //__________________________________
    //  initialize the bullet proofing flags
    t = scinew Task("AMRICE::reflux_BP_zero_CFI_cells",
                    this,
                    &AMRICE::reflux_BP_zero_CFI_cells);
    sched->addTask(t, coarseLevel->eachPatch(), ice_matls);

    scheduleReflux_applyCorrection(coarseLevel, sched); // apply correction

    //__________________________________
    // check the bullet proofing flags
    t = scinew Task("AMRICE::reflux_BP_count_CFI_cells",
                    this,
                    &AMRICE::reflux_BP_count_CFI_cells);

    std::string desc2 = "applyRefluxCorrection";
    t1                = scinew Task("AMRICE::reflux_BP_check_CFI_cells",
                     this,
                     &AMRICE::reflux_BP_check_CFI_cells,
                     desc2);

    sched->addTask(t, patch_set, ice_matls);
    sched->addTask(t1, patch_set, ice_matls);
  }
}

/*___________________________________________________________________
 Function~  AMRICE::Coarsen--
_____________________________________________________________________*/
void
AMRICE::coarsen(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                DataWarehouse*,
                DataWarehouse* new_dw)
{

  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();
  cout_doing << d_myworld->myRank() << " Doing coarsen \t\t\t\t\t AMRICE L-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex();

  //  bool dbg_onOff = cout_dbg.active();      // is cout_dbg switch on or off
  //  (turn off for threaded scheduler)

  for (int p = 0; p < patches->size(); p++) {
    const Patch* coarsePatch = patches->get(p);
    cout_doing << "  patch " << coarsePatch->getID() << endl;

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      CCVariable<double> mass_adv, eng_adv, sp_vol_adv;
      CCVariable<Vector> mom_adv;

      new_dw->getModifiable(
        mass_adv, d_ice_labels->mass_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        sp_vol_adv, d_ice_labels->sp_vol_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        eng_adv, d_ice_labels->eng_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        mom_adv, d_ice_labels->mom_advLabel, indx, coarsePatch);

      // coarsen
      bool computesAve = false;
      AMRCoarsenRefine::fineToCoarseOperator<double>(
        mass_adv,
        computesAve,
        d_ice_labels->mass_advLabel,
        indx,
        new_dw,
        coarsePatch,
        coarseLevel,
        fineLevel);

      AMRCoarsenRefine::fineToCoarseOperator<double>(
        sp_vol_adv,
        computesAve,
        d_ice_labels->sp_vol_advLabel,
        indx,
        new_dw,
        coarsePatch,
        coarseLevel,
        fineLevel);

      AMRCoarsenRefine::fineToCoarseOperator<double>(eng_adv,
                                                     computesAve,
                                                     d_ice_labels->eng_advLabel,
                                                     indx,
                                                     new_dw,
                                                     coarsePatch,
                                                     coarseLevel,
                                                     fineLevel);

      AMRCoarsenRefine::fineToCoarseOperator<Vector>(mom_adv,
                                                     computesAve,
                                                     d_ice_labels->mom_advLabel,
                                                     indx,
                                                     new_dw,
                                                     coarsePatch,
                                                     coarseLevel,
                                                     fineLevel);

      //__________________________________
      // pressure
      if (indx == 0) {
        // pressure
        CCVariable<double> press_CC;
        new_dw->getModifiable(
          press_CC, d_ice_labels->press_CCLabel, 0, coarsePatch);
        computesAve = true;

        AMRCoarsenRefine::fineToCoarseOperator<double>(
          press_CC,
          computesAve,
          d_ice_labels->press_CCLabel,
          0,
          new_dw,
          coarsePatch,
          coarseLevel,
          fineLevel);
      }

      //__________________________________
      // Model with transported variables.
      if (d_models.size()) {
        for (auto& model : d_models) {
          FluidsBasedModel* fb_model =
            dynamic_cast<FluidsBasedModel*>(model.get());

          if (fb_model && fb_model->d_transVars.size()) {
            std::vector<TransportedVariable*>::iterator t_iter;

            for (t_iter = fb_model->d_transVars.begin();
                 t_iter != fb_model->d_transVars.end();
                 t_iter++) {
              TransportedVariable* tvar = *t_iter;

              if (tvar->matls->contains(indx)) {
                CCVariable<double> q_CC_adv;

                new_dw->getModifiable(
                  q_CC_adv, tvar->var_adv, indx, coarsePatch);
                computesAve = false;

                AMRCoarsenRefine::fineToCoarseOperator<double>(q_CC_adv,
                                                               computesAve,
                                                               tvar->var_adv,
                                                               indx,
                                                               new_dw,
                                                               coarsePatch,
                                                               coarseLevel,
                                                               fineLevel);
              }
            }
          }
        }
      }
    }
  } // course patch loop
  //  cout_dbg.setActive(dbg_onOff);  // reset on/off switch for cout_dbg (turn
  //  off for tsanitizer warnings)
}
/*___________________________________________________________________
 Function~  AMRICE::scheduleReflux_computeCorrectionFluxes--
_____________________________________________________________________*/
void
AMRICE::scheduleReflux_computeCorrectionFluxes(const LevelP& coarseLevel,
                                               SchedulerP& sched)
{
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  cout_doing << d_myworld->myRank()
             << " AMRICE::scheduleReflux_computeCorrectionFluxes\tL-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
             << endl;

  Task* task = scinew Task("AMRICE::reflux_computeCorrectionFluxes",
                           this,
                           &AMRICE::reflux_computeCorrectionFluxes);

  Ghost::GhostType gx = Ghost::AroundFacesX;
  Ghost::GhostType gy = Ghost::AroundFacesY;
  Ghost::GhostType gz = Ghost::AroundFacesZ;
  Task::SearchTG OldTG =
    Task::SearchTG::OldTG; // possibly search oldTG for computes

  //__________________________________
  // Fluxes from the fine level
  // MASS
  task->needs(Task::NewDW,
                 d_ice_labels->mass_X_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gx,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mass_Y_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gy,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mass_Z_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gz,
                 1,
                 OldTG);

  // MOMENTUM
  task->needs(Task::NewDW,
                 d_ice_labels->mom_X_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gx,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mom_Y_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gy,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->mom_Z_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gz,
                 1,
                 OldTG);

  // INT_ENG
  task->needs(Task::NewDW,
                 d_ice_labels->int_eng_X_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gx,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->int_eng_Y_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gy,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->int_eng_Z_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gz,
                 1,
                 OldTG);

  // SPECIFIC VOLUME
  task->needs(Task::NewDW,
                 d_ice_labels->sp_vol_X_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gx,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->sp_vol_Y_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gy,
                 1,
                 OldTG);

  task->needs(Task::NewDW,
                 d_ice_labels->sp_vol_Z_FC_fluxLabel,
                 0,
                 Task::FineLevel,
                 0,
                 Task::NormalDomain,
                 gz,
                 1,
                 OldTG);

  //__________________________________
  // Model with reflux variables.
  if (d_models.size()) {
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());

      if (fb_model) {
        for (auto& rvar : fb_model->d_refluxVars) {

          task->needs(Task::NewDW,
                         rvar->var_X_FC_flux,
                         0,
                         Task::FineLevel,
                         0,
                         Task::NormalDomain,
                         gx,
                         1,
                         OldTG);

          task->needs(Task::NewDW,
                         rvar->var_Y_FC_flux,
                         0,
                         Task::FineLevel,
                         0,
                         Task::NormalDomain,
                         gy,
                         1,
                         OldTG);

          task->needs(Task::NewDW,
                         rvar->var_Z_FC_flux,
                         0,
                         Task::FineLevel,
                         0,
                         Task::NormalDomain,
                         gz,
                         1,
                         OldTG);

          task->computes(rvar->var_X_FC_corr);
          task->computes(rvar->var_Y_FC_corr);
          task->computes(rvar->var_Z_FC_corr);
        }
      }
    }
  }

  task->computes(d_ice_labels->mass_X_FC_corrLabel);
  task->computes(d_ice_labels->mass_Y_FC_corrLabel);
  task->computes(d_ice_labels->mass_Z_FC_corrLabel);

  task->computes(d_ice_labels->mom_X_FC_corrLabel);
  task->computes(d_ice_labels->mom_Y_FC_corrLabel);
  task->computes(d_ice_labels->mom_Z_FC_corrLabel);

  task->computes(d_ice_labels->int_eng_X_FC_corrLabel);
  task->computes(d_ice_labels->int_eng_Y_FC_corrLabel);
  task->computes(d_ice_labels->int_eng_Z_FC_corrLabel);

  task->computes(d_ice_labels->sp_vol_X_FC_corrLabel);
  task->computes(d_ice_labels->sp_vol_Y_FC_corrLabel);
  task->computes(d_ice_labels->sp_vol_Z_FC_corrLabel);

  sched->addTask(
    task, coarseLevel->eachPatch(), d_materialManager->allMaterials("ICE"));
}
/*___________________________________________________________________
 Function~  AMRICE::Reflux_computeCorrectionFluxesFluxes--
_____________________________________________________________________*/
void
AMRICE::reflux_computeCorrectionFluxes(const ProcessorGroup*,
                                       const PatchSubset* coarsePatches,
                                       const MaterialSubset* matls,
                                       DataWarehouse*,
                                       DataWarehouse* new_dw)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  cout_doing << d_myworld->myRank()
             << " Doing reflux_computeCorrectionFluxes \t\t\t AMRICE L-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex();

  //  bool dbg_onOff = cout_dbg.active();      // is cout_dbg switch on or off
  //  (turn off for threaded scheduler)

  //__________________________________
  for (int c_p = 0; c_p < coarsePatches->size(); c_p++) {
    const Patch* coarsePatch = coarsePatches->get(c_p);

    cout_doing << "  coarsePatches " << *coarsePatches << endl;

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      Level::selectType finePatches;
      coarsePatch->getOtherLevelPatches(
        1,
        finePatches,
        1); // get with a ghost cell to make sure you get all patches

      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        //__________________________________
        //   compute the correction
        // one_zero:  used to increment the CFI counter.
        if (finePatch->hasCoarseFaces()) {
          cout_doing << d_myworld->myRank() << "  coarsePatch "
                     << coarsePatch->getID() << " finepatch "
                     << finePatch->getID() << endl;

          int one_zero = 1;
          refluxOperator_computeCorrectionFluxes<double>("mass",
                                                         indx,
                                                         coarsePatch,
                                                         finePatch,
                                                         coarseLevel,
                                                         fineLevel,
                                                         new_dw,
                                                         one_zero);

          one_zero = 0;
          refluxOperator_computeCorrectionFluxes<double>("sp_vol",
                                                         indx,
                                                         coarsePatch,
                                                         finePatch,
                                                         coarseLevel,
                                                         fineLevel,
                                                         new_dw,
                                                         one_zero);

          refluxOperator_computeCorrectionFluxes<Vector>("mom",
                                                         indx,
                                                         coarsePatch,
                                                         finePatch,
                                                         coarseLevel,
                                                         fineLevel,
                                                         new_dw,
                                                         one_zero);

          refluxOperator_computeCorrectionFluxes<double>("int_eng",
                                                         indx,
                                                         coarsePatch,
                                                         finePatch,
                                                         coarseLevel,
                                                         fineLevel,
                                                         new_dw,
                                                         one_zero);

          //__________________________________
          // Model with reflux variables.
          if (d_models.size()) {
            for (auto& model : d_models) {
              FluidsBasedModel* fb_model =
                dynamic_cast<FluidsBasedModel*>(model.get());

              if (fb_model) {
                for (auto& rvar : fb_model->d_refluxVars) {

                  if (rvar->matls->contains(indx)) {
                    std::string var_name = rvar->var->getName();
                    refluxOperator_computeCorrectionFluxes<double>(var_name,
                                                                   indx,
                                                                   coarsePatch,
                                                                   finePatch,
                                                                   coarseLevel,
                                                                   fineLevel,
                                                                   new_dw,
                                                                   one_zero);
                  }
                }
              }
            }
          } // model
        }
      } // finePatch loop
    }   // matl loop
  }     // coarse patch loop
  //  cout_dbg.setActive(dbg_onOff);  // reset on/off switch for cout_dbg, (turn
  //  off for tsanitizer warnings)
}

/*___________________________________________________________________
 Function~  ICE::refluxCoarseLevelIterator--
 Purpose:  returns the iterator and face-centered offset that the coarse
           level uses to do refluxing.  THIS IS COMPILCATED AND CONFUSING
_____________________________________________________________________*/
void
AMRICE::refluxCoarseLevelIterator(Patch::FaceType patchFace,
                                  const Patch* coarsePatch,
                                  const Patch* finePatch,
                                  const Level* fineLevel,
                                  CellIterator& iter,
                                  IntVector& coarse_FC_offset,
                                  bool& isRight_CP_FP_pair,
                                  const std::string& whichTask)
{
  Patch::FaceIteratorType IFC = Patch::InteriorFaceCells;
  CellIterator f_iter         = finePatch->getFaceIterator(patchFace, IFC);

  ASSERT(whichTask == "computeRefluxCorrection" ||
         whichTask == "applyRefluxCorrection");

  // find the intersection of the fine patch face iterator and underlying coarse
  // patch
  IntVector dir       = finePatch->getFaceAxes(patchFace); // face axes
  int p_dir           = dir[0];                            // normal direction
  IntVector f_lo_face = f_iter.begin(); // fineLevel face indices
  IntVector f_hi_face = f_iter.end();

  IntVector l = fineLevel->mapCellToCoarser(f_lo_face);
  IntVector h = fineLevel->mapCellToCoarser(f_hi_face);

  //__________________________________
  // Offset for the coarse level iterator
  // shift l & h,  1 cell for x+, y+, z+ finePatchfaces
  // shift l & h, -1 cell for x-, y-, z- finePatchfaces
  //
  // if(minus face){
  //    if(refinement ratio == 1)
  //      shift h
  //    else
  //      no shift for h
  // }
  //
  // Offset for the coarse_face center (c_FC_offset)
  // shift 1 cell   for x-, y-, z- finePatch faces
  // shift 0 cells  for x+, y+, z+ finePatchFaces

  std::string name = finePatch->getFaceName(patchFace);
  IntVector offset = finePatch->faceDirection(patchFace);
  coarse_FC_offset = IntVector(0, 0, 0);

  if (name == "xminus" || name == "yminus" || name == "zminus") {
    coarse_FC_offset = -offset;
    l += offset;

    IntVector rr(fineLevel->getRefinementRatio());
    if (rr[p_dir] == 1) {
      h += offset;
    }
  }
  if (name == "xplus" || name == "yplus" || name == "zplus") {
    l += offset;
    h += offset;
  }

  IntVector coarse_Lo = coarsePatch->getExtraCellLowIndex();
  IntVector coarse_Hi = coarsePatch->getExtraCellHighIndex();
  int y               = dir[1]; // tangential directions
  int z               = dir[2];

/*`==========TESTING==========*/
#if 0
  if(finePatch->getID() == 583){
    cout << "\nrefluxCoarseLevelIterator " << name << " " << whichTask
         << "\n before      " << l << " " << h
         << "\n finePatch   " << *finePatch
         << "\n coarsePatch " << *coarsePatch
         << "\n coarseLevel " << coarsePatch->getLevel()->getIndex()
         << "\n coarse_Lo   " << coarse_Lo << " coarse_Hi " << coarse_Hi << endl;
  }
#endif
  /*===========TESTING==========`*/

  l[y] = Max(l[y], coarse_Lo[y]); // intersection
  l[z] = Max(l[z], coarse_Lo[z]); // only the transerse directions
  h[y] = Min(h[y], coarse_Hi[y]);
  h[z] = Min(h[z], coarse_Hi[z]);

  iter = CellIterator(l, h);

  // Does this iterator exceed the boundaries of the underlying coarse patch?
  // If your computing/applying the reflux correction the conditional is
  // different To understand why you need paper a pencil.  Draw a fine patch
  // over a coarse patch and let the CFI coinside with the boundary  between 2
  // coarse level patches.

  IntVector one(1, 1, 1);         // subtract of 1 (h).
  IntVector h_minusOne = h - one; // h -1 is what we're truely interested in.

  isRight_CP_FP_pair = false;
  if (whichTask == "computeRefluxCorrection" &&
      coarsePatch->containsCell(l + coarse_FC_offset) &&
      coarsePatch->containsCell(h_minusOne + coarse_FC_offset) &&
      l[y] != h[y] && l[z] != h[z]) {
    isRight_CP_FP_pair = true;
  }
  if (whichTask == "applyRefluxCorrection" && coarsePatch->containsCell(l) &&
      coarsePatch->containsCell(h_minusOne) && l[y] != h[y] && l[z] != h[z]) {
    isRight_CP_FP_pair = true;
  }

  /*`==========TESTING==========*/
#if 0
  if(finePatch->getID() == 583){
    cout << " after " << l << " " << h
         << " coarse_FC_offset " << coarse_FC_offset
         << " isRight_CP_FP_pair " << isRight_CP_FP_pair
         << "\ncoarsePatch->containsCell(l)                             " << coarsePatch->containsCell(l)
         << "\ncoarsePatch->containsCell(h_minusOne)                    " << coarsePatch->containsCell(h_minusOne)
         << "\ncoarsePatch->containsCell(l + coarse_FC_offset)          "<< coarsePatch->containsCell(l + coarse_FC_offset)
         << "\ncoarsePatch->containsCell(h_minusOne + coarse_FC_offset) " << coarsePatch->containsCell(h_minusOne + coarse_FC_offset)
         << " \nl[y] != h[y] && l[z] != h[z]                            " << (l[y] != h[y] && l[z] != h[z])<< endl;
  }
#endif
  /*===========TESTING==========`*/

  //____ B U L L E T   P R O O F I N G----
  if (isRight_CP_FP_pair) {
    IntVector diff = Abs(l - h);
    if ((l.x() >= h.x() || l.y() >= h.y() || l.z() >= h.z()) ||
        diff[p_dir] > 1) {
      std::ostringstream warn;
      warn << "AMRICE:refluxCoarseLevelIterator : " << l << " " << h << " "
           << name << "\n  Error:Either l >= h OR l - h > 1"
           << "\n finelevel   " << fineLevel->getIndex() << "\n finePatch   "
           << *finePatch << "\n CoarsePatch " << *coarsePatch
           << "\n coarseLo    " << coarse_Lo << " coarseHi " << coarse_Hi
           << "\n offset      " << offset << " unmodified Iterator "
           << f_iter.begin() << " " << f_iter.end();
      finePatch->printPatchBCs(warn);

      throw InternalError(warn.str(), __FILE__, __LINE__);
    }
  }
}

/*___________________________________________________________________
 Function~  AMRICE::scheduleReflux_applyCorrection--
_____________________________________________________________________*/
void
AMRICE::scheduleReflux_applyCorrection(const LevelP& coarseLevel,
                                       SchedulerP& sched)
{
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  cout_doing << d_myworld->myRank()
             << " AMRICE::scheduleReflux_applyCorrectionFluxes\t\tL-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
             << endl;

  Task* task = scinew Task("AMRICE::reflux_applyCorrectionFluxes",
                           this,
                           &AMRICE::reflux_applyCorrectionFluxes);

  Ghost::GhostType gac = Ghost::AroundCells;

  //__________________________________
  // Correction fluxes  from the coarse level
  // MASS
  task->needs(Task::NewDW, d_ice_labels->mass_X_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->mass_Y_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->mass_Z_FC_corrLabel, gac, 1);
  // MOMENTUM
  task->needs(Task::NewDW, d_ice_labels->mom_X_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->mom_Y_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->mom_Z_FC_corrLabel, gac, 1);
  // INT_ENG
  task->needs(Task::NewDW, d_ice_labels->int_eng_X_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->int_eng_Y_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->int_eng_Z_FC_corrLabel, gac, 1);
  // SPECIFIC VOLUME
  task->needs(Task::NewDW, d_ice_labels->sp_vol_X_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->sp_vol_Y_FC_corrLabel, gac, 1);
  task->needs(Task::NewDW, d_ice_labels->sp_vol_Z_FC_corrLabel, gac, 1);

  //__________________________________
  // Model with reflux variables.
  if (d_models.size()) {
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());

      if (fb_model) {
        for (auto& rvar : fb_model->d_refluxVars) {

          task->needs(Task::NewDW, rvar->var_X_FC_corr, gac, 1);
          task->needs(Task::NewDW, rvar->var_Y_FC_corr, gac, 1);
          task->needs(Task::NewDW, rvar->var_Z_FC_corr, gac, 1);
          task->modifies(rvar->var_adv);
        }
      }
    }
  }

  task->modifies(d_ice_labels->mass_advLabel);
  task->modifies(d_ice_labels->sp_vol_advLabel);
  task->modifies(d_ice_labels->eng_advLabel);
  task->modifies(d_ice_labels->mom_advLabel);

  sched->addTask(
    task, coarseLevel->eachPatch(), d_materialManager->allMaterials("ICE"));
}
/*___________________________________________________________________
 Function~  AMRICE::Reflux_applyCorrectionFluxes--
_____________________________________________________________________*/
void
AMRICE::reflux_applyCorrectionFluxes(const ProcessorGroup*,
                                     const PatchSubset* coarsePatches,
                                     const MaterialSubset* matls,
                                     DataWarehouse*,
                                     DataWarehouse* new_dw)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  cout_doing << d_myworld->myRank()
             << " Doing reflux_applyCorrectionFluxes \t\t\t AMRICE L-"
             << fineLevel->getIndex() << "->" << coarseLevel->getIndex();

  //  bool dbg_onOff = cout_dbg.active();      // is cout_dbg switch on or off
  //  (turn off for threaded scheduler)

  for (int c_p = 0; c_p < coarsePatches->size(); c_p++) {
    const Patch* coarsePatch = coarsePatches->get(c_p);
    cout_doing << "  patch " << coarsePatch->getID() << endl;

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);
      CCVariable<double> mass_adv, eng_adv, sp_vol_adv;
      CCVariable<Vector> mom_adv;

      // Ghost::GhostType  gn  = Ghost::None;
      new_dw->getModifiable(
        mass_adv, d_ice_labels->mass_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        sp_vol_adv, d_ice_labels->sp_vol_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        eng_adv, d_ice_labels->eng_advLabel, indx, coarsePatch);
      new_dw->getModifiable(
        mom_adv, d_ice_labels->mom_advLabel, indx, coarsePatch);

      Level::selectType finePatches;
      coarsePatch->getOtherLevelPatches(
        1,
        finePatches,
        1); // get with a ghost cell to make sure you get all patches

      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];
        // cout_doing << d_myworld->myRank() << "  coarsePatch " <<
        // coarsePatch->getID() <<" finepatch " << finePatch->getID()<< endl;
        //__________________________________
        //  Apply the correction
        //  one_zero:  used to increment the CFI counter.
        if (finePatch->hasCoarseFaces()) {
          int one_zero = 1;
          refluxOperator_applyCorrectionFluxes<double>(mass_adv,
                                                       "mass",
                                                       indx,
                                                       coarsePatch,
                                                       finePatch,
                                                       coarseLevel,
                                                       fineLevel,
                                                       new_dw,
                                                       one_zero);

          one_zero = 0;
          refluxOperator_applyCorrectionFluxes<double>(sp_vol_adv,
                                                       "sp_vol",
                                                       indx,
                                                       coarsePatch,
                                                       finePatch,
                                                       coarseLevel,
                                                       fineLevel,
                                                       new_dw,
                                                       one_zero);

          refluxOperator_applyCorrectionFluxes<Vector>(mom_adv,
                                                       "mom",
                                                       indx,
                                                       coarsePatch,
                                                       finePatch,
                                                       coarseLevel,
                                                       fineLevel,
                                                       new_dw,
                                                       one_zero);

          refluxOperator_applyCorrectionFluxes<double>(eng_adv,
                                                       "int_eng",
                                                       indx,
                                                       coarsePatch,
                                                       finePatch,
                                                       coarseLevel,
                                                       fineLevel,
                                                       new_dw,
                                                       one_zero);

          //__________________________________
          // Model with reflux variables.
          if (d_models.size()) {
            for (auto& model : d_models) {
              FluidsBasedModel* fb_model =
                dynamic_cast<FluidsBasedModel*>(model.get());

              if (fb_model) {
                for (auto& rvar : fb_model->d_refluxVars) {

                  if (rvar->matls->contains(indx)) {
                    CCVariable<double> q_CC_adv;
                    std::string var_name = rvar->var->getName();

                    new_dw->getModifiable(
                      q_CC_adv, rvar->var_adv, indx, coarsePatch);

                    refluxOperator_applyCorrectionFluxes<double>(q_CC_adv,
                                                                 var_name,
                                                                 indx,
                                                                 coarsePatch,
                                                                 finePatch,
                                                                 coarseLevel,
                                                                 fineLevel,
                                                                 new_dw,
                                                                 one_zero);
                  }
                }
              }
            }
          }
        } // patch has a coarseFineInterface
      }   // finePatch loop
    }     // matl loop
  }       // course patch loop

  //  cout_dbg.setActive(dbg_onOff);  // reset on/off switch for cout_dbg, (turn
  //  off for tsanitizer warnings)
}

/*_____________________________________________________________________
 Function~  AMRICE::reflux_BP_zero_CFI_cells--(bulletproofing)
 Purpose~   Initialze coarse fine interface "marks"
______________________________________________________________________*/
void
AMRICE::reflux_BP_zero_CFI_cells(const ProcessorGroup*,
                                 const PatchSubset* coarsePatches,
                                 const MaterialSubset*,
                                 DataWarehouse*,
                                 DataWarehouse*)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  for (int c_p = 0; c_p < coarsePatches->size(); c_p++) {
    const Patch* coarsePatch = coarsePatches->get(c_p);

    Level::selectType finePatches;
    coarsePatch->getOtherLevelPatches(1, finePatches, 1);

    cout_doing << d_myworld->myRank()
               << " Doing reflux_BP_zero_CFI_cells \t\t\t AMRICE L-"
               << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
               << endl;

    for (size_t p = 0; p < finePatches.size(); p++) {
      const Patch* finePatch = finePatches[p];

      std::vector<Patch::FaceType> cf;
      finePatch->getCoarseFaces(cf);
      std::vector<Patch::FaceType>::const_iterator iter;
      for (iter = cf.begin(); iter != cf.end(); ++iter) {
        Patch::FaceType patchFace = *iter;

        setFaceMark(0, finePatch, patchFace, 0);
        setFaceMark(1, finePatch, patchFace, 0);
      }
    } // finePatch loop
  }   // coarsepatch loop
}

/*_____________________________________________________________________
 Function~  AMRICE::reflux_BP_count_CFI_cells--(bulletproofing)
 Purpose~   count up the number of coarse fine interface cells and save
            that number
______________________________________________________________________*/
void
AMRICE::reflux_BP_count_CFI_cells(const ProcessorGroup*,
                                  const PatchSubset* coarsePatches,
                                  const MaterialSubset*,
                                  DataWarehouse*,
                                  DataWarehouse*)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  for (int c_p = 0; c_p < coarsePatches->size(); c_p++) {
    const Patch* coarsePatch = coarsePatches->get(c_p);

    Level::selectType finePatches;
    coarsePatch->getOtherLevelPatches(1, finePatches, 1);

    cout_doing << d_myworld->myRank()
               << " Doing reflux_BP_count_CFI_cells \t\t\t AMRICE L-"
               << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
               << endl;

    for (size_t p = 0; p < finePatches.size(); p++) {
      const Patch* finePatch = finePatches[p];

      std::vector<Patch::FaceType> cf;
      finePatch->getCoarseFaces(cf);
      std::vector<Patch::FaceType>::const_iterator iter;
      for (iter = cf.begin(); iter != cf.end(); ++iter) {
        Patch::FaceType patchFace = *iter;

        bool isRight_CP_FP_pair = false;
        CellIterator f_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
        fineLevel_CFI_Iterator(
          patchFace, coarsePatch, finePatch, f_iter, isRight_CP_FP_pair);

        if (isRight_CP_FP_pair) {

          int n_CFI_cells = 0;
          int count       = getFaceMark(1, finePatch, patchFace);
          for (; !f_iter.done(); f_iter++) {
            n_CFI_cells += 1;
          }

          // divide the number of cells
          IntVector rr  = finePatch->getLevel()->getRefinementRatio();
          IntVector dir = finePatch->getFaceAxes(patchFace);
          int y         = dir[1];
          int z         = dir[2];
          count += n_CFI_cells / (rr[y] * rr[z]);

          setFaceMark(1, finePatch, patchFace, count);
        } // right cp_fp_pair
      }   // face loop
    }     // finePatch loop
  }       // coarsepatch loop
}

/*___________________________________________________________________
 Function~  AMRICE::reflux_BP_check_CFI_cells--  (bulletproofing)
 Purpose~   Check if each coarse fine interface cell was "touched"
            during refluxing
_____________________________________________________________________*/
void
AMRICE::reflux_BP_check_CFI_cells(const ProcessorGroup*,
                                  const PatchSubset* coarsePatches,
                                  const MaterialSubset*,
                                  DataWarehouse*,
                                  DataWarehouse* new_dw,
                                  std::string description)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  for (int c_p = 0; c_p < coarsePatches->size(); c_p++) {
    const Patch* coarsePatch = coarsePatches->get(c_p);

    Level::selectType finePatches;
    coarsePatch->getOtherLevelPatches(1, finePatches, 1);

    cout_doing << d_myworld->myRank()
               << " Doing reflux_BP_check_CFI_cells \t\t\t AMRICE L-"
               << fineLevel->getIndex() << "->" << coarseLevel->getIndex()
               << endl;

    for (size_t p = 0; p < finePatches.size(); p++) {
      const Patch* finePatch = finePatches[p];

      if (finePatch->hasCoarseFaces()) {

        std::vector<Patch::FaceType> cf;
        finePatch->getCoarseFaces(cf);
        std::vector<Patch::FaceType>::const_iterator iter;
        for (iter = cf.begin(); iter != cf.end(); ++iter) {
          Patch::FaceType patchFace = *iter;

          // This makes sure that the processor that "touched" the cell is also
          //  going to check it.  Each processor can have a different instance
          //  of a patch.
          IntVector dummy;
          bool isRight_CP_FP_pair;
          CellIterator dummy_iter(IntVector(0, 0, 0), IntVector(0, 0, 0));
          refluxCoarseLevelIterator(patchFace,
                                    coarsePatch,
                                    finePatch,
                                    fineLevel,
                                    dummy_iter,
                                    dummy,
                                    isRight_CP_FP_pair,
                                    description);

          if (isRight_CP_FP_pair) {

            bool isRight_CP_FP_pair = false;
            CellIterator f_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
            fineLevel_CFI_Iterator(
              patchFace, coarsePatch, finePatch, f_iter, isRight_CP_FP_pair);

            if (isRight_CP_FP_pair) {

              int n_ice_matls = d_materialManager->getNumMaterials("ICE");
              int n_touched_cells =
                (getFaceMark(0, finePatch, patchFace)) / n_ice_matls;
              int n_CFI_cells = getFaceMark(1, finePatch, patchFace);
              //__________________________________
              // If the number of "marked" cells/numICEMatls != n_CFI_cells
              // ignore if a recompute time step has already been requested
              bool rts = new_dw->recomputeTimestep();

              if (n_touched_cells != n_CFI_cells && !rts) {
                std::ostringstream warn;
                warn << d_myworld->myRank() << " AMRICE:refluxing_"
                     << description
                     << " \n CFI face: " << finePatch->getFaceName(patchFace)
                     << " cells were 'touched' " << n_touched_cells << " times"
                     << " it should have been 'touched' " << n_CFI_cells
                     << " times "
                     << "\n patch " << *finePatch << "\n finePatchLevel "
                     << finePatch->getLevel()->getIndex() << endl;
                // cout << warn.str() << endl;
                throw InternalError(warn.str(), __FILE__, __LINE__);
              }
            }
          }
        } // face iter
      }   // has CFI
      clearFaceMarks(0, finePatch);
      clearFaceMarks(1, finePatch);
    } // //finePatches

  } // coarsePatches
}

/*_____________________________________________________________________
 Function~  AMRICE::scheduleInitialErrorEstimate--
______________________________________________________________________*/
void
AMRICE::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                     SchedulerP& sched)
{
  scheduleErrorEstimate(coarseLevel, sched);
}

/*_____________________________________________________________________
 Function~  AMRICE::scheduleErrorEstimate--
______________________________________________________________________*/
void
AMRICE::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  cout_doing << d_myworld->myRank() << " AMRICE::scheduleErrorEstimate \t\t\tL-"
             << coarseLevel->getIndex() << '\n';

  // This method is called at both initialization and otherwise. At
  // initialization the old DW, i.e. DW(0) will not exist. As such,
  // get the time step from the new DW, i.e. DW(1).  Otherwise for a
  // normal time step get the time step from the new DW.
  timeStep_vartype timeStepVar(0);
  if (sched->get_dw(0) && sched->get_dw(0)->exists(getTimestepLabel())) {
    sched->get_dw(0)->get(timeStepVar, getTimestepLabel());
  } else if (sched->get_dw(1) && sched->get_dw(1)->exists(getTimestepLabel())) {
    sched->get_dw(1)->get(timeStepVar, getTimestepLabel());
  }

  bool initial = (timeStepVar == 0); // during initialization

  Task* t =
    scinew Task("AMRICE::errorEstimate", this, &AMRICE::errorEstimate, initial);

  Ghost::GhostType gac          = Ghost::AroundCells;
  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  Task::SearchTG OldTG =
    Task::SearchTG::OldTG; // search the oldTG for the computes

  const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");
  const MaterialSet* all_matls = d_materialManager->allMaterials();

  const MaterialSubset* matls_sub;
  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
  const MaterialSubset* all_matls_sub = all_matls->getUnion();

  // Only require ice_matls during initialization, we don't have *_CC
  // variables for mpm_matls at this point
  if (initial) {
    matls_sub = ice_matls_sub;
  } else {
    matls_sub = all_matls_sub;
  }

  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, matls_sub, gac, 1, OldTG);
  t->needs(
    Task::NewDW, d_ice_labels->temperature_CCLabel, matls_sub, gac, 1, OldTG);
  t->needs(
    Task::NewDW, d_ice_labels->velocity_CCLabel, matls_sub, gac, 1, OldTG);
  t->needs(
    Task::NewDW, d_ice_labels->vol_frac_CCLabel, matls_sub, gac, 1, OldTG);
  t->needs(Task::NewDW,
              d_ice_labels->press_CCLabel,
              d_press_matl,
              oims,
              gac,
              1,
              OldTG);

  t->computes(d_ice_labels->mag_grad_rho_CCLabel);
  t->computes(d_ice_labels->mag_grad_temperature_CCLabel);
  t->computes(d_ice_labels->mag_div_velocity_CCLabel);
  t->computes(d_ice_labels->mag_grad_vol_frac_CCLabel);
  t->computes(d_ice_labels->mag_grad_press_CCLabel, d_press_matl);

  t->modifies(d_regridder->getRefineFlagLabel(),
              d_regridder->refineFlagMaterials(),
              oims);
  t->modifies(d_regridder->getRefinePatchFlagLabel(),
              d_regridder->refineFlagMaterials(),
              oims);

  sched->addTask(
    t, coarseLevel->eachPatch(), d_materialManager->allMaterials());

  //__________________________________
  // Models
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model) {
      fb_model->scheduleErrorEstimate(coarseLevel, sched);
    }

    HEChemModel* hec_model = dynamic_cast<HEChemModel*>(model.get());
    if (hec_model) {
      hec_model->scheduleErrorEstimate(coarseLevel, sched);
    }
  }
}

/*_____________________________________________________________________
 Function~  AMRICE::set_refinementFlags
______________________________________________________________________*/
void
AMRICE::set_refineFlags(constCCVariable<double>& mag_grad_q_CC,
                        double threshold,
                        CCVariable<int>& refineFlag,
                        PerPatch<PatchFlagP>& refinePatchFlag,
                        const Patch* patch)
{
  PatchFlag* refinePatch = refinePatchFlag.get().get_rep();
  for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
    IntVector c = *iter;

    if (mag_grad_q_CC[c] > threshold) {
      refineFlag[c] = true;
      refinePatch->set();
    }
  }
}
/*_____________________________________________________________________
 Function~  AMRICE::errorEstimate--
______________________________________________________________________*/
void
AMRICE::errorEstimate(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* /*matls*/,
                      DataWarehouse*,
                      DataWarehouse* new_dw,
                      bool initial)
{
  const Level* level = getLevel(patches);
  cout_doing << d_myworld->myRank()
             << " Doing errorEstimate \t\t\t\t\t AMRICE L-"
             << level->getIndex();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    cout_doing << " patch " << patch->getID() << endl;
    Ghost::GhostType gac             = Ghost::AroundCells;
    const VarLabel* refineFlagLabel  = d_regridder->getRefineFlagLabel();
    const VarLabel* refinePatchLabel = d_regridder->getRefinePatchFlagLabel();

    CCVariable<int> refineFlag;
    new_dw->getModifiable(refineFlag, refineFlagLabel, 0, patch);

    PerPatch<PatchFlagP> refinePatchFlag;
    new_dw->get(refinePatchFlag, refinePatchLabel, 0, patch);

    //__________________________________
    //  PRESSURE
    constCCVariable<double> press_CC;
    CCVariable<double> mag_grad_press_CC;

    new_dw->get(press_CC, d_ice_labels->press_CCLabel, 0, patch, gac, 1);
    new_dw->allocateAndPut(
      mag_grad_press_CC, d_ice_labels->mag_grad_press_CCLabel, 0, patch);
    mag_grad_press_CC.initialize(0.0);

    compute_Mag_gradient(press_CC, mag_grad_press_CC, patch);

    //__________________________________
    //  initialize mag_grad for all matls
    int numAllMatls = d_materialManager->getNumMaterials();
    std::vector<CCVariable<double>> mag_grad_rho_CC(numAllMatls);
    std::vector<CCVariable<double>> mag_grad_temp_CC(numAllMatls);
    std::vector<CCVariable<double>> mag_grad_vol_frac_CC(numAllMatls);
    std::vector<CCVariable<double>> mag_div_vel_CC(numAllMatls);

    for (int m = 0; m < numAllMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      new_dw->allocateAndPut(
        mag_grad_rho_CC[indx], d_ice_labels->mag_grad_rho_CCLabel, indx, patch);
      new_dw->allocateAndPut(mag_grad_temp_CC[indx],
                             d_ice_labels->mag_grad_temperature_CCLabel,
                             indx,
                             patch);
      new_dw->allocateAndPut(mag_div_vel_CC[indx],
                             d_ice_labels->mag_div_velocity_CCLabel,
                             indx,
                             patch);
      new_dw->allocateAndPut(mag_grad_vol_frac_CC[indx],
                             d_ice_labels->mag_grad_vol_frac_CCLabel,
                             indx,
                             patch);

      mag_grad_rho_CC[indx].initialize(0.0);
      mag_grad_temp_CC[indx].initialize(0.0);
      mag_div_vel_CC[indx].initialize(0.0);
      mag_grad_vol_frac_CC[indx].initialize(0.0);
    } // matls

    //__________________________________
    //  RHO, TEMP, VEL_CC, VOL_FRAC
    // During initialization only compute for ICE matls
    int numMatls = 0;
    if (initial) {
      numMatls = d_materialManager->getNumMaterials("ICE");
    } else {
      numMatls = d_materialManager->getNumMaterials();
    }

    for (int m = 0; m < numMatls; m++) {

      Material* matl;
      if (initial) {
        matl = (ICEMaterial*)d_materialManager->getMaterial("ICE", m);
      } else {
        matl = d_materialManager->getMaterial(m);
      }

      int indx = matl->getDWIndex();
      constCCVariable<double> rho_CC, temp_CC, vol_frac_CC;
      constCCVariable<Vector> vel_CC;

      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, gac, 1);
      new_dw->get(
        temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, gac, 1);
      new_dw->get(vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, gac, 1);
      new_dw->get(
        vol_frac_CC, d_ice_labels->vol_frac_CCLabel, indx, patch, gac, 1);

      //__________________________________
      // compute the magnitude of the gradient/divergence
      compute_Mag_gradient(rho_CC, mag_grad_rho_CC[indx], patch);

      compute_Mag_gradient(temp_CC, mag_grad_temp_CC[indx], patch);

      compute_Mag_gradient(vol_frac_CC, mag_grad_vol_frac_CC[indx], patch);

      compute_Mag_Divergence(vel_CC, mag_div_vel_CC[indx], patch);
    } // matls

    //__________________________________
    // Only set the refinement flags for certain materials
    for (int i = 0; i < (int)d_thresholdVars.size(); i++) {
      thresholdVar data         = d_thresholdVars[i];
      std::string name          = data.name;
      int matl                  = data.matl;
      double thresholdValue     = data.value;
      VarLabel* mag_grad_qLabel = VarLabel::find("mag_grad_" + name);

      if (mag_grad_qLabel == nullptr) { // bulletproofing
        throw InternalError("AMRICE::errorEstimate: label(mag_grad_" + name +
                              ") not found.",
                            __FILE__,
                            __LINE__);
      }
      constCCVariable<double> mag_grad_q_CC;
      new_dw->get(mag_grad_q_CC, mag_grad_qLabel, matl, patch, Ghost::None, 0);

      set_refineFlags(
        mag_grad_q_CC, thresholdValue, refineFlag, refinePatchFlag, patch);
    }
  } // patches
}
