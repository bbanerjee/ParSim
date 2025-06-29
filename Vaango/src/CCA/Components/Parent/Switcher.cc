/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Parent/Switcher.h>

#include <CCA/Components/Parent/ComponentFactory.h>
#include <CCA/Components/ProblemSpecification/ProblemSpecReader.h>
#include <CCA/Components/Solvers/SolverFactory.h>
#include <CCA/Components/SwitchingCriteria/None.h>
#include <CCA/Components/SwitchingCriteria/SwitchingCriteriaFactory.h>

#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>
#include <CCA/Ports/SolverInterface.h>
#include <CCA/Ports/SwitchingCriteria.h>

#include <Core/Containers/ConsecutiveRangeSet.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>

#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/GridP.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <Core/OS/Dir.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Util/StringUtil.h>

#include <sci_defs/uintah_defs.h>

using namespace Uintah;

// export SCI_DEBUG="SWITCHER:+"
DebugStream Switcher::switcher_dbg("SWITCHER",
                                   "Switcher",
                                   "Switcher debug stream",
                                   false);

#define ALL_LEVELS 99

// ToDo:
// - test carry over and init vars
// - test in parallel
// - test restarting capability
// - fix so each subcomponent filebase name is used for uda.
// - test different components (mpmice, impm) for subcomponents

//__________________________________
// In the constructor read the master ups file
// For each subcomponent in the ups file:
//     -
Switcher::Switcher(const ProcessorGroup* myworld,
                   const MaterialManagerP& mat_manager,
                   ProblemSpecP& master_ups,
                   [[maybe_unused]] const std::string& uda)
  : SimulationCommon(myworld, mat_manager)
  , d_master_ups(master_ups)
{
  proc0cout << "-----------------------------Switcher::Switcher top"
            << std::endl;

  d_switch_label =
    VarLabel::create("switchFlag", max_vartype::getTypeDescription());

  std::set<std::string> simComponents;

  ProblemSpecP sim_block = master_ups->findBlock("SimulationComponent");

  //  loop over the subcomponents
  ProblemSpecP child = sim_block->findBlock("subcomponent");
  for (; child != nullptr; child = child->findNextBlock("subcomponent")) {

    //  Read in subcomponent ups file and store the filename
    string input_file("");
    if (!child->get("input_file", input_file)) {
      throw ProblemSetupException(
        "Need 'input_file' for subcomponent", __FILE__, __LINE__);
    }

    proc0cout << "Input file:\t\t" << input_file << std::endl;

    d_in_file.push_back(input_file);
    ProblemSpecP subCompUps = ProblemSpecReader().readInputFile(input_file);

    // get the component name from the input file, and the uda arg is not needed
    // for normal simulations...
    std::string sim_comp;
    ProblemSpecP sim_ps = subCompUps->findBlock("SimulationComponent");
    sim_ps->getAttribute("type", sim_comp);
    simComponents.insert(sim_comp);

    // create simulation port and attach it switcher component
    std::unique_ptr<UintahParallelComponent> comp =
      ComponentFactory::create(subCompUps, myworld, d_materialManager, "");

    SimulationInterface* sim = dynamic_cast<SimulationInterface*>(comp.get());
    attachPort("simulator", sim);

    // create solver  port and attach it to the switcher component
    std::shared_ptr<SolverInterface> solver =
      SolverFactory::create(subCompUps, myworld);
    attachPort("sub_solver", solver.get());
    comp->attachPort("solver", solver.get());

    // create switching criteria port and attach it switcher component
    std::shared_ptr<SwitchingCriteria> switch_criteria =
      SwitchingCriteriaFactory::create(child, myworld);
    if (switch_criteria) {
      switch_criteria->setSwitchLabel(d_switch_label);
      attachPort("switch_criteria", switch_criteria.get());
      comp->attachPort("switch_criteria", switch_criteria.get());
    }

    // Get the variables that will need to be initialized by this subcomponent
    std::unique_ptr<initVars> initVar = std::make_unique<initVars>();
    for (ProblemSpecP var = child->findBlock("init"); var != nullptr;
         var              = var->findNextBlock("init")) {

      std::map<std::string, string> attributes;
      var->getAttributes(attributes);

      // matlsetNames
      std::string matls = attributes["matls"];
      initVar->matlSetNames.push_back(matls);

      // levels
      std::stringstream s_level(attributes["levels"]);
      int levels = ALL_LEVELS;
      s_level >> levels;
      initVar->levels.push_back(levels);

      // variable name
      std::string varName = attributes["var"];
      initVar->varNames.push_back(varName);
    }

    d_initVars.at(d_numComponents) = std::move(initVar);
    d_numComponents++;

    proc0cout << "\n";
  } // loop over subcomponents

  //__________________________________
  // Bulletproofing:
  if (simComponents.count("mpm") && simComponents.count("rmpmice")) {
    throw ProblemSetupException("Switcher: The simulation subComponents "
                                "rmpmice and mpm cannot be used together",
                                __FILE__,
                                __LINE__);
  }

  //__________________________________
  // Bulletproofing:
  // Make sure that a switching criteria was specified.  For n subcomponents,
  // there should be n-1 switching critiera specified.
  int num_switch_criteria = 0;
  for (size_t i = 0; i < d_numComponents; i++) {
    UintahParallelComponent* comp =
      dynamic_cast<UintahParallelComponent*>(getPort("sim", i));
    SwitchingCriteria* sw =
      dynamic_cast<SwitchingCriteria*>(comp->getPort("switch_criteria"));
    if (sw) {
      num_switch_criteria++;
    }
  }

  if (static_cast<size_t>(num_switch_criteria) != d_numComponents - 1) {
    throw ProblemSetupException("Do not have enough switching criteria "
                                "specified for the number of components.",
                                __FILE__,
                                __LINE__);
  }

  //__________________________________
  // Add the "None" SwitchCriteria to the last component, so the switchFlag
  // label is computed in the last stage.
  UintahParallelComponent* last_comp = dynamic_cast<UintahParallelComponent*>(
    getPort("simulator", d_numComponents - 1));

  SwitchingCriteria* none_switch_criteria = scinew None();

  // Attaching to switcher so that the switcher can delete it
  attachPort("switch_criteria", none_switch_criteria);
  last_comp->attachPort("switch_criteria", none_switch_criteria);

  //__________________________________
  // Get the vars that will need to be carried over
  for (ProblemSpecP var = sim_block->findBlock("carry_over"); var != 0;
       var              = var->findNextBlock("carry_over")) {
    std::map<std::string, string> attributes;
    var->getAttributes(attributes);
    string name  = attributes["var"];
    string matls = attributes["matls"];
    string level = attributes["level"];

    if (name != "") {
      d_carryOverVars.push_back(name);
    }

    MaterialSubset* carry_over_matls = 0;
    if (matls != "") {
      carry_over_matls                   = scinew MaterialSubset;
      ConsecutiveRangeSet crs            = matls;
      ConsecutiveRangeSet::iterator iter = crs.begin();

      for (; iter != crs.end(); iter++) {
        carry_over_matls->add(*iter);
      }
      carry_over_matls->addReference();
    }

    d_carryOverVarMatls.push_back(carry_over_matls);
    if (level == "finest") {
      d_carryOverFinestLevelOnly.push_back(true);
    } else {
      d_carryOverFinestLevelOnly.push_back(false);
    }
  } // loop over

  d_computedVars.clear();

  proc0cout << "Number of components " << d_numComponents << endl;
  proc0cout << "-----------------------------Switcher::Switcher bottom"
            << std::endl;
}
//______________________________________________________________________
//
Switcher::~Switcher()
{

  switcher_dbg << d_myworld->myRank() << " Switcher::~Switcher" << endl;

  for (unsigned i = 0; i < d_carryOverVarMatls.size(); i++) {
    if (d_carryOverVarMatls[i] && d_carryOverVarMatls[i]->removeReference()) {
      delete d_carryOverVarMatls[i];
    }
  }
  d_carryOverVarMatls.clear();
}

//______________________________________________________________________
// Setup the first component
void
Switcher::problemSetup([[maybe_unused]] const ProblemSpecP& params,
                       const ProblemSpecP& restart_prob_spec,
                       GridP& grid,
                       [[maybe_unused]] const std::string& input_ups_dir)
{
  switcher_dbg << "Doing ProblemSetup \t\t\t\tSwitcher" << std::endl;
  if (restart_prob_spec) {
    readSwitcherState(restart_prob_spec, d_materialManager);
  }

  switchSimulator(restart_prob_spec, grid);

  // init Variables:
  //   - determine the label from the string names
  //   - determine the MaterialSet from the string matlSetName
  //   - store this info to be used later
  for (auto& [comp, initVar_p] : d_initVars) {
    proc0cout << " init Variables:  component: " << comp << std::endl;

    // Find the varLabel
    std::vector<std::string>& varNames = initVar_p->varNames;
    std::vector<VarLabel*> varLabels   = initVar_p->varLabels;

    for (auto& varName : varNames) {
      VarLabel* label = VarLabel::find(varName);

      if (!label) {
        std::string error =
          "ERROR: Switcher: Cannot find init VarLabel" + varName;
        throw ProblemSetupException(error, __FILE__, __LINE__);
      }

      varLabels.push_back(label);

      // so the variable is not scrubbed from the data warehouse
      d_scheduler->overrideVariableBehavior(
        varName, false, false, true, false, false);
    }

    d_initVars[comp]->varLabels = varLabels;
  }

  // Carry over labels
  for (auto& varName : d_carryOverVars) {
    VarLabel* label = VarLabel::find(varName);
    if (label) {
      d_carryOverVarLabels.push_back(label);

      // so variable is not scrubbed from the data warehouse
      d_scheduler->overrideVariableBehavior(
        varName, false, false, true, false, false);
    } else {
      std::string error =
        "ERROR: Switcher: Cannot find carry_over VarLabel" + varName;
      throw ProblemSetupException(error, __FILE__, __LINE__);
    }
  }
}

void
Switcher::readSwitcherState(const ProblemSpecP& spec,
                            MaterialManagerP& mat_manager)
{
  ProblemSpecP ps = static_cast<ProblemSpecP>(spec);

  int tmp;
  ps->get("switcherComponentIndex", tmp);
  d_componentIndex = tmp;

  ps->get("switcherState", tmp);
  d_switchState = static_cast<switchState>(tmp);

  int numMatls = 0;
  ps->get("switcherCarryOverMatls", numMatls);

  if (numMatls != 0) {
    MaterialSet* new_matls = scinew MaterialSet;
    new_matls->addReference();
    new_matls->createEmptySubsets(1);

    for (int i = 0; i < numMatls; i++) {
      new_matls->getSubset(0)->add(i);
    }

    mat_manager->setOriginalMatlsFromRestart(new_matls);
  }

  proc0cout << "  Switcher RESTART: component index = " << d_componentIndex
            << std::endl;
}

void
Switcher::switchSimulator(const ProblemSpecP& restart_prob_spec,
                          const GridP& grid)
{
  // Get the initial simulation component and initialize the need components
  proc0cout << "\n------------ Switching to application (" << d_componentIndex
            << ") \n";
  proc0cout << "  Reading input file: " << d_in_file[d_componentIndex] << "\n";

  // Read the ups file for the first subcomponent.
  ProblemSpecP subCompUps =
    ProblemSpecReader().readInputFile(d_in_file[d_componentIndex]);

  UintahParallelComponent* simComp = dynamic_cast<UintahParallelComponent*>(
    getPort("simulator", d_componentIndex));

  d_simulator = dynamic_cast<SimulationInterface*>(simComp);
  d_simulator->setComponents(this);

  // Send the subcomponent's UPS file to it's sim interface.
  d_simulator->problemSetup(
    subCompUps, restart_prob_spec, const_cast<GridP&>(grid));

  // Send the subcomponent's UPS file to the data archiver to get the
  // output and checkpointing parameters.
  d_output->problemSetup(subCompUps, restart_prob_spec, d_materialManager);

  // Read in the grid adaptivity flag from the subcomponent's UPS file.
  if (d_regridder) {
    d_regridder->switchInitialize(subCompUps);
  }

  // Send the subcomponent's UPS file to the switcher's simulation
  // time.  Note this goes into the switcher not the subcomponent.
  proc0cout << "  Reading the <Time> block from: "
            << Uintah::basename(subCompUps->getFile()) << "\n";

  problemSetupDeltaT(subCompUps);
}

//______________________________________________________________________
//
void
Switcher::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleInitialize");
  d_simulator->scheduleInitialize(level, sched);
}

//______________________________________________________________________
//
void
Switcher::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleRestartInitialize");
  d_restarting = true;
  d_simulator->scheduleRestartInitialize(level, sched);
}

//______________________________________________________________________
//
void
Switcher::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleComputeStableTimestep");
  d_simulator->scheduleComputeStableTimestep(level, sched);
}

//______________________________________________________________________
//
void
Switcher::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleTimeAdvance");
  d_simulator->scheduleTimeAdvance(level, sched);
}

//______________________________________________________________________
//
void
Switcher::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleFinalizeTimestep");

  d_simulator->scheduleFinalizeTimestep(level, sched);

  scheduleSwitchTest(level, sched);

  // compute variables that are required from the old_dw for the next
  // subcomponent
  scheduleInitNewVars(level, sched);

  scheduleSwitchInitialization(level, sched);

  // carry over vars that will be needed by a future component
  scheduleCarryOverVars(level, sched);
}

void
Switcher::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleSwitchTest");

  d_simulator->scheduleSwitchTest(level, sched); // generates switch test data;

  Task* t = scinew Task("Switcher::switchTest", this, &Switcher::switchTest);

  t->setType(Task::OncePerProc);

  // the component is responsible for determining when it is to switch.
  t->needs(Task::NewDW, d_switch_label);
  sched->addTask(t,
                 d_loadBalancer->getPerProcessorPatchSet(level),
                 d_materialManager->allMaterials());
}

//  Set the flag if switch criteria has been satisfied.
void
Switcher::switchTest(const ProcessorGroup*,
                     [[maybe_unused]] const PatchSubset* patches,
                     [[maybe_unused]] const MaterialSubset* matls,
                     [[maybe_unused]] DataWarehouse* old_dw,
                     DataWarehouse* new_dw)
{
  max_vartype switch_condition;
  new_dw->get(switch_condition, d_switch_label, 0);

  if (switch_condition) {
    // actually PERFORM the switch during the next needRecompile; set back to
    // idle then
    d_switchState = switchState::switching;
  } else {
    d_switchState = switchState::idle;
  }
}

void
Switcher::scheduleInitNewVars(const LevelP& level, SchedulerP& sched)
{
  unsigned int nextComp_indx = d_componentIndex + 1;

  if (nextComp_indx >= d_numComponents) {
    return;
  }

  printSchedule(level, switcher_dbg, "Switcher::scheduleInitNewVars");

  Task* t = scinew Task("Switcher::initNewVars", this, &Switcher::initNewVars);

  initVars* initVar = d_initVars.find(nextComp_indx)->second.get();

  std::vector<const MaterialSet*> matlSet;

  for (unsigned i = 0; i < initVar->varLabels.size(); i++) {

    VarLabel* label = initVar->varLabels[i];

    // Find the MaterialSet for this variable
    // and put that set in the global structure
    const MaterialSet* matls;

    std::string nextComp_matls = initVar->matlSetNames[i];
    if (nextComp_matls == "ice_matls") {
      matls = d_materialManager->allMaterials("ICE");
    } else if (nextComp_matls == "mpm_matls") {
      matls = d_materialManager->allMaterials("MPM");
    } else if (nextComp_matls == "all_matls") {
      matls = d_materialManager->allMaterials();
    } else {
      throw ProblemSetupException("Bad material set", __FILE__, __LINE__);
    }

    matlSet.push_back(matls);
    proc0cout << "init Variable  " << initVar->varNames[i]
              << " \t matls: " << nextComp_matls << " levels "
              << initVar->levels[i] << std::endl;

    const MaterialSubset* matl_ss = matls->getUnion();

    t->computes(label, matl_ss);
  }

  d_initVars[nextComp_indx]->matls = matlSet;

  t->needs(Task::NewDW, d_switch_label);
  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials());
}

//  This only get executed if a switching components has been called for.
void
Switcher::initNewVars(const ProcessorGroup*,
                      const PatchSubset* patches,
                      [[maybe_unused]] const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw)
{
  max_vartype switch_condition;
  new_dw->get(switch_condition, d_switch_label, 0);

  if (!switch_condition) {
    return;
  }

  switcher_dbg << "__________________________________" << std::endl
               << "initNewVars \t\t\t\tSwitcher" << std::endl;

  // loop over the init vars, initialize them and put them in the new_dw
  initVars* initVar = d_initVars.find(d_componentIndex + 1)->second.get();

  for (unsigned i = 0; i < initVar->varLabels.size(); i++) {

    VarLabel* l                 = initVar->varLabels[i];
    const MaterialSubset* matls = initVar->matls[i]->getUnion();

    //__________________________________
    // initialize a variable on this level?
    const Level* level = getLevel(patches);
    int numLevels      = level->getGrid()->numLevels();
    int L_indx         = getLevel(patches)->getIndex();
    int relative_indx  = L_indx - numLevels;
    int init_Levels    = initVar->levels[i];

    switcher_dbg << "    varName: " << l->getName() << " \t\t matls "
                 << initVar->matlSetNames[i] << " level " << init_Levels
                 << std::endl;

    bool onThisLevel = false;

    if (init_Levels == L_indx ||        // user can specify: a level,
        init_Levels == ALL_LEVELS ||    // nothing,
        init_Levels == relative_indx) { // or a relative indx, -1, -2
      onThisLevel = true;
    }

    if (onThisLevel == false) {
      continue;
    }

    // Bulletproofing
    if (l->typeDescription()->getType() ==
          TypeDescription::Type::ParticleVariable &&
        relative_indx != -1) {
      std::ostringstream warn;
      warn << " \nERROR: switcher: subcomponent: init var: (" << l->getName()
           << ") \n particle variables can only be initialized on the finest "
              "level \n"
           << " of a multilevel grid.  Add levels=\"-1\" to that variable"
           << std::endl;
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    // initialization section
    for (int m = 0; m < matls->size(); m++) {
      const int indx = matls->get(m);

      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);

        switcher_dbg << "    indx: " << indx << " patch " << *patch << " "
                     << l->getName() << std::endl;

        switch (l->typeDescription()->getType()) {

          //__________________________________
          //
          case TypeDescription::Type::CCVariable:
            switch (l->typeDescription()->getSubType()->getType()) {
              case TypeDescription::Type::double_type: {
                CCVariable<double> q;
                new_dw->allocateAndPut(q, l, indx, patch);
                q.initialize(0);
                break;
              }
              case TypeDescription::Type::Vector: {
                CCVariable<Vector> q;
                new_dw->allocateAndPut(q, l, indx, patch);
                q.initialize(Vector(0, 0, 0));
                break;
              }
              default:
                throw InternalError(
                  "ERROR:Switcher::initNewVars Unknown CCVariable type",
                  __FILE__,
                  __LINE__);
            }
            break;
          //__________________________________
          //
          case TypeDescription::Type::NCVariable:
            switch (l->typeDescription()->getSubType()->getType()) {
              case TypeDescription::Type::double_type: {
                NCVariable<double> q;
                new_dw->allocateAndPut(q, l, indx, patch);
                q.initialize(0);
                break;
              }
              case TypeDescription::Type::Vector: {
                NCVariable<Vector> q;
                new_dw->allocateAndPut(q, l, indx, patch);
                q.initialize(Vector(0, 0, 0));
                break;
              }
              default:
                throw InternalError(
                  "ERROR:Switcher::initNewVars Unknown NCVariable type",
                  __FILE__,
                  __LINE__);
            }
            break;
          //__________________________________
          //
          case TypeDescription::Type::ParticleVariable: {

            ParticleSubset* pset = old_dw->getParticleSubset(indx, patch);
            switch (l->typeDescription()->getSubType()->getType()) {
              case TypeDescription::Type::int_type: {
                ParticleVariable<int> q;
                new_dw->allocateAndPut(q, l, pset);

                for (ParticleSubset::iterator iter = pset->begin();
                     iter != pset->end();
                     iter++) {
                  q[*iter] = 0;
                }

                break;
              }
              case TypeDescription::Type::double_type: {
                ParticleVariable<double> q;
                new_dw->allocateAndPut(q, l, pset);

                for (ParticleSubset::iterator iter = pset->begin();
                     iter != pset->end();
                     iter++) {
                  q[*iter] = 0;
                }
                break;
              }
              case TypeDescription::Type::Vector: {
                ParticleVariable<Vector> q;
                new_dw->allocateAndPut(q, l, pset);

                for (ParticleSubset::iterator iter = pset->begin();
                     iter != pset->end();
                     iter++) {
                  q[*iter] = Vector(0, 0, 0);
                }
                break;
              }
              case TypeDescription::Type::Matrix3: {
                ParticleVariable<Matrix3> q;
                new_dw->allocateAndPut(q, l, pset);
                for (ParticleSubset::iterator iter = pset->begin();
                     iter != pset->end();
                     iter++) {
                  q[*iter].Identity();
                }
                break;
              }
              default:
                throw InternalError(
                  "ERROR:Switcher::initNewVars Unknown particle type",
                  __FILE__,
                  __LINE__);
            }
            break;
          }
          default:
            throw InternalError(
              "ERROR:Switcher::initNewVars Unknown Variable type",
              __FILE__,
              __LINE__);
        }
      } // patch loop
    }   // matl loop
  }     // varlabel loop
  switcher_dbg << "__________________________________" << std::endl;
}

void
Switcher::scheduleSwitchInitialization(const LevelP& level, SchedulerP& sched)
{
  if (d_doSwitching[level->getIndex()]) {
    printSchedule(
      level, switcher_dbg, "Switcher::scheduleSwitchInitialization");
    d_simulator->switchInitialize(level, sched);
  }
}

void
Switcher::scheduleCarryOverVars(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, switcher_dbg, "Switcher::scheduleCarryOverVars");
  int L_indx = level->getIndex();

  if (d_computedVars.size() == 0) {
    // get the set of computed vars like this, because by scheduling a
    // carry-over var, we add to the compute list
    d_computedVars = sched->getComputedVars();
  }

  if (d_doSwitching[L_indx] || d_restarting) {
    // clear and reset carry-over db
    if (L_indx >= (int)d_doCarryOverVarPerLevel.size()) {
      d_doCarryOverVarPerLevel.resize(L_indx + 1);
    }
    d_doCarryOverVarPerLevel[L_indx].clear();

    // rebuild carry-over database

    // mark each var as carry over if it's not in the computed list
    for (unsigned i = 0; i < d_carryOverVarLabels.size(); i++) {

      bool do_on_this_level = !d_carryOverFinestLevelOnly[i] ||
                              L_indx == level->getGrid()->numLevels() - 1;

      bool no_computes =
        d_computedVars.find(d_carryOverVarLabels[i]) == d_computedVars.end();

      bool trueFalse = (do_on_this_level && no_computes);
      d_doCarryOverVarPerLevel[L_indx].push_back(trueFalse);
    }
  }

  Task* t =
    scinew Task("Switcher::carryOverVars", this, &Switcher::carryOverVars);

  // schedule the vars to be carried over (if this happens before a switch,
  // don't do it)
  if (L_indx < (int)d_doCarryOverVarPerLevel.size()) {

    for (unsigned int i = 0; i < d_carryOverVarLabels.size(); i++) {

      if (d_doCarryOverVarPerLevel[L_indx][i]) {

        VarLabel* var         = d_carryOverVarLabels[i];
        MaterialSubset* matls = d_carryOverVarMatls[i];

        t->needs(Task::OldDW, var, matls, Ghost::None, 0);
        t->computes(var, matls);

        if (d_myworld->myRank() == 0) {
          if (matls) {
            std::cout << d_myworld->myRank() << "  Carry over " << *var
                      << "\t\tmatls: " << *matls << " on level " << L_indx
                      << std::endl;
          } else {
            std::cout << d_myworld->myRank() << "  Carry over " << *var
                      << "\t\tAll matls on level " << L_indx << "\n";
          }
        }
      }
    }
  }
  sched->addTask(
    t, level->eachPatch(), d_materialManager->originalAllMaterials());
}

void
Switcher::carryOverVars(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* matls,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();

  if (L_indx < (int)d_doCarryOverVarPerLevel.size()) {

    for (unsigned int i = 0; i < d_carryOverVarLabels.size(); i++) {

      if (d_doCarryOverVarPerLevel[L_indx][i]) {

        const VarLabel* label = d_carryOverVarLabels[i];
        const MaterialSubset* xfer_matls =
          d_carryOverVarMatls[i] == 0 ? matls : d_carryOverVarMatls[i];

        //__________________________________
        //  reduction variables
        if (label->typeDescription()->isReductionVariable()) {

          switch (label->typeDescription()->getSubType()->getType()) {
            case Uintah::TypeDescription::Type::double_type: {
              ReductionVariable<double, Reductions::Max<double>> var_d;
              old_dw->get(var_d, label);
              new_dw->put(var_d, label);
            } break;
            default:
              throw InternalError("ERROR:Switcher::carryOverVars - Unknown "
                                  "reduction variable type",
                                  __FILE__,
                                  __LINE__);
          }
        } else { // all grid variables
          new_dw->transferFrom(old_dw, label, patches, xfer_matls);
        }
      }
    }
  }
}

//______________________________________________________________________
//  This is where the actual component switching takes place.
bool
Switcher::needRecompile(const GridP& grid)
{
  switcher_dbg << "  Doing Switcher::needRecompile " << std::endl;

  d_restarting = true;
  d_doSwitching.resize(grid->numLevels());

  for (int i = 0; i < grid->numLevels(); i++) {
    d_doSwitching[i] = (d_switchState == switchState::switching);
  }

  if (d_switchState == switchState::switching) {
    d_switchState = switchState::idle;
    d_computedVars.clear();
    d_componentIndex++;

    d_simulator->setupForSwitching();
    d_materialManager->clearMaterials();

    // Reseting the GeometryPieceFactory only (I believe) will ever need to be
    // done by the Switcher component...
    GeometryPieceFactory::resetFactory();

    switchSimulator(nullptr, grid);

    // Each application has their own maximum initial delta T
    // specified.  On a switch from one application to the next, delT
    // needs to be adjusted to the value specified in the input file.
    proc0cout << "Switching the next delT from " << d_delT << " to "
              << d_delTInitialMax << std::endl;

    setDelT(d_delTInitialMax);

    d_materialManager->finalizeMaterials();

    proc0cout << "__________________________________\n\n";

    d_output->setSwitchState(true);

    return true;
  } else {
    d_output->setSwitchState(false);
    return false;
  }
}

void
Switcher::outputProblemSpec(ProblemSpecP& ps)
{
  ps->appendElement("switcherComponentIndex", (int)d_componentIndex);
  ps->appendElement("switcherState", (int)d_switchState);
  ps->appendElement(
    "switcherCarryOverMatls",
    d_materialManager->originalAllMaterials()->getUnion()->size());
  d_simulator->outputProblemSpec(ps);
}

double
Switcher::recomputeDelT(double dt)
{
  return d_simulator->recomputeDelT(dt);
}

//______________________________________________________________________
//     AMR
void
Switcher::scheduleRefineInterface(const LevelP& fineLevel,
                                  SchedulerP& sched,
                                  bool needCoarseOld,
                                  bool needCoarseNew)
{
  d_simulator->scheduleRefineInterface(
    fineLevel, sched, needCoarseOld, needCoarseNew);
}
//______________________________________________________________________
//
void
Switcher::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  d_simulator->scheduleRefine(patches, sched);
}
//______________________________________________________________________
//
void
Switcher::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
  d_simulator->scheduleCoarsen(coarseLevel, sched);
}

//______________________________________________________________________
//
void
Switcher::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                       SchedulerP& sched)
{
  d_simulator->scheduleInitialErrorEstimate(coarseLevel, sched);
}
//______________________________________________________________________
//
void
Switcher::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  d_simulator->scheduleErrorEstimate(coarseLevel, sched);
}
