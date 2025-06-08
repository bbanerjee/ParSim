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

#include <CCA/Components/ICE/ICE.h>

#include <CCA/Components/ICE/impAMRICE.h>

#include <CCA/Components/ICE/Advection/AdvectionFactory.h>
#include <CCA/Components/ICE/TurbulenceModel/TurbulenceFactory.h>
#include <CCA/Components/ICE/WallShearStressModel/WallShearStressFactory.h>
#include <CCA/Components/Models/MultiMatlExchange/ExchangeFactory.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Components/ICE/Core/ConservationTest.h>
#include <CCA/Components/ICE/Core/Diffusion.h>
#include <CCA/Components/ICE/CustomBCs/C_BC_driver.h>
#include <CCA/Components/ICE/EOS/EquationOfState.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/ICE/SpecificHeatModel/SpecificHeat.h>
#include <CCA/Components/ICE/TurbulenceModel/Turbulence.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/Models/FluidsBased/FluidsBasedModel.h>
#include <CCA/Components/Models/HEChem/HEChemModel.h>
#include <CCA/Components/Models/ModelFactory.h>
#include <CCA/Components/Models/ParticleBased/TracerParticles.h>
#include <CCA/Components/Models/SolidReactionModel/SolidReactionModel.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <Core/Grid/BoundaryConditions/BCUtils.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <Core/Math/Expon.h>
#include <Core/Math/FastMatrix.h>

#include <Core/Parallel/MasterLock.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>

#include <Core/Util/DebugStream.h>

#include <sci_defs/hypre_defs.h>

#include <Eigen/Dense>

#include <cfloat>
#include <iostream>
#include <sstream>
#include <vector>

#ifdef HAVE_HYPRE
#include <CCA/Components/Solvers/HypreSolver.h>
#endif

#define SET_CFI_BC 0

using namespace Uintah;
using namespace Uintah::ExchangeModels;

extern Uintah::MasterLock cerrLock;

//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "ICE_NORMAL_COUT:+,ICE_DOING_COUT:+"
//  ICE_NORMAL_COUT:  dumps out during problemSetup
//  ICE_DOING_COUT:   dumps when tasks are scheduled and performed
//  default is OFF
static DebugStream cout_norm("ICE_NORMAL_COUT", false);
static DebugStream ds_EqPress("DBG_EqPress", false);

Dout m_ice_tasks("ICE_tasks",
                 "ICE",
                 "Print task scheduling & execution",
                 false);

ICE::ICE(const ProcessorGroup* myworld, const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
{
  d_ice_labels = std::make_unique<ICELabel>();

#ifdef HAVE_HYPRE
  hypre_solver_label =
    VarLabel::create("hypre_solver_label",
                     SoleVariable<hypre_solver_structP>::getTypeDescription());
#endif

  d_conservationTest        = std::make_unique<conservationTest_flags>();
  d_conservationTest->onOff = false;

  d_customInitialize_basket = std::make_unique<customInitialize_basket>();
  d_BC_globalVars = std::make_unique<CustomBCDriver::customBC_globalVars>();
  d_BC_globalVars->lodi     = std::make_unique<Lodi_globalVars>();
  d_BC_globalVars->slip     = std::make_unique<slip_globalVars>();
  d_BC_globalVars->mms      = std::make_unique<mms_globalVars>();
  d_BC_globalVars->sine     = std::make_unique<sine_globalVars>();
  d_BC_globalVars->inletVel = std::make_unique<inletVel_globalVars>();
}

ICE::~ICE()
{
  cout_norm << d_myworld->myRank() << " Doing: ICE destructor " << std::endl;

#ifdef HAVE_HYPRE
  VarLabel::destroy(hypre_solver_label);
#endif

  if (d_press_matl && d_press_matl->removeReference()) {
    delete d_press_matl;
  }
  if (d_press_matlSet && d_press_matlSet->removeReference()) {
    delete d_press_matlSet;
  }

  //__________________________________
  // MODELS
  for (auto& model : d_models) {
    model->releaseComponents();
  }
}

double
ICE::recomputeDelT(double delT)
{
  return delT * 0.75;
}

/* _____________________________________________________________________
   Function~  ICE::problemSetup--
   _____________________________________________________________________*/
void
ICE::problemSetup(const ProblemSpecP& prob_spec,
                  const ProblemSpecP& restart_prob_spec,
                  GridP& grid,
                  [[maybe_unused]] const std::string& input_ups_dir)
{
  DOUTR(m_ice_tasks, " ICE::problemSetup ");

  d_press_matl = scinew MaterialSubset();
  d_press_matl->add(0);
  d_press_matl->addReference();

  d_press_matlSet = scinew MaterialSet();
  d_press_matlSet->add(0);
  d_press_matlSet->addReference();

  d_ref_press               = 0.0;
  ProblemSpecP phys_cons_ps = prob_spec->findBlock("PhysicalConstants");
  if (phys_cons_ps) {
    phys_cons_ps->require("reference_pressure", d_ref_press);
    phys_cons_ps->require("gravity", d_gravity);
  } else {
    throw ProblemSetupException(
      "\n Could not find the <PhysicalConstants> section in the input file.  "
      "This section contains <gravity> and <reference pressure> \n"
      " This pressure is used during the problem intialization and when\n"
      " the pressure gradient is interpolated to the MPM particles \n"
      " you must have it for all MPMICE and multimaterial ICE problems\n",
      __FILE__,
      __LINE__);
  }

  //__________________________________
  // Pull out from CFD-ICE section
  ProblemSpecP cfd_ps = prob_spec->findBlock("CFD");
  if (!cfd_ps) {
    throw ProblemSetupException(
      "\n Could not find the <CFD> section in the input file\n",
      __FILE__,
      __LINE__);
  }
  cfd_ps->require("cfl", d_CFL);

  ProblemSpecP cfd_ice_ps = cfd_ps->findBlock("ICE");
  if (!cfd_ice_ps) {
    throw ProblemSetupException(
      "\n Could not find the <CFD> <ICE> section in the input file\n",
      __FILE__,
      __LINE__);
  }
  cfd_ice_ps->get("max_iteration_equilibration", d_max_iter_equilibration);
  cfd_ice_ps->get("ClampSpecificVolume", d_clampSpecificVolume);
  cfd_ice_ps->get("applyHydrostaticPressure", d_applyHydrostaticPress);

  d_advector = AdvectionFactory::create(
    cfd_ice_ps, d_useCompatibleFluxes, d_OrderOfAdvection);

  //__________________________________
  //  Pull out add heat section
  ProblemSpecP add_heat_ps = cfd_ice_ps->findBlock("ADD_HEAT");
  if (add_heat_ps) {
    d_add_heat = true;
    add_heat_ps->require("add_heat_matls", d_add_heat_matls);
    add_heat_ps->require("add_heat_coeff", d_add_heat_coeff);
    add_heat_ps->require("add_heat_t_start", d_add_heat_t_start);
    add_heat_ps->require("add_heat_t_final", d_add_heat_t_final);
  }

  //__________________________________
  //  custom Initialization
  customInitialization_problemSetup(
    cfd_ice_ps, d_customInitialize_basket.get(), grid);

  //__________________________________
  // Pull out implicit solver parameters
  ProblemSpecP impSolver = cfd_ice_ps->findBlock("ImplicitSolver");
  if (impSolver) {
    d_delT_speedSoundKnob = 0.5; // default value when running implicit
    d_solver->readParameters(impSolver, "implicitPressure");
    d_solver->getParameters()->setSolveOnExtraCells(false);
    d_solver->getParameters()->setRecomputeTimestepOnFailure(true);
    impSolver->require("max_outer_iterations", d_max_iter_implicit);
    impSolver->require("outer_iteration_tolerance", d_outer_iter_tolerance);
    impSolver->getWithDefault(
      "iters_before_timestep_recompute", d_iters_before_timestep_recompute, 5);
    d_impICE = true;

    d_subsched = d_scheduler->createSubScheduler();
    d_subsched->initialize(3, 1);
    d_subsched->clearMappings();
    d_subsched->mapDataWarehouse(Task::ParentOldDW, 0);
    d_subsched->mapDataWarehouse(Task::ParentNewDW, 1);
    d_subsched->mapDataWarehouse(Task::OldDW, 2);
    d_subsched->mapDataWarehouse(Task::NewDW, 3);

    d_recompileSubsched = true;

    //__________________________________
    // bulletproofing
    double tol;
    ProblemSpecP p = impSolver->findBlock("Parameters");
    if (!p) {
      throw ProblemSetupException(
        "\n Could not find the <Parameters> section in the input file\n",
        __FILE__,
        __LINE__);
    }
    p->get("tolerance", tol);
    if (tol >= d_outer_iter_tolerance) {
      std::ostringstream msg;
      msg << "\n ERROR: implicit pressure: The <outer_iteration_tolerance>"
          << " must be greater than the solver tolerance <tolerance> \n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }

    if (isAMR() && d_solver->getName() != "hypreamr") {
      std::ostringstream msg;
      msg << "\n ERROR: " << d_solver->getName()
          << " cannot be used with an AMR grid \n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  }

  //__________________________________
  // Pull out TimestepControl data
  ProblemSpecP tsc_ps = cfd_ice_ps->findBlock("TimestepControl");
  if (tsc_ps) {
    tsc_ps->require("Scheme_for_delT_calc", d_delT_scheme);
    tsc_ps->require("knob_for_speedSound", d_delT_speedSoundKnob);
    tsc_ps->get("knob_for_diffusion", d_delT_diffusionKnob);

    if (d_delT_scheme != "conservative" && d_delT_scheme != "aggressive") {
      std::string warn = "ERROR:\n Scheme_for_delT_calc:  must specify either "
                         "aggressive or conservative";
      throw ProblemSetupException(warn, __FILE__, __LINE__);
    }

    if (d_delT_speedSoundKnob < 0.0 || d_delT_speedSoundKnob > 1.0) {
      std::string warn =
        "ERROR:\n knob_for_speedSound:  must be between 0 and 1";
      throw ProblemSetupException(warn, __FILE__, __LINE__);
    }
  }

  //  User defined pressure gradient
  ProblemSpecP fpg_ps = cfd_ice_ps->findBlock("fixedPressureGradient");
  if (fpg_ps) {
    double x, y, z;
    fpg_ps->getWithDefault("x_dir", x, d_EVIL_NUM);
    fpg_ps->getWithDefault("y_dir", y, d_EVIL_NUM);
    fpg_ps->getWithDefault("z_dir", z, d_EVIL_NUM);
    d_fixedPressGrad = Vector(x, y, z);
  }

  //__________________________________
  // Pull out Initial Conditions
  ProblemSpecP matProp_ps = nullptr;
  bool isRestart{ false };
  if (prob_spec->findBlockWithOutAttribute("MaterialProperties")) {
    matProp_ps = prob_spec->findBlockWithOutAttribute("MaterialProperties");
  } else if (restart_prob_spec) {
    isRestart = true;
    matProp_ps =
      restart_prob_spec->findBlockWithOutAttribute("MaterialProperties");
  }

  if (!matProp_ps) {
    throw ProblemSetupException(
      "\n Could not find the <MaterialProperties> section in the input file\n",
      __FILE__,
      __LINE__);
  }

  ProblemSpecP ice_ps = matProp_ps->findBlock("ICE");
  if (!ice_ps) {
    throw ProblemSetupException(
      "\n Could not find the <ICE> <MaterialProperties> section in the input "
      "file\n",
      __FILE__,
      __LINE__);
  }

  for (ProblemSpecP mat_ps = ice_ps->findBlock("material"); mat_ps != 0;
       mat_ps              = mat_ps->findNextBlock("material")) {
    std::string index("");
    mat_ps->getAttribute("index", index);

    // create a new ICE material and initalize it
    std::shared_ptr<ICEMaterial> mat =
      std::make_shared<ICEMaterial>(mat_ps, d_materialManager, isRestart);
    if (index != "") {
      d_materialManager->registerMaterial("ICE", mat, std::stoi(index));
    } else {
      d_materialManager->registerMaterial("ICE", mat);
    }
    if (mat->isSurroundingMatl()) {
      d_surroundingMatl_indx = mat->getDWIndex();
    }
  }

  //_________________________________
  // Exchange Coefficients
  proc0cout << "numMatls " << d_materialManager->getNumMaterials() << std::endl;
  d_exchModel =
    ExchangeFactory::create(matProp_ps, d_materialManager, d_with_mpm);
  d_exchModel->problemSetup(matProp_ps);

  // Set up turbulence models - needs to be done after materials are initialized
  d_turbulence = TurbulenceFactory::create(cfd_ice_ps, d_materialManager);

  // Set up wall shaer stress model
  d_WallShearStressModel =
    WallShearStressFactory::create(cfd_ice_ps, d_materialManager);

  //__________________________________
  //  conservationTest
  if (d_output->isLabelSaved("TotalMass")) {
    d_conservationTest->mass  = true;
    d_conservationTest->onOff = true;
  }
  if (d_output->isLabelSaved("TotalMomentum")) {
    d_conservationTest->momentum = true;
    d_conservationTest->onOff    = true;
  }
  if (d_output->isLabelSaved("TotalIntEng") ||
      d_output->isLabelSaved("KineticEnergy")) {
    d_conservationTest->energy = true;
    d_conservationTest->onOff  = true;
  }
  if (d_output->isLabelSaved("eng_exch_error") ||
      d_output->isLabelSaved("mom_exch_error")) {
    d_conservationTest->exchange = true;
    d_conservationTest->onOff    = true;
  }

  //__________________________________
  // WARNINGS
  if (d_impICE && getDelTMaxIncrease() > 10) {
    proc0cout << "\n \n W A R N I N G: " << std::endl;
    proc0cout << " When running implicit ICE you should specify " << std::endl;
    proc0cout << " \t \t <max_delt_increase> to ~2.0 " << std::endl;
    proc0cout << " to a) prevent rapid fluctuations in the timestep and "
              << std::endl;
    proc0cout << "    b) to prevent outflux Vol > cell volume \n \n"
              << std::endl;
  }

  //__________________________________
  //  Custom BC setup
  d_BC_globalVars->d_gravity             = d_gravity;
  d_BC_globalVars->materialManager       = d_materialManager;
  d_BC_globalVars->applyHydrostaticPress = d_applyHydrostaticPress;

  IntVector periodic = grid->getLevel(0)->getPeriodicBoundaries();
  bool periodicDomain =
    (periodic.x() == 1 && periodic.y() == 1 && periodic.z() == 1);

  if (periodicDomain == false) {

    d_BC_globalVars->usingLodi = read_LODI_BC_inputs(
      prob_spec, d_materialManager, d_BC_globalVars->lodi.get());
    d_BC_globalVars->usingMicroSlipBCs =
      read_MicroSlip_BC_inputs(prob_spec, d_BC_globalVars->slip.get());
    d_BC_globalVars->using_MMS_BCs =
      read_MMS_BC_inputs(prob_spec, d_BC_globalVars->mms.get());
    d_BC_globalVars->using_Sine_BCs =
      read_Sine_BC_inputs(prob_spec, d_BC_globalVars->sine.get());
    d_BC_globalVars->using_inletVel_BCs = read_inletVel_BC_inputs(
      prob_spec, d_materialManager, d_BC_globalVars->inletVel.get(), grid);

    //__________________________________
    //  boundary condition warnings
    BC_bulletproofing(prob_spec, d_materialManager, grid);
  }

  //__________________________________
  //  Load Model info.
  // If we are doing a restart, then use the "timestep.xml"
  ProblemSpecP orig_or_restart_ps = 0;
  if (prob_spec->findBlockWithOutAttribute("MaterialProperties")) {
    orig_or_restart_ps = prob_spec;
  } else if (restart_prob_spec) {
    orig_or_restart_ps = restart_prob_spec;
  }
  if (!orig_or_restart_ps) {
    throw ProblemSetupException(
      "\n Could not find the <MaterialProperties> section in the input file\n",
      __FILE__,
      __LINE__);
  }

  ProblemSpecP model_spec = orig_or_restart_ps->findBlock("Models");
  if (model_spec) {

    // Clean up the old models. NOTE This may not be neccessary as the
    // ProblemSpec is only called once EXCEPT when doing switching.
    // Also each application manages its own models. As such, they
    // probably do not need to be regenerated after a switch.
    for (auto& model : d_models) {
      model->releaseComponents();
    }

    d_models = ModelFactory::makeModels(
      d_myworld, d_materialManager, orig_or_restart_ps, prob_spec);

    // Problem setup for each model
    for (auto& model : d_models) {
      model->setComponents(dynamic_cast<SimulationInterface*>(this));
      model->problemSetup(grid, isRestart);

      // Model with transported variables. Bullet proofing each
      // transported variable must have a boundary condition.
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model && fb_model->d_transVars.size()) {
        for (auto& tvar : fb_model->d_transVars) {
          std::string q_CC_name = tvar->var->getName();
          is_BC_specified(prob_spec, q_CC_name, tvar->matls);
        }
      }

      HEChemModel* hec_model = dynamic_cast<HEChemModel*>(model.get());
      if (hec_model) {
        hec_model->setParticleGhostLayer(particle_ghost_type,
                                         particle_ghost_layer);
      }
    }
  }

  //__________________________________
  //  Set up data analysis modules
  if (!d_with_mpm) {
    d_analysisModules =
      AnalysisModuleFactory::create(d_myworld, d_materialManager, prob_spec);

    std::vector<std::vector<const VarLabel*>> dummy;
    for (auto& am : d_analysisModules) {
      am->problemSetup(prob_spec, restart_prob_spec, grid, dummy, dummy);
    }
  } // mpm
}

/*______________________________________________________________________
  Function~  ICE::outputProblemSpec--
  Purpose~   outputs material state
  _____________________________________________________________________*/
void
ICE::outputProblemSpec(ProblemSpecP& root_ps)
{
  DOUTR(m_ice_tasks, " ICE::outputProblemSpec ");

  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP mat_ps = root->findBlockWithOutAttribute("MaterialProperties");
  if (mat_ps == nullptr) {
    mat_ps = root->appendChild("MaterialProperties");
  }

  ProblemSpecP ice_ps = mat_ps->appendChild("ICE");
  for (size_t i = 0; i < d_materialManager->getNumMaterials("ICE"); i++) {
    ICEMaterial* mat =
      static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", i));
    mat->outputProblemSpec(ice_ps);
  }

  // Output the exchange model
  d_exchModel->outputProblemSpec(mat_ps);

  // Models
  ProblemSpecP models_ps = root->appendChild("Models");
  for (auto& model : d_models) {
    model->outputProblemSpec(models_ps);
  }

  // Output data analysis modules
  if (!d_with_mpm) {
    for (auto& am : d_analysisModules) {
      am->outputProblemSpec(root);
    }
  }
}

/* _____________________________________________________________________
   Function~  ICE::scheduleInitialize--
   Notes:     This task actually schedules several tasks.
   _____________________________________________________________________*/
void
ICE::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, m_ice_tasks, " ICE::scheduleInitialize");

  Task* t =
    scinew Task("ICE::actuallyInitialize", this, &ICE::actuallyInitialize);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  t->computes(d_ice_labels->velocity_CCLabel);
  t->computes(d_ice_labels->rho_CCLabel);
  t->computes(d_ice_labels->temperature_CCLabel);
  t->computes(d_ice_labels->specificVolume_CCLabel);
  t->computes(d_ice_labels->vol_frac_CCLabel);
  t->computes(d_ice_labels->rho_micro_CCLabel);
  t->computes(d_ice_labels->speedSound_CCLabel);
  t->computes(d_ice_labels->thermalCondLabel);
  t->computes(d_ice_labels->viscosityLabel);
  t->computes(d_ice_labels->gammaLabel);
  t->computes(d_ice_labels->specific_heatLabel);
  t->computes(d_ice_labels->press_CCLabel, d_press_matl, oims);

  const auto ice_matls = d_materialManager->allMaterials("ICE");
  sched->addTask(t, level->eachPatch(), ice_matls);

  if (d_impICE) {
    d_solver->scheduleInitialize(level, sched, ice_matls);
  }

  //  Exchange Model
  d_exchModel->sched_initialize(sched, level);

  //  Wall shear stress model initialization
  if (d_WallShearStressModel) {
    d_WallShearStressModel->sched_Initialize(sched, level, ice_matls);
  }

  // Models Initialization
  for (auto& model : d_models) {
    model->scheduleInitialize(sched, level);
  }

  // DataAnalysis
  for (auto& am : d_analysisModules) {
    am->scheduleInitialize(sched, level);
  }

  // Make adjustments to the hydrostatic pressure
  // and temperature fields.  You need to do this
  // after the models have initialized the flowfield
  Vector grav              = getGravity();
  const auto ice_matls_sub = ice_matls->getUnion();

  if (grav.length() > 0 && d_applyHydrostaticPress) {

    printSchedule(level, m_ice_tasks, "ICE::initializeSubTask_hydrostaticAdj");

    Task* t2 = scinew Task("ICE::initializeSubTask_hydrostaticAdj",
                           this,
                           &ICE::initializeSubTask_hydrostaticAdj);

    t2->needs(Task::NewDW, d_ice_labels->gammaLabel, ice_matls_sub, d_gn);
    t2->needs(
      Task::NewDW, d_ice_labels->specific_heatLabel, ice_matls_sub, d_gn);

    t2->modifies(d_ice_labels->rho_micro_CCLabel);
    t2->modifies(d_ice_labels->temperature_CCLabel);
    t2->modifies(d_ice_labels->press_CCLabel, d_press_matl, oims);

    sched->addTask(t2, level->eachPatch(), ice_matls);
  }
}

/*!
 * Task:      ICE::scheduleRestartInitialize
 * Purpose:   Set variables that are normally set during the initialization
 *            phase, but get wiped clean when you restart
 */
void
ICE::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
  if (d_impICE) {
    const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");
    d_solver->scheduleRestartInitialize(level, sched, ice_matls);
  }

  // Models Initialization
  for (auto& model : d_models) {
    model->scheduleRestartInitialize(sched, level);
  }

  // dataAnalysis
  for (auto& am : d_analysisModules) {
    am->scheduleRestartInitialize(sched, level);
  }

  // ICE: Material specific flags
  size_t numMatls = d_materialManager->getNumMaterials("ICE");
  for (size_t m = 0; m < numMatls; m++) {
    ICEMaterial* ice_matl =
      static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));

    if (ice_matl->isSurroundingMatl()) {
      d_surroundingMatl_indx = ice_matl->getDWIndex();
    }

    if (ice_matl->getViscosity() > 0.0) {
      d_viscousFlow = true;
    }
  }

  // bulletproofing
  Vector grav = getGravity();
  if (grav.length() > 0.0 && d_surroundingMatl_indx == -9 &&
      d_applyHydrostaticPress) {
    throw ProblemSetupException(
      "ERROR ICE::restartInitialize \n"
      "You must have \n"
      "       <isSurroundingMatl> true </isSurroundingMatl> \n "
      "specified inside the ICE material that is the background matl\n",
      __FILE__,
      __LINE__);
  }
}

/* _____________________________________________________________________
   Function~  ICE::scheduleComputeStableTimestep--
   _____________________________________________________________________*/
void
ICE::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, m_ice_tasks, " ICE::scheduleComputeStableTimestep");

  Task* t = scinew Task("ICE::actuallyComputeStableTimestep",
                        this,
                        &ICE::actuallyComputeStableTimestep);

  const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");
  Task::SearchTG OldTG         = Task::SearchTG::OldTG;

  t->needs(Task::NewDW, d_ice_labels->velocity_CCLabel, d_gac, 1, OldTG);
  t->needs(Task::NewDW, d_ice_labels->speedSound_CCLabel, d_gac, 1, OldTG);
  t->needs(Task::NewDW, d_ice_labels->thermalCondLabel, d_gn, 0, OldTG);
  t->needs(Task::NewDW, d_ice_labels->gammaLabel, d_gn, 0, OldTG);
  t->needs(Task::NewDW, d_ice_labels->specific_heatLabel, d_gn, 0, OldTG);
  t->needs(
    Task::NewDW, d_ice_labels->specificVolume_CCLabel, d_gn, 0, OldTG);
  t->needs(Task::NewDW, d_ice_labels->viscosityLabel, d_gn, 0, OldTG);

  t->computes(d_ice_labels->delTLabel, level.get_rep());
  sched->addTask(t, level->eachPatch(), ice_matls);

  //  If model needs to further restrict the timestep
  for (auto& model : d_models) {
    model->scheduleComputeStableTimestep(sched, level);
  }
}

/*!
 *  Function~  ICE::scheduleTimeAdvance--
 */
void
ICE::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level,
                m_ice_tasks,
                " --------------------------------------------------------");
  printSchedule(level, m_ice_tasks, " ICE::scheduleTimeAdvance");

  // for AMR, we need to reset the initial Delt otherwise some unsuspecting
  // level will get the init delt when it didn't compute delt on L0.
  const PatchSet* patches      = level->eachPatch();
  const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");
  const MaterialSet* mpm_matls = d_materialManager->allMaterials("MPM");
  const MaterialSet* all_matls = d_materialManager->allMaterials();

  MaterialSubset* one_matl = d_press_matl;

  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
  const MaterialSubset* mpm_matls_sub = mpm_matls->getUnion();

  if (d_turbulence) {
    // The turblence model is also called directly from
    // accumlateMomentumSourceSinks.
    d_turbulence->scheduleComputeVariance(sched, patches, ice_matls);
  }

  scheduleMaxMach_on_Lodi_BC_Faces(sched, level, ice_matls);

  scheduleComputeThermoTransportProperties(sched, level, ice_matls);

  scheduleComputePressure(sched, patches, d_press_matl, all_matls);

  scheduleComputeTempFC(
    sched, patches, ice_matls_sub, mpm_matls_sub, all_matls);

  scheduleComputeModelSources(sched, level, all_matls);

  scheduleUpdateVolumeFraction(sched, level, d_press_matl, all_matls);

  scheduleComputeVel_FC(
    sched, patches, ice_matls_sub, mpm_matls_sub, d_press_matl, all_matls);

  d_exchModel->sched_PreExchangeTasks(
    sched, patches, ice_matls_sub, mpm_matls_sub, all_matls);

  d_exchModel->sched_AddExch_VelFC(sched,
                                   patches,
                                   ice_matls_sub,
                                   mpm_matls_sub,
                                   all_matls,
                                   d_BC_globalVars.get(),
                                   false);

  if (d_impICE) { //  I M P L I C I T

    scheduleSetupRHS(sched, patches, one_matl, all_matls, false, "computes");

    scheduleCompute_maxRHS(sched, level, one_matl, all_matls);

    scheduleImplicitPressureSolve(sched,
                                  level,
                                  patches,
                                  one_matl,
                                  d_press_matl,
                                  ice_matls_sub,
                                  mpm_matls_sub,
                                  all_matls);

    scheduleComputeDel_P(
      sched, level, patches, one_matl, d_press_matl, all_matls);
  }

  if (!d_impICE) { //  E X P L I C I T
    scheduleComputeDelPressAndUpdatePressCC(
      sched, patches, d_press_matl, ice_matls_sub, mpm_matls_sub, all_matls);
  }

  scheduleComputePressFC(sched, patches, d_press_matl, all_matls);

  scheduleVelTau_CC(sched, patches, ice_matls);

  scheduleViscousShearStress(sched, patches, ice_matls);

  scheduleAccumulateMomentumSourceSinks(
    sched, patches, d_press_matl, ice_matls_sub, mpm_matls_sub, all_matls);

  scheduleAccumulateEnergySourceSinks(
    sched, patches, ice_matls_sub, mpm_matls_sub, d_press_matl, all_matls);

  scheduleComputeLagrangianValues(sched, patches, all_matls);

  d_exchModel->sched_AddExch_Vel_Temp_CC(sched,
                                         patches,
                                         ice_matls_sub,
                                         mpm_matls_sub,
                                         all_matls,
                                         d_BC_globalVars.get());

  scheduleComputeLagrangianSpecificVolume(
    sched, patches, ice_matls_sub, mpm_matls_sub, d_press_matl, all_matls);

  scheduleComputeLagrangian_Transported_Vars(sched, patches, all_matls);

  scheduleAdvectAndAdvanceInTime(sched, patches, ice_matls_sub, all_matls);

  scheduleConservedtoPrimitive_Vars(
    sched, patches, ice_matls_sub, all_matls, "afterAdvection");

#if 0
  scheduleComputeTaskGraphIndex(          sched, level );
#endif
}

/*!
 *  Function~  ICE::scheduleFinalizeTimestep--
 *  This is called after scheduleTimeAdvance and the scheduleCoarsen
 */
void
ICE::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, m_ice_tasks, " ICE::scheduleFinalizeTimestep");

  const PatchSet* patches             = level->eachPatch();
  const MaterialSet* ice_matls        = d_materialManager->allMaterials("ICE");
  const MaterialSet* all_matls        = d_materialManager->allMaterials();
  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();

  scheduleConservedtoPrimitive_Vars(
    sched, patches, ice_matls_sub, all_matls, "finalizeTimestep");

  //  If a particle model is used schedule relocate
  for (auto& model : d_models) {
    ParticleModel* pModel = dynamic_cast<ParticleModel*>(model.get());
    if (pModel) {

      std::vector<std::vector<const VarLabel*>> old_labels;
      std::vector<std::vector<const VarLabel*>> new_labels;

      MaterialSubset* new_mss = scinew MaterialSubset();
      new_mss->addReference();

      old_labels.push_back(pModel->d_oldLabels);
      new_labels.push_back(pModel->d_newLabels);
      new_mss->addSubset(pModel->d_matl_mss);

      //__________________________________
      //  create a new material set containing the
      //  the updated matlSubset.
      MaterialSet* newMatlSet = scinew MaterialSet();
      newMatlSet->addSubset(new_mss);
      newMatlSet->addReference();

      sched->scheduleParticleRelocation(level,
                                        pModel->pXLabel_preReloc,
                                        old_labels,
                                        pModel->pXLabel,
                                        new_labels,
                                        pModel->pIDLabel,
                                        newMatlSet);

      if (newMatlSet && newMatlSet->removeReference()) {
        delete newMatlSet;
      }
    }
  }
}

/*!
 *  Task:     ICE::scheduleAnalysis--
 *  Purpose:  This is called after ALL other tasks have completed.
 */
void
ICE::scheduleAnalysis(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, m_ice_tasks, " ICE::scheduleAnalysis");

  const PatchSet* patches       = level->eachPatch();
  const MaterialSet* ice_matls  = d_materialManager->allMaterials("ICE");
  const MaterialSet* all_matls  = d_materialManager->allMaterials();
  const MaterialSubset* ice_mss = ice_matls->getUnion();

  // On the fly analysis
  for (auto& am : d_analysisModules) {
    am->scheduleDoAnalysis(sched, level);
  }

  scheduleTestConservation(sched, patches, ice_mss, all_matls);
}

/*!
 *  Function~  ICE::scheduleComputeThermoTransportProperties--
 */
void
ICE::scheduleComputeThermoTransportProperties(SchedulerP& sched,
                                              const LevelP& level,
                                              const MaterialSet* ice_matls)
{
  printSchedule(
    level, m_ice_tasks, " ICE::schedulecomputeThermoTransportProperties");

  Task* t = scinew Task("ICE::computeThermoTransportProperties",
                        this,
                        &ICE::computeThermoTransportProperties);

  t->needs(Task::OldDW,
              d_ice_labels->temperature_CCLabel,
              ice_matls->getUnion(),
              d_gn,
              0);

  t->computes(d_ice_labels->viscosityLabel);
  t->computes(d_ice_labels->thermalCondLabel);
  t->computes(d_ice_labels->gammaLabel);
  t->computes(d_ice_labels->specific_heatLabel);

  sched->addTask(t, level->eachPatch(), ice_matls);

  //  Each model *can* modify the properties
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model && fb_model->computesThermoTransportProps()) {
      fb_model->scheduleModifyThermoTransportProperties(
        sched, level, ice_matls);
    }
  }
}

/*!
 *  Function~  ICE::scheduleComputePressure--
 *  Purpose:  Compute the equilibration pressure step
 */
void
ICE::scheduleComputePressure(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSubset* press_matl,
                             const MaterialSet* ice_matls)
{
  printSchedule(
    patches, m_ice_tasks, " ICE::scheduleComputeEquilibrationPressure");

  Task* t = nullptr;
  if (d_materialManager->getNumMaterials() == 1) {
    t = scinew Task("ICE::computeEquilPressure_1_matl",
                    this,
                    &ICE::computeEquilPressure_1_matl);
  } else {
    t = scinew Task("ICE::computeEquilibrationPressure",
                    this,
                    &ICE::computeEquilibrationPressure);
  }

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  t->needs(Task::OldDW, d_ice_labels->timeStepLabel);
  t->needs(Task::OldDW, d_ice_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(Task::OldDW, d_ice_labels->press_CCLabel, press_matl, oims, d_gn);
  t->needs(Task::OldDW, d_ice_labels->rho_CCLabel, d_gn);
  t->needs(Task::OldDW, d_ice_labels->temperature_CCLabel, d_gn);
  t->needs(Task::OldDW, d_ice_labels->specificVolume_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->gammaLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->specific_heatLabel, d_gn);

  t->computes(d_ice_labels->f_theta_CCLabel);
  t->computes(d_ice_labels->speedSound_CCLabel);
  t->computes(d_ice_labels->vol_frac_CCLabel);
  t->computes(d_ice_labels->specificVolume_CCLabel);
  t->computes(d_ice_labels->rho_CCLabel);
  t->computes(d_ice_labels->compressibilityLabel);
  t->computes(d_ice_labels->sumKappaLabel, press_matl, oims);
  t->computes(d_ice_labels->press_equil_CCLabel, press_matl, oims);
  //  initialized for implicit
  t->computes(d_ice_labels->sum_imp_delPLabel, press_matl, oims);

  computesRequires_CustomBCs(t,
                             "EqPress",
                             d_ice_labels.get(),
                             ice_matls->getUnion(),
                             d_BC_globalVars.get());

  sched->addTask(t, patches, ice_matls);
}

/*!
 *  Function~  ICE::scheduleComputeTempFC--
 */
void
ICE::scheduleComputeTempFC(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSubset* ice_matls,
                           const MaterialSubset* mpm_matls,
                           const MaterialSet* all_matls)
{
  if (d_models.size() > 0) {
    printSchedule(patches, m_ice_tasks, " ICE::computeTempFC");

    Task* t = scinew Task("ICE::computeTempFC", this, &ICE::computeTempFC);

    t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, /*all_matls*/ d_gac, 1);
    t->needs(
      Task::OldDW, d_ice_labels->temperature_CCLabel, ice_matls, d_gac, 1);

    if (mpm_matls) {
      t->needs(
        Task::NewDW, d_ice_labels->temperature_CCLabel, mpm_matls, d_gac, 1);
    }

    t->computes(d_ice_labels->TempX_FCLabel);
    t->computes(d_ice_labels->TempY_FCLabel);
    t->computes(d_ice_labels->TempZ_FCLabel);

    sched->addTask(t, patches, all_matls);
  }
}

/*!
 *  Function~  ICE::scheduleComputeVel_FC--
 */
void
ICE::scheduleComputeVel_FC(SchedulerP& sched,
                           const PatchSet* patches,
                           const MaterialSubset* ice_matls,
                           const MaterialSubset* mpm_matls,
                           const MaterialSubset* press_matl,
                           const MaterialSet* all_matls)
{
  printSchedule(patches, m_ice_tasks, " ICE::scheduleComputeVel_FC");

  Task* t = scinew Task("ICE::computeVel_FC", this, &ICE::computeVel_FC);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(
    Task::NewDW, d_ice_labels->press_equil_CCLabel, press_matl, oims, d_gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->specificVolume_CCLabel, /*all_matls*/ d_gac, 1);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, /*all_matls*/ d_gac, 1);
  t->needs(Task::OldDW, d_ice_labels->velocity_CCLabel, ice_matls, d_gac, 1);

  if (mpm_matls) {
    t->needs(
      Task::NewDW, d_ice_labels->velocity_CCLabel, mpm_matls, d_gac, 1);
  }

  t->computes(d_ice_labels->uvel_FCLabel);
  t->computes(d_ice_labels->vvel_FCLabel);
  t->computes(d_ice_labels->wvel_FCLabel);
  t->computes(d_ice_labels->grad_P_XFCLabel);
  t->computes(d_ice_labels->grad_P_YFCLabel);
  t->computes(d_ice_labels->grad_P_ZFCLabel);

  sched->addTask(t, patches, all_matls);
}

/*!
 *  Task:     ICE::scheduleComputeModelSources--
 *  Purpose:  Compute any sources/sinks of mass, momentum, energy or specific
 * volume
 */
void
ICE::scheduleComputeModelSources(SchedulerP& sched,
                                 const LevelP& level,
                                 const MaterialSet* matls)
{
  if (d_models.size() != 0) {

    printSchedule(level, m_ice_tasks, " ICE::scheduleComputeModelSources");

    Task* task =
      scinew Task("ICE::zeroModelSources", this, &ICE::zeroModelSources);

    task->computes(d_ice_labels->modelMass_srcLabel);
    task->computes(d_ice_labels->modelMom_srcLabel);
    task->computes(d_ice_labels->modelEng_srcLabel);
    task->computes(d_ice_labels->modelVol_srcLabel);

    // Model with transported variables.
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model && fb_model->d_transVars.size()) {
        for (auto& tvar : fb_model->d_transVars) {
          if (tvar->src) {
            task->computes(tvar->src, tvar->matls);
          }
        }
      }
    }

    sched->addTask(task, level->eachPatch(), matls);

    //  Models *can* compute their resources
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model) {
        fb_model->scheduleComputeModelSources(sched, level);
      }

      HEChemModel* hec_model = dynamic_cast<HEChemModel*>(model.get());
      if (hec_model) {
        hec_model->scheduleComputeModelSources(sched, level);
      }

      SolidReactionModel* sr_model =
        dynamic_cast<SolidReactionModel*>(model.get());
      if (sr_model) {
        sr_model->scheduleComputeModelSources(sched, level);
      }

      ParticleModel* p_model = dynamic_cast<ParticleModel*>(model.get());
      if (p_model) {
        p_model->scheduleComputeModelSources(sched, level);
      }
    }
  }
}

/*!
 *  Function~  ICE::scheduleUpdateVolumeFraction--
 */
void
ICE::scheduleUpdateVolumeFraction(SchedulerP& sched,
                                  const LevelP& level,
                                  const MaterialSubset* press_matl,
                                  const MaterialSet* matls)
{
  if (d_models.size() != 0) {
    printSchedule(level, m_ice_tasks, " ICE::scheduleUpdateVolumeFraction");

    Task* task = scinew Task(
      "ICE::updateVolumeFraction", this, &ICE::updateVolumeFraction);

    task->needs(Task::NewDW, d_ice_labels->specificVolume_CCLabel, d_gn);
    task->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gn);
    task->needs(Task::NewDW, d_ice_labels->modelVol_srcLabel, d_gn);
    task->needs(Task::NewDW, d_ice_labels->compressibilityLabel, d_gn);
    task->modifies(d_ice_labels->sumKappaLabel, press_matl);
    task->modifies(d_ice_labels->vol_frac_CCLabel);
    task->modifies(d_ice_labels->f_theta_CCLabel);

    sched->addTask(task, level->eachPatch(), matls);
  }
}

/*!
 *  Function~  ICE::scheduleComputeDelPressAndUpdatePressCC--
 */
void
ICE::scheduleComputeDelPressAndUpdatePressCC(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSubset* press_matl,
  const MaterialSubset* ice_matls,
  const MaterialSubset* /*mpm_matls*/,
  const MaterialSet* matls)
{
  printSchedule(
    patches, m_ice_tasks, " ICE::scheduleComputeDelPressAndUpdatePressCC");
  Task* task = scinew Task("ICE::computeDelPressAndUpdatePressCC",
                           this,
                           &ICE::computeDelPressAndUpdatePressCC);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  task->needs(Task::OldDW, d_ice_labels->timeStepLabel);
  task->needs(Task::OldDW, d_ice_labels->simulationTimeLabel);
  task->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  task->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->uvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->vvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->wvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->specificVolume_CCLabel, d_gn);
  task->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gn);
  task->needs(Task::NewDW, d_ice_labels->speedSound_CCLabel, d_gn);
  task->needs(
    Task::NewDW, d_ice_labels->sumKappaLabel, press_matl, oims, d_gn);
  task->needs(
    Task::NewDW, d_ice_labels->press_equil_CCLabel, press_matl, oims, d_gn);
  //__________________________________
  if (d_models.size() > 0) {
    task->needs(Task::NewDW, d_ice_labels->modelMass_srcLabel, d_gn);
  }

  computesRequires_CustomBCs(task,
                             "update_press_CC",
                             d_ice_labels.get(),
                             ice_matls,
                             d_BC_globalVars.get());

  task->computes(d_ice_labels->press_CCLabel, press_matl, oims);
  task->computes(d_ice_labels->delP_DilatateLabel, press_matl, oims);
  task->computes(d_ice_labels->delP_MassXLabel, press_matl, oims);
  task->computes(d_ice_labels->term2Label, press_matl, oims);
  task->computes(d_ice_labels->sum_rho_CCLabel, press_matl, oims);
  task->computes(d_ice_labels->vol_fracX_FCLabel);
  task->computes(d_ice_labels->vol_fracY_FCLabel);
  task->computes(d_ice_labels->vol_fracZ_FCLabel);

  task->computes(VarLabel::find(abortTimestep_name));
  task->computes(VarLabel::find(recomputeTimestep_name));

  sched->addTask(task, patches, matls);
}

/*!
 *  Function~  ICE::scheduleComputePressFC--
 */
void
ICE::scheduleComputePressFC(SchedulerP& sched,
                            const PatchSet* patches,
                            const MaterialSubset* press_matl,
                            const MaterialSet* matls)
{
  printSchedule(patches, m_ice_tasks, " ICE::scheduleComputePressFC");

  Task* task = scinew Task("ICE::computePressFC", this, &ICE::computePressFC);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  task->needs(
    Task::NewDW, d_ice_labels->press_CCLabel, press_matl, oims, d_gac, 1);
  task->needs(
    Task::NewDW, d_ice_labels->sum_rho_CCLabel, press_matl, oims, d_gac, 1);

  task->computes(d_ice_labels->pressX_FCLabel, press_matl, oims);
  task->computes(d_ice_labels->pressY_FCLabel, press_matl, oims);
  task->computes(d_ice_labels->pressZ_FCLabel, press_matl, oims);

  sched->addTask(task, patches, matls);
}

/*!
 *  Task:    scheduleVelTau_CC
 *  Purpose: Modify the vel_CC in the extra cells so that it behaves
 *           vel_FC[FC] = (vel_CC(c) + vel_CC(ec) )/2
 */
void
ICE::scheduleVelTau_CC(SchedulerP& sched,
                       const PatchSet* patches,
                       const MaterialSet* ice_matls)
{
  if (!d_viscousFlow) {
    return;
  }
  printSchedule(patches, m_ice_tasks, " ICE::scheduleVelTau_CC");

  Task* t = scinew Task("ICE::VelTau_CC", this, &ICE::VelTau_CC);

  t->needs(Task::OldDW, d_ice_labels->velocity_CCLabel, d_gn, 0);
  t->computes(d_ice_labels->velTau_CCLabel);

  sched->addTask(t, patches, ice_matls);
}

/*!
 *  Task:    sscheduleViscousShearStress
 *  Purpose: Compute the turbulent viscosity and the new total viscosity
 */
void
ICE::scheduleViscousShearStress(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* ice_matls)
{
  printSchedule(patches, m_ice_tasks, " ICE::scheduleViscousShearStress");

  Task* t =
    scinew Task("ICE::viscousShearStress", this, &ICE::viscousShearStress);

  if (d_viscousFlow) {
    t->needs(Task::NewDW, d_ice_labels->viscosityLabel, d_gac, 2);
    t->needs(Task::NewDW, d_ice_labels->velTau_CCLabel, d_gac, 2);
    t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gac, 2);
    t->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gac, 2);

    t->computes(d_ice_labels->tau_X_FCLabel);
    t->computes(d_ice_labels->tau_Y_FCLabel);
    t->computes(d_ice_labels->tau_Z_FCLabel);
  }

  if (d_turbulence) {
    t->needs(Task::NewDW, d_ice_labels->uvel_FCMELabel, d_gac, 3);
    t->needs(Task::NewDW, d_ice_labels->vvel_FCMELabel, d_gac, 3);
    t->needs(Task::NewDW, d_ice_labels->wvel_FCMELabel, d_gac, 3);
    t->computes(d_ice_labels->turb_viscosity_CCLabel);
    t->computes(d_ice_labels->total_viscosity_CCLabel);
  }

  if (d_WallShearStressModel) {
    d_WallShearStressModel->sched_AddComputeRequires(t, ice_matls->getUnion());
  }

  //__________________________________
  // bulletproofing
  if ((d_viscousFlow == 0.0 && d_turbulence) ||
      (d_viscousFlow == 0.0 && d_WallShearStressModel)) {
    string warn = "\nERROR:ICE:viscousShearStress\n The viscosity can't be 0 "
                  "when using a turbulence model or a wall shear stress model";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  t->computes(d_ice_labels->viscous_src_CCLabel);
  sched->addTask(t, patches, ice_matls);
}

/* _____________________________________________________________________
   Function~  ICE::scheduleAccumulateMomentumSourceSinks--
   _____________________________________________________________________*/
void
ICE::scheduleAccumulateMomentumSourceSinks(
  SchedulerP& sched,
  const PatchSet* patches,
  const MaterialSubset* press_matl,
  const MaterialSubset* ice_matls,
  const MaterialSubset* /*mpm_matls_sub*/,
  const MaterialSet* matls)
{
  printSchedule(
    patches, m_ice_tasks, " ICE::scheduleAccumulateMomentumSourceSinks");

  Task* t = scinew Task("ICE::accumulateMomentumSourceSinks",
                        this,
                        &ICE::accumulateMomentumSourceSinks);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(
    Task::NewDW, d_ice_labels->pressX_FCLabel, press_matl, oims, d_gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->pressY_FCLabel, press_matl, oims, d_gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->pressZ_FCLabel, press_matl, oims, d_gac, 1);
  t->needs(
    Task::NewDW, d_ice_labels->viscous_src_CCLabel, ice_matls, d_gn, 0);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gn, 0);
  t->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gn, 0);

  t->computes(d_ice_labels->mom_source_CCLabel);
  sched->addTask(t, patches, matls);
}

/*!
 *  Function~  ICE::scheduleAccumulateEnergySourceSinks--
 */
void
ICE::scheduleAccumulateEnergySourceSinks(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSubset* ice_matls,
                                         [[maybe_unused]] const MaterialSubset* mpm_matls,
                                         const MaterialSubset* press_matl,
                                         const MaterialSet* matls)

{
  printSchedule(
    patches, m_ice_tasks, " ICE::scheduleAccumulateEnergySourceSinks");

  Task* t = scinew Task("ICE::accumulateEnergySourceSinks",
                        this,
                        &ICE::accumulateEnergySourceSinks);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  t->needs(Task::OldDW, d_ice_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(Task::NewDW, d_ice_labels->press_CCLabel, press_matl, oims, d_gn);
  t->needs(
    Task::NewDW, d_ice_labels->delP_DilatateLabel, press_matl, oims, d_gn);
  t->needs(Task::NewDW, d_ice_labels->compressibilityLabel, d_gn);
  t->needs(
    Task::OldDW, d_ice_labels->temperature_CCLabel, ice_matls, d_gac, 1);
  t->needs(Task::NewDW, d_ice_labels->thermalCondLabel, ice_matls, d_gac, 1);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gac, 1);
  t->needs(Task::NewDW, d_ice_labels->specificVolume_CCLabel, d_gac, 1);
  t->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gac, 1);

  if (d_with_mpm) {
    t->needs(Task::NewDW, d_ice_labels->TMV_CCLabel, press_matl, oims, d_gn);
  }

  t->computes(d_ice_labels->int_eng_source_CCLabel);
  t->computes(d_ice_labels->heatCond_src_CCLabel);

  sched->addTask(t, patches, matls);
}

/*!
 *  Function~  ICE:: scheduleComputeLagrangianValues--
 *  Purpose:   compute the lagrangian momentum, internal energy and mass
 *  Note:      Only loop over ICE materials
 */
void
ICE::scheduleComputeLagrangianValues(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* ice_matls)
{
  printSchedule(patches, m_ice_tasks, " ICE::scheduleComputeLagrangianValues");

  Task* t = scinew Task(
    "ICE::computeLagrangianValues", this, &ICE::computeLagrangianValues);

  t->needs(Task::NewDW, d_ice_labels->specific_heatLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gn);
  t->needs(Task::OldDW, d_ice_labels->velocity_CCLabel, d_gn);
  t->needs(Task::OldDW, d_ice_labels->temperature_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->mom_source_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->int_eng_source_CCLabel, d_gn);

  if (d_models.size() > 0) {
    t->needs(Task::NewDW, d_ice_labels->modelMass_srcLabel, d_gn);
    t->needs(Task::NewDW, d_ice_labels->modelMom_srcLabel, d_gn);
    t->needs(Task::NewDW, d_ice_labels->modelEng_srcLabel, d_gn);
  }

  t->computes(d_ice_labels->mom_L_CCLabel);
  t->computes(d_ice_labels->int_eng_L_CCLabel);
  t->computes(d_ice_labels->mass_L_CCLabel);

  sched->addTask(t, patches, ice_matls);
}

/*!
 *  Function~  ICE:: scheduleComputeLagrangianSpecificVolume--
 */
void
ICE::scheduleComputeLagrangianSpecificVolume(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSubset* ice_matls,
                                             const MaterialSubset* mpm_matls,
                                             const MaterialSubset* press_matl,
                                             const MaterialSet* matls)
{
  printSchedule(
    patches, m_ice_tasks, " ICE::scheduleComputeLagrangianSpecificVolume");

  Task* t = scinew Task("ICE::computeLagrangianSpecificVolume",
                        this,
                        &ICE::computeLagrangianSpecificVolume);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.

  t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->needs(Task::OldDW, d_ice_labels->temperature_CCLabel, ice_matls, d_gn);
  t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->specificVolume_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->Tdot_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->f_theta_CCLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->compressibilityLabel, d_gn);
  t->needs(Task::NewDW, d_ice_labels->specific_heatLabel, ice_matls, d_gn);
  t->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gac, 1);

  if (mpm_matls) {
    t->needs(
      Task::NewDW, d_ice_labels->temperature_CCLabel, mpm_matls, d_gn);
  }

  t->needs(
    Task::NewDW, d_ice_labels->delP_DilatateLabel, press_matl, oims, d_gn);
  t->needs(Task::NewDW, d_ice_labels->press_CCLabel, press_matl, oims, d_gn);
  if (d_with_mpm) {
    t->needs(Task::NewDW, d_ice_labels->TMV_CCLabel, press_matl, oims, d_gn);
  }

  if (d_models.size() > 0) {
    t->needs(Task::NewDW, d_ice_labels->modelVol_srcLabel, d_gn);
  }

  t->computes(d_ice_labels->sp_vol_L_CCLabel);
  t->computes(d_ice_labels->sp_vol_src_CCLabel);

  t->computes(VarLabel::find(abortTimestep_name));
  t->computes(VarLabel::find(recomputeTimestep_name));

  sched->addTask(t, patches, matls);
}

/*!
 *  Function~  ICE:: scheduleComputeTransportedLagrangianValues--
 *  Purpose:   For each transported variable compute the lagrangian value
 *             q_L_CC = (q_old + q_src) * mass_L
 *  Note:      Be care
 */
void
ICE::scheduleComputeLagrangian_Transported_Vars(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  bool haveTransportVars = false;
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model && fb_model->d_transVars.size()) {
      haveTransportVars = true;
      break;
    }
  }

  if (haveTransportVars) {
    printSchedule(
      patches, m_ice_tasks, "ICE::scheduleComputeLagrangian_Transported_Vars");

    Task* t = scinew Task("ICE::computeLagrangian_Transported_Vars",
                          this,
                          &ICE::computeLagrangian_Transported_Vars);

    t->needs(Task::OldDW, d_ice_labels->timeStepLabel);
    t->needs(Task::NewDW, d_ice_labels->mass_L_CCLabel, d_gn);

    // computes and requires for each transported variable
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model) {
        for (auto& tvar : fb_model->d_transVars) {
          // require q_old
          t->needs(Task::OldDW, tvar->var, tvar->matls, d_gn, 0);
          if (tvar->src) { // require q_src
            t->needs(Task::NewDW, tvar->src, tvar->matls, d_gn, 0);
          }
          t->computes(tvar->var_Lagrangian, tvar->matls);
        }
      }
    }

    sched->addTask(t, patches, matls);
  }
}

/*!
 *  Function~  ICE::scheduleMaxMach_on_Lodi_BC_Faces--
 *  Purpose    compute the reducton variable max_mach_<face>
 *             on Lodi boundary faces
 */
void
ICE::scheduleMaxMach_on_Lodi_BC_Faces(SchedulerP& sched,
                                      const LevelP& level,
                                      const MaterialSet* ice_matls)
{
  if (d_BC_globalVars->usingLodi) {
    printSchedule(level, m_ice_tasks, " ICE::scheduleMaxMach_on_Lodi_BC_Faces");

    Task* task = scinew Task(
      "ICE::maxMach_on_Lodi_BC_Faces", this, &ICE::maxMach_on_Lodi_BC_Faces);

    task->needs(Task::OldDW, d_ice_labels->velocity_CCLabel, d_gn);
    task->needs(Task::OldDW, d_ice_labels->speedSound_CCLabel, d_gn);

    // loop over the Lodi face and add computes for maxMach
    for (auto& face : d_BC_globalVars->lodi->LodiFaces) {
      VarLabel* V_Label = getMaxMach_face_VarLabel(face);
      task->computes(V_Label, ice_matls->getUnion());
    }

    sched->addTask(task, level->eachPatch(), ice_matls);
  }
}

/*!
 *  Function~  ICE::computesRequires_AMR_Refluxing--
 */
void
ICE::computesRequires_AMR_Refluxing(Task* task, [[maybe_unused]] const MaterialSet* ice_matls)
{
  DOUTR(m_ice_tasks, "      computesRequires_AMR_Refluxing\n");

  task->computes(d_ice_labels->mass_X_FC_fluxLabel);
  task->computes(d_ice_labels->mass_Y_FC_fluxLabel);
  task->computes(d_ice_labels->mass_Z_FC_fluxLabel);

  task->computes(d_ice_labels->mom_X_FC_fluxLabel);
  task->computes(d_ice_labels->mom_Y_FC_fluxLabel);
  task->computes(d_ice_labels->mom_Z_FC_fluxLabel);

  task->computes(d_ice_labels->sp_vol_X_FC_fluxLabel);
  task->computes(d_ice_labels->sp_vol_Y_FC_fluxLabel);
  task->computes(d_ice_labels->sp_vol_Z_FC_fluxLabel);

  task->computes(d_ice_labels->int_eng_X_FC_fluxLabel);
  task->computes(d_ice_labels->int_eng_Y_FC_fluxLabel);
  task->computes(d_ice_labels->int_eng_Z_FC_fluxLabel);

  // MODELS
  // DON'T require reflux vars from the OldDW.  Since we only require it
  // between subcycles, we don't want to schedule it.  Otherwise it will
  // just cause excess TG work.  The data will all get to the right place.

  // Models with reflux variables
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model) {
      for (auto& rvar : fb_model->d_refluxVars) {
        task->computes(rvar->var_X_FC_flux);
        task->computes(rvar->var_Y_FC_flux);
        task->computes(rvar->var_Z_FC_flux);
      }
    }
  }
}

/*!
 *  Function~  ICE::scheduleAdvectAndAdvanceInTime--
 */
void
ICE::scheduleAdvectAndAdvanceInTime(SchedulerP& sched,
                                    const PatchSet* patches,
                                    [[maybe_unused]] const MaterialSubset* ice_matlsub,
                                    const MaterialSet* ice_matls)
{
  printSchedule(patches, m_ice_tasks, " ICE::scheduleAdvectAndAdvanceInTime");

  Task* task = scinew Task(
    "ICE::advectAndAdvanceInTime", this, &ICE::advectAndAdvanceInTime);

  task->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  task->needs(Task::NewDW, d_ice_labels->uvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->vvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->wvel_FCMELabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->mom_L_ME_CCLabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->mass_L_CCLabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->eng_L_ME_CCLabel, d_gac, 2);
  task->needs(Task::NewDW, d_ice_labels->sp_vol_L_CCLabel, d_gac, 2);

  if (d_doRefluxing) {
    computesRequires_AMR_Refluxing(task, ice_matls);
  }

  task->computes(d_ice_labels->mass_advLabel);
  task->computes(d_ice_labels->mom_advLabel);
  task->computes(d_ice_labels->eng_advLabel);
  task->computes(d_ice_labels->sp_vol_advLabel);

  // Model Variables.
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model) {
      for (auto& tvar : fb_model->d_transVars) {
        task->needs(
          Task::NewDW, tvar->var_Lagrangian, tvar->matls, d_gac, 2);
        task->computes(tvar->var_adv, tvar->matls);
      }
    }
  }
  task->computes(VarLabel::find(abortTimestep_name));
  task->computes(VarLabel::find(recomputeTimestep_name));

  sched->addTask(task, patches, ice_matls);
}

/*!
 *  Function~  ICE::scheduleConservedtoPrimitive_Vars--
 */
void
ICE::scheduleConservedtoPrimitive_Vars(SchedulerP& sched,
                                       const PatchSet* patch_set,
                                       const MaterialSubset* ice_matlsub,
                                       const MaterialSet* ice_matls,
                                       const string& where)
{
  ASSERT(where == "afterAdvection" || where == "finalizeTimestep");

  int levelIndex = getLevel(patch_set)->getIndex();
  int numLevels  = getLevel(patch_set)->getGrid()->numLevels();

  // single level problems we only need to perform this task once
  // immediately after advecton
  if (numLevels == 1 && where == "finalizeTimestep") {
    return;
  }

  // On the finest level we only need to perform this task once
  // immediately after advecton
  if (levelIndex + 1 == numLevels && where == "finalizeTimestep") {
    return;
  }

  // which taskgraph to search for computes
  Task::SearchTG whichTG = Task::SearchTG::NewTG;
  if (where == "finalizeTimestep") {
    whichTG = Task::SearchTG::OldTG;
  }

  std::string name = "ICE::conservedtoPrimitive_Vars:" + where;
  printSchedule(patch_set, m_ice_tasks, name);

  Task* task = scinew Task(name, this, &ICE::conservedtoPrimitive_Vars);

  task->needs(Task::OldDW, d_ice_labels->timeStepLabel);
  task->needs(Task::OldDW, d_ice_labels->simulationTimeLabel);
  task->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patch_set));
  task->needs(Task::NewDW, d_ice_labels->mass_advLabel, d_gn, 0);
  task->needs(Task::NewDW, d_ice_labels->mom_advLabel, d_gn, 0);
  task->needs(Task::NewDW, d_ice_labels->eng_advLabel, d_gn, 0);
  task->needs(Task::NewDW, d_ice_labels->sp_vol_advLabel, d_gn, 0);
  task->needs(
    Task::NewDW, d_ice_labels->specific_heatLabel, d_gn, 0, whichTG);
  task->needs(
    Task::NewDW, d_ice_labels->speedSound_CCLabel, d_gn, 0, whichTG);
  task->needs(Task::NewDW, d_ice_labels->vol_frac_CCLabel, d_gn, 0, whichTG);
  task->needs(Task::NewDW, d_ice_labels->gammaLabel, d_gn, 0, whichTG);

  computesRequires_CustomBCs(
    task, "Advection", d_ice_labels.get(), ice_matlsub, d_BC_globalVars.get());

  task->modifies(d_ice_labels->rho_CCLabel, whichTG);
  task->modifies(d_ice_labels->specificVolume_CCLabel, whichTG);

  if (where == "afterAdvection") {
    task->computes(d_ice_labels->temperature_CCLabel);
    task->computes(d_ice_labels->velocity_CCLabel);
    task->computes(d_ice_labels->machLabel);
  }

  if (where == "finalizeTimestep") {
    task->modifies(d_ice_labels->temperature_CCLabel, whichTG);
    task->modifies(d_ice_labels->velocity_CCLabel, whichTG);
    task->modifies(d_ice_labels->machLabel, whichTG);
  }

  // Model Variables.
  for (auto& model : d_models) {
    FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
    if (fb_model) {
      for (auto& tvar : fb_model->d_transVars) {
        task->needs(Task::NewDW, tvar->var_adv, tvar->matls, d_gn, 0);
        if (where == "afterAdvection") {
          task->computes(tvar->var, tvar->matls);
        }
        if (where == "finalizeTimestep") {
          task->modifies(tvar->var, tvar->matls, whichTG);
        }
      }
    }
  }

  sched->addTask(task, patch_set, ice_matls);
}

/*!
 *  Function~  ICE::scheduleTestConservation--
 *  Purpose:   Compute the KineticEnergy, total internal energy, total mass and
 error in the energy and momentum exchange
 */
void
ICE::scheduleTestConservation(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSubset* ice_matls,
                              const MaterialSet* all_matls)
{
  int levelIndex = getLevel(patches)->getIndex();
  if (d_conservationTest->onOff && levelIndex == 0) {
    printSchedule(patches, m_ice_tasks, " ICE::scheduleTestConservation");

    Task* t =
      scinew Task("ICE::TestConservation", this, &ICE::TestConservation);

    t->needs(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
    t->needs(Task::NewDW, d_ice_labels->rho_CCLabel, ice_matls, d_gn);
    t->needs(Task::NewDW, d_ice_labels->velocity_CCLabel, ice_matls, d_gn);
    t->needs(
      Task::NewDW, d_ice_labels->temperature_CCLabel, ice_matls, d_gn);
    t->needs(Task::NewDW, d_ice_labels->specific_heatLabel, ice_matls, d_gn);
    t->needs(Task::NewDW, d_ice_labels->uvel_FCMELabel, ice_matls, d_gn);
    t->needs(Task::NewDW, d_ice_labels->vvel_FCMELabel, ice_matls, d_gn);
    t->needs(Task::NewDW, d_ice_labels->wvel_FCMELabel, ice_matls, d_gn);

    // A L L  M A T L S
    t->needs(Task::NewDW, d_ice_labels->mom_L_CCLabel, d_gn);
    t->needs(Task::NewDW, d_ice_labels->int_eng_L_CCLabel, d_gn);
    t->needs(Task::NewDW, d_ice_labels->mom_L_ME_CCLabel, d_gn);
    t->needs(Task::NewDW, d_ice_labels->eng_L_ME_CCLabel, d_gn);

    //  Create reductionMatlSubSet that includes all ice matls
    // and the global matl
    const MaterialSubset* global_mss = t->getGlobalMatlSubset();
    MaterialSubset* reduction_mss    = scinew MaterialSubset();
    reduction_mss->add(global_mss->get(0));
    size_t numICEmatls = d_materialManager->getNumMaterials("ICE");
    if (numICEmatls > 1) { // ignore for single matl problems
      for (size_t m = 0; m < numICEmatls; m++) {
        reduction_mss->add(ice_matls->get(m));
      }
    }
    reduction_mss->addReference();

    if (d_conservationTest->exchange) {
      t->computes(d_ice_labels->mom_exch_errorLabel, global_mss);
      t->computes(d_ice_labels->eng_exch_errorLabel, global_mss);
    }
    if (d_conservationTest->mass) {
      t->computes(d_ice_labels->TotalMassLabel, reduction_mss);
    }
    if (d_conservationTest->energy) {
      t->computes(d_ice_labels->KineticEnergyLabel, reduction_mss);
      t->computes(d_ice_labels->TotalIntEngLabel, reduction_mss);
    }
    if (d_conservationTest->momentum) {
      t->computes(d_ice_labels->TotalMomentumLabel, reduction_mss);
    }

    sched->addTask(t, patches, all_matls);

    if (reduction_mss && reduction_mss->removeReference()) {
      delete reduction_mss;
    }
  }

  //  Each model *can* test conservation
  if (d_models.size() != 0) {
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model) {
        fb_model->scheduleTestConservation(sched, patches);
      }
    }
  }
}

/*!
 *  ICE::scheduleComputeTaskGraphIndex--
 */
void
ICE::scheduleComputeTaskGraphIndex([[maybe_unused]] SchedulerP& sched, [[maybe_unused]] const LevelP& level)
{
}

/*!
 *  Function~  ICE::actuallyComputeStableTimestep--
 *  Purpose:  Compute the stable timestep size.  There are two algorithms to
 *            do it
 */
void
ICE::actuallyComputeStableTimestep(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* /*matls*/,
                                   DataWarehouse* /*old_dw*/,
                                   DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::actuallyComputeStableTimestep");

    Vector dx   = patch->dCell();
    double delX = dx.x();
    double delY = dx.y();
    double delZ = dx.z();
    double delt_CFL;
    double delt_diff;
    double delt;
    double inv_sum_invDelx_sqr =
      1.0 / (1.0 / (delX * delX) + 1.0 / (delY * delY) + 1.0 / (delZ * delZ));

    constCCVariable<double> speedSound, sp_vol_CC, thermalCond, viscosity;
    constCCVariable<double> cv, gamma;
    constCCVariable<Vector> vel_CC;

    IntVector badCell(0, 0, 0);
    delt_CFL  = 1000.0;
    delt_diff = 1000;
    delt      = 1000;

    for (size_t m = 0; m < d_materialManager->getNumMaterials("ICE"); m++) {
      ICEMaterial* ice_matl =
        dynamic_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      new_dw->get(
        speedSound, d_ice_labels->speedSound_CCLabel, indx, patch, d_gac, 1);
      new_dw->get(
        vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gac, 1);
      new_dw->get(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        viscosity, d_ice_labels->viscosityLabel, indx, patch, d_gn, 0);
      new_dw->get(
        thermalCond, d_ice_labels->thermalCondLabel, indx, patch, d_gn, 0);
      new_dw->get(gamma, d_ice_labels->gammaLabel, indx, patch, d_gn, 0);
      new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);

      if (d_delT_scheme == "aggressive") { //      A G G R E S S I V E
        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c        = *iter;
          double speed_Sound = d_delT_speedSoundKnob * speedSound[c];

          double A = d_CFL * delX /
                     (speed_Sound + std::abs(vel_CC[c].x()) + d_SMALL_NUM);
          double B = d_CFL * delY /
                     (speed_Sound + std::abs(vel_CC[c].y()) + d_SMALL_NUM);
          double C = d_CFL * delZ /
                     (speed_Sound + std::abs(vel_CC[c].z()) + d_SMALL_NUM);

          delt_CFL = std::min(A, delt_CFL);
          delt_CFL = std::min(B, delt_CFL);
          delt_CFL = std::min(C, delt_CFL);

          if (A < 1e-20 || B < 1e-20 || C < 1e-20) {
            if (badCell == IntVector(0, 0, 0)) {
              badCell = c;
            }
            std::cout << d_myworld->myRank() << " Bad cell " << c << " ("
                      << patch->getID() << "-" << level->getIndex()
                      << "): " << vel_CC[c] << std::endl;
          }
        }
        //      std::cout << " Aggressive delT Based on currant number "<<
        //      delt_CFL << std::endl;
        //__________________________________
        // stability constraint due to diffusion
        //  I C E  O N L Y
        double thermalCond_test = ice_matl->getThermalConductivity();
        double viscosity_test   = ice_matl->getViscosity();

        if (thermalCond_test != 0 || viscosity_test != 0) {

          for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            double cp   = cv[c] * gamma[c];

            double Mod_thermalCond = d_delT_diffusionKnob * thermalCond[c];
            double Mod_viscosity   = d_delT_diffusionKnob * viscosity[c];

            double inv_thermalDiffusivity =
              cp / (sp_vol_CC[c] * Mod_thermalCond + d_SMALL_NUM);
            double kinematicViscosity = Mod_viscosity * sp_vol_CC[c];
            double inv_diffusionCoeff = std::min(
              inv_thermalDiffusivity, 1.0 / (d_SMALL_NUM + kinematicViscosity));

            double A = d_CFL * 0.5 * inv_sum_invDelx_sqr * inv_diffusionCoeff;

            delt_diff = std::min(A, delt_diff);
            if (delt_diff < 1e-20 && badCell == IntVector(0, 0, 0)) {
              badCell = c;
            }
          }
        } //
        //      std::cout << "delT based on diffusion  "<< delt_diff<<endl;
        delt = std::min(delt_CFL, delt_diff);
      } // aggressive Timestep

      if (d_delT_scheme == "conservative") { //      C O N S E R V A T I V E
        //__________________________________
        // Use a characteristic velocity
        // to compute a sweptvolume. The
        // swept volume can't exceed the cell volume
        std::vector<IntVector> adj_offset(3);
        adj_offset[0] = IntVector(1, 0, 0); // X
        adj_offset[1] = IntVector(0, 1, 0); // Y
        adj_offset[2] = IntVector(0, 0, 1); // Z

        Vector faceArea;
        faceArea[0] = dx.y() * dx.z(); // X
        faceArea[1] = dx.x() * dx.z(); // Y
        faceArea[2] = dx.x() * dx.y(); // Z

        double vol  = dx.x() * dx.y() * dx.z();
        Vector grav = getGravity();
        double grav_vel =
          Sqrt(dx.x() * std::abs(grav.x()) + dx.y() * std::abs(grav.y()) +
               dx.z() * std::abs(grav.z()));

        double dx_length = dx.length();

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          double sumSwept_Vol = 0.0;
          IntVector c         = *iter;
          double cp           = cv[c] * gamma[c];

          for (int dir = 0; dir < 3; dir++) { // loop over all three directions
            IntVector L = c - adj_offset[dir];
            IntVector R = c + adj_offset[dir];

            double vel_R = vel_CC[R][dir];
            double vel_C = vel_CC[c][dir];
            double vel_L = vel_CC[L][dir];

            double vel_FC_R = 0.5 * (vel_R + vel_C);
            double vel_FC_L = 0.5 * (vel_L + vel_C);

            double c_L        = speedSound[L];
            double c_R        = speedSound[R];
            double speedSound = std::max(c_L, c_R);

            double relative_vel = fabs(vel_R - vel_L);

            double thermalDiffusivity = thermalCond[c] * sp_vol_CC[c] / cp;
            double diffusion_vel =
              std::max(thermalDiffusivity, viscosity[c]) / dx_length;

            double characteristicVel_R =
              vel_FC_R + d_delT_speedSoundKnob * speedSound + relative_vel +
              grav_vel + diffusion_vel;
            double characteristicVel_L =
              vel_FC_L - d_delT_speedSoundKnob * speedSound - relative_vel -
              grav_vel - diffusion_vel;

            double sweptVol_R = characteristicVel_R * faceArea[dir];
            double sweptVol_L = characteristicVel_L * -faceArea[dir];

            // only compute outflow volumes
            sweptVol_R = std::max(0.0, sweptVol_R);
            sweptVol_L = std::max(0.0, sweptVol_L);
            sumSwept_Vol += sweptVol_R + sweptVol_L;
          } // dir loop

          double delt_tmp = d_CFL * vol / (sumSwept_Vol + d_SMALL_NUM);
          delt            = std::min(delt, delt_tmp);

          if (delt < 1e-20 && badCell == IntVector(0, 0, 0)) {
            badCell = c;
          }
        } // iter loop
        //      std::cout << " Conservative delT based on swept volumes "<<
        //      delt<<endl;
      }
    } // matl loop

    const Level* level = getLevel(patches);
    //__________________________________
    //  Bullet proofing
    if (delt < 1e-20) {
      std::ostringstream warn;
      warn << "ERROR ICE:(L-" << level->getIndex()
           << "):ComputeStableTimestep: delT < 1e-20 on cell " << badCell;
      throw InvalidValue(warn.str(), __FILE__, __LINE__);
    }
    new_dw->put(delt_vartype(delt), d_ice_labels->delTLabel, level);
  } // patch loop
}

/*!
 *  Function~  ICE::actuallyInitialize--
 *  Purpose~   Initialize CC variables and the pressure
 *             Note that rho_micro, sp_vol, temp and velocity must be defined
 *             everywhere in the domain
 */
void
ICE::actuallyInitialize(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* /*matls*/,
                        DataWarehouse*,
                        DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  new_dw->get(timeStep, d_ice_labels->timeStepLabel);

  bool isNotInitialTimestep = (timeStep > 0);

  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();

  // find max index of all the ICE matls
  // you could have a 1 matl problem with a starting indx of 2
  int max_indx = -100;

  size_t numICEMatls = d_materialManager->getNumMaterials("ICE");
  for (size_t m = 0; m < numICEMatls; m++) {
    ICEMaterial* ice_matl =
      static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
    int indx = ice_matl->getDWIndex();
    max_indx = std::max(max_indx, indx);
  }
  max_indx += 1;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::actuallyInitialize");

    size_t numMatls    = d_materialManager->getNumMaterials("ICE");
    size_t numALLMatls = d_materialManager->getNumMaterials();

    Vector grav = getGravity();

    std::vector<constCCVariable<double>> placeHolder(0);
    std::vector<CCVariable<double>> rho_micro(max_indx);
    std::vector<CCVariable<double>> sp_vol_CC(max_indx);
    std::vector<CCVariable<double>> rho_CC(max_indx);
    std::vector<CCVariable<double>> Temp_CC(max_indx);
    std::vector<CCVariable<double>> speedSound(max_indx);
    std::vector<CCVariable<double>> vol_frac_CC(max_indx);
    std::vector<CCVariable<Vector>> vel_CC(max_indx);
    std::vector<CCVariable<double>> cv(max_indx);
    std::vector<CCVariable<double>> gamma(max_indx);

    CCVariable<double> press_CC, vol_frac_sum;

    new_dw->allocateAndPut(press_CC, d_ice_labels->press_CCLabel, 0, patch);
    new_dw->allocateTemporary(vol_frac_sum, patch);
    vol_frac_sum.initialize(0.0);

    //  Thermo and transport properties
    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      CCVariable<double> viscosity, thermalCond;

      new_dw->allocateAndPut(
        viscosity, d_ice_labels->viscosityLabel, indx, patch);
      new_dw->allocateAndPut(
        thermalCond, d_ice_labels->thermalCondLabel, indx, patch);
      new_dw->allocateAndPut(
        cv[indx], d_ice_labels->specific_heatLabel, indx, patch);
      new_dw->allocateAndPut(
        gamma[indx], d_ice_labels->gammaLabel, indx, patch);

      gamma[indx].initialize(ice_matl->getGamma());
      cv[indx].initialize(ice_matl->getSpecificHeat());
      viscosity.initialize(ice_matl->getViscosity());
      thermalCond.initialize(ice_matl->getThermalConductivity());
    }

    // --------bulletproofing
    if (grav.length() > 0.0 && d_surroundingMatl_indx == -9) {
      throw ProblemSetupException(
        "ERROR ICE::actuallyInitialize \n"
        "You must have \n"
        "       <isSurroundingMatl> true </isSurroundingMatl> \n "
        "specified inside the ICE material that is the background matl\n",
        __FILE__,
        __LINE__);
    }

    //__________________________________
    // Note:
    // The press_CC isn't material dependent even though we loop over numMatls
    // below. This is done so we don't need additional machinery to grab the
    // pressure inside a geom_object
    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      new_dw->allocateAndPut(
        rho_micro[indx], d_ice_labels->rho_micro_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_CC[indx], d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        rho_CC[indx], d_ice_labels->rho_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        Temp_CC[indx], d_ice_labels->temperature_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        speedSound[indx], d_ice_labels->speedSound_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        vol_frac_CC[indx], d_ice_labels->vol_frac_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        vel_CC[indx], d_ice_labels->velocity_CCLabel, indx, patch);
    }

    double p_ref = getRefPress();
    press_CC.initialize(p_ref);

    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      ice_matl->initializeCells(rho_micro[indx],
                                rho_CC[indx],
                                Temp_CC[indx],
                                speedSound[indx],
                                vol_frac_CC[indx],
                                vel_CC[indx],
                                press_CC,
                                numALLMatls,
                                patch,
                                new_dw);

      // if specified, overide the initialization
      customInitialization(patch,
                           rho_CC[indx],
                           Temp_CC[indx],
                           vel_CC[indx],
                           press_CC,
                           ice_matl,
                           d_customInitialize_basket.get());

      setBC(rho_CC[indx],
            "Density",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimestep);
      setBC(rho_micro[indx],
            "Density",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimestep);
      setBC(Temp_CC[indx],
            "Temperature",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimestep);
      setBC(speedSound[indx],
            "zeroNeumann",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimestep);
      setBC(vel_CC[indx],
            "Velocity",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimestep);
      setBC(press_CC,
            rho_micro,
            placeHolder,
            d_surroundingMatl_indx,
            "rho_micro",
            "Pressure",
            patch,
            d_materialManager,
            0,
            new_dw,
            isNotInitialTimestep);

      SpecificHeat* cvModel = ice_matl->getSpecificHeatModel();
      if (cvModel != nullptr) {
        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c    = *iter;
          gamma[indx][c] = cvModel->getGamma(Temp_CC[indx][c]);
          cv[indx][c]    = cvModel->getSpecificHeat(Temp_CC[indx][c]);
        }
      }

      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c        = *iter;
        sp_vol_CC[indx][c] = 1.0 / rho_micro[indx][c];

        // needed for LODI BCs
        vol_frac_CC[indx][c] = rho_CC[indx][c] * sp_vol_CC[indx][c];
        vol_frac_sum[c] += vol_frac_CC[indx][c];

        double dp_drho, dp_de, c_2, press_tmp;
        ice_matl->getEOS()->computePressEOS(rho_micro[indx][c],
                                            gamma[indx][c],
                                            cv[indx][c],
                                            Temp_CC[indx][c],
                                            press_tmp,
                                            dp_drho,
                                            dp_de);

        if (!d_customInitialize_basket->doesComputePressure) {
          press_CC[c] = press_tmp;
        }

        c_2 = dp_drho +
              dp_de * press_CC[c] / (rho_micro[indx][c] * rho_micro[indx][c]);
        speedSound[indx][c] = std::sqrt(c_2);
      }

      //____ B U L L E T   P R O O F I N G----
      IntVector neg_cell;
      std::ostringstream warn, base;
      base << "ERROR ICE:(L-" << L_indx << "):actuallyInitialize, mat " << indx
           << " cell ";

      if (!areAllValuesPositive(press_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << base.str() << neg_cell << " position: " << pt
             << " press_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(rho_CC[indx], neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << base.str() << neg_cell << " position: " << pt
             << " rho_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(Temp_CC[indx], neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << base.str() << neg_cell << " position: " << pt
             << " Temp_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(sp_vol_CC[indx], neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << base.str() << neg_cell << " position: " << pt
             << " sp_vol_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(speedSound[indx], neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << base.str() << neg_cell << " position: " << pt
             << ", speedSound is invalid (" << speedSound[indx][neg_cell] << ")"
             << "  Verify that at least one geom_object covers the entire "
                "computational domain.";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
    } // numMatls

    // make sure volume fractions sum to 1
    if (!d_with_mpm) {
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        if (vol_frac_sum[c] > 1.0 + 1e-10 || vol_frac_sum[c] < 1.0 - 1e-10) {
          Point pt = patch->getCellPosition(c);
          std::ostringstream warn, base;
          base << "ERROR ICE:(L-" << L_indx << "):actuallyInitialize";
          warn << base.str() << "Cell: " << c << " position: " << pt
               << " Volume fractions did not sum to 1. Sum=" << vol_frac_sum[c]
               << "\n";
          throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
        }
      }
    }
  } // patch loop
}

/*!
 *  Function~  ICE::initialize_hydrostaticAdj
 *  Purpose~   adjust the pressure and temperature fields after both
 *  ICE and the models have initialized the fields
 */
void
ICE::initializeSubTask_hydrostaticAdj(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* /*ice_matls*/,
                                      DataWarehouse* /*old_dw*/,
                                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::initializeSubTask_hydrostaticAdj");

    // adjust the pressure field
    CCVariable<double> rho_micro, press_CC;
    new_dw->getModifiable(press_CC, d_ice_labels->press_CCLabel, 0, patch);
    new_dw->getModifiable(rho_micro,
                          d_ice_labels->rho_micro_CCLabel,
                          d_surroundingMatl_indx,
                          patch);

    hydrostaticPressureAdjustment(patch, rho_micro, press_CC);
  }
}

/*!
 *  Function~  ICE::computeThermoTransportProperties
 *  Purpose~
 */
void
ICE::computeThermoTransportProperties(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset* /*ice_matls*/,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::computeThermoTransportProperties");

    size_t numMatls = d_materialManager->getNumMaterials("ICE");
    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      constCCVariable<double> temp_CC;
      old_dw->get(
        temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);

      CCVariable<double> viscosity, thermalCond, gamma, cv;
      new_dw->allocateAndPut(
        thermalCond, d_ice_labels->thermalCondLabel, indx, patch);
      new_dw->allocateAndPut(
        viscosity, d_ice_labels->viscosityLabel, indx, patch);
      new_dw->allocateAndPut(cv, d_ice_labels->specific_heatLabel, indx, patch);
      new_dw->allocateAndPut(gamma, d_ice_labels->gammaLabel, indx, patch);

      viscosity.initialize(ice_matl->getViscosity());
      thermalCond.initialize(ice_matl->getThermalConductivity());
      gamma.initialize(ice_matl->getGamma());
      cv.initialize(ice_matl->getSpecificHeat());

      SpecificHeat* cvModel = ice_matl->getSpecificHeatModel();
      if (cvModel != nullptr) {
        // loop through cells and compute pointwise
        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          cv[c]       = cvModel->getSpecificHeat(temp_CC[c]);
          gamma[c]    = cvModel->getGamma(temp_CC[c]);
        }
      }
    }
  }
}

/*!
 *  Function~  ICE::computeEquilibrationPressure--
 *  Purpose~   Find the equilibration pressure
 *  Reference: Flow of Interpenetrating Material Phases, J. Comp, Phys
 *  18, 440-464, 1975, see the equilibration section
 *
 *  Steps
 *  ----------------
 *  - Compute rho_micro_CC, SpeedSound, vol_frac
 *
 *  For each cell
 *  _ WHILE LOOP(convergence, max_iterations)
 *  - compute the pressure and dp_drho from the EOS of each material.
 *  - Compute delta Pressure
 *  - Compute delta volume fraction and update the
 *  volume fraction and the celldensity.
 *  - Test for convergence of delta pressure and delta volume fraction
 *  - END WHILE LOOP
 *  - bulletproofing
 *  end
 *
 *  Note:  The nomenclature follows the reference.
 */
void
ICE::computeEquilibrationPressure(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* /*matls*/,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimestep = (timeStep > 0);

  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computeEquilibrationPressure");

    constexpr double converg_coeff{ 15.0 };
    constexpr double convergence_crit{ converg_coeff * DBL_EPSILON };

    double sum{ 0.0 }, tmp{ 0.0 };

    size_t numMatls = d_materialManager->getNumMaterials("ICE");
    static int n_passes;
    n_passes++;

    std::vector<double> press_eos(numMatls);
    std::vector<double> dp_drho(numMatls), dp_de(numMatls);

    std::vector<CCVariable<double>> vol_frac(numMatls);
    std::vector<CCVariable<double>> rho_micro(numMatls);
    std::vector<CCVariable<double>> rho_CC_new(numMatls);
    std::vector<CCVariable<double>> sp_vol_new(numMatls);
    std::vector<CCVariable<double>> speedSound(numMatls);
    std::vector<CCVariable<double>> speedSound_new(numMatls);
    std::vector<CCVariable<double>> f_theta(numMatls);
    std::vector<CCVariable<double>> kappa(numMatls);

    std::vector<constCCVariable<double>> Temp(numMatls);
    std::vector<constCCVariable<double>> rho_CC(numMatls);
    std::vector<constCCVariable<double>> sp_vol_CC(numMatls);
    std::vector<constCCVariable<double>> cv(numMatls);
    std::vector<constCCVariable<double>> gamma(numMatls);
    std::vector<constCCVariable<double>> placeHolder(0);

    constCCVariable<double> press;
    CCVariable<int> n_iters_equil_press;
    CCVariable<double> press_new, sumKappa, sum_imp_delP;

    old_dw->get(press, d_ice_labels->press_CCLabel, 0, patch, d_gn, 0);
    new_dw->allocateAndPut(
      press_new, d_ice_labels->press_equil_CCLabel, 0, patch);
    new_dw->allocateAndPut(sumKappa, d_ice_labels->sumKappaLabel, 0, patch);
    new_dw->allocateAndPut(
      sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch);

    sum_imp_delP.initialize(0.0); //-- initialize for implicit pressure

    for (size_t m = 0; m < numMatls; m++) {
      int indx = d_materialManager->getMaterial("ICE", m)->getDWIndex();

      old_dw->get(
        Temp[m], d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
      old_dw->get(rho_CC[m], d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
      old_dw->get(sp_vol_CC[m],
                  d_ice_labels->specificVolume_CCLabel,
                  indx,
                  patch,
                  d_gn,
                  0);
      new_dw->get(
        cv[m], d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);
      new_dw->get(gamma[m], d_ice_labels->gammaLabel, indx, patch, d_gn, 0);

      new_dw->allocateTemporary(rho_micro[m], patch);
      new_dw->allocateAndPut(
        vol_frac[m], d_ice_labels->vol_frac_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        rho_CC_new[m], d_ice_labels->rho_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_new[m], d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        f_theta[m], d_ice_labels->f_theta_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        kappa[m], d_ice_labels->compressibilityLabel, indx, patch);
      new_dw->allocateAndPut(
        speedSound_new[m], d_ice_labels->speedSound_CCLabel, indx, patch);
    }

    press_new.copyData(press);

    //__________________________________
    // Compute rho_micro, volfrac
    for (size_t m = 0; m < numMatls; m++) {
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c     = *iter;
        rho_micro[m][c] = 1.0 / sp_vol_CC[m][c];
        vol_frac[m][c]  = rho_CC[m][c] * sp_vol_CC[m][c];
      }
    }

    //______________________________________________________________________
    // Done with preliminary calcs, now loop over every cell
    int count{ 0 }, test_max_iter{ 0 };
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      double delPress{ 0.0 };
      bool converged{ false };
      count = 0;
      std::vector<EqPress_dbg> dbgEqPress;

      while (count < d_max_iter_equilibration && converged == false) {
        count++;

        // evaluate press_eos at cell i,j,k
        for (size_t m = 0; m < numMatls; m++) {
          ICEMaterial* ice_matl =
            static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
          ice_matl->getEOS()->computePressEOS(rho_micro[m][c],
                                              gamma[m][c],
                                              cv[m][c],
                                              Temp[m][c],
                                              press_eos[m],
                                              dp_drho[m],
                                              dp_de[m]);
        }

        // - compute delPress
        // - update press_CC
        double A = 0., B = 0., C = 0.;
        for (size_t m = 0; m < numMatls; m++) {
          double Q     = press_new[c] - press_eos[m];
          double div_y = (vol_frac[m][c] * vol_frac[m][c]) /
                         (dp_drho[m] * rho_CC[m][c] + d_SMALL_NUM);
          A += vol_frac[m][c];
          B += Q * div_y;
          C += div_y;
        }
        double vol_frac_not_close_packed = 1.0;
        delPress = (A - vol_frac_not_close_packed - B) / C;

        press_new[c] += delPress;

        // backout rho_micro_CC at this new pressure
        for (size_t m = 0; m < numMatls; m++) {
          ICEMaterial* ice_matl =
            static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
          rho_micro[m][c] = ice_matl->getEOS()->computeRhoMicro(
            press_new[c], gamma[m][c], cv[m][c], Temp[m][c], rho_micro[m][c]);

          double div = 1. / rho_micro[m][c];

          // - updated volume fractions
          vol_frac[m][c] = rho_CC[m][c] * div;
        }

        // - Test for convergence
        //  If sum of vol_frac_CC ~= vol_frac_not_close_packed then converged
        sum = 0.0;
        for (size_t m = 0; m < numMatls; m++) {
          sum += vol_frac[m][c];
        }
        if (std::abs(sum - 1.0) < convergence_crit) {
          converged = true;

          // Find the speed of sound based on converged solution
          for (size_t m = 0; m < numMatls; m++) {
            ICEMaterial* ice_matl = static_cast<ICEMaterial*>(
              d_materialManager->getMaterial("ICE", m));
            ice_matl->getEOS()->computePressEOS(rho_micro[m][c],
                                                gamma[m][c],
                                                cv[m][c],
                                                Temp[m][c],
                                                press_eos[m],
                                                dp_drho[m],
                                                dp_de[m]);

            tmp = dp_drho[m] +
                  dp_de[m] * press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]);
            speedSound_new[m][c] = std::sqrt(tmp);
          }
        }

        // Save iteration data for output in case of crash
        if (ds_EqPress.active()) {
          EqPress_dbg dbg;
          dbg.delPress   = delPress;
          dbg.press_new  = press_new[c];
          dbg.sumVolFrac = sum;
          dbg.count      = count;

          for (size_t m = 0; m < numMatls; m++) {
            EqPress_dbgMatl dmatl;
            dmatl.press_eos = press_eos[m];
            dmatl.volFrac   = vol_frac[m][c];
            dmatl.rhoMicro  = rho_micro[m][c];
            dmatl.rho_CC    = rho_CC[m][c];
            dmatl.temp_CC   = Temp[m][c];
            dmatl.mat       = m;
            dbg.matl.push_back(dmatl);
          }
          dbgEqPress.push_back(dbg);
        }
      } // end of converged

      test_max_iter = std::max(test_max_iter, count);

      //__________________________________
      //      BULLET PROOFING
      // ignore BP if a timestep restart has already been requested
      bool rts = new_dw->recomputeTimestep();

      std::string message;
      bool allTestsPassed{ true };
      if (test_max_iter == d_max_iter_equilibration && !rts) {
        allTestsPassed = false;
        message += "Max. iterations reached ";
      }

      for (size_t m = 0; m < numMatls; m++) {
        if ((vol_frac[m][c] > 0.0) || (vol_frac[m][c] < 1.0)) {
          message += " ( vol_frac[m][c] > 0.0 ) ||( vol_frac[m][c] < 1.0) ";
        }
      }

      if (std::abs(sum - 1.0) > convergence_crit && !rts) {
        allTestsPassed = false;
        message += " sum (volumeFractions) != 1 ";
      }

      if (press_new[c] < 0.0 && !rts) {
        allTestsPassed = false;
        message += " Computed pressure is < 0 ";
      }

      for (size_t m = 0; m < numMatls; m++) {
        if ((rho_micro[m][c] < 0.0 || vol_frac[m][c] < 0.0) && !rts) {
          allTestsPassed = false;
          message += " rho_micro < 0 || vol_frac < 0";
        }
      }

      if (allTestsPassed != true) { // throw an exception of there's a problem
        Point pt = patch->getCellPosition(c);
        std::ostringstream warn;
        warn << "\nICE::ComputeEquilibrationPressure: Cell " << c
             << " position: " << pt << ", L-" << L_indx << "\n"
             << message
             << "\nThis usually means that something much deeper has gone "
                "wrong with the simulation. "
             << "\nCompute equilibration pressure task is rarely the problem. "
             << "For more debugging information set the environmental "
                "variable:  \n"
             << "   SCI_DEBUG DBG_EqPress:+\n\n";

        warn << "INPUTS: \n";
        for (size_t m = 0; m < numMatls; m++) {
          warn << "\n matl: " << m << "\n"
               << "   rho_CC:     " << rho_CC[m][c] << "\n"
               << "   Temperature:   " << Temp[m][c] << "\n";
        }
        if (ds_EqPress.active()) {
          warn << "\nDetails on iterations " << std::endl;
          for (auto& d : dbgEqPress) {
            warn << "Iteration:   " << d.count
                 << "  press_new:   " << d.press_new
                 << "  sumVolFrac:  " << d.sumVolFrac
                 << "  delPress:    " << d.delPress << "\n";
            for (size_t m = 0; m < numMatls; m++) {
              warn << "  matl: " << d.matl[m].mat
                   << "  press_eos:  " << d.matl[m].press_eos
                   << "  volFrac:    " << d.matl[m].volFrac
                   << "  rhoMicro:   " << d.matl[m].rhoMicro
                   << "  rho_CC:     " << d.matl[m].rho_CC
                   << "  Temp:       " << d.matl[m].temp_CC << "\n";
            }
          }
        }
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
    } // end of cell interator

    cout_norm << "max. iterations in any cell " << test_max_iter << " on patch "
              << patch->getID() << endl;

    // carry rho_cc forward
    // MPMICE computes rho_CC_new
    // therefore need the machinery here
    for (size_t m = 0; m < numMatls; m++) {
      rho_CC_new[m].copyData(rho_CC[m]);
    }

    // - update Boundary conditions
    std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
      std::make_unique<CustomBCDriver::customBC_localVars>();

    preprocess_CustomBCs("EqPress",
                         old_dw,
                         new_dw,
                         d_ice_labels.get(),
                         patch,
                         999,
                         d_BC_globalVars.get(),
                         BC_localVars.get());

    setBC(press_new,
          rho_micro,
          placeHolder,
          d_surroundingMatl_indx,
          "rho_micro",
          "Pressure",
          patch,
          d_materialManager,
          0,
          new_dw,
          d_BC_globalVars.get(),
          BC_localVars.get(),
          isNotInitialTimestep);

    delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());

    // compute sp_vol_CC
    // - Set BCs on rhoMicro. using press_CC
    // - backout sp_vol_new
    for (size_t m = 0; m < numMatls; m++) {
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c      = *iter;
        sp_vol_new[m][c] = 1.0 / rho_micro[m][c];
      }

      int indx = d_materialManager->getMaterial("ICE", m)->getDWIndex();
      setSpecificVolBC(sp_vol_new[m],
                       "SpecificVol",
                       false,
                       rho_CC[m],
                       vol_frac[m],
                       patch,
                       d_materialManager,
                       indx);
    }

    //  compute f_theta
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      sumKappa[c] = 0.0;
      for (size_t m = 0; m < numMatls; m++) {
        kappa[m][c] =
          sp_vol_new[m][c] / (speedSound_new[m][c] * speedSound_new[m][c]);
        sumKappa[c] += vol_frac[m][c] * kappa[m][c];
      }
      for (size_t m = 0; m < numMatls; m++) {
        f_theta[m][c] = vol_frac[m][c] * kappa[m][c] / sumKappa[c];
      }
    }
  } // patch loop
}

/*!
 *  Function~  ICE::computeEquilPressure_1_matl--
 *  Purpose~   Simple EOS evaluation
 */
void
ICE::computeEquilPressure_1_matl(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 [[maybe_unused]] const MaterialSubset* matls,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimestep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computeEquilPressure_1_matl");

    CCVariable<double> vol_frac, sp_vol_new;
    CCVariable<double> speedSound, f_theta, kappa;
    CCVariable<double> press_eq, sumKappa, sum_imp_delP, rho_CC_new;
    constCCVariable<double> Temp, rho_CC, sp_vol_CC, cv, gamma;
    std::vector<CCVariable<double>> rho_micro(1);

    ICEMaterial* ice_matl =
      static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", 0));
    int indx = ice_matl->getDWIndex();

    old_dw->get(Temp, d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
    old_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
    old_dw->get(
      sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gn, 0);
    new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);
    new_dw->get(gamma, d_ice_labels->gammaLabel, indx, patch, d_gn, 0);

    new_dw->allocateTemporary(rho_micro[0], patch);

    new_dw->allocateAndPut(
      press_eq, d_ice_labels->press_equil_CCLabel, 0, patch);
    new_dw->allocateAndPut(sumKappa, d_ice_labels->sumKappaLabel, 0, patch);
    new_dw->allocateAndPut(
      sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch);
    new_dw->allocateAndPut(
      kappa, d_ice_labels->compressibilityLabel, indx, patch);
    new_dw->allocateAndPut(
      vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch);
    new_dw->allocateAndPut(
      sp_vol_new, d_ice_labels->specificVolume_CCLabel, indx, patch);
    new_dw->allocateAndPut(f_theta, d_ice_labels->f_theta_CCLabel, indx, patch);
    new_dw->allocateAndPut(
      speedSound, d_ice_labels->speedSound_CCLabel, indx, patch);
    new_dw->allocateAndPut(rho_CC_new, d_ice_labels->rho_CCLabel, indx, patch);

    sum_imp_delP.initialize(0.0);

    //  Main loop
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      vol_frac[c]     = 1.0;
      rho_micro[0][c] = rho_CC[c];
      rho_CC_new[c]   = rho_CC[c];
      sp_vol_new[c]   = 1.0 / rho_CC[c];

      double dp_drho, dp_de, c_2;

      // evaluate EOS
      ice_matl->getEOS()->computePressEOS(
        rho_micro[0][c], gamma[c], cv[c], Temp[c], press_eq[c], dp_drho, dp_de);

      c_2 = dp_drho + dp_de * press_eq[c] / (rho_micro[0][c] * rho_micro[0][c]);
      speedSound[c] = std::sqrt(c_2);

      //  compute f_theta
      kappa[c]    = sp_vol_new[c] / (speedSound[c] * speedSound[c]);
      sumKappa[c] = kappa[c];
      f_theta[c]  = 1.0;
    }

    // - apply Boundary conditions
    std::vector<constCCVariable<double>> placeHolder(0);
    std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
      std::make_unique<CustomBCDriver::customBC_localVars>();

    preprocess_CustomBCs("EqPress",
                         old_dw,
                         new_dw,
                         d_ice_labels.get(),
                         patch,
                         999,
                         d_BC_globalVars.get(),
                         BC_localVars.get());

    setBC(press_eq,
          rho_micro,
          placeHolder,
          d_surroundingMatl_indx,
          "rho_micro",
          "Pressure",
          patch,
          d_materialManager,
          0,
          new_dw,
          d_BC_globalVars.get(),
          BC_localVars.get(),
          isNotInitialTimestep);

    delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());

  } // patch loop
}

/*!
 *  Function~  ICE::computeTempFace--
 *  Purpose~   Compute the face centered Temperatures.  This is used by
 *             the HE combustion model
 */
template<class T>
void
ICE::computeTempFace(CellIterator it,
                     IntVector adj_offset,
                     constCCVariable<double>& rho_CC,
                     constCCVariable<double>& Temp_CC,
                     T& Temp_FC)
{
  for (; !it.done(); it++) {
    IntVector R = *it;
    IntVector L = R + adj_offset;

    double rho_FC = rho_CC[L] + rho_CC[R];
    ASSERT(rho_FC > 0.0);
    //__________________________________
    // interpolation to the face
    //  based on continuity of heat flux
    double term1 = (rho_CC[L] * Temp_CC[L] + rho_CC[R] * Temp_CC[R]) / (rho_FC);
    Temp_FC[R]   = term1;
  }
}

/*!
 *  Task:    ICE::computeTempFC
 *  Purpose: This function calculates the face centered temeratures on all faces
 *           in the domain.
 */
void
ICE::computeTempFC(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* /*matls*/,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computeTemp_FCVel");

    size_t numMatls = d_materialManager->getNumMaterials();

    // Compute the face centered Temperatures
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      constCCVariable<double> rho_CC, Temp_CC;

      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      if (ice_matl) {
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 1);
        old_dw->get(
          Temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gac, 1);
      } else {
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 1);
        new_dw->get(
          Temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gac, 1);
      }

      SFCXVariable<double> TempX_FC;
      SFCYVariable<double> TempY_FC;
      SFCZVariable<double> TempZ_FC;

      new_dw->allocateAndPut(
        TempX_FC, d_ice_labels->TempX_FCLabel, indx, patch);
      new_dw->allocateAndPut(
        TempY_FC, d_ice_labels->TempY_FCLabel, indx, patch);
      new_dw->allocateAndPut(
        TempZ_FC, d_ice_labels->TempZ_FCLabel, indx, patch);

      IntVector lowIndex(patch->getExtraSFCXLowIndex());
      TempX_FC.initialize(0.0, lowIndex, patch->getExtraSFCXHighIndex());
      TempY_FC.initialize(0.0, lowIndex, patch->getExtraSFCYHighIndex());
      TempZ_FC.initialize(0.0, lowIndex, patch->getExtraSFCZHighIndex());

      std::vector<IntVector> adj_offset(3);
      adj_offset[0] = IntVector(-1, 0, 0); // X faces
      adj_offset[1] = IntVector(0, -1, 0); // Y faces
      adj_offset[2] = IntVector(0, 0, -1); // Z faces

      CellIterator XFC_iterator = patch->getSFCXIterator();
      CellIterator YFC_iterator = patch->getSFCYIterator();
      CellIterator ZFC_iterator = patch->getSFCZIterator();

      int levelIndex = getLevel(patches)->getIndex();
      if (levelIndex > 0) { // Finer levels need to hit the ghost cells
        IntVector l, h;
        l            = patch->getExtraCellIterator().begin();
        h            = patch->getExtraCellIterator().end();
        XFC_iterator = CellIterator(l + IntVector(1, 0, 0), h);
        YFC_iterator = CellIterator(l + IntVector(0, 1, 0), h);
        ZFC_iterator = CellIterator(l + IntVector(0, 0, 1), h);
      }

      //  Compute the temperature on each face
      //  Currently on used by HEChemistry
      computeTempFace<SFCXVariable<double>>(
        XFC_iterator, adj_offset[0], rho_CC, Temp_CC, TempX_FC);

      computeTempFace<SFCYVariable<double>>(
        YFC_iterator, adj_offset[1], rho_CC, Temp_CC, TempY_FC);

      computeTempFace<SFCZVariable<double>>(
        ZFC_iterator, adj_offset[2], rho_CC, Temp_CC, TempZ_FC);

    } // matls loop
  }   // patch loop
}

/*!
 *  Function~  ICE::computeVelFace--
 *  Purpose~   Compute the face centered velocities minus the exchange
 *             contribution.
 */
template<class T>
void
ICE::computeVelFace(int dir,
                    CellIterator it,
                    IntVector adj_offset,
                    double dx,
                    double delT,
                    double gravity,
                    constCCVariable<double>& rho_CC,
                    constCCVariable<double>& sp_vol_CC,
                    constCCVariable<Vector>& vel_CC,
                    constCCVariable<double>& press_CC,
                    T& vel_FC,
                    T& grad_P_FC,
                    bool include_acc)
{
  double inv_dx = 1.0 / dx;

  double one_or_zero = 1.;
  if (!include_acc) {
    one_or_zero = 0.0;
  }

  for (; !it.done(); it++) {
    IntVector R = *it;
    IntVector L = R + adj_offset;

    double rho_FC = rho_CC[L] + rho_CC[R];
#if SCI_ASSERTION_LEVEL >= 2
    if (rho_FC <= 0.0) {
      std::cout << d_myworld->myRank() << " rho_fc <= 0: " << rho_FC
                << " with L= " << L << " (" << rho_CC[L] << ") R= " << R << " ("
                << rho_CC[R] << ")\n";
    }
#endif
    ASSERT(rho_FC > 0.0);

    //__________________________________
    // interpolation to the face
    double term1 =
      (rho_CC[L] * vel_CC[L][dir] + rho_CC[R] * vel_CC[R][dir]) / (rho_FC);
    //__________________________________
    // pressure term
    double sp_vol_brack =
      2. * (sp_vol_CC[L] * sp_vol_CC[R]) / (sp_vol_CC[L] + sp_vol_CC[R]);

    grad_P_FC[R] = (press_CC[R] - press_CC[L]) * inv_dx;
    double term2 = delT * sp_vol_brack * grad_P_FC[R];

    //__________________________________
    // gravity term
    double term3 = delT * gravity;

    vel_FC[R] = term1 - one_or_zero * term2 + one_or_zero * term3;
  }
}

/*!
 *  Task:    ICE::computeVel_FC
 *  Purpose: This function calculates the face centered velocities on all faces
 *           in the domain.
 */
void
ICE::computeVel_FC(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* /*matls*/,
                   DataWarehouse* old_dw,
                   DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computeVel_FCVel");

    Vector dx      = patch->dCell();
    Vector gravity = getGravity();

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel, getLevel(patches));

    constCCVariable<double> press_CC;
    new_dw->get(
      press_CC, d_ice_labels->press_equil_CCLabel, 0, patch, d_gac, 1);

    // Compute the face centered velocities
    size_t numMatls = d_materialManager->getNumMaterials();
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      int indx              = matl->getDWIndex();
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);

      constCCVariable<double> rho_CC, sp_vol_CC;
      constCCVariable<Vector> vel_CC;
      if (ice_matl) {
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 1);
        old_dw->get(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gac, 1);
      } else {
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 1);
        new_dw->get(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gac, 1);
      }
      new_dw->get(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gac, 1);

      SFCXVariable<double> uvel_FC, grad_P_XFC;
      SFCYVariable<double> vvel_FC, grad_P_YFC;
      SFCZVariable<double> wvel_FC, grad_P_ZFC;

      new_dw->allocateAndPut(uvel_FC, d_ice_labels->uvel_FCLabel, indx, patch);
      new_dw->allocateAndPut(vvel_FC, d_ice_labels->vvel_FCLabel, indx, patch);
      new_dw->allocateAndPut(wvel_FC, d_ice_labels->wvel_FCLabel, indx, patch);

      // debugging variables
      new_dw->allocateAndPut(
        grad_P_XFC, d_ice_labels->grad_P_XFCLabel, indx, patch);
      new_dw->allocateAndPut(
        grad_P_YFC, d_ice_labels->grad_P_YFCLabel, indx, patch);
      new_dw->allocateAndPut(
        grad_P_ZFC, d_ice_labels->grad_P_ZFCLabel, indx, patch);

      IntVector lowIndex(patch->getExtraSFCXLowIndex());
      uvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCXHighIndex());
      vvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCYHighIndex());
      wvel_FC.initialize(0.0, lowIndex, patch->getExtraSFCZHighIndex());

      grad_P_XFC.initialize(0.0);
      grad_P_YFC.initialize(0.0);
      grad_P_ZFC.initialize(0.0);

      std::vector<IntVector> adj_offset(3);
      adj_offset[0] = IntVector(-1, 0, 0); // X faces
      adj_offset[1] = IntVector(0, -1, 0); // Y faces
      adj_offset[2] = IntVector(0, 0, -1); // Z faces

      CellIterator XFC_iterator = patch->getSFCXIterator();
      CellIterator YFC_iterator = patch->getSFCYIterator();
      CellIterator ZFC_iterator = patch->getSFCZIterator();

      bool include_acc = true;
      if (mpm_matl && d_with_rigid_mpm) {
        include_acc = false;
      }

      //__________________________________
      //  Compute vel_FC for each face
      computeVelFace<SFCXVariable<double>>(0,
                                           XFC_iterator,
                                           adj_offset[0],
                                           dx[0],
                                           delT,
                                           gravity[0],
                                           rho_CC,
                                           sp_vol_CC,
                                           vel_CC,
                                           press_CC,
                                           uvel_FC,
                                           grad_P_XFC,
                                           include_acc);

      computeVelFace<SFCYVariable<double>>(1,
                                           YFC_iterator,
                                           adj_offset[1],
                                           dx[1],
                                           delT,
                                           gravity[1],
                                           rho_CC,
                                           sp_vol_CC,
                                           vel_CC,
                                           press_CC,
                                           vvel_FC,
                                           grad_P_YFC,
                                           include_acc);

      computeVelFace<SFCZVariable<double>>(2,
                                           ZFC_iterator,
                                           adj_offset[2],
                                           dx[2],
                                           delT,
                                           gravity[2],
                                           rho_CC,
                                           sp_vol_CC,
                                           vel_CC,
                                           press_CC,
                                           wvel_FC,
                                           grad_P_ZFC,
                                           include_acc);

    } // matls loop
  }   // patch loop
}

/*!
 *  Function~  ICE::updateVelFace--
 *  - tack on delP to the face centered velocity
 */
template<class T>
void
ICE::updateVelFace([[maybe_unused]] int dir,
                   CellIterator it,
                   IntVector adj_offset,
                   double dx,
                   double delT,
                   constCCVariable<double>& sp_vol_CC,
                   constCCVariable<double>& imp_delP,
                   T& vel_FC,
                   T& grad_dp_FC)
{
  double inv_dx = 1.0 / dx;

  for (; !it.done(); it++) {
    IntVector R = *it;
    IntVector L = R + adj_offset;

    //__________________________________
    // pressure term
    double sp_vol_brack =
      2. * (sp_vol_CC[L] * sp_vol_CC[R]) / (sp_vol_CC[L] + sp_vol_CC[R]);

    grad_dp_FC[R] = (imp_delP[R] - imp_delP[L]) * inv_dx;
    double term2  = delT * sp_vol_brack * grad_dp_FC[R];

    vel_FC[R] -= term2;
  }
}

/*!
 *  Task:  ICE::updateVel_FC--
 *  Purpose: Update the face centered velocities after the pressure solve
 */
void
ICE::updateVel_FC(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* /*matls*/,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw,
                  bool recursion)
{
  const Level* level = getLevel(patches);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::updateVel_FCVel");

    Vector dx = patch->dCell();

    // define parent data warehouse
    DataWarehouse* pNewDW;
    DataWarehouse* pOldDW;
    if (recursion) {
      pNewDW = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      pOldDW = new_dw->getOtherDataWarehouse(Task::ParentOldDW);
    } else {
      pNewDW = new_dw;
      pOldDW = old_dw;
    }

    delt_vartype delT;
    pOldDW->get(delT, d_ice_labels->delTLabel, level);

    constCCVariable<double> imp_delP;
    new_dw->get(imp_delP, d_ice_labels->imp_delPLabel, 0, patch, d_gac, 1);

    size_t numMatls = d_materialManager->getNumMaterials();
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      constCCVariable<double> sp_vol_CC;
      pNewDW->get(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gac, 1);

      constSFCXVariable<double> uvel_FC_old;
      constSFCYVariable<double> vvel_FC_old;
      constSFCZVariable<double> wvel_FC_old;

      old_dw->get(
        uvel_FC_old, d_ice_labels->uvel_FCLabel, indx, patch, d_gn, 0);
      old_dw->get(
        vvel_FC_old, d_ice_labels->vvel_FCLabel, indx, patch, d_gn, 0);
      old_dw->get(
        wvel_FC_old, d_ice_labels->wvel_FCLabel, indx, patch, d_gn, 0);

      SFCXVariable<double> uvel_FC, grad_dp_XFC;
      SFCYVariable<double> vvel_FC, grad_dp_YFC;
      SFCZVariable<double> wvel_FC, grad_dp_ZFC;

      new_dw->allocateAndPut(uvel_FC, d_ice_labels->uvel_FCLabel, indx, patch);
      new_dw->allocateAndPut(vvel_FC, d_ice_labels->vvel_FCLabel, indx, patch);
      new_dw->allocateAndPut(wvel_FC, d_ice_labels->wvel_FCLabel, indx, patch);

      new_dw->allocateAndPut(
        grad_dp_XFC, d_ice_labels->grad_dp_XFCLabel, indx, patch);
      new_dw->allocateAndPut(
        grad_dp_YFC, d_ice_labels->grad_dp_YFCLabel, indx, patch);
      new_dw->allocateAndPut(
        grad_dp_ZFC, d_ice_labels->grad_dp_ZFCLabel, indx, patch);

      uvel_FC.copy(uvel_FC_old);
      vvel_FC.copy(vvel_FC_old);
      wvel_FC.copy(wvel_FC_old);

      std::vector<IntVector> adj_offset(3);
      adj_offset[0] = IntVector(-1, 0, 0); // X faces
      adj_offset[1] = IntVector(0, -1, 0); // Y faces
      adj_offset[2] = IntVector(0, 0, -1); // Z faces

      CellIterator XFC_iterator = patch->getSFCXIterator();
      CellIterator YFC_iterator = patch->getSFCYIterator();
      CellIterator ZFC_iterator = patch->getSFCZIterator();

      updateVelFace<SFCXVariable<double>>(0,
                                          XFC_iterator,
                                          adj_offset[0],
                                          dx[0],
                                          delT,
                                          sp_vol_CC,
                                          imp_delP,
                                          uvel_FC,
                                          grad_dp_XFC);

      updateVelFace<SFCYVariable<double>>(1,
                                          YFC_iterator,
                                          adj_offset[1],
                                          dx[1],
                                          delT,
                                          sp_vol_CC,
                                          imp_delP,
                                          vvel_FC,
                                          grad_dp_YFC);

      updateVelFace<SFCZVariable<double>>(2,
                                          ZFC_iterator,
                                          adj_offset[2],
                                          dx[2],
                                          delT,
                                          sp_vol_CC,
                                          imp_delP,
                                          wvel_FC,
                                          grad_dp_ZFC);

    } // matls loop
  }   // patch loop
}

/*!
 *  Function~  ICE::computeDelPressAndUpdatePressCC--
 *  Purpose~   This function calculates the change in pressure explicitly.
 *  Note:      Units of delp_Dilatate and delP_MassX are [Pa]
 *  Reference:  Multimaterial Formalism eq. 1.5
 */
void
ICE::computeDelPressAndUpdatePressCC(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* /*matls*/,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimestep = (timeStep > 0);

  const Level* level = getLevel(patches);
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::computeDelPressAndUpdatePressCC");

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel, level);

    Vector dx      = patch->dCell();
    double inv_vol = 1.0 / (dx.x() * dx.y() * dx.z());

    std::unique_ptr<Advector> advector =
      d_advector->clone(new_dw, patch, isRegridTimestep());

    CCVariable<double> q_advected;
    CCVariable<double> delP_Dilatate;
    CCVariable<double> delP_MassX;
    CCVariable<double> sum_rho_CC;
    CCVariable<double> press_CC;
    CCVariable<double> term1, term2;
    constCCVariable<double> sumKappa, press_equil;

    size_t numMatls = d_materialManager->getNumMaterials();
    std::vector<CCVariable<double>> placeHolder(0);
    std::vector<constCCVariable<double>> sp_vol_CC(numMatls);

    const IntVector gc(1, 1, 1);

    new_dw->get(sumKappa, d_ice_labels->sumKappaLabel, 0, patch, d_gn, 0);
    new_dw->get(
      press_equil, d_ice_labels->press_equil_CCLabel, 0, patch, d_gn, 0);

    new_dw->allocateAndPut(press_CC, d_ice_labels->press_CCLabel, 0, patch);
    new_dw->allocateAndPut(
      delP_Dilatate, d_ice_labels->delP_DilatateLabel, 0, patch);
    new_dw->allocateAndPut(delP_MassX, d_ice_labels->delP_MassXLabel, 0, patch);
    new_dw->allocateAndPut(term2, d_ice_labels->term2Label, 0, patch);
    new_dw->allocateAndPut(sum_rho_CC, d_ice_labels->sum_rho_CCLabel, 0, patch);

    new_dw->allocateTemporary(q_advected, patch);
    new_dw->allocateTemporary(term1, patch);

    term1.initialize(0.);
    term2.initialize(0.);
    sum_rho_CC.initialize(0.0);
    delP_Dilatate.initialize(0.0);
    delP_MassX.initialize(0.0);

    for (size_t m = 0; m < numMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      constCCVariable<double> speedSound;
      constCCVariable<double> vol_frac;
      constCCVariable<double> rho_CC;
      constSFCXVariable<double> uvel_FC;
      constSFCYVariable<double> vvel_FC;
      constSFCZVariable<double> wvel_FC;

      new_dw->get(uvel_FC, d_ice_labels->uvel_FCMELabel, indx, patch, d_gac, 2);
      new_dw->get(vvel_FC, d_ice_labels->vvel_FCMELabel, indx, patch, d_gac, 2);
      new_dw->get(wvel_FC, d_ice_labels->wvel_FCMELabel, indx, patch, d_gac, 2);
      new_dw->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, d_gac, 2);
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(sp_vol_CC[m],
                  d_ice_labels->specificVolume_CCLabel,
                  indx,
                  patch,
                  d_gn,
                  0);
      new_dw->get(
        speedSound, d_ice_labels->speedSound_CCLabel, indx, patch, d_gn, 0);

      SFCXVariable<double> vol_fracX_FC;
      SFCYVariable<double> vol_fracY_FC;
      SFCZVariable<double> vol_fracZ_FC;

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

      // Advection preprocessing
      // - divide vol_frac_cc/vol
      bool bulletProof_test = true;
      std::unique_ptr<advectVarBasket> varBasket =
        std::make_unique<advectVarBasket>();
      advector->inFluxOutFluxVolume(uvel_FC,
                                    vvel_FC,
                                    wvel_FC,
                                    delT,
                                    patch,
                                    indx,
                                    bulletProof_test,
                                    new_dw,
                                    varBasket.get());
      //__________________________________
      //   advect vol_frac
      // common variables that get passed into the advection operators
      varBasket->doRefluxing = false; // don't need to reflux here
      advector->advectQ(vol_frac,
                        patch,
                        q_advected,
                        varBasket.get(),
                        vol_fracX_FC,
                        vol_fracY_FC,
                        vol_fracZ_FC,
                        new_dw);

      varBasket.reset(nullptr);

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        term2[c] -= q_advected[c];
      }

      //   term1 contribution from models
      if (d_models.size() > 0) {
        constCCVariable<double> modelMass_src, modelVol_src;
        new_dw->get(modelMass_src,
                    d_ice_labels->modelMass_srcLabel,
                    indx,
                    patch,
                    d_gn,
                    0);

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          term1[c] += modelMass_src[c] * (sp_vol_CC[m][c] * inv_vol);
        }
      }

      //  compute sum_rho_CC used by press_FC
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        sum_rho_CC[c] += rho_CC[c];
      }
    } // matl loop

    advector.reset(nullptr);

    //__________________________________
    //  add delP to press_equil
    //  AMR:  hit the extra cells, BC aren't set an you need a valid pressure
    //  there
    // THIS COULD BE TROUBLE
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      press_CC[c] = press_equil[c];
    }

    for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c         = *iter;
      double inv_sumKappa = 1.0 / sumKappa[c];
      delP_MassX[c]       = term1[c] * inv_sumKappa;
      delP_Dilatate[c]    = -term2[c] * inv_sumKappa;
      press_CC[c]         = press_equil[c] + delP_MassX[c] + delP_Dilatate[c];
      press_CC[c]         = std::max(1.0e-12, press_CC[c]); // CLAMP
    }

    //__________________________________
    //  set boundary conditions
    std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
      std::make_unique<CustomBCDriver::customBC_localVars>();
    preprocess_CustomBCs("update_press_CC",
                         old_dw,
                         new_dw,
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
#if SET_CFI_BC
    set_CFI_BC<double>(press_CC, patch);
#endif
    delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());
  } // patch loop
}

/*!
 *  Function~  ICE::computePressFC--
 *  Purpose~ This function calculates the face centered pressure on each of the
 *  cell faces for every cell in the computational domain and a single
 *  layer of ghost cells.
 */
template<class T>
void
ICE::computePressFace(CellIterator iter,
                      IntVector adj_offset,
                      constCCVariable<double>& sum_rho,
                      constCCVariable<double>& press_CC,
                      T& press_FC)
{
  for (; !iter.done(); iter++) {
    IntVector R = *iter;
    IntVector L = R + adj_offset;

    press_FC[R] = (press_CC[R] * sum_rho[L] + press_CC[L] * sum_rho[R]) /
                  (sum_rho[R] + sum_rho[L]);
  }
}

/*!
 * Task:  ICE::computePressFC--
 * Purpose:  This function calculates the face centered pressure on each of the
 *     cell faces for every cell in the computational domain and a single
 *     layer of ghost cells.
 */
void
ICE::computePressFC(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* /*matls*/,
                    DataWarehouse*,
                    DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computePressFC");

    constCCVariable<double> press_CC;
    constCCVariable<double> sum_rho_CC;
    new_dw->get(press_CC, d_ice_labels->press_CCLabel, 0, patch, d_gac, 1);
    new_dw->get(sum_rho_CC, d_ice_labels->sum_rho_CCLabel, 0, patch, d_gac, 1);

    SFCXVariable<double> pressX_FC;
    SFCYVariable<double> pressY_FC;
    SFCZVariable<double> pressZ_FC;
    new_dw->allocateAndPut(pressX_FC, d_ice_labels->pressX_FCLabel, 0, patch);
    new_dw->allocateAndPut(pressY_FC, d_ice_labels->pressY_FCLabel, 0, patch);
    new_dw->allocateAndPut(pressZ_FC, d_ice_labels->pressZ_FCLabel, 0, patch);

    std::vector<IntVector> adj_offset(3);
    adj_offset[0] = IntVector(-1, 0, 0); // X faces
    adj_offset[1] = IntVector(0, -1, 0); // Y faces
    adj_offset[2] = IntVector(0, 0, -1); // Z faces

    //__________________________________
    //  For each face compute the pressure
    computePressFace<SFCXVariable<double>>(
      patch->getSFCXIterator(), adj_offset[0], sum_rho_CC, press_CC, pressX_FC);

    computePressFace<SFCYVariable<double>>(
      patch->getSFCYIterator(), adj_offset[1], sum_rho_CC, press_CC, pressY_FC);

    computePressFace<SFCZVariable<double>>(
      patch->getSFCZIterator(), adj_offset[2], sum_rho_CC, press_CC, pressZ_FC);

  } // patch loop
}

/*!
 *  Function~  ICE::zeroModelSources
 *  Purpose~   This function initializes the mass exchange quantities to zero.
 *             These quantities are subsequently modified by the models
 */
void
ICE::zeroModelSources(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* /*old_dw*/,
                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::zeroModelSources");

    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      CCVariable<double> mass_src, energy_src, vol_src;
      CCVariable<Vector> mom_src;

      new_dw->allocateAndPut(
        mass_src, d_ice_labels->modelMass_srcLabel, matl, patch);
      new_dw->allocateAndPut(
        energy_src, d_ice_labels->modelEng_srcLabel, matl, patch);
      new_dw->allocateAndPut(
        mom_src, d_ice_labels->modelMom_srcLabel, matl, patch);
      new_dw->allocateAndPut(
        vol_src, d_ice_labels->modelVol_srcLabel, matl, patch);

      energy_src.initialize(0.0);
      mass_src.initialize(0.0);
      vol_src.initialize(0.0);
      mom_src.initialize(Vector(0.0, 0.0, 0.0));
    }

    // Model with transported variables.
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model) {
        for (auto& tvar : fb_model->d_transVars) {
          for (int m = 0; m < tvar->matls->size(); m++) {
            int matl = tvar->matls->get(m);

            if (tvar->src) {
              CCVariable<double> model_src;
              new_dw->allocateAndPut(model_src, tvar->src, matl, patch);
              model_src.initialize(0.0);
            }
          } // matl loop
        }   // tvar loop
      }     // fb_model if
    }       // model loop
  }         // patches loop
}

/*!
 *  Function~  ICE::updateVolumeFraction
 *  Purpose~   Update the volume fraction to reflect the mass exchange done
 *             by models
 */
void
ICE::updateVolumeFraction(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* /*old_dw*/,
                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "Doing updateVolumeFraction");

    size_t numALLMatls = d_materialManager->getNumMaterials();

    CCVariable<double> sumKappa;
    new_dw->getModifiable(sumKappa, d_ice_labels->sumKappaLabel, 0, patch);

    std::vector<CCVariable<double>> vol_frac(numALLMatls);
    std::vector<CCVariable<double>> f_theta(numALLMatls);

    std::vector<constCCVariable<double>> rho_CC(numALLMatls);
    std::vector<constCCVariable<double>> sp_vol(numALLMatls);
    std::vector<constCCVariable<double>> modVolSrc(numALLMatls);
    std::vector<constCCVariable<double>> kappa(numALLMatls);

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    for (int m = 0; m < matls->size(); m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      new_dw->getModifiable(
        vol_frac[m], d_ice_labels->vol_frac_CCLabel, indx, patch);
      new_dw->getModifiable(
        f_theta[m], d_ice_labels->f_theta_CCLabel, indx, patch);
      new_dw->get(rho_CC[m], d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        sp_vol[m], d_ice_labels->specificVolume_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        modVolSrc[m], d_ice_labels->modelVol_srcLabel, indx, patch, d_gn, 0);
      new_dw->get(
        kappa[m], d_ice_labels->compressibilityLabel, indx, patch, d_gn, 0);
    }

    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      IntVector c      = *iter;
      double total_vol = 0.;
      for (int m = 0; m < matls->size(); m++) {
        total_vol += (rho_CC[m][c] * vol) * sp_vol[m][c];
      }

      double sumKappa_tmp = 0.0;
      for (int m = 0; m < matls->size(); m++) {
        double new_vol = vol_frac[m][c] * total_vol + modVolSrc[m][c];
        vol_frac[m][c] = std::max(new_vol / total_vol, 1.e-100);
        sumKappa_tmp += vol_frac[m][c] * kappa[m][c];
      }
      sumKappa[c] = sumKappa_tmp;

      for (int m = 0; m < matls->size(); m++) {
        f_theta[m][c] = vol_frac[m][c] * kappa[m][c] / sumKappa[c];
      }
    }
  }
}

/*!
 *    See comments in ICE::VelTau_CC()
 */
void
ICE::computeVelTau_CCFace(const Patch* patch,
                          const Patch::FaceType face,
                          constCCVariable<Vector>& vel_CC,
                          CCVariable<Vector>& velTau_CC)
{
  CellIterator iterLimits =
    patch->getFaceIterator(face, Patch::ExtraMinusEdgeCells);
  IntVector oneCell = patch->faceDirection(face);

  for (auto iter = iterLimits; !iter.done(); iter++) {
    IntVector c  = *iter;       // extra cell index
    IntVector in = c - oneCell; // interior cell index

    velTau_CC[c] = 2. * vel_CC[c] - vel_CC[in];
  }
}

/*!
 *   Modify the vel_CC in the extra cells so that it behaves
 *      vel_FC[FC] = (vel_CC(c) + vel_CC(ec) )/2        (1)
 *
 *   Note that at the edge of the domain we assume that vel_FC = vel_CC[ec]
 *   so (1) becomes:
 *      vel_CC[ec]    = (vel_CC(c) + velTau_CC(ec) )/2
 *             or
 *      velTau_CC[ec] = (2 * vel_CC(ec) - velTau_CC(c) );
 *
 *   You need this so the viscous shear stress terms tau_xy = tau_yx.
 *
 *                |           |
 *     ___________|___________|_______________
 *                |           |
 *          o     |     o     |      o         Vel_CC
 *                |     c     |
 *                |           |
 *                |           |
 *     ===========|====FC=====|==============      Edge of computational Domain
 *                |           |
 *                |           |
 *          *     |     o     |      o          Vel_CC in extraCell
 *                |    ec     |
 *                |           |
 *     ___________|___________|_______________
 *
 *    A fundamental assumption is that the boundary conditions
 *   have been vel_CC[ec] in the old_dw.
 */
void
ICE::VelTau_CC(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* /*matls*/,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, m_ice_tasks, "Doing VelTau_CC");

    //__________________________________
    //  ICE matl loop
    size_t numMatls = d_materialManager->getNumMaterials("ICE");

    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      constCCVariable<Vector> vel_CC;
      CCVariable<Vector> velTau_CC;

      old_dw->get(vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gn, 0);
      new_dw->allocateAndPut(
        velTau_CC, d_ice_labels->velTau_CCLabel, indx, patch);

      velTau_CC.copyData(vel_CC); // copy interior values over

      if (patch->hasBoundaryFaces()) {

        // Iterate over the faces encompassing the domain
        vector<Patch::FaceType> bf;
        patch->getBoundaryFaces(bf);

        for (const auto& face : bf) {
          //__________________________________
          //           X faces
          if (face == Patch::xminus || face == Patch::xplus) {
            computeVelTau_CCFace(patch, face, vel_CC, velTau_CC);
            continue;
          }

          //__________________________________
          //           Y faces
          if (face == Patch::yminus || face == Patch::yplus) {
            computeVelTau_CCFace(patch, face, vel_CC, velTau_CC);
            continue;
          }

          //__________________________________
          //           Z faces
          if (face == Patch::zminus || face == Patch::zplus) {
            computeVelTau_CCFace(patch, face, vel_CC, velTau_CC);
            continue;
          }

        } // face loop
      }   // has boundary face
    }     // matl loop
  }       // patch loop
}

/*!
 *  Task:     ICE::viscousShearStress
 *  Purpose:  Compute the viscous shear stress terms on each cell faces
 */
void
ICE::viscousShearStress(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* /*matls*/,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, m_ice_tasks, "ICE::viscousShearStress");

    IntVector right, left, top, bottom, front, back;

    Vector dx    = patch->dCell();
    double areaX = dx.y() * dx.z();
    double areaY = dx.x() * dx.z();
    double areaZ = dx.x() * dx.y();

    //__________________________________
    //  Matl loop
    size_t numMatls = d_materialManager->getNumMaterials("ICE");

    for (size_t m = 0; m < numMatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      CCVariable<Vector> viscous_src;
      new_dw->allocateAndPut(
        viscous_src, d_ice_labels->viscous_src_CCLabel, indx, patch);
      viscous_src.initialize(Vector(0., 0., 0.));

      //__________________________________
      // Compute Viscous diffusion for this matl
      if (d_viscousFlow) {
        constCCVariable<double> vol_frac;
        constCCVariable<double> rho_CC;

        new_dw->get(
          vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, d_gac, 2);
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 2);

        SFCXVariable<Vector> tau_X_FC;
        SFCYVariable<Vector> tau_Y_FC;
        SFCZVariable<Vector> tau_Z_FC;

        SFCXVariable<Vector> Ttau_X_FC;
        SFCYVariable<Vector> Ttau_Y_FC;
        SFCZVariable<Vector> Ttau_Z_FC;

        new_dw->allocateAndPut(
          tau_X_FC, d_ice_labels->tau_X_FCLabel, indx, patch);
        new_dw->allocateAndPut(
          tau_Y_FC, d_ice_labels->tau_Y_FCLabel, indx, patch);
        new_dw->allocateAndPut(
          tau_Z_FC, d_ice_labels->tau_Z_FCLabel, indx, patch);

        tau_X_FC.initialize(Vector(0.0)); // DEFAULT VALUE
        tau_Y_FC.initialize(Vector(0.0));
        tau_Z_FC.initialize(Vector(0.0));

        //__________________________________
        //  compute the shear stress terms
        double viscosity_test = ice_matl->getViscosity();

        if (viscosity_test != 0.0) {
          CCVariable<double> tot_viscosity; // total viscosity
          CCVariable<double> viscosity;     // total *OR* molecular viscosity;
          constCCVariable<double> molecularVis; // molecular viscosity
          constCCVariable<Vector> velTau_CC;

          new_dw->get(
            molecularVis, d_ice_labels->viscosityLabel, indx, patch, d_gac, 2);
          new_dw->get(
            velTau_CC, d_ice_labels->velTau_CCLabel, indx, patch, d_gac, 2);

          // don't alter the original value
          new_dw->allocateTemporary(tot_viscosity, patch, d_gac, 2);
          new_dw->allocateTemporary(viscosity, patch, d_gac, 2);
          viscosity.copyData(molecularVis);

          // Use temporary arrays to eliminate the communication of shear stress
          // components across the network.  Normally you would compute them in
          // a separate task and then require them with ghostCells to compute
          // the divergence.
          SFCXVariable<Vector> Ttau_X_FC;
          SFCYVariable<Vector> Ttau_Y_FC;
          SFCZVariable<Vector> Ttau_Z_FC;

          new_dw->allocateTemporary(Ttau_X_FC, patch, d_gac, 1);
          new_dw->allocateTemporary(Ttau_Y_FC, patch, d_gac, 1);
          new_dw->allocateTemporary(Ttau_Z_FC, patch, d_gac, 1);

          Vector evilNum(-9e30);
          Ttau_X_FC.initialize(evilNum);
          Ttau_Y_FC.initialize(evilNum);
          Ttau_Z_FC.initialize(evilNum);

          // turbulence model
          if (d_turbulence) {
            d_turbulence->callTurb(new_dw,
                                   patch,
                                   velTau_CC,
                                   rho_CC,
                                   vol_frac,
                                   indx,
                                   d_ice_labels.get(),
                                   d_materialManager,
                                   molecularVis,
                                   tot_viscosity);
            viscosity.copyData(tot_viscosity);
          }

          computeTauComponents(patch,
                               vol_frac,
                               velTau_CC,
                               viscosity,
                               Ttau_X_FC,
                               Ttau_Y_FC,
                               Ttau_Z_FC);

          // wall model
          if (d_WallShearStressModel) {
            d_WallShearStressModel->computeWallShearStresses(old_dw,
                                                             new_dw,
                                                             patch,
                                                             indx,
                                                             vol_frac,
                                                             velTau_CC,
                                                             molecularVis,
                                                             rho_CC,
                                                             Ttau_X_FC,
                                                             Ttau_Y_FC,
                                                             Ttau_Z_FC);
          }

          for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            right       = c + IntVector(1, 0, 0);
            left        = c + IntVector(0, 0, 0);
            top         = c + IntVector(0, 1, 0);
            bottom      = c + IntVector(0, 0, 0);
            front       = c + IntVector(0, 0, 1);
            back        = c + IntVector(0, 0, 0);

            viscous_src[c].x(
              (Ttau_X_FC[right].x() - Ttau_X_FC[left].x()) * areaX +
              (Ttau_Y_FC[top].x() - Ttau_Y_FC[bottom].x()) * areaY +
              (Ttau_Z_FC[front].x() - Ttau_Z_FC[back].x()) * areaZ);

            viscous_src[c].y(
              (Ttau_X_FC[right].y() - Ttau_X_FC[left].y()) * areaX +
              (Ttau_Y_FC[top].y() - Ttau_Y_FC[bottom].y()) * areaY +
              (Ttau_Z_FC[front].y() - Ttau_Z_FC[back].y()) * areaZ);

            viscous_src[c].z(
              (Ttau_X_FC[right].z() - Ttau_X_FC[left].z()) * areaX +
              (Ttau_Y_FC[top].z() - Ttau_Y_FC[bottom].z()) * areaY +
              (Ttau_Z_FC[front].z() - Ttau_Z_FC[back].z()) * areaZ);
          }
          // copy the temporary data
          tau_X_FC.copyPatch(
            Ttau_X_FC, tau_X_FC.getLowIndex(), tau_X_FC.getHighIndex());
          tau_Y_FC.copyPatch(
            Ttau_Y_FC, tau_Y_FC.getLowIndex(), tau_Y_FC.getHighIndex());
          tau_Z_FC.copyPatch(
            Ttau_Z_FC, tau_Z_FC.getLowIndex(), tau_Z_FC.getHighIndex());

        } // hasViscosity
      }   // ice_matl
    }     // matl loop
  }       // patch loop
}

/*!
 *  Function~  ICE::accumulateMomentumSourceSinks--
 *  Purpose~   This function accumulates all of the sources/sinks of momentum
 */
void
ICE::accumulateMomentumSourceSinks(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* /*matls*/,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::accumulateMomentumSourceSinks");

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel, level);

    Vector dx      = patch->dCell();
    Vector gravity = getGravity();
    double vol     = dx.x() * dx.y() * dx.z();

    double areaX = dx.y() * dx.z();
    double areaY = dx.x() * dx.z();
    double areaZ = dx.x() * dx.y();

    constSFCXVariable<double> pressX_FC;
    constSFCYVariable<double> pressY_FC;
    constSFCZVariable<double> pressZ_FC;

    new_dw->get(pressX_FC, d_ice_labels->pressX_FCLabel, 0, patch, d_gac, 1);
    new_dw->get(pressY_FC, d_ice_labels->pressY_FCLabel, 0, patch, d_gac, 1);
    new_dw->get(pressZ_FC, d_ice_labels->pressZ_FCLabel, 0, patch, d_gac, 1);

    IntVector right, left, top, bottom, front, back;

    //  Matl loop
    size_t numMatls = d_materialManager->getNumMaterials();
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      int indx = matl->getDWIndex();

      constCCVariable<double> vol_frac;
      constCCVariable<double> rho_CC;
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, d_gn, 0);

      CCVariable<Vector> mom_source;
      new_dw->allocateAndPut(
        mom_source, d_ice_labels->mom_source_CCLabel, indx, patch);
      mom_source.initialize(Vector(0., 0., 0.));

      //  accumulate sources MPM and ICE matls
      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;

        right  = c + IntVector(1, 0, 0);
        left   = c;
        top    = c + IntVector(0, 1, 0);
        bottom = c;
        front  = c + IntVector(0, 0, 1);
        back   = c;

        double press_src_X = (pressX_FC[right] - pressX_FC[left]) * vol_frac[c];
        double press_src_Y = (pressY_FC[top] - pressY_FC[bottom]) * vol_frac[c];
        double press_src_Z = (pressZ_FC[front] - pressZ_FC[back]) * vol_frac[c];

        mom_source[c].x(-press_src_X * areaX);
        mom_source[c].y(-press_src_Y * areaY);
        mom_source[c].z(-press_src_Z * areaZ);
      }

      //  Add user defined pressure gradient
      Vector evilNumV(d_EVIL_NUM);

      if (d_fixedPressGrad.length() != evilNumV.length()) {

        Vector oneZero;
        oneZero.x((d_fixedPressGrad.x() == d_EVIL_NUM) ? 1.0 : 0.0);
        oneZero.y((d_fixedPressGrad.y() == d_EVIL_NUM) ? 1.0 : 0.0);
        oneZero.z((d_fixedPressGrad.z() == d_EVIL_NUM) ? 1.0 : 0.0);

        double src_X    = d_fixedPressGrad.x() * dx.x() * areaX;
        double src_Y    = d_fixedPressGrad.y() * dx.y() * areaY;
        double src_Z    = d_fixedPressGrad.z() * dx.z() * areaZ;
        Vector pressSrc = Vector(src_X, src_Y, src_Z);
        Vector one(1.0, 1.0, 1.0);

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;

          Vector fixedPressSrc = pressSrc * vol_frac[c];
          mom_source[c] = mom_source[c] + (one - oneZero) * fixedPressSrc;
        }
      }

      //  ICE Matls
      if (ice_matl) {
        constCCVariable<Vector> viscous_src;
        new_dw->get(
          viscous_src, d_ice_labels->viscous_src_CCLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          double mass = rho_CC[c] * vol;

          mom_source[c] = (mom_source[c] + viscous_src[c] + mass * gravity);
        }
      }

      //  All Matls
      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        mom_source[*iter] *= delT;
      }

    } // matls loop
  }   // patches
}

/*!
 *  Function~  ICE::accumulateEnergySourceSinks--
 *  Purpose~   This function accumulates all of the sources/sinks of energy
 *             Currently the kinetic energy isn't included.
 */
void
ICE::accumulateEnergySourceSinks(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset* /*matls*/,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, d_ice_labels->simulationTimeLabel);
  double simTime = simTimeVar;

  delt_vartype delT;
  old_dw->get(delT, d_ice_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::accumulateEnergySourceSinks");

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    constCCVariable<double> sp_vol_CC;
    constCCVariable<double> kappa;
    constCCVariable<double> vol_frac;
    constCCVariable<double> press_CC;
    constCCVariable<double> delP_Dilatate;
    constCCVariable<double> matl_press;
    constCCVariable<double> rho_CC;
    constCCVariable<double> TMV_CC;

    new_dw->get(press_CC, d_ice_labels->press_CCLabel, 0, patch, d_gn, 0);
    new_dw->get(
      delP_Dilatate, d_ice_labels->delP_DilatateLabel, 0, patch, d_gn, 0);

    if (d_with_mpm) {
      new_dw->get(TMV_CC, d_ice_labels->TMV_CCLabel, 0, patch, d_gn, 0);
    } else {
      CCVariable<double> TMV_create;
      new_dw->allocateTemporary(TMV_create, patch);
      TMV_create.initialize(vol);
      TMV_CC = TMV_create; // reference created data
    }

    size_t numMatls = d_materialManager->getNumMaterials();
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      int indx = matl->getDWIndex();

      new_dw->get(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gac, 1);
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gac, 1);
      new_dw->get(
        kappa, d_ice_labels->compressibilityLabel, indx, patch, d_gn, 0);
      new_dw->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, d_gac, 1);

      CCVariable<double> int_eng_source;
      CCVariable<double> heatCond_src;
      new_dw->allocateAndPut(
        int_eng_source, d_ice_labels->int_eng_source_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        heatCond_src, d_ice_labels->heatCond_src_CCLabel, indx, patch);
      int_eng_source.initialize(0.0);
      heatCond_src.initialize(0.0);

      //__________________________________
      //  Source due to conduction ICE only
      if (ice_matl) {
        double thermalCond_test = ice_matl->getThermalConductivity();
        if (thermalCond_test != 0.0) {
          constCCVariable<double> Temp_CC;
          constCCVariable<double> thermalCond;
          new_dw->get(
            thermalCond, d_ice_labels->thermalCondLabel, indx, patch, d_gac, 1);
          old_dw->get(
            Temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gac, 1);

          bool use_vol_frac = true; // include vol_frac in diffusion calc.
          scalarDiffusionOperator(new_dw,
                                  patch,
                                  use_vol_frac,
                                  Temp_CC,
                                  vol_frac,
                                  heatCond_src,
                                  thermalCond,
                                  delT);
        }
      }

      //__________________________________
      //   Compute source from volume dilatation
      //   Exclude contribution from delP_MassX
      if (ice_matl->getIncludeFlowWork()) {
        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          double A    = TMV_CC[c] * vol_frac[c] * kappa[c] * press_CC[c];
          int_eng_source[c] += A * delP_Dilatate[c] + heatCond_src[c];
        }
      }

      //__________________________________
      //  User specified source/sink
      if (d_add_heat && simTime >= d_add_heat_t_start &&
          simTime <= d_add_heat_t_final) {

        for (size_t i = 0; i < d_add_heat_matls.size(); i++) {
          if (static_cast<int>(m) == d_add_heat_matls[i]) {
            for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
              IntVector c = *iter;
              if (vol_frac[c] > 0.001) {
                int_eng_source[c] +=
                  d_add_heat_coeff[i] * delT * rho_CC[c] * vol;
              }
            } // iter loop
          }   // if right matl
        }
      } // if add heat

    } // matl loop
  }   // patch loop
}

/*!
 *  Function~  ICE::computeLagrangianValues--
 *  Computes lagrangian mass momentum and energy
 *  Note:    Only loop over ICE materials, mom_L, massL and int_eng_L
 *           for MPM is computed in computeLagrangianValuesMPM()
 */
void
ICE::computeLagrangianValues(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* /*matls*/,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::computeLagrangianValues");

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    //__________________________________
    //  Compute the Lagrangian quantities
    size_t numALLMatls = d_materialManager->getNumMaterials();
    for (size_t m = 0; m < numALLMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      CCVariable<Vector> mom_L;
      CCVariable<double> int_eng_L;
      CCVariable<double> mass_L;
      double tiny_rho = 1.e-12;

      if (ice_matl) { //  I C E
        constCCVariable<double> rho_CC, temp_CC, cv, int_eng_source;
        constCCVariable<Vector> vel_CC, mom_source, mom_comb;
        tiny_rho = ice_matl->getTinyRho();

        new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);
        new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
        old_dw->get(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gn, 0);
        old_dw->get(
          temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(
          mom_source, d_ice_labels->mom_source_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(int_eng_source,
                    d_ice_labels->int_eng_source_CCLabel,
                    indx,
                    patch,
                    d_gn,
                    0);

        new_dw->allocateAndPut(mom_L, d_ice_labels->mom_L_CCLabel, indx, patch);
        new_dw->allocateAndPut(
          int_eng_L, d_ice_labels->int_eng_L_CCLabel, indx, patch);
        new_dw->allocateAndPut(
          mass_L, d_ice_labels->mass_L_CCLabel, indx, patch);

        //  NO mass exchange
        if (d_models.size() == 0) {
          for (auto iter = patch->getExtraCellIterator(); !iter.done();
               iter++) {
            IntVector c  = *iter;
            double mass  = rho_CC[c] * vol;
            mass_L[c]    = mass;
            mom_L[c]     = vel_CC[c] * mass + mom_source[c];
            int_eng_L[c] = mass * cv[c] * temp_CC[c] + int_eng_source[c];
          }
        }

        //__________________________________
        //      MODEL - B A S E D   E X C H A N G E
        //  WITH "model-based" mass exchange
        // Note that the mass exchange can't completely
        // eliminate all the mass, momentum and internal E
        // If it does then we'll get erroneous vel, and temps
        // after advection.  Thus there is always a mininum amount
        if (d_models.size() > 0) {
          constCCVariable<double> modelMass_src;
          constCCVariable<double> modelEng_src;
          constCCVariable<Vector> modelMom_src;
          new_dw->get(modelMass_src,
                      d_ice_labels->modelMass_srcLabel,
                      indx,
                      patch,
                      d_gn,
                      0);
          new_dw->get(modelMom_src,
                      d_ice_labels->modelMom_srcLabel,
                      indx,
                      patch,
                      d_gn,
                      0);
          new_dw->get(modelEng_src,
                      d_ice_labels->modelEng_srcLabel,
                      indx,
                      patch,
                      d_gn,
                      0);

          double massGain = 0.;
          for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
            IntVector c = *iter;
            massGain += modelMass_src[c];
          }

          for (auto iter = patch->getExtraCellIterator(); !iter.done();
               iter++) {
            IntVector c = *iter;

            //  must have a minimum mass
            double mass     = rho_CC[c] * vol;
            double min_mass = tiny_rho * vol;

            mass_L[c] = std::max((mass + modelMass_src[c]), min_mass);

            //  must have a minimum momentum
            for (int dir = 0; dir < 3; dir++) { // loop over all three directons
              double min_mom_L = vel_CC[c][dir] * min_mass;
              double mom_L_tmp = vel_CC[c][dir] * mass + modelMom_src[c][dir];

              // Preserve the original sign on momemtum
              // Use d_SMALL_NUMs to avoid nans when mom_L_temp = 0.0
              double plus_minus_one =
                (mom_L_tmp + d_SMALL_NUM) / (fabs(mom_L_tmp + d_SMALL_NUM));

              mom_L[c][dir] =
                mom_source[c][dir] +
                plus_minus_one * std::max(fabs(mom_L_tmp), min_mom_L);
            }

            // must have a minimum int_eng
            double min_int_eng = min_mass * cv[c] * temp_CC[c];
            double int_eng_tmp = mass * cv[c] * temp_CC[c];

            //  Glossary:
            //  int_eng_tmp    = the amount of internal energy for this
            //                   matl in this cell coming into this task
            //  int_eng_source = thermodynamic work = f(delP_Dilatation)
            //  modelEng_src   = enthalpy of reaction gained by the
            //                   product gas, PLUS (OR, MINUS) the
            //                   internal energy of the reactant
            //                   material that was liberated in the
            //                   reaction
            // min_int_eng     = a small amount of internal energy to keep
            //                   the equilibration pressure from going nuts

            int_eng_L[c] = int_eng_tmp + int_eng_source[c] + modelEng_src[c];

            int_eng_L[c] = std::max(int_eng_L[c], min_int_eng);
          }
#if 0
          if(massGain > 0.0){
            std::cout << "Mass gained by the models this timestep = " 
                 << massGain << "\t L-" <<level->getIndex()<<endl;
          }
#endif
        } //  if (models.size() > 0)

        //____ B U L L E T   P R O O F I N G----
        // catch negative internal energies
        // ignore BP if timestep restart has already been requested
        IntVector neg_cell;
        bool rts = new_dw->recomputeTimestep();
        if (!areAllValuesPositive(int_eng_L, neg_cell) && !rts) {
          std::ostringstream warn;
          int idx = getLevel(patches)->getIndex();
          warn << "ICE:(L-" << idx << "):computeLagrangianValues, mat " << indx
               << " cell " << neg_cell
               << " Negative int_eng_L: " << int_eng_L[neg_cell] << "\n";
          throw InvalidValue(warn.str(), __FILE__, __LINE__);
        }
      } // if (ice_matl)
    }   // end numALLMatl loop
  }     // patch loop
}

/*!
 *  Function~  ICE::computeLagrangianSpecificVolume--
 */
void
ICE::computeLagrangianSpecificVolume(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset* /*matls*/,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::computeLagrangianSpecificVolume");

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel, level);

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    size_t numALLMatls = d_materialManager->getNumMaterials();
    std::vector<constCCVariable<double>> Tdot(numALLMatls);
    std::vector<constCCVariable<double>> vol_frac(numALLMatls);
    std::vector<constCCVariable<double>> Temp_CC(numALLMatls);
    std::vector<CCVariable<double>> alpha(numALLMatls);
    std::vector<double> if_mpm_matl_ignore(numALLMatls);

    constCCVariable<double> rho_CC, f_theta, sp_vol_CC, cv;
    constCCVariable<double> delP, P;
    constCCVariable<double> TMV_CC;

    CCVariable<double> sum_therm_exp;
    new_dw->allocateTemporary(sum_therm_exp, patch);

    new_dw->get(delP, d_ice_labels->delP_DilatateLabel, 0, patch, d_gn, 0);
    new_dw->get(P, d_ice_labels->press_CCLabel, 0, patch, d_gn, 0);
    sum_therm_exp.initialize(0.);

    if (d_with_mpm) {
      new_dw->get(TMV_CC, d_ice_labels->TMV_CCLabel, 0, patch, d_gn, 0);
    } else {
      CCVariable<double> TMV_create;
      new_dw->allocateTemporary(TMV_create, patch);
      TMV_create.initialize(vol);
      TMV_CC = TMV_create; // reference created data
    }

    for (size_t m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      int indx = matl->getDWIndex();

      new_dw->get(Tdot[m], d_ice_labels->Tdot_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        vol_frac[m], d_ice_labels->vol_frac_CCLabel, indx, patch, d_gac, 1);
      new_dw->allocateTemporary(alpha[m], patch);

      if (ice_matl) {
        old_dw->get(
          Temp_CC[m], d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
      }
      if (mpm_matl) {
        new_dw->get(
          Temp_CC[m], d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
      }
    }

    //__________________________________
    // Sum of thermal expansion
    // ignore contributions from mpm_matls
    // UNTIL we have temperature dependent EOS's for the solids
    for (size_t m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      int indx = matl->getDWIndex();

      if (ice_matl) {
        if_mpm_matl_ignore[m] = 1.0;
        new_dw->get(sp_vol_CC,
                    d_ice_labels->specificVolume_CCLabel,
                    indx,
                    patch,
                    d_gn,
                    0);
        new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          alpha[m][c] = ice_matl->getEOS()->getAlpha(
            Temp_CC[m][c], sp_vol_CC[c], P[c], cv[c]);
          sum_therm_exp[c] += vol_frac[m][c] * alpha[m][c] * Tdot[m][c];
        }
      } else {
        if_mpm_matl_ignore[m] = 0.0;
        alpha[m].initialize(0.0);
      }
    }

    //__________________________________
    for (size_t m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);

      int indx = matl->getDWIndex();

      CCVariable<double> sp_vol_L, sp_vol_src;
      constCCVariable<double> kappa;
      new_dw->allocateAndPut(
        sp_vol_L, d_ice_labels->sp_vol_L_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_src, d_ice_labels->sp_vol_src_CCLabel, indx, patch);
      sp_vol_src.initialize(0.);
      double tiny_rho = 1.e-12;

      if (ice_matl) {
        tiny_rho = ice_matl->getTinyRho();
      }

      new_dw->get(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(f_theta, d_ice_labels->f_theta_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        kappa, d_ice_labels->compressibilityLabel, indx, patch, d_gn, 0);

      //__________________________________
      //  compute sp_vol_L * mass
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        sp_vol_L[c] = (rho_CC[c] * vol) * sp_vol_CC[c];
      }

      //__________________________________
      //   Contributions from models
      constCCVariable<double> Modelsp_vol_src;
      if (d_models.size() > 0) {
        new_dw->get(Modelsp_vol_src,
                    d_ice_labels->modelVol_srcLabel,
                    indx,
                    patch,
                    d_gn,
                    0);
        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          sp_vol_L[c] += Modelsp_vol_src[c];
        }
      }

      //__________________________________
      //  add the sources to sp_vol_L
      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        //__________________________________
        //  term1
        //        double term1 = -vol_frac[m][c] * kappa[c] * TMV_CC[c] *
        //        delP[c];
        double term1 = -vol_frac[m][c] * kappa[c] * vol * delP[c];
        //        double term2 = delT * TMV_CC[c] *
        double term2 = delT * vol *
                       (vol_frac[m][c] * alpha[m][c] * Tdot[m][c] -
                        f_theta[c] * sum_therm_exp[c]);

        // This is actually mass * sp_vol
        double src = term1 + if_mpm_matl_ignore[m] * term2;
        sp_vol_L[c] += src;
        sp_vol_src[c] = src / (rho_CC[c] * vol);
      }

      if (d_clampSpecificVolume) {
        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          /*`==========TESTING==========*/
          sp_vol_L[c] = std::max(sp_vol_L[c], tiny_rho * vol * sp_vol_CC[c]);
          /*==========TESTING==========`*/
        }
      }

      //__________________________________
      // Apply boundary conditions
      setSpecificVolBC(sp_vol_L,
                       "SpecificVol",
                       true,
                       rho_CC,
                       vol_frac[m],
                       patch,
                       d_materialManager,
                       indx);

      //____ B U L L E T   P R O O F I N G----
      // ignore BP if timestep restart has already been requested
      IntVector neg_cell;
      bool rts = new_dw->recomputeTimestep();

      if (!areAllValuesPositive(sp_vol_L, neg_cell) && !rts) {
        std::cout << "\nICE:WARNING......Negative specific Volume" << std::endl;
        std::cout << "cell              " << neg_cell << " level "
                  << level->getIndex() << std::endl;
        std::cout << "matl              " << indx << std::endl;
        std::cout << "sum_thermal_exp   " << sum_therm_exp[neg_cell]
                  << std::endl;
        std::cout << "sp_vol_src        " << sp_vol_src[neg_cell] << std::endl;
        std::cout << "mass sp_vol_L     " << sp_vol_L[neg_cell] << std::endl;
        std::cout << "mass sp_vol_L_old "
                  << (rho_CC[neg_cell] * vol * sp_vol_CC[neg_cell])
                  << std::endl;
        std::cout << "-----------------------------------" << endl;
        //         std::ostringstream warn;
        //        int L = level->getIndex();
        //        warn<<"ERROR ICE:("<<L<<"):computeLagrangianSpecificVolumeRF,
        //        mat "<<indx
        //            << " cell " <<neg_cell << " sp_vol_L is negative\n";
        //        throw InvalidValue(warn.str(), __FILE__, __LINE__);
        new_dw->put(bool_or_vartype(true), VarLabel::find(abortTimestep_name));
        new_dw->put(bool_or_vartype(true),
                    VarLabel::find(recomputeTimestep_name));
      }
    } // end numALLMatl loop
  }   // patch loop
}

/*!
 *  Function~  ICE::computeLagrangian_Transported_Vars--
 *  Purpose:  Multiply the passive scalars by mass_L_CC
 */
void
ICE::computeLagrangian_Transported_Vars(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* /*matls*/,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimestep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, m_ice_tasks, "ICE::computeLagrangian_Transported_Vars");

    // get mass_L for all ice matls
    size_t numMatls = d_materialManager->getNumMaterials("ICE");
    std::vector<constCCVariable<double>> mass_L(numMatls);
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = matl->getDWIndex();

      new_dw->get(
        mass_L[m], d_ice_labels->mass_L_CCLabel, indx, patch, d_gn, 0);
    }

    //__________________________________
    //  hit all the transported variables
    for (auto& model : d_models) {
      FluidsBasedModel* fb_model = dynamic_cast<FluidsBasedModel*>(model.get());
      if (fb_model) {
        for (auto& tvar : fb_model->d_transVars) {
          for (size_t m = 0; m < numMatls; m++) {
            ICEMaterial* matl = static_cast<ICEMaterial*>(
              d_materialManager->getMaterial("ICE", m));
            int indx = matl->getDWIndex();

            if (tvar->matls->contains(indx)) {
              constCCVariable<double> q_old, q_src;
              CCVariable<double> q_L_CC;
              old_dw->get(q_old, tvar->var, indx, patch, d_gn, 0);
              new_dw->allocateAndPut(q_L_CC, tvar->var_Lagrangian, indx, patch);

              // initialize q_L to q_old
              q_L_CC.copyData(q_old);

              // If there's a source tack it on.
              if (tvar->src) {
                new_dw->get(q_src, tvar->src, indx, patch, d_gn, 0);
                for (auto iter = patch->getCellIterator(); !iter.done();
                     iter++) {
                  IntVector c = *iter;
                  q_L_CC[c] += q_src[c]; // with source
                }
              }

              // Set boundary conditions on q_L_CC
              // must use var Labelname not var_Lagrangian
              std::string Labelname = tvar->var->getName();
              setBC(q_L_CC,
                    Labelname,
                    patch,
                    d_materialManager,
                    indx,
                    new_dw,
                    isNotInitialTimestep);

              // multiply by mass so advection is conserved
              for (auto iter = patch->getExtraCellIterator(); !iter.done();
                   iter++) {
                IntVector c = *iter;
                q_L_CC[c] *= mass_L[m][c];
              }

            } // tvar matl
          }   // ice matl loop
        }     // tvar loop
      }       // if fb_model
    }         // model loop
  }           // patch loop
}

/*!
 *  Function~  ICE::maxMach_on_Lodi_BC_Faces
 *  Purpose~   Find the max mach Number on all lodi faces
 */
void
ICE::maxMach_on_Lodi_BC_Faces(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* /*matls*/,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::::maxMach_on_Lodi_BC_Faces");

    size_t numAllMatls = d_materialManager->getNumMaterials();
    std::vector<constCCVariable<Vector>> vel_CC(numAllMatls);
    std::vector<constCCVariable<double>> speedSound(numAllMatls);

    for (size_t m = 0; m < numAllMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      if (ice_matl) {
        old_dw->get(
          vel_CC[m], d_ice_labels->velocity_CCLabel, indx, patch, d_gn, 0);
        old_dw->get(speedSound[m],
                    d_ice_labels->speedSound_CCLabel,
                    indx,
                    patch,
                    d_gn,
                    0);
      }
    }

    //__________________________________
    // Work on the lodi faces for each patch.
    // Every patch has to compute a maxMach
    // even if it isn't on a boundary.  We
    // can't do reduction variables with patch subsets yet.
    std::vector<Patch::FaceType>::iterator f;

    for (auto& face : d_BC_globalVars->lodi->LodiFaces) {

      //__________________________________
      // compute maxMach number on this lodi face
      // only ICE matls
      double maxMach = 0.0;
      if (is_LODI_face(patch, face, d_materialManager)) {

        for (size_t m = 0; m < numAllMatls; m++) {
          Material* matl = d_materialManager->getMaterial(m);
          int indx       = matl->getDWIndex();

          ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
          if (ice_matl) {
            Patch::FaceIteratorType MEC = Patch::ExtraMinusEdgeCells;

            for (auto iter = patch->getFaceIterator(face, MEC); !iter.done();
                 iter++) {
              IntVector c = *iter;
              maxMach = Max(maxMach, vel_CC[m][c].length() / speedSound[m][c]);
            }

            VarLabel* V_Label = getMaxMach_face_VarLabel(face);
            new_dw->put(max_vartype(maxMach), V_Label, level, indx);
          } // icematl
        }   // matl loop

      } // is lodi Face
    }   // boundaryFaces
  }     // patches
}

/*!
 *  Function~  ICE::advectAndAdvanceInTime--
 *  Purpose~  This task calculates the The cell-centered, time n+1, mass,
 *            momentum internal energy, sp_vol
 */
void
ICE::advectAndAdvanceInTime(const ProcessorGroup* /*pg*/,
                            const PatchSubset* patches,
                            const MaterialSubset* /*matls*/,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  // the advection calculations care about the position of the old dw subcycle
  double AMR_subCycleProgressVar = getSubCycleProgress(old_dw);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    std::ostringstream msg;
    msg << "ICE::advectAndAdvanceInTime: progressVar "
        << AMR_subCycleProgressVar;
    printTask(patches, patch, m_ice_tasks, msg.str());

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel, level);

    std::unique_ptr<Advector> advector =
      d_advector->clone(new_dw, patch, isRegridTimestep());

    CCVariable<double> q_advected;
    CCVariable<Vector> qV_advected;
    new_dw->allocateTemporary(q_advected, patch);
    new_dw->allocateTemporary(qV_advected, patch);

    size_t numMatls = d_materialManager->getNumMaterials("ICE");
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = matl->getDWIndex();

      CCVariable<double> mass_adv, int_eng_adv, sp_vol_adv;
      CCVariable<Vector> mom_adv;

      constCCVariable<double> int_eng_L_ME, mass_L, sp_vol_L;
      constCCVariable<Vector> mom_L_ME;
      constSFCXVariable<double> uvel_FC;
      constSFCYVariable<double> vvel_FC;
      constSFCZVariable<double> wvel_FC;

      new_dw->get(uvel_FC, d_ice_labels->uvel_FCMELabel, indx, patch, d_gac, 2);
      new_dw->get(vvel_FC, d_ice_labels->vvel_FCMELabel, indx, patch, d_gac, 2);
      new_dw->get(wvel_FC, d_ice_labels->wvel_FCMELabel, indx, patch, d_gac, 2);

      new_dw->get(mass_L, d_ice_labels->mass_L_CCLabel, indx, patch, d_gac, 2);
      new_dw->get(
        mom_L_ME, d_ice_labels->mom_L_ME_CCLabel, indx, patch, d_gac, 2);
      new_dw->get(
        sp_vol_L, d_ice_labels->sp_vol_L_CCLabel, indx, patch, d_gac, 2);
      new_dw->get(
        int_eng_L_ME, d_ice_labels->eng_L_ME_CCLabel, indx, patch, d_gac, 2);

      new_dw->allocateAndPut(
        mass_adv, d_ice_labels->mass_advLabel, indx, patch);
      new_dw->allocateAndPut(mom_adv, d_ice_labels->mom_advLabel, indx, patch);
      new_dw->allocateAndPut(
        int_eng_adv, d_ice_labels->eng_advLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_adv, d_ice_labels->sp_vol_advLabel, indx, patch);

      mass_adv.initialize(0.0);
      mom_adv.initialize(Vector(0.0, 0.0, 0.0));
      int_eng_adv.initialize(0.0);
      sp_vol_adv.initialize(0.0);
      q_advected.initialize(0.0);
      qV_advected.initialize(Vector(0.0, 0.0, 0.0));

      //__________________________________
      // common variables that get passed into the advection operators
      std::unique_ptr<advectVarBasket> varBasket =
        std::make_unique<advectVarBasket>();
      varBasket->new_dw                  = new_dw;
      varBasket->old_dw                  = old_dw;
      varBasket->indx                    = indx;
      varBasket->patch                   = patch;
      varBasket->level                   = level;
      varBasket->lb                      = d_ice_labels.get();
      varBasket->doRefluxing             = d_doRefluxing;
      varBasket->useCompatibleFluxes     = d_useCompatibleFluxes;
      varBasket->AMR_subCycleProgressVar = AMR_subCycleProgressVar;

      //__________________________________
      //   Advection preprocessing
      bool bulletProof_test = true;
      advector->inFluxOutFluxVolume(uvel_FC,
                                    vvel_FC,
                                    wvel_FC,
                                    delT,
                                    patch,
                                    indx,
                                    bulletProof_test,
                                    new_dw,
                                    varBasket.get());
      //__________________________________
      // mass
      advector->advectMass(mass_L, q_advected, varBasket.get());

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        mass_adv[c] = (mass_L[c] + q_advected[c]);
      }

      //__________________________________
      // momentum
      varBasket->is_Q_massSpecific = true;
      varBasket->desc              = "mom";
      advector->advectQ(mom_L_ME, mass_L, qV_advected, varBasket.get());

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        mom_adv[c]  = (mom_L_ME[c] + qV_advected[c]);
      }

      //__________________________________
      // internal energy
      varBasket->is_Q_massSpecific = true;
      varBasket->desc              = "int_eng";
      advector->advectQ(int_eng_L_ME, mass_L, q_advected, varBasket.get());

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c    = *iter;
        int_eng_adv[c] = (int_eng_L_ME[c] + q_advected[c]);
      }

      //__________________________________
      // sp_vol[m] * mass
      varBasket->is_Q_massSpecific = true;
      varBasket->desc              = "sp_vol";
      advector->advectQ(sp_vol_L, mass_L, q_advected, varBasket.get());

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c   = *iter;
        sp_vol_adv[c] = (sp_vol_L[c] + q_advected[c]);
      }

      //__________________________________
      // Advect model variables
      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());
        if (fb_model) {
          for (auto& tvar : fb_model->d_transVars) {
            if (tvar->matls->contains(indx)) {
              std::string Labelname = tvar->var->getName();

              CCVariable<double> q_adv;
              constCCVariable<double> q_L_CC;

              new_dw->allocateAndPut(q_adv, tvar->var_adv, indx, patch);
              new_dw->get(q_L_CC, tvar->var_Lagrangian, indx, patch, d_gac, 2);

              q_adv.initialize(d_EVIL_NUM);

              varBasket->desc              = Labelname;
              varBasket->is_Q_massSpecific = true;

              advector->advectQ(q_L_CC, mass_L, q_advected, varBasket.get());

              for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                q_adv[c]    = (q_L_CC[c] + q_advected[c]);
              }

            } // if tvar contains
          }   // tvar loop
        }     // if fb_model
      }       // models loop
    }         // ice_matls loop
  }           // patch loop
}

/*!
 *  Function~  ICE::conservedtoPrimitive_Vars
 *  Purpose~ This task computes the primitive variables (rho,T,vel,sp_vol,...)
 *           at time n+1, from the conserved variables mass, momentum, energy...
 */
void
ICE::conservedtoPrimitive_Vars(const ProcessorGroup* /*pg*/,
                               const PatchSubset* patches,
                               const MaterialSubset* /*matls*/,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimestep = (timeStep > 0);

  const Level* level = getLevel(patches);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::conservedtoPrimitive_Vars");

    Vector dx     = patch->dCell();
    double invvol = 1.0 / (dx.x() * dx.y() * dx.z());

    size_t numMatls = d_materialManager->getNumMaterials("ICE");
    for (size_t m = 0; m < numMatls; m++) {
      Material* matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = matl->getDWIndex();

      CCVariable<double> rho_CC, temp_CC, sp_vol_CC, mach;
      CCVariable<Vector> vel_CC;

      constCCVariable<double> int_eng_adv, mass_adv, sp_vol_adv, speedSound, cv;
      constCCVariable<double> gamma, placeHolder, vol_frac;
      constCCVariable<Vector> mom_adv;

      new_dw->get(gamma, d_ice_labels->gammaLabel, indx, patch, d_gn, 0);
      new_dw->get(
        speedSound, d_ice_labels->speedSound_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, d_gn, 0);
      new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);

      new_dw->get(mass_adv, d_ice_labels->mass_advLabel, indx, patch, d_gn, 0);
      new_dw->get(mom_adv, d_ice_labels->mom_advLabel, indx, patch, d_gn, 0);
      new_dw->get(
        sp_vol_adv, d_ice_labels->sp_vol_advLabel, indx, patch, d_gn, 0);
      new_dw->get(
        int_eng_adv, d_ice_labels->eng_advLabel, indx, patch, d_gn, 0);

      new_dw->getModifiable(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->getModifiable(rho_CC, d_ice_labels->rho_CCLabel, indx, patch);

      new_dw->allocateAndPut(
        temp_CC, d_ice_labels->temperature_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        vel_CC, d_ice_labels->velocity_CCLabel, indx, patch);
      new_dw->allocateAndPut(mach, d_ice_labels->machLabel, indx, patch);

      rho_CC.initialize(-d_EVIL_NUM);
      temp_CC.initialize(-d_EVIL_NUM);
      vel_CC.initialize(Vector(0.0, 0.0, 0.0));

      //__________________________________
      // Backout primitive quantities from
      // the conserved ones.
      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c         = *iter;
        double inv_mass_adv = 1.0 / mass_adv[c];
        rho_CC[c]           = mass_adv[c] * invvol;
        vel_CC[c]           = mom_adv[c] * inv_mass_adv;
        sp_vol_CC[c]        = sp_vol_adv[c] * inv_mass_adv;
      }

      //__________________________________
      // model variables
      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());
        if (fb_model) {
          for (auto& tvar : fb_model->d_transVars) {
            if (tvar->matls->contains(indx)) {
              std::string Labelname = tvar->var->getName();

              CCVariable<double> q_CC;
              constCCVariable<double> q_adv;

              new_dw->allocateAndPut(q_CC, tvar->var, indx, patch);
              new_dw->get(q_adv, tvar->var_adv, indx, patch, d_gn, 0);
              q_CC.initialize(0.0);

              for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
                IntVector c = *iter;
                q_CC[c]     = q_adv[c] / mass_adv[c];
              }

              //  Set Boundary Conditions
              setBC(q_CC,
                    Labelname,
                    patch,
                    d_materialManager,
                    indx,
                    new_dw,
                    isNotInitialTimestep);
            }
          }
        }
      }

      //__________________________________
      // A model *can* compute the specific heat
      CCVariable<double> cv_new;
      new_dw->allocateTemporary(cv_new, patch, d_gn, 0);
      cv_new.copyData(cv);

      for (auto& model : d_models) {
        FluidsBasedModel* fb_model =
          dynamic_cast<FluidsBasedModel*>(model.get());
        if (fb_model) {
          fb_model->computeSpecificHeat(cv_new, patch, new_dw, indx);
        }
      }

      //__________________________________
      // Backout primitive quantities from
      // the conserved ones.
      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        temp_CC[c]  = int_eng_adv[c] / (mass_adv[c] * cv_new[c]);
      }

      //__________________________________
      // set the boundary conditions
      std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
        std::make_unique<CustomBCDriver::customBC_localVars>();
      preprocess_CustomBCs("Advection",
                           old_dw,
                           new_dw,
                           d_ice_labels.get(),
                           patch,
                           indx,
                           d_BC_globalVars.get(),
                           BC_localVars.get());

      setBC(rho_CC,
            "Density",
            placeHolder,
            placeHolder,
            patch,
            d_materialManager,
            indx,
            new_dw,
            d_BC_globalVars.get(),
            BC_localVars.get(),
            isNotInitialTimestep);
      setBC(vel_CC,
            "Velocity",
            patch,
            d_materialManager,
            indx,
            new_dw,
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
            new_dw,
            d_BC_globalVars.get(),
            BC_localVars.get(),
            isNotInitialTimestep);

      setSpecificVolBC(sp_vol_CC,
                       "SpecificVol",
                       false,
                       rho_CC,
                       vol_frac,
                       patch,
                       d_materialManager,
                       indx);

      delete_CustomBCs(d_BC_globalVars.get(), BC_localVars.get());

      //__________________________________
      // Compute Auxilary quantities
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        mach[c]     = vel_CC[c].length() / speedSound[c];
      }

      //____ B U L L E T   P R O O F I N G----
      // ignore BP if timestep restart has already been requested
      IntVector neg_cell;
      bool rts = new_dw->recomputeTimestep();

      std::ostringstream base, warn;
      base << "ERROR ICE:(L-" << level->getIndex()
           << "):conservedtoPrimitive_Vars, mat " << indx << " cell ";
      if (!areAllValuesPositive(rho_CC, neg_cell) && !rts) {
        warn << base.str() << neg_cell << " negative rho_CC\n ";
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(temp_CC, neg_cell) && !rts) {
        warn << base.str() << neg_cell << " negative temp_CC\n ";
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(sp_vol_CC, neg_cell) && !rts) {
        warn << base.str() << neg_cell << " negative sp_vol_CC\n ";
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
    } // ice_matls loop
  }   // patch loop
}

/*!
 * Function:  TestConservation--
 * Purpose:   Test for conservation of mass, momentum, energy.
 * Test to see if the exchange process is conserving
 */
void
ICE::TestConservation(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* /*matls*/,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  delt_vartype delT;
  old_dw->get(delT, d_ice_labels->delTLabel, level);

  size_t numICEmatls = d_materialManager->getNumMaterials("ICE");
  std::map<int, double> total_mass;
  std::map<int, double> total_KE;
  std::map<int, double> total_int_eng;
  std::map<int, Vector> total_mom;

  double allMatls_totalMass   = 0.0;
  double allMatls_totalKE     = 0.0;
  double allMatls_totalIntEng = 0.0;
  Vector allMatls_totalMom(0.0, 0.0, 0.0);

  Vector mom_exch_error(0, 0, 0);
  double eng_exch_error = 0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, m_ice_tasks, "ICE::TestConservation");

    Vector dx       = patch->dCell();
    double cell_vol = dx.x() * dx.y() * dx.z();

    //__________________________________
    // get face centered velocities to
    // to compute what's being fluxed through the domain
    std::vector<constSFCXVariable<double>> uvel_FC(numICEmatls);
    std::vector<constSFCYVariable<double>> vvel_FC(numICEmatls);
    std::vector<constSFCZVariable<double>> wvel_FC(numICEmatls);

    for (size_t m = 0; m < numICEmatls; m++) {
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      new_dw->get(
        uvel_FC[m], d_ice_labels->uvel_FCMELabel, indx, patch, d_gn, 0);
      new_dw->get(
        vvel_FC[m], d_ice_labels->vvel_FCMELabel, indx, patch, d_gn, 0);
      new_dw->get(
        wvel_FC[m], d_ice_labels->wvel_FCMELabel, indx, patch, d_gn, 0);
    }

    //__________________________________
    // conservation of mass
    constCCVariable<double> rho_CC;
    std::vector<CCVariable<double>> mass(numICEmatls);
    for (size_t m = 0; m < numICEmatls; m++) {

      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      new_dw->allocateTemporary(mass[m], patch);
      new_dw->get(rho_CC, d_ice_labels->rho_CCLabel, indx, patch, d_gn, 0);

      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        mass[m][c]  = rho_CC[c] * cell_vol;
      }
    }

    if (d_conservationTest->mass) {
      for (size_t m = 0; m < numICEmatls; m++) {
        ICEMaterial* ice_matl =
          static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
        int indx         = ice_matl->getDWIndex();
        total_mass[indx] = 0.0;
        double mat_mass  = 0;

        conservationTest<double>(
          patch, delT, mass[m], uvel_FC[m], vvel_FC[m], wvel_FC[m], mat_mass);

        total_mass[indx] += mat_mass;
        allMatls_totalMass += mat_mass;
      }
    }

    //__________________________________
    // conservation of momentum
    if (d_conservationTest->momentum) {
      CCVariable<Vector> mom;
      constCCVariable<Vector> vel_CC;
      new_dw->allocateTemporary(mom, patch);

      for (size_t m = 0; m < numICEmatls; m++) {

        ICEMaterial* ice_matl =
          static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
        int indx = ice_matl->getDWIndex();

        new_dw->get(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          mom[c]      = mass[m][c] * vel_CC[c];
        }

        Vector mat_mom(0, 0, 0);
        conservationTest<Vector>(
          patch, delT, mom, uvel_FC[m], vvel_FC[m], wvel_FC[m], mat_mom);

        total_mom[indx] += mat_mom;
        allMatls_totalMom += mat_mom;
      }
    }

    //__________________________________
    // conservation of internal_energy
    if (d_conservationTest->energy) {
      CCVariable<double> int_eng;
      constCCVariable<double> temp_CC;
      constCCVariable<double> cv;

      new_dw->allocateTemporary(int_eng, patch);

      for (size_t m = 0; m < numICEmatls; m++) {

        ICEMaterial* ice_matl =
          static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
        int indx = ice_matl->getDWIndex();

        new_dw->get(
          temp_CC, d_ice_labels->temperature_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(cv, d_ice_labels->specific_heatLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          int_eng[c]  = mass[m][c] * cv[c] * temp_CC[c];
        }

        double mat_int_eng(0);
        conservationTest<double>(patch,
                                 delT,
                                 int_eng,
                                 uvel_FC[m],
                                 vvel_FC[m],
                                 wvel_FC[m],
                                 mat_int_eng);

        total_int_eng[indx] += mat_int_eng;
        allMatls_totalIntEng += mat_int_eng;
      }
    }

    //__________________________________
    // conservation of kinetic_energy
    if (d_conservationTest->energy) {
      CCVariable<double> KE;
      constCCVariable<Vector> vel_CC;
      new_dw->allocateTemporary(KE, patch);

      for (size_t m = 0; m < numICEmatls; m++) {

        ICEMaterial* ice_matl =
          static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
        int indx = ice_matl->getDWIndex();

        new_dw->get(
          vel_CC, d_ice_labels->velocity_CCLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c    = *iter;
          double vel_mag = vel_CC[c].length();
          KE[c]          = 0.5 * mass[m][c] * vel_mag * vel_mag;
        }

        double mat_KE(0);
        conservationTest<double>(
          patch, delT, KE, uvel_FC[m], vvel_FC[m], wvel_FC[m], mat_KE);

        total_KE[indx] += mat_KE;
        allMatls_totalKE += mat_KE;
      }
    }

    //__________________________________
    // conservation during the exchange process
    if (d_conservationTest->exchange) {
      Vector sum_mom_L_CC     = Vector(0.0, 0.0, 0.0);
      Vector sum_mom_L_ME_CC  = Vector(0.0, 0.0, 0.0);
      double sum_int_eng_L_CC = 0.0;
      double sum_eng_L_ME_CC  = 0.0;

      size_t numALLmatls = d_materialManager->getNumMaterials();
      for (size_t m = 0; m < numALLmatls; m++) {
        Material* matl = d_materialManager->getMaterial(m);
        int indx       = matl->getDWIndex();

        constCCVariable<double> int_eng_L_CC, eng_L_ME_CC;
        constCCVariable<Vector> mom_L_CC, mom_L_ME_CC;

        new_dw->get(
          mom_L_CC, d_ice_labels->mom_L_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(
          int_eng_L_CC, d_ice_labels->int_eng_L_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(
          mom_L_ME_CC, d_ice_labels->mom_L_ME_CCLabel, indx, patch, d_gn, 0);
        new_dw->get(
          eng_L_ME_CC, d_ice_labels->eng_L_ME_CCLabel, indx, patch, d_gn, 0);

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          sum_mom_L_CC += mom_L_CC[c];
          sum_mom_L_ME_CC += mom_L_ME_CC[c];
          sum_int_eng_L_CC += int_eng_L_CC[c];
          sum_eng_L_ME_CC += eng_L_ME_CC[c];
        }
      }

      mom_exch_error = sum_mom_L_CC - sum_mom_L_ME_CC;
      eng_exch_error = sum_int_eng_L_CC - sum_eng_L_ME_CC;
    }
  } // patch loop

  //  Put variables into the DW.
  const MaterialSubset* matls =
    d_materialManager->allMaterials("ICE")->getUnion();
  if (d_conservationTest->mass) {
    new_dw->put(sum_vartype(allMatls_totalMass),
                d_ice_labels->TotalMassLabel,
                nullptr,
                -1);
    new_dw->put_sum_vartype(total_mass, d_ice_labels->TotalMassLabel, matls);
  }

  if (d_conservationTest->exchange) {
    new_dw->put(sumvec_vartype(mom_exch_error),
                d_ice_labels->mom_exch_errorLabel);
    new_dw->put(sum_vartype(eng_exch_error), d_ice_labels->eng_exch_errorLabel);
  }

  if (d_conservationTest->energy) {
    new_dw->put(sum_vartype(allMatls_totalKE),
                d_ice_labels->KineticEnergyLabel,
                nullptr,
                -1);
    new_dw->put(sum_vartype(allMatls_totalIntEng),
                d_ice_labels->TotalIntEngLabel,
                nullptr,
                -1);

    new_dw->put_sum_vartype(total_KE, d_ice_labels->KineticEnergyLabel, matls);
    new_dw->put_sum_vartype(
      total_int_eng, d_ice_labels->TotalIntEngLabel, matls);
  }

  if (d_conservationTest->momentum) {
    new_dw->put(sumvec_vartype(allMatls_totalMom),
                d_ice_labels->TotalMomentumLabel,
                nullptr,
                -1);
    new_dw->put_sum_vartype(total_mom, d_ice_labels->TotalMomentumLabel, matls);
  }
}

/*!
 *  Function:  hydrostaticPressureAdjustment--
 *  Notes:     press_hydro = rho_micro_CC[SURROUNDING_MAT] * grav *
 * some_distance
 */
void
ICE::hydrostaticPressureAdjustment(const Patch* patch,
                                   const CCVariable<double>& rho_micro_CC,
                                   CCVariable<double>& press_CC)
{
  Vector gravity = getGravity();

  // find the upper and lower point of the domain.
  const Level* level = patch->getLevel();
  GridP grid         = level->getGrid();
  BBox b;
  grid->getSpatialRange(b);
  Vector gridMin = b.min().asVector();
  Vector dx_L0   = grid->getLevel(0)->dCell();

  // Pressure reference point is assumed to be
  // at CELL-CENTER of cell 0,0,0
  Vector press_ref_pt = gridMin + 1.5 * dx_L0;

  // Which direction is the gravitational vector pointing
  int dir = -9;
  for (int i = 0; i < 3; i++) {
    if (gravity[i] != 0.0) {
      dir = i;
    }
  }

  //__________________________________
  // Tack on the hydrostatic pressure adjustment
  for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
    IntVector c            = *iter;
    Point here             = level->getCellPosition(c);
    Vector dist_from_p_ref = (here.asVector() - press_ref_pt);

    double press_hydro = rho_micro_CC[c] * gravity[dir] * dist_from_p_ref[dir];
    press_CC[c] += press_hydro;
  }
}

/*!
 *  Function~  ICE::upwindCell--
 *  purpose:   Find the upwind cell in each direction  This is a knock off
 *             of Bucky's logic
 */
IntVector
ICE::upwindCell_X(const IntVector& c,
                  const double& var,
                  double is_logical_R_face)
{
  double plus_minus_half = 0.5 * (var + d_SMALL_NUM) / fabs(var + d_SMALL_NUM);
  int one_or_zero        = int(-0.5 - plus_minus_half + is_logical_R_face);
  IntVector tmp          = c + IntVector(one_or_zero, 0, 0);
  return tmp;
}

IntVector
ICE::upwindCell_Y(const IntVector& c,
                  const double& var,
                  double is_logical_R_face)
{
  double plus_minus_half = 0.5 * (var + d_SMALL_NUM) / fabs(var + d_SMALL_NUM);
  int one_or_zero        = int(-0.5 - plus_minus_half + is_logical_R_face);
  IntVector tmp          = c + IntVector(0, one_or_zero, 0);
  return tmp;
}

IntVector
ICE::upwindCell_Z(const IntVector& c,
                  const double& var,
                  double is_logical_R_face)
{
  double plus_minus_half = 0.5 * (var + d_SMALL_NUM) / fabs(var + d_SMALL_NUM);
  int one_or_zero        = int(-0.5 - plus_minus_half + is_logical_R_face);
  IntVector tmp          = c + IntVector(0, 0, one_or_zero);
  return tmp;
}

//_____________________________________________________________________
//   Stub functions for AMR
void
ICE::addRefineDependencies(Task*, const VarLabel*, int, int)
{
}

void
ICE::refineBoundaries(const Patch*,
                      CCVariable<double>&,
                      DataWarehouse*,
                      const VarLabel*,
                      int,
                      double)
{
  throw InternalError(
    "trying to do AMR iwth the non-AMR component!", __FILE__, __LINE__);
}

void
ICE::refineBoundaries(const Patch*,
                      CCVariable<Vector>&,
                      DataWarehouse*,
                      const VarLabel*,
                      int,
                      double)
{
  throw InternalError(
    "trying to do AMR iwth the non-AMR component!", __FILE__, __LINE__);
}

void
ICE::refineBoundaries(const Patch*,
                      SFCXVariable<double>&,
                      DataWarehouse*,
                      const VarLabel*,
                      int,
                      double)
{
  throw InternalError(
    "trying to do AMR iwth the non-AMR component!", __FILE__, __LINE__);
}

void
ICE::refineBoundaries(const Patch*,
                      SFCYVariable<double>&,
                      DataWarehouse*,
                      const VarLabel*,
                      int,
                      double)
{
  throw InternalError(
    "trying to do AMR iwth the non-AMR component!", __FILE__, __LINE__);
}

void
ICE::refineBoundaries(const Patch*,
                      SFCZVariable<double>&,
                      DataWarehouse*,
                      const VarLabel*,
                      int,
                      double)
{
  throw InternalError(
    "trying to do AMR iwth the non-AMR component!", __FILE__, __LINE__);
}

/*______________________________________________________________________
  S C H E M A T I C   D I A G R A M S

  q_outflux(TOP)

  |    (I/O)flux_EF(TOP_BK)
  |
  (I/O)flux_CF(TOP_L_BK)       _________|___________
  /___/_____|_______/__/|   (I/O)flux_CF(TOP_R_BK)
  /   /      |      /  | |
  /   /       |     /  /| |
  (I/O)flux_EF(TOP_L)      /   /             /  / |/|
  /___/_____________/__/ ------ (I/O)flux_EF(TOP_R)
  _/__ /_____________/__/| /| |
  |   |             |  | |/ | |   (I/O)flux_EF(BCK_R)
  | + |      +      | +| /  | |
  |---|----------------|/|  |/|
  |   |             |  | | /| /  (I/O)flux_CF(BOT_R_BK)
  (I/O)flux(LEFT_FR)    | + |     i,j,k   | +| |/ /
  |   |             |  |/| /   (I/O)flux_EF(BOT_R)
  |---|----------------| |/
  (I/O)flux_CF(BOT_L_FR)| + |      +      | +|/    (I/O)flux_CF(BOT_R_FR)
  ----------------------
  (I/O)flux_EF(BOT_FR)



  (TOP)
  ______________________              ______________________  _
  |   |             |  |              |   |             |  |  |  delY_top
  | + |      +      | +|              | + |      +      | +|  |
  |---|----------------|  --ytop      |---|----------------|  -
  |   |             |  |              |   |             |  |
  | + |     i,j,k   | +| (RIGHT)      | + |     i,j,k   | +|
  |   |             |  |              |   |             |  |
  |---|----------------|  --y0        |---|----------------|  -
  | + |      +      | +|              | + |      +      | +|  | delY_bottom
  ----------------------              ----------------------  -
  |             |                 |---|             |--|
  x0            xright              delX_left         delX_right

  (BACK)
  ______________________              ______________________  _
  |   |             |  |              |   |             |  |  |  delZ_back
  | + |      +      | +|              | + |      +      | +|  |
  |---|----------------|  --z0        |---|----------------|  -
  |   |             |  |              |   |             |  |
  | + |     i,j,k   | +| (RIGHT)      | + |     i,j,k   | +|
  |   |             |  |              |   |             |  |
  |---|----------------|  --z_frt     |---|----------------|  -
  | + |      +      | +|              | + |      +      | +|  | delZ_front
  ----------------------              ----------------------  -
  |             |                 |---|             |--|
  x0            xright              delX_left         delX_right

  ______________________________________________________________________*/
