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

// MPMICE.cc
#include <CCA/Components/MPMICE/MPMICE.h>

#include <CCA/Components/MPMICE/Core/MPMICELabel.h>

#include <CCA/Components/ICE/AMRICE.h>
#include <CCA/Components/ICE/Core/BoundaryCond.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <CCA/Components/ICE/EOS/EquationOfState.h>
#include <CCA/Components/ICE/ICE.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/HeatConduction/HeatConductionTasks.h>
#include <CCA/Components/MPM/RigidMPM.h>
#include <CCA/Components/MPM/SerialMPM.h>
#include <CCA/Components/MPM/ShellMPM.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>

#include <CCA/Components/Models/ParticleBased/ParticleModel.h>

#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ConvergenceFailure.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/AMR_CoarsenRefine.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/Utils.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/MiscMath.h>
#include <Core/Util/DebugStream.h>

#include <cfloat>
#include <cstdio>
#include <vector>

#include <errno.h>
#include <fenv.h>
#include <iomanip>

using namespace Uintah;

//__________________________________
//  To turn on normal output
//  setenv SCI_DEBUG "MPMICE_NORMAL_COUT:+,MPMICE_DOING_COUT".....
//  MPMICE_NORMAL_COUT:  dumps out during problemSetup
//  MPMICE_DOING_COUT:   dumps when tasks are scheduled and performed
//  default is OFF

static DebugStream cout_norm("MPMICE_NORMAL_COUT", false);
static DebugStream cout_doing("MPMICE_DOING_COUT", false);
static DebugStream ds_EqPress("DBG_EqPress", false);

MPMICE::MPMICE(const ProcessorGroup* myworld,
               const MaterialManagerP& mat_manager,
               const MPMType& mpmtype,
               bool doAMR)
  : SimulationCommon(myworld, mat_manager)
{
  d_mpmice_labels = std::make_unique<MPMICELabel>();

  switch (mpmtype) {
    case MPMType::RIGID_MPMICE:
      d_mpm      = std::make_unique<RigidMPM>(myworld, mat_manager);
      d_rigidMPM = true;
      break;
    case MPMType::SHELL_MPMICE:
      d_mpm = std::make_unique<ShellMPM>(myworld, mat_manager);
      break;
    default:
      d_mpm = std::make_unique<SerialMPM>(myworld, mat_manager);
      break;
  }

  // Don't do AMRICE with MPMICE for now...
  if (doAMR) {
    d_ice = std::make_unique<AMRICE>(myworld, mat_manager);
  } else {
    d_ice = std::make_unique<ICE>(myworld, mat_manager);
  }

  d_ice_labels = d_ice->d_ice_labels.get();
  d_mpm_labels = d_mpm->d_mpm_labels.get();

  d_SMALL_NUM = d_ice->d_SMALL_NUM;

  // Note, within MPMICE, d_TINY_RHO is only applied  to MPM materials, for
  // which its value is hardcoded, unlike the situation for ice materials
  d_TINY_RHO = 1.e-12;
}

MPMICE::~MPMICE()
{
  d_mpm->releaseComponents();
  d_ice->releaseComponents();
  for (auto& am : d_analysisModules) {
    am->releaseComponents();
  }
}

// For recomputing timesteps
double
MPMICE::recomputeDelT(double delT)
{
  return delT / 2.0;
}

void
MPMICE::problemSetup(const ProblemSpecP& prob_spec,
                     const ProblemSpecP& restart_prob_spec,
                     GridP& grid)
{
  cout_doing << "Doing MPMICE::problemSetup " << std::endl;

  //  M P M
  d_mpm->setComponents(this);
  dynamic_cast<SimulationCommon*>(d_mpm.get())->problemSetup(prob_spec);

  d_mpm->setWithICE();
  d_mpm->problemSetup(prob_spec, restart_prob_spec, grid);

  d_8or27 = d_mpm->d_mpm_flags->d_8or27;
  if (d_8or27 == 8) {
    d_num_ghost_nodes = 1;
  } else if (d_8or27 == 27) {
    d_num_ghost_nodes = 2;
  }

  Ghost::GhostType gp;
  int ngc_p;
  d_mpm->getParticleGhostLayer(gp, ngc_p);

  ProblemSpecP mpm_ps = prob_spec->findBlock("MPM");
  if (!mpm_ps) {
    mpm_ps = restart_prob_spec->findBlock("MPM");
  }
  mpm_ps->get("testForNegTemps_mpm", d_testForNegTemps_mpm);

  //  I C E
  d_ice->setComponents(this);
  dynamic_cast<SimulationCommon*>(d_ice.get())->problemSetup(prob_spec);

  d_ice->setWithMPM();

  // Communicate the particle ghost from MPM to ICE. Used only in the
  // HEChem/Unsteady_Burn model.
  d_ice->setParticleGhostLayer(gp, ngc_p);

  if (d_rigidMPM) {
    d_ice->setWithRigidMPM();
  }

  d_ice->problemSetup(prob_spec, restart_prob_spec, grid);

  //  Component switching
  d_switchCriteria =
    dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));

  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(
      prob_spec, restart_prob_spec, d_materialManager);
  }

  // Validate input
  if (isAMR() && !isLockstepAMR()) {
    std::ostringstream msg;
    msg << "\n ERROR: You must add \n"
        << " <useLockStep> true </useLockStep> \n"
        << " inside of the <AMR> section for MPMICE and AMR. \n";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (cout_norm.active()) {
    cout_norm << "Done with problemSetup \t\t\t MPMICE" << endl;
    cout_norm << "--------------------------------\n" << endl;
  }

  // Flags specific to MPMICE
  d_useSimpleEquilibrationPressure = false;
  double converg_coeff             = 100.;
  double convergence_crit          = converg_coeff * DBL_EPSILON;
  d_convergence_tolerance          = convergence_crit;
  ProblemSpecP mpmice_ps           = prob_spec->findBlock("MPMICE");
  if (mpmice_ps) {
    mpmice_ps->get("use_simple_equilibration_algorithm",
                   d_useSimpleEquilibrationPressure);
    mpmice_ps->getWithDefault("vol_frac_convergence_tolerance",
                              d_convergence_tolerance,
                              convergence_crit);
  } else {
    mpmice_ps = restart_prob_spec->findBlock("MPMICE");
    if (mpmice_ps) {
      mpmice_ps->get("use_simple_equilibration_algorithm",
                     d_useSimpleEquilibrationPressure);
      mpmice_ps->get("vol_frac_convergence_tolerance", d_convergence_tolerance);
    }
  }

  //__________________________________
  //  Set up data analysis modules
  d_analysisModules =
    AnalysisModuleFactory::create(d_myworld, d_materialManager, prob_spec);
  for (auto& am : d_analysisModules) {
    am->setComponents(dynamic_cast<SimulationInterface*>(this));
    am->problemSetup(prob_spec,
                     restart_prob_spec,
                     grid,
                     d_mpm->d_particleState,
                     d_mpm->d_particleState_preReloc);
  }
}

void
MPMICE::outputProblemSpec(ProblemSpecP& root_ps)
{
  d_mpm->outputProblemSpec(root_ps);
  d_ice->outputProblemSpec(root_ps);

  ProblemSpecP mpm_ps = root_ps->findBlock("MPM");
  mpm_ps->appendElement("testForNegTemps_mpm", d_testForNegTemps_mpm);

  // Global flags required by mpmice
  ProblemSpecP root      = root_ps->getRootNode();
  ProblemSpecP mpmice_ps = root->appendChild("MPMICE");
  mpmice_ps->appendElement("use_simple_equilibration_algorithm",
                           d_useSimpleEquilibrationPressure);
  mpmice_ps->appendElement("vol_frac_convergence_tolerance",
                           d_convergence_tolerance);

  for (auto& am : d_analysisModules) {
    am->outputProblemSpec(root_ps);
  }
}

void
MPMICE::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, cout_doing, "MPMICE::scheduleInitialize");

  d_mpm->scheduleInitialize(level, sched);
  d_ice->scheduleInitialize(level, sched);

  //__________________________________
  //  What isn't initialized in either ice or mpm
  Task* t = scinew Task(
    "MPMICE::actuallyInitialize", this, &MPMICE::actuallyInitialize);

  // Get the material subsets
  const MaterialSubset* ice_matls =
    d_materialManager->allMaterials("ICE")->getUnion();
  const MaterialSubset* mpm_matls =
    d_materialManager->allMaterials("MPM")->getUnion();

  t->requires(Task::NewDW, d_ice_labels->timeStepLabel);

  // These values are calculated for ICE materials in
  // d_ice->actuallyInitialize(...)
  //  so they are only needed for MPM
  t->computes(d_mpmice_labels->velocity_CCLabel, mpm_matls);
  t->computes(d_ice_labels->rho_CCLabel, mpm_matls);
  t->computes(d_ice_labels->temperature_CCLabel, mpm_matls);
  t->computes(d_ice_labels->specificVolume_CCLabel, mpm_matls);
  t->computes(d_ice_labels->speedSound_CCLabel, mpm_matls);
  t->computes(d_mpm_labels->heatRate_CCLabel, mpm_matls);

  // This is computed in d_ice->actuallyInitalize(...), and it is needed in
  //  MPMICE's actuallyInitialize()
  t->requires(
    Task::NewDW, d_ice_labels->vol_frac_CCLabel, ice_matls, Ghost::None, 0);

  if (d_switchCriteria) {
    d_switchCriteria->scheduleInitialize(level, sched);
  }

  // DataAnalysis
  for (auto& am : d_analysisModules) {
    am->scheduleInitialize(sched, level);
  }

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials());
}

void
MPMICE::actuallyInitialize(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  new_dw->get(timeStep, VarLabel::find(timeStep_name));
  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing actuallyInitialize ");

    // Output material indices
    if (patch->getID() == 0) {
      std::cout << "Materials Indicies:   MPM ["
                << *(d_materialManager->allMaterials("MPM")) << "] "
                << "ICE[" << *(d_materialManager->allMaterials("ICE")) << "]"
                << std::endl;

      std::cout << "Material Names:";
      size_t numAllMatls = d_materialManager->getNumMaterials();
      for (size_t m = 0; m < numAllMatls; m++) {
        Material* matl = d_materialManager->getMaterial(m);
        std::cout << " " << matl->getDWIndex() << ") " << matl->getName();
      }
      std::cout << "\n";
    }

    // Sum variable for testing that the volume fractions sum to 1
    CCVariable<double> vol_frac_sum;
    CCVariable<double> vol_frac_sum_mpm;
    new_dw->allocateTemporary(vol_frac_sum, patch);
    new_dw->allocateTemporary(vol_frac_sum_mpm, patch);
    vol_frac_sum.initialize(0.0);
    vol_frac_sum_mpm.initialize(0.0);

    //  Initialize CCVaribles for MPM Materials
    //  Even if mass = 0 in a cell you still need
    //  CC Variables defined.
    double junk{ -9.0 }, tmp{ 0.0 };
    double p_ref = d_ice->getRefPress();

    size_t numMPM_matls = d_materialManager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPM_matls; m++) {

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();

      CCVariable<double> sp_vol_CC, rho_CC, Temp_CC, speedSound;
      CCVariable<Vector> vel_CC;

      new_dw->allocateAndPut(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->allocateAndPut(rho_CC, d_ice_labels->rho_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        speedSound, d_ice_labels->speedSound_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        Temp_CC, d_mpmice_labels->temperature_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        vel_CC, d_mpmice_labels->velocity_CCLabel, indx, patch);

      CCVariable<double> rho_micro, vol_frac_CC;

      new_dw->allocateTemporary(rho_micro, patch);
      new_dw->allocateTemporary(vol_frac_CC, patch);

      CCVariable<double> heatFlux;
      new_dw->allocateAndPut(
        heatFlux, d_mpm_labels->heatRate_CCLabel, indx, patch);
      heatFlux.initialize(0.0);

      mpm_matl->initializeCCVariables(
        rho_micro, rho_CC, Temp_CC, vel_CC, vol_frac_CC, patch);

      setBC(rho_CC,
            "Density",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(rho_micro,
            "Density",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(Temp_CC,
            "Temperature",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(vel_CC,
            "Velocity",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);

      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c  = *iter;
        sp_vol_CC[c] = 1.0 / rho_micro[c];

        mpm_matl->getConstitutiveModel()->computePressEOSCM(
          rho_micro[c], junk, p_ref, junk, tmp, mpm_matl, Temp_CC[c]);
        speedSound[c] = std::sqrt(tmp);

        // sum volume fraction
        vol_frac_sum_mpm[c] += vol_frac_CC[c];
      }

      //__________________________________
      //    B U L L E T   P R O O F I N G
      IntVector neg_cell;
      std::ostringstream warn;
      if (!areAllValuesPositive(rho_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << indx << " cell "
             << neg_cell << " position: " << pt << " rho_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(Temp_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << indx << " cell "
             << neg_cell << " position: " << pt << " Temp_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(sp_vol_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << indx << " cell "
             << neg_cell << " position: " << pt << " sp_vol_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesNumbers(speedSound, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << indx << " cell "
             << neg_cell << " position: " << pt << " speedSound is nan\n";
        warn << "speedSound = " << speedSound[neg_cell]
             << " sp_vol_CC = " << sp_vol_CC[neg_cell]
             << " rho_micro = " << rho_micro[neg_cell]
             << " Temp_CC = " << Temp_CC[neg_cell] << std::endl;
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
    } // num_MPM_matls loop

    //___________________________________
    //   B U L L E T  P R O O F I N G
    // Verify volume fractions sum to 1.0
    // Loop through ICE materials to get their contribution to volume fraction
    size_t numICE_matls = d_materialManager->getNumMaterials("ICE");
    for (size_t m = 0; m < numICE_matls; m++) {
      constCCVariable<double> vol_frac;
      ICEMaterial* ice_matl =
        static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
      int indx = ice_matl->getDWIndex();

      // Get the Volume Fraction computed in ICE's actuallyInitialize(...)
      new_dw->get(
        vol_frac, d_ice_labels->vol_frac_CCLabel, indx, patch, Ghost::None, 0);

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        vol_frac_sum[c] += vol_frac[c];
      }
    } // num_ICE_matls loop

    double errorThresholdTop    = 1.0e0 + 1.0e-10;
    double errorThresholdBottom = 1.0e0 - 1.0e-10;
    bool wrongVolFrac           = false;
    for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;

      double tot_vol_frac = vol_frac_sum[c] + vol_frac_sum_mpm[c];
      if (!(tot_vol_frac <= errorThresholdTop &&
            tot_vol_frac >= errorThresholdBottom)) {

        wrongVolFrac = true;

        //--------------------------------------------------------------------------------------------------
        // Hack for coupling MPMICE and Smooth[Cyl/Sphere]GeomPiece:
        //   Check if the cell volume fracs sum to less than 1.  If that is
        //   true, then increase the volume fraction of the ICE material that
        //   occupies the largest volume fraction of the cell until the total
        //   volume fraction equals 1. Other cell centered quantities will also
        //   have to be changed accordingly.
        // Author: Biswajit Banerjee
        // Date: 9/26/2013
        //--------------------------------------------------------------------------------------------------
        std::vector<double> iceMaterialVolFrac;
        for (size_t mm = 0; mm < numICE_matls; ++mm) {
          constCCVariable<double> vol_frac;
          int materialIndex =
            d_materialManager->getMaterial("ICE", mm)->getDWIndex();
          new_dw->get(vol_frac,
                      d_ice_labels->vol_frac_CCLabel,
                      materialIndex,
                      patch,
                      Ghost::None,
                      0);
          iceMaterialVolFrac.push_back(vol_frac[c]);
        }

        for (unsigned int ii = 0; ii < iceMaterialVolFrac.size(); ++ii) {
          std::cout << "ICE material " << ii
                    << " vol frac = " << iceMaterialVolFrac[ii] << std::endl;
        }
        std::cout << "Max vol frac index = "
                  << std::max_element(iceMaterialVolFrac.begin(),
                                      iceMaterialVolFrac.end()) -
                       iceMaterialVolFrac.begin()
                  << std::endl;

        CCVariable<double> vol_frac_upd;
        int vol_frac_index = (int)(std::max_element(iceMaterialVolFrac.begin(),
                                                    iceMaterialVolFrac.end()) -
                                   iceMaterialVolFrac.begin());
        int ice_material_index =
          d_materialManager->getMaterial("ICE", vol_frac_index)->getDWIndex();
        new_dw->getModifiable(vol_frac_upd,
                              d_ice_labels->vol_frac_CCLabel,
                              ice_material_index,
                              patch,
                              Ghost::None,
                              0);
        double vol_frac_correction = 1.0 - tot_vol_frac;
        vol_frac_upd[c] += vol_frac_correction;
        std::cout << "Vol frac correction = " << vol_frac_correction
                  << " new vol frac = " << vol_frac_upd[c];
        std::cout << std::endl;

        // Now that the volume fraction has been corrected, correct the density
        constCCVariable<double> rho_micro;
        CCVariable<double> rho_CC_upd;
        new_dw->get(rho_micro,
                    d_ice_labels->rho_micro_CCLabel,
                    ice_material_index,
                    patch,
                    Ghost::None,
                    0);
        new_dw->getModifiable(rho_CC_upd,
                              d_ice_labels->rho_CCLabel,
                              ice_material_index,
                              patch,
                              Ghost::None,
                              0);
        std::cout << " old density = " << rho_CC_upd[c];
        rho_CC_upd[c] = rho_micro[c] * std::abs(vol_frac_upd[c]);
        std::cout << " new density = " << rho_CC_upd[c] << std::endl;

      } // end if check for volume fraction
    }   // end cell iterator for volume fraction

    // Redo the sum computation
    if (wrongVolFrac) {

      vol_frac_sum.initialize(0.0);
      int numICE_matls = d_materialManager->getNumMaterials("ICE");
      for (int m = 0; m < numICE_matls; m++) {
        constCCVariable<double> vol_frac;
        ICEMaterial* ice_matl =
          static_cast<ICEMaterial*>(d_materialManager->getMaterial("ICE", m));
        int indx = ice_matl->getDWIndex();

        // Get the Volume Fraction computed in ICE's actuallyInitialize(...)
        new_dw->get(vol_frac,
                    d_ice_labels->vol_frac_CCLabel,
                    indx,
                    patch,
                    Ghost::None,
                    0);

        for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          vol_frac_sum[c] += vol_frac[c];
        }
      } // num_ICE_matls loop

      for (auto iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;

        Point pt = patch->getCellPosition(c);
        if (!(vol_frac_sum[c] <= errorThresholdTop &&
              vol_frac_sum[c] >= errorThresholdBottom)) {
          std::ostringstream warn;
          warn << "ERROR MPMICE::actuallyInitialize cell " << c
               << " position: " << pt << "\n\n"
               << "volume fraction (" << std::setprecision(13)
               << vol_frac_sum[c] << ") does not sum to 1.0 +- 1e-10.\n"
               << "Verify that this region of the domain contains at least 1 "
                  "geometry object.  If you're using the optional\n"
               << "'volumeFraction' tags verify that they're correctly "
                  "specified.\n";
          throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
        }
      } // cell iterator for volume fraction
    }   // end if wrong volume frac
  }     // Patch loop
}

void
MPMICE::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, cout_doing, "MPMICE::scheduleInitialize");

  d_mpm->scheduleRestartInitialize(level, sched);
  d_ice->scheduleRestartInitialize(level, sched);

  for (auto& am : d_analysisModules) {
    am->scheduleRestartInitialize(sched, level);
  }
}

void
MPMICE::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  // Schedule computing the ICE stable timestep
  d_ice->scheduleComputeStableTimestep(level, sched);
  // MPM stable timestep is a by product of the CM
}

void
MPMICE::scheduleTimeAdvance(const LevelP& inlevel, SchedulerP& sched)
{
  // Only do scheduling on level 0 for lockstep AMR
  if (inlevel->getIndex() > 0 && isLockstepAMR()) {
    return;
  }

  // If we have a finer level, then assume that we are doing multilevel MPMICE
  // Otherwise, it is plain-ole MPMICE
  do_mlmpmice = false;
  if (inlevel->hasFinerLevel()) {
    do_mlmpmice = true;
  }

  const LevelP& mpm_level =
    do_mlmpmice
      ? inlevel->getGrid()->getLevel(inlevel->getGrid()->numLevels() - 1)
      : inlevel;

  const PatchSet* mpm_patches  = mpm_level->eachPatch();
  const MaterialSet* ice_matls = d_materialManager->allMaterials("ICE");
  const MaterialSet* mpm_matls = d_materialManager->allMaterials("MPM");
  const MaterialSet* all_matls = d_materialManager->allMaterials();
  MaterialSubset* press_matl   = d_ice->d_press_matl;
  MaterialSubset* one_matl     = d_ice->d_press_matl;

  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();
  const MaterialSubset* mpm_matls_sub = mpm_matls->getUnion();

  cout_doing
    << "---------------------------------------------------------Level ";
  if (do_mlmpmice) {
    cout_doing << inlevel->getIndex() << " (ICE) " << mpm_level->getIndex()
               << " (MPM)" << std::endl;
    ;
  } else {
    cout_doing << inlevel->getIndex() << std::endl;
  }

  // Scheduling
  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level = inlevel->getGrid()->getLevel(l);
    d_ice->scheduleComputeThermoTransportProperties(
      sched, ice_level, ice_matls);

    d_ice->scheduleMaxMach_on_Lodi_BC_Faces(sched, ice_level, ice_matls);
  }

  // diagnostic task
  // d_mpm->scheduleTotalParticleCount(sched, mpm_patches, mpm_matls);

  d_mpm->scheduleComputeParticleBodyForce(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleApplyExternalLoads(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeCurrentParticleSize(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleInterpolateParticlesToGrid(sched, mpm_patches, mpm_matls);
  d_mpm->d_heatConductionTasks->scheduleComputeHeatExchange(
    sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeNormals(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleMomentumExchangeInterpolated(sched, mpm_patches, mpm_matls);

  // schedule the interpolation of mass and volume to the cell centers
  scheduleInterpolateNCToCC_0(sched, mpm_patches, one_matl, mpm_matls);

  // do coarsens in reverse order, and before the other tasks
  if (do_mlmpmice) {
    for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
      const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();
      scheduleCoarsenCC_0(sched, ice_patches, mpm_matls);
      scheduleCoarsenNCMass(sched, ice_patches, mpm_matls);
    }
  }

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
    scheduleComputePressure(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, press_matl, all_matls);

    d_ice->scheduleComputeTempFC(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, all_matls);
    d_ice->scheduleComputeModelSources(sched, ice_level, all_matls);
    d_ice->scheduleUpdateVolumeFraction(
      sched, ice_level, press_matl, all_matls);
    d_ice->scheduleComputeVel_FC(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, press_matl, all_matls);

    d_ice->d_exchModel->sched_PreExchangeTasks(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, all_matls);

    d_ice->d_exchModel->sched_AddExch_VelFC(sched,
                                            ice_patches,
                                            ice_matls_sub,
                                            mpm_matls_sub,
                                            all_matls,
                                            d_ice->d_BC_globalVars.get(),
                                            false);
  }

  if (d_ice->d_impICE) { //  I M P L I C I T, won't work with AMR yet
    // we should use the AMR multi-level pressure solve
    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
      const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();

      d_ice->scheduleSetupRHS(
        sched, ice_patches, one_matl, all_matls, false, "computes");
      d_ice->scheduleCompute_maxRHS(sched, ice_level, one_matl, all_matls);
      d_ice->scheduleImplicitPressureSolve(sched,
                                           ice_level,
                                           ice_patches,
                                           one_matl,
                                           press_matl,
                                           ice_matls_sub,
                                           mpm_matls_sub,
                                           all_matls);
      d_ice->scheduleComputeDel_P(
        sched, ice_level, ice_patches, one_matl, press_matl, all_matls);
    }
  } //  IMPLICIT AND EXPLICIT

  if (!(d_ice->d_impICE)) { //  E X P L I C I T
    for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
      const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();
      d_ice->scheduleComputeDelPressAndUpdatePressCC(sched,
                                                     ice_patches,
                                                     press_matl,
                                                     ice_matls_sub,
                                                     mpm_matls_sub,
                                                     all_matls);
    }
  }

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
    d_ice->scheduleComputePressFC(sched, ice_patches, press_matl, all_matls);
    d_ice->scheduleVelTau_CC(sched, ice_patches, ice_matls);
    d_ice->scheduleViscousShearStress(sched, ice_patches, ice_matls);
    d_ice->scheduleAccumulateMomentumSourceSinks(
      sched, ice_patches, press_matl, ice_matls_sub, mpm_matls_sub, all_matls);
    d_ice->scheduleAccumulateEnergySourceSinks(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, press_matl, all_matls);
  }

  if (!d_rigidMPM) {
    scheduleInterpolatePressCCToPressNC(
      sched, mpm_patches, press_matl, mpm_matls);
    scheduleInterpolatePAndGradP(
      sched, mpm_patches, press_matl, one_matl, mpm_matls_sub, mpm_matls);
  }

  d_mpm->scheduleComputeInternalForce(sched, mpm_patches, mpm_matls);
  d_mpm->d_heatConductionTasks->scheduleComputeInternalHeatRate(
    sched, mpm_patches, mpm_matls);
  d_mpm->d_heatConductionTasks->scheduleComputeNodalHeatFlux(
    sched, mpm_patches, mpm_matls);
  d_mpm->d_heatConductionTasks->scheduleSolveHeatEquations(
    sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeAndIntegrateAcceleration(sched, mpm_patches, mpm_matls);
  d_mpm->d_heatConductionTasks->scheduleIntegrateTemperatureRate(
    sched, mpm_patches, mpm_matls);

  scheduleComputeLagrangianValuesMPM(sched, mpm_patches, one_matl, mpm_matls);

  // do coarsens in reverse order, and before the other tasks
  if (do_mlmpmice) {
    for (int l = inlevel->getGrid()->numLevels() - 2; l >= 0; l--) {
      const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
      const PatchSet* ice_patches = ice_level->eachPatch();
      scheduleCoarsenLagrangianValuesMPM(sched, ice_patches, mpm_matls);
    }
  }

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
    d_ice->scheduleComputeLagrangianValues(sched, ice_patches, ice_matls);
    d_ice->d_exchModel->sched_AddExch_Vel_Temp_CC(sched,
                                                  ice_patches,
                                                  ice_matls_sub,
                                                  mpm_matls_sub,
                                                  all_matls,
                                                  d_ice->d_BC_globalVars.get());
    d_ice->scheduleComputeLagrangianSpecificVolume(
      sched, ice_patches, ice_matls_sub, mpm_matls_sub, press_matl, all_matls);
    d_ice->scheduleComputeLagrangian_Transported_Vars(
      sched, ice_patches, ice_matls);
  }

  scheduleComputeCCVelAndTempRates(sched, mpm_patches, mpm_matls);
  scheduleInterpolateCCToNC(sched, mpm_patches, mpm_matls);

  d_mpm->scheduleMomentumExchangeIntegrated(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleSetGridBoundaryConditions(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleInterpolateToParticlesAndUpdate(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeDeformationGradient(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeStressTensor(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleComputeBasicDamage(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleUpdateErosionParameter(sched, mpm_patches, mpm_matls);
  d_mpm->scheduleFinalParticleUpdate(sched, mpm_patches, mpm_matls);
  // d_mpm->scheduleApplyExternalLoads(          sched, mpm_patches, mpm_matls);

  for (int l = 0; l < inlevel->getGrid()->numLevels(); l++) {
    const LevelP& ice_level     = inlevel->getGrid()->getLevel(l);
    const PatchSet* ice_patches = ice_level->eachPatch();
    d_ice->scheduleAdvectAndAdvanceInTime(
      sched, ice_patches, ice_matls_sub, ice_matls);
    d_ice->scheduleConservedtoPrimitive_Vars(
      sched, ice_patches, ice_matls_sub, ice_matls, "afterAdvection");
  }
} // end scheduleTimeAdvance()

void
MPMICE::scheduleInterpolateNCToCC_0(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSubset* one_matl,
                                    const MaterialSet* mpm_matls)
{
  if (d_mpm->d_mpm_flags->doMPMOnLevel(
        getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {

    printSchedule(patches, cout_doing, "MPMICE::scheduleInterpolateNCToCC_0");

    /* interpolateNCToCC */
    Task* t = scinew Task(
      "MPMICE::interpolateNCToCC_0", this, &MPMICE::interpolateNCToCC_0);

    const MaterialSubset* mss = mpm_matls->getUnion();
    t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::AroundCells, 1);
    t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::AroundCells, 1);
    t->requires(
      Task::NewDW, d_mpm_labels->gVelocityBCLabel, Ghost::AroundCells, 1);
    t->requires(
      Task::NewDW, d_mpm_labels->gTemperatureLabel, Ghost::AroundCells, 1);
    t->requires(
      Task::NewDW, d_mpm_labels->gSpecificVolumeLabel, Ghost::AroundCells, 1);
    t->requires(Task::OldDW,
                d_mpm_labels->NC_CCweightLabel,
                one_matl,
                Ghost::AroundCells,
                1);
    t->requires(
      Task::OldDW, d_ice_labels->specificVolume_CCLabel, Ghost::None, 0);
    t->requires(
      Task::OldDW, d_mpmice_labels->temperature_CCLabel, Ghost::None, 0);

    t->computes(d_mpmice_labels->cMassLabel);
    t->computes(d_mpmice_labels->velocity_CCLabel);
    t->computes(d_mpmice_labels->temperature_CCLabel);
    t->computes(d_ice_labels->specificVolume_CCLabel, mss);
    t->computes(d_ice_labels->rho_CCLabel, mss);

    sched->addTask(t, patches, mpm_matls);
  }
}

void
MPMICE::interpolateNCToCC_0(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, VarLabel::find(timeStep_name));
  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing interpolateNCToCC_0");

    Vector dx       = patch->dCell();
    double cell_vol = dx.x() * dx.y() * dx.z();

    Ghost::GhostType gac = Ghost::AroundCells;
    Ghost::GhostType gn  = Ghost::None;

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);

    size_t numMatls = d_materialManager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();

      // Create arrays for the grid data
      constNCVariable<double> gMass, gVolume, gTemperature, gSp_vol;
      constNCVariable<Vector> gVelocity;
      CCVariable<double> cmass, Temp_CC, sp_vol_CC, rho_CC;
      CCVariable<Vector> vel_CC;
      constCCVariable<double> Temp_CC_ice, sp_vol_CC_ice;

      new_dw->allocateAndPut(cmass, d_mpmice_labels->cMassLabel, indx, patch);
      new_dw->allocateAndPut(
        vel_CC, d_mpmice_labels->velocity_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        Temp_CC, d_mpmice_labels->temperature_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->allocateAndPut(rho_CC, d_ice_labels->rho_CCLabel, indx, patch);

      double very_small_mass = d_TINY_RHO * cell_vol;
      cmass.initialize(very_small_mass);

      new_dw->get(gMass, d_mpm_labels->gMassLabel, indx, patch, gac, 1);
      new_dw->get(gVolume, d_mpm_labels->gVolumeLabel, indx, patch, gac, 1);
      new_dw->get(
        gVelocity, d_mpm_labels->gVelocityBCLabel, indx, patch, gac, 1);
      new_dw->get(
        gTemperature, d_mpm_labels->gTemperatureLabel, indx, patch, gac, 1);
      new_dw->get(
        gSp_vol, d_mpm_labels->gSpecificVolumeLabel, indx, patch, gac, 1);
      old_dw->get(sp_vol_CC_ice,
                  d_ice_labels->specificVolume_CCLabel,
                  indx,
                  patch,
                  gn,
                  0);
      old_dw->get(
        Temp_CC_ice, d_mpmice_labels->temperature_CCLabel, indx, patch, gn, 0);

      IntVector nodeIdx[8];

      //  compute CC Variables
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        patch->findNodesFromCell(*iter, nodeIdx);

        double Temp_CC_mpm = 0.0;
        double sp_vol_mpm  = 0.0;
        Vector vel_CC_mpm  = Vector(0.0, 0.0, 0.0);

        for (int in = 0; in < 8; in++) {
          double NC_CCw_mass = NC_CCweight[nodeIdx[in]] * gMass[nodeIdx[in]];
          cmass[c] += NC_CCw_mass;
          sp_vol_mpm += gSp_vol[nodeIdx[in]] * NC_CCw_mass;
          vel_CC_mpm += gVelocity[nodeIdx[in]] * NC_CCw_mass;
          Temp_CC_mpm += gTemperature[nodeIdx[in]] * NC_CCw_mass;
        }
        double inv_cmass = 1.0 / cmass[c];
        vel_CC_mpm *= inv_cmass;
        Temp_CC_mpm *= inv_cmass;
        sp_vol_mpm *= inv_cmass;

        //__________________________________
        // set *_CC = to either vel/Temp_CC_mpm or some safe values
        // depending upon if there is cmass.  You need
        // a well defined vel/temp_CC even if there isn't any mass
        // If you change this you must also change
        // MPMICE::computeLagrangianValuesMPM
        double one_or_zero = (cmass[c] - very_small_mass) / cmass[c];

        Temp_CC[c] =
          (1.0 - one_or_zero) * Temp_CC_ice[c] + one_or_zero * Temp_CC_mpm;
        sp_vol_CC[c] =
          (1.0 - one_or_zero) * sp_vol_CC_ice[c] + one_or_zero * sp_vol_mpm;

        vel_CC[c] = vel_CC_mpm;
        rho_CC[c] = cmass[c] / cell_vol;
      }

      //  Set BC's
      setBC(Temp_CC,
            "Temperature",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(rho_CC,
            "Density",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(vel_CC,
            "Velocity",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      //  Set if symmetric Boundary conditions
      setBC(cmass,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(sp_vol_CC,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);

      //---- B U L L E T   P R O O F I N G------
      // ignore BP if timestep restart has already been requested
      IntVector neg_cell;
      std::ostringstream warn;
      bool rts = new_dw->recomputeTimeStep();

      int L = getLevel(patches)->getIndex();
      if (d_testForNegTemps_mpm) {
        if (!areAllValuesPositive(Temp_CC, neg_cell) && !rts) {
          warn << "ERROR MPMICE:(" << L << "):interpolateNCToCC_0, mat " << indx
               << " cell " << neg_cell << " Temp_CC " << Temp_CC[neg_cell]
               << "\n ";
          throw InvalidValue(warn.str(), __FILE__, __LINE__);
        }
      }
      if (!areAllValuesPositive(rho_CC, neg_cell) && !rts) {
        warn << "ERROR MPMICE:(" << L << "):interpolateNCToCC_0, mat " << indx
             << " cell " << neg_cell << " rho_CC " << rho_CC[neg_cell] << "\n ";
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(sp_vol_CC, neg_cell) && !rts) {
        warn << "ERROR MPMICE:(" << L << "):interpolateNCToCC_0, mat " << indx
             << " cell " << neg_cell << " sp_vol_CC " << sp_vol_CC[neg_cell]
             << "\n ";
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
    }
  } // patches
}

void
MPMICE::scheduleCoarsenCC_0(SchedulerP& sched,
                            const PatchSet* patches,
                            const MaterialSet* mpm_matls)
{
  printSchedule(patches, cout_doing, "MPMICE::scheduleCoarsenCC_0");

  bool modifies = false;

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_mpmice_labels->cMassLabel,
                            1.9531e-15,
                            modifies,
                            "sum");

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_mpmice_labels->temperature_CCLabel,
                            0.,
                            modifies,
                            "massWeighted");

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_mpmice_labels->velocity_CCLabel,
                            Vector(0, 0, 0),
                            modifies,
                            "massWeighted");

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->specificVolume_CCLabel,
                            0.8479864471,
                            modifies,
                            "massWeighted");

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->rho_CCLabel,
                            1.e-12,
                            modifies,
                            "std");
}

template<typename T>
void
MPMICE::scheduleCoarsenVariableCC(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls,
                                  const VarLabel* variable,
                                  T defaultValue,
                                  bool modifies,
                                  const string& coarsenMethod)
{
  auto func = &MPMICE::coarsenVariableCC<T>;
  std::ostringstream taskName;

  taskName << "MPMICE::coarsenVariableCC(" << variable->getName()
           << (modifies ? " modified" : "") << ")";

  Task* t = scinew Task(taskName.str().c_str(),
                        this,
                        func,
                        variable,
                        defaultValue,
                        modifies,
                        coarsenMethod);

  Ghost::GhostType gn         = Ghost::None;
  Task::MaterialDomainSpec ND = Task::NormalDomain;

  t->requires(Task::OldDW, d_ice_labels->timeStepLabel);
  t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND, gn, 0);

  if (coarsenMethod == "massWeighted") {
    t->requires(Task::NewDW,
                d_mpmice_labels->cMassLabel,
                0,
                Task::FineLevel,
                0,
                ND,
                gn,
                0);
  }

  if (modifies) {
    t->modifies(variable);
  } else {
    t->computes(variable);
  }
  sched->addTask(t, patches, matls);
}

template<typename T>
void
MPMICE::coarsenVariableCC(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw,
                          const VarLabel* variable,
                          T defaultValue,
                          bool modifies,
                          string coarsenMethod)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimeStep = (timeStep > 0);

  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  IntVector refineRatio(fineLevel->getRefinementRatio());
  double ratio = 1. / (refineRatio.x() * refineRatio.y() * refineRatio.z());

  for (int p = 0; p < patches->size(); p++) {
    const Patch* coarsePatch = patches->get(p);
    std::ostringstream message;
    message << "Doing CoarsenVariableCC (" << variable->getName() << ")\t\t\t";
    printTask(patches, coarsePatch, cout_doing, message.str());

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      CCVariable<T> coarse_q_CC;
      if (modifies) {
        new_dw->getModifiable(coarse_q_CC, variable, indx, coarsePatch);
      } else {
        new_dw->allocateAndPut(coarse_q_CC, variable, indx, coarsePatch);
      }
      coarse_q_CC.initialize(defaultValue);

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);
      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        IntVector cl, ch, fl, fh;
        getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);
        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }

        constCCVariable<T> fine_q_CC;
        new_dw->getRegion(fine_q_CC, variable, indx, fineLevel, fl, fh, false);

        //__________________________________
        //  call the coarsening function
        ASSERT((coarsenMethod == "std" || coarsenMethod == "sum" ||
                coarsenMethod == "massWeighted"));
        if (coarsenMethod == "std") {
          AMRCoarsenRefine::coarsenDriver_std(cl,
                                              ch,
                                              fl,
                                              fh,
                                              refineRatio,
                                              ratio,
                                              coarseLevel,
                                              fine_q_CC,
                                              coarse_q_CC);
        }
        if (coarsenMethod == "sum") {
          ratio = 1.0;
          AMRCoarsenRefine::coarsenDriver_std(cl,
                                              ch,
                                              fl,
                                              fh,
                                              refineRatio,
                                              ratio,
                                              coarseLevel,
                                              fine_q_CC,
                                              coarse_q_CC);
        }
        if (coarsenMethod == "massWeighted") {
          constCCVariable<double> cMass;
          new_dw->getRegion(
            cMass, d_mpmice_labels->cMassLabel, indx, fineLevel, fl, fh, false);

          AMRCoarsenRefine::coarsenDriver_massWeighted(cl,
                                                       ch,
                                                       fl,
                                                       fh,
                                                       refineRatio,
                                                       coarseLevel,
                                                       cMass,
                                                       fine_q_CC,
                                                       coarse_q_CC);
        }
      } // fine patches

      // Set BCs on coarsened data.  This sucks--Steve
      if (variable->getName() == "temp_CC") {
        setBC(coarse_q_CC,
              "Temperature",
              coarsePatch,
              d_materialManager,
              indx,
              new_dw,
              isNotInitialTimeStep);
      } else if (variable->getName() == "rho_CC") {
        setBC(coarse_q_CC,
              "Density",
              coarsePatch,
              d_materialManager,
              indx,
              new_dw,
              isNotInitialTimeStep);
      } else if (variable->getName() == "vel_CC") {
        setBC(coarse_q_CC,
              "Velocity",
              coarsePatch,
              d_materialManager,
              indx,
              new_dw,
              isNotInitialTimeStep);
      } else if (variable->getName() == "c.mass" ||
                 variable->getName() == "sp_vol_CC" ||
                 variable->getName() == "mom_L_CC" ||
                 variable->getName() == "int_eng_L_CC") {
        setBC(coarse_q_CC,
              "set_if_sym_BC",
              coarsePatch,
              d_materialManager,
              indx,
              new_dw,
              isNotInitialTimeStep);
      }
    } // matls
  }   // coarse level
}

void
MPMICE::scheduleCoarsenNCMass(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* mpm_matls)
{
  printSchedule(patches, cout_doing, "MPMICE::scheduleCoarsenNCMass");

  bool modifies = false;

  scheduleCoarsenVariableNC(sched,
                            patches,
                            mpm_matls,
                            d_mpm_labels->gMassLabel,
                            1.e-200,
                            modifies,
                            "sum");
}

template<typename T>
void
MPMICE::scheduleCoarsenVariableNC(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls,
                                  const VarLabel* variable,
                                  T defaultValue,
                                  bool modifies,
                                  string coarsenMethod)
{
  auto func = &MPMICE::coarsenVariableNC<T>;

  std::ostringstream taskName;

  taskName << "MPMICE::coarsenVariableNC(" << variable->getName()
           << (modifies ? " modified" : "") << ")";

  Task* t = scinew Task(taskName.str().c_str(),
                        this,
                        func,
                        variable,
                        defaultValue,
                        modifies,
                        coarsenMethod);

  // Ghost::GhostType  gn = Ghost::None;
  Ghost::GhostType gan        = Ghost::AroundNodes;
  Task::MaterialDomainSpec ND = Task::NormalDomain;

  const LevelP fineLevel = getLevel(patches)->getFinerLevel();
  IntVector refineRatio(fineLevel->getRefinementRatio());
  int ghost = std::max(refineRatio.x(), refineRatio.y());
  ghost     = std::max(ghost, refineRatio.z());

  t->requires(Task::NewDW, variable, 0, Task::FineLevel, 0, ND, gan, ghost);

  if (modifies) {
    t->modifies(variable);
  } else {
    t->computes(variable);
  }
  sched->addTask(t, patches, matls);
}

template<typename T>
void
MPMICE::coarsenVariableNC(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset* matls,
                          DataWarehouse*,
                          DataWarehouse* new_dw,
                          const VarLabel* variable,
                          T defaultValue,
                          bool modifies,
                          string coarsenMethod)
{
  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();

  IntVector refineRatio(fineLevel->getRefinementRatio());
  double ratio = 1. / (refineRatio.x() * refineRatio.y() * refineRatio.z());

  for (int p = 0; p < patches->size(); p++) {
    const Patch* coarsePatch = patches->get(p);
    std::ostringstream message;
    message << "Doing CoarsenVariableNC (" << variable->getName() << ")\t\t\t";
    printTask(patches, coarsePatch, cout_doing, message.str());

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      NCVariable<T> coarse_q_NC;
      if (modifies) {
        new_dw->getModifiable(coarse_q_NC, variable, indx, coarsePatch);
      } else {
        new_dw->allocateAndPut(coarse_q_NC, variable, indx, coarsePatch);
      }
      coarse_q_NC.initialize(defaultValue);

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);
      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        IntVector cl, ch, fl, fh;

        IntVector padding(
          refineRatio.x() / 2, refineRatio.y() / 2, refineRatio.z() / 2);
        getFineLevelRangeNodes(coarsePatch, finePatch, cl, ch, fl, fh, padding);

        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }

        constNCVariable<T> fine_q_NC;
        new_dw->getRegion(fine_q_NC, variable, indx, fineLevel, fl, fh, false);

        //__________________________________
        //  call the coarsening function
        ASSERT(coarsenMethod == "sum");
        if (coarsenMethod == "sum") {
          coarsenDriver_stdNC(
            cl, ch, refineRatio, ratio, coarseLevel, fine_q_NC, coarse_q_NC);
        }
      } // fine patches
    }   // matls
  }     // coarse level
}

template<typename T>
void
MPMICE::coarsenDriver_stdNC(IntVector cl,
                            IntVector ch,
                            IntVector refinementRatio,
                            double ratio,
                            const Level* coarseLevel,
                            constNCVariable<T>& fine_q_NC,
                            NCVariable<T>& coarse_q_NC)
{
  T zero(0.0);

  // iterate over coarse level cells
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  Vector DX              = coarseLevel->dCell();
  IntVector range(
    refinementRatio.x() / 2, refinementRatio.y() / 2, refinementRatio.z() / 2);

  IntVector varLow  = fine_q_NC.getLowIndex();
  IntVector varHigh = fine_q_NC.getHighIndex();

  for (NodeIterator iter(cl, ch); !iter.done(); iter++) {
    IntVector c        = *iter;
    IntVector fineNode = coarseLevel->mapNodeToFiner(c);
    Point coarseLoc    = coarseLevel->getNodePosition(c);

    IntVector start = Max(fineNode - range, varLow);
    IntVector end   = Min(fineNode + range, varHigh);

    // for each coarse level cell iterate over the fine level cells
    T q_NC_tmp(zero);

    for (NodeIterator inner(start, end); !inner.done(); inner++) {
      IntVector fc(*inner);
      Point fineLoc  = fineLevel->getNodePosition(fc);
      Vector C2F     = fineLoc - coarseLoc;
      Vector Vweight = C2F / DX;
      double weight  = (1. - std::abs(Vweight.x())) *
                      (1. - std::abs(Vweight.y())) *
                      (1. - std::abs(Vweight.z()));
      q_NC_tmp += fine_q_NC[fc] * weight;
    }
    coarse_q_NC[c] = q_NC_tmp;
  }
}

/*_____________________________________________________________________
  Function~  MPMICE::scheduleComputePressure--
  Note:  Handles both Rate and Equilibration form of solution technique
  _____________________________________________________________________*/
void
MPMICE::scheduleComputePressure(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSubset* ice_matls,
                                const MaterialSubset* mpm_matls,
                                const MaterialSubset* press_matl,
                                const MaterialSet* all_matls)
{
  printSchedule(
    patches, cout_doing, "MPMICE::scheduleComputeEquilibrationPressure");

  Task* t = scinew Task("MPMICE::computeEquilibrationPressure",
                        this,
                        &MPMICE::computeEquilibrationPressure,
                        press_matl);

  t->requires(Task::OldDW, d_ice_labels->timeStepLabel);
  t->requires(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));

  // I C E
  Ghost::GhostType gn = Ghost::None;

  t->requires(Task::OldDW, d_ice_labels->temperature_CCLabel, ice_matls, gn);
  t->requires(Task::OldDW, d_ice_labels->rho_CCLabel, ice_matls, gn);
  t->requires(Task::OldDW, d_ice_labels->specificVolume_CCLabel, ice_matls, gn);
  t->requires(Task::NewDW, d_ice_labels->specific_heatLabel, ice_matls, gn);
  t->requires(Task::NewDW, d_ice_labels->gammaLabel, ice_matls, gn);

  // M P M
  t->requires(Task::NewDW, d_mpmice_labels->temperature_CCLabel, mpm_matls, gn);
  t->requires(Task::NewDW, d_ice_labels->rho_CCLabel, mpm_matls, gn);
  t->requires(Task::NewDW, d_ice_labels->specificVolume_CCLabel, mpm_matls, gn);

  t->requires(Task::OldDW, d_ice_labels->press_CCLabel, press_matl, gn);
  t->requires(Task::OldDW, d_ice_labels->velocity_CCLabel, ice_matls, gn);
  t->requires(Task::NewDW, d_mpmice_labels->velocity_CCLabel, mpm_matls, gn);

  computesRequires_CustomBCs(t,
                             "EqPress",
                             d_ice->d_ice_labels.get(),
                             ice_matls,
                             d_ice->d_BC_globalVars.get());

  //  A L L _ M A T L S
  t->computes(d_ice_labels->f_theta_CCLabel);
  t->computes(d_ice_labels->compressibilityLabel, ice_matls);
  t->computes(d_ice_labels->compressibilityLabel, mpm_matls);

  t->computes(d_ice_labels->speedSound_CCLabel);
  t->computes(d_ice_labels->vol_frac_CCLabel);
  t->computes(d_ice_labels->sumKappaLabel, press_matl);
  t->computes(d_ice_labels->TMV_CCLabel, press_matl);
  t->computes(d_ice_labels->press_equil_CCLabel, press_matl);
  t->computes(d_ice_labels->sum_imp_delPLabel,
              press_matl); // needed by implicit ICE

  t->modifies(d_ice_labels->specificVolume_CCLabel, mpm_matls);
  t->computes(d_ice_labels->specificVolume_CCLabel, ice_matls);
  t->computes(d_ice_labels->rho_CCLabel, ice_matls);

  sched->addTask(t, patches, all_matls);
}

/* ---------------------------------------------------------------------
   Function~  MPMICE::computeEquilibrationPressure--
   Reference: Flow of Interpenetrating Material Phases, J. Comp, Phys
   18, 440-464, 1975, see the equilibration section

   Steps
   ----------------
   - Compute rho_micro_CC, SpeedSound, vol_frac for ALL matls

   For each cell
   _ WHILE LOOP(convergence, max_iterations)
   - compute the pressure and dp_drho from the EOS of each material.
   - Compute delta Pressure
   - Compute delta volume fraction and update the
   volume fraction and the celldensity.
   - Test for convergence of delta pressure and delta volume fraction
   - END WHILE LOOP
   - bulletproofing
   end

   Note:  The nomenclature follows the reference.
   This is identical to  ICE::computeEquilibrationPressure except
   we now include EOS for MPM matls.
   _____________________________________________________________________*/
void
MPMICE::computeEquilibrationPressure(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw,
                                     const MaterialSubset* press_matl)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, d_ice_labels->timeStepLabel);
  bool isNotInitialTimeStep = (timeStep > 0);

  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeEquilibrationPressure");

    // double    converg_coeff = 100.;
    // double    convergence_crit = converg_coeff * DBL_EPSILON;
    double convergence_crit = d_convergence_tolerance;
    double c_2;
    double press_ref = d_ice->getRefPress();
    int numICEMatls  = d_materialManager->getNumMaterials("ICE");
    int numMPMMatls  = d_materialManager->getNumMaterials("MPM");
    int numALLMatls  = numICEMatls + numMPMMatls;

    Vector dx       = patch->dCell();
    double cell_vol = dx.x() * dx.y() * dx.z();

    std::vector<double> press_eos(numALLMatls);
    std::vector<double> dp_drho(numALLMatls), dp_de(numALLMatls);
    std::vector<double> mat_volume(numALLMatls);

    std::vector<CCVariable<double>> vol_frac(numALLMatls);
    std::vector<CCVariable<double>> rho_micro(numALLMatls);
    std::vector<CCVariable<double>> rho_CC_new(numALLMatls);
    std::vector<CCVariable<double>> speedSound(numALLMatls);
    std::vector<CCVariable<double>> sp_vol_new(numALLMatls);
    std::vector<CCVariable<double>> f_theta(numALLMatls);
    std::vector<CCVariable<double>> kappa(numALLMatls);

    std::vector<constCCVariable<double>> placeHolder(0);
    std::vector<constCCVariable<double>> cv(numALLMatls);
    std::vector<constCCVariable<double>> gamma(numALLMatls);
    std::vector<constCCVariable<double>> sp_vol_CC(numALLMatls);
    std::vector<constCCVariable<double>> Temp(numALLMatls);
    std::vector<constCCVariable<double>> rho_CC_old(numALLMatls);
    std::vector<constCCVariable<double>> mass_CC(numALLMatls);
    std::vector<constCCVariable<Vector>> vel_CC(numALLMatls);

    constCCVariable<double> press;
    CCVariable<double> press_new, delPress_tmp, sumKappa, TMV_CC;
    CCVariable<double> sum_imp_delP;
    CCVariable<int> nIterations;

    Ghost::GhostType gn = Ghost::None;

    //__________________________________
    old_dw->get(press, d_ice_labels->press_CCLabel, 0, patch, gn, 0);
    new_dw->allocateAndPut(
      press_new, d_ice_labels->press_equil_CCLabel, 0, patch);
    new_dw->allocateAndPut(TMV_CC, d_ice_labels->TMV_CCLabel, 0, patch);
    new_dw->allocateAndPut(sumKappa, d_ice_labels->sumKappaLabel, 0, patch);
    new_dw->allocateAndPut(
      sum_imp_delP, d_ice_labels->sum_imp_delPLabel, 0, patch);
    new_dw->allocateAndPut(
      nIterations, d_ice_labels->eq_press_itersLabel, 0, patch);
    new_dw->allocateTemporary(delPress_tmp, patch);

    sum_imp_delP.initialize(0.0);
    nIterations.initialize(0.0);

    std::vector<MPMMaterial*> mpm_matl(numALLMatls);
    std::vector<ICEMaterial*> ice_matl(numALLMatls);
    for (int m = 0; m < numALLMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      ice_matl[m]    = dynamic_cast<ICEMaterial*>(matl);
      mpm_matl[m]    = dynamic_cast<MPMMaterial*>(matl);
    }

    for (int m = 0; m < numALLMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();

      if (ice_matl[m]) { // I C E
        old_dw->get(
          Temp[m], d_ice_labels->temperature_CCLabel, indx, patch, gn, 0);
        old_dw->get(
          rho_CC_old[m], d_ice_labels->rho_CCLabel, indx, patch, gn, 0);
        old_dw->get(sp_vol_CC[m],
                    d_ice_labels->specificVolume_CCLabel,
                    indx,
                    patch,
                    gn,
                    0);
        old_dw->get(
          vel_CC[m], d_ice_labels->velocity_CCLabel, indx, patch, gn, 0);
        new_dw->get(
          cv[m], d_ice_labels->specific_heatLabel, indx, patch, gn, 0);
        new_dw->get(gamma[m], d_ice_labels->gammaLabel, indx, patch, gn, 0);

        new_dw->allocateAndPut(
          rho_CC_new[m], d_ice_labels->rho_CCLabel, indx, patch);
      }

      if (mpm_matl[m]) { // M P M
        new_dw->get(
          Temp[m], d_mpmice_labels->temperature_CCLabel, indx, patch, gn, 0);
        new_dw->get(
          vel_CC[m], d_mpmice_labels->velocity_CCLabel, indx, patch, gn, 0);
        new_dw->get(sp_vol_CC[m],
                    d_ice_labels->specificVolume_CCLabel,
                    indx,
                    patch,
                    gn,
                    0);
        new_dw->get(
          rho_CC_old[m], d_ice_labels->rho_CCLabel, indx, patch, gn, 0);

        new_dw->allocateTemporary(rho_CC_new[m], patch);
      }

      new_dw->allocateTemporary(rho_micro[m], patch);
      new_dw->allocateAndPut(
        vol_frac[m], d_ice_labels->vol_frac_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        f_theta[m], d_ice_labels->f_theta_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        speedSound[m], d_ice_labels->speedSound_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        sp_vol_new[m], d_ice_labels->specificVolume_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        kappa[m], d_ice_labels->compressibilityLabel, indx, patch);
    }

    press_new.copyData(press);

    //__________________________________
    // Compute rho_micro, speedSound, volfrac, rho_CC
    // see Docs/MPMICE.txt for explanation of why we ONlY
    // use eos evaulations for rho_micro_mpm
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {

      const IntVector& c   = *iter;
      double total_mat_vol = 0.0;

      for (int m = 0; m < numALLMatls; m++) {
        if (ice_matl[m]) { // I C E
          rho_micro[m][c] = 1.0 / sp_vol_CC[m][c];
        } else if (mpm_matl[m]) { //  M P M
          rho_micro[m][c] =
            mpm_matl[m]->getConstitutiveModel()->computeRhoMicroCM(
              press_new[c],
              press_ref,
              mpm_matl[m],
              Temp[m][c],
              1.0 / sp_vol_CC[m][c]);
        }
        mat_volume[m] = (rho_CC_old[m][c] * cell_vol) / rho_micro[m][c];
        total_mat_vol += mat_volume[m];
      } // numAllMatls loop

      TMV_CC[c] = total_mat_vol;

      for (int m = 0; m < numALLMatls; m++) {
        vol_frac[m][c]   = mat_volume[m] / total_mat_vol;
        rho_CC_new[m][c] = vol_frac[m][c] * rho_micro[m][c];
      }
    } // cell iterator

    //______________________________________________________________________
    // Done with preliminary calcs, now loop over every cell
    int count{ 0 }, test_max_iter{ 0 };
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      const IntVector& c = *iter;
      double delPress    = 0.;
      bool converged     = false;
      double sum{ 0.0 };
      count = 0;
      std::vector<EqPress_dbg> dbgEqPress;

      if (d_useSimpleEquilibrationPressure) {

        // Simplified newton iterations
        double p_new = press_new[c];
        while (count < d_ice->d_max_iter_equilibration && converged == false) {

          ++count;

          double A = 0.0;
          double B = 0.0;

          for (int m = 0; m < numALLMatls; m++) {

            double rho_m_CC = rho_CC_new[m][c];
            double rho_m_0  = rho_m_CC;

            // Compute rho_0(p) and Compute drho_0/dp = 1/(dp/drho_0)
            if (ice_matl[m]) {
              rho_m_0 = ice_matl[m]->getEOS()->computeRhoMicro(
                p_new, gamma[m][c], cv[m][c], Temp[m][c], rho_m_0);
              ice_matl[m]->getEOS()->computePressEOS(rho_m_0,
                                                     gamma[m][c],
                                                     cv[m][c],
                                                     Temp[m][c],
                                                     press_eos[m],
                                                     dp_drho[m],
                                                     dp_de[m]);
            } else if (mpm_matl[m]) {
              rho_m_0 = mpm_matl[m]->getConstitutiveModel()->computeRhoMicroCM(
                p_new, press_ref, mpm_matl[m], Temp[m][c], rho_m_0);
              mpm_matl[m]->getConstitutiveModel()->computePressEOSCM(
                rho_m_0,
                press_eos[m],
                press_ref,
                dp_drho[m],
                c_2,
                mpm_matl[m],
                Temp[m][c]);
            }

            double rho_CC_over_rho_0 = rho_m_CC / rho_m_0;

            double rho_CC_over_rho_0_sq = rho_CC_over_rho_0 / rho_m_0;
            double drho_dp              = 1.0 / dp_drho[m];

            // Update saved data
            rho_micro[m][c] = rho_m_0;
            vol_frac[m][c]  = rho_CC_over_rho_0;

            // Compute sum of rho_CC/rho_0
            A += rho_CC_over_rho_0;

            // Compute sum of rho_CC/rho_0^2*drho_0/dp
            B += rho_CC_over_rho_0_sq * drho_dp;
          }

          // Update pressure
          double p_upd = p_new + (A - 1.0) / B;
          p_new        = (p_upd > 0) ? p_upd : p_new;

          // Check that the new pressure is not significantly differenet from
          // the old
          delPress = p_new - press_new[c];
          if (std::abs(p_new - press_new[c]) / std::abs(p_new) < 1.0e-6) {

            // Update the volume fractions
            sum = 0;
            for (int m = 0; m < numALLMatls; m++) {
              if (ice_matl[m]) {
                rho_micro[m][c] = ice_matl[m]->getEOS()->computeRhoMicro(
                  p_new, gamma[m][c], cv[m][c], Temp[m][c], rho_micro[m][c]);
              } else if (mpm_matl[m]) {
                rho_micro[m][c] =
                  mpm_matl[m]->getConstitutiveModel()->computeRhoMicroCM(
                    p_new, press_ref, mpm_matl[m], Temp[m][c], rho_micro[m][c]);
              }
              vol_frac[m][c] = rho_CC_new[m][c] / rho_micro[m][c];
              sum += vol_frac[m][c];
            }

            //__________________________________
            // - Test for convergence
            //  If sum of vol_frac_CC ~= 1.0 then converged
            if (std::abs(sum - 1.0) < convergence_crit) {
              converged = true;
              //__________________________________
              // Find the speed of sound based on the converged solution
              for (int m = 0; m < numALLMatls; m++) {
                if (ice_matl[m]) {
                  ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c],
                                                         gamma[m][c],
                                                         cv[m][c],
                                                         Temp[m][c],
                                                         press_eos[m],
                                                         dp_drho[m],
                                                         dp_de[m]);

                  c_2 = dp_drho[m] +
                        dp_de[m] *
                          (press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]));

                } else if (mpm_matl[m]) {
                  mpm_matl[m]->getConstitutiveModel()->computePressEOSCM(
                    rho_micro[m][c],
                    press_eos[m],
                    press_ref,
                    dp_drho[m],
                    c_2,
                    mpm_matl[m],
                    Temp[m][c]);
                }

                speedSound[m][c] = sqrt(c_2); // Isentropic speed of sound

                //____ BB : B U L L E T   P R O O F I N G----
                // catch inf and nan speed sound values
                if (std::isnan(speedSound[m][c]) || c_2 == 0.0) {
                  std::ostringstream warn;
                  warn << "ERROR MPMICE::computeEquilPressure, MPM mat= " << m
                       << " cell= " << c << " p = " << p_new
                       << " sound speed is imaginary.\n";
                  warn << "speedSound = " << speedSound[m][c]
                       << " c_2 = " << c_2 << " press_eos = " << press_eos[m]
                       << " dp_drho = " << dp_drho[m] << " dp_de = " << dp_de[m]
                       << " rho_micro = " << rho_micro[m][c]
                       << " Temp = " << Temp[m][c] << std::endl;
                  throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
                }
              }
            }
          }

          // Update press_new
          press_new[c] = p_new;
        }

      } else { // Use the old equilibration algorithm

        while (count < d_ice->d_max_iter_equilibration && converged == false) {
          count++;
          //__________________________________
          // evaluate press_eos at cell i,j,k
          for (int m = 0; m < numALLMatls; m++) {
            if (ice_matl[m]) { // ICE
              ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c],
                                                     gamma[m][c],
                                                     cv[m][c],
                                                     Temp[m][c],
                                                     press_eos[m],
                                                     dp_drho[m],
                                                     dp_de[m]);
            } else if (mpm_matl[m]) { // MPM
              mpm_matl[m]->getConstitutiveModel()->computePressEOSCM(
                rho_micro[m][c],
                press_eos[m],
                press_ref,
                dp_drho[m],
                c_2,
                mpm_matl[m],
                Temp[m][c]);
            }
          }

          //__________________________________
          // - compute delPress
          // - update press_CC
          double A = 0., B = 0., C = 0.;

          for (int m = 0; m < numALLMatls; m++) {
            double Q     = press_new[c] - press_eos[m];
            double inv_y = (vol_frac[m][c] * vol_frac[m][c]) /
                           (dp_drho[m] * rho_CC_new[m][c] + d_SMALL_NUM);

            A += vol_frac[m][c];
            B += Q * inv_y;
            C += inv_y;
          }
          double vol_frac_not_close_packed = 1.;
          delPress = (A - vol_frac_not_close_packed - B) / C;

          press_new[c] += delPress;

          if (press_new[c] < convergence_crit) {
            press_new[c] = fabs(delPress);
          }

          //__________________________________
          // backout rho_micro_CC at this new pressure
          // - compute the updated volume fractions
          sum = 0;
          for (int m = 0; m < numALLMatls; m++) {
            if (ice_matl[m]) {
              rho_micro[m][c] =
                ice_matl[m]->getEOS()->computeRhoMicro(press_new[c],
                                                       gamma[m][c],
                                                       cv[m][c],
                                                       Temp[m][c],
                                                       rho_micro[m][c]);
            }
            if (mpm_matl[m]) {
              rho_micro[m][c] =
                mpm_matl[m]->getConstitutiveModel()->computeRhoMicroCM(
                  press_new[c],
                  press_ref,
                  mpm_matl[m],
                  Temp[m][c],
                  rho_micro[m][c]);
            }
            vol_frac[m][c] = rho_CC_new[m][c] / rho_micro[m][c];
            sum += vol_frac[m][c];
          }

          //__________________________________
          // - Test for convergence
          //  If sum of vol_frac_CC ~= 1.0 then converged
          if (std::abs(sum - vol_frac_not_close_packed) < convergence_crit) {
            converged = true;
            //__________________________________
            // Find the speed of sound based on the converged solution
            feclearexcept(FE_ALL_EXCEPT);
            for (int m = 0; m < numALLMatls; m++) {
              if (ice_matl[m]) {
                ice_matl[m]->getEOS()->computePressEOS(rho_micro[m][c],
                                                       gamma[m][c],
                                                       cv[m][c],
                                                       Temp[m][c],
                                                       press_eos[m],
                                                       dp_drho[m],
                                                       dp_de[m]);

                c_2 = dp_drho[m] +
                      dp_de[m] *
                        (press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]));

                speedSound[m][c] = sqrt(c_2); // Isentropic speed of sound

                //____ BB : B U L L E T   P R O O F I N G----
                // catch inf and nan speed sound values
                if (fetestexcept(FE_INVALID) != 0 || c_2 == 0.0) {
                  std::ostringstream warn;
                  warn << "ERROR MPMICE::computeEquilPressure, ICE mat= " << m
                       << " cell= " << c << " sound speed is imaginary.\n";
                  warn << "speedSound = " << speedSound[m][c]
                       << " c_2 = " << c_2 << " press_eos = " << press_eos[m]
                       << " dp_drho = " << dp_drho[m] << " dp_de = " << dp_de[m]
                       << " rho_micro = " << rho_micro[m][c]
                       << " Temp = " << Temp[m][c] << std::endl;
                  throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
                }

              } else if (mpm_matl[m]) {
                mpm_matl[m]->getConstitutiveModel()->computePressEOSCM(
                  rho_micro[m][c],
                  press_eos[m],
                  press_ref,
                  dp_drho[m],
                  c_2,
                  mpm_matl[m],
                  Temp[m][c]);

                speedSound[m][c] = sqrt(c_2); // Isentropic speed of sound

                //____ BB : B U L L E T   P R O O F I N G----
                // catch inf and nan speed sound values
                if (fetestexcept(FE_INVALID) != 0 || c_2 == 0.0) {
                  std::ostringstream warn;
                  warn << "ERROR MPMICE::computeEquilPressure, MPM mat= " << m
                       << " cell= " << c << " sound speed is imaginary.\n";
                  warn << "speedSound = " << speedSound[m][c]
                       << " c_2 = " << c_2 << " press_eos = " << press_eos[m]
                       << " dp_drho = " << dp_drho[m] << " dp_de = " << dp_de[m]
                       << " rho_micro = " << rho_micro[m][c]
                       << " Temp = " << Temp[m][c] << std::endl;
                  throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
                }
              }
            }
          }

          // Save iteration data for output in case of crash
          if (ds_EqPress.active()) {
            EqPress_dbg dbg;
            dbg.delPress   = delPress;
            dbg.press_new  = press_new[c];
            dbg.sumVolFrac = sum;
            dbg.count      = count;

            for (int m = 0; m < numALLMatls; m++) {
              EqPress_dbgMatl dmatl;
              dmatl.press_eos = press_eos[m];
              dmatl.volFrac   = vol_frac[m][c];
              dmatl.rhoMicro  = rho_micro[m][c];
              dmatl.rho_CC    = rho_CC_new[m][c];
              dmatl.temp_CC   = Temp[m][c];
              dmatl.mat       = m;
              dbg.matl.push_back(dmatl);
            }
            dbgEqPress.push_back(dbg);
          }

        } // end of converged
      }   // end if old algorithm for equilibration

      delPress_tmp[c] = delPress;

      //__________________________________
      // If the pressure solution has stalled out
      //  then try a binary search
      int binaryPressCount = 0;
      if (count >= d_ice->d_max_iter_equilibration) {
        // int lev = patch->getLevel()->getIndex();
        // cout << "WARNING:MPMICE:ComputeEquilibrationPressure "
        //      << " Cell : " << c << " on level " << lev << " having a
        //      difficult time converging. \n"
        //     << " Now performing a binary pressure search " << std::endl;

        binaryPressureSearch(Temp,
                             rho_micro,
                             vol_frac,
                             rho_CC_new,
                             speedSound,
                             dp_drho,
                             dp_de,
                             press_eos,
                             press,
                             press_new,
                             press_ref,
                             cv,
                             gamma,
                             convergence_crit,
                             numALLMatls,
                             count,
                             sum,
                             c);
        binaryPressCount = count;
      }

      nIterations[c] = count + binaryPressCount;
      test_max_iter  = std::max(test_max_iter, count);

      //__________________________________
      //      BULLET PROOFING
      // ignore BP if timestep restart has already been requested
      bool rts = new_dw->recomputeTimeStep();

      std::string message;
      bool allTestsPassed = true;
      if (test_max_iter == d_ice->d_max_iter_equilibration && !rts) {
        allTestsPassed = false;
        message += "Max. iterations reached ";
      }

      for (int m = 0; m < numALLMatls; m++) {
        ASSERT((vol_frac[m][c] > 0.0) || (vol_frac[m][c] < 1.0));
      }

      if (fabs(sum - 1.0) > convergence_crit && !rts) {
        allTestsPassed = false;
        message += " sum (volumeFractions) != 1 ";
      }

      if (press_new[c] < 0.0 && !rts) {
        allTestsPassed = false;
        message += " Computed pressure is < 0 ";
      }

      for (int m = 0; m < numALLMatls; m++) {
        if ((rho_micro[m][c] < 0.0 || vol_frac[m][c] < 0.0) && !rts) {
          allTestsPassed = false;
          message += " rho_micro < 0 || vol_frac < 0";
        }
      }
      if (allTestsPassed != true) { // throw an exception of there's a problem
        std::ostringstream warn;
        warn << "\nMPMICE::ComputeEquilibrationPressure: Cell " << c << ", L-"
             << L_indx << "\n"
             << message
             << "\nThis usually means that something much deeper has gone "
                "wrong with the simulation. "
             << "\nCompute equilibration pressure task is rarely the problem. "
             << "For more debugging information set the environmental "
                "variable:  \n"
             << "   SCI_DEBUG DBG_EqPress:+\n\n";

        warn << "INPUTS: \n";
        for (int m = 0; m < numALLMatls; m++) {
          warn << "\n matl: " << m << "\n"
               << "   rho_CC:     " << rho_CC_new[m][c] << "\n"
               << "   Temperature:   " << Temp[m][c] << "\n";
        }
        if (ds_EqPress.active()) {
          warn << "\nDetails on iterations " << std::endl;
          std::vector<EqPress_dbg>::iterator dbg_iter;
          for (dbg_iter = dbgEqPress.begin(); dbg_iter != dbgEqPress.end();
               dbg_iter++) {
            EqPress_dbg& d = *dbg_iter;
            warn << "Iteration:   " << d.count
                 << "  press_new:   " << d.press_new
                 << "  sumVolFrac:  " << d.sumVolFrac
                 << "  delPress:    " << d.delPress << "\n";
            for (int m = 0; m < numALLMatls; m++) {
              warn << "  matl: " << d.matl[m].mat
                   << "  press_eos:  " << d.matl[m].press_eos
                   << "  volFrac:    " << d.matl[m].volFrac
                   << "  rhoMicro:   " << d.matl[m].rhoMicro
                   << "  rho_CC:     " << d.matl[m].rho_CC
                   << "  Temp:       " << d.matl[m].temp_CC << "\n";
            }
          }
        }
      } // all testsPassed
    }   // end of cell interator
    if (cout_norm.active()) {
      cout_norm << "max number of iterations in any cell \t" << test_max_iter
                << endl;
    }

    //__________________________________
    // Now change how rho_CC is defined to
    // rho_CC = mass/cell_volume  NOT mass/mat_volume
    for (int m = 0; m < numALLMatls; m++) {
      if (ice_matl[m]) {
        rho_CC_new[m].copyData(rho_CC_old[m]);
      }
    }

    //__________________________________
    // - update Boundary conditions
    //   Don't set Lodi bcs, we already compute Press
    //   in all the extra cells.
    // - make copy of press for implicit calc.
    std::unique_ptr<CustomBCDriver::customBC_localVars> BC_localVars =
      std::make_unique<CustomBCDriver::customBC_localVars>();
    preprocess_CustomBCs("EqPress",
                         old_dw,
                         new_dw,
                         d_ice->d_ice_labels.get(),
                         patch,
                         999,
                         d_ice->d_BC_globalVars.get(),
                         BC_localVars.get());

    setBC(press_new,
          rho_micro,
          placeHolder,
          d_ice->d_surroundingMatl_indx,
          "rho_micro",
          "Pressure",
          patch,
          d_materialManager,
          0,
          new_dw,
          d_ice->d_BC_globalVars.get(),
          BC_localVars.get(),
          isNotInitialTimeStep);

    delete_CustomBCs(d_ice->d_BC_globalVars.get(), BC_localVars.get());

    //__________________________________
    // compute sp_vol_CC
    for (int m = 0; m < numALLMatls; m++) {
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        const IntVector& c = *iter;
        sp_vol_new[m][c]   = 1.0 / rho_micro[m][c];
      }

      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      setSpecificVolBC(sp_vol_new[m],
                       "SpecificVol",
                       false,
                       rho_CC_new[m],
                       vol_frac[m],
                       patch,
                       d_materialManager,
                       indx);
    }
    //__________________________________
    //  compute f_theta
    for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
      const IntVector& c = *iter;
      sumKappa[c]        = 0.0;
      for (int m = 0; m < numALLMatls; m++) {
        kappa[m][c] = sp_vol_new[m][c] / (speedSound[m][c] * speedSound[m][c]);
        sumKappa[c] += vol_frac[m][c] * kappa[m][c];
      }
      for (int m = 0; m < numALLMatls; m++) {
        f_theta[m][c] = vol_frac[m][c] * kappa[m][c] / sumKappa[c];
      }
    }

    //____ BB : B U L L E T   P R O O F I N G----
    for (int m = 0; m < numALLMatls; m++) {
      Material* matl = d_materialManager->getMaterial(m);
      int indx       = matl->getDWIndex();
      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        if (std::isnan(kappa[m][*iter]) || std::isinf(kappa[m][*iter])) {
          std::ostringstream warn;
          warn << "MPMICE:(L-" << m << "):computeEquilibrationPressure, mat "
               << indx << " cell " << *iter << " kappa = " << kappa[m][*iter]
               << " vol_frac = " << vol_frac[m][*iter]
               << " sp_vol_new = " << sp_vol_new[m][*iter]
               << " speedSound = " << speedSound[m][*iter] << std::endl;
          throw InvalidValue(warn.str(), __FILE__, __LINE__);
        }
      }
    }
  } // patches
}

/* ---------------------------------------------------------------------
   Function~  MPMICE::binaryPressureSearch--
   Purpose:   When the technique for find the equilibration pressure
   craps out then try this method.
   Reference:  See Numerical Methods by Robert W. Hornbeck.
   _____________________________________________________________________*/
void
MPMICE::binaryPressureSearch(std::vector<constCCVariable<double>>& Temp,
                             std::vector<CCVariable<double>>& rho_micro,
                             std::vector<CCVariable<double>>& vol_frac,
                             std::vector<CCVariable<double>>& rho_CC_new,
                             std::vector<CCVariable<double>>& speedSound,
                             std::vector<double>& dp_drho,
                             std::vector<double>& dp_de,
                             std::vector<double>& press_eos,
                             constCCVariable<double>& press,
                             CCVariable<double>& press_new,
                             double press_ref,
                             std::vector<constCCVariable<double>>& cv,
                             std::vector<constCCVariable<double>>& gamma,
                             double convergence_crit,
                             int numALLMatls,
                             int& count,
                             double& sum,
                             IntVector c)
{
  // Start over for this cell using a binary search
  //  std::cout << " cell " << c << " Starting binary pressure search "<<
  //  std::endl;
  count          = 0;
  bool converged = false;
  double c_2;
  double Pleft = 0., Pright = 0., Ptemp = 0., Pm = 0.;
  double rhoMicroR = 0., rhoMicroL = 0.;
  std::vector<double> vfR(numALLMatls);
  std::vector<double> vfL(numALLMatls);
  Pm              = press[c];
  double residual = 1.0;

  while (count < d_ice->d_max_iter_equilibration && converged == false) {
    count++;
    sum = 0.;
    for (int m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
      if (ice_matl) { // ICE
        rho_micro[m][c] = ice_matl->getEOS()->computeRhoMicro(
          Pm, gamma[m][c], cv[m][c], Temp[m][c], rho_micro[m][c]);
      }
      if (mpm_matl) { // MPM
        rho_micro[m][c] = mpm_matl->getConstitutiveModel()->computeRhoMicroCM(
          Pm, press_ref, mpm_matl, Temp[m][c], rho_micro[m][c]);
      }
      vol_frac[m][c] = rho_CC_new[m][c] / rho_micro[m][c];
      sum += vol_frac[m][c];
    } // loop over matls

    residual = 1. - sum;

    if (fabs(residual) <= convergence_crit) {
      converged    = true;
      press_new[c] = Pm;
      //__________________________________
      // Find the speed of sound at ijk
      // needed by eos and the the explicit
      // del pressure function
      // feclearexcept(FE_ALL_EXCEPT);
      for (int m = 0; m < numALLMatls; m++) {
        Material* matl        = d_materialManager->getMaterial(m);
        ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
        MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
        if (ice_matl) { // ICE
          ice_matl->getEOS()->computePressEOS(rho_micro[m][c],
                                              gamma[m][c],
                                              cv[m][c],
                                              Temp[m][c],
                                              press_eos[m],
                                              dp_drho[m],
                                              dp_de[m]);
          c_2 = dp_drho[m] +
                dp_de[m] * (press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]));
        }
        if (mpm_matl) { // MPM
          mpm_matl->getConstitutiveModel()->computePressEOSCM(rho_micro[m][c],
                                                              press_eos[m],
                                                              press_ref,
                                                              dp_drho[m],
                                                              c_2,
                                                              mpm_matl,
                                                              Temp[m][c]);
        }
        speedSound[m][c] = sqrt(c_2); // Isentropic speed of sound

        //____ BB : B U L L E T   P R O O F I N G----
        // catch inf and nan speed sound values
        // if (fetestexcept(FE_INVALID) != 0 || c_2 == 0.0) {
        //   std::ostringstream warn;
        //  warn<<"ERROR MPMICE::binaryPressSearch, mat= "<< m << " cell= "
        //      << c << " sound speed is imaginary.\n";
        //  warn << "speedSound = " << speedSound[m][c] << " c_2 = " << c_2
        //       << " press_eos = " << press_eos[m]
        //       << " dp_drho = " << dp_drho[m]
        //       << " dp_de = " << dp_de[m]
        //       << " rho_micro = " << rho_micro[m][c] << " Temp = " <<
        //       Temp[m][c] << std::endl;
        //  throw ProblemSetupException(warn.str(), __FILE__, __LINE__ );
        //}
      }
    }
    // Initial guess
    if (count == 1) {
      Pleft  = DBL_EPSILON;
      Pright = 1.0 / DBL_EPSILON;
      Ptemp  = .5 * (Pleft + Pright);
    }

    double sumR = 0., sumL = 0.;
    for (int m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);

      if (ice_matl) { //  ICE
        rhoMicroR = ice_matl->getEOS()->computeRhoMicro(
          Pright, gamma[m][c], cv[m][c], Temp[m][c], rho_micro[m][c]);
        rhoMicroL = ice_matl->getEOS()->computeRhoMicro(
          Pleft, gamma[m][c], cv[m][c], Temp[m][c], rho_micro[m][c]);
      }
      if (mpm_matl) { //  MPM
        rhoMicroR = mpm_matl->getConstitutiveModel()->computeRhoMicroCM(
          Pright, press_ref, mpm_matl, Temp[m][c], rho_micro[m][c]);
        rhoMicroL = mpm_matl->getConstitutiveModel()->computeRhoMicroCM(
          Pleft, press_ref, mpm_matl, Temp[m][c], rho_micro[m][c]);
      }
      vfR[m] = rho_CC_new[m][c] / rhoMicroR;
      vfL[m] = rho_CC_new[m][c] / rhoMicroL;
      sumR += vfR[m];
      sumL += vfL[m];

      //      std::cout << "matl: " << m << " vol_frac_L: " << vfL[m] << "
      //      vol_frac_R: " << vfR[m]
      //           << " rho_CC: " << rho_CC_new[m][c] << " rho_micro_L: " <<
      //           rhoMicroL << " rhoMicroR: " << rhoMicroR << std::endl;
    } // all matls

    //    std::cout << "Pm = " << Pm << "\t P_L: " << Pleft << "\t P_R: " <<
    //    Pright << "\t 1.-sum " << residual << " \t sumR: " << sumR << " \t
    //    sumL " << sumL << std::endl;

    //__________________________________
    //  come up with a new guess
    double prod = (1. - sumR) * (1. - sumL);
    if (prod < 0.) {
      Ptemp = Pleft;
      Pleft = .5 * (Pleft + Pright);
    } else {
      Pright = Pleft;
      Pleft  = Ptemp;
      Pleft  = 0.5 * (Pleft + Pright);
    }
    Pm = .5 * (Pleft + Pright);
    //   std::cout << setprecision(15);

  } // end of converged

#ifdef D_STRICT
  if (count >= d_ice->d_max_iter_equilibration) {
    std::ostringstream desc;
    desc << "**ERROR** Binary pressure search failed to converge in cell" << c
         << std::endl;
    throw ConvergenceFailure(desc.str(),
                             d_ice->d_max_iter_equilibration,
                             fabs(residual),
                             convergence_crit,
                             __FILE__,
                             __LINE__);
  }
#else
  if (count >= d_ice->d_max_iter_equilibration) {
    std::cout
      << "**WARNING** Binary pressure search failed to converge in cell " << c
      << " after " << d_ice->d_max_iter_equilibration
      << " iterations.  Final residual is " << fabs(residual)
      << " and convergence tolerance is " << convergence_crit << std::endl;
    std::cout << "  Continuing with unconverged value of pressure = " << Pm
              << std::endl;
    press_new[c] = Pm;
    //__________________________________
    // Find the speed of sound at ijk
    // needed by eos and the the explicit
    // del pressure function
    for (int m = 0; m < numALLMatls; m++) {
      Material* matl        = d_materialManager->getMaterial(m);
      ICEMaterial* ice_matl = dynamic_cast<ICEMaterial*>(matl);
      MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial*>(matl);
      if (ice_matl) { // ICE
        ice_matl->getEOS()->computePressEOS(rho_micro[m][c],
                                            gamma[m][c],
                                            cv[m][c],
                                            Temp[m][c],
                                            press_eos[m],
                                            dp_drho[m],
                                            dp_de[m]);
        c_2 = dp_drho[m] +
              dp_de[m] * (press_eos[m] / (rho_micro[m][c] * rho_micro[m][c]));
      }
      if (mpm_matl) { // MPM
        mpm_matl->getConstitutiveModel()->computePressEOSCM(rho_micro[m][c],
                                                            press_eos[m],
                                                            press_ref,
                                                            dp_drho[m],
                                                            c_2,
                                                            mpm_matl,
                                                            Temp[m][c]);
      }
      speedSound[m][c] = sqrt(c_2); // Isentropic speed of sound
      if (std::isnan(vol_frac[m][c])) {
        std::ostringstream warn;
        warn << "    Material " << m << " vol. frac = " << vol_frac[m][c]
             << " rho = " << rho_micro[m][c] << " press = " << press_eos[m]
             << " dp_drho = " << dp_drho[m] << " c^2 = " << c_2 << std::endl;
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }
    }
  }
#endif
}

void
MPMICE::scheduleInterpolatePressCCToPressNC(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSubset* press_matl,
                                            const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();
  if (!d_mpm->d_mpm_flags->doMPMOnLevel(L_indx,
                                        level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(
    patches, cout_doing, "MPMICE::scheduleInterpolatePressCCToPressNC");

  Task* t = scinew Task("MPMICE::interpolatePressCCToPressNC",
                        this,
                        &MPMICE::interpolatePressCCToPressNC);

  Ghost::GhostType gac = Ghost::AroundCells;
  t->requires(Task::NewDW, d_ice_labels->press_CCLabel, press_matl, gac, 1);
  t->computes(d_mpmice_labels->press_NCLabel, press_matl);

  sched->addTask(t, patches, matls);
}

void
MPMICE::interpolatePressCCToPressNC(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse*,
                                    DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing interpolatePressCCToPressNC");

    constCCVariable<double> pressCC;
    NCVariable<double> pressNC;

    Ghost::GhostType gac = Ghost::AroundCells;
    new_dw->get(pressCC, d_ice_labels->press_CCLabel, 0, patch, gac, 1);
    new_dw->allocateAndPut(pressNC, d_mpmice_labels->press_NCLabel, 0, patch);
    pressNC.initialize(0.0);

    IntVector cIdx[8];
    // Interpolate CC pressure to nodes
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      patch->findCellsFromNode(*iter, cIdx);
      for (int in = 0; in < 8; in++) {
        pressNC[*iter] += .125 * pressCC[cIdx[in]];
      }
    }

    // Apply grid boundary conditions to the pressure before storing the data
    string inter_type = d_mpm->d_mpm_flags->d_interpolatorType;
    MPMBoundCond bc;
    bc.setBoundaryCondition(patch, 0, "Pressure", pressNC, inter_type);
  }
}

void
MPMICE::scheduleInterpolatePAndGradP(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSubset* press_matl,
                                     const MaterialSubset* /*one_matl*/,
                                     const MaterialSubset* mpm_matl,
                                     const MaterialSet* all_matls)
{
  if (!d_mpm->d_mpm_flags->doMPMOnLevel(
        getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPMICE::scheduleInterpolatePAndGradP");

  Task* t = scinew Task(
    "MPMICE::interpolatePAndGradP", this, &MPMICE::interpolatePAndGradP);
  Ghost::GhostType gac = Ghost::AroundCells;

  t->requires(Task::NewDW,
              d_mpmice_labels->press_NCLabel,
              press_matl,
              gac,
              d_num_ghost_nodes);
  t->requires(Task::NewDW, d_mpmice_labels->cMassLabel, mpm_matl, gac, 1);
  t->requires(Task::OldDW, d_mpm_labels->pXLabel, mpm_matl, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pSizeLabel, mpm_matl, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pDefGradLabel, mpm_matl, Ghost::None);

  t->computes(d_mpm_labels->pPressureLabel, mpm_matl);
  sched->addTask(t, patches, all_matls);
}

void
MPMICE::interpolatePAndGradP(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(
      patches, patch, cout_doing, "Doing interpolatePressureToParticles");

    auto interpolator = d_mpm->d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    double p_ref = d_ice->getRefPress();
    constNCVariable<double> pressNC;
    Ghost::GhostType gac = Ghost::AroundCells;
    new_dw->get(pressNC,
                d_mpmice_labels->press_NCLabel,
                0,
                patch,
                gac,
                d_num_ghost_nodes);

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(indx, patch);
      ParticleVariable<double> pPressure;
      constParticleVariable<Point> px;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
      new_dw->allocateAndPut(pPressure, d_mpm_labels->pPressureLabel, pset);

      //__________________________________
      // Interpolate NC pressure to particles
      for (auto idx : *pset) {
        double press = 0.;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(
          px[idx], ni, S, pSize[idx], pDefGrad[idx]);

        for (size_t k = 0; k < ni.size(); k++) {
          press += pressNC[ni[k]] * S[k];
        }
        pPressure[idx] = press - p_ref;
      }
    } // numMPMMatls
  }   // patches
}

void
MPMICE::scheduleComputeLagrangianValuesMPM(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSubset* one_matl,
                                           const MaterialSet* mpm_matls)
{
  if (d_mpm->d_mpm_flags->doMPMOnLevel(
        getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {

    printSchedule(
      patches, cout_doing, "MPMICE::scheduleComputeLagrangianValuesMPM");

    Task* t = scinew Task("MPMICE::computeLagrangianValuesMPM",
                          this,
                          &MPMICE::computeLagrangianValuesMPM);

    const MaterialSubset* mss = mpm_matls->getUnion();
    Ghost::GhostType gac      = Ghost::AroundCells;
    Ghost::GhostType gn       = Ghost::None;
    t->requires(Task::NewDW, d_mpm_labels->gVelocityStarLabel, mss, gac, 1);
    t->requires(Task::NewDW, d_mpm_labels->gMassLabel, gac, 1);
    t->requires(Task::NewDW, d_mpm_labels->gTemperatureStarLabel, gac, 1);
    t->requires(Task::OldDW, d_mpm_labels->NC_CCweightLabel, one_matl, gac, 1);
    t->requires(Task::NewDW, d_mpmice_labels->cMassLabel, gn);
    t->requires(Task::NewDW, d_ice_labels->int_eng_source_CCLabel, gn);
    t->requires(Task::NewDW, d_ice_labels->mom_source_CCLabel, gn);

    t->requires(Task::NewDW, d_mpmice_labels->temperature_CCLabel, gn);
    t->requires(Task::NewDW, d_mpmice_labels->velocity_CCLabel, gn);

    if (d_ice->d_models.size() > 0 && !do_mlmpmice) {
      t->requires(Task::NewDW, d_ice_labels->modelMass_srcLabel, gn);
      t->requires(Task::NewDW, d_ice_labels->modelMom_srcLabel, gn);
      t->requires(Task::NewDW, d_ice_labels->modelEng_srcLabel, gn);
    }

    t->computes(d_ice_labels->mass_L_CCLabel);
    t->computes(d_ice_labels->mom_L_CCLabel);
    t->computes(d_ice_labels->int_eng_L_CCLabel);

    sched->addTask(t, patches, mpm_matls);
  }
}

void
MPMICE::computeLagrangianValuesMPM(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, VarLabel::find(timeStep_name));
  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeLagrangianValuesMPM");

    int numMatls           = d_materialManager->getNumMaterials("MPM");
    Vector dx              = patch->dCell();
    double cellVol         = dx.x() * dx.y() * dx.z();
    double very_small_mass = d_TINY_RHO * cellVol;
    Ghost::GhostType gn    = Ghost::None;
    Ghost::GhostType gac   = Ghost::AroundCells;

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();

      // Create arrays for the grid data
      constNCVariable<double> gMass, gVolume, gtempstar;
      constNCVariable<Vector> gVelocity;
      CCVariable<Vector> cmomentum;
      CCVariable<double> int_eng_L, mass_L;
      constCCVariable<double> cmass, Temp_CC_sur, int_eng_src;
      constCCVariable<Vector> vel_CC_sur, mom_source;
      new_dw->get(gMass, d_mpm_labels->gMassLabel, indx, patch, gac, 1);
      new_dw->get(
        gVelocity, d_mpm_labels->gVelocityStarLabel, indx, patch, gac, 1);
      new_dw->get(
        gtempstar, d_mpm_labels->gTemperatureStarLabel, indx, patch, gac, 1);
      new_dw->get(cmass, d_mpmice_labels->cMassLabel, indx, patch, gn, 0);
      new_dw->get(
        Temp_CC_sur, d_mpmice_labels->temperature_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        vel_CC_sur, d_mpmice_labels->velocity_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        mom_source, d_ice_labels->mom_source_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        int_eng_src, d_ice_labels->int_eng_source_CCLabel, indx, patch, gn, 0);

      new_dw->allocateAndPut(mass_L, d_ice_labels->mass_L_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        cmomentum, d_ice_labels->mom_L_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        int_eng_L, d_ice_labels->int_eng_L_CCLabel, indx, patch);

      cmomentum.initialize(Vector(0.0, 0.0, 0.0));
      int_eng_L.initialize(0.);
      mass_L.initialize(0.);
      double cv = mpm_matl->getSpecificHeat();

      IntVector nodeIdx[8];
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        patch->findNodesFromCell(c, nodeIdx);
        double int_eng_L_mpm = 0.0;
        double int_eng_L_sur = cmass[c] * Temp_CC_sur[c] * cv;
        Vector cmomentum_mpm = Vector(0.0, 0.0, 0.0);
        Vector cmomentum_sur = vel_CC_sur[c] * cmass[c];

        for (int in = 0; in < 8; in++) {
          double NC_CCw_mass = NC_CCweight[nodeIdx[in]] * gMass[nodeIdx[in]];
          cmomentum_mpm += gVelocity[nodeIdx[in]] * NC_CCw_mass;
          int_eng_L_mpm += gtempstar[nodeIdx[in]] * cv * NC_CCw_mass;
        }

        // Biswajit: 10/3/2013: Hack to make the bubble simulation go through
        // TODO: Check the internal energy rate calculation
        int_eng_L_mpm = (int_eng_L_mpm > 0.0) ? int_eng_L_mpm : 0.0;
        int_eng_L_mpm += int_eng_src[c];
        if (!d_rigidMPM) {
          cmomentum_mpm += mom_source[c];
        }

        //__________________________________
        // set cmomentum/int_eng_L to either
        // what's calculated from mpm or
        // the surrounding value.
        // If you change this you must also change MPMICE::interpolateNCToCC_0
        double one_or_zero = (cmass[c] - very_small_mass) / cmass[c];
        cmomentum[c] =
          (1.0 - one_or_zero) * cmomentum_sur + one_or_zero * cmomentum_mpm;
        int_eng_L[c] =
          (1.0 - one_or_zero) * int_eng_L_sur + one_or_zero * int_eng_L_mpm;
      }
      //__________________________________
      //  NO REACTION
      if (d_ice->d_models.size() == 0 || do_mlmpmice) {
        for (CellIterator iter = patch->getExtraCellIterator(); !iter.done();
             iter++) {
          IntVector c = *iter;
          mass_L[c]   = cmass[c];
        }
      }
      //__________________________________
      //   M O D E L   B A S E D   E X C H A N G E
      // The reaction can't completely eliminate all the mass, momentum and
      // internal E. If it does then we'll get erroneous vel, and temps in
      // CCMomExchange.  If the mass goes to min_mass then cmomentum and
      // int_eng_L need to be scaled by min_mass to avoid inf temp and vel_CC
      if (d_ice->d_models.size() > 0 && !do_mlmpmice) {
        constCCVariable<double> modelMass_src;
        constCCVariable<double> modelEng_src;
        constCCVariable<Vector> modelMom_src;
        new_dw->get(
          modelMass_src, d_ice_labels->modelMass_srcLabel, indx, patch, gn, 0);
        new_dw->get(
          modelMom_src, d_ice_labels->modelMom_srcLabel, indx, patch, gn, 0);
        new_dw->get(
          modelEng_src, d_ice_labels->modelEng_srcLabel, indx, patch, gn, 0);

        for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          //  must have a minimum mass
          double min_mass  = very_small_mass;
          double inv_cmass = 1.0 / cmass[c];
          mass_L[c]        = std::max((cmass[c] + modelMass_src[c]), min_mass);

          //  must have a minimum momentum
          for (int dir = 0; dir < 3; dir++) { // loop over all three directons
            double min_mom_L = min_mass * cmomentum[c][dir] * inv_cmass;
            double mom_L_tmp = cmomentum[c][dir] + modelMom_src[c][dir];

            // Preserve the original sign on momemtum
            // Use d_SMALL_NUMs to avoid nans when mom_L_temp = 0.0
            double plus_minus_one =
              (mom_L_tmp + d_SMALL_NUM) / (fabs(mom_L_tmp + d_SMALL_NUM));

            mom_L_tmp = (mom_L_tmp / mass_L[c]) * (cmass[c] + modelMass_src[c]);

            cmomentum[c][dir] =
              plus_minus_one * std::max(fabs(mom_L_tmp), fabs(min_mom_L));
          }
          // must have a minimum int_eng
          double min_int_eng = min_mass * int_eng_L[c] * inv_cmass;
          int_eng_L[c] =
            (int_eng_L[c] / mass_L[c]) * (cmass[c] + modelMass_src[c]);
          int_eng_L[c] =
            std::max((int_eng_L[c] + modelEng_src[c]), min_int_eng);
        }
      } // if(model.size() >0)

      //__________________________________
      //  Set Boundary conditions
      setBC(cmomentum,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(int_eng_L,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);

      //---- B U L L E T   P R O O F I N G------
      // ignore BP if timestep restart has already been requested
      IntVector neg_cell;
      std::ostringstream warn;
      bool rts = new_dw->recomputeTimeStep();

      if (d_testForNegTemps_mpm) {
        if (!areAllValuesPositive(int_eng_L, neg_cell) && !rts) {
          int L = getLevel(patches)->getIndex();
          warn << "ERROR MPMICE:(" << L << "):computeLagrangianValuesMPM, mat "
               << indx << " cell " << neg_cell << " int_eng_L_CC "
               << int_eng_L[neg_cell] << "\n ";
          throw InvalidValue(warn.str(), __FILE__, __LINE__);
        }
      }
    } // numMatls
  }   // patches
}

void
MPMICE::scheduleCoarsenLagrangianValuesMPM(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* mpm_matls)
{
  printSchedule(
    patches, cout_doing, "MPMICE:scheduleCoarsenLagrangianValues mpm_matls");

  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->rho_CCLabel,
                            1e-12,
                            true,
                            "std"); // modifies
  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->mass_L_CCLabel,
                            1.9531e-15,
                            false,
                            "sum");
  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->mom_L_CCLabel,
                            Vector(0, 0, 0),
                            false,
                            "sum");
  scheduleCoarsenVariableCC(sched,
                            patches,
                            mpm_matls,
                            d_ice_labels->int_eng_L_CCLabel,
                            0.0,
                            false,
                            "sum");
}

void
MPMICE::scheduleComputeCCVelAndTempRates(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* mpm_matls)
{
  printSchedule(
    patches, cout_doing, "MPMICE::scheduleComputeCCVelAndTempRates");

  Task* t = scinew Task("MPMICE::computeCCVelAndTempRates",
                        this,
                        &MPMICE::computeCCVelAndTempRates);

  Ghost::GhostType gn = Ghost::None;

  t->requires(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->requires(Task::NewDW, d_ice_labels->mass_L_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->mom_L_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->int_eng_L_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->mom_L_ME_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->eng_L_ME_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->int_eng_source_CCLabel, gn);
  t->requires(Task::NewDW, d_ice_labels->mom_source_CCLabel, gn);
  t->requires(Task::OldDW, d_mpm_labels->heatRate_CCLabel, gn);

  t->computes(d_ice_labels->dTdt_CCLabel);
  t->computes(d_ice_labels->dVdt_CCLabel);
  t->computes(d_mpm_labels->heatRate_CCLabel);

  sched->addTask(t, patches, mpm_matls);
}

void
MPMICE::computeCCVelAndTempRates(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, VarLabel::find(timeStep_name));
  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeCCVelAndTempRates");

    //__________________________________
    // This is where I interpolate the CC
    // changes to NCs for the MPMMatls
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();
      CCVariable<double> dTdt_CC, heatRate;
      CCVariable<Vector> dVdt_CC;

      constCCVariable<double> mass_L_CC, old_heatRate;
      constCCVariable<Vector> mom_L_ME_CC, old_mom_L_CC, mom_source;
      constCCVariable<double> eng_L_ME_CC, old_int_eng_L_CC, int_eng_src;

      double cv = mpm_matl->getSpecificHeat();

      Ghost::GhostType gn = Ghost::None;
      new_dw->get(
        old_mom_L_CC, d_ice_labels->mom_L_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        old_int_eng_L_CC, d_ice_labels->int_eng_L_CCLabel, indx, patch, gn, 0);
      new_dw->get(mass_L_CC, d_ice_labels->mass_L_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        mom_L_ME_CC, d_ice_labels->mom_L_ME_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        eng_L_ME_CC, d_ice_labels->eng_L_ME_CCLabel, indx, patch, gn, 0);
      old_dw->get(
        old_heatRate, d_mpm_labels->heatRate_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        mom_source, d_ice_labels->mom_source_CCLabel, indx, patch, gn, 0);
      new_dw->get(
        int_eng_src, d_ice_labels->int_eng_source_CCLabel, indx, patch, gn, 0);

      new_dw->allocateAndPut(dTdt_CC, d_ice_labels->dTdt_CCLabel, indx, patch);
      new_dw->allocateAndPut(dVdt_CC, d_ice_labels->dVdt_CCLabel, indx, patch);
      new_dw->allocateAndPut(
        heatRate, d_mpm_labels->heatRate_CCLabel, indx, patch);

      dTdt_CC.initialize(0.0);
      dVdt_CC.initialize(Vector(0.0));
      //__________________________________
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        if (!d_rigidMPM) {
          dVdt_CC[c] = (mom_L_ME_CC[c] - (old_mom_L_CC[c] - mom_source[c])) /
                       (mass_L_CC[c] * delT);
        }
        dTdt_CC[c] = (eng_L_ME_CC[c] - (old_int_eng_L_CC[c] - int_eng_src[c])) /
                     (mass_L_CC[c] * cv * delT);
        double heatRte = (eng_L_ME_CC[c] - old_int_eng_L_CC[c]) / delT;
        heatRate[c]    = .05 * heatRte + .95 * old_heatRate[c];
      }
      setBC(dTdt_CC,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
      setBC(dVdt_CC,
            "set_if_sym_BC",
            patch,
            d_materialManager,
            indx,
            new_dw,
            isNotInitialTimeStep);
    }
  } // patches
}

void
MPMICE::scheduleInterpolateCCToNC(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* mpm_matls)
{
  if (!d_mpm->d_mpm_flags->doMPMOnLevel(
        getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPMICE::scheduleInterpolateCCToNC");

  Task* t =
    scinew Task("MPMICE::interpolateCCToNC", this, &MPMICE::interpolateCCToNC);
  const MaterialSubset* mss = mpm_matls->getUnion();
  Ghost::GhostType gan      = Ghost::AroundNodes;
  Ghost::GhostType gac      = Ghost::AroundCells;

  t->requires(Task::OldDW, d_ice_labels->delTLabel, getLevel(patches));
  t->requires(Task::NewDW, d_ice_labels->dVdt_CCLabel, gan, 1);
  t->requires(Task::NewDW, d_ice_labels->dTdt_CCLabel, gan, 1);

  if (d_ice->d_models.size() > 0) {
    t->requires(Task::NewDW, d_mpmice_labels->cMassLabel, gac, 1);
    t->requires(Task::NewDW, d_ice_labels->modelMass_srcLabel, gac, 1);
  }

  t->modifies(d_mpm_labels->gVelocityStarLabel, mss);
  t->modifies(d_mpm_labels->gAccelerationLabel, mss);
  t->computes(d_mpm_labels->massBurnFractionLabel, mss);
  t->computes(d_mpm_labels->dTdt_NCLabel);

  sched->addTask(t, patches, mpm_matls);
}

void
MPMICE::interpolateCCToNC(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing interpolateCCToNC");

    //__________________________________
    // This is where I interpolate the CC
    // changes to NCs for the MPMMatls
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    delt_vartype delT;
    old_dw->get(delT, d_ice_labels->delTLabel);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();
      NCVariable<Vector> gacceleration, gVelocity;
      NCVariable<double> dTdt_NC, massBurnFraction;

      constCCVariable<double> dTdt_CC;
      constCCVariable<Vector> dVdt_CC;

      new_dw->getModifiable(
        gVelocity, d_mpm_labels->gVelocityStarLabel, indx, patch);
      new_dw->getModifiable(
        gacceleration, d_mpm_labels->gAccelerationLabel, indx, patch);

      Ghost::GhostType gan = Ghost::AroundNodes;
      new_dw->get(dTdt_CC, d_ice_labels->dTdt_CCLabel, indx, patch, gan, 1);
      new_dw->get(dVdt_CC, d_ice_labels->dVdt_CCLabel, indx, patch, gan, 1);

      new_dw->allocateAndPut(
        massBurnFraction, d_mpm_labels->massBurnFractionLabel, indx, patch);
      new_dw->allocateAndPut(dTdt_NC, d_mpm_labels->dTdt_NCLabel, indx, patch);

      dTdt_NC.initialize(0.0);
      massBurnFraction.initialize(0.);
      IntVector cIdx[8];
      //__________________________________
      //
      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        patch->findCellsFromNode(*iter, cIdx);
        for (int in = 0; in < 8; in++) {
          gVelocity[*iter] += dVdt_CC[cIdx[in]] * delT * .125;
          gacceleration[*iter] += dVdt_CC[cIdx[in]] * .125;
          dTdt_NC[*iter] += dTdt_CC[cIdx[in]] * .125;
        }
      }
      //__________________________________
      //  inter-material phase transformation
      if (d_ice->d_models.size() > 0) {
        constCCVariable<double> modelMass_src, mass_CC;
        Ghost::GhostType gac = Ghost::AroundCells;
        new_dw->get(
          modelMass_src, d_ice_labels->modelMass_srcLabel, indx, patch, gac, 1);
        new_dw->get(mass_CC, d_mpmice_labels->cMassLabel, indx, patch, gac, 1);

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          patch->findCellsFromNode(*iter, cIdx);
          for (int in = 0; in < 8; in++) {
            massBurnFraction[*iter] +=
              (std::abs(modelMass_src[cIdx[in]]) / mass_CC[cIdx[in]]) * .125;
          }
        }
      } // if(models >0 )
    }   // mpmMatls
  }     // patches
}

/* _____________________________________________________________________
   MPMICE::scheduleFinalizeTimestep--
   This task called at the very bottom of the timestep,
   after scheduleTimeAdvance and the scheduleCoarsen.

   This is scheduled on every level.
   _____________________________________________________________________*/
void
MPMICE::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{
  cout_doing << "----------------------------" << endl;
  cout_doing << d_myworld->myRank()
             << " MPMICE::scheduleFinalizeTimestep\t\t\t\tL-"
             << level->getIndex() << std::endl;

  const PatchSet* ice_patches         = level->eachPatch();
  const MaterialSet* ice_matls        = d_materialManager->allMaterials("ICE");
  const MaterialSet* all_matls        = d_materialManager->allMaterials();
  const MaterialSet* mpm_matls        = d_materialManager->allMaterials("MPM");
  const MaterialSubset* ice_matls_sub = ice_matls->getUnion();

  d_ice->scheduleConservedtoPrimitive_Vars(
    sched, ice_patches, ice_matls_sub, ice_matls, "finalizeTimestep");

  d_ice->scheduleTestConservation(sched, ice_patches, ice_matls_sub, all_matls);

  for (auto& am : d_analysisModules) {
    am->scheduleDoAnalysis_preReloc(sched, level);
  }

  // only do on finest level until we get AMR MPM
  if (level->getIndex() == level->getGrid()->numLevels() - 1) {
    MPMICE::scheduleParticleRelocation(sched, level, mpm_matls);
  }

  for (auto& am : d_analysisModules) {
    am->scheduleDoAnalysis(sched, level);
  }
  cout_doing << "---------------------------------------------------------"
             << std::endl;
}

/*!
 *  Purpose:  If needed update the labels and materialSubsets that are passed
 *            to Relocate and schedule particle relocate
 */
void
MPMICE::scheduleParticleRelocation(SchedulerP& sched,
                                   const LevelP& level,
                                   const MaterialSet* mpm_matls)
{
  //__________________________________
  //  Unmodified labels and matls subset
  std::vector<std::vector<const VarLabel*>> old_labels;
  std::vector<std::vector<const VarLabel*>> new_labels;

  old_labels = d_mpm->d_particleState_preReloc;
  new_labels = d_mpm->d_particleState;

  const MaterialSubset* old_mss = mpm_matls->getSubset(0);
  MaterialSubset* new_mss       = scinew MaterialSubset();
  new_mss->addReference();
  new_mss->addSubset(old_mss);

  //__________________________________
  // Update VarLabels and matl Subsets
  for (auto& model : d_ice->d_models) {
    ParticleModel* pModel = dynamic_cast<ParticleModel*>(model.get());
    if (pModel) {
      old_labels.push_back(pModel->d_oldLabels);
      new_labels.push_back(pModel->d_newLabels);
      new_mss->addSubset(pModel->d_matl_mss);
    }
  }

  //__________________________________
  //  create a new material set containing the
  //  the updated matlSubset.
  MaterialSet* newMatlSet = scinew MaterialSet();
  newMatlSet->addSubset(new_mss);
  newMatlSet->addReference();

  sched->scheduleParticleRelocation(level,
                                    d_mpm_labels->pXLabel_preReloc,
                                    old_labels,
                                    d_mpm_labels->pXLabel,
                                    new_labels,
                                    d_mpm_labels->pParticleIDLabel,
                                    newMatlSet);

  if (newMatlSet && newMatlSet->removeReference()) {
    delete newMatlSet;
  }
}

void
MPMICE::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  if (d_switchCriteria) {
    d_switchCriteria->scheduleSwitchTest(level, sched);
  }
}

void
MPMICE::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                     SchedulerP& sched)
{
  d_ice->scheduleInitialErrorEstimate(coarseLevel, sched);
  d_mpm->scheduleInitialErrorEstimate(coarseLevel, sched);
}

void
MPMICE::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  d_ice->scheduleErrorEstimate(coarseLevel, sched);
  d_mpm->scheduleErrorEstimate(coarseLevel, sched);
}

void
MPMICE::scheduleRefineInterface(const LevelP& fineLevel,
                                SchedulerP& scheduler,
                                bool needOld,
                                bool needNew)
{
  d_ice->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);
  d_mpm->scheduleRefineInterface(fineLevel, scheduler, needOld, needNew);

  if (fineLevel->getIndex() > 0 && scheduler->isCopyDataTimestep() &&
      d_mpm->d_mpm_flags->doMPMOnLevel(fineLevel->getIndex(),
                                       fineLevel->getGrid()->numLevels())) {
    cout_doing << d_myworld->myRank()
               << " MPMICE::scheduleRefineInterface \t\t\tL-"
               << fineLevel->getIndex() << std::endl;

    Task* task = scinew Task("MPMICE::refineCoarseFineInterface",
                             this,
                             &MPMICE::refineCoarseFineInterface);

    const MaterialSet* all_matls   = d_materialManager->allMaterials();
    const MaterialSubset* one_matl = d_ice->d_press_matl;

    task->modifies(d_mpm_labels->NC_CCweightLabel, one_matl);

    scheduler->addTask(task, fineLevel->eachPatch(), all_matls);
  }
}

void
MPMICE::refineCoarseFineInterface(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset*,
                                  DataWarehouse* fine_old_dw,
                                  DataWarehouse* fine_new_dw)
{
  // This isn't actually refining anything, it is simply reinitializing
  // NC_CCweight after regridding on all levels finer than 0 because
  // copyData doesn't copy extra cell data.
  const Level* level = getLevel(patches);
  if (level->getIndex() > 0) {
    cout_doing << d_myworld->myRank() << " Doing refineCoarseFineInterface"
               << "\t\t\t MPMICE L-" << level->getIndex()
               << " Patches: " << *patches << std::endl;

    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      //__________________________________
      // NC_CCweight
      NCVariable<double> NC_CCweight;
      fine_new_dw->getModifiable(
        NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch);
      //__________________________________
      // - Initialize NC_CCweight = 0.125
      // - Find the walls with symmetry BC and double NC_CCweight
      NC_CCweight.initialize(0.125);

      std::vector<Patch::FaceType> bf;
      patch->getBoundaryFaces(bf);
      for (auto& face : bf) {
        int mat_id = 0;
        if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {

          for (CellIterator iter =
                 patch->getFaceIterator(face, Patch::FaceNodes);
               !iter.done();
               iter++) {
            NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
          } // cell iterator
        }   // if symmetry
      }     // for patch faces
    }       // for patches
  }         // if level
}

void
MPMICE::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  d_ice->scheduleRefine(patches, sched);
  d_mpm->scheduleRefine(patches, sched);

  printSchedule(patches, cout_doing, "MPMICE::scheduleRefine");

  Task* task = scinew Task("MPMICE::refine", this, &MPMICE::refine);

  task->computes(d_mpm_labels->heatRate_CCLabel);
  task->computes(d_ice_labels->specificVolume_CCLabel);
  task->computes(d_mpmice_labels->velocity_CCLabel);
  task->computes(d_ice_labels->temperature_CCLabel);

  sched->addTask(task, patches, d_materialManager->allMaterials("MPM"));
}

void
MPMICE::refine(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* /*matls*/,
               DataWarehouse* old_dw,
               DataWarehouse* new_dw)
{
  timeStep_vartype timeStep;
  old_dw->get(timeStep, VarLabel::find(timeStep_name));
  bool isNotInitialTimeStep = (timeStep > 0);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing refine");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      cout_doing << d_myworld->myRank() << " Doing refine on patch "
                 << patch->getID() << " material # = " << dwi << std::endl;

      // for now, create 0 heat flux
      CCVariable<double> heatFlux;
      new_dw->allocateAndPut(
        heatFlux, d_mpm_labels->heatRate_CCLabel, dwi, patch);
      heatFlux.initialize(0.0);

      CCVariable<double> rho_micro, sp_vol_CC, rho_CC, Temp_CC, vol_frac_CC;
      CCVariable<Vector> vel_CC;

      new_dw->allocateTemporary(rho_micro, patch);
      new_dw->allocateTemporary(rho_CC, patch);
      new_dw->allocateTemporary(vol_frac_CC, patch);

      new_dw->allocateAndPut(
        sp_vol_CC, d_ice_labels->specificVolume_CCLabel, dwi, patch);
      new_dw->allocateAndPut(
        Temp_CC, d_mpmice_labels->temperature_CCLabel, dwi, patch);
      new_dw->allocateAndPut(
        vel_CC, d_mpmice_labels->velocity_CCLabel, dwi, patch);

      mpm_matl->initializeDummyCCVariables(
        rho_micro, rho_CC, Temp_CC, vel_CC, vol_frac_CC, patch);
      //__________________________________
      //  Set boundary conditions
      setBC(rho_micro,
            "Density",
            patch,
            d_materialManager,
            dwi,
            new_dw,
            isNotInitialTimeStep);
      setBC(Temp_CC,
            "Temperature",
            patch,
            d_materialManager,
            dwi,
            new_dw,
            isNotInitialTimeStep);
      setBC(vel_CC,
            "Velocity",
            patch,
            d_materialManager,
            dwi,
            new_dw,
            isNotInitialTimeStep);
      for (auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        sp_vol_CC[*iter] = 1.0 / rho_micro[*iter];
      }

      //__________________________________
      //    B U L L E T   P R O O F I N G
      IntVector neg_cell;
      std::ostringstream warn;
      if (!areAllValuesPositive(rho_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << dwi << " cell "
             << neg_cell << " position: " << pt << " rho_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(Temp_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << dwi << " cell "
             << neg_cell << " position: " << pt << " Temp_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      if (!areAllValuesPositive(sp_vol_CC, neg_cell)) {
        Point pt = patch->getCellPosition(neg_cell);
        warn << "ERROR MPMICE::actuallyInitialize, mat " << dwi << " cell "
             << neg_cell << " poistion: " << pt << " sp_vol_CC is negative\n";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
    } // mpmMatls
  }   // patches
}

void
MPMICE::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
  d_ice->scheduleCoarsen(coarseLevel, sched);
  d_mpm->scheduleCoarsen(coarseLevel, sched);
}

/* For old RigidMPM code*/
void
MPMICE::scheduleRefinePressCC(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSubset* press_matl,
                              const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "MPMICE::scheduleRefinePressCC");

  MaterialSet* press_matls = scinew MaterialSet();
  press_matls->add(0);
  press_matls->addReference();

  scheduleRefineVariableCC(
    sched, patches, press_matls, d_ice_labels->press_CCLabel);
  if (press_matls->removeReference()) {
    delete press_matls;
  }
}

void
MPMICE::scheduleRefineCC(SchedulerP& sched,
                         const PatchSet* patches,
                         const MaterialSet* mpm_matls)
{
  if (!d_mpm->d_mpm_flags->doMPMOnLevel(
        getLevel(patches)->getIndex(),
        getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPMICE::scheduleRefineCC");
  scheduleRefineVariableCC(
    sched, patches, mpm_matls, d_ice_labels->dTdt_CCLabel);
  scheduleRefineVariableCC(
    sched, patches, mpm_matls, d_ice_labels->dVdt_CCLabel);
}

void
MPMICE::scheduleRefineVariableCC(SchedulerP& sched,
                                 const PatchSet* patches,
                                 const MaterialSet* matls,
                                 const VarLabel* variable)
{
  std::ostringstream taskName;
  taskName << "MPMICE::refineVariable(" << variable->getName() << ")";
  Task* t;

  switch (variable->typeDescription()->getSubType()->getType()) {
    case TypeDescription::Type::double_type:
      {
      auto func = &MPMICE::refineVariableCC<double>;
      t         = scinew Task(taskName.str().c_str(), this, func, variable);
      }
      break;
    case TypeDescription::Type::Vector:
      {
      auto func = &MPMICE::refineVariableCC<Vector>;
      t         = scinew Task(taskName.str().c_str(), this, func, variable);
      }
      break;
    default:
      throw InternalError(
        "Unknown variable type for refine", __FILE__, __LINE__);
  }

  Ghost::GhostType gac = Ghost::AroundCells;
  t->requires(
    Task::NewDW, variable, 0, Task::CoarseLevel, 0, Task::NormalDomain, gac, 1);
  t->computes(variable);
  sched->addTask(t, patches, matls);
}

template<typename T>
void
MPMICE::refineVariableCC(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* matls,
                         DataWarehouse*,
                         DataWarehouse* new_dw,
                         const VarLabel* variable)
{
  const Level* fineLevel   = getLevel(patches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
  IntVector refineRatio(fineLevel->getRefinementRatio());

  for (int p = 0; p < patches->size(); p++) {
    const Patch* finePatch = patches->get(p);
    std::ostringstream message;
    message << "Doing refineVariableCC (" << variable->getName() << ")\t\t\t";
    printTask(patches, finePatch, cout_doing, message.str());

    // region of fine space that will correspond to the coarse we need to get
    IntVector cl, ch, fl, fh;
    IntVector bl(0, 0, 0); // boundary layer cells
    int nGhostCells           = 1;
    bool returnExclusiveRange = true;

    getCoarseLevelRange(finePatch,
                        coarseLevel,
                        cl,
                        ch,
                        fl,
                        fh,
                        bl,
                        nGhostCells,
                        returnExclusiveRange);

    for (int m = 0; m < matls->size(); m++) {
      int indx = matls->get(m);

      CCVariable<T> fine_q_CC;
      new_dw->allocateAndPut(fine_q_CC, variable, indx, finePatch);

      constCCVariable<T> coarse_q_CC;

      new_dw->getRegion(
        coarse_q_CC, variable, indx, coarseLevel, cl, ch, false);

      // Only interpolate over the intersection of the fine and coarse patches
      // coarse cell
      //      linearInterpolation<T>(coarse_q_CC, coarseLevel, fineLevel,
      //                             refineRatio, lo, hi, fine_q_CC);

      piecewiseConstantInterpolation<T>(
        coarse_q_CC, fineLevel, fl, fh, fine_q_CC);
    }
  }
}
