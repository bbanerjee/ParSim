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

#include <CCA/Components/MPM/ImpMPM.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/MPM/Core/ImpMPMFlags.h>
#include <CCA/Components/MPM/Core/ImpMPMLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/MPMUtils.h>

#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConduction.h>
#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConductionTasks.h>

#include <CCA/Components/MPM/ImpMPMSolvers/PetscSolver.h>
#include <CCA/Components/MPM/ImpMPMSolvers/SimpleSolver.h>

#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>

#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/PerPatchVars.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Math/FastMatrix.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>

using namespace Uintah;

using NCdoubleArray            = std::vector<NCVariable<double>>;
using NCVectorArray            = std::vector<NCVariable<Vector>>;
using NCMatrix3Array           = std::vector<NCVariable<Matrix3>>;
using CellParticleTempMap      = std::multimap<IntVector, ParticleTempShape>;
using CellParticleTempPair     = std::pair<IntVector, ParticleTempShape>;
using CellParticleTempMapArray = std::vector<CellParticleTempMap>;

//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "IMPM:+,IMPM_debug:+".....
//  bash     : export SCI_DEBUG="IMPM:+,IMPM_debug:+"
//  default is OFF
static DebugStream cout_doing("IMPM", false);
static DebugStream cout_dbg("IMPM_Debug", false);

ImpMPM::ImpMPM(const ProcessorGroup* myworld,
               const MaterialManagerP& mat_manager)
  : SimulationCommon(myworld, mat_manager)
  , MPMCommon(d_materialManager)
{
  d_mpm_labels       = std::make_unique<MPMLabel>();
  d_impmpmLabels    = std::make_unique<ImpMPMLabel>();
  d_mpm_flags        = std::make_unique<ImpMPMFlags>(myworld);

  d_oneMaterial = scinew MaterialSubset();
  d_oneMaterial->add(0);
  d_oneMaterial->addReference();
}

ImpMPM::~ImpMPM()
{
  if (d_perprocPatches && d_perprocPatches->removeReference()) {
    delete d_perprocPatches;
    std::cout << "Freeing patches!!\n";
  }

  if (d_oneMaterial && d_oneMaterial->removeReference()) {
    delete d_oneMaterial;
  }

  if (d_loadCurveIndex && d_loadCurveIndex->removeReference()) {
    delete d_loadCurveIndex;
  }

  MPMPhysicalBCFactory::clean();
}

void
ImpMPM::problemSetup(const ProblemSpecP& prob_spec,
                     const ProblemSpecP& restart_prob_spec,
                     GridP& grid)
{
  cout_doing << " Doing ImpMPM::problemSetup " << std::endl;

  d_scheduler->setPositionVar(d_mpm_labels->pXLabel);

  ProblemSpecP mpm_ps         = 0;
  ProblemSpecP restart_mat_ps = 0;

  bool isRestart = false;

  ProblemSpecP prob_spec_mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");
  if (prob_spec_mat_ps) {
    restart_mat_ps = prob_spec;
  } else if (restart_prob_spec) {
    isRestart      = true;
    restart_mat_ps = restart_prob_spec;
  } else {
    restart_mat_ps = prob_spec;
  }

  ProblemSpecP mpm_soln_ps = restart_mat_ps->findBlock("MPM");

  std::string integrator_type;
  if (mpm_soln_ps) {

    // Read all MPM flags (look in MPMFlags.cc)
    d_mpm_flags->readMPMFlags(restart_mat_ps, d_output);

    if (d_mpm_flags->d_integratorType != "implicit") {
      throw ProblemSetupException("Can't use explicit integration with -impm",
                                  __FILE__,
                                  __LINE__);
    }

    // convert text representation of face into FaceType
    for (auto face : d_mpm_flags->d_boundaryTractionFaceStrings) {
      auto faceType = Patch::invalidFace;
      for (auto ft = Patch::startFace; ft <= Patch::endFace;
           ft      = Patch::nextFace(ft)) {
        if (Patch::getFaceName(ft) == face) {
          faceType = ft;
          break;
        }
      }
      if (faceType != Patch::invalidFace) {
        d_boundaryTractionFaces.push_back(faceType);
      } else {
        std::cerr << "warning: ignoring unknown face '" << face << "'\n";
      }
    }
  }

  // read in AMR flags from the main ups file
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  if (amr_ps) {
    ProblemSpecP mpm_amr_ps = amr_ps->findBlock("MPM");
    mpm_amr_ps->getWithDefault("min_grid_level", d_mpm_flags->d_minGridLevel, 0);
    mpm_amr_ps->getWithDefault("max_grid_level",
                               d_mpm_flags->d_maxGridLevel,
                               1000);
  }

  if (d_mpm_flags->d_8or27 == 8) {
    d_numGhostParticles = 1;
    d_numGhostNodes     = 1;
  } else if (d_mpm_flags->d_8or27 == 27 || d_mpm_flags->d_8or27 == 64) {
    d_numGhostParticles = 2;
    d_numGhostNodes     = 2;
  }

  // Search for the MaterialProperties block and then get the MPM section
  ProblemSpecP mat_ps =
    restart_mat_ps->findBlockWithOutAttribute("MaterialProperties");
  ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");
  ProblemSpecP contact_ps = mpm_mat_ps->findBlock("contact");

  if (contact_ps) {
    contact_ps->getWithDefault("type", d_contactType, "null");
  }

  if (d_contactType == "rigid") {
    d_rigidContact = true;
    Vector defaultDir(1, 1, 1);
    contact_ps->getWithDefault("direction", d_contactDirections, defaultDir);
    contact_ps->getWithDefault("stop_time",
                               d_contactStopTime,
                               std::numeric_limits<double>::max());
    contact_ps->getWithDefault("velocity_after_stop",
                               d_velocityAfterContactStop,
                               Vector(0, 0, 0));
  }

  setParticleGhostLayer(Ghost::AroundNodes, 1);

  MPMPhysicalBCFactory::create(restart_mat_ps, grid, d_mpm_flags.get());
  if (MPMPhysicalBCFactory::mpmPhysicalBCs.size() == 0 &&
      d_mpm_flags->d_useLoadCurves) {
    throw ProblemSetupException("No load curve in ups, d_useLoadCurve==true?",
                                __FILE__,
                                __LINE__);
  }

  materialProblemSetup(restart_mat_ps, d_mpm_flags.get(), isRestart);

  if (d_mpm_flags->d_solverType == "petsc") {
    d_solver = std::make_unique<MPMPetscSolver>();
  } else {
    d_solver = std::make_unique<SimpleSolver>();
  }
  d_solver->initialize();

  d_defGradComputer =
    std::make_unique<DeformationGradientComputer>(d_materialManager,
                                                  d_mpm_labels.get(),
                                                  d_mpm_flags.get());

  // setup sub scheduler
  d_subsched = d_scheduler->createSubScheduler();
  d_subsched->initialize(3, 1);
  d_subsched->clearMappings();
  d_subsched->mapDataWarehouse(Task::ParentOldDW, 0);
  d_subsched->mapDataWarehouse(Task::ParentNewDW, 1);
  d_subsched->mapDataWarehouse(Task::OldDW, 2);
  d_subsched->mapDataWarehouse(Task::NewDW, 3);

  d_recompileSubsched = true;

  d_heatConductionTasks =
    std::make_unique<ImplicitHeatConductionTasks>(restart_mat_ps,
                                                  d_materialManager,
                                                  d_mpm_labels.get(),
                                                  d_mpm_flags.get());

  d_switchCriteria =
    dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));
  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(restart_mat_ps,
                                   restart_prob_spec,
                                   d_materialManager);
  }

  // Pull out from Time section
  ProblemSpecP time_ps = restart_mat_ps->findBlock("Time");
  if (time_ps) {
    time_ps->get("delt_init", d_initialDt);
  }
}

void
ImpMPM::outputProblemSpec(ProblemSpecP& root_ps)
{
  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP flags_ps = root->appendChild("MPM");
  d_mpm_flags->outputProblemSpec(flags_ps);

  ProblemSpecP mat_ps = root->findBlockWithOutAttribute("MaterialProperties");
  if (mat_ps == nullptr) {
    mat_ps = root->appendChild("MaterialProperties");
  }

  ProblemSpecP mpm_ps = mat_ps->appendChild("MPM");
  for (size_t i = 0; i < d_materialManager->getNumMaterials("MPM"); i++) {
    MPMMaterial* mat =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", i));
    ProblemSpecP cm_ps = mat->outputProblemSpec(mpm_ps);
  }

  ProblemSpecP contact_ps = mpm_ps->appendChild("contact");
  contact_ps->appendElement("type", d_contactType);
  contact_ps->appendElement("direction", d_contactDirections);
  contact_ps->appendElement("stop_time", d_contactStopTime);
  contact_ps->appendElement("velocity_after_stop", d_velocityAfterContactStop);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps   = physical_bc_ps->appendChild("MPM");
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    particleBC->outputProblemSpec(mpm_ph_bc_ps);
  }
}

void
ImpMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("ImpMPM::actuallyInitialize",
                        this,
                        &ImpMPM::actuallyInitialize);

  const PatchSet* patches = level->eachPatch();

  t->computes(d_mpm_labels->partCountLabel);
  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_mpm_labels->pDispLabel);
  t->computes(d_mpm_labels->pMassLabel);
  t->computes(d_mpm_labels->pVolumeLabel);
  t->computes(d_mpm_labels->pFiberDirLabel);
  t->computes(d_mpm_labels->pVelocityLabel);
  t->computes(d_mpm_labels->pAccelerationLabel);
  t->computes(d_mpm_labels->pExternalForceLabel);
  t->computes(d_mpm_labels->pTemperatureLabel);
  t->computes(d_mpm_labels->pTempPreviousLabel);
  t->computes(d_mpm_labels->pSizeLabel);
  t->computes(d_mpm_labels->pParticleIDLabel);
  t->computes(d_mpm_labels->pStressLabel);
  t->computes(d_mpm_labels->pRefinedLabel);
  t->computes(d_mpm_labels->pCellNAPIDLabel);

  t->computes(d_mpm_labels->pCoriolisImportanceLabel);
  t->computes(d_mpm_labels->pBodyForceAccLabel);

  t->computes(d_mpm_labels->pExternalHeatFluxLabel);
  t->computes(d_mpm_labels->heatRate_CCLabel);

  t->computes(d_mpm_labels->delTLabel, level.get_rep());

  size_t numMPM = d_materialManager->getNumMaterials("MPM");
  for (size_t m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  if (d_mpm_flags->d_artificialViscosity) {
    t->computes(d_mpm_labels->p_qLabel); //  only used for imp -> exp transition
  }

  if (d_mpm_flags->d_useLoadCurves) {
    t->computes(d_mpm_labels->pLoadCurveIDLabel);
  }

  // For friction contact
  t->computes(d_mpm_labels->pSurfLabel);

  if (!d_mpm_flags->d_doGridReset) {
    t->computes(d_mpm_labels->gDisplacementLabel);
  }

  t->computes(d_mpm_labels->NC_CCweightLabel, d_oneMaterial);
  if (!d_mpm_flags->d_tempSolve) {
    t->computes(d_mpm_labels->gTemperatureLabel, d_oneMaterial);
  }

  d_perprocPatches = d_loadBalancer->getPerProcessorPatchSet(level);
  sched->addTask(t, d_perprocPatches, d_materialManager->allMaterials("MPM"));

  t = scinew Task("ImpMPM::printParticleCount",
                  this,
                  &ImpMPM::printParticleCount);
  t->requires(Task::NewDW, d_mpm_labels->partCountLabel);
  t->setType(Task::OncePerProc);
  sched->addTask(t, d_perprocPatches, d_materialManager->allMaterials("MPM"));

  if (d_mpm_flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    d_heatConductionTasks->scheduleInitializeHeatFluxBCs(level, sched);

    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
  }
}

void
ImpMPM::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
}

void
ImpMPM::scheduleSwitchInitialization(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  if (d_mpm_flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    if (d_myworld->myRank() == 0) {
      std::cout
        << " \n--------------------------------------------------------------"
        << std::endl;
      std::cout << " ImpMPM: the heat flux BC cannot be applied on the timestep"
                << std::endl;
      std::cout
        << " immediately after a component switch.  The computes/requires "
        << std::endl;
      std::cout << " cannot be met and one pseudo timestep must take place"
                << endl;
      std::cout
        << " ---------------------------------------------------------------\n"
        << std::endl;
    }

    d_heatConductionTasks->scheduleInitializeHeatFluxBCs(level, sched);

    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
  }
}

void
ImpMPM::actuallyInitialize(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  particleIndex totalParticles = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing IMPM::actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpm_labels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    CCVariable<double> heatFlux;
    new_dw->allocateAndPut(heatFlux, d_mpm_labels->heatRate_CCLabel, 0, patch);
    heatFlux.initialize(1.0);

    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", matl));
      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);
      mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                         mpm_matl,
                                                         new_dw);
      if (!d_mpm_flags->d_doGridReset) {
        int matID = mpm_matl->getDWIndex();
        NCVariable<Vector> gDisplacement;
        new_dw->allocateAndPut(gDisplacement,
                               d_mpm_labels->gDisplacementLabel,
                               matID,
                               patch);
        gDisplacement.initialize(Vector(0.));
      }
    }

    std::string interp_type = d_mpm_flags->d_interpolatorType;
    if ((interp_type == "gimp" || interp_type == "3rdorderBS" ||
         interp_type == "cpdi" || interp_type == "cpti" ||
         interp_type == "cpgimp")) {
      proc0cout << "__________________________________\n"
                << "WARNING: Use of GIMP/3rdorderBS/cpdi/cpgimp with Implicit "
                   "MPM is untested and may not work at this time.\n\n";
    }

    IntVector num_extra_cells = patch->getExtraCells();
    IntVector periodic        = patch->getLevel()->getPeriodicBoundaries();

    if (interp_type == "linear" && num_extra_cells != IntVector(0, 0, 0)) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>linear</interpolator> \n"
          << " you should also use <extraCells>[0,0,0]</extraCells> \n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    } else if ((interp_type == "gimp" || interp_type == "3rdorderBS" ||
                interp_type == "cpdi" || interp_type == "cpti" ||
                interp_type == "cpgimp") &&
               (num_extra_cells + periodic) != IntVector(1, 1, 1)) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>gimp</interpolator> \n"
          << " or <interpolator>3rdorderBS</interpolator> \n"
          << " or <interpolator>cpdi</interpolator> \n"
          << " or <interpolator>cpti</interpolator> \n"
          << " or <interpolator>cpgimp</interpolator> \n"
          << " you must also use extraCells and/or periodicBCs such that\n"
          << " the sum of the two is [1,1,1].\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }

    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight,
                           d_mpm_labels->NC_CCweightLabel,
                           0,
                           patch);
    NC_CCweight.initialize(0.125);
    for (auto face = Patch::startFace; face <= Patch::endFace;
         face      = Patch::nextFace(face)) {
      int mat_id = 0;
      if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {
        for (auto iter = patch->getFaceIterator(face, Patch::FaceNodes);
             !iter.done();
             iter++) {
          NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
        }
      }
    }
    if (!d_mpm_flags->d_tempSolve) {
      NCVariable<double> gTemperature;
      new_dw->allocateAndPut(gTemperature,
                             d_mpm_labels->gTemperatureLabel,
                             0,
                             patch);
      gTemperature.initialize(0.);
    }
  }
  new_dw->put(sumlong_vartype(totalParticles), d_mpm_labels->partCountLabel);
}

void
ImpMPM::printParticleCount(const ProcessorGroup* pg,
                           const PatchSubset*,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  if (pg->myRank() == 0) {
    sumlong_vartype pcount;
    new_dw->get(pcount, d_mpm_labels->partCountLabel);
    std::cerr << "Created " << (long)pcount << " total particles\n";
  }
}

void
ImpMPM::scheduleInitializePressureBCs(const LevelP& level, SchedulerP& sched)
{
  d_loadCurveIndex = scinew MaterialSubset();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int nofPressureBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "Pressure") {
      d_loadCurveIndex->add(nofPressureBCs++);
    }
  }
  if (nofPressureBCs > 0) {

    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task("ImpMPM::countMaterialPointsPerLoadCurve",
                          this,
                          &ImpMPM::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t,
                   level->eachPatch(),
                   d_materialManager->allMaterials("MPM"));

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("ImpMPM::initializePressureBC",
                    this,
                    &ImpMPM::initializePressureBC);
    t->requires(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpm_labels->pExternalForceLabel);
    sched->addTask(t,
                   level->eachPatch(),
                   d_materialManager->allMaterials("MPM"));
  }

  d_loadCurveIndex->removeReference();
}

void
ImpMPM::countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse*,
                                        DataWarehouse* new_dw)
{
  // Find the number of pressure BCs in the problem
  int nofPressureBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "Pressure") {
      nofPressureBCs++;

      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        int numPts         = 0;
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == (nofPressureBCs)) {
              ++numPts;
            }
          }
        } // matl loop
        new_dw->put(sumlong_vartype(numPts),
                    d_mpm_labels->materialPointsPerLoadCurveLabel,
                    0,
                    nofPressureBCs - 1);
      } // patch loop
    }
  }
}

void
ImpMPM::initializePressureBC(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse*,
                             DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;

  if (cout_dbg.active()) {
    cout_dbg << "Current Time (Initialize Pressure BC) = " << time << std::endl;
  }

  // Calculate the force vector at each particle
  int nofPressureBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "Pressure") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart,
                  d_mpm_labels->materialPointsPerLoadCurveLabel,
                  0,
                  nofPressureBCs++);

      // Save the material points per load curve in the PressureBC object
      PressureBC* pbc = dynamic_cast<PressureBC*>(particleBC.get());
      pbc->numMaterialPoints(numPart);

      if (cout_dbg.active()) {
        cout_dbg << "    Load Curve = " << nofPressureBCs
                 << " Num Particles = " << numPart << std::endl;
      }

      // Calculate the force per particle at t = 0.0
      double forcePerPart = pbc->forcePerParticle(time);

      // Loop through the patches and calculate the force vector
      // at each particle
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<Point> pX;
          constParticleVariable<int> pLoadCurveID;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX, d_mpm_labels->pXLabel, pset);
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);
          new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

          constParticleVariable<Vector> pDisp;
          new_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(pExternalForce,
                                d_mpm_labels->pExternalForceLabel,
                                pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == nofPressureBCs) {
              pExternalForce[idx] = pbc->getForceVector(pX[idx],
                                                        pDisp[idx],
                                                        forcePerPart,
                                                        time,
                                                        pDefGrad[idx]);
            }
          }

        } // matl loop
      }   // patch loop
    }
  }
}

void
ImpMPM::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  printSchedule(level, cout_doing, "IMPM::scheduleComputeStableTimestep");

  Task* t = scinew Task("ImpMPM::actuallyComputeStableTimestep",
                        this,
                        &ImpMPM::actuallyComputeStableTimestep);

  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  if (d_mpm_flags->doMPMOnLevel(level->getIndex(),
                               level->getGrid()->numLevels())) {
    t->requires(Task::OldDW, d_mpm_labels->delTLabel);
    t->requires(Task::NewDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  }

  // compute a delT on all levels even levels where impm is not running
  t->computes(d_mpm_labels->delTLabel, level.get_rep());

  sched->addTask(t, level->eachPatch(), matls);
}

void
ImpMPM::actuallyComputeStableTimestep(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  // compute a delT on all levels even levels where impm is not running
  const Level* level = getLevel(patches);
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    new_dw->put(delt_vartype(999), d_mpm_labels->delTLabel, level);
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::actuallyComputeStableTimestep");

    if (d_numIterations == 0) {
      new_dw->put(delt_vartype(d_initialDt),
                  d_mpm_labels->delTLabel,
                  patch->getLevel());
    } else {
      Vector dx = patch->dCell();
      delt_vartype old_delT;
      old_dw->get(old_delT, d_mpm_labels->delTLabel, patch->getLevel());

      int numMPMMatls = d_materialManager->getNumMaterials("MPM");
      for (int m = 0; m < numMPMMatls; m++) {
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
        int dwindex = mpm_matl->getDWIndex();

        ParticleSubset* pset = new_dw->getParticleSubset(dwindex, patch);

        constParticleVariable<Vector> pVelocity;
        new_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);

        Vector ParticleSpeed(1.e-12, 1.e-12, 1.e-12);

        for (auto idx : *pset) {
          ParticleSpeed =
            Vector(Max(std::abs(pVelocity[idx].x()), ParticleSpeed.x()),
                   Max(std::abs(pVelocity[idx].y()), ParticleSpeed.y()),
                   Max(std::abs(pVelocity[idx].z()), ParticleSpeed.z()));
        }
        ParticleSpeed   = dx / ParticleSpeed;
        double delT_new = .8 * ParticleSpeed.minComponent();

        delT_new = old_delT;

        double old_dt = old_delT;
        if (d_numIterations <= d_mpm_flags->d_numItersToIncreaseDelT) {
          old_dt = d_mpm_flags->d_delTIncreaseFactor * old_delT;
        }
        if (d_numIterations >= d_mpm_flags->d_numItersToDecreaseDelT) {
          old_dt = d_mpm_flags->d_delTDecreaseFactor * old_delT;
        }
        delT_new = std::min(delT_new, old_dt);

        new_dw->put(delt_vartype(delT_new),
                    d_mpm_labels->delTLabel,
                    patch->getLevel());
      }
    }
  }
}

void
ImpMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  d_perprocPatches = d_loadBalancer->getPerProcessorPatchSet(level);
  d_perprocPatches->addReference();

  scheduleComputeParticleBodyForce(sched, d_perprocPatches, matls);

  scheduleApplyExternalLoads(sched, d_perprocPatches, matls);

  scheduleInterpolateParticlesToGrid(sched,
                                     d_perprocPatches,
                                     d_oneMaterial,
                                     matls);

  d_heatConductionTasks->scheduleProjectHeatSource(sched,
                                                   d_perprocPatches,
                                                   d_oneMaterial,
                                                   matls);

  scheduleFindSurfaceParticles(sched, d_perprocPatches, matls);

  scheduleDestroyMatrix(sched, d_perprocPatches, matls, false);
  scheduleCreateMatrix(sched, d_perprocPatches, matls);

  d_heatConductionTasks->scheduleDestroyAndCreateMatrix(sched,
                                                        d_perprocPatches,
                                                        matls);

  scheduleApplyBoundaryConditions(sched, d_perprocPatches, matls);

  d_heatConductionTasks->scheduleApplyBoundaryConditions(sched,
                                                         d_perprocPatches,
                                                         matls);

  scheduleComputeContact(sched, d_perprocPatches, matls);
  scheduleFindFixedDOF(sched, d_perprocPatches, matls);

  d_heatConductionTasks->scheduleSolve(sched, d_perprocPatches, matls);

  scheduleIterate(sched, level, d_perprocPatches, matls);

  if (!d_mpm_flags->d_doGridReset) {
    scheduleUpdateTotalDisplacement(sched, d_perprocPatches, matls);
  }

  scheduleComputeDeformationGradient(sched, d_perprocPatches, matls);
  scheduleComputeStressTensor(sched, d_perprocPatches, matls);
  scheduleComputeAcceleration(sched, d_perprocPatches, matls);
  scheduleInterpolateToParticlesAndUpdate(sched, d_perprocPatches, matls);
  scheduleInterpolateStressToGrid(sched, d_perprocPatches, matls);

  sched->scheduleParticleRelocation(level,
                                    d_mpm_labels->pXLabel_preReloc,
                                    d_particleState_preReloc,
                                    d_mpm_labels->pXLabel,
                                    d_particleState,
                                    d_mpm_labels->pParticleIDLabel,
                                    matls);
}

void
ImpMPM::scheduleComputeParticleBodyForce(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeParticleBodyForce");
  Task* t = scinew Task("IMPM::computeParticleBodyForce",
                        this,
                        &ImpMPM::computeParticleBodyForce);

  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  t->computes(d_mpm_labels->pBodyForceAccLabel_preReloc);
  t->computes(d_mpm_labels->pCoriolisImportanceLabel_preReloc);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeParticleBodyForce(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  // Get the MPM flags and make local copies
  Uintah::Point rotation_center = d_mpm_flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis  = d_mpm_flags->d_coordRotationAxis;
  double rotation_speed         = d_mpm_flags->d_coordRotationSpeed;
  // Uintah::Point body_ref_point = d_mpm_flags->d_coord_rotation_body_ref_point;

  // Compute angular velocity vector (omega)
  Uintah::Vector omega = rotation_axis * rotation_speed;

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleBodyForce");

    // Loop thru materials
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {

      // Get the material ID
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      new_dw->allocateAndPut(pBodyForceAcc,
                             d_mpm_labels->pBodyForceAccLabel_preReloc,
                             pset);

      // Create space for particle coriolis importance
      ParticleVariable<double> pCoriolisImportance;
      new_dw->allocateAndPut(pCoriolisImportance,
                             d_mpm_labels->pCoriolisImportanceLabel_preReloc,
                             pset);

      // Don't do much if coord rotation is off
      if (!d_mpm_flags->d_useCoordRotation) {

        // Iterate over the particles
        for (auto pidx : *pset) {

          // Compute the body force acceleration (g)
          pBodyForceAcc[pidx] = d_mpm_flags->d_gravity;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[pidx] = 0.0;
        } // particle loop

      } else { // Use coordinate rotation

        // Get the particle data
        constParticleVariable<Point> pPosition;
        old_dw->get(pPosition, d_mpm_labels->pXLabel, pset);

        constParticleVariable<Vector> pVelocity;
        old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);

        for (auto pidx : *pset) {

          // Compute the Coriolis term (omega x v)
          Vector coriolis_accel = Uintah::Cross(omega, pVelocity[pidx]) * 2.0;

          // Compute the centrifugal term (omega x omega x r)
          // Simplified version where body ref point is not needed
          Vector rVec              = pPosition[pidx] - rotation_center;
          Vector omega_x_r         = Uintah::Cross(omega, rVec);
          Vector centrifugal_accel = Uintah::Cross(omega, omega_x_r);

          // Compute the body force acceleration (g - omega x omega x r - 2
          // omega x v)
          pBodyForceAcc[pidx] =
            d_mpm_flags->d_gravity - centrifugal_accel - coriolis_accel;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[pidx] =
            coriolis_accel.length() /
            (centrifugal_accel.length() + coriolis_accel.length());

        } // particle loop
      }   // end if coordinate rotation
    }     // matl loop
  }       // patch loop
}

void
ImpMPM::scheduleApplyExternalLoads(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)

{
  printSchedule(patches, cout_doing, "IMPM::scheduleApplyExternalLoads");
  Task* t =
    scinew Task("IMPM::applyExternalLoads", this, &ImpMPM::applyExternalLoads);

  t->requires(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->requires(Task::OldDW, d_mpm_labels->pExternalForceLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pExternalHeatFluxLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pDispLabel, Ghost::None);
  t->computes(d_mpm_labels->pExtForceLabel_preReloc);
  t->computes(d_mpm_labels->pExternalHeatRateLabel);
  t->computes(d_mpm_labels->pExternalHeatFluxLabel_preReloc);
  if (d_mpm_flags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->pLoadCurveIDLabel_preReloc);
  }

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::applyExternalLoads(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  // Get the current time
  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
  double time = simTimeVar;

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << std::endl;
  }

  // Calculate the force vector at each particle for each bc
  std::vector<double> forceMagPerPart;
  std::vector<PressureBC*> pbcP;
  std::vector<double> heatFluxMagPerPart;
  std::vector<HeatFluxBC*> hfbcP;
  if (d_mpm_flags->d_useLoadCurves) {

    // Currently, only one load curve at a time is supported, but
    // I've left the infrastructure in place to go to multiple
    for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
      auto bcType = particleBC->getType();
      if (bcType == "Pressure") {

        // std::cerr << "Pressure BCs is being supported in ImpMPM" <<
        // std::endl;
        PressureBC* pbc = dynamic_cast<PressureBC*>(particleBC.get());
        pbcP.push_back(pbc);
        forceMagPerPart.push_back(pbc->forcePerParticle(time));

      } else if (bcType == "HeatFlux") {

        HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(particleBC.get());
#if 0
        std::cout << *hfbc << std::endl;
        std::cout << "hfbc type = " << hfbc->getType() << std::endl;
        std::cout << "surface area = " << hfbc->getSurfaceArea() << std::endl;
        std::cout << "heat flux = " << hfbc->heatflux(time) << std::endl;
        std::cout << "flux per particle = " << hfbc->fluxPerParticle(time) << std::endl;
#endif
        hfbcP.push_back(hfbc);
        heatFluxMagPerPart.push_back(hfbc->fluxPerParticle(time));
      }
    }
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing applyExternalLoads");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Point> pX;
      constParticleVariable<double> pExternalHeatFlux;
      constParticleVariable<Vector> pDisp, pExternalForce;
      constParticleVariable<Matrix3> pDefGrad;
      old_dw->get(pX, d_mpm_labels->pXLabel, pset);
      old_dw->get(pExternalHeatFlux, d_mpm_labels->pExternalHeatFluxLabel, pset);
      old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);
      old_dw->get(pExternalForce, d_mpm_labels->pExternalForceLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

      ParticleVariable<double> pExternalHeatFlux_new, pExtHeatRate;
      ParticleVariable<Vector> pExternalForce_new;
      new_dw->allocateAndPut(pExternalHeatFlux_new,
                             d_mpm_labels->pExternalHeatFluxLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pExtHeatRate,
                             d_mpm_labels->pExternalHeatRateLabel,
                             pset);
      new_dw->allocateAndPut(pExternalForce_new,
                             d_mpm_labels->pExtForceLabel_preReloc,
                             pset);

      for (auto idx : *pset) {
        pExternalHeatFlux_new[idx] = 0.;
        pExtHeatRate[idx]          = 0.0;
#if 0
        // Prescribe an external heat rate to some particles
        if(pX[idx].x()*pX[idx].x() + pX[idx].y()*pX[idx].y() > 0.0562*0.0562 ||
           pX[idx].z()>.0562 || pX[idx].z()<-.0562){
          pExtHeatRate[idx]=0.001;
        }
#endif
      }

      if (d_mpm_flags->d_useLoadCurves) {
        bool do_PressureBCs = false;
        for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          auto bcType = particleBC->getType();
          if (bcType == "Pressure") {
            do_PressureBCs = true;
          }
        }

        // Get the load curve data
        constParticleVariable<int> pLoadCurveID;
        old_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

        if (do_PressureBCs) {
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = pExternalForce[idx];
            } else {
              PressureBC* pbc         = pbcP[loadCurveID];
              double force            = forceMagPerPart[loadCurveID];
              pExternalForce_new[idx] = pbc->getForceVector(pX[idx],
                                                            pDisp[idx],
                                                            force,
                                                            time,
                                                            pDefGrad[idx]);
            }
          }
        } else {
          for (auto idx : *pset) {
            pExternalForce_new[idx] =
              pExternalForce[idx] * d_mpm_flags->d_forceIncrementFactor;
          }
        } // end do_PressureBCs

        if (!heatFluxMagPerPart.empty()) {
          // double mag = heatFluxMagPerPart[0];
          // std::cout << "heat flux mag = " << mag << std::endl;
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalHeatFlux_new[idx] = 0.;
            } else {
              //              pExternalHeatFlux_new[idx] = mag;
              pExternalHeatFlux_new[idx] = pExternalHeatFlux[idx];
            }
          }
        }

        // Recycle the loadCurveIDs
        ParticleVariable<int> pLoadCurveID_new;
        new_dw->allocateAndPut(pLoadCurveID_new,
                               d_mpm_labels->pLoadCurveIDLabel_preReloc,
                               pset);
        pLoadCurveID_new.copyData(pLoadCurveID);

      } else { // not use pLoadCurve

        for (auto idx : *pset) {
          pExternalForce_new[idx] =
            pExternalForce[idx] * d_mpm_flags->d_forceIncrementFactor;
          pExternalHeatFlux_new[idx] = pExternalHeatFlux[idx];
        }
      }
    } // matl loop

    printTask(patches, patch, cout_doing, "Completed applyExternalLoads");

  } // patch loop
}

void
ImpMPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSubset* oneMaterial,
                                           const MaterialSet* matls)
{
  printSchedule(patches,
                cout_doing,
                "IMPM::scheduleInterpolateParticlesToGrid");
  Task* t = scinew Task("ImpMPM::interpolateParticlesToGrid",
                        this,
                        &ImpMPM::interpolateParticlesToGrid);

  t->requires(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, d_mpm_labels->pVolumeLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW,
              d_mpm_labels->pTemperatureLabel,
              Ghost::AroundNodes,
              1);
  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW,
              d_mpm_labels->pAccelerationLabel,
              Ghost::AroundNodes,
              1);
  t->requires(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::AroundNodes, 1);

  t->requires(Task::NewDW,
              d_mpm_labels->pBodyForceAccLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pExtForceLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pExternalHeatRateLabel,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pExternalHeatFluxLabel_preReloc,
              Ghost::AroundNodes,
              1);

  if (!d_mpm_flags->d_doGridReset) {
    t->requires(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
    t->computes(d_mpm_labels->gDisplacementLabel);
  }

  t->requires(Task::OldDW,
              d_mpm_labels->NC_CCweightLabel,
              oneMaterial,
              Ghost::AroundCells,
              1);

  if (!d_mpm_flags->d_tempSolve) {
    t->requires(Task::OldDW,
                d_mpm_labels->gTemperatureLabel,
                oneMaterial,
                Ghost::None,
                0);
  }

  t->computes(d_mpm_labels->gMassLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpm_labels->gVolumeLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);

  t->computes(d_mpm_labels->gMassLabel);
  t->computes(d_mpm_labels->gMassAllLabel);
  t->computes(d_mpm_labels->gVolumeLabel);
  t->computes(d_mpm_labels->gVelocityOldLabel);
  t->computes(d_mpm_labels->gVelocityLabel);
  t->computes(d_mpm_labels->dispNewLabel);
  t->computes(d_mpm_labels->gAccelerationLabel);
  t->computes(d_mpm_labels->gBodyForceLabel);
  t->computes(d_mpm_labels->gExternalForceLabel);
  t->computes(d_mpm_labels->gInternalForceLabel);
  t->computes(d_mpm_labels->TotalMassLabel);
  t->computes(d_mpm_labels->gTemperatureLabel, oneMaterial);
  t->computes(d_mpm_labels->gExternalHeatRateLabel);
  t->computes(d_mpm_labels->gExternalHeatFluxLabel);
  t->computes(d_mpm_labels->NC_CCweightLabel, oneMaterial);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);

  d_heatConductionTasks->scheduleProjectHeatSource(sched,
                                                   patches,
                                                   oneMaterial,
                                                   matls);
}

void
ImpMPM::interpolateParticlesToGrid(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;

  static int timestep = 0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing interpolateParticlesToGrid");

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    int i_size        = interpolator->size();
    std::vector<IntVector> ni(i_size);
    std::vector<double> S(i_size);

    NCVariable<double> gMass_global, gVolume_global;
    new_dw->allocateAndPut(gMass_global,
                           d_mpm_labels->gMassLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gVolume_global,
                           d_mpm_labels->gVolumeLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    gMass_global.initialize(d_SMALL_NUM_MPM);
    gVolume_global.initialize(d_SMALL_NUM_MPM);

    NCVariable<double> gMass_sum, gVolume_sum, gSpecificHeat;
    NCVariable<Vector> gVelocity_sum_old, gAcceleration_sum, gBodyForce_sum,
      gExternalForce_sum;
    new_dw->allocateTemporary(gMass_sum, patch, gnone, 0);
    new_dw->allocateTemporary(gVolume_sum, patch, gnone, 0);
    new_dw->allocateTemporary(gSpecificHeat, patch, gnone, 0);
    new_dw->allocateTemporary(gVelocity_sum_old, patch, gnone, 0);
    new_dw->allocateTemporary(gAcceleration_sum, patch, gnone, 0);
    new_dw->allocateTemporary(gBodyForce_sum, patch, gnone, 0);
    new_dw->allocateTemporary(gExternalForce_sum, patch, gnone, 0);
    gMass_sum.initialize(0.);
    gVolume_sum.initialize(0.);
    gSpecificHeat.initialize(0.);
    gVelocity_sum_old.initialize(Vector(0., 0., 0.));
    gAcceleration_sum.initialize(Vector(0., 0., 0.));
    gBodyForce_sum.initialize(Vector(0., 0., 0.));
    gExternalForce_sum.initialize(Vector(0., 0., 0.));

    constNCVariable<double> gTemperature_old;
    NCVariable<double> gTemperature;
    bool switching_to_implicit_from_explicit = false;

    if (!d_mpm_flags->d_tempSolve) {
      if (old_dw->exists(d_mpm_labels->gTemperatureLabel, 0, patch)) {
        old_dw->get(gTemperature_old,
                    d_mpm_labels->gTemperatureLabel,
                    0,
                    patch,
                    gnone,
                    0);
        switching_to_implicit_from_explicit = true;
      }
    }
    new_dw->allocateAndPut(gTemperature,
                           d_mpm_labels->gTemperatureLabel,
                           0,
                           patch);
    gTemperature.initialize(0.0);

    // carry forward interpolation weight
    IntVector low = patch->getExtraNodeLowIndex();
    IntVector hi  = patch->getExtraNodeHighIndex();
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_copy;
    old_dw->get(NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gac, 1);
    new_dw->allocateAndPut(NC_CCweight_copy,
                           d_mpm_labels->NC_CCweightLabel,
                           0,
                           patch);
    NC_CCweight_copy.copyPatch(NC_CCweight, low, hi);

    int numMatls = d_materialManager->getNumMaterials("MPM");
    NCdoubleArray gMass(numMatls), gVolume(numMatls),
      gExternalHeatRate(numMatls), gExternalHeatFlux(numMatls),
      gMass_all(numMatls);
    NCVectorArray gVelocity_old(numMatls), gVelocity(numMatls),
      gAcceleration(numMatls), dispNew(numMatls), gBodyForce(numMatls),
      gExternalForce(numMatls), gInternalForce(numMatls);

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      double Cp            = mpm_matl->getSpecificHeat();
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       1,
                                                       d_mpm_labels->pXLabel);

      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume, pTemperature,
        pExternalHeatRate, pExternalHeatFlux;
      constParticleVariable<Vector> pVelocity, pAcceleration, pBodyForceAcc,
        pExternalForce;
      constParticleVariable<Matrix3> pSize, pDefGrad;

      old_dw->get(pX, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pVolume, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
      old_dw->get(pAcceleration, d_mpm_labels->pAccelerationLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
      new_dw->get(pBodyForceAcc,
                  d_mpm_labels->pBodyForceAccLabel_preReloc,
                  pset);
      new_dw->get(pExternalForce, d_mpm_labels->pExtForceLabel_preReloc, pset);
      new_dw->get(pExternalHeatRate, d_mpm_labels->pExternalHeatRateLabel, pset);
      new_dw->get(pExternalHeatFlux,
                  d_mpm_labels->pExternalHeatFluxLabel_preReloc,
                  pset);

      new_dw->allocateAndPut(gMass[m], d_mpm_labels->gMassLabel, matID, patch);
      new_dw->allocateAndPut(gMass_all[m],
                             d_mpm_labels->gMassAllLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gVolume[m],
                             d_mpm_labels->gVolumeLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gVelocity_old[m],
                             d_mpm_labels->gVelocityOldLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gVelocity[m],
                             d_mpm_labels->gVelocityLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(dispNew[m],
                             d_mpm_labels->dispNewLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gAcceleration[m],
                             d_mpm_labels->gAccelerationLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gBodyForce[m],
                             d_mpm_labels->gBodyForceLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gExternalForce[m],
                             d_mpm_labels->gExternalForceLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gInternalForce[m],
                             d_mpm_labels->gInternalForceLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gExternalHeatRate[m],
                             d_mpm_labels->gExternalHeatRateLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gExternalHeatFlux[m],
                             d_mpm_labels->gExternalHeatFluxLabel,
                             matID,
                             patch);

      if (!d_mpm_flags->d_doGridReset) {
        constNCVariable<Vector> gDisplacement;
        NCVariable<Vector> gDisplacement_new;
        old_dw->get(gDisplacement,
                    d_mpm_labels->gDisplacementLabel,
                    matID,
                    patch,
                    gnone,
                    0);
        new_dw->allocateAndPut(gDisplacement_new,
                               d_mpm_labels->gDisplacementLabel,
                               matID,
                               patch);
        gDisplacement_new.copyData(gDisplacement);
      }

      gMass[m].initialize(d_SMALL_NUM_MPM);
      gVolume[m].initialize(0);
      gExternalHeatRate[m].initialize(0.0);
      gExternalHeatFlux[m].initialize(0.0);
      gVelocity_old[m].initialize(Vector(0, 0, 0));
      gVelocity[m].initialize(Vector(0, 0, 0));
      dispNew[m].initialize(Vector(0, 0, 0));
      gAcceleration[m].initialize(Vector(0, 0, 0));
      gBodyForce[m].initialize(Vector(0, 0, 0));
      gExternalForce[m].initialize(Vector(0, 0, 0));
      gInternalForce[m].initialize(Vector(0, 0, 0));

      double totalMass = 0;
      for (auto idx : *pset) {

        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad[idx]);

        auto pMassAcc = pAcceleration[idx] * pMass[idx];
        auto pMom     = pVelocity[idx] * pMass[idx];
        totalMass += pMass[idx];

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        for (int k = 0; k < i_size; k++) {
          auto node = ni[k];
          if (patch->containsNode(node)) {
            gMass[m][node] += pMass[idx] * S[k];
            gMass_global[node] += pMass[idx] * S[k];
            gVolume[m][node] += pVolume[idx] * S[k];
            gVolume_global[node] += pVolume[idx] * S[k];
            gSpecificHeat[node] += Cp * S[k];
            gExternalHeatRate[m][node] += pExternalHeatRate[idx] * S[k];
            gExternalHeatFlux[m][node] += pExternalHeatFlux[idx] * S[k];
            if (!d_mpm_flags->d_tempSolve) {
              gTemperature[node] += pTemperature[idx] * pMass[idx] * S[k];
            }
            gBodyForce[m][node] += pBodyForceAcc[idx] * pMass[idx] * S[k];
            gExternalForce[m][node] += pExternalForce[idx] * S[k];
            gVelocity_old[m][node] += pMom * S[k];
            gAcceleration[m][node] += pMassAcc * S[k];
          }
        }
      }

      if (mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          gVelocity_old[m][c] /= gMass[m][c];
          gAcceleration[m][c] /= gMass[m][c];
        }
      } else {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          gMass_sum[c] += gMass[m][c];
          gVolume_sum[c] += gVolume[m][c];
          gBodyForce_sum[c] += gBodyForce[m][c];
          gExternalForce_sum[c] += gExternalForce[m][c];
          gVelocity_sum_old[c] += gVelocity_old[m][c];
          gAcceleration_sum[c] += gAcceleration[m][c];
        }
      }

      new_dw->put(sum_vartype(totalMass), d_mpm_labels->TotalMassLabel);
    } // End loop over materials

    // Single material temperature field
    if (!d_mpm_flags->d_tempSolve) {

      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        gTemperature[c] /= gMass_global[c];
      }

    } else {

      // This actually solves for the grid temperatures assuming a linear
      // form of the temperature field.  Uses the particle temperatures
      // as known points within the cell.  For number of particles less
      // than the number of grid cells, take the transpose of the matrix
      // and multiply the original matrix by the transpose.  Also multiply
      // the right hand side by the transpose.  This will yield a matrix
      // that is 8 x 8 and will yield the nodal temperatures even when
      // the number of particles is less than the 8 (number of grid nodes).

      std::vector<IntVector> ni_cell(interpolator->size());

      CellParticleTempMap cell_map;
      CellParticleTempMapArray sparse_cell_map(7);

      for (int m = 0; m < numMatls; m++) {

        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
        int matID            = mpm_matl->getDWIndex();
        ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                         patch,
                                                         Ghost::AroundNodes,
                                                         1,
                                                         d_mpm_labels->pXLabel);

        constParticleVariable<double> pTemperature;
        constParticleVariable<Point> pX;
        constParticleVariable<Matrix3> pSize, pDefGrad;

        old_dw->get(pX, d_mpm_labels->pXLabel, pset);
        old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
        old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
        old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

        for (auto idx : *pset) {

          interpolator->findCellAndWeights(pX[idx],
                                           ni_cell,
                                           S,
                                           pSize[idx],
                                           pDefGrad[idx]);

          ParticleTempShape ptshape;
          ptshape.pTemperature  = pTemperature[idx];
          ptshape.cellNodes     = ni_cell;
          ptshape.shapeFnValues = S;

          IntVector cellID = ni_cell[0];
          cell_map.insert(CellParticleTempPair(cellID, ptshape));
        }
      }

#ifdef debug
      std::cout << "size of cell_map before = " << cell_map.size() << std::endl;
#endif
      for (auto iter = cell_map.begin(); iter != cell_map.end();
           iter      = cell_map.upper_bound(iter->first)) {
#ifdef debug
        std::cout << "cell = " << iter->first
                  << " temp = " << iter->second.particleTemps
                  << " count = " << cell_map.count(iter->first) << std::endl;
#endif

        if (cell_map.count(iter->first) < 8) {
#ifdef debug
          std::cout << "Inserting cell " << iter->first
                    << " into sparse_cell_map\n";
#endif
          auto& smap       = sparse_cell_map[cell_map.count(iter->first) - 1];
          auto eq_range    = cell_map.equal_range(iter->first);
          IntVector cellID = iter->first;
          smap.insert(eq_range.first, eq_range.second);
          cell_map.erase(eq_range.first, eq_range.second);
        }
      }
#ifdef debug
      std::cout << "size of cell_map after = " << cell_map.size() << std::endl;
      for (int i = 0; i < 7; i++) {
        std::cout << "size of sparse_cell_map[" << i
                  << "] after = " << sparse_cell_map[i].size() << std::endl;
      }
#endif

      // Process all of the cells with 8 particles in them
      FastMatrix A(8, 8);
      double B[8];
#ifdef debug
      std::cout << "Working on cells with 8 particles" << std::endl;
#endif
      for (auto iter = cell_map.begin(); iter != cell_map.end();
           iter      = cell_map.upper_bound(iter->first)) {
#ifdef debug
        std::cout << "working on cell " << iter->first << std::endl;
#endif

        auto eq_range = cell_map.equal_range(iter->first);
        int count     = 0;

        ParticleTempShape ptshape;
        for (auto it = eq_range.first; it != eq_range.second; ++it) {
          ptshape  = it->second;
          B[count] = ptshape.pTemperature;
          for (int j = 0; j < 8; j++) {
            A(count, j) = ptshape.shapeFnValues[j];
          }
          count++;
        }

        A.destructiveSolve(B);
        A.zero();
        for (int j = 0; j < 8; j++) {
          if (patch->containsNode(ptshape.cellNodes[j])) {
            gTemperature[ptshape.cellNodes[j]] = B[j];
#ifdef debug
            std::cout << "gTemperature[" << ptshape.cellNodes[j]
                      << "] = " << gTemperature[ptshape.cellNodes[j]]
                      << std::endl;
#endif
          }
        }
      }

      // Work on the cells that have fewer than 8 particles in them
      for (int i = 6; i >= 0; i--) {
#ifdef debug
        std::cout << "Working on cells with " << i + 1 << " particles"
                  << std::endl;
#endif
        auto& smap = sparse_cell_map[i];

        FastMatrix A(8, 8);
        double B[8];
        for (auto it = smap.begin(); it != smap.end();
             it      = smap.upper_bound(it->first)) {
#ifdef debug
          std::cout << "working on cell " << it->first << std::endl;
#endif

          auto eq_range = smap.equal_range(it->first);
          int count     = 0;
          A.zero();
          for (int i = 0; i < 8; i++) {
            B[i] = 0.;
          }
          ParticleTempShape ptshape;
          for (auto it = eq_range.first; it != eq_range.second; ++it) {
            ptshape  = it->second;
            B[count] = ptshape.pTemperature;
            for (int j = 0; j < 8; j++) {
              A(count, j) = ptshape.shapeFnValues[j];
            }
            count++;
          }

          FastMatrix A_t(8, 8);
          A_t.transpose(A);
          double A_tB[8];
          A_t.multiply(B, A_tB);
          FastMatrix A_tA(8, 8);
          A_tA.multiply(A_t, A);

          for (int i = 0; i < 8; i++) {
            if (patch->containsNode(ptshape.cellNodes[i])) {
              if (gTemperature[ptshape.cellNodes[i]] != 0.0) {
#ifdef debug
                std::cout << "i = " << i << " setting gTemperature["
                          << ptshape.cellNodes[i]
                          << "]=" << gTemperature[ptshape.cellNodes[i]]
                          << std::endl;
#endif
                for (int j = 0; j < 8; j++) {
                  A_tA(i, j) = 0.;
                }
                A_tA(i, i) = 1.0;
                A_tB[i]    = gTemperature[ptshape.cellNodes[i]];
              }
            }
          }

          A_tA.destructiveSolve(A_tB);
          for (int j = 0; j < 8; j++) {
            if (patch->containsNode(ptshape.cellNodes[j])) {
              gTemperature[ptshape.cellNodes[j]] = A_tB[j];
#ifdef debug
              std::cout << "gTemperature[" << ptshape.cellNodes[j]
                        << "] = " << gTemperature[ptshape.cellNodes[j]]
                        << std::endl;
#endif
            }
          }
        }
      }
    }
    if (!d_mpm_flags->d_interpolateParticleTempToGridEveryStep) {
      if (timestep > 0) {
        if (!switching_to_implicit_from_explicit) {
          for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
            IntVector c     = *iter;
            gTemperature[c] = gTemperature_old[c];
          }
        }
      }
    }

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector c          = *iter;
          gMass_all[m][c]      = gMass_sum[c];
          gVolume[m][c]        = gVolume_sum[c];
          gBodyForce[m][c]     = gBodyForce_sum[c];
          gExternalForce[m][c] = gExternalForce_sum[c];
          gVelocity_old[m][c] = gVelocity_sum_old[c] / (gMass_sum[c] + 1.e-200);
          gAcceleration[m][c] = gAcceleration_sum[c] / (gMass_sum[c] + 1.e-200);
        }
      }
    } // End loop over materials
  }   // End loop over patches
  timestep++;
}

/*!----------------------------------------------------------------------
 * scheduleFindSurfaceParticles
 *-----------------------------------------------------------------------*/
void
ImpMPM::scheduleFindSurfaceParticles(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFindSurfaceParticles");

  Task* t = scinew Task("ImpMPM::findSurfaceParticles",
                        this,
                        &ImpMPM::findSurfaceParticles);

  t->requires(Task::OldDW,
              d_mpm_labels->pSurfLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->computes(d_mpm_labels->pSurfLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
ImpMPM::findSurfaceParticles(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  auto numMPMMatls = d_materialManager->getNumMaterials("MPM");

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing findSurfaceParticles");

    for (size_t mat = 0; mat < numMPMMatls; mat++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", mat));
      int matID = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<double> pSurf_old;
      ParticleVariable<double> pSurf;

      old_dw->get(pSurf_old, d_mpm_labels->pSurfLabel, pset);
      new_dw->allocateAndPut(pSurf, d_mpm_labels->pSurfLabel_preReloc, pset);

      // For now carry forward the particle surface data
      for (auto particle : *pset) {
        pSurf[particle] = pSurf_old[particle];
      }
    } // matl loop
  }   // patches
}

void
ImpMPM::scheduleDestroyMatrix(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* matls,
                              bool recursion)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleDestroyMatrix");
  Task* t = scinew Task("ImpMPM::destroyMatrix",
                        this,
                        &ImpMPM::destroyMatrix,
                        recursion);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::destroyMatrix(const ProcessorGroup*,
                      const PatchSubset* /*patches*/,
                      const MaterialSubset*,
                      DataWarehouse* /* old_dw */,
                      DataWarehouse* /* new_dw */,
                      bool recursion)
{
  if (cout_doing.active()) {
    cout_doing << "Doing destroyMatrix \t\t\t\t\t IMPM"
               << "\n";
  }

  d_solver->destroyMatrix(recursion);
}

void
ImpMPM::scheduleCreateMatrix(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleCreateMatrix");
  Task* t = scinew Task("ImpMPM::createMatrix", this, &ImpMPM::createMatrix);

  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::createMatrix(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw)
{
  std::map<int, int> dof_diag;
  d_solver->createLocalToGlobalMapping(UintahParallelComponent::d_myworld,
                                       d_perprocPatches,
                                       patches,
                                       3,
                                       d_mpm_flags->d_8or27);
  int global_offset = 0;
  int numMatls      = d_materialManager->getNumMaterials("MPM");
  int n8or27        = d_mpm_flags->d_8or27;

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::createMatrix");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, n8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;

    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    // set the global offset if this is the first patch
    if (pp == 0) {
      global_offset = l2g[lowIndex];
    }

    CCVariable<int> visited;
    new_dw->allocateTemporary(visited, patch, Ghost::AroundCells, 1);
    visited.initialize(0);

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       1,
                                                       d_mpm_labels->pXLabel);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpm_labels->pXLabel, pset);

      for (auto idx : *pset) {
        IntVector cell, ni[27];
        patch->findCell(pX[idx], cell);
        if (visited[cell] == 0) {
          visited[cell] = 1;
          patch->findNodesFromCell(cell, ni);
          std::vector<int> dof(0);
          int l2g_node_num;
          for (int k = 0; k < n8or27; k++) {
            if (patch->containsNode(ni[k])) {
              // subtract global offset in order to map into array correctly
              l2g_node_num = l2g[ni[k]] - global_offset;
              dof.push_back(l2g_node_num);
              dof.push_back(l2g_node_num + 1);
              dof.push_back(l2g_node_num + 2);
            }
          }
          for (auto dofi : dof) {
            for (auto jj = 0u; jj < dof.size(); ++jj) {
              dof_diag[dofi] += 1;
            }
          }
        }
      }
    }
  }
  d_solver->createMatrix(UintahParallelComponent::d_myworld, dof_diag);
}

void
ImpMPM::scheduleApplyBoundaryConditions(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleApplyBoundaryConditions");
  Task* t = scinew Task("ImpMPM::applyBoundaryCondition",
                        this,
                        &ImpMPM::applyBoundaryConditions);

  t->modifies(d_mpm_labels->gVelocityOldLabel);
  t->modifies(d_mpm_labels->gAccelerationLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::applyBoundaryConditions(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* /*old_dw*/,
                                DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing applyBoundaryConditions");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;

    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    // Apply grid boundary conditions to the velocity before storing the data
    IntVector offset = IntVector(0, 0, 0);
    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      NCVariable<Vector> gAcceleration, gVelocity_old;
      new_dw->getModifiable(gVelocity_old,
                            d_mpm_labels->gVelocityOldLabel,
                            matID,
                            patch);
      new_dw->getModifiable(gAcceleration,
                            d_mpm_labels->gAccelerationLabel,
                            matID,
                            patch);

      for (auto face = Patch::startFace; face <= Patch::endFace;
           face      = Patch::nextFace(face)) {

        if (patch->getBCType(face) == Patch::None) {

          int numChildren =
            patch->getBCDataArray(face)->getNumberChildren(matID);
          for (int child = 0; child < numChildren; child++) {

            Iterator nbound_ptr;
            Iterator nu; // not used;

            BoundCondBaseSP vel_bcs = patch->getArrayBCValues(face,
                                                             matID,
                                                             "Velocity",
                                                             nu,
                                                             nbound_ptr,
                                                             child);

            auto bc = std::dynamic_pointer_cast<BoundCond<Vector>>(vel_bcs);

            if (bc != 0) {
              if (bc->getBCType() == "Dirichlet") {
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] = bc->getValue();
                  gAcceleration[*nbound_ptr] = bc->getValue();
                }
                IntVector l, h;
                patch->getFaceNodes(face, 0, l, h);
                for (NodeIterator it(l, h); !it.done(); it++) {
                  IntVector n      = *it;
                  int l2g_node_num = l2g[n];
                  d_solver->d_DOF.insert(l2g_node_num);
                  d_solver->d_DOF.insert(l2g_node_num + 1);
                  d_solver->d_DOF.insert(l2g_node_num + 2);
                }
              }
            }

            BoundCondBaseSP sym_bcs = patch->getArrayBCValues(face,
                                                             matID,
                                                             "Symmetric",
                                                             nu,
                                                             nbound_ptr,
                                                             child);
            auto sbc = std::dynamic_pointer_cast<BoundCond<NoValue>>(sym_bcs);

            if (sbc != 0) {
              if (face == Patch::xplus || face == Patch::xminus) {
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] =
                    Vector(0.,
                           gVelocity_old[*nbound_ptr].y(),
                           gVelocity_old[*nbound_ptr].z());
                  gAcceleration[*nbound_ptr] =
                    Vector(0.,
                           gAcceleration[*nbound_ptr].y(),
                           gAcceleration[*nbound_ptr].z());
                }
              }
              if (face == Patch::yplus || face == Patch::yminus) {
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] =
                    Vector(gVelocity_old[*nbound_ptr].x(),
                           0.,
                           gVelocity_old[*nbound_ptr].z());
                  gAcceleration[*nbound_ptr] =
                    Vector(gAcceleration[*nbound_ptr].x(),
                           0.,
                           gAcceleration[*nbound_ptr].z());
                }
              }
              if (face == Patch::zplus || face == Patch::zminus) {
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] =
                    Vector(gVelocity_old[*nbound_ptr].x(),
                           gVelocity_old[*nbound_ptr].y(),
                           0.);
                  gAcceleration[*nbound_ptr] =
                    Vector(gAcceleration[*nbound_ptr].x(),
                           gAcceleration[*nbound_ptr].y(),
                           0.);
                }
              }
              IntVector l, h;
              patch->getFaceNodes(face, 0, l, h);
              for (NodeIterator it(l, h); !it.done(); it++) {
                IntVector n = *it;
                // The DOF is an IntVector which is initially (0,0,0).
                // Inserting a 1 into any of the components indicates that
                // the component should be inserted into the DOF array.
                IntVector DOF(0, 0, 0);
                if (face == Patch::xminus || face == Patch::xplus) {
                  DOF = IntVector(std::max(DOF.x(), 1),
                                  std::max(DOF.y(), 0),
                                  std::max(DOF.z(), 0));
                }
                if (face == Patch::yminus || face == Patch::yplus) {
                  DOF = IntVector(std::max(DOF.x(), 0),
                                  std::max(DOF.y(), 1),
                                  std::max(DOF.z(), 0));
                }
                if (face == Patch::zminus || face == Patch::zplus) {
                  DOF = IntVector(std::max(DOF.x(), 0),
                                  std::max(DOF.y(), 0),
                                  std::max(DOF.z(), 1));
                }

                int l2g_node_num = l2g[n];
                if (DOF.x()) {
                  d_solver->d_DOF.insert(l2g_node_num);
                }
                if (DOF.y()) {
                  d_solver->d_DOF.insert(l2g_node_num + 1);
                }
                if (DOF.z()) {
                  d_solver->d_DOF.insert(l2g_node_num + 2);
                }
              }
            } // endif (sbc)
          }   // endfor child
        }     // endif patchBC
      }       // endfor face
    }         // endfor mpmmat
  }           // endfor patch
}

void
ImpMPM::scheduleComputeContact(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeContact");
  Task* t =
    scinew Task("ImpMPM::computeContact", this, &ImpMPM::computeContact);

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  if (d_rigidContact) {
    t->requires(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpm_labels->gVelocityOldLabel, Ghost::None);
    t->modifies(d_mpm_labels->dispNewLabel);
    if (!d_mpm_flags->d_doGridReset) {
      t->requires(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
      t->modifies(d_mpm_labels->gDisplacementLabel);
    }
  }
  t->computes(d_mpm_labels->gContactLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeContact(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "Doing computeContact on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n";
    }

    delt_vartype dt;

    int numMatls = d_materialManager->getNumMaterials("MPM");
    std::vector<NCVariable<int>> contact(numMatls);
    for (int n = 0; n < numMatls; n++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", n));
      int matID = mpm_matl->getDWIndex();
      new_dw->allocateAndPut(contact[n],
                             d_mpm_labels->gContactLabel,
                             matID,
                             patch);
      contact[n].initialize(0);
      // std::cout << "matID = " << matID << " n = " << n << "\n";
    }

    if (d_rigidContact) {
      // std::cout << "Rigid = " << std::boolalpha << d_rigidContact << "\n";
      constNCVariable<Vector> vel_rigid;
      constNCVariable<double> mass_rigid;
      int numMatls = d_materialManager->getNumMaterials("MPM");
      for (int n = 0; n < numMatls; n++) {
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", n));
        if (mpm_matl->getIsRigid()) {
          int matID = mpm_matl->getDWIndex();
          new_dw->get(vel_rigid,
                      d_mpm_labels->gVelocityOldLabel,
                      matID,
                      patch,
                      Ghost::None,
                      0);
          new_dw->get(mass_rigid,
                      d_mpm_labels->gMassLabel,
                      matID,
                      patch,
                      Ghost::None,
                      0);
          // std::cout << "matID = " << matID << " n = " << n << "\n";
        }
      }

      // Get and modify non-rigid data
      for (int m = 0; m < numMatls; m++) {
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
        int matID = mpm_matl->getDWIndex();
        NCVariable<Vector> dispNew;
        new_dw->getModifiable(dispNew, d_mpm_labels->dispNewLabel, matID, patch);

        delt_vartype dt;
        old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          if (!compare(mass_rigid[node], 0.0)) {
            // std::cout << "node = " << node << "\n";
            dispNew[node] =
              Vector(vel_rigid[node].x() * dt * d_contactDirections.x(),
                     vel_rigid[node].y() * dt * d_contactDirections.y(),
                     vel_rigid[node].z() * dt * d_contactDirections.z());
            contact[m][node] = 2;
          }
        }

        if (!d_mpm_flags->d_doGridReset) {
          NCVariable<Vector> gDisplacement_new;
          constNCVariable<Vector> gDisplacement_old;
          new_dw->get(gDisplacement_old,
                      d_mpm_labels->gDisplacementLabel,
                      matID,
                      patch,
                      Ghost::None,
                      0);
          new_dw->getModifiable(gDisplacement_new,
                                d_mpm_labels->gDisplacementLabel,
                                matID,
                                patch);
          for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
            IntVector node          = *iter;
            gDisplacement_new[node] = gDisplacement_old[node] + dispNew[node];
          }
        }
      } // endfor numMatls
    }   // endif rigid_body
  }     // endfor patches
}

void
ImpMPM::scheduleFindFixedDOF(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFindFixedDOF");
  Task* t = scinew Task("ImpMPM::findFixedDOF", this, &ImpMPM::findFixedDOF);

  t->requires(Task::NewDW, d_mpm_labels->gMassAllLabel, Ghost::None, 0);
  t->requires(Task::NewDW, d_mpm_labels->gContactLabel, Ghost::None, 0);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::findFixedDOF(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse*,
                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing ImpMPM::findFixedDOF");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    bool firstTimeThrough = true;
    int numMatls          = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid() && firstTimeThrough) {
        firstTimeThrough = false;
        int matID        = mpm_matl->getDWIndex();
        constNCVariable<double> mass;
        constNCVariable<int> contact;
        new_dw
          ->get(mass, d_mpm_labels->gMassAllLabel, matID, patch, Ghost::None, 0);
        new_dw->get(contact,
                    d_mpm_labels->gContactLabel,
                    matID,
                    patch,
                    Ghost::None,
                    0);

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node   = *iter;
          int l2g_node_num = l2g[node];

          // Just look on the grid to see if the gMass is 0 and then remove that
          if (compare(mass[node], 0.)) {
            d_solver->d_DOF.insert(l2g_node_num);
            d_solver->d_DOF.insert(l2g_node_num + 1);
            d_solver->d_DOF.insert(l2g_node_num + 2);
          }
          if (contact[node] == 2) { // Rigid Contact imposed on these nodes
            for (int i = 0; i < 3; i++) {
              if (d_contactDirections[i] == 1) {
                d_solver->d_DOF.insert(l2g_node_num +
                                       i); // specifically, these DOFs
              }
            }
          } // contact ==2
        }   // node iterator
      }     // if not rigid
    }       // loop over matls
  }         // patches
}

void
ImpMPM::scheduleIterate(SchedulerP& sched,
                        const LevelP& level,
                        const PatchSet* patches,
                        const MaterialSet*)
{
  d_recompileSubsched = true;
  printSchedule(patches, cout_doing, "IMPM::scheduleIterate");
  Task* task = scinew Task("ImpMPM::iterate",
                           this,
                           &ImpMPM::iterate,
                           level,
                           sched.get_rep());

  task->hasSubScheduler();

  task->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None, 0);
  task->requires(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None, 0);
  task->requires(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None, 0);
  task->requires(Task::OldDW, d_mpm_labels->pVolumeLabel, Ghost::None, 0);
  task->requires(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None, 0);

  task->modifies(d_mpm_labels->dispNewLabel);
  task->modifies(d_mpm_labels->gVelocityLabel);
  task->modifies(d_mpm_labels->gInternalForceLabel);
  if (!d_mpm_flags->d_doGridReset) {
    task->requires(Task::OldDW,
                   d_mpm_labels->gDisplacementLabel,
                   Ghost::None,
                   0);
    task->modifies(d_mpm_labels->gDisplacementLabel);
  }

  task->requires(Task::NewDW, d_mpm_labels->gVelocityOldLabel, Ghost::None, 0);
  task->requires(Task::NewDW, d_mpm_labels->gMassAllLabel, Ghost::None, 0);
  task->requires(Task::NewDW, d_mpm_labels->gBodyForceLabel, Ghost::None, 0);
  task->requires(Task::NewDW, d_mpm_labels->gExternalForceLabel, Ghost::None, 0);
  task->requires(Task::NewDW, d_mpm_labels->gAccelerationLabel, Ghost::None, 0);
  task->requires(Task::NewDW, d_mpm_labels->gContactLabel, Ghost::None, 0);

  if (d_mpm_flags->d_doMechanics) {
    int numMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

      d_defGradComputer->addComputesAndRequires(task,
                                                mpm_matl,
                                                patches,
                                                true,
                                                false);

      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
      cm->addComputesAndRequires(task, mpm_matl, patches, true, false);
    }
  }

  task->requires(Task::OldDW, d_mpm_labels->delTLabel);

  task->setType(Task::OncePerProc);
  sched->addTask(task, patches, d_materialManager->allMaterials());
}

void
ImpMPM::iterate(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset*,
                DataWarehouse* old_dw,
                DataWarehouse* new_dw,
                LevelP level,
                Scheduler* sched)
{
  Ghost::GhostType gnone = Ghost::None;

  auto old_dw_scrubmode = old_dw->setScrubbing(DataWarehouse::ScrubNone);

  GridP grid = level->getGrid();

  /*
  std::cout << __FILE__ << ":" << __LINE__ << "\n";
  std::cout << "old_dw = " << old_dw << "\n";
  old_dw->print();
  */

  d_subsched->setParentDWs(old_dw, new_dw);
  d_subsched->advanceDataWarehouse(grid);
  d_subsched->setInitTimestep(true);
  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  if (d_recompileSubsched) {
    int numOldDW = 3, numNewDW = 1;
    d_subsched->initialize(numOldDW, numNewDW);

    // This task only zeros out the stiffness matrix it doesn't free any memory.
    scheduleDestroyMatrix(d_subsched, d_perprocPatches, matls, true);

    if (d_mpm_flags->d_doMechanics) {
      scheduleComputeDeformationGradient(d_subsched,
                                         d_perprocPatches,
                                         matls,
                                         true);
      scheduleComputeStressTensor(d_subsched, d_perprocPatches, matls, true);
      scheduleFormStiffnessMatrix(d_subsched, d_perprocPatches, matls);
      scheduleComputeInternalForce(d_subsched, d_perprocPatches, matls);
      scheduleFormQ(d_subsched, d_perprocPatches, matls);
      scheduleSolveForDuCG(d_subsched, d_perprocPatches, matls);
    }

    scheduleGetDisplacementIncrement(d_subsched, d_perprocPatches, matls);
    scheduleUpdateGridKinematics(d_subsched, d_perprocPatches, matls);
    scheduleCheckConvergence(d_subsched, level, d_perprocPatches, matls);

    d_subsched->compile();
    d_recompileSubsched = false;
  }

  int count     = 0;
  bool dispInc  = false;
  bool dispIncQ = false;
  sum_vartype dispIncQNorm, dispIncNorm, dispIncQNorm0, dispIncNormMax;

  // Get all of the required particle data that is in the old_dw and put it
  // in the subscheduler's  new_dw.  Then once dw is advanced, subscheduler
  // will be pulling data out of the old_dw.
  auto subsched_parent_old_dw = d_subsched->get_dw(0);
  auto subsched_new_dw        = d_subsched->get_dw(3);
  // std::cout << "subsched: parent_old_dw = " << subsched_parent_old_dw <<
  // "\n"; std::cout << "subsched: new_dw = " << subsched_new_dw << "\n";

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::iterate-----------------------");

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset =
        subsched_parent_old_dw->getParticleSubset(matID, patch);

      delt_vartype dt;
      old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());

      subsched_new_dw->put(sum_vartype(0.0), d_mpm_labels->dispIncQNorm0);
      subsched_new_dw->put(sum_vartype(0.0), d_mpm_labels->dispIncNormMax);

      // New data to be stored in the subscheduler
      NCVariable<Vector> dispNew, dispNew_sub;
      new_dw->getModifiable(dispNew, d_mpm_labels->dispNewLabel, matID, patch);
      subsched_new_dw->allocateAndPut(dispNew_sub,
                                      d_mpm_labels->dispNewLabel,
                                      matID,
                                      patch);
      dispNew_sub.copyData(dispNew);

      if (!d_mpm_flags->d_doGridReset) {
        NCVariable<Vector> gDisplacement, gDisplacement_new;
        new_dw->getModifiable(gDisplacement,
                              d_mpm_labels->gDisplacementLabel,
                              matID,
                              patch);
        subsched_new_dw->allocateAndPut(gDisplacement_new,
                                        d_mpm_labels->gDisplacementLabel,
                                        matID,
                                        patch);
        gDisplacement_new.copyData(gDisplacement);
      }

      subsched_new_dw->saveParticleSubset(pset, matID, patch);

      // These variables are ultimately retrieved from the subschedulers
      // old datawarehouse after the advancement of the data warehouse.
      double new_dt;
      new_dt = dt;
      subsched_new_dw->put(delt_vartype(new_dt), d_mpm_labels->delTLabel);
    }
  }

  subsched_new_dw->finalize();
  d_subsched->advanceDataWarehouse(grid);
  d_subsched->setInitTimestep(false);

  d_numIterations = 0;
  while (!(dispInc && dispIncQ)) {
    proc0cout << "    Beginning Iteration = " << count << "\n";

    auto subsched_old_dw = d_subsched->get_dw(2);
    subsched_new_dw      = d_subsched->get_dw(3);

    /*
    std::cout << __FILE__ << ":" << __LINE__ << "\n";
    std::cout << "Before scrub\n";
    subsched_old_dw->print();
    std::cout << "subsched: parent_old_dw = " << d_subsched->get_dw(0) << "\n";
    std::cout << "subsched: parent_new_dw = " << d_subsched->get_dw(1) << "\n";
    std::cout << "subsched: old_dw = " << d_subsched->get_dw(2) << "\n";
    std::cout << "subsched: new_dw = " << d_subsched->get_dw(3) << "\n";
    */

    count++;
    subsched_old_dw->setScrubbing(DataWarehouse::ScrubComplete);
    subsched_new_dw->setScrubbing(DataWarehouse::ScrubNone);

    d_subsched->execute(); // THIS ACTUALLY GETS THE WORK DONE
    subsched_new_dw->get(dispIncNorm, d_mpm_labels->dispIncNorm);
    subsched_new_dw->get(dispIncQNorm, d_mpm_labels->dispIncQNorm);
    subsched_new_dw->get(dispIncNormMax, d_mpm_labels->dispIncNormMax);
    subsched_new_dw->get(dispIncQNorm0, d_mpm_labels->dispIncQNorm0);

    double frac_Norm  = dispIncNorm / (dispIncNormMax + 1.e-100);
    double frac_QNorm = dispIncQNorm / (dispIncQNorm0 + 1.e-100);

    if (UintahParallelComponent::d_myworld->myRank() == 0) {
      std::cerr << "  dispIncNorm/dispIncNormMax = " << frac_Norm << "\n";
      std::cerr << "  dispIncQNorm/dispIncQNorm0 = " << frac_QNorm << "\n";
    }
    if ((frac_Norm <= d_mpm_flags->d_convCritDisp) ||
        (dispIncNormMax <= d_mpm_flags->d_convCritDisp)) {
      dispInc = true;
    }
    if ((frac_QNorm <= d_mpm_flags->d_convCritEnergy) ||
        (dispIncQNorm0 <= d_mpm_flags->d_convCritEnergy)) {
      dispIncQ = true;
    }

    // Check to see if the residual is likely a nan, if so, we'll restart.
    bool restart_nan = false;
    if ((std::isnan(dispIncQNorm / dispIncQNorm0) ||
         std::isnan(dispIncNorm / dispIncNormMax)) &&
        dispIncQNorm0 != 0.) {
      restart_nan = true;
      if (UintahParallelComponent::d_myworld->myRank() == 0) {
        std::cerr << "Restarting due to a nan residual\n";
      }
    }

    bool restart_neg_residual = false;
    if (dispIncQNorm / (dispIncQNorm0 + 1e-100) < 0. ||
        dispIncNorm / (dispIncNormMax + 1e-100) < 0.) {
      restart_neg_residual = true;
      if (UintahParallelComponent::d_myworld->myRank() == 0) {
        std::cerr << "Restarting due to a negative residual\n";
      }
    }

    bool restart_num_iters = false;
    if (count > d_mpm_flags->d_maxNumIterations) {
      restart_num_iters = true;
      if (UintahParallelComponent::d_myworld->myRank() == 0) {
        std::cerr << "Restarting due to exceeding max number of iterations\n";
      }
    }

    if (restart_nan || restart_neg_residual || restart_num_iters) {
      new_dw->put(bool_or_vartype(true), VarLabel::find(abortTimeStep_name));
      new_dw->put(bool_or_vartype(true),
                  VarLabel::find(recomputeTimeStep_name));
      return;
    }

    d_subsched->advanceDataWarehouse(grid);

  } // endwhile

  d_numIterations = count;

  // Move the particle data from subscheduler to scheduler.
  auto subsched_old_dw = d_subsched->get_dw(2);
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing << "  Getting the recursive data on patch " << patch->getID()
                 << "\t\t\t IMPM"
                 << "\n"
                 << "\n";
    }

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      // Needed in computeAcceleration
      constNCVariable<Vector> gVelocity, dispNew, gInternalForce;
      subsched_old_dw
        ->get(gVelocity, d_mpm_labels->gVelocityLabel, matID, patch, gnone, 0);
      subsched_old_dw
        ->get(dispNew, d_mpm_labels->dispNewLabel, matID, patch, gnone, 0);
      if (d_mpm_flags->d_doMechanics) {
        subsched_old_dw->get(gInternalForce,
                             d_mpm_labels->gInternalForceLabel,
                             matID,
                             patch,
                             gnone,
                             0);
      }

      NCVariable<Vector> gVelocity_new, dispNew_new, gInternalForce_new;
      new_dw->getModifiable(gVelocity_new,
                            d_mpm_labels->gVelocityLabel,
                            matID,
                            patch);
      new_dw->getModifiable(dispNew_new,
                            d_mpm_labels->dispNewLabel,
                            matID,
                            patch);
      new_dw->getModifiable(gInternalForce_new,
                            d_mpm_labels->gInternalForceLabel,
                            matID,
                            patch);
      gVelocity_new.copyData(gVelocity);
      dispNew_new.copyData(dispNew);
      if (d_mpm_flags->d_doMechanics) {
        gInternalForce_new.copyData(gInternalForce);
      }
    }
  }
  old_dw->setScrubbing(old_dw_scrubmode);
}

void
ImpMPM::scheduleComputeDeformationGradient(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls,
                                           bool recursion)
{
  printSchedule(patches,
                cout_doing,
                "ImpMPM::scheduleComputeDeformationGradient");
  Task* t = scinew Task("ImpMPM::computeDeformationGradient",
                        this,
                        &ImpMPM::computeDeformationGradient,
                        recursion);

  // std::cout << "OldDW = " << Task::OldDW << __FILE__ << __LINE__ << "\n";
  // std::cout << "ParentOldDW = " << Task::ParentOldDW << __FILE__ << __LINE__
  // << "\n"; t->requires(Task::ParentOldDW, d_mpm_labels->delTLabel);
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    d_defGradComputer->addComputesAndRequires(t,
                                              mpm_matl,
                                              patches,
                                              recursion,
                                              true);
  }

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeDeformationGradient(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   bool recursion)
{
  printTask(patches,
            patches->get(0),
            cout_doing,
            "Doing ImpMPM::computeDeformationGradient with recursion");
  d_defGradComputer->computeDeformationGradient(patches,
                                                old_dw,
                                                new_dw,
                                                recursion);
}

void
ImpMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls,
                                    bool recursion)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeStressTensor");
  Task* t = scinew Task("ImpMPM::computeStressTensor",
                        this,
                        &ImpMPM::computeStressTensorImplicit,
                        recursion);

  t->requires(Task::ParentOldDW, d_mpm_labels->delTLabel);
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches, recursion, true);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeStressTensorImplicit(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw,
                                    bool recursion)
{
  if (cout_doing.active()) {
    cout_doing << "Doing computeStressTensor (wrapper) "
               << "\t\t\t IMPM"
               << "\n";
  }

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    ImplicitCM* cmi       = dynamic_cast<ImplicitCM*>(cm);
    if (cmi) {
      cmi->computeStressTensorImplicit(patches,
                                       mpm_matl,
                                       old_dw,
                                       new_dw,
                                       d_solver.get(),
                                       recursion);
    }
  }
}

void
ImpMPM::scheduleComputeDeformationGradient(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  printSchedule(patches,
                cout_doing,
                "IMPM::scheduleComputeDeformationGradientImplicit");
  Task* t      = scinew Task("ImpMPM::computeDeformationGradientImplicit",
                        this,
                        &ImpMPM::computeDeformationGradient);
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    d_defGradComputer->addComputesAndRequires(t, mpm_matl, patches);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeDeformationGradient(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  printTask(patches,
            patches->get(0),
            cout_doing,
            "Doing IMPM::computeDeformationGradient");
  d_defGradComputer->computeDeformationGradient(patches, old_dw, new_dw);
}

void
ImpMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  int numMatls = d_materialManager->getNumMaterials("MPM");
  printSchedule(patches,
                cout_doing,
                "IMPM::scheduleComputeStressTensorImplicit");
  Task* t = scinew Task("ImpMPM::computeStressTensorImplicit",
                        this,
                        &ImpMPM::computeStressTensorImplicit);

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeStressTensorImplicit(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  if (cout_doing.active()) {
    cout_doing << "Doing computeStressTensorImplicit (wrapper)"
               << "\t\t IMPM"
               << "\n";
  }

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->computeStressTensorImplicit(patches, mpm_matl, old_dw, new_dw);
  }
}

void
ImpMPM::scheduleFormStiffnessMatrix(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFormStiffnessMatrix");
  Task* t = scinew Task("ImpMPM::formStiffnessMatrix",
                        this,
                        &ImpMPM::formStiffnessMatrix);

  t->requires(Task::ParentOldDW, d_mpm_labels->delTLabel);
  t->requires(Task::ParentNewDW, d_mpm_labels->gMassAllLabel, Ghost::None);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::formStiffnessMatrix(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* /*old_dw*/,
                            DataWarehouse* new_dw)
{
  if (!d_mpm_flags->d_dynamic) {
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::formStiffnessMatrix");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);

    bool firstTimeThrough = true;
    int numMatls          = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid() && firstTimeThrough) {

        firstTimeThrough = false;
        int matID        = mpm_matl->getDWIndex();
        d_solver->copyL2G(l2g, patch);

        DataWarehouse* parent_old_dw =
          new_dw->getOtherDataWarehouse(Task::ParentOldDW);
        DataWarehouse* parent_new_dw =
          new_dw->getOtherDataWarehouse(Task::ParentNewDW);

        delt_vartype dt;
        constNCVariable<double> gMass;

        parent_old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());
        parent_new_dw->get(gMass,
                           d_mpm_labels->gMassAllLabel,
                           matID,
                           patch,
                           Ghost::None,
                           0);

        double v[1];
        int dof[3];
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node   = *iter;
          int l2g_node_num = l2g[node];
          v[0]             = gMass[node] * (4. / (dt * dt));
          dof[0]           = l2g_node_num;
          dof[1]           = l2g_node_num + 1;
          dof[2]           = l2g_node_num + 2;

          d_solver->fillMatrix(1, &dof[0], 1, &dof[0], v);
          d_solver->fillMatrix(1, &dof[1], 1, &dof[1], v);
          d_solver->fillMatrix(1, &dof[2], 1, &dof[2], v);
        } // end node iterator

      } // endif
    }   // endfor matls
  }
}

void
ImpMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeInternalForce");
  Task* t = scinew Task("ImpMPM::computeInternalForce",
                        this,
                        &ImpMPM::computeInternalForce);

  t->requires(Task::ParentOldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);
  t->requires(Task::ParentOldDW,
              d_mpm_labels->pSizeLabel,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pDefGradLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pStressLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pVolumeLabel_preReloc,
              Ghost::AroundNodes,
              1);

  t->computes(d_mpm_labels->gInternalForceLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeInternalForce(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* /*old_dw*/,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::computeInternalForce");

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    int n8or27      = d_mpm_flags->d_8or27;

    NCVectorArray gInternalForce(numMPMMatls);
    NCVariable<Vector> gInternalForce_sum;
    new_dw->allocateTemporary(gInternalForce_sum, patch, Ghost::None, 0);
    gInternalForce_sum.initialize(Vector(0, 0, 0));

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      new_dw->allocateAndPut(gInternalForce[m],
                             d_mpm_labels->gInternalForceLabel,
                             matID,
                             patch);
      gInternalForce[m].initialize(Vector(0, 0, 0));

      if (!mpm_matl->getIsRigid()) {

        DataWarehouse* parent_old_dw =
          new_dw->getOtherDataWarehouse(Task::ParentOldDW);
        ParticleSubset* pset =
          parent_old_dw->getParticleSubset(matID,
                                           patch,
                                           Ghost::AroundNodes,
                                           1,
                                           d_mpm_labels->pXLabel);

        constParticleVariable<Point> pX;
        constParticleVariable<double> pVolume;
        constParticleVariable<Matrix3> pStress;
        constParticleVariable<Matrix3> pSize;
        constParticleVariable<Matrix3> pDefGrad;

        parent_old_dw->get(pX, d_mpm_labels->pXLabel, pset);
        parent_old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
        new_dw->get(pVolume, d_mpm_labels->pVolumeLabel_preReloc, pset);
        new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel_preReloc, pset);
        new_dw->get(pStress, d_mpm_labels->pStressLabel_preReloc, pset);

        Matrix3 pStressVol;

        for (auto idx : *pset) {

          interpolator->findCellAndShapeDerivatives(pX[idx],
                                                    ni,
                                                    d_S,
                                                    pSize[idx],
                                                    pDefGrad[idx]);
          pStressVol = pStress[idx] * pVolume[idx];
          for (int k = 0; k < n8or27; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              Vector div(d_S[k].x() * oodx[0],
                         d_S[k].y() * oodx[1],
                         d_S[k].z() * oodx[2]);
              gInternalForce[m][node] -= (div * pStress[idx]) * pVolume[idx];
            }
          }
        }

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gInternalForce_sum[node] += gInternalForce[m][node];
        }
      } // if matl isn't rigid
    }   // matls

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node          = *iter;
          gInternalForce[m][node] = gInternalForce_sum[node];
        }
      }
    } // matls
  }   // patches
}

void
ImpMPM::scheduleFormQ(SchedulerP& sched,
                      const PatchSet* patches,
                      const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleFormQ");
  Task* t = scinew Task("ImpMPM::formQ", this, &ImpMPM::formQ);

  Ghost::GhostType gnone = Ghost::None;

  t->requires(Task::ParentOldDW, d_mpm_labels->delTLabel);
  t->requires(Task::ParentNewDW, d_mpm_labels->gMassAllLabel, gnone, 0);
  t->requires(Task::ParentNewDW, d_mpm_labels->gBodyForceLabel, gnone, 0);
  t->requires(Task::ParentNewDW, d_mpm_labels->gExternalForceLabel, gnone, 0);
  t->requires(Task::ParentNewDW, d_mpm_labels->gVelocityOldLabel, gnone, 0);
  t->requires(Task::ParentNewDW, d_mpm_labels->gAccelerationLabel, gnone, 0);
  t->requires(Task::OldDW, d_mpm_labels->dispNewLabel, gnone, 0);
  t->requires(Task::NewDW, d_mpm_labels->gInternalForceLabel, gnone, 0);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::formQ(const ProcessorGroup*,
              const PatchSubset* patches,
              const MaterialSubset*,
              DataWarehouse* old_dw,
              DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::formQ");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    bool firstTimeThrough = true;
    int numMatls          = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid() && firstTimeThrough) {
        firstTimeThrough = false;
        int matID        = mpm_matl->getDWIndex();

        DataWarehouse* parent_new_dw =
          new_dw->getOtherDataWarehouse(Task::ParentNewDW);
        DataWarehouse* parent_old_dw =
          new_dw->getOtherDataWarehouse(Task::ParentOldDW);

        delt_vartype dt;
        constNCVariable<double> gMass;
        constNCVariable<Vector> gBodyForce, gExternalForce, gInternalForce;
        constNCVariable<Vector> dispNew, gVelocity, gAcceleration;

        parent_old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());
        parent_new_dw->get(gBodyForce,
                           d_mpm_labels->gBodyForceLabel,
                           matID,
                           patch,
                           gnone,
                           0);
        parent_new_dw->get(gExternalForce,
                           d_mpm_labels->gExternalForceLabel,
                           matID,
                           patch,
                           gnone,
                           0);
        parent_new_dw->get(gVelocity,
                           d_mpm_labels->gVelocityOldLabel,
                           matID,
                           patch,
                           gnone,
                           0);
        parent_new_dw->get(gAcceleration,
                           d_mpm_labels->gAccelerationLabel,
                           matID,
                           patch,
                           gnone,
                           0);
        parent_new_dw
          ->get(gMass, d_mpm_labels->gMassAllLabel, matID, patch, gnone, 0);

        old_dw->get(dispNew, d_mpm_labels->dispNewLabel, matID, patch, gnone, 0);
        new_dw->get(gInternalForce,
                    d_mpm_labels->gInternalForceLabel,
                    matID,
                    patch,
                    gnone,
                    0);

        double fodts = 4. / (dt * dt);
        double fodt  = 4. / dt;

        double Q = 0.;

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node   = *iter;
          int l2g_node_num = l2g[node];

          Vector force =
            gExternalForce[node] + gInternalForce[node] + gBodyForce[node];
          double v[3];
          v[0] = force.x();
          v[1] = force.y();
          v[2] = force.z();

          // temp2 = M*a^(k-1)(t+dt)
          if (d_mpm_flags->d_dynamic) {
            Vector damping = (dispNew[node] * fodts - gVelocity[node] * fodt -
                              gAcceleration[node]) *
                             gMass[node];
            v[0] -= damping.x();
            v[1] -= damping.y();
            v[2] -= damping.z();
          }
          d_solver->fillVector(l2g_node_num, double(v[0]));
          d_solver->fillVector(l2g_node_num + 1, double(v[1]));
          d_solver->fillVector(l2g_node_num + 2, double(v[2]));
          Q += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        }
        if (std::isnan(Q)) {
          std::cout << "RHS contains a nan, restarting timestep" << std::endl;
          new_dw->put(bool_or_vartype(true),
                      VarLabel::find(abortTimeStep_name));
          new_dw->put(bool_or_vartype(true),
                      VarLabel::find(recomputeTimeStep_name));
          return;
        }
      } // first time through non-rigid
    }   // matls
  }     // patches
}

void
ImpMPM::scheduleSolveForDuCG(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleSolveForDuCG");
  Task* t = scinew Task("ImpMPM::solveForDuCG", this, &ImpMPM::solveForDuCG);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::solveForDuCG(const ProcessorGroup* /*pg*/,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse*,
                     DataWarehouse* new_dw)

{
  if (cout_doing.active()) {
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      cout_doing << "Doing solveForDuCG on patch " << patch->getID()
                 << "\t\t\t\t IMPM"
                 << "\n";
    }
  }

  // if a recompute time step has already been called for don't do the solve
  if (!new_dw->recomputeTimeStep()) {
    d_solver->assembleVector();
    d_solver->removeFixedDOF();
    std::vector<double> guess;
    d_solver->solve(guess);
  } else {
    std::cout << "skipping solve, timestep has already called for a restart"
              << std::endl;
  }
}

void
ImpMPM::scheduleGetDisplacementIncrement(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleGetDisplacementIncrement");
  Task* t = scinew Task("ImpMPM::getDisplacementIncrement",
                        this,
                        &ImpMPM::getDisplacementIncrement);

  t->computes(d_mpm_labels->dispIncLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::getDisplacementIncrement(const ProcessorGroup* /*pg*/,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse*,
                                 DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::getDisplacementIncrement");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    std::vector<double> x;
    int begin = d_solver->getSolution(x);

    int numMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      NCVariable<Vector> dispInc;
      new_dw->allocateAndPut(dispInc, d_mpm_labels->dispIncLabel, matID, patch);
      dispInc.initialize(Vector(0.));

      if (d_mpm_flags->d_doMechanics) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node   = *iter;
          int l2g_node_num = l2g[node] - begin;
          dispInc[node] =
            Vector(x[l2g_node_num], x[l2g_node_num + 1], x[l2g_node_num + 2]);
        }
      }
    }
  }
}

void
ImpMPM::scheduleUpdateGridKinematics(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleUpdateGridKinematics");
  Task* t = scinew Task("ImpMPM::updateGridKinematics",
                        this,
                        &ImpMPM::updateGridKinematics);

  t->requires(Task::ParentOldDW, d_mpm_labels->delTLabel);
  t->requires(Task::ParentNewDW,
              d_mpm_labels->gVelocityOldLabel,
              Ghost::None,
              0);
  t->requires(Task::ParentNewDW, d_mpm_labels->gContactLabel, Ghost::None, 0);
  t->requires(Task::OldDW, d_mpm_labels->dispNewLabel, Ghost::None, 0);
  t->requires(Task::NewDW, d_mpm_labels->dispIncLabel, Ghost::None, 0);
  if (!d_mpm_flags->d_doGridReset) {
    t->requires(Task::ParentOldDW,
                d_mpm_labels->gDisplacementLabel,
                Ghost::None,
                0);
    t->computes(d_mpm_labels->gDisplacementLabel);
  }
  t->computes(d_mpm_labels->dispNewLabel);
  t->computes(d_mpm_labels->gVelocityLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::updateGridKinematics(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::updateGridKinematics");

    int matID_rigid = -99;
    int numMatls    = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (mpm_matl->getIsRigid()) {
        matID_rigid = mpm_matl->getDWIndex();
      }
    }

    constNCVariable<Vector> gVelocity_rigid;
    if (d_rigidContact) {
      DataWarehouse* parent_new_dw =
        new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      parent_new_dw->get(gVelocity_rigid,
                         d_mpm_labels->gVelocityOldLabel,
                         matID_rigid,
                         patch,
                         gnone,
                         0);
    }

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      delt_vartype dt;
      constNCVariable<int> gContact;
      constNCVariable<Vector> dispInc, dispNew_old, gVelocity_old;
      NCVariable<Vector> dispNew, gVelocity;

      DataWarehouse* parent_new_dw =
        new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      DataWarehouse* parent_old_dw =
        new_dw->getOtherDataWarehouse(Task::ParentOldDW);

      parent_old_dw->get(dt, d_mpm_labels->delTLabel, patch->getLevel());
      parent_new_dw->get(gVelocity_old,
                         d_mpm_labels->gVelocityOldLabel,
                         matID,
                         patch,
                         gnone,
                         0);
      parent_new_dw
        ->get(gContact, d_mpm_labels->gContactLabel, matID, patch, gnone, 0);
      old_dw
        ->get(dispNew_old, d_mpm_labels->dispNewLabel, matID, patch, gnone, 0);
      new_dw->get(dispInc, d_mpm_labels->dispIncLabel, matID, patch, gnone, 0);

      new_dw->allocateAndPut(dispNew, d_mpm_labels->dispNewLabel, matID, patch);
      new_dw->allocateAndPut(gVelocity,
                             d_mpm_labels->gVelocityLabel,
                             matID,
                             patch);
      dispNew.copyData(dispNew_old);

      double oneifdyn = 0.;
      if (d_mpm_flags->d_dynamic) {
        oneifdyn = 1.;
      }

      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          dispNew[node] += dispInc[node];
          gVelocity[node] =
            dispNew[node] * (2. / dt) - oneifdyn * gVelocity_old[node];
        }
      }

      if (d_rigidContact) { // overwrite some of the values computed above
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          if (gContact[node] == 2) {
            dispNew[node] = Vector(
              (1. - d_contactDirections.x()) * dispNew[node].x() +
                d_contactDirections.x() * gVelocity_rigid[node].x() * dt,
              (1. - d_contactDirections.y()) * dispNew[node].y() +
                d_contactDirections.y() * gVelocity_rigid[node].y() * dt,
              (1. - d_contactDirections.z()) * dispNew[node].z() +
                d_contactDirections.z() * gVelocity_rigid[node].z() * dt);

            gVelocity[node] =
              dispNew[node] * (2. / dt) - oneifdyn * gVelocity_old[node];
          } // if contact == 2
        }   // for
      }     // if d_rigidContact

      if (!d_mpm_flags->d_doGridReset) {
        constNCVariable<Vector> gDisplacement_old;
        NCVariable<Vector> gDisplacement;
        parent_old_dw->get(gDisplacement_old,
                           d_mpm_labels->gDisplacementLabel,
                           matID,
                           patch,
                           gnone,
                           0);
        new_dw->allocateAndPut(gDisplacement,
                               d_mpm_labels->gDisplacementLabel,
                               matID,
                               patch);
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node      = *iter;
          gDisplacement[node] = gDisplacement_old[node] + dispNew[node];
        }
      }
    } // matls
  }
}

void
ImpMPM::scheduleCheckConvergence(SchedulerP& sched,
                                 const LevelP& /* level */,
                                 const PatchSet* patches,
                                 const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleCheckConvergence");
  Task* t =
    scinew Task("ImpMPM::checkConvergence", this, &ImpMPM::checkConvergence);

  t->requires(Task::OldDW, d_mpm_labels->dispIncQNorm0);
  t->requires(Task::OldDW, d_mpm_labels->dispIncNormMax);
  t->requires(Task::NewDW, d_mpm_labels->dispIncLabel, Ghost::None, 0);

  t->computes(d_mpm_labels->dispIncNormMax);
  t->computes(d_mpm_labels->dispIncQNorm0);
  t->computes(d_mpm_labels->dispIncNorm);
  t->computes(d_mpm_labels->dispIncQNorm);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::checkConvergence(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw)
{
  int global_offset = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::checkConvergence");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, d_mpm_flags->d_8or27);
    auto lowIndex  = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);
    d_solver->copyL2G(l2g, patch);

    int matID = 0;
    constNCVariable<Vector> dispInc;
    new_dw
      ->get(dispInc, d_mpm_labels->dispIncLabel, matID, patch, Ghost::None, 0);

    double dispIncNorm  = 0.;
    double dispIncQNorm = 0.;
    std::vector<double> getQ;

    d_solver->getRHS(getQ);
    if (p == 0) {
      global_offset = l2g[lowIndex];
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node   = *iter;
      int l2g_node_num = l2g[node] - global_offset;
      dispIncNorm += Dot(dispInc[node], dispInc[node]);
      dispIncQNorm += dispInc[node].x() * getQ[l2g_node_num] +
                      dispInc[node].y() * getQ[l2g_node_num + 1] +
                      dispInc[node].z() * getQ[l2g_node_num + 2];
    }

    // We are computing both dispIncQNorm0 and dispIncNormMax (max residuals)
    // We are computing both dispIncQNorm and dispIncNorm (current residuals)
    double dispIncQNorm0, dispIncNormMax;
    sum_vartype dispIncQNorm0_var, dispIncNormMax_var;
    old_dw->get(dispIncQNorm0_var, d_mpm_labels->dispIncQNorm0);
    old_dw->get(dispIncNormMax_var, d_mpm_labels->dispIncNormMax);

    dispIncQNorm0  = dispIncQNorm0_var;
    dispIncNormMax = dispIncNormMax_var;

    bool first_iteration = false;
    if (compare(dispIncQNorm0, 0.)) {
      first_iteration = true;
      dispIncQNorm0   = dispIncQNorm;
    }

    if (dispIncNorm > dispIncNormMax) {
      dispIncNormMax = dispIncNorm;
    }

    // The following is being done because the denominator in the
    // convergence criteria is carried forward each iteration.  Since
    // every patch puts this into the sum_vartype, the value is multiplied
    // by the number of patches.  Predividing by numPatches fixes this.
    int numPatches = patch->getLevel()->numPatches();
    if (!first_iteration) {
      dispIncQNorm0 /= ((double)numPatches);
      if (dispIncNormMax != dispIncNorm) {
        dispIncNormMax /= ((double)numPatches);
      }
    }

    new_dw->put(sum_vartype(dispIncNorm), d_mpm_labels->dispIncNorm);
    new_dw->put(sum_vartype(dispIncQNorm), d_mpm_labels->dispIncQNorm);
    new_dw->put(sum_vartype(dispIncNormMax), d_mpm_labels->dispIncNormMax);
    new_dw->put(sum_vartype(dispIncQNorm0), d_mpm_labels->dispIncQNorm0);
  } // End of loop over patches
}

void
ImpMPM::scheduleUpdateTotalDisplacement(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (!d_mpm_flags->d_doGridReset) {
    printSchedule(patches, cout_doing, "IMPM::scheduleUpdateTotalDisplacement");
    Task* t = scinew Task("ImpMPM::updateTotalDisplacement",
                          this,
                          &ImpMPM::updateTotalDisplacement);

    t->requires(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpm_labels->dispNewLabel, Ghost::None);
    t->modifies(d_mpm_labels->gDisplacementLabel);

    t->setType(Task::OncePerProc);
    sched->addTask(t, patches, matls);
  }
}

void
ImpMPM::updateTotalDisplacement(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::updateTotalDisplacement");

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      constNCVariable<Vector> dispNew, gDisplacement_old;
      NCVariable<Vector> gDisplacement;
      old_dw->get(gDisplacement_old,
                  d_mpm_labels->gDisplacementLabel,
                  matID,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(gDisplacement,
                            d_mpm_labels->gDisplacementLabel,
                            matID,
                            patch);
      new_dw->get(dispNew, d_mpm_labels->dispNewLabel, matID, patch, gnone, 0);

      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node      = *iter;
        gDisplacement[node] = gDisplacement_old[node] + dispNew[node];
      }
    }
  }
}

void
ImpMPM::scheduleComputeAcceleration(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeAcceleration");
  Task* t = scinew Task("ImpMPM::computeAcceleration",
                        this,
                        &ImpMPM::computeAcceleration);

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);
  t->requires(Task::NewDW, d_mpm_labels->gVelocityOldLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->dispNewLabel, Ghost::None);
  if (!d_mpm_flags->d_doGridReset) {
    t->requires(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
    t->modifies(d_mpm_labels->gDisplacementLabel);
  }
  t->modifies(d_mpm_labels->gAccelerationLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeAcceleration(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  if (!d_mpm_flags->d_dynamic) {
    return;
  }
  Ghost::GhostType gnone = Ghost::None;

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::computeAcceleration");

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      constNCVariable<Vector> gVelocity, dispNew;
      NCVariable<Vector> gAcceleration;

      new_dw->get(dispNew, d_mpm_labels->dispNewLabel, matID, patch, gnone, 0);
      new_dw->get(gVelocity,
                  d_mpm_labels->gVelocityOldLabel,
                  matID,
                  patch,
                  gnone,
                  0);
      new_dw->getModifiable(gAcceleration,
                            d_mpm_labels->gAccelerationLabel,
                            matID,
                            patch);

      double fodts = 4. / (delT * delT);
      double fodt  = 4. / (delT);

      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        gAcceleration[node] =
          dispNew[node] * fodts - gVelocity[node] * fodt - gAcceleration[node];
      }
    }
  }
}

void
ImpMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  printSchedule(patches,
                cout_doing,
                "IMPM::scheduleInterpolateToParticlesAndUpdate");
  Task* t = scinew Task("ImpMPM::interpolateToParticlesAndUpdate",
                        this,
                        &ImpMPM::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_mpm_labels->delTLabel);

  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pParticleIDLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pAccelerationLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pTemperatureLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pTempPreviousLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pDispLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);

  t->requires(Task::NewDW, d_mpm_labels->pVolumeLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_mpm_labels->pDefGradLabel_preReloc, Ghost::None);

  t->requires(Task::NewDW,
              d_mpm_labels->gAccelerationLabel,
              Ghost::AroundCells,
              1);
  t->requires(Task::NewDW, d_mpm_labels->dispNewLabel, Ghost::AroundCells, 1);
  t->requires(Task::NewDW,
              d_mpm_labels->gTemperatureRateLabel,
              d_oneMaterial,
              Ghost::AroundCells,
              1);

  t->computes(d_mpm_labels->pVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pAccelerationLabel_preReloc);
  t->computes(d_mpm_labels->pXLabel_preReloc);
  t->computes(d_mpm_labels->pXXLabel);
  t->computes(d_mpm_labels->pParticleIDLabel_preReloc);
  t->computes(d_mpm_labels->pMassLabel_preReloc);
  t->computes(d_mpm_labels->pTemperatureLabel_preReloc);
  t->computes(d_mpm_labels->pDispLabel_preReloc);
  t->computes(d_mpm_labels->pSizeLabel_preReloc);
  t->computes(d_mpm_labels->pTempPreviousLabel_preReloc);

  if (d_mpm_flags->d_artificialViscosity) {
    t->requires(Task::OldDW, d_mpm_labels->p_qLabel, Ghost::None);
    t->computes(d_mpm_labels->p_qLabel_preReloc);
  }

  t->computes(d_mpm_labels->KineticEnergyLabel);
  t->computes(d_mpm_labels->CenterOfMassPositionLabel);
  t->computes(d_mpm_labels->TotalMomentumLabel);
  t->computes(d_mpm_labels->ThermalEnergyLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  Ghost::GhostType gac = Ghost::AroundCells;

  delt_vartype delT;
  old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::interpolateToParticlesAndUpdate");

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and displacement to the particles to update their
    // velocity and position respectively
    Vector disp(0.0, 0.0, 0.0);
    Vector acc(0.0, 0.0, 0.0);

    // DON'T MOVE THESE!!!
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke             = 0;
    double thermal_energy = 0.0;
    int numMPMMatls       = d_materialManager->getNumMaterials("MPM");
    int n8or27            = d_mpm_flags->d_8or27;

    double move_particles = 1.;
    if (!d_mpm_flags->d_doGridReset) {
      move_particles = 0.;
    }

    constNCVariable<double> gTemperatureRate;
    new_dw->get(gTemperatureRate,
                d_mpm_labels->gTemperatureRateLabel,
                0,
                patch,
                gac,
                1);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      double Cp = mpm_matl->getSpecificHeat();

      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume_new, pTemp, pq;
      constParticleVariable<Vector> pVelocity, pAcceleration, pDisp;
      constParticleVariable<Matrix3> pSize, pDefGrad;

      ParticleVariable<double> pMass_new, pTemp_new, pq_new;
      ParticleVariable<double> pTempPre_new;
      ParticleVariable<Point> pX_new, pXx;
      ParticleVariable<Vector> pVelocity_new, pAcceleration_new, pDisp_new;
      ParticleVariable<Matrix3> pSize_new;

      constNCVariable<Vector> dispNew, gAcceleration;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      ParticleSubset* delete_particles = scinew ParticleSubset(0, matID, patch);

      old_dw->get(pX, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pAcceleration, d_mpm_labels->pAccelerationLabel, pset);
      old_dw->get(pTemp, d_mpm_labels->pTemperatureLabel, pset);
      old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->get(pVolume_new, d_mpm_labels->pVolumeLabel_preReloc, pset);
      new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel_preReloc, pset);

      new_dw->allocateAndPut(pMass_new, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pTemp_new,
                             d_mpm_labels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempPre_new,
                             d_mpm_labels->pTempPreviousLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pVelocity_new,
                             d_mpm_labels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pAcceleration_new,
                             d_mpm_labels->pAccelerationLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pX_new, d_mpm_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pXx, d_mpm_labels->pXXLabel, pset);
      new_dw->allocateAndPut(pDisp_new, d_mpm_labels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(pSize_new, d_mpm_labels->pSizeLabel_preReloc, pset);

      new_dw->get(dispNew, d_mpm_labels->dispNewLabel, matID, patch, gac, 1);
      new_dw->get(gAcceleration,
                  d_mpm_labels->gAccelerationLabel,
                  matID,
                  patch,
                  gac,
                  1);

      if (d_mpm_flags->d_artificialViscosity) {
        old_dw->get(pq, d_mpm_labels->p_qLabel, pset);
        new_dw->allocateAndPut(pq_new, d_mpm_labels->p_qLabel_preReloc, pset);
        pq_new.copyData(pq);
      }
      pSize_new.copyData(pSize);
      pTempPre_new.copyData(pTemp);

      for (auto idx : *pset) {

        interpolator->findCellAndWeightsAndShapeDerivatives(pX[idx],
                                                            ni,
                                                            S,
                                                            d_S,
                                                            pSize[idx],
                                                            pDefGrad[idx]);

        disp            = Vector(0.0, 0.0, 0.0);
        acc             = Vector(0.0, 0.0, 0.0);
        double tempRate = 0.;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < n8or27; k++) {
          auto node = ni[k];
          disp += dispNew[node] * S[k];
          acc += gAcceleration[node] * S[k];
          tempRate += gTemperatureRate[node] * S[k];
        }

        // Update the particle's position and velocity
        pX_new[idx]    = pX[idx] + disp * move_particles;
        pDisp_new[idx] = pDisp[idx] + disp;
        pVelocity_new[idx] =
          pVelocity[idx] + (pAcceleration[idx] + acc) * (.5 * delT);

        // pXx is only useful if we're not in normal grid resetting mode.
        pXx[idx] = pX[idx] + pDisp_new[idx];

        pAcceleration_new[idx] = acc;
        pMass_new[idx]         = pMass[idx];
        pTemp_new[idx]         = pTemp[idx] + tempRate * delT;

        if (pMass_new[idx] <= 0.0) {
          delete_particles->addParticle(idx);
          pVelocity_new[idx] = Vector(0., 0., 0);
          pX_new[idx]        = pX[idx];
        }

        thermal_energy += pTemp_new[idx] * pMass[idx] * Cp;
        // Thermal energy due to temperature flux (spatially varying part).
        // thermal_energy2 += potential_energy* pVolume_new[idx];
        ke += .5 * pMass[idx] * pVelocity_new[idx].length2();
        CMX = CMX + (pX_new[idx] * pMass[idx]).asVector();
        totalMom += pVelocity_new[idx] * pMass[idx];
      }

      if (mpm_matl->getIsRigid()) {
        const double tcurr = simTimeVar;
        if (tcurr >= d_contactStopTime) {
          for (auto idx : *pset) {
            pVelocity_new[idx] = d_velocityAfterContactStop;
          }
        }
      }

      new_dw->deleteParticles(delete_particles);

      constParticleVariable<long64> pids;
      ParticleVariable<long64> pids_new;
      old_dw->get(pids, d_mpm_labels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pids_new,
                             d_mpm_labels->pParticleIDLabel_preReloc,
                             pset);
      pids_new.copyData(pids);
    }

    // DON'T MOVE THESE!!!
    new_dw->put(sum_vartype(ke), d_mpm_labels->KineticEnergyLabel);
    new_dw->put(sumvec_vartype(CMX), d_mpm_labels->CenterOfMassPositionLabel);
    new_dw->put(sumvec_vartype(totalMom), d_mpm_labels->TotalMomentumLabel);
    new_dw->put(sum_vartype(thermal_energy), d_mpm_labels->ThermalEnergyLabel);
  }
}

void
ImpMPM::scheduleInterpolateStressToGrid(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  // This task is done for visualization only
  Task* t = scinew Task("ImpMPM::interpolateStressToGrid",
                        this,
                        &ImpMPM::interpolateStressToGrid);

  t->requires(Task::OldDW, d_mpm_labels->pXLabel, Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::AroundNodes, 1);
  t->requires(Task::NewDW,
              d_mpm_labels->pVolumeLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW,
              d_mpm_labels->pStressLabel_preReloc,
              Ghost::AroundNodes,
              1);
  t->requires(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW,
              d_mpm_labels->gVolumeLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain,
              Ghost::None);
  t->requires(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::AroundNodes, 1);

  t->modifies(d_mpm_labels->gInternalForceLabel);
  t->computes(d_mpm_labels->gStressForSavingLabel);

  if (!d_boundaryTractionFaces.empty()) {
    for (auto face : d_boundaryTractionFaces) {
      t->computes(d_mpm_labels->BndyForceLabel[face]);       // node based
      t->computes(d_mpm_labels->BndyContactAreaLabel[face]); // node based
    }
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::interpolateStressToGrid(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  // node based forces
  Vector bndyForce[6];
  double bndyArea[6];
  for (int iface = 0; iface < 6; iface++) {
    bndyForce[iface] = Vector(0.);
    bndyArea[iface]  = 0.;
  }
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing ImpMPM::interpolateStressToGrid");

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    int numMatls = d_materialManager->getNumMaterials("MPM");
    int n8or27   = d_mpm_flags->d_8or27;

    constNCVariable<double> gVolume_sum;
    new_dw->get(gVolume_sum,
                d_mpm_labels->gVolumeLabel,
                d_materialManager->getAllInOneMaterial()->get(0),
                patch,
                Ghost::None,
                0);

    NCVariable<Vector> gInternalForce_sum;
    NCVariable<Matrix3> gStress_sum;
    new_dw->allocateTemporary(gInternalForce_sum, patch, Ghost::None, 0);
    new_dw->allocateTemporary(gStress_sum, patch, Ghost::None, 0);
    gInternalForce_sum.initialize(Vector(0, 0, 0));
    gStress_sum.initialize(Matrix3(0.));

    std::vector<constNCVariable<double>> gVolume(numMatls);
    NCMatrix3Array gStress(numMatls);
    NCVectorArray gInternalForce(numMatls);

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      new_dw->get(gVolume[m],
                  d_mpm_labels->gVolumeLabel,
                  matID,
                  patch,
                  Ghost::None,
                  0);
      new_dw->getModifiable(gInternalForce[m],
                            d_mpm_labels->gInternalForceLabel,
                            matID,
                            patch);
      new_dw->allocateAndPut(gStress[m],
                             d_mpm_labels->gStressForSavingLabel,
                             matID,
                             patch);
      gInternalForce[m].initialize(Vector(0.));
      gStress[m].initialize(Matrix3(0.));

      if (!mpm_matl->getIsRigid()) {

        constParticleVariable<Point> pX;
        constParticleVariable<double> pVolume;
        constParticleVariable<Matrix3> pSize, pStress, pDefGrad;

        ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                         patch,
                                                         Ghost::AroundNodes,
                                                         1,
                                                         d_mpm_labels->pXLabel);
        old_dw->get(pX, d_mpm_labels->pXLabel, pset);
        new_dw->get(pVolume, d_mpm_labels->pVolumeLabel_preReloc, pset);
        old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
        old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
        new_dw->get(pStress, d_mpm_labels->pStressLabel_preReloc, pset);

        Matrix3 pStressVol;
        for (auto idx : *pset) {

          interpolator->findCellAndWeightsAndShapeDerivatives(pX[idx],
                                                              ni,
                                                              S,
                                                              d_S,
                                                              pSize[idx],
                                                              pDefGrad[idx]);
          pStressVol = pStress[idx] * pVolume[idx];
          for (int k = 0; k < n8or27; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              gStress[m][node] += pStressVol * S[k];
              Vector div(d_S[k].x() * oodx[0],
                         d_S[k].y() * oodx[1],
                         d_S[k].z() * oodx[2]);
              gInternalForce[m][node] -= div * pStressVol;
            }
          }
        }

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gStress_sum[node] += (gStress[m][node]);
          gInternalForce_sum[node] += gInternalForce[m][node];
        }
      } // if matl isn't rigid
    }   // Loop over matls

    // gStress will be normalized by gVolume (same for all matls)
    for (int m = 0; m < numMatls; m++) {
      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node   = *iter;
        gStress[m][node] = gStress_sum[node] / (gVolume[m][node] + 1.e-200);
      }
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node          = *iter;
          gInternalForce[m][node] = gInternalForce_sum[node];
        }
      }
    } // Loop over matls

    // Fill in the value for the all in one material
    // gStress_sum will be normalized by gVolume_global
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      gStress_sum[node] /= (gVolume_sum[node] + 1.e-200);
    }

    // save boundary forces before apply symmetry boundary condition.
    bool did_it_already = false;
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (!did_it_already && !mpm_matl->getIsRigid()) {
        did_it_already = true;
        for (auto face : d_boundaryTractionFaces) {

          // Check if the face is on an external boundary
          if (patch->getBCType(face) == Patch::Neighbor) {
            continue;
          }

          // We are on the boundary, i.e. not on an interior patch
          // boundary, and also on the correct side,
          // so do the traction accumulation . . .
          // loop nodes to find forces
          IntVector projlow, projhigh;
          patch->getFaceNodes(face, 0, projlow, projhigh);

          for (int i = projlow.x(); i < projhigh.x(); i++) {
            for (int j = projlow.y(); j < projhigh.y(); j++) {
              for (int k = projlow.z(); k < projhigh.z(); k++) {
                IntVector node(i, j, k);
                // flip sign so that pushing on boundary gives positive force
                bndyForce[face] -= gInternalForce[m][node];

                double celldepth = dx[face / 2];
                bndyArea[face] += gVolume[m][node] / celldepth;
              }
            }
          }
        } // faces
      }   // if
    }     // matls
  }

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (auto face : d_boundaryTractionFaces) {
    new_dw->put(sumvec_vartype(bndyForce[face]),
                d_mpm_labels->BndyForceLabel[face]);
    new_dw->put(sum_vartype(bndyArea[face]),
                d_mpm_labels->BndyContactAreaLabel[face]);
  }
}

void
ImpMPM::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  printSchedule(patches, cout_doing, "ImpMPM::scheduleRefine");
  Task* t = scinew Task("ImpMPM::refine", this, &ImpMPM::refine);

  t->computes(d_mpm_labels->partCountLabel);
  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_mpm_labels->pMassLabel);
  t->computes(d_mpm_labels->pVolumeLabel);
  t->computes(d_mpm_labels->pDispLabel);
  t->computes(d_mpm_labels->pVelocityLabel);
  t->computes(d_mpm_labels->pAccelerationLabel);
  t->computes(d_mpm_labels->pExternalForceLabel);
  t->computes(d_mpm_labels->pTemperatureLabel);
  t->computes(d_mpm_labels->pTempPreviousLabel);
  t->computes(d_mpm_labels->pSizeLabel);
  t->computes(d_mpm_labels->pParticleIDLabel);
  t->computes(d_mpm_labels->pDefGradLabel);
  t->computes(d_mpm_labels->pStressLabel);
  t->computes(d_mpm_labels->pCellNAPIDLabel);
  t->computes(d_mpm_labels->delTLabel, getLevel(patches));

  t->computes(d_mpm_labels->pExternalHeatFluxLabel);

  t->computes(d_mpm_labels->heatRate_CCLabel);
  if (!d_mpm_flags->d_doGridReset) {
    t->computes(d_mpm_labels->gDisplacementLabel);
  }

  if (d_mpm_flags->d_useLoadCurves) {
    // Computes the load curve ID associated with each particle
    t->computes(d_mpm_labels->pLoadCurveIDLabel);
  }

  t->computes(d_mpm_labels->NC_CCweightLabel, d_oneMaterial);

  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));

  Level* level = const_cast<Level*>(getLevel(patches));

  if (d_mpm_flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    d_heatConductionTasks->scheduleInitializeHeatFluxBCs(level, sched);
  }
}

void
ImpMPM::scheduleRefineInterface(const LevelP& /*fineLevel*/,
                                SchedulerP& /*scheduler*/,
                                bool,
                                bool)
{
  // do nothing for now
}

void
ImpMPM::scheduleCoarsen(const LevelP& /*coarseLevel*/, SchedulerP& /*sched*/)
{
  // do nothing for now
}
//______________________________________________________________________
// Schedule to mark flags for AMR regridding
void
ImpMPM::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  // main way is to count particles, but for now we only want particles on
  // the finest level.  Thus to schedule cells for regridding during the
  // execution, we'll coarsen the flagged cells (see coarsen).
  printSchedule(coarseLevel, cout_doing, "ImpMPM::scheduleErrorEstimate");

  if (cout_doing.active()) {
    cout_doing << "ImpMPM::scheduleErrorEstimate on level "
               << coarseLevel->getIndex() << '\n';
  }

  // The simulation controller should not schedule it every time step
  Task* task = scinew Task("errorEstimate", this, &ImpMPM::errorEstimate);

  // if the finest level, compute flagged cells
  if (coarseLevel->getIndex() == coarseLevel->getGrid()->numLevels() - 1) {
    task->requires(Task::NewDW, d_mpm_labels->pXLabel, Ghost::AroundCells, 0);
  } else {
    task->requires(Task::NewDW,
                   d_regridder->getRefineFlagLabel(),
                   0,
                   Task::FineLevel,
                   d_regridder->refineFlagMaterials(),
                   Task::NormalDomain,
                   Ghost::None,
                   0);
  }
  task->modifies(d_regridder->getRefineFlagLabel(),
                 d_regridder->refineFlagMaterials());
  task->modifies(d_regridder->getRefinePatchFlagLabel(),
                 d_regridder->refineFlagMaterials());
  sched->addTask(task,
                 coarseLevel->eachPatch(),
                 d_materialManager->allMaterials("MPM"));
}
//______________________________________________________________________
// Schedule to mark initial flags for AMR regridding
void
ImpMPM::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                     SchedulerP& sched)
{
  scheduleErrorEstimate(coarseLevel, sched);
}

void
ImpMPM::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  if (d_switchCriteria) {
    d_switchCriteria->scheduleSwitchTest(level, sched);
  }
}

void
ImpMPM::initialErrorEstimate(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* /*matls*/,
                             DataWarehouse*,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing initialErrorEstimate");

    CCVariable<int> refineFlag;
    PerPatch<PatchFlagP> refinePatchFlag;
    new_dw->getModifiable(refineFlag,
                          d_regridder->getRefineFlagLabel(),
                          0,
                          patch);
    new_dw->get(refinePatchFlag,
                d_regridder->getRefinePatchFlagLabel(),
                0,
                patch);

    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
      constParticleVariable<Point> pX;
      new_dw->get(pX, d_mpm_labels->pXLabel, pset);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        refineFlag[patch->getLevel()->getCellIndex(pX[*iter])] = true;
        refinePatch->set();
      }
    }
  }
}

void
ImpMPM::errorEstimate(const ProcessorGroup* group,
                      const PatchSubset* patches,
                      const MaterialSubset* matls,
                      DataWarehouse* old_dw,
                      DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  if (level->getIndex() == level->getGrid()->numLevels() - 1) {
    // on finest level, we do the same thing as initialErrorEstimate, so call it
    initialErrorEstimate(group, patches, matls, old_dw, new_dw);
  } else {
    // coarsen the errorflag.
    const Level* fineLevel = level->getFinerLevel().get_rep();

    for (int p = 0; p < patches->size(); p++) {
      const Patch* coarsePatch = patches->get(p);

      printTask(patches,
                coarsePatch,
                cout_doing,
                "Doing IMPM::errorEstimate\t\t\t\t\t");

      CCVariable<int> refineFlag;
      PerPatch<PatchFlagP> refinePatchFlag;

      new_dw->getModifiable(refineFlag,
                            d_regridder->getRefineFlagLabel(),
                            0,
                            coarsePatch);
      new_dw->get(refinePatchFlag,
                  d_regridder->getRefinePatchFlagLabel(),
                  0,
                  coarsePatch);

      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);

      // coarsen the fineLevel flag
      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        IntVector cl, ch, fl, fh;
        getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);
        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }
        constCCVariable<int> fineErrorFlag;
        new_dw->getRegion(fineErrorFlag,
                          d_regridder->getRefineFlagLabel(),
                          0,
                          fineLevel,
                          fl,
                          fh,
                          false);

        //__________________________________
        // if the fine level flag has been set
        // then set the corrsponding coarse level flag
        for (CellIterator iter(cl, ch); !iter.done(); iter++) {
          IntVector fineStart(level->mapCellToFiner(*iter));

          for (CellIterator inside(IntVector(0, 0, 0),
                                   fineLevel->getRefinementRatio());
               !inside.done();
               inside++) {

            if (fineErrorFlag[fineStart + *inside]) {
              refineFlag[*iter] = 1;
              refinePatch->set();
            }
          }
        } // coarse patch iterator
      }   // fine patch loop
    }     // coarse patch loop
  }
}

void
ImpMPM::refine(const ProcessorGroup*,
               const PatchSubset* patches,
               const MaterialSubset* /*matls*/,
               DataWarehouse*,
               DataWarehouse* new_dw)
{
  // just create a particle subset if one doesn't exist
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing refine");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      if (cout_doing.active()) {
        cout_doing << "Doing refine on patch " << patch->getID()
                   << " material # = " << dwi << std::endl;
      }

      // this is a new patch, so create empty particle variables.
      if (!new_dw->haveParticleSubset(dwi, patch)) {
        ParticleSubset* pset = new_dw->createParticleSubset(0, dwi, patch);

        // Create arrays for the particle data
        ParticleVariable<Point> pX;
        ParticleVariable<double> pMass, pVolume, pTemperature;
        ParticleVariable<Vector> pVelocity, pExternalForce, pDisp;
        ParticleVariable<Matrix3> pSize;
        ParticleVariable<double> pTempPrev;
        ParticleVariable<int> pLoadCurve;
        ParticleVariable<long64> pID;
        ParticleVariable<Matrix3> pDefGrad, pStress;

        new_dw->allocateAndPut(pX, d_mpm_labels->pXLabel, pset);
        new_dw->allocateAndPut(pMass, d_mpm_labels->pMassLabel, pset);
        new_dw->allocateAndPut(pVolume, d_mpm_labels->pVolumeLabel, pset);
        new_dw->allocateAndPut(pVelocity, d_mpm_labels->pVelocityLabel, pset);
        new_dw->allocateAndPut(pTemperature,
                               d_mpm_labels->pTemperatureLabel,
                               pset);
        new_dw->allocateAndPut(pTempPrev,
                               d_mpm_labels->pTempPreviousLabel,
                               pset);
        new_dw->allocateAndPut(pExternalForce,
                               d_mpm_labels->pExternalForceLabel,
                               pset);
        new_dw->allocateAndPut(pID, d_mpm_labels->pParticleIDLabel, pset);
        new_dw->allocateAndPut(pDisp, d_mpm_labels->pDispLabel, pset);
        if (d_mpm_flags->d_useLoadCurves) {
          new_dw->allocateAndPut(pLoadCurve,
                                 d_mpm_labels->pLoadCurveIDLabel,
                                 pset);
        }
        new_dw->allocateAndPut(pSize, d_mpm_labels->pSizeLabel, pset);

        mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                           mpm_matl,
                                                           new_dw);
      }
    }
    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight,
                           d_mpm_labels->NC_CCweightLabel,
                           0,
                           patch);
    NC_CCweight.initialize(0.125);
    for (Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
         face                 = Patch::nextFace(face)) {
      int mat_id = 0;
      if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {

        for (CellIterator iter = patch->getFaceIterator(face, Patch::FaceNodes);
             !iter.done();
             iter++) {
          NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
        }
      }
    }
  }
} // end refine()
