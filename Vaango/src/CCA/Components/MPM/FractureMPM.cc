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

#include <CCA/Components/MPM/FractureMPM.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>
#include <CCA/Components/MPM/HeatConduction/HeatConductionTasks.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
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
#include <Core/Math/MinMax.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>

#undef KUMAR
// #define KUMAR

using namespace Uintah;

static DebugStream cout_doing("MPM", false);
static DebugStream cout_dbg("FractureMPM", false);
static DebugStream cout_convert("MPMConv", false);
static DebugStream cout_heat("MPMHeat", false);
static DebugStream amr_doing("AMRMPM", false);

static Vector
face_norm(Patch::FaceType f)
{
  switch (f) {
    case Patch::xminus:
      return Vector(-1, 0, 0);
    case Patch::xplus:
      return Vector(1, 0, 0);
    case Patch::yminus:
      return Vector(0, -1, 0);
    case Patch::yplus:
      return Vector(0, 1, 0);
    case Patch::zminus:
      return Vector(0, 0, -1);
    case Patch::zplus:
      return Vector(0, 0, 1);
    default:
      return Vector(0, 0, 0); // oops !
  }
}

FractureMPM::FractureMPM(const ProcessorGroup* myworld,
                         const MaterialManagerP& mat_manager)
  : SerialMPM(myworld, mat_manager)
{
  crackModel = nullptr;
}

FractureMPM::~FractureMPM()
{
  delete crackModel;
}

void
FractureMPM::problemSetup(const ProblemSpecP& prob_spec,
                          const ProblemSpecP& restart_prob_spec,
                          GridP& grid)
{
  SerialMPM::problemSetup(prob_spec, restart_prob_spec, grid);

  // for FractureMPM
  crackModel     = scinew Crack(prob_spec,
                            d_materialManager,
                            d_output,
                            d_mpmLabels.get(),
                            d_mpmFlags.get());
}

void
FractureMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  Task* t = scinew Task("FractureMPM::actuallyInitialize",
                        this,
                        &FractureMPM::actuallyInitialize);

  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_mpmLabels->partCountLabel);
  t->computes(d_mpmLabels->pXLabel);
  t->computes(d_mpmLabels->pDispLabel);
  t->computes(d_mpmLabels->pMassLabel);
  t->computes(d_mpmLabels->pFiberDirLabel);
  t->computes(d_mpmLabels->pVolumeLabel);
  t->computes(d_mpmLabels->pTemperatureLabel);
  t->computes(d_mpmLabels->pTempPreviousLabel); // for thermal stress
  t->computes(d_mpmLabels->pdTdtLabel);
  t->computes(d_mpmLabels->pVelocityLabel);
  t->computes(d_mpmLabels->pExternalForceLabel);
  t->computes(d_mpmLabels->pParticleIDLabel);
  t->computes(d_mpmLabels->pDefGradLabel);
  t->computes(d_mpmLabels->pStressLabel);
  t->computes(d_mpmLabels->pSizeLabel);
  t->computes(d_mpmLabels->pDispGradsLabel);
  t->computes(d_mpmLabels->pStrainEnergyDensityLabel);
  t->computes(d_mpmLabels->delTLabel, level.get_rep());
  t->computes(d_mpmLabels->pCellNAPIDLabel, zeroth_matl);

  // Debugging Scalar
  if (d_mpmFlags->d_withColor) {
    t->computes(d_mpmLabels->pColorLabel);
  }

  if (d_mpmFlags->d_useLoadCurves) {
    // Computes the load curve ID associated with each particle
    t->computes(d_mpmLabels->pLoadCurveIDLabel);
  }

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    // Computes accumulated strain energy
    t->computes(d_mpmLabels->AccStrainEnergyLabel);
  }

  // artificial damping coeff initialized to 0.0
  if (cout_dbg.active()) {
    cout_doing << "Artificial Damping Coeff = "
               << d_mpmFlags->d_artificialDampCoeff
               << " 8 or 27 = " << d_mpmFlags->d_8or27 << std::endl;
  }

  int numMPM              = d_materialManager->getNumMaterials("MPM");
  const PatchSet* patches = level->eachPatch();
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));

  schedulePrintParticleCount(level, sched);

  // for FractureMPM: Descritize crack plane into triangular elements
  t = scinew Task("Crack:CrackDiscretization",
                  crackModel,
                  &Crack::CrackDiscretization);
  crackModel->addComputesAndRequiresCrackDiscretization(
    t,
    level->eachPatch(),
    d_materialManager->allMaterials("MPM"));
  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) {
    delete zeroth_matl; // shouln't happen, but...
  }

  if (d_mpmFlags->d_useLoadCurves) {
    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
  }
}

void
FractureMPM::scheduleInitializeAddedMaterial(const LevelP& level,
                                             SchedulerP& sched)
{
  if (cout_doing.active()) {
    cout_doing << "Doing FractureMPM::scheduleInitializeAddedMaterial "
               << std::endl;
  }

  Task* t = scinew Task("FractureMPM::actuallyInitializeAddedMaterial",
                        this,
                        &FractureMPM::actuallyInitializeAddedMaterial);

  int numALLMatls          = d_materialManager->getNumMaterials();
  int numMPMMatls          = d_materialManager->getNumMaterials("MPM");
  MaterialSubset* add_matl = scinew MaterialSubset();
  std::cout << "Added Material = " << numALLMatls - 1 << std::endl;
  add_matl->add(numALLMatls - 1);
  add_matl->addReference();

  t->computes(d_mpmLabels->partCountLabel, add_matl);
  t->computes(d_mpmLabels->pXLabel, add_matl);
  t->computes(d_mpmLabels->pDispLabel, add_matl);
  t->computes(d_mpmLabels->pMassLabel, add_matl);
  t->computes(d_mpmLabels->pVolumeLabel, add_matl);
  t->computes(d_mpmLabels->pTemperatureLabel, add_matl);
  t->computes(d_mpmLabels->pTempPreviousLabel, add_matl); // for thermal stress
  t->computes(d_mpmLabels->pdTdtLabel, add_matl);
  t->computes(d_mpmLabels->pVelocityLabel, add_matl);
  t->computes(d_mpmLabels->pExternalForceLabel, add_matl);
  t->computes(d_mpmLabels->pParticleIDLabel, add_matl);
  t->computes(d_mpmLabels->pDefGradLabel, add_matl);
  t->computes(d_mpmLabels->pStressLabel, add_matl);
  t->computes(d_mpmLabels->pSizeLabel, add_matl);

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    // Computes accumulated strain energy
    t->computes(d_mpmLabels->AccStrainEnergyLabel);
  }

  const PatchSet* patches = level->eachPatch();

  MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
    d_materialManager->getMaterial("MPM", numMPMMatls - 1));
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  cm->addInitialComputesAndRequires(t, mpm_matl, patches);

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));

  // The task will have a reference to add_matl
  if (add_matl->removeReference()) {
    delete add_matl; // shouln't happen, but...
  }
}

void
FractureMPM::scheduleInitializePressureBCs(const LevelP& level,
                                           SchedulerP& sched)
{
  MaterialSubset* loadCurveIndex = scinew MaterialSubset();
  int nofPressureBCs             = 0;
  for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
       ii++) {
    string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
    if (bcs_type == "Pressure") {
      loadCurveIndex->add(nofPressureBCs++);
    }
  }
  if (nofPressureBCs > 0) {

    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task("FractureMPM::countMaterialPointsPerLoadCurve",
                          this,
                          &FractureMPM::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpmLabels->materialPointsPerLoadCurveLabel,
                loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t,
                   level->eachPatch(),
                   d_materialManager->allMaterials("MPM"));

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("FractureMPM::initializePressureBC",
                    this,
                    &FractureMPM::initializePressureBC);
    t->requires(Task::OldDW, d_mpmLabels->simulationTimeLabel);
    t->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpmLabels->materialPointsPerLoadCurveLabel,
                loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpmLabels->pExternalForceLabel);
    sched->addTask(t,
                   level->eachPatch(),
                   d_materialManager->allMaterials("MPM"));
  }
}

void
FractureMPM::scheduleComputeStableTimestep(const LevelP&, SchedulerP&)
{
  // Nothing to do here - delt is computed as a by-product of the
  // consitutive model
}

void
FractureMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  const PatchSet* patches  = level->eachPatch();
  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleParticleVelocityField(sched, patches, matls); // for FractureMPM
  scheduleInterpolateParticlesToGrid(sched, patches, matls);
  scheduleAdjustCrackContactInterpolated(sched,
                                         patches,
                                         matls); // for FractureMPM
  scheduleMomentumExchangeInterpolated(sched, patches, matls);
  scheduleComputeContactArea(sched, patches, matls);
  scheduleComputeInternalForce(sched, patches, matls);
  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  scheduleAdjustCrackContactIntegrated(sched, patches, matls); // for
                                                               // FractureMPM
  scheduleMomentumExchangeIntegrated(sched, patches, matls);
  scheduleSetGridBoundaryConditions(sched, patches, matls);
  scheduleComputeStressTensor(sched, patches, matls);

  d_heatConductionTasks->scheduleCompute(sched, patches, matls);

  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  scheduleCalculateFractureParameters(sched, patches, matls); // for FractureMPM
  scheduleDoCrackPropagation(sched, patches, matls);          // for FractureMPM
  scheduleMoveCracks(sched, patches, matls);                  // for FractureMPM
  scheduleUpdateCrackFront(sched, patches, matls);            // for FractureMPM

  sched->scheduleParticleRelocation(level,
                                    d_mpmLabels->pXLabel_preReloc,
                                    d_particleState_preReloc,
                                    d_mpmLabels->pXLabel,
                                    d_particleState,
                                    d_mpmLabels->pParticleIDLabel,
                                    matls);
}

void
FractureMPM::scheduleApplyExternalLoads(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  /*
   * applyExternalLoads
   *   in(p.externalForce, p.externalheatrate)
   *   out(p.externalForceNew, p.externalheatrateNew) */
  Task* t = scinew Task("FractureMPM::applyExternalLoads",
                        this,
                        &FractureMPM::applyExternalLoads);

  t->requires(Task::OldDW, d_mpmLabels->simulationTimeLabel);
  t->requires(Task::OldDW, d_mpmLabels->pExternalForceLabel, Ghost::None);
  t->computes(d_mpmLabels->pExtForceLabel_preReloc);
  if (d_mpmFlags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
    t->requires(Task::OldDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpmLabels->pLoadCurveIDLabel_preReloc);
  }

  //  t->computes(Task::OldDW, d_mpmLabels->pExternalHeatRateLabel_preReloc);

  sched->addTask(t, patches, matls);
}

// Determine velocity field for each particle-node pair
void
FractureMPM::scheduleParticleVelocityField(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::ParticleVelocityField",
                        crackModel,
                        &Crack::ParticleVelocityField);

  crackModel->addComputesAndRequiresParticleVelocityField(t, patches, matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  /* interpolateParticlesToGrid
   *   in(P.MASS, P.VELOCITY, P.NAT_X)
   *   operation(interpolate the P.MASS and P.VEL to the grid
   *             using P.NAT_X and some shape function evaluations)
   *   out(G.MASS, G.VELOCITY) */

  Task* t = scinew Task("FractureMPM::interpolateParticlesToGrid",
                        this,
                        &FractureMPM::interpolateParticlesToGrid);

  t->requires(Task::OldDW,
              d_mpmLabels->pMassLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pVolumeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pVelocityLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pTemperatureLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pExtForceLabel_preReloc,
              Ghost::AroundNodes,
              d_numGhostParticles);
  // t->requires(Task::OldDW, d_mpmLabels->pExternalHeatRateLabel,
  // Ghost::AroundNodes,d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);

  t->computes(d_mpmLabels->gMassLabel);
  t->computes(d_mpmLabels->gMassLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpmLabels->gTemperatureLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpmLabels->gVolumeLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpmLabels->gVelocityLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpmLabels->gSp_volLabel);
  t->computes(d_mpmLabels->gVolumeLabel);
  t->computes(d_mpmLabels->gVelocityLabel);
  t->computes(d_mpmLabels->gExternalForceLabel);
  t->computes(d_mpmLabels->gTemperatureLabel);
  t->computes(d_mpmLabels->gTemperatureNoBCLabel);
  t->computes(d_mpmLabels->gTemperatureRateLabel);
  t->computes(d_mpmLabels->gExternalHeatRateLabel);
  t->computes(d_mpmLabels->gNumNearParticlesLabel);
  t->computes(d_mpmLabels->TotalMassLabel);

  // for FractureMPM
  t->requires(Task::OldDW,
              d_mpmLabels->pDispLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->computes(d_mpmLabels->GMassLabel);
  t->computes(d_mpmLabels->GSp_volLabel);
  t->computes(d_mpmLabels->GVolumeLabel);
  t->computes(d_mpmLabels->GVelocityLabel);
  t->computes(d_mpmLabels->GExternalForceLabel);
  t->computes(d_mpmLabels->GTemperatureLabel);
  t->computes(d_mpmLabels->GTemperatureNoBCLabel);
  t->computes(d_mpmLabels->GExternalHeatRateLabel);
  t->computes(d_mpmLabels->gDisplacementLabel);
  t->computes(d_mpmLabels->GDisplacementLabel);

  sched->addTask(t, patches, matls);
}

// Check crack contact and adjust grid velocity field
void
FractureMPM::scheduleAdjustCrackContactInterpolated(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::AdjustCrackContactInterpolated",
                        crackModel,
                        &Crack::AdjustCrackContactInterpolated);

  crackModel->addComputesAndRequiresAdjustCrackContactInterpolated(t,
                                                                   patches,
                                                                   matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleMomentumExchangeInterpolated(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  contactModel->addComputesAndRequires(sched,
                                       patches,
                                       matls,
                                       d_mpmLabels->gVelocityLabel);
}

/////////////////////////////////////////////////////////////////////////
/*!  **WARNING** In addition to the stresses and deformations, the internal
 *               heat rate in the particles (pdTdtLabel)
 *               is computed here */
/////////////////////////////////////////////////////////////////////////
void
FractureMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  // for thermal stress analysis
  scheduleComputeParticleTempFromGrid(sched, patches, matls);

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("FractureMPM::computeStressTensor",
                        this,
                        &FractureMPM::computeStressTensor);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    t->computes(d_mpmLabels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_mpmLabels->delTLabel, getLevel(patches));
  t->computes(d_mpmLabels->StrainEnergyLabel);

  sched->addTask(t, patches, matls);

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }
  if (d_mpmFlags->d_artificialViscosity) {
    scheduleComputeArtificialViscosity(sched, patches, matls);
  }
}

// Compute particle temperature by interpolating grid temperature
// for thermal stress analysis
void
FractureMPM::scheduleComputeParticleTempFromGrid(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  Task* t = scinew Task("FractureMPM::computeParticleTempFromGrid",
                        this,
                        &FractureMPM::computeParticleTempFromGrid);
  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->computes(d_mpmLabels->pTempCurrentLabel);
  sched->addTask(t, patches, matls);
}

// Compute the accumulated strain energy
void
FractureMPM::scheduleComputeAccStrainEnergy(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls)
{
  Task* t = scinew Task("FractureMPM::computeAccStrainEnergy",
                        this,
                        &FractureMPM::computeAccStrainEnergy);
  t->requires(Task::OldDW, d_mpmLabels->AccStrainEnergyLabel);
  t->requires(Task::NewDW, d_mpmLabels->StrainEnergyLabel);
  t->computes(d_mpmLabels->AccStrainEnergyLabel);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleComputeArtificialViscosity(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  Task* t = scinew Task("FractureMPM::computeArtificialViscosity",
                        this,
                        &FractureMPM::computeArtificialViscosity);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->pVolumeLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes); // for FractureMPM
  t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, Ghost::None);
  t->computes(d_mpmLabels->p_qLabel);

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleComputeContactArea(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  /*
   * computeContactArea */
  if (d_boundaryTractionFaces.size() > 0) {
    Task* t = scinew Task("FractureMPM::computeContactArea",
                          this,
                          &FractureMPM::computeContactArea);

    t->requires(Task::NewDW, d_mpmLabels->gVolumeLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpmLabels->GVolumeLabel,
                Ghost::None); // for FractureMPM
    for (auto& face : d_boundaryTractionFaces) {
      int iface = (int)(face);
      t->computes(d_mpmLabels->BndyContactCellAreaLabel[iface]);
    }
    sched->addTask(t, patches, matls);
  }
}

void
FractureMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  /*
   * computeInternalForce
   *   in(P.CONMOD, P.NAT_X, P.VOLUME)
   *   operation(evaluate the divergence of the stress (stored in
   *   P.CONMOD) using P.NAT_X and the gradients of the
   *   shape functions)
   * out(G.F_INTERNAL) */

  Task* t = scinew Task("FractureMPM::computeInternalForce",
                        this,
                        &FractureMPM::computeInternalForce);

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW,
              d_mpmLabels->gMassLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain,
              Ghost::None);
  t->requires(Task::OldDW,
              d_mpmLabels->pStressLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pVolumeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pMassLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);

  // for FractureMPM
  t->requires(Task::NewDW,
              d_mpmLabels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW, d_mpmLabels->GMassLabel, Ghost::None);
  t->computes(d_mpmLabels->GInternalForceLabel);
  t->computes(d_mpmLabels->TotalVolumeDeformedLabel);

  if (d_mpmFlags->d_withICE) {
    t->requires(Task::NewDW,
                d_mpmLabels->pPressureLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }

  if (d_mpmFlags->d_artificialViscosity) {
    t->requires(Task::NewDW,
                d_mpmLabels->p_qLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }

  t->computes(d_mpmLabels->gInternalForceLabel);

  for (std::list<Patch::FaceType>::const_iterator ftit(
         d_boundaryTractionFaces.begin());
       ftit != d_boundaryTractionFaces.end();
       ftit++) {
    int iface = (int)(*ftit);
    t->requires(Task::NewDW, d_mpmLabels->BndyContactCellAreaLabel[iface]);
    t->computes(d_mpmLabels->BndyForceLabel[iface]);
    t->computes(d_mpmLabels->BndyContactAreaLabel[iface]);
    t->computes(d_mpmLabels->BndyTractionLabel[iface]);
  }

  t->computes(d_mpmLabels->gStressForSavingLabel);
  t->computes(d_mpmLabels->gStressForSavingLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleComputeAndIntegrateAcceleration");

  Task* t = scinew Task("MPM::computeAndIntegrateAcceleration",
                        this,
                        &FractureMPM::computeAndIntegrateAcceleration);

  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gExternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  t->requires(Task::NewDW, d_mpmLabels->GMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->GVelocityLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->GInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->GExternalForceLabel, Ghost::None);

  t->computes(d_mpmLabels->gVelocityStarLabel);
  t->computes(d_mpmLabels->gAccelerationLabel);

  t->computes(d_mpmLabels->GAccelerationLabel);
  t->computes(d_mpmLabels->GVelocityStarLabel);

  sched->addTask(t, patches, matls);
}

// Check crack contact and adjust nodal velocities and accelerations
void
FractureMPM::scheduleAdjustCrackContactIntegrated(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::AdjustCrackContactIntegrated",
                        crackModel,
                        &Crack::AdjustCrackContactIntegrated);

  crackModel->addComputesAndRequiresAdjustCrackContactIntegrated(t,
                                                                 patches,
                                                                 matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleMomentumExchangeIntegrated(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  contactModel->addComputesAndRequires(sched,
                                       patches,
                                       matls,
                                       d_mpmLabels->gVelocityStarLabel);
}

void
FractureMPM::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)

{
  Task* t = scinew Task("FractureMPM::setGridBoundaryConditions",
                        this,
                        &FractureMPM::setGridBoundaryConditions);

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  t->modifies(d_mpmLabels->gAccelerationLabel, mss);
  t->modifies(d_mpmLabels->gVelocityStarLabel, mss);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  // for FractureMPM
  t->modifies(d_mpmLabels->GAccelerationLabel, mss);
  t->modifies(d_mpmLabels->GVelocityStarLabel, mss);
  t->requires(Task::NewDW, d_mpmLabels->GVelocityLabel, Ghost::None);

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls)

{
  /*
   * interpolateToParticlesAndUpdate
   *   in(G.ACCELERATION, G.VELOCITY_STAR, P.NAT_X)
   *   operation(interpolate acceleration and v* to particles and
   *   integrate these to get new particle velocity and position)
   * out(P.VELOCITY, P.X, P.NAT_X) */

  Task* t = scinew Task("FractureMPM::interpolateToParticlesAndUpdate",
                        this,
                        &FractureMPM::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  t->requires(Task::NewDW,
              d_mpmLabels->gAccelerationLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureRateLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureNoBCLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->frictionalWorkLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pTemperatureLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pVelocityLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pDispLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t->modifies(d_mpmLabels->pVolumeLabel_preReloc);
  // for thermal stress analysis
  t->requires(Task::NewDW, d_mpmLabels->pTempCurrentLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel_preReloc, Ghost::None);

  // for FractureMPM
  t->requires(Task::NewDW,
              d_mpmLabels->GAccelerationLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GTemperatureRateLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->GTemperatureNoBCLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->requires(Task::NewDW, d_mpmLabels->pgCodeLabel, Ghost::None);

  if (d_mpmFlags->d_withICE) {
    t->requires(Task::NewDW,
                d_mpmLabels->dTdt_NCLabel,
                Ghost::AroundCells,
                d_numGhostNodes);
    t->requires(Task::NewDW,
                d_mpmLabels->massBurnFractionLabel,
                Ghost::AroundCells,
                d_numGhostNodes);
  }

  t->computes(d_mpmLabels->pDispLabel_preReloc);
  t->computes(d_mpmLabels->pVelocityLabel_preReloc);
  t->computes(d_mpmLabels->pXLabel_preReloc);
  t->computes(d_mpmLabels->pParticleIDLabel_preReloc);
  t->computes(d_mpmLabels->pTemperatureLabel_preReloc);
  t->computes(d_mpmLabels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpmLabels->pMassLabel_preReloc);
  t->computes(d_mpmLabels->pSizeLabel_preReloc);
  t->computes(d_mpmLabels->pXXLabel);
  t->computes(d_mpmLabels->pKineticEnergyDensityLabel); // for FractureMPM

  t->computes(d_mpmLabels->KineticEnergyLabel);
  t->computes(d_mpmLabels->ThermalEnergyLabel);
  t->computes(d_mpmLabels->CenterOfMassPositionLabel);
  t->computes(d_mpmLabels->TotalMomentumLabel);

  // debugging scalar
  if (d_mpmFlags->d_withColor) {
    t->requires(Task::OldDW, d_mpmLabels->pColorLabel, Ghost::None);
    t->computes(d_mpmLabels->pColorLabel_preReloc);
  }

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleCalculateFractureParameters(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  // Get nodal solutions
  Task* t = scinew Task("Crack::GetNodalSolutions",
                        crackModel,
                        &Crack::GetNodalSolutions);
  crackModel->addComputesAndRequiresGetNodalSolutions(t, patches, matls);
  sched->addTask(t, patches, matls);

  // cfnset & cfsset
  t = scinew Task("Crack::CrackFrontNodeSubset",
                  crackModel,
                  &Crack::CrackFrontNodeSubset);
  crackModel->addComputesAndRequiresCrackFrontNodeSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Compute fracture parameters (J, K,...)
  t = scinew Task("Crack::CalculateFractureParameters",
                  crackModel,
                  &Crack::CalculateFractureParameters);
  crackModel->addComputesAndRequiresCalculateFractureParameters(t,
                                                                patches,
                                                                matls);
  sched->addTask(t, patches, matls);
}

// Do crack propgation
void
FractureMPM::scheduleDoCrackPropagation(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  // Propagate crack-front points
  Task* t = scinew Task("Crack::PropagateCrackFrontPoints",
                        crackModel,
                        &Crack::PropagateCrackFrontPoints);
  crackModel->addComputesAndRequiresPropagateCrackFrontPoints(t,
                                                              patches,
                                                              matls);
  sched->addTask(t, patches, matls);

  // Construct the new crack-front elems and new crack-front segments.
  // The new crack-front is temporary, and will be updated after moving cracks
  t = scinew Task("Crack::ConstructNewCrackFrontElems",
                  crackModel,
                  &Crack::ConstructNewCrackFrontElems);
  crackModel->addComputesAndRequiresConstructNewCrackFrontElems(t,
                                                                patches,
                                                                matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleMoveCracks(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls)
{
  // Set up cpset -- crack node subset in each patch
  Task* t = scinew Task("Crack::CrackPointSubset",
                        crackModel,
                        &Crack::CrackPointSubset);
  crackModel->addComputesAndRequiresCrackPointSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Move crack points
  t = scinew Task("Crack::MoveCracks", crackModel, &Crack::MoveCracks);
  crackModel->addComputesAndRequiresMoveCracks(t, patches, matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleUpdateCrackFront(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  // Set up cfnset & cfsset -- subset for the temporary crack-front
  // and crack-front segment subset for each patch
  Task* t = scinew Task("Crack::CrackFrontNodeSubset",
                        crackModel,
                        &Crack::CrackFrontNodeSubset);
  crackModel->addComputesAndRequiresCrackFrontNodeSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Recollect crack-front segments, discarding the dead segments,
  // calculating normals, indexes and so on
  t = scinew Task("Crack::RecollectCrackFrontSegments",
                  crackModel,
                  &Crack::RecollectCrackFrontSegments);
  crackModel->addComputesAndRequiresRecollectCrackFrontSegments(t,
                                                                patches,
                                                                matls);
  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  Task* task = scinew Task("FractureMPM::refine", this, &FractureMPM::refine);
  sched->addTask(task, patches, d_materialManager->allMaterials("MPM"));
  // do nothing for now
}

void
FractureMPM::scheduleRefineInterface(const LevelP& /*fineLevel*/,
                                     SchedulerP& /*scheduler*/,
                                     bool,
                                     bool)
{
  // do nothing for now
}

void
FractureMPM::scheduleCoarsen(const LevelP& /*coarseLevel*/,
                             SchedulerP& /*sched*/)
{
  // do nothing for now
}

/// Schedule to mark flags for AMR regridding
void
FractureMPM::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  // main way is to count particles, but for now we only want particles on
  // the finest level.  Thus to schedule cells for regridding during the
  // execution, we'll coarsen the flagged cells (see coarsen).

  if (cout_doing.active()) {
    cout_doing << "FractureMPM::scheduleErrorEstimate on level "
               << coarseLevel->getIndex() << '\n';
  }

  // Estimate error - this should probably be in it's own schedule,
  // and the simulation controller should not schedule it every time step
  Task* task = scinew Task("errorEstimate", this, &FractureMPM::errorEstimate);

  // if the finest level, compute flagged cells
  if (coarseLevel->getIndex() == coarseLevel->getGrid()->numLevels() - 1) {
    task->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::AroundCells, 0);
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

/// Schedule to mark initial flags for AMR regridding
void
FractureMPM::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                          SchedulerP& sched)
{

  if (cout_doing.active()) {
    cout_doing << "FractureMPM::scheduleErrorEstimate on level "
               << coarseLevel->getIndex() << '\n';
  }

  // Estimate error - this should probably be in it's own schedule,
  // and the simulation controller should not schedule it every time step
  Task* task =
    scinew Task("errorEstimate", this, &FractureMPM::initialErrorEstimate);
  task->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::AroundCells, 0);

  task->modifies(d_regridder->getRefineFlagLabel(),
                 d_regridder->refineFlagMaterials());
  task->modifies(d_regridder->getRefinePatchFlagLabel(),
                 d_regridder->refineFlagMaterials());
  sched->addTask(task,
                 coarseLevel->eachPatch(),
                 d_materialManager->allMaterials("MPM"));
}

void
FractureMPM::computeAccStrainEnergy(const ProcessorGroup*,
                                    const PatchSubset*,
                                    const MaterialSubset*,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  // Get the totalStrainEnergy from the old datawarehouse
  max_vartype accStrainEnergy;
  old_dw->get(accStrainEnergy, d_mpmLabels->AccStrainEnergyLabel);

  // Get the incremental strain energy from the new datawarehouse
  sum_vartype incStrainEnergy;
  new_dw->get(incStrainEnergy, d_mpmLabels->StrainEnergyLabel);

  // Add the two a put into new dw
  double totalStrainEnergy = (double)accStrainEnergy + (double)incStrainEnergy;
  new_dw->put(max_vartype(totalStrainEnergy),
              d_mpmLabels->AccStrainEnergyLabel);
}

// Calculate the number of material points per load curve
void
FractureMPM::countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset*,
                                             DataWarehouse*,
                                             DataWarehouse* new_dw)
{
  // Find the number of pressure BCs in the problem
  int nofPressureBCs = 0;
  for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
       ii++) {
    string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
    if (bcs_type == "Pressure") {
      nofPressureBCs++;

      // Loop through the patches and count
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        int numPts         = 0;
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(static_cast<MPMMaterial*>(
              d_materialManager->getMaterial("MPM", m)));
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);

          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            if (pLoadCurveID[idx] == (nofPressureBCs)) {
              ++numPts;
            }
          }
        } // matl loop
        new_dw->put(sumlong_vartype(numPts),
                    d_mpmLabels->materialPointsPerLoadCurveLabel,
                    0,
                    nofPressureBCs - 1);
      } // patch loop
    }
  }
}

// Calculate the number of material points per load curve
void
FractureMPM::initializePressureBC(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset*,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  // Get the current time
  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, d_mpmLabels->simulationTimeLabel);
  double time = simTimeVar;

  if (cout_dbg.active()) {
    cout_dbg << "Current Time (Initialize Pressure BC) = " << time << std::endl;
  }

  // Calculate the force vector at each particle
  int nofPressureBCs = 0;
  for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
       ii++) {
    string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
    if (bcs_type == "Pressure") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart,
                  d_mpmLabels->materialPointsPerLoadCurveLabel,
                  0,
                  nofPressureBCs++);

      // Save the material points per load curve in the PressureBC object
      PressureBC* pbc = dynamic_cast<PressureBC*>(
        MPMPhysicalBCFactory::mpmPhysicalBCs[ii].get());
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
            static_cast<MPMMaterial*>(static_cast<MPMMaterial*>(
              d_materialManager->getMaterial("MPM", m)));
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<Point> px;
          constParticleVariable<int> pLoadCurveID;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(px, d_mpmLabels->pXLabel, pset);
          new_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);
          new_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);

          constParticleVariable<Vector> pDisp;
          new_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);
          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(pExternalForce,
                                d_mpmLabels->pExternalForceLabel,
                                pset);

          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            if (pLoadCurveID[idx] == nofPressureBCs) {
              pExternalForce[idx] = pbc->getForceVector(px[idx],
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
FractureMPM::actuallyInitialize(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* matls,
                                DataWarehouse*,
                                DataWarehouse* new_dw)
{
  particleIndex totalParticles = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing actuallyInitialize on patch " << patch->getID()
                 << "\t\t\t MPM" << std::endl;
    }

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    for (int m = 0; m < matls->size(); m++) {
      // cerrLock.lock();
      // NOT_FINISHED("not quite right - mapping of matls, use matls->get()");
      // cerrLock.unlock();
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                         mpm_matl,
                                                         new_dw);
      // scalar used for debugging
      if (d_mpmFlags->d_withColor) {
        ParticleVariable<double> pcolor;
        int index            = mpm_matl->getDWIndex();
        ParticleSubset* pset = new_dw->getParticleSubset(index, patch);
        setParticleDefault<double>(pcolor,
                                   d_mpmLabels->pColorLabel,
                                   pset,
                                   new_dw,
                                   0.0);
      }
    }
  }

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    // Initialize the accumulated strain energy
    new_dw->put(max_vartype(0.0), d_mpmLabels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_mpmLabels->partCountLabel);
}

void
FractureMPM::actuallyInitializeAddedMaterial(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset* matls,
                                             DataWarehouse*,
                                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing actuallyInitializeAddedMaterial on patch "
                 << patch->getID() << "\t\t\t MPM" << std::endl;
    }

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    std::cout << "num MPM Matls = " << numMPMMatls << std::endl;
    CCVariable<short int> cellNAPID;
    int m                 = numMPMMatls - 1;
    MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
    new_dw->unfinalize();
    mpm_matl->createParticles(cellNAPID, patch, new_dw);

    mpm_matl->getConstitutiveModel()->initializeCMData(patch, mpm_matl, new_dw);
    new_dw->refinalize();
  }
}

void
FractureMPM::actuallyComputeStableTimestep(const ProcessorGroup*,
                                           const PatchSubset*,
                                           const MaterialSubset*,
                                           DataWarehouse*,
                                           DataWarehouse*)
{
}

void
FractureMPM::interpolateParticlesToGrid(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing interpolateParticlesToGrid on patch "
                 << patch->getID() << "\t\t MPM" << std::endl;
    }

    int numMatls      = d_materialManager->getNumMaterials("MPM");
    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    NCVariable<double> gMassglobal, gtempglobal, gVolumeglobal;
    NCVariable<Vector> gvelglobal;
    new_dw->allocateAndPut(gMassglobal,
                           d_mpmLabels->gMassLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gtempglobal,
                           d_mpmLabels->gTemperatureLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gVolumeglobal,
                           d_mpmLabels->gVolumeLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gvelglobal,
                           d_mpmLabels->gVelocityLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    gMassglobal.initialize(d_SMALL_NUM_MPM);
    gVolumeglobal.initialize(d_SMALL_NUM_MPM);
    gtempglobal.initialize(0.0);
    gvelglobal.initialize(Vector(0.0));

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Point> px;
      constParticleVariable<double> pmass, pvolume, pTemperature;
      constParticleVariable<Vector> pvelocity, pexternalforce, pdisp;
      constParticleVariable<Matrix3> psize;
      constParticleVariable<Matrix3> pDeformationMeasure;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(px, d_mpmLabels->pXLabel, pset);
      old_dw->get(pmass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pvolume, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pvelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      old_dw->get(psize, d_mpmLabels->pSizeLabel, pset);
      new_dw->get(pexternalforce, d_mpmLabels->pExtForceLabel_preReloc, pset);
      old_dw->get(pDeformationMeasure, d_mpmLabels->pDefGradLabel, pset);

      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpmLabels->pgCodeLabel, pset);
      old_dw->get(pdisp, d_mpmLabels->pDispLabel, pset);

      // Create arrays for the grid data
      NCVariable<double> gMass;
      NCVariable<double> gVolume;
      NCVariable<Vector> gVelocity;
      NCVariable<Vector> gexternalforce;
      NCVariable<double> gexternalheatrate;
      NCVariable<double> gTemperature;
      NCVariable<double> gSp_vol;
      NCVariable<double> gTemperatureNoBC;
      NCVariable<double> gTemperatureRate;
      NCVariable<double> gnumnearparticles;

      new_dw->allocateAndPut(gMass, d_mpmLabels->gMassLabel, dwi, patch);
      new_dw->allocateAndPut(gSp_vol, d_mpmLabels->gSp_volLabel, dwi, patch);
      new_dw->allocateAndPut(gVolume, d_mpmLabels->gVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(gVelocity,
                             d_mpmLabels->gVelocityLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperature,
                             d_mpmLabels->gTemperatureLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperatureNoBC,
                             d_mpmLabels->gTemperatureNoBCLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperatureRate,
                             d_mpmLabels->gTemperatureRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gexternalforce,
                             d_mpmLabels->gExternalForceLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gexternalheatrate,
                             d_mpmLabels->gExternalHeatRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gnumnearparticles,
                             d_mpmLabels->gNumNearParticlesLabel,
                             dwi,
                             patch);

      gMass.initialize(d_SMALL_NUM_MPM);
      gVolume.initialize(d_SMALL_NUM_MPM);
      gVelocity.initialize(Vector(0, 0, 0));
      gexternalforce.initialize(Vector(0, 0, 0));
      gTemperature.initialize(0);
      gTemperatureNoBC.initialize(0);
      gTemperatureRate.initialize(0);
      gexternalheatrate.initialize(0);
      gnumnearparticles.initialize(0.);
      gSp_vol.initialize(0.);

      // for FractureMPM
      NCVariable<double> Gmass;
      NCVariable<double> Gvolume;
      NCVariable<Vector> Gvelocity;
      NCVariable<Vector> Gexternalforce;
      NCVariable<double> Gexternalheatrate;
      NCVariable<double> GTemperature;
      NCVariable<double> GSp_vol;
      NCVariable<double> GTemperatureNoBC;
      NCVariable<Vector> gdisplacement;
      NCVariable<Vector> Gdisplacement;

      new_dw->allocateAndPut(Gmass, d_mpmLabels->GMassLabel, dwi, patch);
      new_dw->allocateAndPut(GSp_vol, d_mpmLabels->GSp_volLabel, dwi, patch);
      new_dw->allocateAndPut(Gvolume, d_mpmLabels->GVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(Gvelocity,
                             d_mpmLabels->GVelocityLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(GTemperature,
                             d_mpmLabels->GTemperatureLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(GTemperatureNoBC,
                             d_mpmLabels->GTemperatureNoBCLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gexternalforce,
                             d_mpmLabels->GExternalForceLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gexternalheatrate,
                             d_mpmLabels->GExternalHeatRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gdisplacement,
                             d_mpmLabels->gDisplacementLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gdisplacement,
                             d_mpmLabels->GDisplacementLabel,
                             dwi,
                             patch);

      // initialization
      Gmass.initialize(d_SMALL_NUM_MPM);
      Gvolume.initialize(d_SMALL_NUM_MPM);
      Gvelocity.initialize(Vector(0, 0, 0));
      Gexternalforce.initialize(Vector(0, 0, 0));
      GTemperature.initialize(0);
      GTemperatureNoBC.initialize(0);
      Gexternalheatrate.initialize(0);
      GSp_vol.initialize(0.);
      gdisplacement.initialize(Vector(0, 0, 0));
      Gdisplacement.initialize(Vector(0, 0, 0));

      // Interpolate particle data to Grid data.
      // This currently consists of the particle velocity and mass
      // Need to compute the lumped global mass matrix and velocity
      // Vector from the individual mass matrix and velocity vector
      // GridMass * GridVelocity =  S^T*M_D*ParticleVelocity

      double totalmass = 0;
      Vector total_mom(0.0, 0.0, 0.0);
      Vector pmom;

      double pSp_vol = 1. / mpm_matl->getInitialDensity();
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(px[idx],
                                         ni,
                                         S,
                                         psize[idx],
                                         pDeformationMeasure[idx]);

        pmom = pvelocity[idx] * pmass[idx];
        total_mom += pvelocity[idx] * pmass[idx];

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          if (patch->containsNode(ni[k])) {
            if (pgCode[idx][k] == 1) { // above crack
              gMass[ni[k]] += pmass[idx] * S[k];
              gVelocity[ni[k]] += pmom * S[k];
              gVolume[ni[k]] += pvolume[idx] * S[k];
              gexternalforce[ni[k]] += pexternalforce[idx] * S[k];
              gTemperature[ni[k]] += pTemperature[idx] * pmass[idx] * S[k];
              // gexternalheatrate[ni[k]] += pexternalheatrate[idx]      * S[k];
              gnumnearparticles[ni[k]] += 1.0;
              gSp_vol[ni[k]] += pSp_vol * pmass[idx] * S[k];
              gdisplacement[ni[k]] += pdisp[idx] * pmass[idx] * S[k];
            } else if (pgCode[idx][k] == 2) { // below crack
              Gmass[ni[k]] += pmass[idx] * S[k];
              Gvelocity[ni[k]] += pmom * S[k];
              Gvolume[ni[k]] += pvolume[idx] * S[k];
              Gexternalforce[ni[k]] += pexternalforce[idx] * S[k];
              GTemperature[ni[k]] += pTemperature[idx] * pmass[idx] * S[k];
              // Gexternalheatrate[ni[k]] += pexternalheatrate[idx]      * S[k];
              GSp_vol[ni[k]] += pSp_vol * pmass[idx] * S[k];
              Gdisplacement[ni[k]] += pdisp[idx] * pmass[idx] * S[k];
            }
          }
        } // End of loop over k
      }   // End of loop over iter

      string interp_type = d_mpmFlags->d_interpolatorType;
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        totalmass += (gMass[c] + Gmass[c]);
        gMassglobal[c] += (gMass[c] + Gmass[c]);
        gVolumeglobal[c] += (gVolume[c] + Gvolume[c]);
        gvelglobal[c] += (gVelocity[c] + Gvelocity[c]);
        gtempglobal[c] += (gTemperature[c] + GTemperature[c]);

        // above crack
        gVelocity[c] /= gMass[c];
        gTemperature[c] /= gMass[c];
        gSp_vol[c] /= gMass[c];
        gdisplacement[c] /= gMass[c];
        gTemperatureNoBC[c] = gTemperature[c];

        // below crack
        Gvelocity[c] /= Gmass[c];
        GTemperature[c] /= Gmass[c];
        GSp_vol[c] /= Gmass[c];
        Gdisplacement[c] /= Gmass[c];
        GTemperatureNoBC[c] = GTemperature[c];
      }

      // Apply grid boundary conditions to the velocity before storing the data
      MPMBoundCond bc;
      // above crack
      bc.setBoundaryCondition(patch, dwi, "Velocity", gVelocity, interp_type);
      bc.setBoundaryCondition(patch, dwi, "Symmetric", gVelocity, interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Temperature",
                              gTemperature,
                              interp_type);
      // below crack
      bc.setBoundaryCondition(patch, dwi, "Velocity", Gvelocity, interp_type);
      bc.setBoundaryCondition(patch, dwi, "Symmetric", Gvelocity, interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Temperature",
                              GTemperature,
                              interp_type);

      new_dw->put(sum_vartype(totalmass), d_mpmLabels->TotalMassLabel);

    } // End loop over materials

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gtempglobal[c] /= gMassglobal[c];
      gvelglobal[c] /= gMassglobal[c];
    }
  } // End loop over patches
}

void
FractureMPM::computeStressTensor(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  if (cout_doing.active()) {
    cout_doing << "Doing computeStressTensor:FractureMPM: \n";
  }

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

    if (cout_dbg.active()) {
      cout_dbg << " Patch = " << (patches->get(0))->getID();
      cout_dbg << " Mat = " << m;
    }

    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    if (cout_dbg.active()) {
      cout_dbg << " MPM_Mat = " << mpm_matl;
    }

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cout_dbg.active()) {
      cout_dbg << " CM = " << cm;
    }

    cm->setWorld(UintahParallelComponent::d_myworld);
    cm->computeStressTensor(patches, mpm_matl, old_dw, new_dw);

    if (cout_dbg.active()) {
      cout_dbg << " Exit\n";
    }
  }
}

void
FractureMPM::computeArtificialViscosity(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  double C0 = d_mpmFlags->d_artificialViscCoeff1;
  double C1 = d_mpmFlags->d_artificialViscCoeff2;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing computeArtificialViscosity on patch "
                 << patch->getID() << "\t\t MPM" << std::endl;
    }

    // The following scheme for removing ringing behind a shock comes from:
    // VonNeumann, J.; Richtmyer, R. D. (1950): A method for the numerical
    // calculation of hydrodynamic shocks. J. Appl. Phys., vol. 21, pp. 232.

    int numMatls      = d_materialManager->getNumMaterials("MPM");
    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      constNCVariable<Vector> gVelocity;
      ParticleVariable<double> p_q;
      constParticleVariable<Matrix3> psize;
      constParticleVariable<Point> px;
      constParticleVariable<double> pmass, pvol;
      constParticleVariable<Matrix3> pDeformationMeasure;

      new_dw->get(gVelocity,
                  d_mpmLabels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostNodes);
      old_dw->get(px, d_mpmLabels->pXLabel, pset);
      old_dw->get(pmass, d_mpmLabels->pMassLabel, pset);
      new_dw->get(pvol, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(psize, d_mpmLabels->pSizeLabel, pset);
      new_dw->allocateAndPut(p_q, d_mpmLabels->p_qLabel, pset);
      old_dw->get(pDeformationMeasure, d_mpmLabels->pDefGradLabel, pset);

      // for FractureMPM
      constNCVariable<Vector> Gvelocity;
      new_dw->get(Gvelocity,
                  d_mpmLabels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostNodes);
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpmLabels->pgCodeLabel, pset);

      Matrix3 velGrad;
      Vector dx      = patch->dCell();
      double oodx[3] = { 1. / dx.x(), 1. / dx.y(), 1. / dx.z() };
      double dx_ave  = (dx.x() + dx.y() + dx.z()) / 3.0;

      double K = 1. / mpm_matl->getConstitutiveModel()->getCompressibility();
      double c_dil;

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndShapeDerivatives(px[idx],
                                                  ni,
                                                  d_S,
                                                  psize[idx],
                                                  pDeformationMeasure[idx]);

        // get particle's velocity gradients
        Vector gvel(0., 0., 0.);
        velGrad.set(0.0);
        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          if (pgCode[idx][k] == 1) {
            gvel = gVelocity[ni[k]];
          }
          if (pgCode[idx][k] == 2) {
            gvel = Gvelocity[ni[k]];
          }
          for (int j = 0; j < 3; j++) {
            double d_SXoodx = d_S[k][j] * oodx[j];
            for (int i = 0; i < 3; i++) {
              velGrad(i, j) += gvel[i] * d_SXoodx;
            }
          }
        }

        Matrix3 D = (velGrad + velGrad.Transpose()) * .5;

        double DTrace = D.Trace();
        p_q[idx]      = 0.0;
        if (DTrace < 0.) {
          c_dil    = sqrt(K * pvol[idx] / pmass[idx]);
          p_q[idx] = (C0 * fabs(c_dil * DTrace * dx_ave) +
                      C1 * (DTrace * DTrace * dx_ave * dx_ave)) *
                     (pmass[idx] / pvol[idx]);
        }
      }
    }
  }
}

void
FractureMPM::computeContactArea(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  // six indices for each of the faces
  double bndyCArea[6] = { 0, 0, 0, 0, 0, 0 };

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing computeContactArea on patch " << patch->getID()
                 << "\t\t\t MPM" << std::endl;
    }

    Vector dx      = patch->dCell();
    double cellvol = dx.x() * dx.y() * dx.z();

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      constNCVariable<double> gVolume, Gvolume;

      new_dw
        ->get(gVolume, d_mpmLabels->gVolumeLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(Gvolume,
                  d_mpmLabels->GVolumeLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0); // for FractureMPM

      for (auto& face : d_boundaryTractionFaces) {
        int iface = (int)(face);

        // Check if the face is on an external boundary
        if (patch->getBCType(face) == Patch::Neighbor) {
          continue;
        }

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,
        // so do the traction accumulation . . .
        // loop cells to find boundary areas
        IntVector projlow, projhigh;
        patch->getFaceCells(face, 0, projlow, projhigh);
        // Vector norm = face_norm(face);

        for (int i = projlow.x(); i < projhigh.x(); i++) {
          for (int j = projlow.y(); j < projhigh.y(); j++) {
            for (int k = projlow.z(); k < projhigh.z(); k++) {
              IntVector ijk(i, j, k);
              double nodevol = gVolume[ijk] + Gvolume[ijk];
              if (nodevol > 0) // FIXME: uses node index to get node volume ...
              {
                const double celldepth = dx[iface / 2];
                bndyCArea[iface] += cellvol / celldepth;
              }
            }
          }
        }

      } // faces
    }   // materials
  }     // patches

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (std::list<Patch::FaceType>::const_iterator ftit(
         d_boundaryTractionFaces.begin());
       ftit != d_boundaryTractionFaces.end();
       ftit++) {
    int iface = (int)(*ftit);
    new_dw->put(sum_vartype(bndyCArea[iface]),
                d_mpmLabels->BndyContactCellAreaLabel[iface]);
  }
}

void
FractureMPM::computeInternalForce(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset*,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  // node based forces
  Vector bndyForce[6];
  Vector bndyTraction[6];
  for (int iface = 0; iface < 6; iface++) {
    bndyForce[iface]    = Vector(0.);
    bndyTraction[iface] = Vector(0.);
  }
  double partvoldef = 0.;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing computeInternalForce on patch " << patch->getID()
                 << "\t\t\t MPM" << std::endl;
    }

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0]        = 1.0 / dx.x();
    oodx[1]        = 1.0 / dx.y();
    oodx[2]        = 1.0 / dx.z();
    double cellvol = dx.x() * dx.y() * dx.z();
    Matrix3 Id;
    Id.Identity();

    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    NCVariable<Matrix3> gstressglobal;
    constNCVariable<double> gMassglobal;
    new_dw->get(gMassglobal,
                d_mpmLabels->gMassLabel,
                d_materialManager->getAllInOneMaterial()->get(0),
                patch,
                Ghost::None,
                0);
    new_dw->allocateAndPut(gstressglobal,
                           d_mpmLabels->gStressForSavingLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    // gstressglobal.initialize(Matrix3(0.));

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      // Create arrays for the particle position, volume
      // and the constitutive model
      constParticleVariable<Point> px;
      constParticleVariable<double> pvol, pmass;
      constParticleVariable<double> p_pressure;
      constParticleVariable<double> p_q;
      constParticleVariable<Matrix3> pstress;
      constParticleVariable<Matrix3> psize;
      constParticleVariable<Matrix3> pDeformationMeasure;
      NCVariable<Vector> internalforce;
      NCVariable<Matrix3> gstress;
      constNCVariable<double> gMass;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(px, d_mpmLabels->pXLabel, pset);
      old_dw->get(pmass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pvol, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pstress, d_mpmLabels->pStressLabel, pset);
      old_dw->get(psize, d_mpmLabels->pSizeLabel, pset);
      new_dw->get(gMass, d_mpmLabels->gMassLabel, dwi, patch, Ghost::None, 0);
      old_dw->get(pDeformationMeasure, d_mpmLabels->pDefGradLabel, pset);

      new_dw->allocateAndPut(gstress,
                             d_mpmLabels->gStressForSavingLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(internalforce,
                             d_mpmLabels->gInternalForceLabel,
                             dwi,
                             patch);
      // gstress.initialize(Matrix3(0.));

      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpmLabels->pgCodeLabel, pset);
      constNCVariable<double> Gmass;
      new_dw->get(Gmass, d_mpmLabels->GMassLabel, dwi, patch, Ghost::None, 0);
      NCVariable<Vector> Ginternalforce;
      new_dw->allocateAndPut(Ginternalforce,
                             d_mpmLabels->GInternalForceLabel,
                             dwi,
                             patch);

      if (d_mpmFlags->d_withICE) {
        new_dw->get(p_pressure, d_mpmLabels->pPressureLabel, pset);
      } else {
        ParticleVariable<double> p_pressure_create;
        new_dw->allocateTemporary(p_pressure_create, pset);
        for (ParticleSubset::iterator it = pset->begin(); it != pset->end();
             it++) {
          p_pressure_create[*it] = 0.0;
        }
        p_pressure = p_pressure_create; // reference created data
      }

      if (d_mpmFlags->d_artificialViscosity) {
        old_dw->get(p_q, d_mpmLabels->p_qLabel, pset);
      } else {
        ParticleVariable<double> p_q_create;
        new_dw->allocateTemporary(p_q_create, pset);
        for (ParticleSubset::iterator it = pset->begin(); it != pset->end();
             it++) {
          p_q_create[*it] = 0.0;
        }
        p_q = p_q_create; // reference created data
      }

      internalforce.initialize(Vector(0, 0, 0));
      Ginternalforce.initialize(Vector(0, 0, 0));

      Matrix3 stressmass;
      Matrix3 stresspress;

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(
          px[idx],
          ni,
          S,
          d_S,
          psize[idx],
          pDeformationMeasure[idx]);

        stressmass = pstress[idx] * pmass[idx];
        // stresspress = pstress[idx] + Id*p_pressure[idx];
        stresspress = pstress[idx] + Id * p_pressure[idx] - Id * p_q[idx];
        partvoldef += pvol[idx];

        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          if (patch->containsNode(ni[k])) {
            Vector div(d_S[k].x() * oodx[0],
                       d_S[k].y() * oodx[1],
                       d_S[k].z() * oodx[2]);
            if (pgCode[idx][k] == 1) {
              internalforce[ni[k]] -= (div * stresspress) * pvol[idx];
            } else if (pgCode[idx][k] == 2) {
              Ginternalforce[ni[k]] -= (div * stresspress) * pvol[idx];
            }
            gstress[ni[k]] += stressmass * S[k];
          }
        }
      }

      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        gstressglobal[c] += gstress[c];
        gstress[c] /= (gMass[c] + Gmass[c]); // add in addtional field
      }

      // save boundary forces before apply symmetry boundary condition.
      for (auto& face : d_boundaryTractionFaces) {

        // Check if the face is on an external boundary
        if (patch->getBCType(face) == Patch::Neighbor) {
          continue;
        }

        const int iface = (int)face;

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,
        // so do the traction accumulation . . .
        // loop nodes to find forces
        IntVector projlow, projhigh;
        patch->getFaceNodes(face, 0, projlow, projhigh);
        Vector norm = face_norm(face);

        for (int i = projlow.x(); i < projhigh.x(); i++) {
          for (int j = projlow.y(); j < projhigh.y(); j++) {
            for (int k = projlow.z(); k < projhigh.z(); k++) {
              IntVector ijk(i, j, k);

              // flip sign so that pushing on boundary gives positive force
              bndyForce[iface] -= internalforce[ijk];
            }
          }
        }

        patch->getFaceCells(face, 0, projlow, projhigh);
        for (int i = projlow.x(); i < projhigh.x(); i++) {
          for (int j = projlow.y(); j < projhigh.y(); j++) {
            for (int k = projlow.z(); k < projhigh.z(); k++) {
              IntVector ijk(i, j, k);

              double celldepth =
                dx[iface / 2]; // length in direction perpendicular to boundary
              double dA_c = cellvol / celldepth; // cell based volume

              for (int ic = 0; ic < 3; ic++) {
                for (int jc = 0; jc < 3; jc++) {
                  bndyTraction[iface][ic] +=
                    gstress[ijk](ic, jc) * norm[jc] * dA_c;
                }
              }
            }
          }
        }

      } // faces
      string interp_type = d_mpmFlags->d_interpolatorType;
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Symmetric",
                              internalforce,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Symmetric",
                              Ginternalforce,
                              interp_type);

#ifdef KUMAR
      internalforce.initialize(Vector(0, 0, 0));
      Ginternalforce.initialize(Vector(0, 0, 0));
#endif
    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gstressglobal[c] /= gMassglobal[c];
    }
  }
  new_dw->put(sum_vartype(partvoldef), d_mpmLabels->TotalVolumeDeformedLabel);

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (std::list<Patch::FaceType>::const_iterator ftit(
         d_boundaryTractionFaces.begin());
       ftit != d_boundaryTractionFaces.end();
       ftit++) {
    int iface = (int)(*ftit);
    new_dw->put(sumvec_vartype(bndyForce[iface]),
                d_mpmLabels->BndyForceLabel[iface]);

    sum_vartype bndyContactCellArea_iface;
    new_dw->get(bndyContactCellArea_iface,
                d_mpmLabels->BndyContactCellAreaLabel[iface]);

    if (bndyContactCellArea_iface > 0) {
      bndyTraction[iface] /= bndyContactCellArea_iface;
    }

    new_dw->put(sumvec_vartype(bndyTraction[iface]),
                d_mpmLabels->BndyTractionLabel[iface]);

    double bndyContactArea_iface = bndyContactCellArea_iface;
    if (bndyTraction[iface].length2() > 0) {
      bndyContactArea_iface =
        ::sqrt(bndyForce[iface].length2() / bndyTraction[iface].length2());
    }
    new_dw->put(sum_vartype(bndyContactArea_iface),
                d_mpmLabels->BndyContactAreaLabel[iface]);
  }
}

void
FractureMPM::computeAndIntegrateAcceleration(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset*,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing integrateAcceleration on patch " << patch->getID()
                 << "\t\t\t MPM" << std::endl;
    }

    Vector gravity = d_mpmFlags->d_gravity;
    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      delt_vartype delT;
      old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

      // Get required variables for this patch
      constNCVariable<Vector> velocity;
      constNCVariable<Vector> internalforce;
      constNCVariable<Vector> externalforce;
      constNCVariable<double> mass;

      // for FractureMPM
      constNCVariable<Vector> Gvelocity;
      new_dw->get(internalforce,
                  d_mpmLabels->gInternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(externalforce,
                  d_mpmLabels->gExternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(mass, d_mpmLabels->gMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(velocity,
                  d_mpmLabels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(Gvelocity,
                  d_mpmLabels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      NCVariable<Vector> acceleration;
      NCVariable<Vector> velocity_star;
      NCVariable<Vector> Gvelocity_star;
      new_dw->allocateAndPut(velocity_star,
                             d_mpmLabels->gVelocityStarLabel,
                             dwi,
                             patch);
      velocity_star.initialize(Vector(0.0));
      new_dw->allocateAndPut(Gvelocity_star,
                             d_mpmLabels->GVelocityStarLabel,
                             dwi,
                             patch);
      Gvelocity_star.initialize(Vector(0.0));

      // Create variables for the results
      new_dw->allocateAndPut(acceleration,
                             d_mpmLabels->gAccelerationLabel,
                             dwi,
                             patch);
      acceleration.initialize(Vector(0., 0., 0.));

      // for FractureMPM
      constNCVariable<double> Gmass;
      constNCVariable<Vector> Ginternalforce;
      constNCVariable<Vector> Gexternalforce;
      new_dw->get(Gmass, d_mpmLabels->GMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(Ginternalforce,
                  d_mpmLabels->GInternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(Gexternalforce,
                  d_mpmLabels->GExternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      NCVariable<Vector> Gacceleration;
      new_dw->allocateAndPut(Gacceleration,
                             d_mpmLabels->GAccelerationLabel,
                             dwi,
                             patch);
      Gacceleration.initialize(Vector(0., 0., 0.));

      string interp_type = d_mpmFlags->d_interpolatorType;
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        // above crack
        acceleration[c] =
          (internalforce[c] + externalforce[c]) / mass[c] + gravity;
        // below crack
        Gacceleration[c] =
          (Ginternalforce[c] + Gexternalforce[c]) / Gmass[c] + gravity;
        // above crack
        velocity_star[c] = velocity[c] + acceleration[c] * delT;
        // below crack
        Gvelocity_star[c] = Gvelocity[c] + Gacceleration[c] * delT;
      }
    }
  }
}

void
FractureMPM::setGridBoundaryConditions(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing setGridBoundaryConditions on patch "
                 << patch->getID() << "\t\t MPM" << std::endl;
    }

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    delt_vartype delT;
    old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gacceleration;
      constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(gacceleration,
                            d_mpmLabels->gAccelerationLabel,
                            dwi,
                            patch);
      new_dw->getModifiable(gVelocity_star,
                            d_mpmLabels->gVelocityStarLabel,
                            dwi,
                            patch);
      new_dw->get(gVelocity,
                  d_mpmLabels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      // for FractureMPM
      NCVariable<Vector> Gvelocity_star, Gacceleration;
      constNCVariable<Vector> Gvelocity;
      new_dw->getModifiable(Gacceleration,
                            d_mpmLabels->GAccelerationLabel,
                            dwi,
                            patch);
      new_dw->getModifiable(Gvelocity_star,
                            d_mpmLabels->GVelocityStarLabel,
                            dwi,
                            patch);
      new_dw->get(Gvelocity,
                  d_mpmLabels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      string interp_type = d_mpmFlags->d_interpolatorType;
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Velocity",
                              gVelocity_star,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Velocity",
                              Gvelocity_star,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Symmetric",
                              gVelocity_star,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Symmetric",
                              Gvelocity_star,
                              interp_type);

      // Now recompute acceleration as the difference between the velocity
      // interpolated to the grid (no bcs applied) and the new velocity_star
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c      = *iter;
        gacceleration[c] = (gVelocity_star[c] - gVelocity[c]) / delT;
        Gacceleration[c] = (Gvelocity_star[c] - Gvelocity[c]) / delT;
      }
    } // matl loop
  }   // patch loop
}

void
FractureMPM::applyExternalLoads(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  // Get the current time
  simTime_vartype simTimeVar;
  old_dw->get(simTimeVar, d_mpmLabels->simulationTimeLabel);
  double time = simTimeVar;

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << std::endl;
  }

  // Calculate the force vector at each particle for each pressure bc
  std::vector<double> forcePerPart;
  std::vector<PressureBC*> pbcP;
  if (d_mpmFlags->d_useLoadCurves) {
    for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
         ii++) {
      string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
      if (bcs_type == "Pressure") {

        // Get the material points per load curve
        PressureBC* pbc = dynamic_cast<PressureBC*>(
          MPMPhysicalBCFactory::mpmPhysicalBCs[ii].get());
        pbcP.push_back(pbc);

        // Calculate the force per particle at current time
        forcePerPart.push_back(pbc->forcePerParticle(time));
      }
    }
  }

  // Loop thru patches to update external force vector
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing applyExternalLoads on patch " << patch->getID()
                 << "\t MPM" << std::endl;
    }

    // Place for user defined loading scenarios to be defined,
    // otherwise pExternalForce is just carried forward.

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      if (d_mpmFlags->d_useLoadCurves) {
        bool do_PressureBCs = false;
        for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
             ii++) {
          string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
          if (bcs_type == "Pressure") {
            do_PressureBCs = true;
          }
        }
        if (do_PressureBCs) {
          // Get the particle position data
          constParticleVariable<Point> px;
          old_dw->get(px, d_mpmLabels->pXLabel, pset);

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);

          // Get the defromation gradient
          constParticleVariable<Matrix3> pDefGrad;
          old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);

          constParticleVariable<Vector> pDisp;
          old_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);

          // Get the external force data and allocate new space for
          // external force
          ParticleVariable<Vector> pExternalForce;
          ParticleVariable<Vector> pExternalForce_new;
          old_dw->getModifiable(pExternalForce,
                                d_mpmLabels->pExternalForceLabel,
                                pset);
          new_dw->allocateAndPut(pExternalForce_new,
                                 d_mpmLabels->pExtForceLabel_preReloc,
                                 pset);

          // Iterate over the particles
          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            int loadCurveID   = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = pExternalForce[idx];
            } else {
              PressureBC* pbc         = pbcP[loadCurveID];
              double force            = forcePerPart[loadCurveID];
              pExternalForce_new[idx] = pbc->getForceVector(px[idx],
                                                            pDisp[idx],
                                                            force,
                                                            time,
                                                            pDefGrad[idx]);
            }
          }

          // Recycle the loadCurveIDs
          ParticleVariable<int> pLoadCurveID_new;
          new_dw->allocateAndPut(pLoadCurveID_new,
                                 d_mpmLabels->pLoadCurveIDLabel_preReloc,
                                 pset);
          pLoadCurveID_new.copyData(pLoadCurveID);
        }
      } else { // Carry forward the old pEF, scale by d_forceIncrementFactor
        // Get the external force data and allocate new space for
        // external force and copy the data
        constParticleVariable<Vector> pExternalForce;
        ParticleVariable<Vector> pExternalForce_new;
        old_dw->get(pExternalForce, d_mpmLabels->pExternalForceLabel, pset);
        new_dw->allocateAndPut(pExternalForce_new,
                               d_mpmLabels->pExtForceLabel_preReloc,
                               pset);

        // Iterate over the particles
        ParticleSubset::iterator iter = pset->begin();
        for (; iter != pset->end(); iter++) {
          particleIndex idx = *iter;
          pExternalForce_new[idx] =
            pExternalForce[idx] * d_mpmFlags->d_forceIncrementFactor;
        }
      }
    } // matl loop
  }   // patch loop
}

// for thermal stress analysis
void
FractureMPM::computeParticleTempFromGrid(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset*,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      constNCVariable<double> gTemperature, GTemperature;
      new_dw->get(gTemperature,
                  d_mpmLabels->gTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperature,
                  d_mpmLabels->GTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      constParticleVariable<Point> px;
      constParticleVariable<Matrix3> psize;
      constParticleVariable<Matrix3> pDeformationMeasure;

      old_dw->get(px, d_mpmLabels->pXLabel, pset);
      old_dw->get(psize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDeformationMeasure, d_mpmLabels->pDefGradLabel, pset);

      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpmLabels->pgCodeLabel, pset);

      ParticleVariable<double> pTempCur;
      new_dw->allocateAndPut(pTempCur, d_mpmLabels->pTempCurrentLabel, pset);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        double pTemp      = 0.0;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(
          px[idx],
          ni,
          S,
          d_S,
          psize[idx],
          pDeformationMeasure[idx]);
        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          IntVector node = ni[k];
          if (pgCode[idx][k] == 1) {
            pTemp += gTemperature[node] * S[k];
          } else if (pgCode[idx][k] == 2) {
            pTemp += GTemperature[node] * S[k];
          }
        }
        pTempCur[idx] = pTemp;
      } // End of loop over iter
    }   // End of loop over m
  }     // End of loop over p
}

void
FractureMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                             const PatchSubset* patches,
                                             const MaterialSubset*,
                                             DataWarehouse* old_dw,
                                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (cout_doing.active()) {
      cout_doing << "Doing interpolateToParticlesAndUpdate on patch "
                 << patch->getID() << "\t MPM" << std::endl;
    }

    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke       = 0;
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    delt_vartype delT;
    old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

    /*
    Material* reactant;
    reactant = d_materialManager->getMaterialByName("reactant");

    bool combustion_problem=false;
    int RMI = -99;
    if(reactant != 0){
      RMI = reactant->getDWIndex();
      combustion_problem=true;
    }
    */

    double move_particles = 1.;
    if (!d_mpmFlags->d_doGridReset) {
      move_particles = 0.;
    }

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> px;
      ParticleVariable<Point> pxnew, pxx;
      constParticleVariable<Vector> pvelocity;
      constParticleVariable<Matrix3> psize;
      ParticleVariable<Vector> pvelocitynew;
      ParticleVariable<Matrix3> psizeNew;
      constParticleVariable<double> pmass, pTemperature;
      ParticleVariable<double> pmassNew, pvolume, pTempNew;
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pids_new;
      constParticleVariable<Vector> pdisp;
      ParticleVariable<Vector> pdispnew;
      ParticleVariable<double> pkineticEnergyDensity;
      constParticleVariable<Matrix3> pDeformationMeasure;

      // for thermal stress analysis
      constParticleVariable<double> pTempCurrent;
      ParticleVariable<double> pTempPreNew;

      // Get the arrays of grid data on which the new part. values depend
      constNCVariable<Vector> gVelocity_star, gacceleration;
      constNCVariable<double> gTemperatureRate, gTemperature, gTemperatureNoBC;
      constNCVariable<double> dTdt, massBurnFrac, frictionTempRate;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(px, d_mpmLabels->pXLabel, pset);
      old_dw->get(pdisp, d_mpmLabels->pDispLabel, pset);
      old_dw->get(pmass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pids, d_mpmLabels->pParticleIDLabel, pset);
      old_dw->get(pvelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      new_dw->get(pDeformationMeasure,
                  d_mpmLabels->pDefGradLabel_preReloc,
                  pset);

      // for thermal stress analysis
      new_dw->get(pTempCurrent, d_mpmLabels->pTempCurrentLabel, pset);
      new_dw->getModifiable(pvolume, d_mpmLabels->pVolumeLabel_preReloc, pset);
      new_dw->allocateAndPut(pvelocitynew,
                             d_mpmLabels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pxnew, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pxx, d_mpmLabels->pXXLabel, pset);
      new_dw->allocateAndPut(pdispnew, d_mpmLabels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(pmassNew, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pids_new,
                             d_mpmLabels->pParticleIDLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempNew,
                             d_mpmLabels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pkineticEnergyDensity,
                             d_mpmLabels->pKineticEnergyDensityLabel,
                             pset);
      // for thermal stress analysis
      new_dw->allocateAndPut(pTempPreNew,
                             d_mpmLabels->pTempPreviousLabel_preReloc,
                             pset);

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      pids_new.copyData(pids);
      old_dw->get(psize, d_mpmLabels->pSizeLabel, pset);
      new_dw->allocateAndPut(psizeNew, d_mpmLabels->pSizeLabel_preReloc, pset);
      psizeNew.copyData(psize);

      new_dw->get(gVelocity_star,
                  d_mpmLabels->gVelocityStarLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gacceleration,
                  d_mpmLabels->gAccelerationLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperatureRate,
                  d_mpmLabels->gTemperatureRateLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperature,
                  d_mpmLabels->gTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperatureNoBC,
                  d_mpmLabels->gTemperatureNoBCLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpmLabels->frictionalWorkLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpmLabels->pgCodeLabel, pset);
      constNCVariable<Vector> Gvelocity_star, Gacceleration;
      constNCVariable<double> GTemperatureRate, GTemperature, GTemperatureNoBC;
      new_dw->get(Gvelocity_star,
                  d_mpmLabels->GVelocityStarLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(Gacceleration,
                  d_mpmLabels->GAccelerationLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperatureRate,
                  d_mpmLabels->GTemperatureRateLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperature,
                  d_mpmLabels->GTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperatureNoBC,
                  d_mpmLabels->GTemperatureNoBCLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      if (d_mpmFlags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpmLabels->dTdt_NCLabel,
                    dwi,
                    patch,
                    Ghost::AroundCells,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpmLabels->massBurnFractionLabel,
                    dwi,
                    patch,
                    Ghost::AroundCells,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create,
                                  patch,
                                  Ghost::AroundCells,
                                  d_numGhostParticles);
        new_dw->allocateTemporary(massBurnFrac_create,
                                  patch,
                                  Ghost::AroundCells,
                                  d_numGhostParticles);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);
        dTdt         = dTdt_create;         // reference created data
        massBurnFrac = massBurnFrac_create; // reference created data
      }

      double Cp       = mpm_matl->getSpecificHeat();
      double rho_init = mpm_matl->getInitialDensity();

      /*
      double rho_frac_min = 0.;
      if(m == RMI){
        rho_frac_min = .1;
      }
      */

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(
          px[idx],
          ni,
          S,
          d_S,
          psize[idx],
          pDeformationMeasure[idx]);

        Vector vel(0.0, 0.0, 0.0);
        Vector acc(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0;
        double burnFraction = 0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          IntVector node = ni[k];
          fricTempRate = frictionTempRate[node] * d_mpmFlags->d_addFrictionWork;
          if (pgCode[idx][k] == 1) {
            vel += gVelocity_star[node] * S[k];
            acc += gacceleration[node] * S[k];
            tempRate +=
              (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
            burnFraction += massBurnFrac[node] * S[k];
          } else if (pgCode[idx][k] == 2) {
            vel += Gvelocity_star[node] * S[k];
            acc += Gacceleration[node] * S[k];
            tempRate +=
              (GTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
            burnFraction += massBurnFrac[node] * S[k];
          }
        }

        // Update the particle's position and velocity
        pxnew[idx]        = px[idx] + vel * delT * move_particles;
        pdispnew[idx]     = pdisp[idx] + vel * delT;
        pvelocitynew[idx] = pvelocity[idx] + acc * delT;
        // pxx is only useful if we're not in normal grid resetting mode.
        pxx[idx]         = px[idx] + pdispnew[idx];
        pTempNew[idx]    = pTemperature[idx] + tempRate * delT;
        pTempPreNew[idx] = pTempCurrent[idx]; // for thermal stress

        if (cout_heat.active()) {
          cout_heat << "FractureMPM::Particle = " << idx
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate << " dT = " << (tempRate * delT)
                    << " T_new = " << pTempNew[idx] << std::endl;
        }

        double rho;
        if (pvolume[idx] > 0.) {
          rho = pmass[idx] / pvolume[idx];
        } else {
          rho = rho_init;
        }
        pkineticEnergyDensity[idx] = 0.5 * rho * pvelocitynew[idx].length2();
        pmassNew[idx]              = Max(pmass[idx] * (1. - burnFraction), 0.);
        pvolume[idx]               = pmassNew[idx] / rho;

        thermal_energy += pTemperature[idx] * pmass[idx] * Cp;
        ke += .5 * pmass[idx] * pvelocitynew[idx].length2();
        CMX = CMX + (pxnew[idx] * pmass[idx]).asVector();
        totalMom += pvelocitynew[idx] * pmass[idx];
      }

      // Delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        if (pmassNew[idx] <= d_mpmFlags->d_minPartMass) {
          delset->addParticle(idx);
        }
        if (pvelocitynew[idx].length() > d_mpmFlags->d_maxVel) {
          pvelocitynew[idx] = pvelocity[idx];
        }
      }

      new_dw->deleteParticles(delset);
      //__________________________________
      //  particle debugging label-- carry forward
      if (d_mpmFlags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_mpmLabels->pColorLabel, pset);
        new_dw->allocateAndPut(pColor_new,
                               d_mpmLabels->pColorLabel_preReloc,
                               pset);
        pColor_new.copyData(pColor);
      }
    }

    // DON'T MOVE THESE!!!
    new_dw->put(sum_vartype(ke), d_mpmLabels->KineticEnergyLabel);
    new_dw->put(sum_vartype(thermal_energy), d_mpmLabels->ThermalEnergyLabel);
    new_dw->put(sumvec_vartype(CMX), d_mpmLabels->CenterOfMassPositionLabel);
    new_dw->put(sumvec_vartype(totalMom), d_mpmLabels->TotalMomentumLabel);

    // std::cout << "Solid mass lost this timestep = " << massLost << std::endl;
    // std::cout << "Solid momentum after advection = " << totalMom <<
    // std::endl;

    // std::cout << "THERMAL ENERGY " << thermal_energy << std::endl;
  }
}

void
FractureMPM::initialErrorEstimate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* /*matls*/,
                                  DataWarehouse*,
                                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    if (amr_doing.active()) {
      amr_doing << "Doing FractureMPM::initialErrorEstimate on patch "
                << patch->getID() << std::endl;
    }

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
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
      constParticleVariable<Point> px;
      new_dw->get(px, d_mpmLabels->pXLabel, pset);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        refineFlag[patch->getLevel()->getCellIndex(px[*iter])] = true;
        refinePatch->set();
      }
    }
  }
}

void
FractureMPM::errorEstimate(const ProcessorGroup* group,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  // coarsen the errorflag.

  if (cout_doing.active()) {
    cout_doing << "Doing FractureMPM::errorEstimate" << '\n';
  }

  const Level* level = getLevel(patches);
  if (level->getIndex() == level->getGrid()->numLevels() - 1) {
    // on finest level, we do the same thing as initialErrorEstimate, so call it
    initialErrorEstimate(group, patches, matls, old_dw, new_dw);
  } else {
    const Level* fineLevel = level->getFinerLevel().get_rep();

    for (int p = 0; p < patches->size(); p++) {
      const Patch* coarsePatch = patches->get(p);

      if (amr_doing.active()) {
        amr_doing << "Doing FractureMPM::errorEstimate on patch "
                  << coarsePatch->getID() << std::endl;
      }

      // Find the overlapping regions...

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

      for (size_t i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        // Get the particle data
        constCCVariable<int> fineErrorFlag;
        new_dw->get(fineErrorFlag,
                    d_regridder->getRefineFlagLabel(),
                    0,
                    finePatch,
                    Ghost::None,
                    0);

        IntVector fl(finePatch->getExtraCellLowIndex());
        IntVector fh(finePatch->getExtraCellHighIndex());
        IntVector l(fineLevel->mapCellToCoarser(fl));
        IntVector h(fineLevel->mapCellToCoarser(fh));
        l = Max(l, coarsePatch->getExtraCellLowIndex());
        h = Min(h, coarsePatch->getExtraCellHighIndex());

        for (CellIterator iter(l, h); !iter.done(); iter++) {
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
        }
      } // fine patch loop
    }   // coarse patch loop
  }
}

void
FractureMPM::refine(const ProcessorGroup*,
                    const PatchSubset* patches,
                    const MaterialSubset* /*matls*/,
                    DataWarehouse*,
                    DataWarehouse* new_dw)
{
  // just create a particle subset if one doesn't exist
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();

      if (cout_doing.active()) {
        cout_doing << "Doing refine on patch " << patch->getID()
                   << " material # = " << dwi << std::endl;
      }

      // this is a new patch, so create empty particle variables.
      if (!new_dw->haveParticleSubset(dwi, patch)) {
        ParticleSubset* pset = new_dw->createParticleSubset(0, dwi, patch);

        // Create arrays for the particle data
        ParticleVariable<Point> px;
        ParticleVariable<double> pmass, pvolume, pTemperature;
        ParticleVariable<Vector> pvelocity, pexternalforce, psize, pdisp;
        ParticleVariable<int> pLoadCurve;
        ParticleVariable<long64> pID;
        ParticleVariable<Matrix3> pdeform, pstress;

        new_dw->allocateAndPut(px, d_mpmLabels->pXLabel, pset);
        new_dw->allocateAndPut(pmass, d_mpmLabels->pMassLabel, pset);
        new_dw->allocateAndPut(pvolume, d_mpmLabels->pVolumeLabel, pset);
        new_dw->allocateAndPut(pvelocity, d_mpmLabels->pVelocityLabel, pset);
        new_dw->allocateAndPut(pTemperature,
                               d_mpmLabels->pTemperatureLabel,
                               pset);
        new_dw->allocateAndPut(pexternalforce,
                               d_mpmLabels->pExternalForceLabel,
                               pset);
        new_dw->allocateAndPut(pID, d_mpmLabels->pParticleIDLabel, pset);
        new_dw->allocateAndPut(pdisp, d_mpmLabels->pDispLabel, pset);
        new_dw->allocateAndPut(pdeform, d_mpmLabels->pDefGradLabel, pset);
        new_dw->allocateAndPut(pstress, d_mpmLabels->pStressLabel, pset);
        if (d_mpmFlags->d_useLoadCurves) {
          new_dw->allocateAndPut(pLoadCurve,
                                 d_mpmLabels->pLoadCurveIDLabel,
                                 pset);
        }
        new_dw->allocateAndPut(psize, d_mpmLabels->pSizeLabel, pset);
      }
    }
  }

} // end refine()
