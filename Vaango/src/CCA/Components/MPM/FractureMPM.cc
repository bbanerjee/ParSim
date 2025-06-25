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

// usage: export SCI_DEBUG="MPM:+,FractureMPM:+"
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
}

void
FractureMPM::problemSetup(const ProblemSpecP& prob_spec,
                          const ProblemSpecP& restart_prob_spec,
                          GridP& grid,
                          const std::string& input_ups_dir)
{
  SerialMPM::problemSetup(prob_spec, restart_prob_spec, grid, input_ups_dir);

  // for FractureMPM
  crackModel = std::make_unique<Crack>(prob_spec,
                                       d_materialManager,
                                       d_output,
                                       d_mpm_labels.get(),
                                       d_mpm_flags.get());
}

void
FractureMPM::outputProblemSpec(ProblemSpecP& root_ps)
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

  contactModel->outputProblemSpec(mpm_ps);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps   = physical_bc_ps->appendChild("MPM");
  for (auto& bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    bc->outputProblemSpec(mpm_ph_bc_ps);
  }

  // **TODO** Add crack information
  ProblemSpecP uda_ps = root_ps->findBlock("DataArchiver");
  crackModel->outputProblemSpec(mpm_ps, uda_ps);
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

  t->computes(d_mpm_labels->partCountLabel);
  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_mpm_labels->pFiberDirLabel);
  t->computes(d_mpm_labels->pMassLabel);
  t->computes(d_mpm_labels->pVolumeLabel);
  t->computes(d_mpm_labels->pTemperatureLabel);
  t->computes(d_mpm_labels->pTempPreviousLabel); // for thermal stress
  t->computes(d_mpm_labels->pdTdtLabel);
  t->computes(d_mpm_labels->pDispLabel);
  t->computes(d_mpm_labels->pVelocityLabel);
  t->computes(d_mpm_labels->pAccelerationLabel);
  t->computes(d_mpm_labels->pExternalForceLabel);
  t->computes(d_mpm_labels->pParticleIDLabel);
  t->computes(d_mpm_labels->pStressLabel);
  t->computes(d_mpm_labels->pSizeLabel);
  t->computes(d_mpm_labels->pDispGradsLabel);
  t->computes(d_mpm_labels->pStrainEnergyDensityLabel);
  t->computes(d_mpm_labels->delTLabel, level.get_rep());
  t->computes(d_mpm_labels->pCellNAPIDLabel, zeroth_matl);

  // Debugging Scalar
  if (d_mpm_flags->d_withColor) {
    t->computes(d_mpm_labels->pColorLabel);
  }

  if (d_mpm_flags->d_useLoadCurves) {
    // Computes the load curve ID associated with each particle
    t->computes(d_mpm_labels->pLoadCurveIDLabel);
  }

  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    // Computes accumulated strain energy
    t->computes(d_mpm_labels->AccStrainEnergyLabel);
  }

  // artificial damping coeff initialized to 0.0
  if (cout_dbg.active()) {
    cout_doing << "Artificial Damping Coeff = "
               << d_mpm_flags->d_artificialDampCoeff
               << " 8 or 27 = " << d_mpm_flags->d_8or27 << std::endl;
  }

  const PatchSet* patches = level->eachPatch();
  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    // For velocity gradient and deformation gradient
    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Constitutive model computes and requires
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  // Add initialization of body force and coriolis importance terms
  // These are initialized to zero in ParticleCreator
  t->computes(d_mpm_labels->pCoriolisImportanceLabel);
  t->computes(d_mpm_labels->pBodyForceAccLabel);

  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) {
    delete zeroth_matl; // shouln't happen, but...
  }

  schedulePrintParticleCount(level, sched);

  // for FractureMPM: Discretize crack plane into triangular elements
  t = scinew Task("Crack:CrackDiscretization",
                  crackModel.get(),
                  &Crack::CrackDiscretization);
  crackModel->addComputesAndRequiresCrackDiscretization(
    t,
    level->eachPatch(),
    d_materialManager->allMaterials("MPM"));
  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));

  if (d_mpm_flags->d_useLoadCurves) {
    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
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

    printTask(patches, patch, cout_doing, "FractureMPM::Doing actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpm_labels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      // Initialize deformation gradient
      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

      // Initialize constitutive models
      mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                         mpm_matl,
                                                         new_dw);
      // scalar used for debugging
      if (d_mpm_flags->d_withColor) {
        ParticleVariable<double> pcolor;
        int index            = mpm_matl->getDWIndex();
        ParticleSubset* pset = new_dw->getParticleSubset(index, patch);
        setParticleDefault<double>(pcolor,
                                   d_mpm_labels->pColorLabel,
                                   pset,
                                   new_dw,
                                   0.0);
      }
    }
  }

  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    // Initialize the accumulated strain energy
    new_dw->put(max_vartype(0.0), d_mpm_labels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_mpm_labels->partCountLabel);
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
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  const PatchSet* patches  = level->eachPatch();
  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  // Compute body forces first
  scheduleComputeParticleBodyForce(sched, patches, matls);

  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleParticleVelocityField(sched, patches, matls); // for FractureMPM
  scheduleInterpolateParticlesToGrid(sched, patches, matls);

  // For contact
  scheduleComputeNormals(sched, patches, matls);
  scheduleFindSurfaceParticles(sched, patches, matls);
  scheduleComputeLogisticRegression(sched, patches, matls);

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

  // Schedule compute of the deformation gradient
  scheduleComputeDeformationGradient(sched, patches, matls);

  // Schedule compute of the stress tensor
  scheduleComputeStressTensor(sched, patches, matls);

  d_heatConductionTasks->scheduleCompute(sched, patches, matls);

  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  scheduleCalculateFractureParameters(sched, patches, matls); // for FractureMPM
  scheduleDoCrackPropagation(sched, patches, matls);          // for FractureMPM
  scheduleMoveCracks(sched, patches, matls);                  // for FractureMPM
  scheduleUpdateCrackFront(sched, patches, matls);            // for FractureMPM

  sched->scheduleParticleRelocation(level,
                                    d_mpm_labels->pXLabel_preReloc,
                                    d_particleState_preReloc,
                                    d_mpm_labels->pXLabel,
                                    d_particleState,
                                    d_mpm_labels->pParticleIDLabel,
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

  t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_mpm_labels->pExternalForceLabel, Ghost::None);
  t->computes(d_mpm_labels->pExtForceLabel_preReloc);
  if (d_mpm_flags->d_useLoadCurves) {
    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
    t->needs(Task::OldDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->pLoadCurveIDLabel_preReloc);
  }

  //  t->computes(Task::OldDW, d_mpm_labels->pExternalHeatRateLabel_preReloc);

  sched->addTask(t, patches, matls);
}

// Determine velocity field for each particle-node pair
void
FractureMPM::scheduleParticleVelocityField(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::ParticleVelocityField",
                        crackModel.get(),
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

  t->needs(Task::OldDW,
              d_mpm_labels->pMassLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pVolumeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pVelocityLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pTemperatureLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::NewDW,
              d_mpm_labels->pExtForceLabel_preReloc,
              Ghost::AroundNodes,
              d_numGhostParticles);
  // t->needs(Task::OldDW, d_mpm_labels->pExternalHeatRateLabel,
  // Ghost::AroundNodes,d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);

  t->computes(d_mpm_labels->gMassLabel);
  t->computes(d_mpm_labels->gMassLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpm_labels->gTemperatureLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpm_labels->gVolumeLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpm_labels->gVelocityLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);
  t->computes(d_mpm_labels->gSpecificVolumeLabel);
  t->computes(d_mpm_labels->gVolumeLabel);
  t->computes(d_mpm_labels->gVelocityLabel);
  t->computes(d_mpm_labels->gExternalForceLabel);
  t->computes(d_mpm_labels->gTemperatureLabel);
  t->computes(d_mpm_labels->gTemperatureNoBCLabel);
  t->computes(d_mpm_labels->gTemperatureRateLabel);
  t->computes(d_mpm_labels->gExternalHeatRateLabel);
  t->computes(d_mpm_labels->gNumNearParticlesLabel);
  t->computes(d_mpm_labels->TotalMassLabel);

  // for FractureMPM
  t->needs(Task::OldDW,
              d_mpm_labels->pDispLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::NewDW,
              d_mpm_labels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->computes(d_mpm_labels->GMassLabel);
  t->computes(d_mpm_labels->GSp_volLabel);
  t->computes(d_mpm_labels->GVolumeLabel);
  t->computes(d_mpm_labels->GVelocityLabel);
  t->computes(d_mpm_labels->GExternalForceLabel);
  t->computes(d_mpm_labels->GTemperatureLabel);
  t->computes(d_mpm_labels->GTemperatureNoBCLabel);
  t->computes(d_mpm_labels->GExternalHeatRateLabel);
  t->computes(d_mpm_labels->gDisplacementLabel);
  t->computes(d_mpm_labels->GDisplacementLabel);

  sched->addTask(t, patches, matls);
}

// Check crack contact and adjust grid velocity field
void
FractureMPM::scheduleAdjustCrackContactInterpolated(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::AdjustCrackContactInterpolated",
                        crackModel.get(),
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
                                       d_mpm_labels->gVelocityLabel);
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

  scheduleUnrotateStressAndDeformationRate(sched, patches, matls);

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
    t->computes(d_mpm_labels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_mpm_labels->delTLabel, getLevel(patches));
  t->computes(d_mpm_labels->StrainEnergyLabel);

  sched->addTask(t, patches, matls);

  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }
  if (d_mpm_flags->d_artificialViscosity) {
    scheduleComputeArtificialViscosity(sched, patches, matls);
  }

  scheduleRotateStress(sched, patches, matls);
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
  t->needs(Task::OldDW,
              d_mpm_labels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->computes(d_mpm_labels->pTempCurrentLabel);
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
  t->needs(Task::OldDW, d_mpm_labels->AccStrainEnergyLabel);
  t->needs(Task::NewDW, d_mpm_labels->StrainEnergyLabel);
  t->computes(d_mpm_labels->AccStrainEnergyLabel);
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

  t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->pVolumeLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->needs(Task::NewDW,
              d_mpm_labels->gVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes); // for FractureMPM
  t->needs(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None);
  t->computes(d_mpm_labels->p_qLabel);

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

    t->needs(Task::NewDW, d_mpm_labels->gVolumeLabel, Ghost::None);
    t->needs(Task::NewDW,
                d_mpm_labels->GVolumeLabel,
                Ghost::None); // for FractureMPM
    for (auto& face : d_boundaryTractionFaces) {
      int iface = (int)(face);
      t->computes(d_mpm_labels->BndyContactCellAreaLabel[iface]);
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

  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->needs(Task::NewDW,
              d_mpm_labels->gMassLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain,
              Ghost::None);
  t->needs(Task::OldDW,
              d_mpm_labels->pStressLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pVolumeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pMassLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::OldDW,
              d_mpm_labels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);

  // for FractureMPM
  t->needs(Task::NewDW,
              d_mpm_labels->pgCodeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->needs(Task::NewDW, d_mpm_labels->GMassLabel, Ghost::None);
  t->computes(d_mpm_labels->GInternalForceLabel);
  t->computes(d_mpm_labels->TotalVolumeDeformedLabel);

  if (d_mpm_flags->d_withICE) {
    t->needs(Task::NewDW,
                d_mpm_labels->pPressureLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }

  if (d_mpm_flags->d_artificialViscosity) {
    t->needs(Task::NewDW,
                d_mpm_labels->p_qLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }

  t->computes(d_mpm_labels->gInternalForceLabel);

  for (std::list<Patch::FaceType>::const_iterator ftit(
         d_boundaryTractionFaces.begin());
       ftit != d_boundaryTractionFaces.end();
       ftit++) {
    int iface = (int)(*ftit);
    t->needs(Task::NewDW, d_mpm_labels->BndyContactCellAreaLabel[iface]);
    t->computes(d_mpm_labels->BndyForceLabel[iface]);
    t->computes(d_mpm_labels->BndyContactAreaLabel[iface]);
    t->computes(d_mpm_labels->BndyTractionLabel[iface]);
  }

  t->computes(d_mpm_labels->gStressForSavingLabel);
  t->computes(d_mpm_labels->gStressForSavingLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain);

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
                                                     const PatchSet* patches,
                                                     const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleComputeAndIntegrateAcceleration");

  Task* t = scinew Task("MPM::computeAndIntegrateAcceleration",
                        this,
                        &FractureMPM::computeAndIntegrateAcceleration);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  t->needs(Task::NewDW, d_mpm_labels->gMassLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gInternalForceLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gExternalForceLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->gVelocityLabel, Ghost::None);

  t->needs(Task::NewDW, d_mpm_labels->GMassLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->GVelocityLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->GInternalForceLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->GExternalForceLabel, Ghost::None);

  t->computes(d_mpm_labels->gVelocityStarLabel);
  t->computes(d_mpm_labels->gAccelerationLabel);

  t->computes(d_mpm_labels->GAccelerationLabel);
  t->computes(d_mpm_labels->GVelocityStarLabel);

  sched->addTask(t, patches, matls);
}

// Check crack contact and adjust nodal velocities and accelerations
void
FractureMPM::scheduleAdjustCrackContactIntegrated(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)
{
  Task* t = scinew Task("Crack::AdjustCrackContactIntegrated",
                        crackModel.get(),
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
                                       d_mpm_labels->gVelocityStarLabel);
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
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  t->modifies(d_mpm_labels->gAccelerationLabel, mss);
  t->modifies(d_mpm_labels->gVelocityStarLabel, mss);
  t->needs(Task::NewDW, d_mpm_labels->gVelocityLabel, Ghost::None);

  // for FractureMPM
  t->modifies(d_mpm_labels->GAccelerationLabel, mss);
  t->modifies(d_mpm_labels->GVelocityStarLabel, mss);
  t->needs(Task::NewDW, d_mpm_labels->GVelocityLabel, Ghost::None);

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

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  t->needs(Task::NewDW,
              d_mpm_labels->gAccelerationLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureRateLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->gTemperatureNoBCLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->frictionalWorkLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->frictionalWorkCrackLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pTemperatureLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pDispLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->modifies(d_mpm_labels->pVolumeLabel_preReloc);
  // for thermal stress analysis
  t->needs(Task::NewDW, d_mpm_labels->pTempCurrentLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel_preReloc, Ghost::None);

  // for FractureMPM
  t->needs(Task::NewDW,
              d_mpm_labels->GAccelerationLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GVelocityStarLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GTemperatureRateLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GTemperatureLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW,
              d_mpm_labels->GTemperatureNoBCLabel,
              Ghost::AroundCells,
              d_numGhostNodes);
  t->needs(Task::NewDW, d_mpm_labels->pgCodeLabel, Ghost::None);

  if (d_mpm_flags->d_withICE) {
    t->needs(Task::NewDW,
                d_mpm_labels->dTdt_NCLabel,
                Ghost::AroundCells,
                d_numGhostNodes);
    t->needs(Task::NewDW,
                d_mpm_labels->massBurnFractionLabel,
                Ghost::AroundCells,
                d_numGhostNodes);
  }

  t->computes(d_mpm_labels->pDispLabel_preReloc);
  t->computes(d_mpm_labels->pVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pAccelerationLabel_preReloc);
  t->computes(d_mpm_labels->pXLabel_preReloc);
  t->computes(d_mpm_labels->pParticleIDLabel_preReloc);
  t->computes(d_mpm_labels->pTemperatureLabel_preReloc);
  t->computes(d_mpm_labels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpm_labels->pMassLabel_preReloc);
  t->computes(d_mpm_labels->pSizeLabel_preReloc);
  t->computes(d_mpm_labels->pXXLabel);
  t->computes(d_mpm_labels->pKineticEnergyDensityLabel); // for FractureMPM

  t->computes(d_mpm_labels->KineticEnergyLabel);
  t->computes(d_mpm_labels->ThermalEnergyLabel);
  t->computes(d_mpm_labels->CenterOfMassPositionLabel);
  t->computes(d_mpm_labels->TotalMomentumLabel);

  // debugging scalar
  if (d_mpm_flags->d_withColor) {
    t->needs(Task::OldDW, d_mpm_labels->pColorLabel, Ghost::None);
    t->computes(d_mpm_labels->pColorLabel_preReloc);
  }

  // Carry forward external heat flux for switch from explicit to implicit
  t->needs(Task::OldDW, d_mpm_labels->pExternalHeatFluxLabel, Ghost::None);
  t->computes(d_mpm_labels->pExternalHeatFluxLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
FractureMPM::scheduleCalculateFractureParameters(SchedulerP& sched,
                                                 const PatchSet* patches,
                                                 const MaterialSet* matls)
{
  // Get nodal solutions
  Task* t = scinew Task("Crack::GetNodalSolutions",
                        crackModel.get(),
                        &Crack::GetNodalSolutions);
  crackModel->addComputesAndRequiresGetNodalSolutions(t, patches, matls);
  sched->addTask(t, patches, matls);

  // cfnset & cfsset
  t = scinew Task("Crack::CrackFrontNodeSubset",
                  crackModel.get(),
                  &Crack::CrackFrontNodeSubset);
  crackModel->addComputesAndRequiresCrackFrontNodeSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Compute fracture parameters (J, K,...)
  t = scinew Task("Crack::CalculateFractureParameters",
                  crackModel.get(),
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
                        crackModel.get(),
                        &Crack::PropagateCrackFrontPoints);
  crackModel->addComputesAndRequiresPropagateCrackFrontPoints(t,
                                                              patches,
                                                              matls);
  sched->addTask(t, patches, matls);

  // Construct the new crack-front elems and new crack-front segments.
  // The new crack-front is temporary, and will be updated after moving cracks
  t = scinew Task("Crack::ConstructNewCrackFrontElems",
                  crackModel.get(),
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
                        crackModel.get(),
                        &Crack::CrackPointSubset);
  crackModel->addComputesAndRequiresCrackPointSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Move crack points
  t = scinew Task("Crack::MoveCracks", crackModel.get(), &Crack::MoveCracks);
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
                        crackModel.get(),
                        &Crack::CrackFrontNodeSubset);
  crackModel->addComputesAndRequiresCrackFrontNodeSubset(t, patches, matls);
  sched->addTask(t, patches, matls);

  // Recollect crack-front segments, discarding the dead segments,
  // calculating normals, indexes and so on
  t = scinew Task("Crack::RecollectCrackFrontSegments",
                  crackModel.get(),
                  &Crack::RecollectCrackFrontSegments);
  crackModel->addComputesAndRequiresRecollectCrackFrontSegments(t,
                                                                patches,
                                                                matls);
  sched->addTask(t, patches, matls);
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
  old_dw->get(accStrainEnergy, d_mpm_labels->AccStrainEnergyLabel);

  // Get the incremental strain energy from the new datawarehouse
  sum_vartype incStrainEnergy;
  new_dw->get(incStrainEnergy, d_mpm_labels->StrainEnergyLabel);

  // Add the two a put into new dw
  double totalStrainEnergy = (double)accStrainEnergy + (double)incStrainEnergy;
  new_dw->put(max_vartype(totalStrainEnergy),
              d_mpm_labels->AccStrainEnergyLabel);
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
    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

    NCVariable<double> gMassglobal, gtempglobal, gVolumeglobal;
    NCVariable<Vector> gvelglobal;
    new_dw->allocateAndPut(gMassglobal,
                           d_mpm_labels->gMassLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gtempglobal,
                           d_mpm_labels->gTemperatureLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gVolumeglobal,
                           d_mpm_labels->gVolumeLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gvelglobal,
                           d_mpm_labels->gVelocityLabel,
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
      constParticleVariable<double> pMass, pVolume, pTemperature;
      constParticleVariable<Vector> pVelocity, pExternalForce, pDisplacement;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDeformationMeasure;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpm_labels->pXLabel);

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pVolume, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->get(pExternalForce, d_mpm_labels->pExtForceLabel_preReloc, pset);
      old_dw->get(pDeformationMeasure, d_mpm_labels->pDefGradLabel, pset);

      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);
      old_dw->get(pDisplacement, d_mpm_labels->pDispLabel, pset);

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

      new_dw->allocateAndPut(gMass, d_mpm_labels->gMassLabel, dwi, patch);
      new_dw->allocateAndPut(gSp_vol, d_mpm_labels->gSpecificVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(gVolume, d_mpm_labels->gVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(gVelocity,
                             d_mpm_labels->gVelocityLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperature,
                             d_mpm_labels->gTemperatureLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperatureNoBC,
                             d_mpm_labels->gTemperatureNoBCLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperatureRate,
                             d_mpm_labels->gTemperatureRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gexternalforce,
                             d_mpm_labels->gExternalForceLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gexternalheatrate,
                             d_mpm_labels->gExternalHeatRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gnumnearparticles,
                             d_mpm_labels->gNumNearParticlesLabel,
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

      new_dw->allocateAndPut(Gmass, d_mpm_labels->GMassLabel, dwi, patch);
      new_dw->allocateAndPut(GSp_vol, d_mpm_labels->GSp_volLabel, dwi, patch);
      new_dw->allocateAndPut(Gvolume, d_mpm_labels->GVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(Gvelocity,
                             d_mpm_labels->GVelocityLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(GTemperature,
                             d_mpm_labels->GTemperatureLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(GTemperatureNoBC,
                             d_mpm_labels->GTemperatureNoBCLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gexternalforce,
                             d_mpm_labels->GExternalForceLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gexternalheatrate,
                             d_mpm_labels->GExternalHeatRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gdisplacement,
                             d_mpm_labels->gDisplacementLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(Gdisplacement,
                             d_mpm_labels->GDisplacementLabel,
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
                                         pSize[idx],
                                         pDeformationMeasure[idx]);

        pmom = pVelocity[idx] * pMass[idx];
        total_mom += pVelocity[idx] * pMass[idx];

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          if (patch->containsNode(ni[k])) {
            if (pgCode[idx][k] == 1) { // above crack
              gMass[ni[k]] += pMass[idx] * S[k];
              gVelocity[ni[k]] += pmom * S[k];
              gVolume[ni[k]] += pVolume[idx] * S[k];
              gexternalforce[ni[k]] += pExternalForce[idx] * S[k];
              gTemperature[ni[k]] += pTemperature[idx] * pMass[idx] * S[k];
              // gexternalheatrate[ni[k]] += pexternalheatrate[idx]      * S[k];
              gnumnearparticles[ni[k]] += 1.0;
              gSp_vol[ni[k]] += pSp_vol * pMass[idx] * S[k];
              gdisplacement[ni[k]] += pDisplacement[idx] * pMass[idx] * S[k];
            } else if (pgCode[idx][k] == 2) { // below crack
              Gmass[ni[k]] += pMass[idx] * S[k];
              Gvelocity[ni[k]] += pmom * S[k];
              Gvolume[ni[k]] += pVolume[idx] * S[k];
              Gexternalforce[ni[k]] += pExternalForce[idx] * S[k];
              GTemperature[ni[k]] += pTemperature[idx] * pMass[idx] * S[k];
              // Gexternalheatrate[ni[k]] += pexternalheatrate[idx]      * S[k];
              GSp_vol[ni[k]] += pSp_vol * pMass[idx] * S[k];
              Gdisplacement[ni[k]] += pDisplacement[idx] * pMass[idx] * S[k];
            }
          }
        } // End of loop over k
      }   // End of loop over iter

      string interp_type = d_mpm_flags->d_interpolatorType;
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

      new_dw->put(sum_vartype(totalmass), d_mpm_labels->TotalMassLabel);

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
  double C0 = d_mpm_flags->d_artificialViscCoeff1;
  double C1 = d_mpm_flags->d_artificialViscCoeff2;
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
    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      constNCVariable<Vector> gVelocity;
      ParticleVariable<double> p_q;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Point> px;
      constParticleVariable<double> pMass, pvol;
      constParticleVariable<Matrix3> pDeformationMeasure;

      new_dw->get(gVelocity,
                  d_mpm_labels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostNodes);
      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      new_dw->get(pvol, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(p_q, d_mpm_labels->p_qLabel, pset);
      old_dw->get(pDeformationMeasure, d_mpm_labels->pDefGradLabel, pset);

      // for FractureMPM
      constNCVariable<Vector> Gvelocity;
      new_dw->get(Gvelocity,
                  d_mpm_labels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostNodes);
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);

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
                                                  pSize[idx],
                                                  pDeformationMeasure[idx]);

        // get particle's velocity gradients
        Vector gvel(0., 0., 0.);
        velGrad.set(0.0);
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
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
          c_dil    = sqrt(K * pvol[idx] / pMass[idx]);
          p_q[idx] = (C0 * fabs(c_dil * DTrace * dx_ave) +
                      C1 * (DTrace * DTrace * dx_ave * dx_ave)) *
                     (pMass[idx] / pvol[idx]);
        }
      }
    }
  }
}

void
FractureMPM::computeContactArea(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                [[maybe_unused]] DataWarehouse* old_dw,
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
        ->get(gVolume, d_mpm_labels->gVolumeLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(Gvolume,
                  d_mpm_labels->GVolumeLabel,
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
                d_mpm_labels->BndyContactCellAreaLabel[iface]);
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

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    NCVariable<Matrix3> gstressglobal;
    constNCVariable<double> gMassglobal;
    new_dw->get(gMassglobal,
                d_mpm_labels->gMassLabel,
                d_materialManager->getAllInOneMaterial()->get(0),
                patch,
                Ghost::None,
                0);
    new_dw->allocateAndPut(gstressglobal,
                           d_mpm_labels->gStressForSavingLabel,
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
      constParticleVariable<double> pvol, pMass;
      constParticleVariable<double> p_pressure;
      constParticleVariable<double> p_q;
      constParticleVariable<Matrix3> pstress;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDeformationMeasure;
      NCVariable<Vector> internalforce;
      NCVariable<Matrix3> gstress;
      constNCVariable<double> gMass;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpm_labels->pXLabel);

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pvol, d_mpm_labels->pVolumeLabel, pset);
      old_dw->get(pstress, d_mpm_labels->pStressLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->get(gMass, d_mpm_labels->gMassLabel, dwi, patch, Ghost::None, 0);
      old_dw->get(pDeformationMeasure, d_mpm_labels->pDefGradLabel, pset);

      new_dw->allocateAndPut(gstress,
                             d_mpm_labels->gStressForSavingLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(internalforce,
                             d_mpm_labels->gInternalForceLabel,
                             dwi,
                             patch);
      // gstress.initialize(Matrix3(0.));

      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);
      constNCVariable<double> Gmass;
      new_dw->get(Gmass, d_mpm_labels->GMassLabel, dwi, patch, Ghost::None, 0);
      NCVariable<Vector> Ginternalforce;
      new_dw->allocateAndPut(Ginternalforce,
                             d_mpm_labels->GInternalForceLabel,
                             dwi,
                             patch);

      if (d_mpm_flags->d_withICE) {
        new_dw->get(p_pressure, d_mpm_labels->pPressureLabel, pset);
      } else {
        ParticleVariable<double> p_pressure_create;
        new_dw->allocateTemporary(p_pressure_create, pset);
        for (ParticleSubset::iterator it = pset->begin(); it != pset->end();
             it++) {
          p_pressure_create[*it] = 0.0;
        }
        p_pressure = p_pressure_create; // reference created data
      }

      if (d_mpm_flags->d_artificialViscosity) {
        old_dw->get(p_q, d_mpm_labels->p_qLabel, pset);
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
          pSize[idx],
          pDeformationMeasure[idx]);

        stressmass = pstress[idx] * pMass[idx];
        // stresspress = pstress[idx] + Id*p_pressure[idx];
        stresspress = pstress[idx] + Id * p_pressure[idx] - Id * p_q[idx];
        partvoldef += pvol[idx];

        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
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
      string interp_type = d_mpm_flags->d_interpolatorType;
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
  new_dw->put(sum_vartype(partvoldef), d_mpm_labels->TotalVolumeDeformedLabel);

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (std::list<Patch::FaceType>::const_iterator ftit(
         d_boundaryTractionFaces.begin());
       ftit != d_boundaryTractionFaces.end();
       ftit++) {
    int iface = (int)(*ftit);
    new_dw->put(sumvec_vartype(bndyForce[iface]),
                d_mpm_labels->BndyForceLabel[iface]);

    sum_vartype bndyContactCellArea_iface;
    new_dw->get(bndyContactCellArea_iface,
                d_mpm_labels->BndyContactCellAreaLabel[iface]);

    if (bndyContactCellArea_iface > 0) {
      bndyTraction[iface] /= bndyContactCellArea_iface;
    }

    new_dw->put(sumvec_vartype(bndyTraction[iface]),
                d_mpm_labels->BndyTractionLabel[iface]);

    double bndyContactArea_iface = bndyContactCellArea_iface;
    if (bndyTraction[iface].length2() > 0) {
      bndyContactArea_iface =
        ::sqrt(bndyForce[iface].length2() / bndyTraction[iface].length2());
    }
    new_dw->put(sum_vartype(bndyContactArea_iface),
                d_mpm_labels->BndyContactAreaLabel[iface]);
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

    Vector gravity = d_mpm_flags->d_gravity;
    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      delt_vartype delT;
      old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

      // Get required variables for this patch
      constNCVariable<Vector> velocity;
      constNCVariable<Vector> internalforce;
      constNCVariable<Vector> externalforce;
      constNCVariable<double> mass;

      // for FractureMPM
      constNCVariable<Vector> Gvelocity;
      new_dw->get(internalforce,
                  d_mpm_labels->gInternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(externalforce,
                  d_mpm_labels->gExternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(mass, d_mpm_labels->gMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(velocity,
                  d_mpm_labels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(Gvelocity,
                  d_mpm_labels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      NCVariable<Vector> acceleration;
      NCVariable<Vector> velocity_star;
      NCVariable<Vector> Gvelocity_star;
      new_dw->allocateAndPut(velocity_star,
                             d_mpm_labels->gVelocityStarLabel,
                             dwi,
                             patch);
      velocity_star.initialize(Vector(0.0));
      new_dw->allocateAndPut(Gvelocity_star,
                             d_mpm_labels->GVelocityStarLabel,
                             dwi,
                             patch);
      Gvelocity_star.initialize(Vector(0.0));

      // Create variables for the results
      new_dw->allocateAndPut(acceleration,
                             d_mpm_labels->gAccelerationLabel,
                             dwi,
                             patch);
      acceleration.initialize(Vector(0., 0., 0.));

      // for FractureMPM
      constNCVariable<double> Gmass;
      constNCVariable<Vector> Ginternalforce;
      constNCVariable<Vector> Gexternalforce;
      new_dw->get(Gmass, d_mpm_labels->GMassLabel, dwi, patch, Ghost::None, 0);
      new_dw->get(Ginternalforce,
                  d_mpm_labels->GInternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      new_dw->get(Gexternalforce,
                  d_mpm_labels->GExternalForceLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      NCVariable<Vector> Gacceleration;
      new_dw->allocateAndPut(Gacceleration,
                             d_mpm_labels->GAccelerationLabel,
                             dwi,
                             patch);
      Gacceleration.initialize(Vector(0., 0., 0.));

      string interp_type = d_mpm_flags->d_interpolatorType;
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
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gacceleration;
      constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(gacceleration,
                            d_mpm_labels->gAccelerationLabel,
                            dwi,
                            patch);
      new_dw->getModifiable(gVelocity_star,
                            d_mpm_labels->gVelocityStarLabel,
                            dwi,
                            patch);
      new_dw->get(gVelocity,
                  d_mpm_labels->gVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);
      // for FractureMPM
      NCVariable<Vector> Gvelocity_star, Gacceleration;
      constNCVariable<Vector> Gvelocity;
      new_dw->getModifiable(Gacceleration,
                            d_mpm_labels->GAccelerationLabel,
                            dwi,
                            patch);
      new_dw->getModifiable(Gvelocity_star,
                            d_mpm_labels->GVelocityStarLabel,
                            dwi,
                            patch);
      new_dw->get(Gvelocity,
                  d_mpm_labels->GVelocityLabel,
                  dwi,
                  patch,
                  Ghost::None,
                  0);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      string interp_type = d_mpm_flags->d_interpolatorType;
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
  old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
  double time = simTimeVar;

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << std::endl;
  }

  // Calculate the force vector at each particle for each pressure bc
  std::vector<double> forcePerPart;
  std::vector<PressureBC*> pbcP;
  if (d_mpm_flags->d_useLoadCurves) {
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

      if (d_mpm_flags->d_useLoadCurves) {
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
          old_dw->get(px, d_mpm_labels->pXLabel, pset);

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          // Get the defromation gradient
          constParticleVariable<Matrix3> pDefGrad;
          old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

          constParticleVariable<Vector> pDisp;
          old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

          // Get the external force data and allocate new space for
          // external force
          ParticleVariable<Vector> pExternalForce;
          ParticleVariable<Vector> pExternalForce_new;
          old_dw->getModifiable(pExternalForce,
                                d_mpm_labels->pExternalForceLabel,
                                pset);
          new_dw->allocateAndPut(pExternalForce_new,
                                 d_mpm_labels->pExtForceLabel_preReloc,
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
                                 d_mpm_labels->pLoadCurveIDLabel_preReloc,
                                 pset);
          pLoadCurveID_new.copyData(pLoadCurveID);
        }
      } else { // Carry forward the old pEF, scale by d_forceIncrementFactor
        // Get the external force data and allocate new space for
        // external force and copy the data
        constParticleVariable<Vector> pExternalForce;
        ParticleVariable<Vector> pExternalForce_new;
        old_dw->get(pExternalForce, d_mpm_labels->pExternalForceLabel, pset);
        new_dw->allocateAndPut(pExternalForce_new,
                               d_mpm_labels->pExtForceLabel_preReloc,
                               pset);

        // Iterate over the particles
        ParticleSubset::iterator iter = pset->begin();
        for (; iter != pset->end(); iter++) {
          particleIndex idx = *iter;
          pExternalForce_new[idx] =
            pExternalForce[idx] * d_mpm_flags->d_forceIncrementFactor;
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

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
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
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperature,
                  d_mpm_labels->GTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      constParticleVariable<Point> px;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDeformationMeasure;

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pDeformationMeasure, d_mpm_labels->pDefGradLabel, pset);

      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);

      ParticleVariable<double> pTempCur;
      new_dw->allocateAndPut(pTempCur, d_mpm_labels->pTempCurrentLabel, pset);

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
          pSize[idx],
          pDeformationMeasure[idx]);
        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
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

    auto interpolator = d_mpm_flags->d_interpolator->clone(patch);
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
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    double move_particles = 1.;
    if (!d_mpm_flags->d_doGridReset) {
      move_particles = 0.;
    }

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m)));
      int dwi = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> px;
      ParticleVariable<Point> pxnew, pxx;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pSizeNew;
      constParticleVariable<double> pMass, pTemperature;
      ParticleVariable<double> pMassNew, pVolume, pTempNew;
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pids_new;
      constParticleVariable<Vector> pDisplacement;
      ParticleVariable<Vector> pDisp_new, pVel_new, pAcc_new;
      ParticleVariable<double> pkineticEnergyDensity;
      constParticleVariable<Matrix3> pDeformationMeasure;

      // for thermal stress analysis
      constParticleVariable<double> pTempCurrent;
      ParticleVariable<double> pTempPreNew;

      // Get the arrays of grid data on which the new part. values depend
      constNCVariable<Vector> gVelocity_star, gacceleration;
      constNCVariable<double> gTemperatureRate, gTemperature, gTemperatureNoBC;
      constNCVariable<double> dTdt, massBurnFrac, frictionTempRate, frictionTempRateCrack;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(px, d_mpm_labels->pXLabel, pset);
      old_dw->get(pDisplacement, d_mpm_labels->pDispLabel, pset);
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pids, d_mpm_labels->pParticleIDLabel, pset);
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
      new_dw->get(pDeformationMeasure,
                  d_mpm_labels->pDefGradLabel_preReloc,
                  pset);

      // for thermal stress analysis
      new_dw->get(pTempCurrent, d_mpm_labels->pTempCurrentLabel, pset);
      new_dw->getModifiable(pVolume, d_mpm_labels->pVolumeLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pDisp_new, d_mpm_labels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pVel_new, d_mpm_labels->pVelocityLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pAcc_new, d_mpm_labels->pAccelerationLabel_preReloc, pset);
      new_dw->allocateAndPut(pxnew, d_mpm_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pxx, d_mpm_labels->pXXLabel, pset);
      new_dw->allocateAndPut(pMassNew, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pids_new,
                             d_mpm_labels->pParticleIDLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempNew,
                             d_mpm_labels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pkineticEnergyDensity,
                             d_mpm_labels->pKineticEnergyDensityLabel,
                             pset);
      // for thermal stress analysis
      new_dw->allocateAndPut(pTempPreNew,
                             d_mpm_labels->pTempPreviousLabel_preReloc,
                             pset);

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      pids_new.copyData(pids);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(pSizeNew, d_mpm_labels->pSizeLabel_preReloc, pset);
      pSizeNew.copyData(pSize);

      // Copy needed for switch from explicit to implicit MPM
      constParticleVariable<double> pExtHeatFlux;
      ParticleVariable<double> pExtHeatFlux_new;
      old_dw->get(pExtHeatFlux, d_mpm_labels->pExternalHeatFluxLabel, pset);
      new_dw->allocateAndPut(
        pExtHeatFlux_new, d_mpm_labels->pExternalHeatFluxLabel_preReloc, pset);
      pExtHeatFlux_new.copyData(pExtHeatFlux);

      new_dw->get(gVelocity_star,
                  d_mpm_labels->gVelocityStarLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gacceleration,
                  d_mpm_labels->gAccelerationLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperatureRate,
                  d_mpm_labels->gTemperatureRateLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperature,
                  d_mpm_labels->gTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gTemperatureNoBC,
                  d_mpm_labels->gTemperatureNoBCLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpm_labels->frictionalWorkLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(frictionTempRateCrack,
                  d_mpm_labels->frictionalWorkCrackLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      // for FractureMPM
      constParticleVariable<Short27> pgCode;
      new_dw->get(pgCode, d_mpm_labels->pgCodeLabel, pset);
      constNCVariable<Vector> Gvelocity_star, Gacceleration;
      constNCVariable<double> GTemperatureRate, GTemperature, GTemperatureNoBC;
      new_dw->get(Gvelocity_star,
                  d_mpm_labels->GVelocityStarLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(Gacceleration,
                  d_mpm_labels->GAccelerationLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperatureRate,
                  d_mpm_labels->GTemperatureRateLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperature,
                  d_mpm_labels->GTemperatureLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(GTemperatureNoBC,
                  d_mpm_labels->GTemperatureNoBCLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      if (d_mpm_flags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpm_labels->dTdt_NCLabel,
                    dwi,
                    patch,
                    Ghost::AroundCells,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpm_labels->massBurnFractionLabel,
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
          pSize[idx],
          pDeformationMeasure[idx]);

        Vector vel(0.0, 0.0, 0.0);
        Vector acc(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0;
        double burnFraction = 0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < d_mpm_flags->d_8or27; k++) {
          IntVector node = ni[k];
          fricTempRate = frictionTempRate[node] * d_mpm_flags->d_addFrictionWork;
          fricTempRate += frictionTempRateCrack[node] * d_mpm_flags->d_addFrictionWork;
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
        pDisp_new[idx]     = pDisplacement[idx] + vel * delT;
        pVel_new[idx] = pVelocity[idx] + acc * delT;
        // pxx is only useful if we're not in normal grid resetting mode.
        pxx[idx]         = px[idx] + pDisp_new[idx];
        pTempNew[idx]    = pTemperature[idx] + tempRate * delT;
        pTempPreNew[idx] = pTempCurrent[idx]; // for thermal stress

        pAcc_new[idx] = acc;

        if (cout_heat.active()) {
          cout_heat << "FractureMPM::Particle = " << idx
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate << " dT = " << (tempRate * delT)
                    << " T_new = " << pTempNew[idx] << std::endl;
        }

        double rho;
        if (pVolume[idx] > 0.) {
          rho = pMass[idx] / pVolume[idx];
        } else {
          rho = rho_init;
        }
        pkineticEnergyDensity[idx] = 0.5 * rho * pVel_new[idx].length2();
        pMassNew[idx]              = Max(pMass[idx] * (1. - burnFraction), 0.);
        pVolume[idx]               = pMassNew[idx] / rho;

        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
        ke += .5 * pMass[idx] * pVel_new[idx].length2();
        CMX = CMX + (pxnew[idx] * pMass[idx]).asVector();
        totalMom += pVel_new[idx] * pMass[idx];
      }

      // Delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        if (pMassNew[idx] <= d_mpm_flags->d_minPartMass) {
          delset->addParticle(idx);
        }
        if (pVel_new[idx].length() > d_mpm_flags->d_maxVel) {
          pVel_new[idx] = pVelocity[idx];
        }
      }

      new_dw->deleteParticles(delset);
      //__________________________________
      //  particle debugging label-- carry forward
      if (d_mpm_flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_mpm_labels->pColorLabel, pset);
        new_dw->allocateAndPut(pColor_new,
                               d_mpm_labels->pColorLabel_preReloc,
                               pset);
        pColor_new.copyData(pColor);
      }
    }

    // DON'T MOVE THESE!!!
    new_dw->put(sum_vartype(ke), d_mpm_labels->KineticEnergyLabel);
    new_dw->put(sum_vartype(thermal_energy), d_mpm_labels->ThermalEnergyLabel);
    new_dw->put(sumvec_vartype(CMX), d_mpm_labels->CenterOfMassPositionLabel);
    new_dw->put(sumvec_vartype(totalMom), d_mpm_labels->TotalMomentumLabel);

    // std::cout << "Solid mass lost this timestep = " << massLost << std::endl;
    // std::cout << "Solid momentum after advection = " << totalMom <<
    // std::endl;

    // std::cout << "THERMAL ENERGY " << thermal_energy << std::endl;
  }
}

