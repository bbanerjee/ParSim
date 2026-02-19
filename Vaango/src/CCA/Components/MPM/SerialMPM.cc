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
#include <CCA/Components/MPM/SerialMPM.h>

#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/CohesiveZone/CohesiveZoneTasks.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Components/MPM/DEM/DEMTasks.h>

#include <CCA/Components/MPM/Core/CZLabel.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/Core/MPMUtils.h>

#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>
#include <CCA/Components/MPM/HeatConduction/HeatConductionTasks.h>

#include <CCA/Components/MPM/MMS/MMS.h>

#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>

#include <CCA/Components/MPM/PhysicalBC/FluxBCModel.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>

#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>

#include <CCA/Components/MPM/ReactionDiffusion/ScalarDiffusionTasks.h>

#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/BBox.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/LocalSDF.h>
#include <Core/GeometryPiece/DynamicSDFGeometry.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/UnknownVariable.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/PerPatchVars.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/CubicPolyRoots.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>

#ifdef __NVCC__
#pragma nv_diag_suppress 20011
#pragma nv_diag_suppress 20013
#pragma nv_diag_suppress 20014
#pragma nv_diag_suppress 20015
#endif
#include <Eigen/Dense>
#ifdef __NVCC__
#pragma nv_diag_default 20011
#pragma nv_diag_default 20013
#pragma nv_diag_default 20014
#pragma nv_diag_default 20015
#endif

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <format>

// #define XPIC2_UPDATE
#define CHECK_PARTICLE_DELETION
// #define TIME_COMPUTE_STRESS
#define CHECK_ISFINITE
// #define DEBUG_WITH_PARTICLE_ID
//  constexpr long64 testParticleID = testParticleID;

using namespace Uintah;

// Debug streams
//__________________________________
//  To turn on debug d_mpm_flags
//  csh/tcsh : setenv SCI_DEBUG "MPM:+,SerialMPM:+".....
//  bash     : export SCI_DEBUG="MPM:+,SerialMPM:+" )
//  default is OFF
static DebugStream cout_doing("MPM", false);
static DebugStream cout_dbg("SerialMPM", false);
static DebugStream cout_convert("MPMConv", false);
static DebugStream cout_heat("MPMHeat", false);
static DebugStream amr_doing("AMRMPM", false);
static DebugStream cout_damage("Damage", false);
static DebugStream cout_dem("DEM", false);
Dout mpm_delt_dbg("MPMdelt", "MPM", "Debug MPM delta t", false);

SerialMPM::SerialMPM(const ProcessorGroup* myworld,
                     const MaterialManagerP& materialManager)
  : SimulationCommon(myworld, materialManager)
  , MPMCommon(d_materialManager)
{
  d_mpm_labels = std::make_unique<MPMLabel>();
  d_czLabels   = std::make_unique<CZLabel>();
  d_mpm_flags  = std::make_unique<MPMFlags>(myworld);
}

SerialMPM::~SerialMPM()
{
  MPMPhysicalBCFactory::clean();
  d_analysisModules.clear();
}

/*!----------------------------------------------------------------------
 * problemSetup
 *-----------------------------------------------------------------------*/
void
SerialMPM::problemSetup(const ProblemSpecP& prob_spec,
                        const ProblemSpecP& restart_prob_spec,
                        GridP& grid,
                        const std::string& input_ups_dir)
{
  cout_doing << "Doing problemSetup\t\t\t\t\t MPM"
             << "\n";
  d_scheduler->setPositionVar(d_mpm_labels->pXLabel);

  SimulationCommon::problemSetup( prob_spec );

  ProblemSpecP restart_mat_ps = nullptr;
  ProblemSpecP prob_spec_mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");

  if (prob_spec_mat_ps) {
    restart_mat_ps = prob_spec;
  } else if (restart_prob_spec) {
    d_isRestart    = true;
    restart_mat_ps = restart_prob_spec;
  } else {
    restart_mat_ps = prob_spec;
  }

  ProblemSpecP mpm_soln_ps = restart_mat_ps->findBlock("MPM");
  if (!mpm_soln_ps) {
    std::ostringstream warn;
    warn << "ERROR:MPM:\n missing MPM section in the input file\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // Read all MPM d_mpm_flags (look in MPMFlags.cc)
  d_mpm_flags->readMPMFlags(restart_mat_ps, d_output);
  if (d_mpm_flags->d_integratorType == "implicit") {
    throw ProblemSetupException(
      "Can't use implicit integration with -mpm", __FILE__, __LINE__);
  }

  // convert text representation of face into FaceType
  for (auto faceStr : d_mpm_flags->d_boundaryTractionFaceStrings) {
    Patch::FaceType face = Patch::invalidFace;
    for (auto ft = Patch::startFace; ft <= Patch::endFace;
         ft      = Patch::nextFace(ft)) {
      if (Patch::getFaceName(ft) == faceStr) {
        face = ft;
      }
    }
    if (face != Patch::invalidFace) {
      d_boundaryTractionFaces.push_back(face);
    } else {
      std::cerr << "warning: ignoring unknown face '" << face << "'"
                << "\n";
    }
  }

  // read in AMR d_mpm_flags from the main ups file
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  if (amr_ps) {
    ProblemSpecP mpm_amr_ps = amr_ps->findBlock("MPM");
    if (!mpm_amr_ps) {
      std::ostringstream warn;
      warn << "ERROR:MPM:\n missing MPM section in the AMR section of the "
              "input file\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    mpm_amr_ps->getWithDefault(
      "min_grid_level", d_mpm_flags->d_minGridLevel, 0);
    mpm_amr_ps->getWithDefault(
      "max_grid_level", d_mpm_flags->d_maxGridLevel, 1000);
  }

  if (d_mpm_flags->d_8or27 == 8) {
    d_numGhostParticles = 1;
    d_numGhostNodes     = 1;
  } else {
    d_numGhostParticles = 2;
    d_numGhostNodes     = 2;
  }

  if (d_mpm_flags->d_prescribeDeformation) {
    readPrescribedDeformations(d_mpm_flags->d_prescribedDeformationFile);
  }
  if (d_mpm_flags->d_insertParticles) {
    readInsertParticlesFile(d_mpm_flags->d_insertParticlesFile);
  }

  setParticleGhostLayer(Ghost::AroundNodes, d_numGhostParticles);

  MPMPhysicalBCFactory::create(restart_mat_ps, grid, d_mpm_flags.get());

  contactModel = ContactFactory::create(UintahParallelComponent::d_myworld,
                                        restart_mat_ps,
                                        d_materialManager,
                                        d_mpm_labels.get(),
                                        d_mpm_flags.get());
  contactModel->setContactMaterialAttributes();

  // Creates MPM material w/ constitutive models and damage models
  materialProblemSetup(
    restart_mat_ps, d_mpm_flags.get(), d_isRestart, input_ups_dir);

  // Cohesize zones
  d_cohesiveZoneTasks = std::make_unique<CohesiveZoneTasks>(restart_mat_ps,
                                                            d_materialManager,
                                                            d_mpm_labels.get(),
                                                            d_czLabels.get(),
                                                            d_mpm_flags.get());

  // Scalar diffusion
  d_diffusionTasks = std::make_unique<ScalarDiffusionTasks>(restart_mat_ps,
                                                            d_materialManager,
                                                            d_mpm_labels.get(),
                                                            d_mpm_flags.get(),
                                                            d_numGhostParticles,
                                                            d_numGhostNodes);

  // Heat conduction
  d_heatConductionTasks = std::make_unique<HeatConductionTasks>(
    restart_mat_ps, d_materialManager, d_mpm_labels.get(), d_mpm_flags.get());

  // DEM
  d_demTasks = std::make_unique<DEMTasks>(restart_mat_ps,
                                          d_materialManager,
                                          d_mpm_labels.get(),
                                          d_mpm_flags.get(),
                                          d_numGhostParticles,
                                          d_numGhostNodes);

  // Create deformation gradient computer
  d_defGradComputer = std::make_unique<DeformationGradientComputer>(
    d_materialManager, d_mpm_labels.get(), d_mpm_flags.get());

  // create analysis modules
  if (!d_mpm_flags->d_withICE) { // mpmice handles this
    d_analysisModules =
      AnalysisModuleFactory::create(d_myworld, d_materialManager, prob_spec);

    for (auto& module_name : d_analysisModules) {
      module_name->setComponents(dynamic_cast<SimulationInterface*>(this));
      module_name->problemSetup(prob_spec,
                           restart_prob_spec,
                           grid,
                           d_particleState,
                           d_particleState_preReloc);
    }
  }

  //  create the switching criteria port
  d_switchCriteria =
    dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));

  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(
      restart_mat_ps, restart_prob_spec, d_materialManager);
  }
}

/*!----------------------------------------------------------------------
 * readPrescribedDeformations
 *-----------------------------------------------------------------------*/
void
SerialMPM::readPrescribedDeformations(std::string filename)
{

  if (filename != "") {
    std::ifstream is(filename.c_str());
    if (!is) {
      throw ProblemSetupException(
        "ERROR Opening prescribed deformation file '" + filename + "'\n",
        __FILE__,
        __LINE__);
    } else {
      std::cout << "Reading prescribed deformations from file:" << filename
                << "\n";
    }
    double t0(-1.e9);
    while (is) {
      double t1, F11, F12, F13, F21, F22, F23, F31, F32, F33, Theta, a1, a2, a3;
      is >> t1 >> F11 >> F12 >> F13 >> F21 >> F22 >> F23 >> F31 >> F32 >> F33 >>
        Theta >> a1 >> a2 >> a3;
      if (is) {
        if (t1 <= t0) {
          throw ProblemSetupException("ERROR: Time in prescribed deformation "
                                      "file is not monotomically increasing",
                                      __FILE__,
                                      __LINE__);
        }
        d_prescribedTimes.push_back(t1);
        d_prescribedF.emplace_back(F11, F12, F13, F21, F22, F23, F31, F32, F33);
        d_prescribedAngle.push_back(Theta);
        d_prescribedRotationAxis.emplace_back(a1, a2, a3);
      }
      t0 = t1;
    }
    if (d_prescribedTimes.size() < 2) {
      throw ProblemSetupException(
        "ERROR: Failed to generate valid deformation profile",
        __FILE__,
        __LINE__);
    }
  }
}

/*!----------------------------------------------------------------------
 * readInsertParticlesFile
 *-----------------------------------------------------------------------*/
void
SerialMPM::readInsertParticlesFile(std::string filename)
{

  if (filename != "") {
    std::ifstream is(filename.c_str());
    if (!is) {
      throw ProblemSetupException("ERROR Opening particle insertion file '" +
                                    filename + "'\n",
                                  __FILE__,
                                  __LINE__);
    }

    double t0(-1.e9);
    while (is) {
      double t1, color, transx, transy, transz, v_new_x, v_new_y, v_new_z;
      is >> t1 >> color >> transx >> transy >> transz >> v_new_x >> v_new_y >>
        v_new_z;
      if (is) {
        if (t1 <= t0) {
          throw ProblemSetupException(
            "ERROR: Time in insertParticleFile is not monotomically increasing",
            __FILE__,
            __LINE__);
        }
        d_IPTimes.push_back(t1);
        d_IPColor.push_back(color);
        d_IPTranslate.emplace_back(transx, transy, transz);
        d_IPVelNew.emplace_back(v_new_x, v_new_y, v_new_z);
      }
      t0 = t1;
    }
  }
}

/*!----------------------------------------------------------------------
 * outputProblemSpec
 *-----------------------------------------------------------------------*/
void
SerialMPM::outputProblemSpec(ProblemSpecP& root_ps)
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

  d_cohesiveZoneTasks->outputProblemSpec(mpm_ps);
  d_heatConductionTasks->outputProblemSpec(mpm_ps);
  d_diffusionTasks->outputProblemSpec(mpm_ps);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps   = physical_bc_ps->appendChild("MPM");
  for (auto& bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    bc->outputProblemSpec(mpm_ph_bc_ps);
  }

  //  output data analysis modules. Mpmice handles this
  if (!d_mpm_flags->d_withICE && d_analysisModules.size() != 0) {
    for (auto& am : d_analysisModules) {
      am->outputProblemSpec(root_ps);
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleInitialize
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                 level->getGrid()->numLevels())) {
    return;
  }
  Task* t = scinew Task(
    "MPM::actuallyInitialize", this, &SerialMPM::actuallyInitialize);

  const PatchSet* patches = level->eachPatch();
  printSchedule(patches, cout_doing, "MPM::scheduleInitialize");

  auto* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_mpm_labels->partCountLabel);
  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_mpm_labels->pFiberDirLabel);
  t->computes(d_mpm_labels->pMassLabel);
  t->computes(d_mpm_labels->pVolumeLabel);
  t->computes(d_mpm_labels->pTemperatureLabel);
  t->computes(d_mpm_labels->pTempPreviousLabel); // for thermal stress analysis
  t->computes(d_mpm_labels->pdTdtLabel);
  t->computes(d_mpm_labels->pDispLabel);
  t->computes(d_mpm_labels->pVelocityLabel);
  t->computes(d_mpm_labels->pAccelerationLabel);
  t->computes(d_mpm_labels->pExternalForceLabel);
  t->computes(d_mpm_labels->pParticleIDLabel);
  t->computes(d_mpm_labels->pStressLabel);
  t->computes(d_mpm_labels->pSizeLabel);
  t->computes(d_mpm_labels->pRefinedLabel);

  // DEM labels
  if (d_mpm_flags->d_enableDEM) {
    t->computes(d_mpm_labels->pX0Label);
    t->computes(d_mpm_labels->pRigidBodyIDLabel);
    t->computes(d_mpm_labels->pAngularVelocityLabel);
    t->computes(d_mpm_labels->pTorqueLabel);
    t->computes(d_mpm_labels->pOrientationLabel);
    t->computes(d_mpm_labels->pRadiusLabel);
    t->computes(d_mpm_labels->pInertiaTensorLabel);
  }

  t->computes(d_mpm_labels->delTLabel, level.get_rep());
  t->computes(d_mpm_labels->pCellNAPIDLabel, zeroth_matl);
  t->computes(d_mpm_labels->NC_CCweightLabel, zeroth_matl);

  if (!d_mpm_flags->d_doGridReset) {
    t->computes(d_mpm_labels->gDisplacementLabel);
  }

  // Debugging Scalar
  if (d_mpm_flags->d_withColor) {
    t->computes(d_mpm_labels->pColorLabel);
  }

  // Computes the load curve ID associated with each particle
  if (d_mpm_flags->d_useLoadCurves) {
    t->computes(d_mpm_labels->pLoadCurveIDLabel);
  }

  // Computes accumulated strain energy
  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    t->computes(d_mpm_labels->AccStrainEnergyLabel);
  }

  if (d_mpm_flags->d_artificialViscosity) {
    t->computes(d_mpm_labels->p_qLabel);
  }

  // artificial damping coeff initialized to 0.0
  if (cout_dbg.active()) {
    cout_dbg << "Artificial Damping Coeff = "
             << d_mpm_flags->d_artificialDampCoeff
             << " 8 or 27 = " << d_mpm_flags->d_8or27 << "\n";
  }

  // Scalar diffusion
  d_diffusionTasks->scheduleInitialize(t);

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    // For velocity gradient and deformation gradient
    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Add constitutive model computes
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Add damage model computes
    if (cout_damage.active()) {
      cout_damage << "Damage::Material = " << m << " MPMMaterial = " << mpm_matl
                  << " Do damage = " << mpm_matl->doBasicDamage() << "\n";
    }
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->addInitialComputesAndRequires(
        t, mpm_matl, patches, d_mpm_labels.get());
    }

    // Scalar diffusion
    d_diffusionTasks->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  // Add initialization of body force and coriolis importance terms
  // These are initialized to zero in ParticleCreator
  t->computes(d_mpm_labels->pCoriolisImportanceLabel);
  t->computes(d_mpm_labels->pBodyForceAccLabel);

  // Needed for switch from explicit to implicit MPM
  t->computes(d_mpm_labels->pExternalHeatFluxLabel);

  // For friction contact
  t->computes(d_mpm_labels->pSurfLabel);

  // Add task to scheduler
  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) {
    delete zeroth_matl; // shouln't happen, but...
  }

  // Print particle count
  schedulePrintParticleCount(level, sched);

  // Compute initial stresses due to body forces and recompute the initial
  // deformation gradient
  if (d_mpm_flags->d_initializeStressFromBodyForce) {
    scheduleInitializeStressAndDefGradFromBodyForce(level, sched);
  }

  // Schedule the initialization of pressure BCs per particle
  if (d_mpm_flags->d_useLoadCurves) {
    if (MPMPhysicalBCFactory::mpmPhysicalBCs.size() > 0) {
      std::string bcType = MPMPhysicalBCFactory::mpmPhysicalBCs[0]->getType();
      if (bcType == "Pressure") {
        scheduleInitializePressureBCs(level, sched);
      } else if (bcType == "Moment") {
        scheduleInitializeMomentBCs(level, sched);
      }
    }
  }

  // Scalar diffusion flux particle BCs
  d_diffusionTasks->scheduleInitializeFluxBCs(level, sched);

  // Cohesive zones
  d_cohesiveZoneTasks->scheduleInitialize(level, sched);

  // Data analysis
  if (d_analysisModules.size() != 0) {
    for (auto& am : d_analysisModules) {
      am->scheduleInitialize(sched, level);
    }
  }

  if (d_mpm_flags->d_deleteGeometryObjects) {
    scheduleDeleteGeometryObjects(level, sched);
  }
}

/*!----------------------------------------------------------------------
 * actuallyInitialize
 *-----------------------------------------------------------------------*/
void
SerialMPM::actuallyInitialize(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse*,
                              DataWarehouse* new_dw)
{
  particleIndex totalParticles = 0;

  const Level* level = getLevel(patches);
  IntVector lowNode, highNode;
  level->findInteriorNodeIndexRange(lowNode, highNode);

  // Determine dimensionality for particle splitting
  // To be recognized as 2D, must be in the x-y plane
  // A 1D problem must be in the x-direction.
  d_mpm_flags->d_ndim = 3;
  if (highNode.z() - lowNode.z() == 2) {
    d_mpm_flags->d_ndim = 2;
    if (highNode.y() - lowNode.y() == 2) {
      d_mpm_flags->d_ndim = 1;
    }
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpm_labels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch);

    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NC_CCweight.initialize(0.125);

    std::vector<Patch::FaceType> boundary_faces;
    patch->getBoundaryFaces(boundary_faces);
    for (auto& face : boundary_faces) {
      int mat_id = 0;

      if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {
        for (CellIterator iter = patch->getFaceIterator(face, Patch::FaceNodes);
             !iter.done();
             iter++) {
          NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
        }
      }
    }

    for (int m = 0; m < matls->size(); m++) {
      cout_dbg << "Before get material\n";
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      cout_dbg << "After get material\n";

      if (!d_mpm_flags->d_doGridReset) {
        NCVariable<Vector> gDisplacement;
        new_dw->allocateAndPut(gDisplacement,
                               d_mpm_labels->gDisplacementLabel,
                               mpm_matl->getDWIndex(),
                               patch);
        gDisplacement.initialize(Vector(0.));
      }

      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      // Initialize deformation gradient
      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

      // Initialize constitutive models
      mpm_matl->getConstitutiveModel()->initializeCMData(
        patch, mpm_matl, new_dw);

      // Initialize basic damage model
      if (mpm_matl->doBasicDamage()) {
        mpm_matl->getBasicDamageModel()->initializeDamageData(
          patch, mpm_matl, new_dw, d_mpm_labels.get());
      }

      // Diffusion
      d_diffusionTasks->actuallyInitialize(patch, mpm_matl, new_dw);

      // DEM
      d_demTasks->actuallyInitialize(patch, mpm_matl, new_dw);
    } // end matls loop
  }   // end patches loop

  // Only allow axisymmetric runs if the grid is one cell
  // thick in the theta dir.
  if (d_mpm_flags->d_axisymmetric) {
    int num_cells_in_theta = (highNode.z() - lowNode.z()) - 1;
    if (num_cells_in_theta > 1) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <axisymmetric>true</axisymmetric> the \n"
          << "grid can only have one cell in the circumferential direction.\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  }

  // Bulletproofing for extra cells/interpolators/periodic BCs
  auto num_extra_cells = level->getExtraCells();
  auto periodic        = level->getPeriodicBoundaries();
  auto interp_type     = d_mpm_flags->d_interpolatorType;
  if (interp_type == "linear" && num_extra_cells != IntVector(0, 0, 0)) {
    if (!d_mpm_flags->d_withICE) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>linear</interpolator> \n"
          << " you should also use <extraCells>[0,0,0]</extraCells> \n"
          << " unless you are running an MPMICE case.\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  } else if (((interp_type == "gimp" || interp_type == "3rdorderBS" ||
               interp_type == "fast_cpdi" || interp_type == "cpti" ||
               interp_type == "cpdi") &&
              ((num_extra_cells + periodic) != IntVector(1, 1, 1) &&
               (!((num_extra_cells + periodic) == IntVector(1, 1, 0) &&
                  d_mpm_flags->d_axisymmetric))))) {
    std::ostringstream msg;
    msg << "\n ERROR: When using <interpolator>gimp</interpolator> \n"
        << " or <interpolator>3rdorderBS</interpolator> \n"
        << " or <interpolator>cpdi</interpolator> \n"
        << " or <interpolator>fast_cpdi</interpolator> \n"
        << " or <interpolator>cpti</interpolator> \n"
        << " you must also use extraCells and/or periodicBCs such\n"
        << " that the sum of the two is [1,1,1].\n"
        << " If using axisymmetry, the sum of the two can be [1,1,0].\n";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    // Initialize the accumulated strain energy
    new_dw->put(max_vartype(0.0), d_mpm_labels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_mpm_labels->partCountLabel);

  // For diffusion
  d_diffusionTasks->actuallyInitializeReductionVars(new_dw);

  cout_dbg << "Completed actuallyinitialize\n";
}

void
SerialMPM::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
  /*`==========TESTING==========*/
  Task* t = scinew Task(
    "SerialMPM::restartInitialize", this, &SerialMPM::restartInitialize);

  const PatchSet* patches = level->eachPatch();

  size_t numMPM = d_materialManager->getNumMaterials("MPM");
  for (size_t m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addReinitializeComputesAndRequires(t, mpm_matl, patches);
  }

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));
  /*`==========TESTING==========*/
}

/*! Task:  SerialMPM::restartInitialize
 *  Purpose:  Modify variables on a restart.  You MUST schedule a
 *            computes<label> in the restartInitializeHACK.
 */
void
SerialMPM::restartInitialize(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             [[maybe_unused]] DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    const std::string msg = "Doing SerialMPM::restartInitializeTask";
    printTask(patches, patch, cout_doing, msg);

    size_t numMatls = d_materialManager->getNumMaterials("MPM");

    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

      /*`==========TESTING==========*/
      cm->reinitializeCMData(patch, mpm_matl, new_dw);
      /*===========TESTING==========`*/
    }
  }
}

void
SerialMPM::scheduleDeleteGeometryObjects(const LevelP& level, SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();

  Task* t = scinew Task(
    "MPM::deleteGeometryObjects", this, &SerialMPM::deleteGeometryObjects);
  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
}

void
SerialMPM::deleteGeometryObjects(const ProcessorGroup*,
                                 [[maybe_unused]] const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse*,
                                 [[maybe_unused]] DataWarehouse* new_dw)
{
  printTask(cout_doing, "Doing MPM::deleteGeometryObjects");

  size_t numMPMMatls = d_materialManager->getNumMaterials("MPM");
  for (size_t m = 0; m < numMPMMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    std::cout << "MPM::Deleting Geometry Objects  matl: "
              << mpm_matl->getDWIndex() << "\n";
    mpm_matl->deleteGeomObjects();
  }

  // The call below is necessary because the GeometryPieceFactory holds on to a
  // pointer to all geom_pieces (so that it can look them up by name during
  // initialization) The pieces are never actually deleted until the factory is
  // destroyed at the end of the program. resetFactory() will rid of the pointer
  // (lookup table) and allow the deletion of the unneeded pieces.
  GeometryPieceFactory::resetFactory();
}

/*!----------------------------------------------------------------------
 * schedulePrintParticleCount
 *-----------------------------------------------------------------------*/
void
SerialMPM::schedulePrintParticleCount(const LevelP& level, SchedulerP& sched)
{
  Task* t = scinew Task(
    "MPM::printParticleCount", this, &SerialMPM::printParticleCount);
  t->needs(Task::NewDW, d_mpm_labels->partCountLabel);
  t->setType(Task::OncePerProc);
  sched->addTask(t,
                 sched->getLoadBalancer()->getPerProcessorPatchSet(level),
                 d_materialManager->allMaterials("MPM"));
}

/*!----------------------------------------------------------------------
 * printParticleCount
 *-----------------------------------------------------------------------*/
void
SerialMPM::printParticleCount(const ProcessorGroup* pg,
                              const PatchSubset*,
                              const MaterialSubset*,
                              DataWarehouse*,
                              DataWarehouse* new_dw)
{
  sumlong_vartype pcount;
  new_dw->get(pcount, d_mpm_labels->partCountLabel);

  if (pg->myRank() == 0) {
    std::cerr << "Created " << (long)pcount << " total particles\n";
  }
}

//  Diagnostic task: compute the total number of particles
void
SerialMPM::scheduleTotalParticleCount(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task(
    "SerialMPM::totalParticleCount", this, &SerialMPM::totalParticleCount);
  t->computes(d_mpm_labels->partCountLabel);
  sched->addTask(t, patches, matls);
}

//  Diagnostic task: compute the total number of particles
void
SerialMPM::totalParticleCount(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    long int totalParticles = 0;
    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
      int numParticles     = pset->end() - pset->begin();

      totalParticles += numParticles;
    }
    new_dw->put(sumlong_vartype(totalParticles), d_mpm_labels->partCountLabel);
  }
}

/*!------------------------------------------------------------------------
 * Schedule the initialization of the stress and deformation gradient
 * based on the body forces (which also have to be computed)
 *------------------------------------------------------------------------*/
void
SerialMPM::scheduleInitializeStressAndDefGradFromBodyForce(const LevelP& level,
                                                           SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();
  printSchedule(
    patches, cout_doing, "MPM::initializeStressAndDefGradFromBodyForce");

  // First compute the body force
  Task* t1 = scinew Task(
    "MPM::initializeBodyForce", this, &SerialMPM::initializeBodyForce);
  t1->needs(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
  t1->modifies(d_mpm_labels->pBodyForceAccLabel);
  sched->addTask(t1, patches, d_materialManager->allMaterials("MPM"));

  // Compute the stress and deformation gradient only for selected
  // constitutive models that have a "initializeWithBodyForce" flag as true.
  // This is because a more general implementation is quite involved and
  // not worth the effort at this time. BB
  Task* t2 = scinew Task("MPM::initializeStressAndDefGradFromBodyForce",
                         this,
                         &SerialMPM::initializeStressAndDefGradFromBodyForce);

  t2->needs(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
  t2->needs(Task::NewDW, d_mpm_labels->pBodyForceAccLabel, Ghost::None);
  t2->modifies(d_mpm_labels->pStressLabel);
  t2->modifies(d_mpm_labels->pDefGradLabel);
  sched->addTask(t2, patches, d_materialManager->allMaterials("MPM"));
}

/*!------------------------------------------------------------------------
 * Actually initialize the body force acceleration
 *-------------------------------------------------------------------------*/
void
SerialMPM::initializeBodyForce(const ProcessorGroup*,
                               const PatchSubset* patches,
                               [[maybe_unused]] const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse* new_dw)
{
  // Get the MPM d_mpm_flags and make local copies
  Uintah::Point rotation_center = d_mpm_flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis  = d_mpm_flags->d_coordRotationAxis;
  double rotation_speed         = d_mpm_flags->d_coordRotationSpeed;

  // Compute angular velocity std::vector (omega)
  Uintah::Vector omega = rotation_axis * rotation_speed;

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing initializeBodyForce");

    // Loop thru materials
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {

      // Get the material ID
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      new_dw->getModifiable(
        pBodyForceAcc, d_mpm_labels->pBodyForceAccLabel, pset);

      // Get the position data
      constParticleVariable<Point> pPosition;
      new_dw->get(pPosition, d_mpm_labels->pXLabel, pset);

      // Iterate over the particles
      for (auto pidx : *pset) {

        // Compute the body force acceleration (g)
        // Just use gravity if rotation is off
        pBodyForceAcc[pidx] = d_mpm_flags->d_gravity;

        // If rotating add centrifugal force
        if (d_mpm_flags->d_useCoordRotation) {

          // Compute the centrifugal term (omega x omega x r)
          // Simplified version where body ref point is not needed
          Vector rVec              = pPosition[pidx] - rotation_center;
          Vector omega_x_r         = Uintah::Cross(omega, rVec);
          Vector centrifugal_accel = Uintah::Cross(omega, omega_x_r);

          // Compute the body force acceleration (g - omega x omega x r)
          pBodyForceAcc[pidx] -= centrifugal_accel;
        } // coord rotation end if

      } // end particle loop
    }   // end matl loop
  }     // end patch loop
}

/*!--------------------------------------------------------------------------
 * Actually initialize the stress and deformation gradient assuming linear
 * elastic behavior after computing the body force acceleration
 *
 * **WARNING** Assumes zero shear stresses and that body forces are aligned
 *             with coordinate directions
 *--------------------------------------------------------------------------*/
void
SerialMPM::initializeStressAndDefGradFromBodyForce(const ProcessorGroup*,
                                                   const PatchSubset* patches,
                                                   const MaterialSubset* matls,
                                                   DataWarehouse*,
                                                   DataWarehouse* new_dw)
{
  // Loop over patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              cout_doing,
              "Doing initializeStressAndDefGradFromBodyForce");

    // Loop over materials
    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

      // Compute the stress and deformation gradient only for selected
      // constitutive models that have a "initializeWithBodyForce" flag as true.
      // A more general implementation is not worth the significant extra
      // effort. BB
      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
      cm->initializeStressAndDefGradFromBodyForce(patch, mpm_matl, new_dw);

    } // end matl loop

  } // end patches loop
}

/*!--------------------------------------------------------------------------
 * Schedule the initialization of the external forces: Pressure
 *---------------------------------------------------------------------------*/
void
SerialMPM::scheduleInitializePressureBCs(const LevelP& level, SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();

  d_loadCurveIndex = scinew MaterialSubset();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int pressureBCId = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    std::string bcType = bc->getType();
    if (bcType == "Pressure") {
      d_loadCurveIndex->add(pressureBCId++);
    }
  }
  if (pressureBCId > 0) {
    printSchedule(patches, cout_doing, "MPM::countMaterialPointsPerLoadCurve");
    printSchedule(patches, cout_doing, "MPM::scheduleInitializePressureBCs");
    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task("MPM::countMaterialPointsPerLoadCurve",
                          this,
                          &SerialMPM::countMaterialPointsPerLoadCurve);
    t->needs(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task(
      "MPM::initializePressureBC", this, &SerialMPM::initializePressureBC);
    t->needs(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pSizeLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pDispLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->needs(Task::NewDW,
                d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpm_labels->pExternalForceLabel);
    if (d_mpm_flags->d_useCBDI) {
      t->computes(d_mpm_labels->pExternalForceCorner1Label);
      t->computes(d_mpm_labels->pExternalForceCorner2Label);
      t->computes(d_mpm_labels->pExternalForceCorner3Label);
      t->computes(d_mpm_labels->pExternalForceCorner4Label);
    }
    sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
  }

  if (d_loadCurveIndex->removeReference()) {
    delete d_loadCurveIndex;
  }
}

/*!----------------------------------------------------------------------
 * countMaterialPointsPerLoadCurve
 *   Calculate the number of material points per load curve
 *-----------------------------------------------------------------------*/
void
SerialMPM::countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset*,
                                           DataWarehouse*,
                                           DataWarehouse* new_dw)
{
  printTask(
    patches, patches->get(0), cout_doing, "countMaterialPointsPerLoadCurve");
  // Find the number of pressure BCs in the problem
  int nofPressureBCs = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    std::string bcType = bc->getType();
    if (bcType == "Pressure" || bcType == "Moment") {
      nofPressureBCs++;

      // Loop through the patches and count
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        int numPts         = 0;
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            if (pLoadCurveID[idx] == (nofPressureBCs)) {
              ++numPts;
            }
          }
        } // matl loop
        new_dw->put(sumlong_vartype(numPts),
                    d_mpm_labels->materialPointsPerLoadCurveLabel,
                    nullptr,
                    nofPressureBCs - 1);
      } // patch loop
    }
  }
}

/*!----------------------------------------------------------------------
 * initializePressureBC
 *-----------------------------------------------------------------------*/
void
SerialMPM::initializePressureBC(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse*,
                                DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;
  printTask(patches, patches->get(0), cout_doing, "Doing initializePressureBC");
  if (cout_dbg.active()) {
    cout_dbg << "Current Time (Initialize Pressure BC) = " << time << "\n";
  }

  // Calculate the force std::vector at each particle
  int pressureBCId = 0;
  int ii           = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {

    std::string bcType = bc->getType();

    if (bcType == "Pressure") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart,
                  d_mpm_labels->materialPointsPerLoadCurveLabel,
                  nullptr,
                  pressureBCId++);

      // Save the material points per load curve in the PressureBC object
      auto* pbc = dynamic_cast<PressureBC*>(bc.get());
      pbc->numMaterialPoints(numPart);

      if (cout_dbg.active()) {
        cout_dbg << "    Load Curve = " << pressureBCId
                 << " Num Particles = " << numPart << "\n";
      }

      // Calculate the force per particle at t = 0.0
      double forcePerPart = pbc->forcePerParticle(time);

      // Loop through the patches and calculate the force std::vector
      // at each particle
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

          constParticleVariable<Point> pX;
          constParticleVariable<Matrix3> pSize;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX, d_mpm_labels->pXLabel, pset);
          new_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
          new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);
          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(
            pExternalForce, d_mpm_labels->pExternalForceLabel, pset);

          constParticleVariable<Vector> pDisp;
          new_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

          ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
            pExternalForceCorner3, pExternalForceCorner4;
          if (d_mpm_flags->d_useCBDI) {
            if (ii == 0) {
              new_dw->allocateAndPut(pExternalForceCorner1,
                                     d_mpm_labels->pExternalForceCorner1Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner2,
                                     d_mpm_labels->pExternalForceCorner2Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner3,
                                     d_mpm_labels->pExternalForceCorner3Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner4,
                                     d_mpm_labels->pExternalForceCorner4Label,
                                     pset);
            } else {
              new_dw->getModifiable(pExternalForceCorner1,
                                    d_mpm_labels->pExternalForceCorner1Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner2,
                                    d_mpm_labels->pExternalForceCorner2Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner3,
                                    d_mpm_labels->pExternalForceCorner3Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner4,
                                    d_mpm_labels->pExternalForceCorner4Label,
                                    pset);
            }
          }

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == pressureBCId) {
              if (d_mpm_flags->d_useCBDI) {
                Vector dxCell = patch->dCell();
                pExternalForce[idx] =
                  pbc->getForceVectorCBDI(pX[idx],
                                          pDisp[idx],
                                          pSize[idx],
                                          pDefGrad[idx],
                                          forcePerPart,
                                          time,
                                          pExternalForceCorner1[idx],
                                          pExternalForceCorner2[idx],
                                          pExternalForceCorner3[idx],
                                          pExternalForceCorner4[idx],
                                          dxCell);
              } else {
                pExternalForce[idx] = pbc->getForceVector(
                  pX[idx], pDisp[idx], forcePerPart, time, pDefGrad[idx]);
              }
            }
          }
        } // matl loop
      }   // patch loop
    }
    ++ii;
  } // bc loop
}

/*!---------------------------------------------------------------------------
 * scheduleInitializeMomentBCs
 *   Schedule the initialization of the external forces: Moments
 *---------------------------------------------------------------------------*/
void
SerialMPM::scheduleInitializeMomentBCs(const LevelP& level, SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();

  d_loadCurveIndex = scinew MaterialSubset();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int nofMomentBCs = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    std::string bcType = bc->getType();
    if (bcType == "Moment") {
      d_loadCurveIndex->add(nofMomentBCs++);
    }
  }
  if (nofMomentBCs > 0) {
    printSchedule(patches, cout_doing, "MPM::countMaterialPointsPerLoadCurve");
    printSchedule(patches, cout_doing, "MPM::scheduleInitializeMomentBCs");
    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task("MPM::countMaterialPointsPerLoadCurve",
                          this,
                          &SerialMPM::countMaterialPointsPerLoadCurve);
    t->needs(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));

    // Create a task that calculates the force to be associated with
    // each particle based on the moment BCs
    t = scinew Task(
      "MPM::initializeMomentBC", this, &SerialMPM::initializeMomentBC);
    t->needs(Task::NewDW, d_mpm_labels->pXLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pSizeLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel, Ghost::None);
    t->needs(Task::NewDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->needs(Task::NewDW,
                d_mpm_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpm_labels->pExternalForceLabel);
    if (d_mpm_flags->d_useCBDI) {
      t->computes(d_mpm_labels->pExternalForceCorner1Label);
      t->computes(d_mpm_labels->pExternalForceCorner2Label);
      t->computes(d_mpm_labels->pExternalForceCorner3Label);
      t->computes(d_mpm_labels->pExternalForceCorner4Label);
    }
    sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
  }

  if (d_loadCurveIndex->removeReference()) {
    delete d_loadCurveIndex;
  }
}

/*!----------------------------------------------------------------------
 * initializeMomentBC
 *-----------------------------------------------------------------------*/
void
SerialMPM::initializeMomentBC(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse*,
                              DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;
  printTask(patches, patches->get(0), cout_doing, "Doing initializeMomentBC");
  if (cout_dbg.active()) {
    cout_dbg << "Current Time (Initialize Moment BC) = " << time << "\n";
  }

  // Calculate the force std::vector at each particle
  int nofMomentBCs = 0;
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    std::string bcType = bc->getType();
    if (bcType == "Moment") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart,
                  d_mpm_labels->materialPointsPerLoadCurveLabel,
                  nullptr,
                  nofMomentBCs++);

      // Save the material points per load curve in the MomentBC object
      auto* pbc = dynamic_cast<MomentBC*>(bc.get());
      pbc->numMaterialPoints(numPart);

      if (cout_dbg.active()) {
        cout_dbg << "    Load Curve = " << nofMomentBCs
                 << " Num Particles = " << numPart << "\n";
      }

      // Calculate the force per particle at t = 0.0
      double forcePerPart = pbc->forcePerParticle(time);

      // Loop through the patches and calculate the force std::vector
      // at each particle
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        int numMPMMatls    = d_materialManager->getNumMaterials("MPM");
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl =
            static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<Point> pX;
          constParticleVariable<Matrix3> pSize;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX, d_mpm_labels->pXLabel, pset);
          new_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
          new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);
          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(
            pExternalForce, d_mpm_labels->pExternalForceLabel, pset);

          ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
            pExternalForceCorner3, pExternalForceCorner4;
          if (d_mpm_flags->d_useCBDI) {
            new_dw->allocateAndPut(pExternalForceCorner1,
                                   d_mpm_labels->pExternalForceCorner1Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner2,
                                   d_mpm_labels->pExternalForceCorner2Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner3,
                                   d_mpm_labels->pExternalForceCorner3Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner4,
                                   d_mpm_labels->pExternalForceCorner4Label,
                                   pset);
          }
          // std::cout << "d_mpm_flags->d_useCBDI: " << d_mpm_flags->d_useCBDI
          // <<
          // "\n";
          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            if (pLoadCurveID[idx] == nofMomentBCs) {
              if (d_mpm_flags->d_useCBDI) {
                Vector dxCell = patch->dCell();
                pExternalForce[idx] =
                  pbc->getForceVectorCBDI(pX[idx],
                                          pSize[idx],
                                          pDefGrad[idx],
                                          forcePerPart,
                                          time,
                                          pExternalForceCorner1[idx],
                                          pExternalForceCorner2[idx],
                                          pExternalForceCorner3[idx],
                                          pExternalForceCorner4[idx],
                                          dxCell);
              } else {
                pExternalForce[idx] = pbc->getForceVector(
                  pX[idx], forcePerPart, time, pDefGrad[idx]);
              }
            }
          }
        } // matl loop
      }   // patch loop
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeStableTimsetp
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  // Nothing to do here - delt is computed as a by-product of the
  // constitutive model
  // However, this task needs to do something in the case that MPM
  // is being run on more than one level.
  cout_doing << UintahParallelComponent::d_myworld->myRank()
             << " MPM::scheduleComputeStableTimestep \t\t\t\tL-"
             << level->getIndex() << "\n";

  Task* t = scinew Task("MPM::actuallyComputeStableTimestep",
                        this,
                        &SerialMPM::actuallyComputeStableTimestep);

  if (d_mpm_flags->d_enableDEM) {
    t->needs(Task::NewDW, d_mpm_labels->pMassLabel, Ghost::None);
  }

  t->computes(d_mpm_labels->delTLabel, level.get_rep());
  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));
}

/*!----------------------------------------------------------------------
 * actuallyComputeStableTimestep
 *-----------------------------------------------------------------------*/
void
SerialMPM::actuallyComputeStableTimestep(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset*,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{
  // The standard MPM CFL is usually computed in the constitutive models.
  // Here we compute the DEM stability limit: dt < sqrt(m/k)
  double dt_dem = 1.0e30;

  int numMPMMatls = d_materialManager->getNumMaterials("MPM");

  if (!old_dw) {
    // Initial timestep: use a conservative estimate based on minPartMass
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      d_demTasks->computeStableTimestep(nullptr, mpm_matl, old_dw, new_dw, dt_dem);
    }
    const Level* level = getLevel(patches);
    new_dw->put(delt_vartype(dt_dem), d_mpm_labels->delTLabel, level);
    DOUT(mpm_delt_dbg, "Init MPM Discrete/Rigid delT = " << dt_dem );
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      d_demTasks->computeStableTimestep(patch, mpm_matl, old_dw, new_dw, dt_dem);
    }
  }

  DOUT(mpm_delt_dbg, "[MPM:ComputeStableTimestep] DEM/Rigid delt = " << dt_dem);

  const Level* level = getLevel(patches);
  new_dw->put(delt_vartype(dt_dem), d_mpm_labels->delTLabel, level);
}

/*!----------------------------------------------------------------------
 * scheduleTimeAdvance
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  MALLOC_TRACE_TAG_SCOPE("SerialMPM::scheduleTimeAdvance()");
  if (!d_mpm_flags->doMPMOnLevel(level->getIndex(),
                                 level->getGrid()->numLevels())) {
    return;
  }

  const PatchSet* patches  = level->eachPatch();
  const MaterialSet* matls = d_materialManager->allMaterials("MPM");

  // Compute body forces first
  scheduleComputeParticleBodyForce(sched, patches, matls);

  // Compute particle size for improved interpolators
  scheduleComputeCurrentParticleSize(sched, patches, matls);

  // Apply mechanical loads (PhysicalBCs)
  scheduleApplyExternalLoads(sched, patches, matls);

  // For discrete element modeling
  if (d_mpm_flags->d_enableDEM) {
    d_demTasks->scheduleComputeDEMForces(sched, patches, matls);
  }

  // For scalar diffusion
  d_diffusionTasks->scheduleApplyExternalScalarFlux(sched, patches, matls);

  // Move data to grid
  scheduleInterpolateParticlesToGrid(sched, patches, matls);

  // For contact
  scheduleComputeNormals(sched, patches, matls);
  scheduleFindSurfaceParticles(sched, patches, matls);
  scheduleComputeLogisticRegression(sched, patches, matls);

  // Exchange momentum on grid
  scheduleMomentumExchangeInterpolated(sched, patches, matls);

  // For scalar diffusion
  d_diffusionTasks->scheduleConcInterpolated(sched, patches, matls);

  // For cohesive zones
  d_cohesiveZoneTasks->scheduleUpdate(sched, patches, matls);

  scheduleComputeContactArea(sched, patches, matls);

  scheduleComputeInternalForce(sched, patches, matls);

  // For scalar diffusion
  d_diffusionTasks->scheduleCompute(sched, patches, matls);

  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);

  // For discrete element modeling
  if (d_mpm_flags->d_enableDEM) {
    d_demTasks->scheduleIntegrateDEMRotation(sched, patches, matls);
  }

  // For scalar diffusion
  d_diffusionTasks->scheduleIntegrate(sched, patches, matls);

  scheduleMomentumExchangeIntegrated(sched, patches, matls);

  scheduleSetGridBoundaryConditions(sched, patches, matls);

  scheduleSetPrescribedMotion(sched, patches, matls);

  // For XPIC(2) computations
  if (d_mpm_flags->d_useXPIC) {
    scheduleComputeXPICVelocities(sched, patches, matls);
  }

  // For explicit heat conduction
  d_heatConductionTasks->scheduleCompute(sched, patches, matls);

  // Update stress last
  if (d_mpm_flags->d_updateStressLast) {
    scheduleReduceVars(sched, patches, matls);
    scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  }

  // Schedule compute of the deformation gradient
  scheduleComputeDeformationGradient(sched, patches, matls);

  // Schedule compute of the stress tensor
  scheduleComputeStressTensor(sched, patches, matls);

  // Create a task for computing damage and updating stress
  scheduleComputeBasicDamage(sched, patches, matls);

  // Schedule update of the erosion parameter
  scheduleUpdateErosionParameter(sched, patches, matls);

  // Schedule task to find rogue particles
  scheduleFindRogueParticles(sched, patches, matls);

  // Schedule task to compute the accumulated strain energy
  if (d_mpm_flags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }

  // Update stress first
  if (!d_mpm_flags->d_updateStressLast) {
    scheduleReduceVars(sched, patches, matls);
    scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  }

  // Final update of particles
  // This is from the Uintah implementation and may be separated out
  // from interpolateToParticlesAndUpdate at some stage
  // scheduleFinalParticleUpdate(sched, patches, matls);

  scheduleInsertParticles(sched, patches, matls);

  if (d_mpm_flags->d_refineParticles) {
    scheduleAddParticles(sched, patches, matls);
  }
  if (d_mpm_flags->d_computeScaleFactor) {
    scheduleComputeParticleScaleFactor(sched, patches, matls);
  }

  if (d_analysisModules.size() != 0) {
    for (auto& am : d_analysisModules) {
      am->scheduleDoAnalysis_preReloc(sched, level);
    }
  }

  scheduleParticleRelocation(sched, level, patches, matls);

  //__________________________________
  //  on the fly analysis
  if (d_analysisModules.size() != 0) {
    for (auto& module_name : d_analysisModules) {
      module_name->scheduleDoAnalysis(sched, level);
    }
  }
}

/*!====================================================================================
 * Method: scheduleComputeParticleBodyForce
 * Purpose: Schedule a task to compute particle body forces
 * Inputs:  p.x
 * Outputs: p.bodyForce
 *====================================================================================*/
void
SerialMPM::scheduleComputeParticleBodyForce(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleBodyForce");

  Task* t = scinew Task("MPM::computeParticleBodyForce",
                        this,
                        &SerialMPM::computeParticleBodyForce);

  t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, Ghost::None);
  // t->computes(d_mpm_labels->pBodyForceAccLabel);
  // t->computes(d_mpm_labels->pCoriolisImportanceLabel);
  t->computes(d_mpm_labels->pBodyForceAccLabel_preReloc);
  t->computes(d_mpm_labels->pCoriolisImportanceLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*!====================================================================================
 * Method: computeParticleBodyForce
 * Purpose: Actually compute particle body forces
 * Inputs:  p.x
 * Outputs: p.bodyForce
 *====================================================================================*/
void
SerialMPM::computeParticleBodyForce(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  // Get the MPM d_mpm_flags and make local copies
  Uintah::Point rotation_center = d_mpm_flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis  = d_mpm_flags->d_coordRotationAxis;
  double rotation_speed         = d_mpm_flags->d_coordRotationSpeed;
  // Uintah::Point body_ref_point =
  // d_mpm_flags->d_coord_rotation_body_ref_point;

  // Compute angular velocity std::vector (omega)
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
      // new_dw->allocateAndPut(pBodyForceAcc, d_mpm_labels->pBodyForceAccLabel,
      // pset);
      new_dw->allocateAndPut(
        pBodyForceAcc, d_mpm_labels->pBodyForceAccLabel_preReloc, pset);

      // Create space for particle coriolis importance
      ParticleVariable<double> pCoriolisImportance;
      // new_dw->allocateAndPut(pCoriolisImportance,
      // d_mpm_labels->pCoriolisImportanceLabel, pset);
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

        // Iterate over the particles
        // std::cout << "Mat id = " << matID << " patch = " << patch << "\n";
        // std::cout << "Particle subset = " << *pset;
        // std::cout << "Num particles = " << pset->numParticles() << "\n";
        for (int pidx : *pset) {
          // std::cout << " Particle # = " << pidx << "\n";
          //  Compute the local "x" std::vector wrt ref point in body
          // Vector xVec = pPosition[pidx].std::vector() - body_ref_point;

          // Compute reference std::vector R wrt rotation center
          // Uintah::Vector Rvec = body_ref_point - rotation_center;

          // Compute the local "r" std::vector with respect to rotation center
          // Vector rVec = Rvec + pPosition[pidx].std::vector();

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

          /*
          //if (pVelocity[pidx].length2() > 0.0) {
          if (pCoriolisImportance[pidx] > 0.7) {
          std::cout << "pidx = " << pidx << " omega = " << omega << " x = " <<
          pPosition[pidx] << " r = " << rVec << " v = " << pVelocity[pidx] <<
          "\n"; std::cout << "\t omega x r = " << omega_x_r << " omega x omega x
          r = " << centrifugal_accel << " omega x v = " << coriolis_accel <<
          "\n" ; std::cout << "\t b = " << pBodyForceAcc[pidx]
          << " cor. imp. = " << pCoriolisImportance[pidx] << "\n";
          }
          */
        } // particle loop
      }   // end if coordinate rotation

      // Copy data for relocation if particles cross patch boundaries
      /*
        ParticleVariable<double> pCoriolisImportance_new;
        new_dw->allocateAndPut(pCoriolisImportance_new,
        d_mpm_labels->pCoriolisImportanceLabel_preReloc, pset);
        pCoriolisImportance_new.copyData(pCoriolisImportance);

        ParticleVariable<Vector> pBodyForceAcc_new;
        new_dw->allocateAndPut(pBodyForceAcc_new,
        d_mpm_labels->pBodyForceAccLabel_preReloc, pset);
        pBodyForceAcc_new.copyData(pBodyForceAcc);
      */

    } // matl loop
  }   // patch loop
}

void
SerialMPM::scheduleComputeCurrentParticleSize(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeCurrentParticleSize");

  Task* t = scinew Task("MPM::computeCurrentParticleSize",
                        this,
                        &SerialMPM::computeCurrentParticleSize);

  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None);

  t->computes(d_mpm_labels->pCurSizeLabel);

  sched->addTask(t, patches, matls);
}

void
SerialMPM::computeCurrentParticleSize(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(
      patches, patch, cout_doing, "Doing MPM::computeCurrentParticleSize");

    size_t numMatls = d_materialManager->getNumMaterials("MPM");

    std::string interp_type = d_mpm_flags->d_interpolatorType;

    for (size_t m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad_old;
      ParticleVariable<Matrix3> pCurSize;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_mpm_labels->pDefGradLabel, pset);
      new_dw->allocateAndPut(pCurSize, d_mpm_labels->pCurSizeLabel, pset);

      if (interp_type == "cpdi" || interp_type == "fast_cpdi" ||
          interp_type == "cpti") {
        if (d_mpm_flags->d_axisymmetric) {
          for (auto idx : *pset) {
            Matrix3 defgrad1 = Matrix3(pDefGrad_old[idx](0, 0),
                                       pDefGrad_old[idx](0, 1),
                                       0.0,
                                       pDefGrad_old[idx](1, 0),
                                       pDefGrad_old[idx](1, 1),
                                       0.0,
                                       0.0,
                                       0.0,
                                       1.0);

            pCurSize[idx] = defgrad1 * pSize[idx];
          }
        } else {
          for (auto idx : *pset) {
            pCurSize[idx] = pDefGrad_old[idx] * pSize[idx];
          }
        }
      } else {
        pCurSize.copyData(pSize);
      }
    }
  }
}

/*====================================================================================*/
// Apply external loads
//*
//* applyExternalLoads
//*   in(p.externalForce)
//*   out(p.externalForceNew) */
/*====================================================================================*/
void
SerialMPM::scheduleApplyExternalLoads(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleApplyExternalLoads");

  Task* t = scinew Task(
    "MPM::applyExternalLoads", this, &SerialMPM::applyExternalLoads);

  t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pDispLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None);
  t->needs(Task::OldDW, d_mpm_labels->pExternalForceLabel, Ghost::None);
  t->computes(d_mpm_labels->pExtForceLabel_preReloc);
  if (d_mpm_flags->d_useLoadCurves) {
    t->needs(Task::OldDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpm_labels->pLoadCurveIDLabel_preReloc);
    if (d_mpm_flags->d_useCBDI) {
      t->computes(d_mpm_labels->pExternalForceCorner1Label);
      t->computes(d_mpm_labels->pExternalForceCorner2Label);
      t->computes(d_mpm_labels->pExternalForceCorner3Label);
      t->computes(d_mpm_labels->pExternalForceCorner4Label);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * addExternalLoads
 *-----------------------------------------------------------------------*/
void
SerialMPM::applyExternalLoads(const ProcessorGroup*,
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
    cout_doing << "Current Time (applyExternalLoads) = " << time << "\n";
  }

  // Calculate the force std::vector at each particle for each pressure bc
  std::vector<double> forcePerPart;
  std::vector<PressureBC*> pbcP;
  std::vector<MomentBC*> pbcM;
  if (d_mpm_flags->d_useLoadCurves) {
    for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
      std::string bcType = bc->getType();
      if (bcType == "Pressure") {
        auto* pbc = dynamic_cast<PressureBC*>(bc.get());
        pbcP.push_back(pbc);

        // Calculate the force per particle at current time
        forcePerPart.push_back(pbc->forcePerParticle(time));
      } else if (bcType == "Moment") {
        auto* pbc = dynamic_cast<MomentBC*>(bc.get());
        pbcM.push_back(pbc);

        // Calculate the moment at current time.
        forcePerPart.push_back(pbc->forcePerParticle(time));
      }
    }
  }

  // Loop thru patches to update external force std::vector
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing applyExternalLoads");

    // Place for user defined loading scenarios to be defined,
    // otherwise pExternalForce is just carried forward.

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Get the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;

      old_dw->get(pX, d_mpm_labels->pXLabel, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, pset);

      ParticleVariable<Vector> pExternalForce_new;
      ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
        pExternalForceCorner3, pExternalForceCorner4;

      new_dw->allocateAndPut(
        pExternalForce_new, d_mpm_labels->pExtForceLabel_preReloc, pset);
      if (d_mpm_flags->d_useCBDI) {
        new_dw->allocateAndPut(pExternalForceCorner1,
                               d_mpm_labels->pExternalForceCorner1Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner2,
                               d_mpm_labels->pExternalForceCorner2Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner3,
                               d_mpm_labels->pExternalForceCorner3Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner4,
                               d_mpm_labels->pExternalForceCorner4Label,
                               pset);
      }
      // std::cout << "applyloads: patch = " << patch << " matID = " << matID
      //           << " numparticles = " << pset->numParticles() << "
      //           d_numGhostParticles = " << d_numGhostParticles << "\n";

      if (d_mpm_flags->d_useLoadCurves) {

        // Get the load curve data
        // Recycle the loadCurveIDs
        constParticleVariable<int> pLoadCurveID;
        old_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

        ParticleVariable<int> pLoadCurveID_new;
        new_dw->allocateAndPut(
          pLoadCurveID_new, d_mpm_labels->pLoadCurveIDLabel_preReloc, pset);
        pLoadCurveID_new.copyData(pLoadCurveID);
        // std::cout << " Recycled load curve ID" << "\n";

        // Check whether it's a presure or moment bc
        bool do_PressureBCs = false;
        bool do_MomentBCs   = false;
        for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          std::string bcType = bc->getType();
          if (bcType == "Pressure") {
            do_PressureBCs = true;
          } else if (bcType == "Moment") {
            do_MomentBCs = true;
          }
        }

        if (do_PressureBCs) {

          // Get the external force data and allocate new space for
          // external force
          constParticleVariable<Vector> pDisp;
          old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

          // Iterate over the particles
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = Vector(0.0, 0.0, 0.0);
              if (d_mpm_flags->d_useCBDI) {
                pExternalForceCorner1[idx] = pX[idx];
                pExternalForceCorner2[idx] = pX[idx];
                pExternalForceCorner3[idx] = pX[idx];
                pExternalForceCorner4[idx] = pX[idx];
              }
            } else {
              PressureBC* pbc = pbcP[loadCurveID];
              double force    = forcePerPart[loadCurveID];

              if (d_mpm_flags->d_useCBDI) {
                Vector dxCell = patch->dCell();
                pExternalForce_new[idx] =
                  pbc->getForceVectorCBDI(pX[idx],
                                          pDisp[idx],
                                          pSize[idx],
                                          pDefGrad[idx],
                                          force,
                                          time,
                                          pExternalForceCorner1[idx],
                                          pExternalForceCorner2[idx],
                                          pExternalForceCorner3[idx],
                                          pExternalForceCorner4[idx],
                                          dxCell);
                /*
                std::cout << "idx = " << idx << "PX = " << pX[idx] << " fext = "
                << pExternalForce_new[idx] << "\n"; std::cout << "corners: \n"
                          << pExternalForceCorner1[idx] << ", "
                          << pExternalForceCorner2[idx] << ", "
                          << pExternalForceCorner3[idx] << ", "
                          << pExternalForceCorner4[idx] << "\n";
                */
              } else {
                pExternalForce_new[idx] = pbc->getForceVector(
                  pX[idx], pDisp[idx], force, time, pDefGrad[idx]);
              }
            }
          }
        } else if (do_MomentBCs) {
          // Get the external force data and allocate new space for
          // external force
          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_mpm_labels->pExternalForceLabel, pset);

          // Iterate over the particles
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = pExternalForce[idx];
            } else {
              MomentBC* pbc = pbcM[loadCurveID];
              double force  = forcePerPart[loadCurveID];

              pExternalForce_new[idx] =
                pbc->getForceVector(pX[idx], force, time, pDefGrad[idx]);
            }
          }
        } else {
          for (auto idx : *pset) {
            pExternalForce_new[idx] = 0.;
          }
        }

        // MMS (compute body force)
        std::string mms_type = d_mpm_flags->d_mmsType;
        if (!mms_type.empty()) {
          MMS MMSObject;
          MMSObject.computeBodyForceForMMS(old_dw,
                                           new_dw,
                                           time,
                                           pset,
                                           d_mpm_labels.get(),
                                           d_mpm_flags.get(),
                                           pExternalForce_new);
        }

      } else { // d_useLoadCurves = False
        // MMS
        std::string mms_type = d_mpm_flags->d_mmsType;
        if (!mms_type.empty()) {
          MMS MMSObject;
          MMSObject.computeExternalForceForMMS(old_dw,
                                               new_dw,
                                               time,
                                               pset,
                                               d_mpm_labels.get(),
                                               d_mpm_flags.get(),
                                               pExternalForce_new);
        } else {
          // Get the external force data and allocate new space for
          // external force and copy the data
          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_mpm_labels->pExternalForceLabel, pset);

          for (auto idx : *pset) {
            pExternalForce_new[idx] =
              pExternalForce[idx] * d_mpm_flags->d_forceIncrementFactor;
          }
        }
      } // end if (d_useLoadCurves)
    }   // matl loop
  }     // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleMomentumExchangeIntegrated
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleMomentumExchangeIntegrated(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* exMomIntegrated
   *   in(G.MASS, G.VELOCITY_STAR, G.ACCELERATION)
   *   operation(peform operations which will cause each of
   *              velocity fields to feel the influence of the
   *              the others according to specific rules)
   *   out(G.VELOCITY_STAR, G.ACCELERATION) */
  printSchedule(patches, cout_doing, "MPM::scheduleMomentumExchangeIntegrated");
  contactModel->addComputesAndRequires(
    sched, patches, matls, d_mpm_labels->gVelocityStarLabel);
}

/*!----------------------------------------------------------------------
 * scheduleSetGridBoundaryConditions
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleSetGridBoundaryConditions");
  Task* t = scinew Task("MPM::setGridBoundaryConditions",
                        this,
                        &SerialMPM::setGridBoundaryConditions);

  const MaterialSubset* mss = matls->getUnion();
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  t->modifies(d_mpm_labels->gAccelerationLabel, mss);
  t->modifies(d_mpm_labels->gVelocityStarLabel, mss);
  t->needs(Task::NewDW, d_mpm_labels->gVelocityLabel, Ghost::None);

  if (!d_mpm_flags->d_doGridReset) {
    t->needs(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
    t->computes(d_mpm_labels->gDisplacementLabel);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * setGridBoundaryConditions
 *-----------------------------------------------------------------------*/
void
SerialMPM::setGridBoundaryConditions(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing setGridBoundaryConditions");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    std::string interp_type = d_mpm_flags->d_interpolatorType;
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gAcceleration;
      constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(
        gAcceleration, d_mpm_labels->gAccelerationLabel, matID, patch);
      new_dw->getModifiable(
        gVelocity_star, d_mpm_labels->gVelocityStarLabel, matID, patch);
      new_dw->get(
        gVelocity, d_mpm_labels->gVelocityLabel, matID, patch, Ghost::None, 0);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      MPMBoundCond bc;
      bc.setBoundaryCondition(
        patch, matID, "Velocity", gVelocity_star, interp_type);
      bc.setBoundaryCondition(
        patch, matID, "Symmetric", gVelocity_star, interp_type);

      // Now recompute acceleration as the difference between the velocity
      // interpolated to the grid (no bcs applied) and the new velocity_star
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c      = *iter;
        gAcceleration[c] = (gVelocity_star[c] - gVelocity[c]) / delT;
      }

      if (!d_mpm_flags->d_doGridReset) {
        NCVariable<Vector> displacement;
        constNCVariable<Vector> displacementOld;
        new_dw->allocateAndPut(
          displacement, d_mpm_labels->gDisplacementLabel, matID, patch);
        old_dw->get(displacementOld,
                    d_mpm_labels->gDisplacementLabel,
                    matID,
                    patch,
                    Ghost::None,
                    0);
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c     = *iter;
          displacement[c] = displacementOld[c] + gVelocity_star[c] * delT;
        }
      } // d_doGridReset
    }   // matl loop
  }     // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleSetPrescribedMotion
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleSetPrescribedMotion(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpm_flags->d_prescribeDeformation) {
    printSchedule(patches, cout_doing, "MPM::scheduleSetPrescribedMotion");

    Task* t = scinew Task(
      "MPM::setPrescribedMotion", this, &SerialMPM::setPrescribedMotion);

    const MaterialSubset* mss = matls->getUnion();
    t->modifies(d_mpm_labels->gAccelerationLabel, mss);
    t->modifies(d_mpm_labels->gVelocityStarLabel, mss);
    t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
    t->needs(Task::OldDW, d_mpm_labels->delTLabel);
    if (!d_mpm_flags->d_doGridReset) {
      t->needs(Task::OldDW, d_mpm_labels->gDisplacementLabel, Ghost::None);
      t->modifies(d_mpm_labels->gDisplacementLabel, mss);
    }

    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * setPrescribedMotion
 *-----------------------------------------------------------------------*/
void
SerialMPM::setPrescribedMotion(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing setPrescribedMotion");

    // Get the current time
    simTime_vartype simTimeVar;
    old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
    double time = simTimeVar;

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matlID = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gAcceleration;

      new_dw->getModifiable(
        gVelocity_star, d_mpm_labels->gVelocityStarLabel, matlID, patch);
      new_dw->getModifiable(
        gAcceleration, d_mpm_labels->gAccelerationLabel, matlID, patch);

      gAcceleration.initialize(Vector(0.0));

      // std::cout << " tval = ";
      // for (auto tval : d_prescribedTimes) {
      //   std::cout << tval << ", ";
      // }
      // std::cout << "\n";

      // Get F and Q from file by interpolating between available times
      auto t_upper_iter = std::upper_bound(
        d_prescribedTimes.begin(), d_prescribedTimes.end(), time);

      double ss{ 0.0 };
      double t_lower{ 0.0 };
      double t_upper{ 0.0 };
      int t_upper_index{ 0 };
      if (t_upper_iter == d_prescribedTimes.begin()) {
        t_lower       = *t_upper_iter;
        t_upper       = *(t_upper_iter + 1);
        ss            = (time - t_lower) / (t_upper - t_lower);
        t_upper_index = 1;
      } else if (t_upper_iter == d_prescribedTimes.end()) {
        t_lower       = *(t_upper_iter - 2);
        t_upper       = *(t_upper_iter - 1);
        ss            = (time - t_lower) / (t_upper - t_lower);
        t_upper_index = t_upper_iter - 1 - d_prescribedTimes.begin();
      } else {
        t_lower       = *(t_upper_iter - 1);
        t_upper       = *t_upper_iter;
        ss            = (time - t_lower) / (t_upper - t_lower);
        t_upper_index = t_upper_iter - d_prescribedTimes.begin();
      }

      // Interpolate to get the deformation gradient at the current time:
      auto F_lower = d_prescribedF[t_upper_index - 1];
      auto F_upper = d_prescribedF[t_upper_index];
      auto Ft      = (1 - ss) * F_lower + ss * F_upper;

      // Calculate the rate of the deformation gradient without the rotation:
      auto Fdot   = (F_upper - F_lower) / (t_upper - t_lower);
      auto Ft_inv = Ft.Inverse();
      auto L      = Fdot * Ft_inv;

      // Now we need to construct the rotation matrix and its time rate:
      // We are only interested in the rotation information at the next
      // specified time since the rotations specified should be relative to the
      // previously specified time. For example if I specify Theta=90 at
      // time=1.0, and Theta = 91 and time=2.0 the total rotation at time=2.0
      // will be 181 degrees.
      const double pi = M_PI; // 3.1415926535897932384626433832795028841972;
      const double degtorad = pi / 180.0;

      auto theta_upper = d_prescribedAngle[t_upper_index];
      auto thetat      = ss * theta_upper * degtorad;

      auto rot_axis_upper = d_prescribedRotationAxis[t_upper_index];
      Matrix3 QQ(thetat, rot_axis_upper);
      auto Qt = QQ.Transpose();

      auto thetadot = theta_upper * degtorad / (t_upper - t_lower);

      // Exact Deformation Update
      /*
       ** TODO ** Add computes for delT for this code to run.
       **
       ** Warning: Tries to access data that is outside bounds
       **/
      if (d_mpm_flags->d_exactDeformation) {
        // Check to see we do not exceed bounds
        int count              = 0;
        auto t_upper_iter_copy = t_upper_iter;
        // auto t_upper_index_copy = t_upper_index;
        while (++t_upper_iter_copy != d_prescribedTimes.end()) {
          // t_upper_index_copy = t_upper_iter_copy - d_prescribedTimes.begin();
          ++count;
          // std::cout << "t_upper_index_copy = " << t_upper_index_copy << "
          // count = " << count << "\n";
          if (count > 1) {
            break;
          }
        }

        // If there are at least two extra data points
        if (count > 1) {

          double t3 = d_prescribedTimes[t_upper_index + 1];
          double t4 = d_prescribedTimes[t_upper_index + 2];
          if (time == 0 && t4 != 0) {

            new_dw->put(delt_vartype(t3 - t_upper),
                        d_mpm_labels->delTLabel,
                        getLevel(patches));

          } else {

            F_lower = d_prescribedF[t_upper_index]; // last prescribed
                                                    // deformation gradient
            F_upper  = d_prescribedF[t_upper_index +
                                    1]; // next prescribed deformation gradient
            Ft       = (1 - ss) * F_lower + ss * F_upper;
            Ft_inv   = Ft.Inverse();
            Fdot     = (F_upper - F_lower) / (t3 - t_upper);
            thetadot = theta_upper * degtorad / (t3 - t_upper);

            double tst = t4 - t3;
            new_dw->put(
              delt_vartype(tst), d_mpm_labels->delTLabel, getLevel(patches));
          }
        }
      }

      // Construct Qdot:
      const double costhetat = cos(thetat);
      const double sinthetat = sin(thetat);
      Matrix3 Ident;
      Ident.Identity();
      Matrix3 aa(rot_axis_upper, rot_axis_upper);
      Matrix3 AA(rot_axis_upper);
      auto Qdot =
        (Ident - aa) * (-sinthetat * thetadot) + AA * costhetat * thetadot;

      // Now we need to compute the total previous rotation:
      Matrix3 R_previous;
      R_previous.Identity();
      for (auto ii = 0; ii < t_upper_index; ii++) {
        auto thetai = d_prescribedAngle[ii] * degtorad;
        auto ai     = d_prescribedRotationAxis[ii];
        Matrix3 Qi(thetai, ai);
        R_previous = Qi * R_previous;
      }

      // Fstar is the deformation gradient with the superimposed rotations
      // included Fdotstar is the rate of the deformation gradient with
      // superimposed rotations included
      auto Fstar          = Qt * R_previous * Ft;
      auto Fdotstar       = Qdot * R_previous * Ft + Qt * R_previous * Fdot;
      auto R_previous_inv = R_previous.Inverse();

      // Update grid velocities
      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector n     = *iter;
        Vector position = patch->getNodePosition(n).asVector();

        // Exact Deformation Update
        if (d_mpm_flags->d_exactDeformation) {
          gVelocity_star[n] = (F_upper * F_lower.Inverse() - Ident) *
                              R_previous_inv * QQ * position / delT;
        } else {
          gVelocity_star[n] =
            Fdotstar * Ft_inv * R_previous_inv * QQ * position;
        }

      } // Node Iterator

      if (!d_mpm_flags->d_doGridReset) {
        NCVariable<Vector> displacement;
        constNCVariable<Vector> displacementOld;
        new_dw->allocateAndPut(
          displacement, d_mpm_labels->gDisplacementLabel, matlID, patch);
        old_dw->get(displacementOld,
                    d_mpm_labels->gDisplacementLabel,
                    matlID,
                    patch,
                    Ghost::None,
                    0);
        for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
          IntVector c     = *iter;
          displacement[c] = displacementOld[c] + gVelocity_star[c] * delT;
        }
      } // d_doGridReset

    } // matl loop
  }   // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleComputeXPICVelocities
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeXPICVelocities(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeXPICVelocities");

  // Particle velocities
  Task* t_part = scinew Task("MPM::computeParticleVelocityXPIC",
                             this,
                             &SerialMPM::computeParticleVelocityXPIC);

  t_part->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
  t_part->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t_part->needs(Task::OldDW, d_mpm_labels->pDefGradLabel, Ghost::None);

  t_part->needs(Task::NewDW,
                   d_mpm_labels->gVelocityLabel,
                   Ghost::AroundCells,
                   d_numGhostNodes);

  t_part->computes(d_mpm_labels->pVelocityXPICLabel);

  sched->addTask(t_part, patches, matls);

  // Grid velocities
  Task* t_grid = scinew Task(
    "MPM::computeGridVelocityXPIC", this, &SerialMPM::computeGridVelocityXPIC);

  t_grid->needs(Task::OldDW,
                   d_mpm_labels->pXLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->needs(Task::OldDW,
                   d_mpm_labels->pMassLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->needs(Task::OldDW,
                   d_mpm_labels->pSizeLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->needs(Task::OldDW,
                   d_mpm_labels->pDefGradLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->needs(Task::NewDW,
                   d_mpm_labels->pVelocityXPICLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);

  t_grid->needs(
    Task::NewDW, d_mpm_labels->gMassLabel, Ghost::AroundCells, d_numGhostNodes);

  t_grid->computes(d_mpm_labels->gVelocityXPICLabel);

  sched->addTask(t_grid, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeParticleVelocityXPIC
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeParticleVelocityXPIC(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset*,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  for (auto patch : *patches) {
    printTask(patches, patch, cout_doing, "Doing computeParticleVelocityXPIC");

    auto interpolator      = d_mpm_flags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    size_t numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPMMatls; m++) {

      MPMMaterial* material =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

      auto mat_id    = material->getDWIndex();
      auto particles = old_dw->getParticleSubset(mat_id, patch);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpm_labels->pXLabel, particles);

      constParticleVariable<Matrix3> pSize, pDefGrad;
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, particles);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, particles);

      constNCVariable<Vector> gVelocity;
      new_dw->get(gVelocity,
                  d_mpm_labels->gVelocityLabel,
                  mat_id,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      ParticleVariable<Vector> pVelocityXPIC;
      new_dw->allocateAndPut(
        pVelocityXPIC, d_mpm_labels->pVelocityXPICLabel, particles);

      for (auto particle : *particles) {

        interpolator->findCellAndWeights(
          pX[particle], ni, S, pSize[particle], pDefGrad[particle]);
        Vector pVelocity(0.0, 0.0, 0.0);
        for (int k = 0; k < numInfluenceNodes; k++) {
          pVelocity += gVelocity[ni[k]] * S[k];
        }
        pVelocityXPIC[particle] = pVelocity;
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * computeGridVelocityXPIC
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeGridVelocityXPIC(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  for (auto patch : *patches) {
    printTask(patches, patch, cout_doing, "Doing computeGridVelocityXPIC");

    auto interpolator      = d_mpm_flags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    size_t numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPMMatls; m++) {

      MPMMaterial* material =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

      auto mat_id    = material->getDWIndex();
      auto particles = old_dw->getParticleSubset(mat_id, patch);

      constParticleVariable<double> pMass;
      old_dw->get(pMass, d_mpm_labels->pMassLabel, particles);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpm_labels->pXLabel, particles);

      constParticleVariable<Matrix3> pSize, pDefGrad;
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, particles);
      old_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel, particles);

      constParticleVariable<Vector> pVelocityXPIC;
      new_dw->get(pVelocityXPIC, d_mpm_labels->pVelocityXPICLabel, particles);

      constNCVariable<double> gMass;
      new_dw->get(gMass,
                  d_mpm_labels->gMassLabel,
                  mat_id,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      NCVariable<Vector> gVelocityXPIC;
      new_dw->allocateAndPut(
        gVelocityXPIC, d_mpm_labels->gVelocityXPICLabel, mat_id, patch);
      gVelocityXPIC.initialize(Vector(0.0, 0.0, 0.0));

      IntVector node;
      for (auto particle : *particles) {

        interpolator->findCellAndWeights(
          pX[particle], ni, S, pSize[particle], pDefGrad[particle]);
        Vector pMomentum = pVelocityXPIC[particle] * pMass[particle];
        for (int k = 0; k < numInfluenceNodes; k++) {
          node = ni[k];
          if (patch->containsNode(node)) {
            gVelocityXPIC[node] += pMomentum * S[k];
          }
        }
      }

      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); ++iter) {
        node = *iter;
        gVelocityXPIC[node] /= gMass[node];
      }
    }
  }
}

//  Schedule the reduction of variables that are computed multiple times in a
//  timestep Use Label->schedReductionTask( true ); to the tell scheduler to
//  perform the reduction. The actual task is inside MPIScheduler.
void
SerialMPM::scheduleReduceVars(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* matls)
{
  if (!d_mpm_flags->d_reductionVars->sumTransmittedForce) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleReduceVars");

  Task* t = scinew Task("MPM::reductionTask", Task::Reduction);

  // Create reductionMatlSubSet that includes all mpm matls
  // and the global matl.
  const MaterialSubset* global_mss = t->getGlobalMatlSubset();
  const MaterialSubset* mpm_mss    = (matls ? matls->getUnion() : nullptr);

  auto* reduction_mss = scinew MaterialSubset();
  reduction_mss->add(global_mss->get(0));

  size_t nMatls = d_materialManager->getNumMaterials("MPM");

  if (nMatls > 1) { // ignore for single matl problems
    for (size_t m = 0; m < nMatls; m++) {
      reduction_mss->add(mpm_mss->get(m));
    }
  }

  reduction_mss->addReference();

  // Tell the scheduler to reduce this variable
  d_mpm_labels->SumTransmittedForceLabel->schedReductionTask(true);
  t->computes(
    d_mpm_labels->SumTransmittedForceLabel, reduction_mss, Task::OutOfDomain);

  sched->addTask(t, patches, matls);

  if (reduction_mss && reduction_mss->removeReference()) {
    delete reduction_mss;
  }
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateToParticlesAndUpdate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(
    patches, cout_doing, "MPM::scheduleInterpolateToParticlesAndUpdate");

  Task* t = scinew Task("MPM::interpolateToParticlesAndUpdate",
                        this,
                        &SerialMPM::interpolateToParticlesAndUpdate);

  t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->needs(
    Task::NewDW, d_mpm_labels->gAccelerationLabel, gac, d_numGhostNodes);
  t->needs(
    Task::NewDW, d_mpm_labels->gVelocityStarLabel, gac, d_numGhostNodes);
  if (d_mpm_flags->d_useXPIC) {
    t->needs(
      Task::NewDW, d_mpm_labels->gVelocityXPICLabel, gac, d_numGhostNodes);
  }
  t->needs(
    Task::NewDW, d_mpm_labels->gTemperatureRateLabel, gac, d_numGhostNodes);
  t->needs(
    Task::NewDW, d_mpm_labels->frictionalWorkLabel, gac, d_numGhostNodes);
  t->needs(Task::OldDW, d_mpm_labels->pXLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pMassLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pTemperatureLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pVelocityLabel, gnone);
  if (d_mpm_flags->d_useXPIC) {
    t->needs(Task::OldDW, d_mpm_labels->pVelocityXPICLabel, gnone);
  }
  t->needs(Task::OldDW, d_mpm_labels->pDispLabel, gnone);
  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, gnone);
  if (d_mpm_flags->d_enableDEM) {
    t->needs(Task::OldDW, d_mpm_labels->pRadiusLabel, gnone);
    t->needs(Task::OldDW, d_mpm_labels->pX0Label, gnone);
    t->needs(Task::OldDW, d_mpm_labels->pRigidBodyIDLabel, gnone);
    t->needs(Task::OldDW, d_mpm_labels->pAngularVelocityLabel, gnone);
  }
  t->needs(Task::NewDW, d_mpm_labels->pdTdtLabel_preReloc, gnone);
  t->needs(Task::NewDW, d_mpm_labels->pLocalizedMPMLabel, gnone);
  t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel_preReloc, gnone);
  t->modifies(d_mpm_labels->pVolumeLabel_preReloc);

  if (d_mpm_flags->d_useLoadCurves) {
    t->needs(Task::OldDW, d_mpm_labels->pLoadCurveIDLabel, Ghost::None);
  }

  if (d_mpm_flags->d_withICE) {
    t->needs(Task::NewDW, d_mpm_labels->dTdt_NCLabel, gac, d_numGhostNodes);
    t->needs(
      Task::NewDW, d_mpm_labels->massBurnFractionLabel, gac, d_numGhostNodes);
  }

  t->computes(d_mpm_labels->pDispLabel_preReloc);
  t->computes(d_mpm_labels->pVelocityLabel_preReloc);
  t->computes(d_mpm_labels->pAccelerationLabel_preReloc);
  t->computes(d_mpm_labels->pXLabel_preReloc);
  if (d_mpm_flags->d_enableDEM) {
    t->computes(d_mpm_labels->pX0Label_preReloc);
  }
  t->computes(d_mpm_labels->pParticleIDLabel_preReloc);
  t->computes(d_mpm_labels->pTemperatureLabel_preReloc);
  t->computes(d_mpm_labels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpm_labels->pMassLabel_preReloc);
  t->computes(d_mpm_labels->pSizeLabel_preReloc);
  t->computes(d_mpm_labels->pXXLabel);

  //__________________________________
  //  reduction variables
  if (d_mpm_flags->d_reductionVars->momentum) {
    t->computes(d_mpm_labels->TotalMomentumLabel);
  }
  if (d_mpm_flags->d_reductionVars->KE) {
    t->computes(d_mpm_labels->KineticEnergyLabel);
  }
  if (d_mpm_flags->d_reductionVars->thermalEnergy) {
    t->computes(d_mpm_labels->ThermalEnergyLabel);
  }
  if (d_mpm_flags->d_reductionVars->centerOfMass) {
    t->computes(d_mpm_labels->CenterOfMassPositionLabel);
  }
  if (d_mpm_flags->d_reductionVars->mass) {
    t->computes(d_mpm_labels->TotalMassLabel);
  }
  if (d_mpm_flags->d_reductionVars->volDeformed) {
    t->computes(d_mpm_labels->TotalVolumeDeformedLabel);
  }

  // debugging scalar
  if (d_mpm_flags->d_withColor) {
    t->needs(Task::OldDW, d_mpm_labels->pColorLabel, Ghost::None);
    t->computes(d_mpm_labels->pColorLabel_preReloc);
  }

  // Carry Forward particle refinement flag
  if (d_mpm_flags->d_refineParticles) {
    t->needs(Task::OldDW, d_mpm_labels->pRefinedLabel, Ghost::None);
    t->computes(d_mpm_labels->pRefinedLabel_preReloc);
  }

  // Carry forward external heat flux for switch from explicit to implicit
  t->needs(Task::OldDW, d_mpm_labels->pExternalHeatFluxLabel, Ghost::None);
  t->computes(d_mpm_labels->pExternalHeatFluxLabel_preReloc);

  auto* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->needs(Task::OldDW, d_mpm_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(d_mpm_labels->NC_CCweightLabel, z_matl);

  sched->addTask(t, patches, matls);

  // The task will have a reference to z_matl
  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}

/*!----------------------------------------------------------------------
 * interpolateToParticlesAndUpdate
 *-----------------------------------------------------------------------*/
void
SerialMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset*,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;
  Ghost::GhostType gac   = Ghost::AroundCells;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(
      patches, patch, cout_doing, "Doing interpolateToParticlesAndUpdate");

    auto interpolator      = d_mpm_flags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    double totalmass      = 0;
    double partvoldef     = 0.;
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke       = 0;
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    const Level* level = getLevel(patches);
    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, level);

    BBox domain;
    level->getSpatialRange(domain);
    Point domain_low  = domain.min();
    Point domain_high = domain.max();

    // bool combustion_problem=false;

    // Material* reactant;
    // int RMI = -99;
    // reactant = d_materialManager->getMaterialByName("reactant");
    // if(reactant != 0){
    //   RMI = reactant->getDWIndex();
    //   combustion_problem=true;
    // }
    double move_particles = 1.;
    if (!d_mpm_flags->d_doGridReset) {
      move_particles = 0.;
    }

    // Copy NC_CCweight (only material 0)
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_new;
    old_dw->get(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(
      NC_CCweight_new, d_mpm_labels->NC_CCweightLabel, 0, patch);
    NC_CCweight_new.copyData(NC_CCweight);

    std::vector<Patch::FaceType> bf;
    patch->getBoundaryFaces(bf);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Copy particle IDs and particle size
      constParticleVariable<long64> pParticleID;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<long64> pParticleID_new;
      ParticleVariable<Matrix3> pSize_new;
      old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(
        pParticleID_new, d_mpm_labels->pParticleIDLabel_preReloc, pset);
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(
        pSize_new, d_mpm_labels->pSizeLabel_preReloc, pset);
      pParticleID_new.copyData(pParticleID);
      pSize_new.copyData(pSize);

      // Copy needed for switch from explicit to implicit MPM
      constParticleVariable<double> pExtHeatFlux;
      ParticleVariable<double> pExtHeatFlux_new;
      old_dw->get(pExtHeatFlux, d_mpm_labels->pExternalHeatFluxLabel, pset);
      new_dw->allocateAndPut(
        pExtHeatFlux_new, d_mpm_labels->pExternalHeatFluxLabel_preReloc, pset);
      pExtHeatFlux_new.copyData(pExtHeatFlux);

      // Get particle variables
      constParticleVariable<int> pLocalized_new;
      new_dw->get(pLocalized_new, d_mpm_labels->pLocalizedMPMLabel, pset);

      constParticleVariable<double> pMass, pTemperature, pdTdt_new;
      old_dw->get(pMass, d_mpm_labels->pMassLabel, pset);
      old_dw->get(pTemperature, d_mpm_labels->pTemperatureLabel, pset);
      new_dw->get(pdTdt_new, d_mpm_labels->pdTdtLabel_preReloc, pset);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpm_labels->pXLabel, pset);

      constParticleVariable<Vector> pVelocity, pDisp;
      old_dw->get(pVelocity, d_mpm_labels->pVelocityLabel, pset);
      old_dw->get(pDisp, d_mpm_labels->pDispLabel, pset);

      constParticleVariable<double> pRadius;
      constParticleVariable<long64> pRigidBodyID;
      constParticleVariable<Vector> pAngularVelocity;
      if (d_mpm_flags->d_enableDEM) {
        old_dw->get(pRadius, d_mpm_labels->pRadiusLabel, pset);
        old_dw->get(pRigidBodyID, d_mpm_labels->pRigidBodyIDLabel, pset);
        old_dw->get(pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, pset);
      }

      constParticleVariable<Vector> pVelocityXPIC;
      if (d_mpm_flags->d_useXPIC) {
        new_dw->get(pVelocityXPIC, d_mpm_labels->pVelocityXPICLabel, pset);
      }

      constParticleVariable<Matrix3> pDefGrad_new;
      new_dw->get(pDefGrad_new, d_mpm_labels->pDefGradLabel_preReloc, pset);

      // Allocate updated particle variables
      ParticleVariable<double> pMass_new, pVolume_new, pTemp_new, pTempPrev_new;
      new_dw->allocateAndPut(
        pMass_new, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->getModifiable(
        pVolume_new, d_mpm_labels->pVolumeLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pTemp_new, d_mpm_labels->pTemperatureLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pTempPrev_new, d_mpm_labels->pTempPreviousLabel_preReloc, pset);

      ParticleVariable<Point> pX_new, pX0_new, pXx;
      new_dw->allocateAndPut(pX_new, d_mpm_labels->pXLabel_preReloc, pset);
      if (d_mpm_flags->d_enableDEM) {
        constParticleVariable<Point> pX0;
        old_dw->get(pX0, d_mpm_labels->pX0Label, pset);
        new_dw->allocateAndPut(pX0_new, d_mpm_labels->pX0Label_preReloc, pset);
        pX0_new.copyData(pX0);
      }
      new_dw->allocateAndPut(pXx, d_mpm_labels->pXXLabel, pset);

      ParticleVariable<Vector> pVelocity_new, pDisp_new, pAcc_new;
      new_dw->allocateAndPut(
        pVelocity_new, d_mpm_labels->pVelocityLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pDisp_new, d_mpm_labels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(
        pAcc_new, d_mpm_labels->pAccelerationLabel_preReloc, pset);

      // Get grid variables
      constNCVariable<double> gTemperatureRate, frictionTempRate;
      new_dw->get(gTemperatureRate,
                  d_mpm_labels->gTemperatureRateLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpm_labels->frictionalWorkLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);

      constNCVariable<double> dTdt, massBurnFrac;
      if (d_mpm_flags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpm_labels->dTdt_NCLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpm_labels->massBurnFractionLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create, patch, gac, d_numGhostParticles);
        new_dw->allocateTemporary(
          massBurnFrac_create, patch, gac, d_numGhostParticles);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);

        dTdt         = dTdt_create;         // reference created data
        massBurnFrac = massBurnFrac_create; // reference created data
      }

      constNCVariable<Vector> gVelocityStar, gAcceleration;
      new_dw->get(gVelocityStar,
                  d_mpm_labels->gVelocityStarLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gAcceleration,
                  d_mpm_labels->gAccelerationLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);

      constNCVariable<Vector> gVelocityXPIC;
      if (d_mpm_flags->d_useXPIC) {
        new_dw->get(gVelocityXPIC,
                    d_mpm_labels->gVelocityXPICLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
      }

      auto* delset = scinew ParticleSubset(0, matID, patch);

      // Map to store master particle kinematics for synchronization
      struct MasterState {
        Point pos;
        Point pos_old;
        Vector vel;
        Vector omega;
      };
      std::map<long64, MasterState> master_states;

      double Cp       = mpm_matl->getSpecificHeat();
      double rho_init = mpm_matl->getInitialDensity();

      double rho_frac_min = 0.;
      /*
      if(m == RMI){
        rho_frac_min = .1;
      }
      */

      // Loop over particles
      for (auto idx : *pset) {

        interpolator->findCellAndWeights(
          pX[idx], ni, S, pSize[idx], pDefGrad_new[idx]);

        Vector velocity(0.0, 0.0, 0.0);
        Vector acceleration(0.0, 0.0, 0.0);
        Vector velocityXPIC(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0.0;
        double burnFraction = 0.0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          velocity += gVelocityStar[node] * S[k];
          acceleration += gAcceleration[node] * S[k];

          if (d_mpm_flags->d_useXPIC) {
            velocityXPIC += gVelocityXPIC[node] * S[k];
          }

#ifdef CHECK_ISFINITE
          if (!std::isfinite(velocity.x()) || !std::isfinite(velocity.y()) ||
              !std::isfinite(velocity.z()) ||
              !std::isfinite(acceleration.x()) ||
              !std::isfinite(acceleration.y()) ||
              !std::isfinite(acceleration.z())) {
            std::cout << "particle ID = " << pParticleID[idx]
                      << " node = " << node << " k = " << k
                      << " S[k] = " << S[k] << " v_g* = "
                      << gVelocityStar[node]
                      //<< " v_g(2) = " << gVelocityXPIC[node]
                      << " a_g* = " << gAcceleration[node] << "\n";
          }
#endif

          fricTempRate =
            frictionTempRate[node] * d_mpm_flags->d_addFrictionWork;
          tempRate +=
            (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
          burnFraction += massBurnFrac[node] * S[k];
        }

        // Update the particle's position and velocity
        if (d_mpm_flags->d_useXPIC) {
          pX_new[idx] = pX[idx] + velocity * delT -
                        0.5 *
                          (acceleration * delT + pVelocity[idx] -
                           2.0 * pVelocityXPIC[idx] + velocityXPIC) *
                          delT;
          pVelocity_new[idx] =
            2.0 * pVelocityXPIC[idx] - velocityXPIC + acceleration * delT;
          pDisp_new[idx] = pDisp[idx] + (pX_new[idx] - pX[idx]);
        } else {
          pX_new[idx]        = pX[idx] + velocity * delT * move_particles;
          pDisp_new[idx]     = pDisp[idx] + velocity * delT;
          pVelocity_new[idx] = pVelocity[idx] + acceleration * delT;
        }

        // Collect master particle states for discrete or rigid materials
        if (d_mpm_flags->d_enableDEM && (mpm_matl->isDiscrete() || mpm_matl->getIsRigid()) && pRadius[idx] > 0) {
          master_states[pRigidBodyID[idx]] = {pX_new[idx], pX[idx], pVelocity_new[idx], pAngularVelocity[idx]};
        }

        // Explicit boundary enforcement for discrete particles
        // In quasi-2D simulations, the domain thickness in one dimension (e.g., Y)
        // might be much smaller than the physical radius of a discrete object.
        // We only enforce symmetry boundary corrections if the domain is thick
        // enough to contain at least one full particle diameter. Otherwise,
        // the correction would push the particle outside the opposite boundary.
        if (d_mpm_flags->d_enableDEM && mpm_matl->isDiscrete() && pMass[idx] > 0) {
          double radius = pRadius[idx];
          for (auto face : bf) {
            if (patch->haveBC(face, matID, "symmetry", "Symmetric")) {
              if (face == Patch::xminus &&
                  (domain_high.x() - domain_low.x() > 2.0 * radius) &&
                  pX_new[idx].x() - radius < domain_low.x()) {
                pVelocity_new[idx] =
                  Vector(0, pVelocity_new[idx].y(), pVelocity_new[idx].z());
                pX_new[idx] = Point(domain_low.x() + radius,
                                    pX_new[idx].y(),
                                    pX_new[idx].z());
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
              if (face == Patch::xplus &&
                  (domain_high.x() - domain_low.x() > 2.0 * radius) &&
                  pX_new[idx].x() + radius > domain_high.x()) {
                pVelocity_new[idx] =
                  Vector(0, pVelocity_new[idx].y(), pVelocity_new[idx].z());
                pX_new[idx] = Point(domain_high.x() - radius,
                                    pX_new[idx].y(),
                                    pX_new[idx].z());
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
              if (face == Patch::yminus &&
                  (domain_high.y() - domain_low.y() > 2.0 * radius) &&
                  pX_new[idx].y() - radius < domain_low.y()) {
                pVelocity_new[idx] =
                  Vector(pVelocity_new[idx].x(), 0, pVelocity_new[idx].z());
                pX_new[idx] = Point(pX_new[idx].x(),
                                    domain_low.y() + radius,
                                    pX_new[idx].z());
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
              if (face == Patch::yplus &&
                  (domain_high.y() - domain_low.y() > 2.0 * radius) &&
                  pX_new[idx].y() + radius > domain_high.y()) {
                pVelocity_new[idx] =
                  Vector(pVelocity_new[idx].x(), 0, pVelocity_new[idx].z());
                pX_new[idx] = Point(pX_new[idx].x(),
                                    domain_high.y() - radius,
                                    pX_new[idx].z());
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
              if (face == Patch::zminus &&
                  (domain_high.z() - domain_low.z() > 2.0 * radius) &&
                  pX_new[idx].z() - radius < domain_low.z()) {
                pVelocity_new[idx] =
                  Vector(pVelocity_new[idx].x(), pVelocity_new[idx].y(), 0);
                pX_new[idx] = Point(pX_new[idx].x(),
                                    pX_new[idx].y(),
                                    domain_low.z() + radius);
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
              if (face == Patch::zplus &&
                  (domain_high.z() - domain_low.z() > 2.0 * radius) &&
                  pX_new[idx].z() + radius > domain_high.z()) {
                pVelocity_new[idx] =
                  Vector(pVelocity_new[idx].x(), pVelocity_new[idx].y(), 0);
                pX_new[idx] = Point(pX_new[idx].x(),
                                    pX_new[idx].y(),
                                    domain_high.z() - radius);
                pDisp_new[idx] = pX_new[idx] - pX[idx];
              }
            }
          }
        }

        pAcc_new[idx] = acceleration;

#ifdef CHECK_ISFINITE
        if (!std::isfinite(pVelocity_new[idx].x()) ||
            !std::isfinite(pVelocity_new[idx].y()) ||
            !std::isfinite(pVelocity_new[idx].z())) {
          std::cout << "particle ID = " << pParticleID[idx]
                    << " v_p = " << pVelocity[idx] << " a_p = " << acceleration
                    << "\n";
        }
#endif

        // pXx is only useful if we're not in normal grid resetting mode.
        pXx[idx]       = pX[idx] + pDisp_new[idx];
        pTemp_new[idx] = pTemperature[idx] + (tempRate + pdTdt_new[idx]) * delT;
        pTempPrev_new[idx] = pTemperature[idx]; // for thermal stress

        // Clamp negative temperatures
        if (pTemp_new[idx] < 0) {
          pTemp_new[idx] = pTemperature[idx];
        }

        if (cout_heat.active()) {
          cout_heat << "MPM::Particle = " << pParticleID[idx]
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate << " dT = " << (tempRate * delT)
                    << " T_new = " << pTemp_new[idx] << "\n";
        }

        double rho;
        if (pVolume_new[idx] > 0.) {
          rho =
            std::max(pMass[idx] / pVolume_new[idx], rho_frac_min * rho_init);
        } else {
          rho = rho_init;
        }

        pMass_new[idx] = Max(pMass[idx] * (1.0 - burnFraction), 0.);
        // std::cout << "m = " << pMass[idx] << " burnFraction = " <<
        // burnFraction
        //           << "m_new = " << pMass_new[idx] << "\n";
        pVolume_new[idx] = pMass_new[idx] / rho;

        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
        ke += .5 * pMass[idx] * pVelocity_new[idx].length2();
        CMX = CMX + (pX_new[idx] * pMass[idx]).asVector();
        totalMom += pVelocity_new[idx] * pMass[idx];
        totalmass += pMass_new[idx];
        partvoldef += pVolume_new[idx];

      } // End loop over particles

      // Synchronize dummy particles with master particle for discrete or rigid materials
      if (d_mpm_flags->d_enableDEM && (mpm_matl->isDiscrete() || mpm_matl->getIsRigid())) {
        for (auto idx : *pset) {
          // If it's a dummy particle (discrete only) OR any particle in a rigid material (except the master)
          if ((mpm_matl->isDiscrete() && pMass[idx] <= 0) || (mpm_matl->getIsRigid() && pRadius[idx] <= 0)) {
            long64 rbID = pRigidBodyID[idx];
            if (master_states.count(rbID)) {
              const auto& master = master_states[rbID];
              Vector arm = pX[idx] - master.pos_old;
              double arm_len = arm.length();
              
              if (master.omega.length2() > 1e-12) {
                 // Predict new position with explicit rotation
                 Vector arm_new = arm + Cross(master.omega, arm) * delT;
                 
                 // Corrector: Enforce rigid body distance constraint
                 if (arm_new.length() > 1e-12) {
                    arm_new = arm_new * (arm_len / arm_new.length());
                 }
                 
                 pX_new[idx] = master.pos + arm_new;
                 pVelocity_new[idx] = (pX_new[idx] - pX[idx]) / delT;
              } else {
                 pX_new[idx] = master.pos + arm;
                 pVelocity_new[idx] = master.vel;
              }
              pDisp_new[idx] = pX_new[idx] - pX[idx];
            }
          }
        }
      }

      // If load curves are being used with VelocityBC then apply
      // these BCs to the boundary particles
      if (d_mpm_flags->d_useLoadCurves) {

        std::vector<VelocityBC*> vbcP;
        bool do_VelocityBCs = false;
        for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          std::string bcType = bc->getType();
          if (bcType == "Velocity") {
            do_VelocityBCs = true;
            auto* vbc      = dynamic_cast<VelocityBC*>(bc.get());
            vbcP.push_back(vbc);
          }
        }

        // std::cout << "do_VelocityBCs = " << do_VelocityBCs << "\n";
        if (do_VelocityBCs) {

          // Get the current time
          simTime_vartype simTimeVar;
          old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
          double time = simTimeVar;

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, d_mpm_labels->pLoadCurveIDLabel, pset);

          // Iterate over the particles
          for (int idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (!(loadCurveID < 0)) {
              VelocityBC* vbc = vbcP[loadCurveID];
              pVelocity_new[idx] =
                vbc->getVelocityVector(pX[idx], pDisp[idx], time);
              pDisp_new[idx] = pDisp[idx] + pVelocity_new[idx] * delT;
              pX_new[idx] =
                pX[idx] + pVelocity_new[idx] * delT * move_particles;
              // std::cout << " Load curve ID = " << loadCurveID
              //           << " V = " << pVelocity_new[idx]
              //           << " U = " << pDisp_new[idx]
              //           << " x = " << pX_new[idx]
              //           << " num = " << pset->numParticles() << "\n";
            }
          }
        }
      }

      // Delete particles that have left the domain
      // This is only needed if extra cells are being used.
      // Also delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (auto idx : *pset) {

        /*
        if (pParticleID[idx] == 562644844544) {
          std::cout << "pID=" << pParticleID[idx] << " pLocalized_new = " <<
        pLocalized_new[idx] << "\n";
        }
        */

        if ((pMass_new[idx] <= d_mpm_flags->d_minPartMass && !mpm_matl->isDiscrete()) ||
            (pTemp_new[idx] < 0.0) || (pLocalized_new[idx] == -999)) {
          if (d_mpm_flags->d_erosionAlgorithm != "none") {
            delset->addParticle(idx);
          }
          if (mpm_matl->isDiscrete()) {
             std::cout << "  DEBUG: DELETING discrete particle ID: " << pParticleID[idx]
                       << " in interpolateToParticlesAndUpdate"
                       << " mass: " << pMass_new[idx] << " temp: " << pTemp_new[idx]
                       << " localized: " << pLocalized_new[idx] << std::endl;
          }
          proc0cout
            << "\n Warning: particle " << pParticleID[idx]
            << " being deleted: low mass or low temperature or localized\n";
          proc0cout << "\t mass = " << pMass_new[idx]
                    << " temperature = " << pTemp_new[idx]
                    << " localized = " << pLocalized_new[idx] << "\n";
#ifdef CHECK_PARTICLE_DELETION
          proc0cout << "In " << __FILE__ << ":" << __LINE__ << "\n";
          proc0cout << "Material = " << m
                    << " Deleted Particle = " << pParticleID_new[idx]
                    << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
                    << " vold = " << pVelocity[idx]
                    << " vnew = " << pVelocity_new[idx]
                    << " massold = " << pMass[idx]
                    << " massnew = " << pMass_new[idx]
                    << " tempold = " << pTemperature[idx]
                    << " tempnew = " << pTemp_new[idx]
                    << " pLocalized = " << pLocalized_new[idx]
                    << " volnew = " << pVolume_new[idx] << "\n";
#endif
        }

        if (pVelocity_new[idx].length() > d_mpm_flags->d_maxVel) {
          if (d_mpm_flags->d_deleteRogueParticles) {
            if (d_mpm_flags->d_erosionAlgorithm != "none") {
              delset->addParticle(idx);
            }
            if (mpm_matl->isDiscrete()) {
               std::cout << "  DEBUG: DELETING discrete particle ID: " << pParticleID[idx]
                         << " hit speed ceiling. Vel: " << pVelocity_new[idx].length() << std::endl;
            }
            std::cout << "\n Warning: particle " << pParticleID[idx]
                      << " hit speed ceiling #1. Deleting particle."
                      << "\n";
          } else {
            if (pVelocity_new[idx].length() >= pVelocity[idx].length()) {
              pVelocity_new[idx] =
                (pVelocity_new[idx] / pVelocity_new[idx].length()) *
                (d_mpm_flags->d_maxVel * .9);
              std::cout << "\n Warning: particle " << pParticleID[idx]
                        << " hit speed ceiling #1. Modifying particle velocity "
                           "accordingly."
                        << "\n";
            }
          }
        }
      }

      /*
      std::cout << "Particles in domain = " << pset->numParticles() << "\n";
      std::cout << "Particles to be deleted = " << delset->numParticles() <<
      "\n";
      */

      new_dw->deleteParticles(delset);

      //__________________________________
      //  particle debugging label-- carry forward
      if (d_mpm_flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_mpm_labels->pColorLabel, pset);
        new_dw->allocateAndPut(
          pColor_new, d_mpm_labels->pColorLabel_preReloc, pset);
        pColor_new.copyData(pColor);
      }
    }

    // DON'T MOVE THESE!!!
    //__________________________________
    //  reduction variables
    if (d_mpm_flags->d_reductionVars->mass) {
      new_dw->put(sum_vartype(totalmass), d_mpm_labels->TotalMassLabel);
    }
    if (d_mpm_flags->d_reductionVars->volDeformed) {
      new_dw->put(sum_vartype(partvoldef),
                  d_mpm_labels->TotalVolumeDeformedLabel);
    }
    if (d_mpm_flags->d_reductionVars->momentum) {
      new_dw->put(sumvec_vartype(totalMom), d_mpm_labels->TotalMomentumLabel);
    }
    if (d_mpm_flags->d_reductionVars->KE) {
      new_dw->put(sum_vartype(ke), d_mpm_labels->KineticEnergyLabel);
    }
    if (d_mpm_flags->d_reductionVars->thermalEnergy) {
      new_dw->put(sum_vartype(thermal_energy),
                  d_mpm_labels->ThermalEnergyLabel);
    }
    if (d_mpm_flags->d_reductionVars->centerOfMass) {
      new_dw->put(sumvec_vartype(CMX), d_mpm_labels->CenterOfMassPositionLabel);
    }

    // std::cout << "Solid mass lost this timestep = " << massLost << "\n";
    // std::cout << "Solid momentum after advection = " << totalMom << "\n";

    // std::cout << "THERMAL ENERGY " << thermal_energy << "\n";

    // delete interpolator;
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeDeformationGradient
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeDeformationGradient(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* Create a task for computing the deformation gradient */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeDeformationGradient");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeDeformationGradient",
                        this,
                        &SerialMPM::computeDeformationGradient);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    // Add requires and computes for vel grad/def grad
    d_defGradComputer->addComputesAndRequires(t, mpm_matl, patches);
    if (d_mpm_flags->d_enableDEM) {
      t->needs(Task::OldDW, d_mpm_labels->pRadiusLabel, Ghost::None);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeDeformationGradient
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeDeformationGradient(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  printTask(
    patches, patches->get(0), cout_doing, "Doing computeDeformationGradient");

  if (cout_doing.active()) {
    cout_doing << "Before compute def grad: old_dw\n";
    old_dw->print();
  }

  // Compute deformation gradient
  d_defGradComputer->computeDeformationGradient(patches, old_dw, new_dw);

  if (cout_doing.active()) {
    cout_doing << "After compute def grad: new_dw\n";
    new_dw->print();
  }

  /*
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
    static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    int matID = mpm_matl->getDWIndex();
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Matrix3> pVelGrad_mid, pDefGrad_mid;
      new_dw->get(pVelGrad_mid,  d_mpm_labels->pVelGradLabel_preReloc, pset);
      new_dw->get(pDefGrad_mid,  d_mpm_labels->pDefGradLabel_preReloc, pset);
      for (auto particle : *pset) {
      std::cout << "After compute vel/def gradients: material = " << m
                << " particle = " << particle << "\n"
                << " L_mid = " << pVelGrad_mid[particle] << "\n"
                << " F_mid = " << pDefGrad_mid[particle] << "\n";
      }
    }
  }
  */
}

/////////////////////////////////////////////////////////////////////////
/*!  **WARNING** In addition to the stresses and deformations, the internal
 *               heat rate in the particles (pdTdtLabel)
 *               is computed here */
/////////////////////////////////////////////////////////////////////////
void
SerialMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* Create a task for computing the stress tensor */
  scheduleUnrotateStressAndDeformationRate(sched, patches, matls);

  /* Create a task for computing the stress tensor */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeStressTensor");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task(
    "MPM::computeStressTensor", this, &SerialMPM::computeStressTensor);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    // Add requires and computes for constitutive model
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    t->computes(d_mpm_labels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_mpm_labels->delTLabel, getLevel(patches));

  if (d_mpm_flags->d_reductionVars->accStrainEnergy ||
      d_mpm_flags->d_reductionVars->strainEnergy) {
    t->computes(d_mpm_labels->StrainEnergyLabel);
  }

  sched->addTask(t, patches, matls);

  scheduleRotateStress(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleUnrotateStressAndDeformationRate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleUnrotateStressAndDeformationRate(SchedulerP& sched,
                                                    const PatchSet* patches,
                                                    const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(
    patches, cout_doing, "MPM::scheduleUnrotateStressAndDeformationRate");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeUnrotatedStressAndDeformationRate",
                        this,
                        &SerialMPM::computeUnrotatedStressAndDeformationRate);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->needs(
        Task::OldDW, d_mpm_labels->pParticleIDLabel, matlset, Ghost::None);
      t->needs(
        Task::OldDW, d_mpm_labels->pPolarDecompRLabel, matlset, Ghost::None);
      t->needs(
        Task::NewDW, d_mpm_labels->pPolarDecompRMidLabel, matlset, Ghost::None);
      t->needs(
        Task::OldDW, d_mpm_labels->pStressLabel, matlset, Ghost::None);
      t->needs(Task::NewDW,
                  d_mpm_labels->pVelGradLabel_preReloc,
                  matlset,
                  Ghost::None);

      t->computes(d_mpm_labels->pDeformRateMidLabel, matlset);
      t->computes(d_mpm_labels->pStressUnrotatedLabel, matlset);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeUnrotatedStressAndDeformationRate
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeUnrotatedStressAndDeformationRate(const ProcessorGroup*,
                                                    const PatchSubset* patches,
                                                    const MaterialSubset*,
                                                    DataWarehouse* old_dw,
                                                    DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing, "Doing computeUnrotate");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch   = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pStress_old, pR_old, pR_mid,
          pVelGrad_mid;
        old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
        old_dw->get(pR_old, d_mpm_labels->pPolarDecompRLabel, pset);
        new_dw->get(pR_mid, d_mpm_labels->pPolarDecompRMidLabel, pset);
        old_dw->get(pStress_old, d_mpm_labels->pStressLabel, pset);
        new_dw->get(pVelGrad_mid, d_mpm_labels->pVelGradLabel_preReloc, pset);

        ParticleVariable<Matrix3> pDeformRate_mid, pStress_old_unrotated;
        new_dw->allocateAndPut(
          pDeformRate_mid, d_mpm_labels->pDeformRateMidLabel, pset);
        new_dw->allocateAndPut(
          pStress_old_unrotated, d_mpm_labels->pStressUnrotatedLabel, pset);

        for (auto particle : *pset) {
          pStress_old_unrotated[particle] =
            (pR_old[particle].Transpose()) *
            (pStress_old[particle] * pR_old[particle]);
          Matrix3 DD =
            (pVelGrad_mid[particle] + pVelGrad_mid[particle].Transpose()) * 0.5;
          pDeformRate_mid[particle] =
            (pR_mid[particle].Transpose()) * (DD * pR_mid[particle]);
        }
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * computeStressTensor
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeStressTensor(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{

  printTask(patches, patches->get(0), cout_doing, "Doing computeStressTensor");

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

    // Compute stress
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cout_dbg.active()) {
      cout_dbg << " CM = " << cm;
    }

    cm->setWorld(UintahParallelComponent::d_myworld);
#ifdef TIME_COMPUTE_STRESS
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    cm->computeStressTensor(patches, mpm_matl, old_dw, new_dw);
#ifdef TIME_COMPUTE_STRESS
    end = std::chrono::system_clock::now();
    std::cout << "Compute stress : Time taken = "
              << std::chrono::duration<double>(end - start).count() << "\n";
#endif

    if (cout_dbg.active()) {
      cout_dbg << " Exit\n";
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleRotateStress
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleRotateStress(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleRotateStress");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task(
    "MPM::computeRotatedStress", this, &SerialMPM::computeRotatedStress);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->needs(
        Task::OldDW, d_mpm_labels->pParticleIDLabel, matlset, Ghost::None);
      t->needs(Task::NewDW,
                  d_mpm_labels->pPolarDecompRLabel_preReloc,
                  matlset,
                  Ghost::None);

      t->modifies(d_mpm_labels->pStressLabel_preReloc, matlset);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeRotatedStress
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeRotatedStress(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing, "Doing computeRotate");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch   = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pR_new;
        ParticleVariable<Matrix3> pStress_new;
        old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
        new_dw->get(pR_new, d_mpm_labels->pPolarDecompRLabel_preReloc, pset);
        new_dw->getModifiable(
          pStress_new, d_mpm_labels->pStressLabel_preReloc, pset);

        for (auto particle : *pset) {
          pStress_new[particle] = (pR_new[particle] * pStress_new[particle]) *
                                  (pR_new[particle].Transpose());
        }
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeBasicDamage
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeBasicDamage(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* Create a task for computing the damage variables */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeBasicDamage");

  size_t numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t         = scinew Task(
    "MPM::computeBasicDamage", this, &SerialMPM::computeBasicDamage);
  for (size_t m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    // Add requires and computes for vel grad/def grad
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* d_basicDamageModel =
        mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addComputesAndRequires(
        t, mpm_matl, patches, d_mpm_labels.get());
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeBasicDamage
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeBasicDamage(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{

  printTask(patches, patches->get(0), cout_doing, "Doing computeBasicDamage");

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

    if (cout_dbg.active()) {
      cout_dbg << " Patch = " << (patches->get(0))->getID() << " Mat = " << m;
    }

    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    if (cout_dbg.active()) {
      cout_dbg << " MPM_Mat = " << mpm_matl;
    }

    // Compute basic damage
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->computeBasicDamage(
        patches, mpm_matl, old_dw, new_dw, d_mpm_labels.get());
      if (cout_dbg.active()) {
        cout_dbg << " Damage model = " << basicDamageModel;
      }
    }

    if (cout_dbg.active()) {
      cout_dbg << " Exit\n";
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleUpdateErosionParameter
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleUpdateErosionParameter(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleUpdateErosionParameter");

  Task* t = scinew Task(
    "MPM::updateErosionParameter", this, &SerialMPM::updateErosionParameter);
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* d_basicDamageModel =
        mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addRequiresLocalizationParameter(
        t, mpm_matl, patches);
    }
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addRequiresDamageParameter(t, mpm_matl, patches);
  }
  /*
  t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, Ghost::None);
  */
  t->computes(d_mpm_labels->pLocalizedMPMLabel);

  if (d_mpm_flags->d_deleteRogueParticles) {
    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
    t->computes(d_mpm_labels->numLocInCellLabel);
    t->computes(d_mpm_labels->numInCellLabel);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * updateErosionParameter
 *-----------------------------------------------------------------------*/
void
SerialMPM::updateErosionParameter(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset*,
                                  DataWarehouse* old_dw,
                                  DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing updateErosionParameter");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: material # = " << m << "\n";
      }

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: mpm_matl* = " << mpm_matl
                 << " matID = " << matID << " pset* = " << pset << "\n";
      }

      // Get the localization info
      ParticleVariable<int> isLocalized;
      new_dw->allocateAndPut(
        isLocalized, d_mpm_labels->pLocalizedMPMLabel, pset);
      for (auto particle : *pset) {
        isLocalized[particle] = 0;
      }

      // Update the localization info from basic damage model
      if (mpm_matl->doBasicDamage()) {
        Vaango::BasicDamageModel* basicDamageModel =
          mpm_matl->getBasicDamageModel();
        basicDamageModel->getLocalizationParameter(
          patch, isLocalized, matID, old_dw, new_dw);
      }

      // Update the localization info from constitutive model
      mpm_matl->getConstitutiveModel()->getDamageParameter(
        patch, isLocalized, matID, old_dw, new_dw);

      /*
      constParticleVariable<long64> pParticleID;;
      old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
      for (auto particle : *pset) {
        if (pParticleID[particle] == 562644844544) {
          std::cout << "pID=" << pParticleID[particle]
                    << " pLocalized_mpm = " << isLocalized[particle] << "\n";
        }
      }
      */

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: Got Damage Parameter"
                 << "\n";
      }

      if (d_mpm_flags->d_deleteRogueParticles) {
        // The following looks for localized particles that are isolated
        // either individually or in small groups
        // Ghost::GhostType  gac = Ghost::AroundCells;
        CCVariable<int> numLocInCell, numInCell;
        new_dw->allocateAndPut(
          numLocInCell, d_mpm_labels->numLocInCellLabel, matID, patch);
        new_dw->allocateAndPut(
          numInCell, d_mpm_labels->numInCellLabel, matID, patch);
        numLocInCell.initialize(0);
        numInCell.initialize(0);

        constParticleVariable<Point> pX;
        old_dw->get(pX, d_mpm_labels->pXLabel, pset);

        // Count the number of localized particles in each cell
        for (auto particle : *pset) {
          IntVector c;
          patch->findCell(pX[particle], c);
          numInCell[c]++;
          if (isLocalized[particle]) {
            numLocInCell[c]++;
          }
        }
      } // if d_deleteRogueParticles

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: Updated Erosion "
                 << "\n";
      }
    }

    if (cout_dbg.active()) {
      cout_dbg << "Done updateErosionParamter on patch " << patch->getID()
               << "\t MPM"
               << "\n";
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleFindRogueParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleFindRogueParticles(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpm_flags->d_deleteRogueParticles) {
    printSchedule(patches, cout_doing, "MPM::scheduleFindRogueParticles");

    Task* t = scinew Task(
      "MPM::findRogueParticles", this, &SerialMPM::findRogueParticles);
    Ghost::GhostType gac = Ghost::AroundCells;
    t->needs(Task::NewDW, d_mpm_labels->numLocInCellLabel, gac, 1);
    t->needs(Task::NewDW, d_mpm_labels->numInCellLabel, gac, 1);
    t->needs(Task::OldDW, d_mpm_labels->pXLabel, Ghost::None);
    t->needs(Task::OldDW, d_mpm_labels->pParticleIDLabel, Ghost::None);
    t->modifies(d_mpm_labels->pLocalizedMPMLabel);

    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * findRogueParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::findRogueParticles(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing findRogueParticles");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // The following looks for localized particles that are isolated
      // either individually or in small groups
      Ghost::GhostType gac = Ghost::AroundCells;
      constCCVariable<int> numLocInCell, numInCell;
      constParticleVariable<Point> pX;

      ParticleVariable<int> isLocalized;
      constParticleVariable<long64> pParticleID;

      new_dw->get(
        numLocInCell, d_mpm_labels->numLocInCellLabel, matID, patch, gac, 1);
      new_dw->get(
        numInCell, d_mpm_labels->numInCellLabel, matID, patch, gac, 1);
      old_dw->get(pX, d_mpm_labels->pXLabel, pset);
      old_dw->get(pParticleID, d_mpm_labels->pParticleIDLabel, pset);
      new_dw->getModifiable(
        isLocalized, d_mpm_labels->pLocalizedMPMLabel, pset);

      // Look at the number of localized particles in the current and
      // surrounding cells
      for (auto particle : *pset) {
        if (isLocalized[particle] == 1) {
          IntVector c;
          patch->findCell(pX[particle], c);
          int totalInCells = 0;
          for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
              for (int k = -1; k < 2; k++) {
                IntVector cell = c + IntVector(i, j, k);
                totalInCells += numInCell[cell];
              }
            }
          }
          // If the localized particles are sufficiently isolated, set
          // a flag for deletion in interpolateToParticlesAndUpdate
          if (numLocInCell[c] <= 3 && totalInCells <= 3) {
            proc0cout << "**WARNING** Particle " << pParticleID[particle]
                      << " is isolated and will be removed.\n"
                      << " cell = " << c
                      << " isLocalized = " << isLocalized[particle]
                      << " numLocIncell = " << numLocInCell[c]
                      << " totalInCells = " << totalInCells << "\n";
            isLocalized[particle] = -999;
          }

        } // if localized

        /*
        if (pParticleID[particle] == 111670263811) {
          std::cout << "pID=" << pParticleID[particle] << " " <<
        isLocalized[particle] << "\n";
        }
        */
      } // particles
    }   // matls
  }     // patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeAccStrainEnergy
 *   Compute the accumulated strain energy
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeAccStrainEnergy(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeAccStrainEnergy");

  Task* t = scinew Task(
    "MPM::computeAccStrainEnergy", this, &SerialMPM::computeAccStrainEnergy);
  t->needs(Task::OldDW, d_mpm_labels->AccStrainEnergyLabel);
  t->needs(Task::NewDW, d_mpm_labels->StrainEnergyLabel);
  t->computes(d_mpm_labels->AccStrainEnergyLabel);
  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeAccStrainEnergy
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeAccStrainEnergy(const ProcessorGroup*,
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

/*!----------------------------------------------------------------------
 * scheduleFinalParticleUpdate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleFinalParticleUpdate(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleFinalParticleUpdate");

  Task* t = scinew Task(
    "MPM::finalParticleUpdate", this, &SerialMPM::finalParticleUpdate);

  t->needs(Task::OldDW, d_mpm_labels->delTLabel);

  Ghost::GhostType gnone = Ghost::None;
  t->needs(Task::NewDW, d_mpm_labels->pdTdtLabel, gnone);
  t->needs(Task::NewDW, d_mpm_labels->pLocalizedMPMLabel_preReloc, gnone);
  t->needs(Task::NewDW, d_mpm_labels->pMassLabel_preReloc, gnone);

  t->modifies(d_mpm_labels->pTemperatureLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * finalParticleUpdate
 *-----------------------------------------------------------------------*/
void
SerialMPM::finalParticleUpdate(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing finalParticleUpdate");

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<int> pLocalized;
      constParticleVariable<double> pdTdt, pMass_new;
      ParticleVariable<double> pTemp_new;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      auto* delset         = scinew ParticleSubset(0, matID, patch);

      new_dw->get(pdTdt, d_mpm_labels->pdTdtLabel, pset);
      new_dw->get(pMass_new, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->get(pLocalized, d_mpm_labels->pLocalizedMPMLabel_preReloc, pset);

      new_dw->getModifiable(
        pTemp_new, d_mpm_labels->pTemperatureLabel_preReloc, pset);

      // Loop over particles
      for (auto idx : *pset) {
        pTemp_new[idx] += pdTdt[idx] * delT;

        if ((pMass_new[idx] <= d_mpm_flags->d_minPartMass && !mpm_matl->isDiscrete()) ||
            pTemp_new[idx] < 0. || (pLocalized[idx] == -999)) {
          if (d_mpm_flags->d_erosionAlgorithm != "none") {
            delset->addParticle(idx);
          }
        }

      } // particles

      new_dw->deleteParticles(delset);

    } // materials
  }   // patches
}

/*!----------------------------------------------------------------------
 * printParticleLabels
 *-----------------------------------------------------------------------*/
void
SerialMPM::printParticleLabels(std::vector<const VarLabel*> labels,
                               DataWarehouse* dw,
                               int matID,
                               const Patch* patch)
{
  for (auto label : labels) {
    if (dw->exists(label, matID, patch)) {
      std::cout << label->getName() << " does exists"
                << "\n";
    } else {
      std::cout << label->getName() << " does NOT exists"
                << "\n";
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleInsertParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInsertParticles(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpm_flags->d_insertParticles) {
    printSchedule(patches, cout_doing, "MPM::scheduleInsertParticles");

    Task* t =
      scinew Task("MPM::insertParticles", this, &SerialMPM::insertParticles);

    t->needs(Task::OldDW, d_mpm_labels->simulationTimeLabel);
    t->needs(Task::OldDW, d_mpm_labels->delTLabel);

    t->modifies(d_mpm_labels->pXLabel_preReloc);
    t->modifies(d_mpm_labels->pVelocityLabel_preReloc);
    t->needs(Task::OldDW, d_mpm_labels->pColorLabel, Ghost::None);

    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * insertParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::insertParticles(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing insertParticles");

    // Get current time and timestep size
    simTime_vartype simTimeVar;
    old_dw->get(simTimeVar, d_mpm_labels->simulationTimeLabel);
    double time = simTimeVar;

    delt_vartype delT;
    old_dw->get(delT, d_mpm_labels->delTLabel, getLevel(patches));

    // activate materials if it is their time
    size_t numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (size_t m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      if (time >= mpm_matl->getActivationTime()) {
        mpm_matl->setActive(true);
      }
    }

    int index = -999;

    for (auto i = 0u; i < d_IPTimes.size(); i++) {
      if (time + delT > d_IPTimes[i] && time <= d_IPTimes[i]) {
        index = i;
        if (index >= 0) {
          size_t numMPMMatls = d_materialManager->getNumMaterials("MPM");
          for (size_t m = 0; m < numMPMMatls; m++) {
            MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(
              d_materialManager->getMaterial("MPM", m));
            int matID            = mpm_matl->getDWIndex();
            ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

            // Get the arrays of particle values to be changed
            ParticleVariable<Point> pX;
            ParticleVariable<Vector> pVelocity;
            constParticleVariable<double> pColor;

            old_dw->get(pColor, d_mpm_labels->pColorLabel, pset);
            new_dw->getModifiable(pX, d_mpm_labels->pXLabel_preReloc, pset);
            new_dw->getModifiable(
              pVelocity, d_mpm_labels->pVelocityLabel_preReloc, pset);

            int numParticles = pset->end() - pset->begin();
            std::cout << "Insertion: Patch " << p << " now contains "
                      << numParticles << " particles\n";
            // Loop over particles here
            for (auto idx : *pset) {
              if (pColor[idx] == d_IPColor[index]) {
                pVelocity[idx] = d_IPVelNew[index];
                pX[idx] += d_IPTranslate[index];
              }
            }
          }
        }
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleAddParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleAddParticles(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleAddParticles");

  Task* t = scinew Task("MPM::addParticles", this, &SerialMPM::addParticles);

  auto* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->modifies(d_mpm_labels->pParticleIDLabel_preReloc);
  t->modifies(d_mpm_labels->pXLabel_preReloc);
  t->modifies(d_mpm_labels->pVolumeLabel_preReloc);
  t->modifies(d_mpm_labels->pVelocityLabel_preReloc);
  t->modifies(d_mpm_labels->pMassLabel_preReloc);
  t->modifies(d_mpm_labels->pSizeLabel_preReloc);
  t->modifies(d_mpm_labels->pDispLabel_preReloc);
  t->modifies(d_mpm_labels->pStressLabel_preReloc);
  t->modifies(d_mpm_labels->pColorLabel_preReloc);
  t->modifies(d_mpm_labels->pLocalizedMPMLabel_preReloc);
  t->modifies(d_mpm_labels->pExtForceLabel_preReloc);
  t->modifies(d_mpm_labels->pTemperatureLabel_preReloc);
  t->modifies(d_mpm_labels->pTempPreviousLabel_preReloc);
  t->modifies(d_mpm_labels->pDefGradLabel_preReloc);
  t->modifies(d_mpm_labels->pRefinedLabel_preReloc);
  t->modifies(d_mpm_labels->pVelGradLabel_preReloc);

  t->needs(
    Task::OldDW, d_mpm_labels->pCellNAPIDLabel, zeroth_matl, Ghost::None);
  t->computes(d_mpm_labels->pCellNAPIDLabel, zeroth_matl);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * addParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::addParticles(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset*,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing addParticles");
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    // Carry forward CellNAPID
    constCCVariable<short int> NAPID;
    CCVariable<short int> NAPID_new;
    Ghost::GhostType gnone = Ghost::None;
    old_dw->get(NAPID, d_mpm_labels->pCellNAPIDLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NAPID_new, d_mpm_labels->pCellNAPIDLabel, 0, patch);
    NAPID_new.copyData(NAPID);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      ParticleVariable<Point> pX;
      ParticleVariable<Matrix3> pF, pSize, pStress, pvelgrad;
      ParticleVariable<long64> pParticleID;
      ParticleVariable<double> pVolume, pMass, ptemp, ptempP, pcolor;
      ParticleVariable<Vector> pVelocity, pextforce, pDisp;
      ParticleVariable<int> pref, ploc;
      new_dw->getModifiable(pX, d_mpm_labels->pXLabel_preReloc, pset);
      new_dw->getModifiable(
        pParticleID, d_mpm_labels->pParticleIDLabel_preReloc, pset);
      new_dw->getModifiable(pMass, d_mpm_labels->pMassLabel_preReloc, pset);
      new_dw->getModifiable(pSize, d_mpm_labels->pSizeLabel_preReloc, pset);
      new_dw->getModifiable(pDisp, d_mpm_labels->pDispLabel_preReloc, pset);
      new_dw->getModifiable(pStress, d_mpm_labels->pStressLabel_preReloc, pset);
      new_dw->getModifiable(pcolor, d_mpm_labels->pColorLabel_preReloc, pset);
      new_dw->getModifiable(pVolume, d_mpm_labels->pVolumeLabel_preReloc, pset);
      new_dw->getModifiable(
        pVelocity, d_mpm_labels->pVelocityLabel_preReloc, pset);
      new_dw->getModifiable(
        pextforce, d_mpm_labels->pExtForceLabel_preReloc, pset);
      new_dw->getModifiable(
        ptemp, d_mpm_labels->pTemperatureLabel_preReloc, pset);
      new_dw->getModifiable(
        ptempP, d_mpm_labels->pTempPreviousLabel_preReloc, pset);
      new_dw->getModifiable(pref, d_mpm_labels->pRefinedLabel_preReloc, pset);
      new_dw->getModifiable(
        ploc, d_mpm_labels->pLocalizedMPMLabel_preReloc, pset);
      new_dw->getModifiable(
        pvelgrad, d_mpm_labels->pVelGradLabel_preReloc, pset);
      new_dw->getModifiable(pF, d_mpm_labels->pDefGradLabel_preReloc, pset);

      int numNewPartNeeded = 0;
      // Put refinement criteria here
      const unsigned int origNParticles = pset->addParticles(0);
      for (unsigned int pp = 0; pp < origNParticles; ++pp) {
        if (pref[pp] == 0 && pStress[pp].Norm() > 1) {
          pref[pp] = 2;
          numNewPartNeeded++;
        }
      }
      numNewPartNeeded *= 8;

      const unsigned int oldNumPar = pset->addParticles(numNewPartNeeded);

      ParticleVariable<Point> pXtmp;
      ParticleVariable<Matrix3> pFtmp, pSizetmp, pstrstmp, pvgradtmp;
      ParticleVariable<long64> pParticleIDtmp;
      ParticleVariable<double> pVoltmp, pMasstmp, ptemptmp, ptempPtmp,
        pcolortmp;
      ParticleVariable<Vector> pveltmp, pextFtmp, pDisptmp;
      ParticleVariable<int> preftmp, ploctmp;
      new_dw->allocateTemporary(pParticleIDtmp, pset);
      new_dw->allocateTemporary(pXtmp, pset);
      new_dw->allocateTemporary(pVoltmp, pset);
      new_dw->allocateTemporary(pveltmp, pset);
      new_dw->allocateTemporary(pextFtmp, pset);
      new_dw->allocateTemporary(ptemptmp, pset);
      new_dw->allocateTemporary(ptempPtmp, pset);
      new_dw->allocateTemporary(pFtmp, pset);
      new_dw->allocateTemporary(pSizetmp, pset);
      new_dw->allocateTemporary(pDisptmp, pset);
      new_dw->allocateTemporary(pstrstmp, pset);
      new_dw->allocateTemporary(pcolortmp, pset);
      new_dw->allocateTemporary(pMasstmp, pset);
      new_dw->allocateTemporary(preftmp, pset);
      new_dw->allocateTemporary(ploctmp, pset);
      new_dw->allocateTemporary(pvgradtmp, pset);

      // copy data from old variables for particle IDs and the position
      // std::vector
      for (unsigned int pp = 0; pp < oldNumPar; ++pp) {
        pParticleIDtmp[pp] = pParticleID[p];
        pXtmp[pp]          = pX[pp];
        pVoltmp[pp]        = pVolume[pp];
        pveltmp[pp]        = pVelocity[pp];
        pextFtmp[pp]       = pextforce[pp];
        ptemptmp[pp]       = ptemp[pp];
        ptempPtmp[pp]      = ptempP[pp];
        pFtmp[pp]          = pF[pp];
        pSizetmp[pp]       = pSize[pp];
        pDisptmp[pp]       = pDisp[pp];
        pstrstmp[pp]       = pStress[pp];
        pcolortmp[pp]      = pcolor[pp];
        pMasstmp[pp]       = pMass[pp];
        preftmp[pp]        = pref[pp];
        ploctmp[pp]        = ploc[pp];
        pvgradtmp[pp]      = pvelgrad[pp];
      }

      Vector dx     = patch->dCell();
      int numRefPar = 0;
      for (unsigned int idx = 0; idx < oldNumPar; ++idx) {
        if (pref[idx] == 2) {
          std::vector<Point> new_part_pos;

          Matrix3 dsize = (pF[idx] * pSize[idx] *
                           Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));

          // Find std::vectors to new particle locations, based on particle size
          // and deformation (patterned after CPDI interpolator code)
          Vector r[4];
          r[0] = Vector(-dsize(0, 0) - dsize(0, 1) + dsize(0, 2),
                        -dsize(1, 0) - dsize(1, 1) + dsize(1, 2),
                        -dsize(2, 0) - dsize(2, 1) + dsize(2, 2)) *
                 0.25;
          r[1] = Vector(dsize(0, 0) - dsize(0, 1) + dsize(0, 2),
                        dsize(1, 0) - dsize(1, 1) + dsize(1, 2),
                        dsize(2, 0) - dsize(2, 1) + dsize(2, 2)) *
                 0.25;
          r[2] = Vector(dsize(0, 0) + dsize(0, 1) + dsize(0, 2),
                        dsize(1, 0) + dsize(1, 1) + dsize(1, 2),
                        dsize(2, 0) + dsize(2, 1) + dsize(2, 2)) *
                 0.25;
          r[3] = Vector(-dsize(0, 0) + dsize(0, 1) + dsize(0, 2),
                        -dsize(1, 0) + dsize(1, 1) + dsize(1, 2),
                        -dsize(2, 0) + dsize(2, 1) + dsize(2, 2)) *
                 0.25;

          new_part_pos.push_back(pX[idx] + r[0]);
          new_part_pos.push_back(pX[idx] + r[1]);
          new_part_pos.push_back(pX[idx] + r[2]);
          new_part_pos.push_back(pX[idx] + r[3]);
          new_part_pos.push_back(pX[idx] - r[0]);
          new_part_pos.push_back(pX[idx] - r[1]);
          new_part_pos.push_back(pX[idx] - r[2]);
          new_part_pos.push_back(pX[idx] - r[3]);

          //        new_part_pos.push_back(pX[idx]+Vector(dxp,dxp,dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(-dxp,-dxp,-dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(dxp,dxp,-dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(dxp,-dxp,dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(-dxp,dxp,dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(dxp,-dxp,-dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(-dxp,-dxp,dxp));
          //        new_part_pos.push_back(pX[idx]+Vector(-dxp,dxp,-dxp));
          std::cout << "new_part_pos = " << new_part_pos[0] << "\n";

          for (int i = 0; i < 8; i++) {
            IntVector c;
            patch->findCell(new_part_pos[i], c);

            long64 cellID = ((long64)c.x() << 16) | ((long64)c.y() << 32) |
                            ((long64)c.z() << 48);

            short int& myCellNAPID = NAPID_new[c];
            int new_index;
            if (i == 0) {
              new_index = idx;
            } else {
              new_index = oldNumPar + 8 * numRefPar + i;
            }
            pParticleIDtmp[new_index] = (cellID | (long64)myCellNAPID);
            pXtmp[new_index]          = new_part_pos[i];
            pVoltmp[new_index]        = .125 * pVolume[idx];
            pMasstmp[new_index]       = .125 * pMass[idx];
            pveltmp[new_index]        = pVelocity[idx];
            pextFtmp[new_index]       = pextforce[idx];
            pFtmp[new_index]          = pF[idx];
            pSizetmp[new_index]       = 0.5 * pSize[idx];
            pDisptmp[new_index]       = pDisp[idx];
            pstrstmp[new_index]       = pStress[idx];
            pcolortmp[new_index]      = pcolor[idx];
            ptemptmp[new_index]       = ptemp[idx];
            ptempPtmp[new_index]      = ptempP[idx];
            preftmp[new_index]        = 1;
            ploctmp[new_index]        = ploc[idx];
            pvgradtmp[new_index]      = pvelgrad[idx];
            NAPID_new[c]++;
          }
          numRefPar++;
        } // if particle flagged for refinement
      }   // for particles

      // put back temporary data
      new_dw->put(
        pParticleIDtmp, d_mpm_labels->pParticleIDLabel_preReloc, true);
      new_dw->put(pXtmp, d_mpm_labels->pXLabel_preReloc, true);
      new_dw->put(pVoltmp, d_mpm_labels->pVolumeLabel_preReloc, true);
      new_dw->put(pveltmp, d_mpm_labels->pVelocityLabel_preReloc, true);
      new_dw->put(pextFtmp, d_mpm_labels->pExtForceLabel_preReloc, true);
      new_dw->put(pMasstmp, d_mpm_labels->pMassLabel_preReloc, true);
      new_dw->put(ptemptmp, d_mpm_labels->pTemperatureLabel_preReloc, true);
      new_dw->put(ptempPtmp, d_mpm_labels->pTempPreviousLabel_preReloc, true);
      new_dw->put(pSizetmp, d_mpm_labels->pSizeLabel_preReloc, true);
      new_dw->put(pDisptmp, d_mpm_labels->pDispLabel_preReloc, true);
      new_dw->put(pstrstmp, d_mpm_labels->pStressLabel_preReloc, true);
      new_dw->put(pcolortmp, d_mpm_labels->pColorLabel_preReloc, true);
      new_dw->put(pFtmp, d_mpm_labels->pDefGradLabel_preReloc, true);
      new_dw->put(preftmp, d_mpm_labels->pRefinedLabel_preReloc, true);
      new_dw->put(ploctmp, d_mpm_labels->pLocalizedMPMLabel_preReloc, true);
      new_dw->put(pvgradtmp, d_mpm_labels->pVelGradLabel_preReloc, true);
    } // for matls
  }   // for patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeParticleScaleFactor
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeParticleScaleFactor(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)

{
  if (!d_mpm_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                 getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleScaleFactor");

  Task* t = scinew Task("MPM::computeParticleScaleFactor",
                        this,
                        &SerialMPM::computeParticleScaleFactor);

  t->needs(Task::OldDW, d_mpm_labels->pSizeLabel, Ghost::None);
  t->needs(Task::NewDW, d_mpm_labels->pDefGradLabel_preReloc, Ghost::None);
  t->computes(d_mpm_labels->pScaleFactorLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeParticleScaleFactor
 *   This task computes the particles initial physical size, to be used
 *   in scaling particles for the deformed particle vis feature
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeParticleScaleFactor(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleScaleFactor");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID            = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pScaleFactor;
      old_dw->get(pSize, d_mpm_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(
        pScaleFactor, d_mpm_labels->pScaleFactorLabel_preReloc, pset);

      if (d_output->isOutputTimestep()) {
        Vector dx = patch->dCell();

        if (d_mpm_flags->d_interpolatorType != "cpdi" &&
            d_mpm_flags->d_interpolatorType != "cpti") {
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pDefGrad, d_mpm_labels->pDefGradLabel_preReloc, pset);
          for (auto pidx : *pset) {
            pScaleFactor[pidx] =
              ((pDefGrad[pidx] * pSize[pidx]) *
               Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));
          }
        } else {
          for (int idx : *pset) {
            pScaleFactor[idx] =
              (pSize[idx] * Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));

          } // for particles
        }
      } // isOutputTimestep
    }   // matls
  }     // patches
}

void
SerialMPM::scheduleParticleRelocation(SchedulerP& sched,
                                      const LevelP& level,
                                      [[maybe_unused]] const PatchSet* patches,
                                      const MaterialSet* matls)
{
  //  Unmodified labels and matls subset
  std::vector<std::vector<const VarLabel*>> old_labels =
    d_particleState_preReloc;
  std::vector<std::vector<const VarLabel*>> new_labels = d_particleState;

  const MaterialSubset* old_mss = matls->getSubset(0);

  auto* new_mss = scinew MaterialSubset();
  new_mss->addReference();
  new_mss->addSubset(old_mss);

  // If needed concatenate the labels and matls that are passed into
  // the ParticleRelocate
  d_cohesiveZoneTasks->scheduleParticleRelocation(
    new_mss, old_labels, new_labels);

  //  create a new material set containing the
  //  the updated matlSubset.
  auto* newMatlSet = scinew MaterialSet();
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

/*!----------------------------------------------------------------------
 * scheduleRefine
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  printSchedule(patches, cout_doing, "MPM::scheduleRefine");
  Task* t = scinew Task("SerialMPM::refine", this, &SerialMPM::refine);

  t->computes(d_mpm_labels->pXLabel);
  t->computes(d_mpm_labels->p_qLabel);
  t->computes(d_mpm_labels->pDispLabel);
  t->computes(d_mpm_labels->pMassLabel);
  t->computes(d_mpm_labels->pVolumeLabel);
  t->computes(d_mpm_labels->pTemperatureLabel);
  t->computes(d_mpm_labels->pTempPreviousLabel); // for thermal stress analysis
  t->computes(d_mpm_labels->pdTdtLabel);
  t->computes(d_mpm_labels->pVelocityLabel);
  t->computes(d_mpm_labels->pExternalForceLabel);
  t->computes(d_mpm_labels->pParticleIDLabel);
  t->computes(d_mpm_labels->pStressLabel);
  t->computes(d_mpm_labels->pSizeLabel);

  // DEM labels
  if (d_mpm_flags->d_enableDEM) {
    t->computes(d_mpm_labels->pX0Label);
    t->computes(d_mpm_labels->pRigidBodyIDLabel);
    t->computes(d_mpm_labels->pAngularVelocityLabel);
    t->computes(d_mpm_labels->pTorqueLabel);
    t->computes(d_mpm_labels->pOrientationLabel);
    t->computes(d_mpm_labels->pRadiusLabel);
    t->computes(d_mpm_labels->pInertiaTensorLabel);
  }

  t->computes(d_mpm_labels->NC_CCweightLabel);
  t->computes(d_mpm_labels->delTLabel, getLevel(patches));

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

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    // For deformation gradient computer
    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    // For constitutive models
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Basic damage model related stuff
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->addInitialComputesAndRequires(
        t, mpm_matl, patches, d_mpm_labels.get());
    }
  }

  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
}

/*!----------------------------------------------------------------------
 * refine
 *-----------------------------------------------------------------------*/
void
SerialMPM::refine(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset* /*matls*/,
                  DataWarehouse*,
                  DataWarehouse* new_dw)
{
  // just create a particle subset if one doesn't exist
  // and initialize NC_CCweights

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing refine");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    // First do NC_CCweight
    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(
      NC_CCweight, d_mpm_labels->NC_CCweightLabel, 0, patch);
    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and
    //   double NC_CCweight
    NC_CCweight.initialize(0.125);
    std::vector<Patch::FaceType>::const_iterator iter;
    std::vector<Patch::FaceType> bf;
    patch->getBoundaryFaces(bf);

    for (iter = bf.begin(); iter != bf.end(); ++iter) {
      Patch::FaceType face = *iter;
      int mat_id           = 0;
      if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {

        for (CellIterator iter = patch->getFaceIterator(face, Patch::FaceNodes);
             !iter.done();
             iter++) {
          NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
        }
      }
    }

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();

      if (cout_doing.active()) {
        cout_doing << "Doing refine on patch " << patch->getID()
                   << " material # = " << matID << "\n";
      }

      // this is a new patch, so create empty particle variables.
      if (!new_dw->haveParticleSubset(matID, patch)) {
        ParticleSubset* pset = new_dw->createParticleSubset(0, matID, patch);

        // Create arrays for the particle data
        ParticleVariable<Point> pX;
        ParticleVariable<double> pMass, pVolume, pTemperature;
        ParticleVariable<Vector> pVelocity, pExternalForce, pDisp;
        ParticleVariable<Matrix3> pSize;
        ParticleVariable<double> pTempPrev, p_q;
        ParticleVariable<int> pLoadCurve;
        ParticleVariable<long64> pID;
        ParticleVariable<Matrix3> pDeform, pStress;

        new_dw->allocateAndPut(pX, d_mpm_labels->pXLabel, pset);
        new_dw->allocateAndPut(p_q, d_mpm_labels->p_qLabel, pset);
        new_dw->allocateAndPut(pMass, d_mpm_labels->pMassLabel, pset);
        new_dw->allocateAndPut(pVolume, d_mpm_labels->pVolumeLabel, pset);
        new_dw->allocateAndPut(pVelocity, d_mpm_labels->pVelocityLabel, pset);
        new_dw->allocateAndPut(
          pTemperature, d_mpm_labels->pTemperatureLabel, pset);
        new_dw->allocateAndPut(
          pTempPrev, d_mpm_labels->pTempPreviousLabel, pset);
        new_dw->allocateAndPut(
          pExternalForce, d_mpm_labels->pExternalForceLabel, pset);
        new_dw->allocateAndPut(pID, d_mpm_labels->pParticleIDLabel, pset);
        new_dw->allocateAndPut(pDisp, d_mpm_labels->pDispLabel, pset);

        // DEM labels
        ParticleVariable<Point> pX0;
        ParticleVariable<long64> pRigidBodyID;
        ParticleVariable<Vector> pAngularVelocity, pTorque;
        ParticleVariable<Matrix3> pOrientation, pInertiaTensor;
        ParticleVariable<double> pRadius;
        if (d_mpm_flags->d_enableDEM) {
          new_dw->allocateAndPut(pX0, d_mpm_labels->pX0Label, pset);
          pX0.copyData(pX);
          new_dw->allocateAndPut(
            pRigidBodyID, d_mpm_labels->pRigidBodyIDLabel, pset);
          new_dw->allocateAndPut(
            pAngularVelocity, d_mpm_labels->pAngularVelocityLabel, pset);
          new_dw->allocateAndPut(pTorque, d_mpm_labels->pTorqueLabel, pset);
          new_dw->allocateAndPut(
            pOrientation, d_mpm_labels->pOrientationLabel, pset);
          new_dw->allocateAndPut(pRadius, d_mpm_labels->pRadiusLabel, pset);
          new_dw->allocateAndPut(
            pInertiaTensor, d_mpm_labels->pInertiaTensorLabel, pset);
        }

        if (d_mpm_flags->d_useLoadCurves) {
          new_dw->allocateAndPut(
            pLoadCurve, d_mpm_labels->pLoadCurveIDLabel, pset);
        }
        new_dw->allocateAndPut(pSize, d_mpm_labels->pSizeLabel, pset);

        // Init deformation gradient
        d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

        // Init constitutive model
        mpm_matl->getConstitutiveModel()->initializeCMData(
          patch, mpm_matl, new_dw);

        // Initialize basic damage models
        if (mpm_matl->doBasicDamage()) {
          mpm_matl->getBasicDamageModel()->initializeDamageData(
            patch, mpm_matl, new_dw, d_mpm_labels.get());
        }

#if 0
        if(d_mpm_flags->d_withColor) {
          ParticleVariable<double> pcolor;
          int index = mpm_matl->getDWIndex();
          ParticleSubset* pset = new_dw->getParticleSubset(index, patch);
          setParticleDefault<double>(pcolor, d_mpm_labels->pColorLabel, pset, new_dw, 0.0);
        }
#endif
      }
    }
  }

} // end refine()

/*!----------------------------------------------------------------------
 * scheduleRefineInterface
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleRefineInterface(const LevelP& /*fineLevel*/,
                                   SchedulerP& /*scheduler*/,
                                   bool,
                                   bool)
{
  //  do nothing for now
}

/*!----------------------------------------------------------------------
 * scheduleCoarsen
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleCoarsen(const LevelP& /*coarseLevel*/, SchedulerP& /*sched*/)
{
  // do nothing for now
}

/*!----------------------------------------------------------------------
 * scheduleErrorEstimate
 *   Schedule to mark d_mpm_flags for AMR regridding
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  // main way is to count particles, but for now we only want particles on
  // the finest level.  Thus to schedule cells for regridding during the
  // execution, we'll coarsen the flagged cells (see coarsen).

  if (amr_doing.active()) {
    amr_doing << "SerialMPM::scheduleErrorEstimate on level "
              << coarseLevel->getIndex() << '\n';
  }

  // The simulation controller should not schedule it every time step
  Task* task = scinew Task("errorEstimate", this, &SerialMPM::errorEstimate);

  // if the finest level, compute flagged cells
  if (coarseLevel->getIndex() == coarseLevel->getGrid()->numLevels() - 1) {
    task->needs(Task::NewDW, d_mpm_labels->pXLabel, Ghost::AroundCells, 0);
  } else {
    task->needs(Task::NewDW,
                   d_regridder->getRefineFlagLabel(),
                   nullptr,
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
  sched->addTask(
    task, coarseLevel->eachPatch(), d_materialManager->allMaterials("MPM"));
}

/*!----------------------------------------------------------------------
 * errorEstimate
 *-----------------------------------------------------------------------*/
void
SerialMPM::errorEstimate(const ProcessorGroup* group,
                         const PatchSubset* coarsePatches,
                         const MaterialSubset* matls,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  if (coarseLevel->getIndex() == coarseLevel->getGrid()->numLevels() - 1) {
    // on finest level, we do the same thing as initialErrorEstimate, so call it
    initialErrorEstimate(group, coarsePatches, matls, old_dw, new_dw);
  } else {
    // coarsen the errorflag.
    const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();

    for (int p = 0; p < coarsePatches->size(); p++) {
      const Patch* coarsePatch = coarsePatches->get(p);
      printTask(coarsePatches, coarsePatch, cout_doing, "Doing errorEstimate");

      CCVariable<int> refineFlag;
      PerPatch<PatchFlagP> refinePatchFlag;

      new_dw->getModifiable(
        refineFlag, d_regridder->getRefineFlagLabel(), 0, coarsePatch);
      new_dw->get(refinePatchFlag,
                  d_regridder->getRefinePatchFlagLabel(),
                  0,
                  coarsePatch);

      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);

      // coarsen the fineLevel flag
      for (auto finePatch : finePatches) {
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
        for (CellIterator iter(fl, fh); !iter.done(); iter++) {

          IntVector coarseCell(fineLevel->mapCellToCoarser(*iter));

          if (fineErrorFlag[*iter]) {
            refineFlag[coarseCell] = 1;
            refinePatch->set();
          }
        }
      } // fine patch loop
    }   // coarse patch loop
  }
}

/*!----------------------------------------------------------------------
 * scheduleInitialErrorEstimate
 *   Schedule to mark initial d_mpm_flags for AMR regridding
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                        SchedulerP& sched)
{
  scheduleErrorEstimate(coarseLevel, sched);
}

/*!----------------------------------------------------------------------
 * initialErrorEstimate
 *-----------------------------------------------------------------------*/
void
SerialMPM::initialErrorEstimate(const ProcessorGroup*,
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
    new_dw->getModifiable(
      refineFlag, d_regridder->getRefineFlagLabel(), 0, patch);
    new_dw->get(
      refinePatchFlag, d_regridder->getRefinePatchFlagLabel(), 0, patch);

    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int matID = mpm_matl->getDWIndex();
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
      constParticleVariable<Point> pX;
      new_dw->get(pX, d_mpm_labels->pXLabel, pset);

      for (int& iter : *pset) {
        refineFlag[patch->getLevel()->getCellIndex(pX[iter])] = true;
        refinePatch->set();
      }
    }
  }
}
/*!----------------------------------------------------------------------
 * scheduleSwitchTest
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  if (d_switchCriteria) {
    d_switchCriteria->scheduleSwitchTest(level, sched);
  }
}

/*!----------------------------------------------------------------------
 * Set particle default
 *-----------------------------------------------------------------------*/
template<typename T>
void
SerialMPM::setParticleDefault(ParticleVariable<T>& pvar,
                              const VarLabel* label,
                              ParticleSubset* pset,
                              DataWarehouse* new_dw,
                              const T& val)
{
  new_dw->allocateAndPut(pvar, label, pset);
  for (auto idx : *pset) {
    pvar[idx] = val;
  }
}

namespace Uintah {
template void
SerialMPM::setParticleDefault<>(ParticleVariable<double>& pvar,
                                const VarLabel* label,
                                ParticleSubset* pset,
                                DataWarehouse* new_dw,
                                const double& val);
}
