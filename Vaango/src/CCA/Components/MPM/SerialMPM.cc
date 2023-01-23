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
#include <CCA/Components/MPM/SerialMPM.h>

#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/HeatConduction/HeatConduction.h>
#include <CCA/Components/MPM/MMS/MMS.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>
#include <CCA/Components/OnTheFlyAnalysis/AnalysisModuleFactory.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
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
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Math/CubicPolyRoots.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <Eigen/Dense>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

// #define XPIC2_UPDATE
#define CHECK_PARTICLE_DELETION
// #define TIME_COMPUTE_STRESS
#define CHECK_ISFINITE
// #define DEBUG_WITH_PARTICLE_ID
//  constexpr long64 testParticleID = testParticleID;

using namespace Uintah;

//__________________________________
//  To turn on debug d_mpmFlags
//  csh/tcsh : setenv SCI_DEBUG "MPM:+,SerialMPM:+".....
//  bash     : export SCI_DEBUG="MPM:+,SerialMPM:+" )
//  default is OFF

static DebugStream cout_doing("MPM", false);
static DebugStream cout_dbg("SerialMPM", false);
static DebugStream cout_convert("MPMConv", false);
static DebugStream cout_heat("MPMHeat", false);
static DebugStream amr_doing("AMRMPM", false);
static DebugStream cout_damage("Damage", false);

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern Mutex cerrLock;

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

SerialMPM::SerialMPM(const ProcessorGroup* myworld,
                     const MaterialManagerP materialManager)
  : SimulationCommon(myworld, materialManager)
  , MPMCommon(materialManager)
{
  d_mpmLabels = std::make_unique<MPMLabel>();
  d_mpmFlags  = std::make_unique<MPMFlags>(myworld);
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
                        MaterialManagerP& mat_manager)
{
  cout_doing << "Doing problemSetup\t\t\t\t\t MPM"
             << "\n";
  d_materialManager = mat_manager;
  dynamic_cast<Scheduler*>(getPort("scheduler"))
    ->setPositionVar(d_mpmLabels->pXLabel);

  d_dataArchiver = dynamic_cast<Output*>(getPort("output"));
  if (!d_dataArchiver) {
    throw InternalError("MPM:couldn't get output port", __FILE__, __LINE__);
  }

  ProblemSpecP restart_mat_ps = nullptr;
  ProblemSpecP prob_spec_mat_ps =
    prob_spec->findBlockWithOutAttribute("MaterialProperties");

  if (prob_spec_mat_ps) {
    restart_mat_ps = prob_spec;
  } else if (restart_prob_spec) {
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

  // Read all MPM d_mpmFlags (look in MPMFlags.cc)
  d_mpmFlags->readMPMFlags(restart_mat_ps, d_dataArchiver);
  if (d_mpmFlags->d_integratorType == "implicit") {
    throw ProblemSetupException("Can't use implicit integration with -mpm",
                                __FILE__,
                                __LINE__);
  }

  // convert text representation of face into FaceType
  for (auto faceStr : d_mpmFlags->d_boundaryTractionFaceStrings) {
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

  // read in AMR d_mpmFlags from the main ups file
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  if (amr_ps) {
    ProblemSpecP mpm_amr_ps = amr_ps->findBlock("MPM");
    if (!mpm_amr_ps) {
      std::ostringstream warn;
      warn << "ERROR:MPM:\n missing MPM section in the AMR section of the "
              "input file\n";
      throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
    }

    mpm_amr_ps->getWithDefault("min_grid_level", d_mpmFlags->d_minGridLevel, 0);
    mpm_amr_ps->getWithDefault("max_grid_level",
                               d_mpmFlags->d_maxGridLevel,
                               1000);
  }

  if (d_mpmFlags->d_8or27 == 8) {
    d_numGhostParticles = 1;
    d_numGhostNodes     = 1;
  } else {
    d_numGhostParticles = 2;
    d_numGhostNodes     = 2;
  }

  if (d_mpmFlags->d_prescribeDeformation) {
    readPrescribedDeformations(d_mpmFlags->d_prescribedDeformationFile);
  }
  if (d_mpmFlags->d_insertParticles) {
    readInsertParticlesFile(d_mpmFlags->d_insertParticlesFile);
  }

  d_materialManager->setParticleGhostLayer(Ghost::AroundNodes,
                                           d_numGhostParticles);

  MPMPhysicalBCFactory::create(restart_mat_ps, grid, d_mpmFlags);

  contactModel = ContactFactory::create(UintahParallelComponent::d_myworld,
                                        restart_mat_ps,
                                        d_materialManager,
                                        d_mpmLabels,
                                        d_mpmFlags);
  thermalContactModel = ThermalContactFactory::create(restart_mat_ps,
                                                      d_materialManager,
                                                      d_mpmLabels,
                                                      d_mpmFlags);

  heatConductionModel =
    scinew HeatConduction(d_materialManager, d_mpmLabels, d_mpmFlags);

  // Creates MPM material w/ constitutive models and damage models
  materialProblemSetup(restart_mat_ps, grid, d_materialManager, d_mpmFlags);

  cohesiveZoneProblemSetup(restart_mat_ps, d_materialManager, d_mpmFlags);

  // Create deformation gradient computer
  d_defGradComputer =
    scinew DeformationGradientComputer(flags, d_materialManager);

  // create analysis modules
  if (!d_mpmFlags->d_withICE) { // mpmice handles this
    d_analysisModules =
      AnalysisModuleFactory::create(prob_spec, d_materialManager, dataArchiver);

    for (auto& module : d_analysisModules) {
      module->problemSetup(prob_spec,
                           restart_prob_spec,
                           grid,
                           d_materialManager);
    }
  }
}

//__________________________________
//  create the switching criteria port
d_switchCriteria = dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));

if (d_switchCriteria) {
  d_switchCriteria->problemSetup(restart_mat_ps,
                                 restart_prob_spec,
                                 d_materialManager);
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
        d_prescribedF.push_back(
          Matrix3(F11, F12, F13, F21, F22, F23, F31, F32, F33));
        d_prescribedAngle.push_back(Theta);
        d_prescribedRotationAxis.push_back(Vector(a1, a2, a3));
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
        d_IPTranslate.push_back(Vector(transx, transy, transz));
        d_IPVelNew.push_back(Vector(v_new_x, v_new_y, v_new_z));
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

  ProblemSpecP d_mpmFlags_ps = root->appendChild("MPM");
  d_mpmFlags->outputProblemSpec(flags_ps);

  ProblemSpecP mat_ps = 0;
  mat_ps              = root->findBlockWithOutAttribute("MaterialProperties");

  if (mat_ps == 0) {
    mat_ps = root->appendChild("MaterialProperties");
  }

  ProblemSpecP mpm_ps = mat_ps->appendChild("MPM");
  for (int i = 0; i < d_materialManager->getNumMaterials("MPM"); i++) {
    MPMMaterial* mat   = d_materialManager->getMaterial("MPM", i);
    ProblemSpecP cm_ps = mat->outputProblemSpec(mpm_ps);
  }

  contactModel->outputProblemSpec(mpm_ps);
  thermalContactModel->outputProblemSpec(mpm_ps);

  for (int i = 0; i < d_materialManager->getNumMatls("CohesizeZone"); i++) {
    CZMaterial* mat    = d_materialManager->getCZMaterial(i);
    ProblemSpecP cm_ps = mat->outputProblemSpec(mpm_ps);
  }

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps   = physical_bc_ps->appendChild("MPM");
  for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    bc->outputProblemSpec(mpm_ph_bc_ps);
  }
}

/*!----------------------------------------------------------------------
 * scheduleInitialize
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }
  Task* t = scinew Task("MPM::actuallyInitialize",
                        this,
                        &SerialMPM::actuallyInitialize);

  const PatchSet* patches = level->eachPatch();
  printSchedule(patches, cout_doing, "MPM::scheduleInitialize");
  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->computes(d_mpmLabels->partCountLabel);
  t->computes(d_mpmLabels->pXLabel);
  t->computes(d_mpmLabels->pFiberDirLabel);
  t->computes(d_mpmLabels->pMassLabel);
  t->computes(d_mpmLabels->pVolumeLabel);
  t->computes(d_mpmLabels->pTemperatureLabel);
  t->computes(d_mpmLabels->pTempPreviousLabel); // for thermal stress analysis
  t->computes(d_mpmLabels->pdTdtLabel);
  t->computes(d_mpmLabels->pDispLabel);
  t->computes(d_mpmLabels->pVelocityLabel);
  t->computes(d_mpmLabels->pAccelerationLabel);
  t->computes(d_mpmLabels->pExternalForceLabel);
  t->computes(d_mpmLabels->pParticleIDLabel);
  t->computes(d_mpmLabels->pStressLabel);
  t->computes(d_mpmLabels->pSizeLabel);
  t->computes(d_mpmLabels->pRefinedLabel);
  t->computes(d_materialManager->get_delt_label(), level.get_rep());
  t->computes(d_mpmLabels->pCellNAPIDLabel, zeroth_matl);
  t->computes(d_mpmLabels->NC_CCweightLabel, zeroth_matl);

  if (!d_mpmFlags->d_doGridReset) {
    t->computes(d_mpmLabels->gDisplacementLabel);
  }

  // Debugging Scalar
  if (d_mpmFlags->d_withColor) {
    t->computes(d_mpmLabels->pColorLabel);
  }

  // Computes the load curve ID associated with each particle
  if (d_mpmFlags->d_useLoadCurves) {
    t->computes(d_mpmLabels->pLoadCurveIDLabel);
  }

  // Computes accumulated strain energy
  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    t->computes(d_mpmLabels->AccStrainEnergyLabel);
  }

  if (d_mpmFlags->d_artificialViscosity) {
    t->computes(d_mpmLabels->p_qLabel);
  }

  // artificial damping coeff initialized to 0.0
  if (cout_dbg.active()) {
    cout_dbg << "Artificial Damping Coeff = "
             << d_mpmFlags->d_artificialDampCoeff
             << " 8 or 27 = " << d_mpmFlags->d_8or27 << "\n";
  }

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    // For velocity gradinet and deformation gradient
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
      basicDamageModel->addInitialComputesAndRequires(t,
                                                      mpm_matl,
                                                      patches,
                                                      d_mpmLabels);
    }
  }

  // Add initialization of body force and coriolis importance terms
  // These are initialized to zero in ParticleCreator
  t->computes(d_mpmLabels->pCoriolisImportanceLabel);
  t->computes(d_mpmLabels->pBodyForceAccLabel);

  // Needed for switch from explicit to implicit MPM
  t->computes(d_mpmLabels->pExternalHeatFluxLabel);

  // For friction contact
  t->computes(d_mpmLabels->pSurfLabel);

  // Add task to scheduler
  sched->addTask(t, patches, d_materialManager->allMPMMaterials());

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference()) {
    delete zeroth_matl; // shouln't happen, but...
  }

  // Print particle count
  schedulePrintParticleCount(level, sched);

  // Compute initial stresses due to body forces and recompute the initial
  // deformation gradient
  if (d_mpmFlags->d_initializeStressFromBodyForce) {
    scheduleInitializeStressAndDefGradFromBodyForce(level, sched);
  }

  // Schedule the initialization of pressure BCs per particle
  if (d_mpmFlags->d_useLoadCurves) {
    if (MPMPhysicalBCFactory::mpmPhysicalBCs.size() > 0) {
      std::string bcType = MPMPhysicalBCFactory::mpmPhysicalBCs[0]->getType();
      if (bcType == "Pressure") {
        scheduleInitializePressureBCs(level, sched);
      } else if (bcType == "Moment") {
        scheduleInitializeMomentBCs(level, sched);
      }
    }
  }

  // dataAnalysis
  if (d_analysisModules.size() != 0) {
    std::vector<AnalysisModule*>::iterator iter;
    for (iter = d_analysisModules.begin(); iter != d_analysisModules.end();
         iter++) {
      AnalysisModule* am = *iter;
      am->scheduleInitialize(sched, level);
    }
  }

  // Cohesive zones
  int numCZM = d_materialManager->getNumMatls("CohesizeZone");
  for (int m = 0; m < numCZM; m++) {
    CZMaterial* cz_matl = d_materialManager->getMaterial("CohesiveZone", m);
    CohesiveZone* ch    = cz_matl->getCohesiveZone();
    ch->scheduleInitialize(level, sched, cz_matl);
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
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight,
                           d_mpmLabels->NC_CCweightLabel,
                           0,
                           patch);

    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
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

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int indx              = mpm_matl->getDWIndex();
      if (!d_mpmFlags->d_doGridReset) {
        NCVariable<Vector> gDisplacement;
        new_dw->allocateAndPut(gDisplacement,
                               d_mpmLabels->gDisplacementLabel,
                               indx,
                               patch);
        gDisplacement.initialize(Vector(0.));
      }

      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

      // Initialize deformation gradient
      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

      // Initialize constitutive models
      cm->initializeCMData(patch, mpm_matl, new_dw);

      // Initialize basic damage model
      if (mpm_matl->doBasicDamage()) {
        mpm_matl->getBasicDamageModel()->initializeDamageData(patch,
                                                              mpm_matl,
                                                              new_dw,
                                                              lb);
      }
    }
    IntVector num_extra_cells = patch->getExtraCells();
    IntVector periodic        = patch->getLevel()->getPeriodicBoundaries();
    std::string interp_type   = d_mpmFlags->d_interpolatorType;
    if (interp_type == "linear" && num_extra_cells != IntVector(0, 0, 0)) {
      if (!d_mpmFlags->d_withICE) {
        std::ostringstream msg;
        msg << "\n ERROR: When using <interpolator>linear</interpolator> \n"
            << " you should also use <extraCells>[0,0,0]</extraCells> \n"
            << " unless you are running an MPMICE or MPMARCHES case.\n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
    } else if (((interp_type == "gimp" || interp_type == "3rdorderBS" ||
                 interp_type == "cpdi") &&
                ((num_extra_cells + periodic) != IntVector(1, 1, 1) &&
                 ((num_extra_cells + periodic) != IntVector(1, 1, 0) &&
                  d_mpmFlags->d_axisymmetric)))) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>gimp</interpolator> \n"
          << " or <interpolator>3rdorderBS</interpolator> \n"
          << " or <interpolator>cpdi</interpolator> \n"
          << " you must also use extraCells and/or periodicBCs such\n"
          << " the sum of the two is [1,1,1].\n"
          << " If using axisymmetry, the sum of the two can be [1,1,0].\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }

    // Only allow axisymmetric runs if the grid is one cell thick in the theta
    // dir.
    if (d_mpmFlags->d_axisymmetric) {
      IntVector patchLowNode  = patch->getNodeLowIndex();
      IntVector patchHighNode = patch->getNodeHighIndex();
      int num_cells_in_theta  = (patchHighNode.z() - patchLowNode.z()) - 1;
      if (num_cells_in_theta > 1) {
        std::ostringstream msg;
        msg << "\n ERROR: When using <axisymmetric>true</axisymmetric> \n"
            << "the grid can only have one cell in the circumferential "
               "direction.\n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
    }
  }

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    // Initialize the accumulated strain energy
    new_dw->put(max_vartype(0.0), d_mpmLabels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_mpmLabels->partCountLabel);
}

/*!----------------------------------------------------------------------
 * schedulePrintParticleCount
 *-----------------------------------------------------------------------*/
void
SerialMPM::schedulePrintParticleCount(const LevelP& level, SchedulerP& sched)
{
  Task* t = scinew Task("MPM::printParticleCount",
                        this,
                        &SerialMPM::printParticleCount);
  t->requires(Task::NewDW, d_mpmLabels->partCountLabel);
  t->setType(Task::OncePerProc);
  sched->addTask(t,
                 sched->getLoadBalancer()->getPerProcessorPatchSet(level),
                 d_materialManager->allMPMMaterials());
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
  new_dw->get(pcount, d_mpmLabels->partCountLabel);

  if (pg->myRank() == 0) {
    std::cerr << "Created " << (long)pcount << " total particles\n";
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
  printSchedule(patches,
                cout_doing,
                "MPM::initializeStressAndDefGradFromBodyForce");

  // First compute the body force
  Task* t1 = scinew Task("MPM::initializeBodyForce",
                         this,
                         &SerialMPM::initializeBodyForce);
  t1->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::None);
  t1->modifies(d_mpmLabels->pBodyForceAccLabel);
  sched->addTask(t1, patches, d_materialManager->allMPMMaterials());

  // Compute the stress and deformation gradient only for selected
  // constitutive models that have a "initializeWithBodyForce" flag as true.
  // This is because a more general implementation is quite involved and
  // not worth the effort at this time. BB
  Task* t2 = scinew Task("MPM::initializeStressAndDefGradFromBodyForce",
                         this,
                         &SerialMPM::initializeStressAndDefGradFromBodyForce);

  t2->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::None);
  t2->requires(Task::NewDW, d_mpmLabels->pBodyForceAccLabel, Ghost::None);
  t2->modifies(d_mpmLabels->pStressLabel);
  t2->modifies(d_mpmLabels->pDefGradLabel);
  sched->addTask(t2, patches, d_materialManager->allMPMMaterials());
}

/*!------------------------------------------------------------------------
 * Actually initialize the body force acceleration
 *-------------------------------------------------------------------------*/
void
SerialMPM::initializeBodyForce(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset* matls,
                               DataWarehouse*,
                               DataWarehouse* new_dw)
{
  // Get the MPM d_mpmFlags and make local copies
  Uintah::Point rotation_center = d_mpmFlags->d_coordRotationCenter;
  Uintah::Vector rotation_axis  = d_mpmFlags->d_coordRotationAxis;
  double rotation_speed         = d_mpmFlags->d_coordRotationSpeed;

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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      new_dw->getModifiable(pBodyForceAcc,
                            d_mpmLabels->pBodyForceAccLabel,
                            pset);

      // Get the position data
      constParticleVariable<Point> pPosition;
      new_dw->get(pPosition, d_mpmLabels->pXLabel, pset);

      // Iterate over the particles
      for (auto pidx : *pset) {

        // Compute the body force acceleration (g)
        // Just use gravity if rotation is off
        pBodyForceAcc[pidx] = d_mpmFlags->d_gravity;

        // If rotating add centrifugal force
        if (d_mpmFlags->d_useCoordRotation) {

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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

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
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpmLabels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_materialManager->allMPMMaterials());

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("MPM::initializePressureBC",
                    this,
                    &SerialMPM::initializePressureBC);
    t->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pSizeLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pDispLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpmLabels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpmLabels->pExternalForceLabel);
    if (d_mpmFlags->d_useCBDI) {
      t->computes(d_mpmLabels->pExternalForceCorner1Label);
      t->computes(d_mpmLabels->pExternalForceCorner2Label);
      t->computes(d_mpmLabels->pExternalForceCorner3Label);
      t->computes(d_mpmLabels->pExternalForceCorner4Label);
    }
    sched->addTask(t, patches, d_materialManager->allMPMMaterials());
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
  printTask(patches,
            patches->get(0),
            cout_doing,
            "countMaterialPointsPerLoadCurve");
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
          MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
          int matID             = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
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
                  d_mpmLabels->materialPointsPerLoadCurveLabel,
                  0,
                  pressureBCId++);

      // Save the material points per load curve in the PressureBC object
      PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
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
          MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
          int matID             = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

          constParticleVariable<Point> pX;
          constParticleVariable<Matrix3> pSize;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX, d_mpmLabels->pXLabel, pset);
          new_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
          new_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);

          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);
          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(pExternalForce,
                                d_mpmLabels->pExternalForceLabel,
                                pset);

          constParticleVariable<Vector> pDisp;
          new_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);

          ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
            pExternalForceCorner3, pExternalForceCorner4;
          if (d_mpmFlags->d_useCBDI) {
            if (ii == 0) {
              new_dw->allocateAndPut(pExternalForceCorner1,
                                     d_mpmLabels->pExternalForceCorner1Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner2,
                                     d_mpmLabels->pExternalForceCorner2Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner3,
                                     d_mpmLabels->pExternalForceCorner3Label,
                                     pset);
              new_dw->allocateAndPut(pExternalForceCorner4,
                                     d_mpmLabels->pExternalForceCorner4Label,
                                     pset);
            } else {
              new_dw->getModifiable(pExternalForceCorner1,
                                    d_mpmLabels->pExternalForceCorner1Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner2,
                                    d_mpmLabels->pExternalForceCorner2Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner3,
                                    d_mpmLabels->pExternalForceCorner3Label,
                                    pset);
              new_dw->getModifiable(pExternalForceCorner4,
                                    d_mpmLabels->pExternalForceCorner4Label,
                                    pset);
            }
          }

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == pressureBCId) {
              if (d_mpmFlags->d_useCBDI) {
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
                pExternalForce[idx] = pbc->getForceVector(pX[idx],
                                                          pDisp[idx],
                                                          forcePerPart,
                                                          time,
                                                          pDefGrad[idx]);
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
 * scheduleInitializeMoemntBCs
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
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpmLabels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_materialManager->allMPMMaterials());

    // Create a task that calculates the force to be associated with
    // each particle based on the moment BCs
    t = scinew Task("MPM::initializeMomentBC",
                    this,
                    &SerialMPM::initializeMomentBC);
    t->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pSizeLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpmLabels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex,
                Task::OutOfDomain,
                Ghost::None);
    t->modifies(d_mpmLabels->pExternalForceLabel);
    if (d_mpmFlags->d_useCBDI) {
      t->computes(d_mpmLabels->pExternalForceCorner1Label);
      t->computes(d_mpmLabels->pExternalForceCorner2Label);
      t->computes(d_mpmLabels->pExternalForceCorner3Label);
      t->computes(d_mpmLabels->pExternalForceCorner4Label);
    }
    sched->addTask(t, patches, d_materialManager->allMPMMaterials());
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
                  d_mpmLabels->materialPointsPerLoadCurveLabel,
                  0,
                  nofMomentBCs++);

      // Save the material points per load curve in the MomentBC object
      MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
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
          MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
          int matID             = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<Point> pX;
          constParticleVariable<Matrix3> pSize;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX, d_mpmLabels->pXLabel, pset);
          new_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
          new_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);
          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(pExternalForce,
                                d_mpmLabels->pExternalForceLabel,
                                pset);

          ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
            pExternalForceCorner3, pExternalForceCorner4;
          if (d_mpmFlags->d_useCBDI) {
            new_dw->allocateAndPut(pExternalForceCorner1,
                                   d_mpmLabels->pExternalForceCorner1Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner2,
                                   d_mpmLabels->pExternalForceCorner2Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner3,
                                   d_mpmLabels->pExternalForceCorner3Label,
                                   pset);
            new_dw->allocateAndPut(pExternalForceCorner4,
                                   d_mpmLabels->pExternalForceCorner4Label,
                                   pset);
          }
          // std::cout << "d_mpmFlags->d_useCBDI: " << d_mpmFlags->d_useCBDI <<
          // "\n";
          ParticleSubset::iterator iter = pset->begin();
          for (; iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            if (pLoadCurveID[idx] == nofMomentBCs) {
              if (d_mpmFlags->d_useCBDI) {
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
                pExternalForce[idx] = pbc->getForceVector(pX[idx],
                                                          forcePerPart,
                                                          time,
                                                          pDefGrad[idx]);
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
  Task* t = 0;
  cout_doing << UintahParallelComponent::d_myworld->myRank()
             << " MPM::scheduleComputeStableTimestep \t\t\t\tL-"
             << level->getIndex() << "\n";

  t = scinew Task("MPM::actuallyComputeStableTimestep",
                  this,
                  &SerialMPM::actuallyComputeStableTimestep);

  const MaterialSet* mpm_matls = d_materialManager->allMPMMaterials();

  t->computes(d_materialManager->get_delt_label(), level.get_rep());
  sched->addTask(t, level->eachPatch(), mpm_matls);
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
  // Put something here to satisfy the need for a reduction operation in
  // the case that there are multiple levels present
  const Level* level = getLevel(patches);
  new_dw->put(delt_vartype(999.0), d_mpmLabels->delTLabel, level);
}

/*!----------------------------------------------------------------------
 * scheduleTimeAdvance
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  MALLOC_TRACE_TAG_SCOPE("SerialMPM::scheduleTimeAdvance()");
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  const PatchSet* patches      = level->eachPatch();
  const MaterialSet* matls     = d_materialManager->allMPMMaterials();
  const MaterialSet* cz_matls  = d_materialManager->allCZMaterials();
  const MaterialSet* all_matls = d_materialManager->allMaterials();

  const MaterialSubset* mpm_matls_sub = matls->getUnion();
  const MaterialSubset* cz_matls_sub  = cz_matls->getUnion();

  // Compute body forces first
  scheduleComputeParticleBodyForce(sched, patches, matls);

  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleInterpolateParticlesToGrid(sched, patches, matls);

  scheduleComputeNormals(sched, patches, matls);
  scheduleFindSurfaceParticles(sched, patches, matls);
  scheduleComputeLogisticRegression(sched, patches, matls);

  scheduleExMomInterpolated(sched, patches, matls);
  if (d_mpmFlags->d_useCohesiveZones) {

    scheduleUpdateCohesiveZones(sched,
                                patches,
                                mpm_matls_sub,
                                cz_matls_sub,
                                all_matls);

    scheduleAddCohesiveZoneForces(sched,
                                  patches,
                                  mpm_matls_sub,
                                  cz_matls_sub,
                                  all_matls);
  }
  scheduleComputeContactArea(sched, patches, matls);
  scheduleComputeInternalForce(sched, patches, matls);

  scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
  scheduleExMomIntegrated(sched, patches, matls);
  scheduleSetGridBoundaryConditions(sched, patches, matls);
  scheduleSetPrescribedMotion(sched, patches, matls);

// For XPIC(2) computations
#ifdef XPIC2_UPDATE
  scheduleComputeXPICVelocities(sched, patches, matls);
#endif

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
  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }

  if (d_mpmFlags->d_doExplicitHeatConduction) {
    scheduleComputeHeatExchange(sched, patches, matls);
    scheduleComputeInternalHeatRate(sched, patches, matls);
    scheduleComputeNodalHeatFlux(sched, patches, matls);
    scheduleSolveHeatEquations(sched, patches, matls);
    scheduleIntegrateTemperatureRate(sched, patches, matls);
  }

  scheduleAddNewParticles(sched, patches, matls);
  scheduleConvertLocalizedParticles(sched, patches, matls);
  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);

  scheduleInsertParticles(sched, patches, matls);
  if (d_mpmFlags->d_refineParticles) {
    scheduleAddParticles(sched, patches, matls);
  }
  if (d_mpmFlags->d_computeScaleFactor) {
    scheduleComputeParticleScaleFactor(sched, patches, matls);
  }

  if (d_analysisModules.size() != 0) {
    std::vector<AnalysisModule*>::iterator iter;
    for (iter = d_analysisModules.begin(); iter != d_analysisModules.end();
         iter++) {
      AnalysisModule* am = *iter;
      am->scheduleDoAnalysis_preReloc(sched, level);
    }
  }

  sched->scheduleParticleRelocation(level,
                                    d_mpmLabels->pXLabel_preReloc,
                                    d_materialManager->d_particleState_preReloc,
                                    d_mpmLabels->pXLabel,
                                    d_materialManager->d_particleState,
                                    d_mpmLabels->pParticleIDLabel,
                                    matls,
                                    1);

  if (d_mpmFlags->d_useCohesiveZones) {
    sched->scheduleParticleRelocation(
      level,
      d_mpmLabels->pXLabel_preReloc,
      d_materialManager->d_cohesiveZoneState_preReloc,
      d_mpmLabels->pXLabel,
      d_materialManager->d_cohesiveZoneState,
      d_mpmLabels->czIDLabel,
      cz_matls,
      2);
  }

  //__________________________________
  //  on the fly analysis
  if (d_analysisModules.size() != 0) {
    std::vector<AnalysisModule*>::iterator iter;
    for (iter = d_analysisModules.begin(); iter != d_analysisModules.end();
         iter++) {
      AnalysisModule* am = *iter;
      am->scheduleDoAnalysis(sched, level);
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleBodyForce");

  Task* t = scinew Task("MPM::computeParticleBodyForce",
                        this,
                        &SerialMPM::computeParticleBodyForce);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pVelocityLabel, Ghost::None);
  // t->computes(d_mpmLabels->pBodyForceAccLabel);
  // t->computes(d_mpmLabels->pCoriolisImportanceLabel);
  t->computes(d_mpmLabels->pBodyForceAccLabel_preReloc);
  t->computes(d_mpmLabels->pCoriolisImportanceLabel_preReloc);

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
  // Get the MPM d_mpmFlags and make local copies
  Uintah::Point rotation_center = d_mpmFlags->d_coordRotationCenter;
  Uintah::Vector rotation_axis  = d_mpmFlags->d_coordRotationAxis;
  double rotation_speed         = d_mpmFlags->d_coordRotationSpeed;
  // Uintah::Point body_ref_point = d_mpmFlags->d_coord_rotation_body_ref_point;

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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      // new_dw->allocateAndPut(pBodyForceAcc, d_mpmLabels->pBodyForceAccLabel,
      // pset);
      new_dw->allocateAndPut(pBodyForceAcc,
                             d_mpmLabels->pBodyForceAccLabel_preReloc,
                             pset);

      // Create space for particle coriolis importance
      ParticleVariable<double> pCoriolisImportance;
      // new_dw->allocateAndPut(pCoriolisImportance,
      // d_mpmLabels->pCoriolisImportanceLabel, pset);
      new_dw->allocateAndPut(pCoriolisImportance,
                             d_mpmLabels->pCoriolisImportanceLabel_preReloc,
                             pset);

      // Don't do much if coord rotation is off
      if (!d_mpmFlags->d_useCoordRotation) {

        // Iterate over the particles
        for (auto pidx : *pset) {

          // Compute the body force acceleration (g)
          pBodyForceAcc[pidx] = d_mpmFlags->d_gravity;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[pidx] = 0.0;
        } // particle loop

      } else { // Use coordinate rotation

        // Get the particle data
        constParticleVariable<Point> pPosition;
        old_dw->get(pPosition, d_mpmLabels->pXLabel, pset);

        constParticleVariable<Vector> pVelocity;
        old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);

        // Iterate over the particles
        // std::cout << "Mat id = " << matID << " patch = " << patch << "\n";
        // std::cout << "Particle subset = " << *pset;
        // std::cout << "Num particles = " << pset->numParticles() << "\n";
        for (auto iter = pset->begin(); iter != pset->end(); iter++) {
          particleIndex pidx = *iter;

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
            d_mpmFlags->d_gravity - centrifugal_accel - coriolis_accel;

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
        d_mpmLabels->pCoriolisImportanceLabel_preReloc, pset);
        pCoriolisImportance_new.copyData(pCoriolisImportance);

        ParticleVariable<Vector> pBodyForceAcc_new;
        new_dw->allocateAndPut(pBodyForceAcc_new,
        d_mpmLabels->pBodyForceAccLabel_preReloc, pset);
        pBodyForceAcc_new.copyData(pBodyForceAcc);
      */

    } // matl loop
  }   // patch loop
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleApplyExternalLoads");

  Task* t = scinew Task("MPM::applyExternalLoads",
                        this,
                        &SerialMPM::applyExternalLoads);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pDispLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pExternalForceLabel, Ghost::None);
  t->computes(d_mpmLabels->pExtForceLabel_preReloc);
  if (d_mpmFlags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_mpmLabels->pLoadCurveIDLabel_preReloc);
    if (d_mpmFlags->d_useCBDI) {
      t->computes(d_mpmLabels->pExternalForceCorner1Label);
      t->computes(d_mpmLabels->pExternalForceCorner2Label);
      t->computes(d_mpmLabels->pExternalForceCorner3Label);
      t->computes(d_mpmLabels->pExternalForceCorner4Label);
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
  double time = d_materialManager->getElapsedTime();

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << "\n";
  }

  // Calculate the force std::vector at each particle for each pressure bc
  std::vector<double> forcePerPart;
  std::vector<PressureBC*> pbcP;
  std::vector<MomentBC*> pbcM;
  if (d_mpmFlags->d_useLoadCurves) {
    for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
      std::string bcType = bc->getType();
      if (bcType == "Pressure") {
        PressureBC* pbc = dynamic_cast<PressureBC*>(bc.get());
        pbcP.push_back(pbc);

        // Calculate the force per particle at current time
        forcePerPart.push_back(pbc->forcePerParticle(time));
      } else if (bcType == "Moment") {
        MomentBC* pbc = dynamic_cast<MomentBC*>(bc.get());
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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      // Get the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);

      ParticleVariable<Vector> pExternalForce_new;
      ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
        pExternalForceCorner3, pExternalForceCorner4;

      new_dw->allocateAndPut(pExternalForce_new,
                             d_mpmLabels->pExtForceLabel_preReloc,
                             pset);
      if (d_mpmFlags->d_useCBDI) {
        new_dw->allocateAndPut(pExternalForceCorner1,
                               d_mpmLabels->pExternalForceCorner1Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner2,
                               d_mpmLabels->pExternalForceCorner2Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner3,
                               d_mpmLabels->pExternalForceCorner3Label,
                               pset);
        new_dw->allocateAndPut(pExternalForceCorner4,
                               d_mpmLabels->pExternalForceCorner4Label,
                               pset);
      }
      // std::cout << "applyloads: patch = " << patch << " matID = " << matID
      //           << " numparticles = " << pset->numParticles() << "
      //           d_numGhostParticles = " << d_numGhostParticles << "\n";

      if (d_mpmFlags->d_useLoadCurves) {

        // Get the load curve data
        // Recycle the loadCurveIDs
        constParticleVariable<int> pLoadCurveID;
        old_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);

        ParticleVariable<int> pLoadCurveID_new;
        new_dw->allocateAndPut(pLoadCurveID_new,
                               d_mpmLabels->pLoadCurveIDLabel_preReloc,
                               pset);
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
          old_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);

          // Iterate over the particles
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = Vector(0.0, 0.0, 0.0);
              if (d_mpmFlags->d_useCBDI) {
                pExternalForceCorner1[idx] = pX[idx];
                pExternalForceCorner2[idx] = pX[idx];
                pExternalForceCorner3[idx] = pX[idx];
                pExternalForceCorner4[idx] = pX[idx];
              }
            } else {
              PressureBC* pbc = pbcP[loadCurveID];
              double force    = forcePerPart[loadCurveID];

              if (d_mpmFlags->d_useCBDI) {
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
                pExternalForce_new[idx] = pbc->getForceVector(pX[idx],
                                                              pDisp[idx],
                                                              force,
                                                              time,
                                                              pDefGrad[idx]);
              }
            }
          }
        } else if (do_MomentBCs) {
          // Get the external force data and allocate new space for
          // external force
          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_mpmLabels->pExternalForceLabel, pset);

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
        std::string mms_type = d_mpmFlags->d_mmsType;
        if (!mms_type.empty()) {
          MMS MMSObject;
          MMSObject.computeBodyForceForMMS(old_dw,
                                           new_dw,
                                           time,
                                           pset,
                                           d_mpmLabels,
                                           d_mpmFlags,
                                           pExternalForce_new);
        }

      } else { // d_useLoadCurves = False
        // MMS
        std::string mms_type = d_mpmFlags->d_mmsType;
        if (!mms_type.empty()) {
          MMS MMSObject;
          MMSObject.computeExternalForceForMMS(old_dw,
                                               new_dw,
                                               time,
                                               pset,
                                               d_mpmLabels,
                                               d_mpmFlags,
                                               pExternalForce_new);
        } else {
          // Get the external force data and allocate new space for
          // external force and copy the data
          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_mpmLabels->pExternalForceLabel, pset);

          for (auto idx : *pset) {
            pExternalForce_new[idx] =
              pExternalForce[idx] * d_mpmFlags->d_forceIncrementFactor;
          }
        }
      } // end if (d_useLoadCurves)
    }   // matl loop
  }     // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateParticlesToGrid
 * interpolateParticlesToGrid
 *   in(P.MASS, P.VELOCITY, P.NAT_X)
 *   operation(interpolate the P.MASS and P.VEL to the grid
 *             using P.NAT_X and some shape function evaluations)
 *   out(G.MASS, G.VELOCITY)
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleInterpolateParticlesToGrid");

  Task* t              = scinew Task("MPM::interpolateParticlesToGrid",
                        this,
                        &SerialMPM::interpolateParticlesToGrid);
  Ghost::GhostType gan = Ghost::AroundNodes;
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pVolumeLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pVelocityLabel,
              gan,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gan, d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pBodyForceAccLabel_preReloc,
              gan,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pExtForceLabel_preReloc,
              gan,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pTemperatureLabel,
              gan,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              gan,
              d_numGhostParticles);
  if (d_mpmFlags->d_useLoadCurves) {
    t->requires(Task::OldDW,
                d_mpmLabels->pLoadCurveIDLabel,
                gan,
                d_numGhostParticles);
    if (d_mpmFlags->d_useCBDI) {
      t->requires(Task::NewDW,
                  d_mpmLabels->pExternalForceCorner1Label,
                  gan,
                  d_numGhostParticles);
      t->requires(Task::NewDW,
                  d_mpmLabels->pExternalForceCorner2Label,
                  gan,
                  d_numGhostParticles);
      t->requires(Task::NewDW,
                  d_mpmLabels->pExternalForceCorner3Label,
                  gan,
                  d_numGhostParticles);
      t->requires(Task::NewDW,
                  d_mpmLabels->pExternalForceCorner4Label,
                  gan,
                  d_numGhostParticles);
    }
  }

#ifdef DEBUG_WITH_PARTICLE_ID
  t->requires(Task::OldDW,
              d_mpmLabels->pParticleIDLabel,
              gan,
              d_numGhostParticles);
#endif

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
  t->computes(d_mpmLabels->gBodyForceLabel);
  t->computes(d_mpmLabels->gExternalForceLabel);
  t->computes(d_mpmLabels->gTemperatureLabel);
  t->computes(d_mpmLabels->gTemperatureNoBCLabel);
  t->computes(d_mpmLabels->gTemperatureRateLabel);
  t->computes(d_mpmLabels->gExternalHeatRateLabel);

  if (d_mpmFlags->d_withICE) {
    t->computes(d_mpmLabels->gVelocityBCLabel);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * interpolateParticlesToGrid
 *-----------------------------------------------------------------------*/
void
SerialMPM::interpolateParticlesToGrid(const ProcessorGroup*,
                                      const PatchSubset* patches,
                                      const MaterialSubset*,
                                      DataWarehouse* old_dw,
                                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing interpolateParticlesToGrid");

    int numMatls             = d_materialManager->getNumMaterials("MPM");
    auto interpolator        = d_mpmFlags->d_interpolator->clone(patch);
    auto linear_interpolator = std::make_unique<LinearInterpolator>(patch);

    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);
    std::string interp_type = d_mpmFlags->d_interpolatorType;

    NCVariable<double> gMassglobal, gTempglobal, gVolumeglobal;
    NCVariable<Vector> gVelglobal;
    new_dw->allocateAndPut(gMassglobal,
                           d_mpmLabels->gMassLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gTempglobal,
                           d_mpmLabels->gTemperatureLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gVolumeglobal,
                           d_mpmLabels->gVolumeLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    new_dw->allocateAndPut(gVelglobal,
                           d_mpmLabels->gVelocityLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);
    gMassglobal.initialize(d_SMALL_NUM_MPM);
    gVolumeglobal.initialize(d_SMALL_NUM_MPM);
    gTempglobal.initialize(0.0);
    gVelglobal.initialize(Vector(0.0));
    Ghost::GhostType gan = Ghost::AroundNodes;
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume, pTemperature;
      constParticleVariable<Vector> pVelocity, pBodyForceAcc, pExternalForce;
      constParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
        pExternalForceCorner3, pExternalForceCorner4;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad_old;

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pVolume, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);
      new_dw->get(pBodyForceAcc,
                  d_mpmLabels->pBodyForceAccLabel_preReloc,
                  pset);
      new_dw->get(pExternalForce, d_mpmLabels->pExtForceLabel_preReloc, pset);

      /*
      std::cout << "patch = " << patch << " matID = " << matID
                << " numparticles = " << pset->numParticles() << "
      d_numGhostParticles = " << d_numGhostParticles << "\n";
      */

      constParticleVariable<int> pLoadCurveID;
      if (d_mpmFlags->d_useLoadCurves) {
        old_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);
        if (d_mpmFlags->d_useCBDI) {
          new_dw->get(pExternalForceCorner1,
                      d_mpmLabels->pExternalForceCorner1Label,
                      pset);

          /*
          for (auto idx : *pset) {
            std::cout << "idx = " << idx << " px = " << pX[idx] << " fext = " <<
          pExternalForce[idx]
                      << " curve id = " << pLoadCurveID[idx] << "\n";
            std::cout << "corner1 = " << pExternalForceCorner1[idx] << "\n";
          }
          */

          new_dw->get(pExternalForceCorner2,
                      d_mpmLabels->pExternalForceCorner2Label,
                      pset);
          new_dw->get(pExternalForceCorner3,
                      d_mpmLabels->pExternalForceCorner3Label,
                      pset);
          new_dw->get(pExternalForceCorner4,
                      d_mpmLabels->pExternalForceCorner4Label,
                      pset);
        }
      }

#ifdef DEBUG_WITH_PARTICLE_ID
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
#endif

      // Create arrays for the grid data
      NCVariable<double> gMass;
      NCVariable<double> gVolume;
      NCVariable<Vector> gVelocity;
      NCVariable<Vector> gBodyForce;
      NCVariable<Vector> gExternalForce;
      NCVariable<double> gExternalheatrate;
      NCVariable<double> gTemperature;
      NCVariable<double> gSp_vol;
      NCVariable<double> gTemperatureNoBC;
      NCVariable<double> gTemperatureRate;
      // NCVariable<double> gnumnearparticles;

      new_dw->allocateAndPut(gMass, d_mpmLabels->gMassLabel, matID, patch);
      new_dw->allocateAndPut(gSp_vol, d_mpmLabels->gSp_volLabel, matID, patch);
      new_dw->allocateAndPut(gVolume, d_mpmLabels->gVolumeLabel, matID, patch);
      new_dw->allocateAndPut(gVelocity,
                             d_mpmLabels->gVelocityLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gTemperature,
                             d_mpmLabels->gTemperatureLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gTemperatureNoBC,
                             d_mpmLabels->gTemperatureNoBCLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gTemperatureRate,
                             d_mpmLabels->gTemperatureRateLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gBodyForce,
                             d_mpmLabels->gBodyForceLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gExternalForce,
                             d_mpmLabels->gExternalForceLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gExternalheatrate,
                             d_mpmLabels->gExternalHeatRateLabel,
                             matID,
                             patch);

      gMass.initialize(d_SMALL_NUM_MPM);
      gVolume.initialize(d_SMALL_NUM_MPM);
      gVelocity.initialize(Vector(0, 0, 0));
      gBodyForce.initialize(Vector(0, 0, 0));
      gExternalForce.initialize(Vector(0, 0, 0));
      gTemperature.initialize(0);
      gTemperatureNoBC.initialize(0);
      gTemperatureRate.initialize(0);
      gExternalheatrate.initialize(0);
      gSp_vol.initialize(0.);
      // gnumnearparticles.initialize(0.);

      // Interpolate particle data to Grid data.
      // This currently consists of the particle velocity and mass
      // Need to compute the lumped global mass matrix and velocity
      // Vector from the individual mass matrix and velocity std::vector
      // GridMass * GridVelocity =  S^T*M_D*ParticleVelocity

      Vector total_mom(0.0, 0.0, 0.0);
      Vector pMom;
      double pSp_vol = 1. / mpm_matl->getInitialDensity();

      // loop over all particles in the patch:
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {
        particleIndex idx = *iter;
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad_old[idx]);
        pMom = pVelocity[idx] * pMass[idx];
        total_mom += pMom;

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        IntVector node;
        for (int k = 0; k < numInfluenceNodes;
             k++) { // Iterates through the nodes which
                    // receive information from the current particle
          node = ni[k];
          if (patch->containsNode(node)) {
            gMass[node] += pMass[idx] * S[k];
            gVelocity[node] += pMom * S[k];
            gVolume[node] += pVolume[idx] * S[k];
            if (!d_mpmFlags->d_useCBDI) {
              gExternalForce[node] += pExternalForce[idx] * S[k];
            }
            gBodyForce[node] += pBodyForceAcc[idx] * pMass[idx] * S[k];
            gTemperature[node] += pTemperature[idx] * pMass[idx] * S[k];
            gSp_vol[node] += pSp_vol * pMass[idx] * S[k];
#ifdef DEBUG_WITH_PARTICLE_ID
            if (pParticleID[idx] == testParticleID) {
              proc0cout << pParticleID[idx] << pExternalForce[idx]
                        << " pMom = " << pMom
                        << " pVelocity = " << pVelocity[idx]
                        << " node = " << node
                        << " gVelocity = " << gVelocity[node] << "\n";
            }
#endif
          }
        }
        if (d_mpmFlags->d_useLoadCurves && d_mpmFlags->d_useCBDI) {
          std::vector<IntVector> niCorner1(linear_interpolator->size());
          std::vector<IntVector> niCorner2(linear_interpolator->size());
          std::vector<IntVector> niCorner3(linear_interpolator->size());
          std::vector<IntVector> niCorner4(linear_interpolator->size());
          std::vector<double> SCorner1(linear_interpolator->size());
          std::vector<double> SCorner2(linear_interpolator->size());
          std::vector<double> SCorner3(linear_interpolator->size());
          std::vector<double> SCorner4(linear_interpolator->size());
          linear_interpolator->findCellAndWeights(pExternalForceCorner1[idx],
                                                  niCorner1,
                                                  SCorner1,
                                                  pSize[idx],
                                                  pDefGrad_old[idx]);
          linear_interpolator->findCellAndWeights(pExternalForceCorner2[idx],
                                                  niCorner2,
                                                  SCorner2,
                                                  pSize[idx],
                                                  pDefGrad_old[idx]);
          linear_interpolator->findCellAndWeights(pExternalForceCorner3[idx],
                                                  niCorner3,
                                                  SCorner3,
                                                  pSize[idx],
                                                  pDefGrad_old[idx]);
          linear_interpolator->findCellAndWeights(pExternalForceCorner4[idx],
                                                  niCorner4,
                                                  SCorner4,
                                                  pSize[idx],
                                                  pDefGrad_old[idx]);
          for (int k = 0; k < 8;
               k++) { // Iterates through the nodes which receive information
                      // from the current particle
            node = niCorner1[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] += pExternalForce[idx] * SCorner1[k];
            }
            node = niCorner2[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] += pExternalForce[idx] * SCorner2[k];
            }
            node = niCorner3[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] += pExternalForce[idx] * SCorner3[k];
            }
            node = niCorner4[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] += pExternalForce[idx] * SCorner4[k];
            }
          }
        }
      } // End of particle loop
      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        gMassglobal[c] += gMass[c];
        gVolumeglobal[c] += gVolume[c];
        gVelglobal[c] += gVelocity[c];
        gVelocity[c] /= gMass[c];
        gTempglobal[c] += gTemperature[c];
        // gBodyForce[c]     /= gMass[c];
        gTemperature[c] /= gMass[c];
        gTemperatureNoBC[c] = gTemperature[c];
        gSp_vol[c] /= gMass[c];
      }

      // Apply boundary conditions to the temperature and velocity (if symmetry)
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch, matID, "Velocity", gVelocity, interp_type);
      bc.setBoundaryCondition(patch,
                              matID,
                              "Temperature",
                              gTemperature,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              matID,
                              "Symmetric",
                              gVelocity,
                              interp_type);

      // If an MPMICE problem, create a velocity with BCs variable for NCToCC_0
      if (d_mpmFlags->d_withICE) {
        NCVariable<Vector> gVelocityWBC;
        new_dw->allocateAndPut(gVelocityWBC,
                               d_mpmLabels->gVelocityBCLabel,
                               matID,
                               patch);
        gVelocityWBC.copyData(gVelocity);
        bc.setBoundaryCondition(patch,
                                matID,
                                "Velocity",
                                gVelocityWBC,
                                interp_type);
        bc.setBoundaryCondition(patch,
                                matID,
                                "Symmetric",
                                gVelocityWBC,
                                interp_type);
      }
    } // End loop over materials

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gTempglobal[c] /= gMassglobal[c];
      gVelglobal[c] /= gMassglobal[c];
    }
    // delete interpolator;
    // delete linear_interpolator;
  } // End loop over patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeNormals: For contact
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeNormals(SchedulerP& sched,
                                  const PatchSet* patches,
                                  const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "MPM::scheduleComputeNormals");

  Task* t =
    scinew Task("MPM::computeNormals", this, &SerialMPM::computeNormals);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();

  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pMassLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pVolumeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pSizeLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pStressLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::AroundNodes, 1);
  // t->requires(Task::OldDW, d_mpmLabels->NC_CCweightLabel, z_matl,
  // Ghost::None);

  t->computes(d_mpmLabels->gSurfNormLabel);
  t->computes(d_mpmLabels->gStressLabel);
  t->computes(d_mpmLabels->gNormTractionLabel);
  t->computes(d_mpmLabels->gPositionLabel);

  sched->addTask(t, patches, matls);

  if (z_matl->removeReference()) {
    delete z_matl; // shouldn't happen, but...
  }
}

//______________________________________________________________________
//
void
SerialMPM::computeNormals(const ProcessorGroup*,
                          const PatchSubset* patches,
                          const MaterialSubset*,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  Ghost::GhostType gan = Ghost::AroundNodes;
  // Ghost::GhostType  gnone = Ghost::None;

  auto numMPMMatls = d_materialManager->getNumMaterials("MPM");
  std::vector<constNCVariable<double>> gMass(numMPMMatls);
  std::vector<NCVariable<Point>> gPosition(numMPMMatls);
  std::vector<NCVariable<Vector>> gvelocity(numMPMMatls);
  std::vector<NCVariable<Vector>> gSurfNorm(numMPMMatls);
  std::vector<NCVariable<double>> gNormTraction(numMPMMatls);
  std::vector<NCVariable<Matrix3>> gStress(numMPMMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing MPM::computeNormals");

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();

    // constNCVariable<double>    NC_CCweight;
    // old_dw->get(NC_CCweight,   d_mpmLabels->NC_CCweightLabel,  0, patch,
    // gnone, 0);

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);
    std::vector<Vector> d_S(numInfluenceNodes);
    std::string interp_type = d_mpmFlags->d_interpolatorType;

    // Find surface normal at each material based on a gradient of nodal mass
    for (auto m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      new_dw->get(gMass[m], d_mpmLabels->gMassLabel, matID, patch, gan, 1);

      new_dw->allocateAndPut(gSurfNorm[m],
                             d_mpmLabels->gSurfNormLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gPosition[m],
                             d_mpmLabels->gPositionLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gStress[m],
                             d_mpmLabels->gStressLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gNormTraction[m],
                             d_mpmLabels->gNormTractionLabel,
                             matID,
                             patch);

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume;
      constParticleVariable<Matrix3> pSize, pStress;
      constParticleVariable<Matrix3> pDefGrad_old;

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pVolume, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pStress, d_mpmLabels->pStressLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);

      gSurfNorm[m].initialize(Vector(0.0, 0.0, 0.0));
      gPosition[m].initialize(Point(0.0, 0.0, 0.0));
      gNormTraction[m].initialize(0.0);
      gStress[m].initialize(Matrix3(0.0));

      if (d_mpmFlags->d_axisymmetric) {
        for (auto idx : *pset) {
          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[idx],
            ni,
            S,
            d_S,
            pSize[idx],
            pDefGrad_old[idx]);
          double rho = pMass[idx] / pVolume[idx];
          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              Vector G(d_S[k].x(), d_S[k].y(), 0.0);
              gSurfNorm[m][node] += rho * G;
              gPosition[m][node] += pX[idx].asVector() * pMass[idx] * S[k];
              gStress[m][node] += pStress[idx] * S[k];
            }
          }
        }
      } else {
        for (auto idx : *pset) {
          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[idx],
            ni,
            S,
            d_S,
            pSize[idx],
            pDefGrad_old[idx]);
          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              Vector grad(d_S[k].x() * oodx[0],
                          d_S[k].y() * oodx[1],
                          d_S[k].z() * oodx[2]);
              gSurfNorm[m][node] += pMass[idx] * grad;
              gPosition[m][node] += pX[idx].asVector() * pMass[idx] * S[k];
              gStress[m][node] += pStress[idx] * S[k];
            }
          }
        }
      } // axisymmetric conditional
    }   // matl loop

    // Make normals collinear by taking an average with the
    // other materials at a node
    if (d_mpmFlags->d_computeCollinearNormals) {
      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        std::vector<Vector> norm_temp(numMPMMatls);
        for (auto m = 0; m < numMPMMatls; m++) {
          norm_temp[m] = Vector(0., 0., 0);
          if (gMass[m][node] > 1.e-200) {
            Vector mWON(0., 0., 0.);
            double mON = 0.0;
            for (auto n = 0; n < numMPMMatls; n++) {
              if (n != m) {
                mWON += gMass[n][node] * gSurfNorm[n][node];
                mON += gMass[n][node];
              }
            } // loop over other matls
            mWON /= (mON + 1.e-100);
            norm_temp[m] = 0.5 * (gSurfNorm[m][node] - mWON);
          } // If node has mass
        }   // Outer loop over materials

        // Now put temporary norm into main array
        for (auto m = 0; m < numMPMMatls; m++) {
          gSurfNorm[m][node] = norm_temp[m];
        } // Outer loop over materials
      }   // Loop over nodes
    }     // if(flags..)

    // Make traditional norms unit length, compute gNormTraction
    for (auto m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              matID,
                              "Symmetric",
                              gSurfNorm[m],
                              interp_type);

      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        double length  = gSurfNorm[m][node].length();
        if (length > 1.0e-15) {
          gSurfNorm[m][node] = gSurfNorm[m][node] / length;
        }
        Vector norm            = gSurfNorm[m][node];
        gNormTraction[m][node] = Dot((norm * gStress[m][node]), norm);
        gPosition[m][node] /= gMass[m][node];
      }
    }
  } // patches
}

/*!----------------------------------------------------------------------
 * scheduleFindSurfaceParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleFindSurfaceParticles(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "SerialMPM::scheduleFindSurfaceParticles");

  Task* t = scinew Task("MPM::findSurfaceParticles",
                        this,
                        &SerialMPM::findSurfaceParticles);

  Ghost::GhostType gp = Ghost::AroundNodes;
  int ngc_p           = d_numGhostParticles;

  t->requires(Task::OldDW, d_mpmLabels->pSurfLabel, gp, ngc_p);
  t->computes(d_mpmLabels->pSurfLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
SerialMPM::findSurfaceParticles(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  auto numMPMMatls = d_materialManager->getNumMaterials("MPM");

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing findSurfaceParticles");

    for (int mat = 0; mat < numMPMMatls; mat++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", mat);
      int matID             = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<double> pSurf_old;
      ParticleVariable<double> pSurf;

      old_dw->get(pSurf_old, d_mpmLabels->pSurfLabel, pset);
      new_dw->allocateAndPut(pSurf, d_mpmLabels->pSurfLabel_preReloc, pset);

      // For now carry forward the particle surface data
      for (auto particle : *pset) {
        pSurf[particle] = pSurf_old[particle];
      }
    } // matl loop
  }   // patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeLogisticRegression
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeLogisticRegression(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (contactModel->useLogisticRegression()) {
    printSchedule(patches,
                  cout_doing,
                  "MPM::scheduleComputeLogisticRegression");

    Task* t = scinew Task("MPM::computeLogisticRegression",
                          this,
                          &SerialMPM::computeLogisticRegression);

    Ghost::GhostType gp = Ghost::AroundNodes;
    int ngc_p           = d_numGhostParticles;

    MaterialSubset* z_matl = scinew MaterialSubset();
    z_matl->add(0);
    z_matl->addReference();

    t->requires(Task::OldDW, d_mpmLabels->pXLabel, gp, ngc_p);
    t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gp, ngc_p);
    t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, gp, ngc_p);
    t->requires(Task::NewDW, d_mpmLabels->pSurfLabel_preReloc, gp, ngc_p);
    t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
    // t->requires(Task::OldDW, d_mpmLabels->NC_CCweightLabel, z_matl,
    // Ghost::None);

    t->computes(d_mpmLabels->gMatlProminenceLabel);
    t->computes(d_mpmLabels->gAlphaMaterialLabel);
    t->computes(d_mpmLabels->gNormAlphaToBetaLabel, z_matl);

    sched->addTask(t, patches, matls);

    if (z_matl->removeReference()) {
      delete z_matl; // shouln't happen, but...
    }
  }
}

void
SerialMPM::computeLogisticRegression(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{

  // As of 5/22/19, this uses John Nairn's and Chad Hammerquist's
  // Logistic Regression method for finding normals used in contact.
  // These "NormAlphaToBeta" are then used to find each material's
  // greatest prominence at a node.  That is, the portion of a
  // particle that projects farthest along the direction of that normal.
  // One material at each multi-material node is identified as the
  // alpha material, this is the material with the most mass at the node.
  // All other materials are beta materials, and the NormAlphaToBeta
  // is perpendicular to the plane separating those materials
  Ghost::GhostType gan   = Ghost::AroundNodes;
  Ghost::GhostType gnone = Ghost::None;

  auto numMPMMatls = d_materialManager->getNumMaterials("MPM");
  std::vector<constNCVariable<double>> gMass(numMPMMatls);
  std::vector<constParticleVariable<Point>> pX(numMPMMatls);

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              cout_doing,
              "Doing MPM::computeLogisticRegression");

    Vector dx = patch->dCell();

    // constNCVariable<double>  NC_CCweight;
    // old_dw->get(NC_CCweight, d_mpmLabels->NC_CCweightLabel,  0, patch, gnone,
    // 0);

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);
    std::vector<Vector> d_S(numInfluenceNodes);
    std::string interp_type = d_mpmFlags->d_interpolatorType;

    // Declare and allocate storage for use in the Logistic Regression
    NCVariable<int> gAlphaMaterial;
    NCVariable<int> gNumMatlsOnNode;
    NCVariable<int> gNumParticlesOnNode;
    NCVariable<Vector> gNormAlphaToBeta;
    std::vector<NCVariable<Int130>> gParticleList(numMPMMatls);

    new_dw->allocateAndPut(gAlphaMaterial,
                           d_mpmLabels->gAlphaMaterialLabel,
                           0,
                           patch);
    new_dw->allocateAndPut(gNormAlphaToBeta,
                           d_mpmLabels->gNormAlphaToBetaLabel,
                           0,
                           patch);
    new_dw->allocateTemporary(gNumMatlsOnNode, patch);
    new_dw->allocateTemporary(gNumParticlesOnNode, patch);
    gAlphaMaterial.initialize(-99);
    gNumMatlsOnNode.initialize(0);
    gNumParticlesOnNode.initialize(0);
    gNormAlphaToBeta.initialize(Vector(-99. - 99. - 99.));

    // Get out the mass first, need access to the mass of all materials
    // at each node.
    // Also, at each multi-material node, we need the position of the
    // particles around that node.  Rather than store the positions, for now,
    // store a list of particle indices for each material at each node and
    // use those to point into the particle set to get the particle positions
    for (int mat = 0; mat < numMPMMatls; mat++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", mat);
      int matID             = mpm_matl->getDWIndex();
      new_dw->get(gMass[mat], d_mpmLabels->gMassLabel, matID, patch, gnone, 0);
      new_dw->allocateTemporary(gParticleList[mat], patch);
    }

    // Here, find out two things:
    // 1.  How many materials have mass on a node
    // 2.  Which material has the most mass on a node.  That is the alpha matl.
    for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      double maxMass = -9.e99;
      for (int mat = 0; mat < numMPMMatls; mat++) {
        if (gMass[mat][node] > 1.e-16) {
          gNumMatlsOnNode[node]++;
          if (gMass[mat][node] > maxMass) {
            // This is the alpha material, all other matls are beta
            gAlphaMaterial[node] = mat;
            maxMass              = gMass[mat][node];
          }
        }
      } // Loop over materials
      if (gNumMatlsOnNode[node] < 2) {
        gAlphaMaterial[node] = -99;
      }
    } // Node Iterator

    // In this section of code, we find the particles that are in the
    // vicinity of a multi-material node and put their indices in a list
    // so we can retrieve their positions later.

    // I hope to improve on this gParticleList later, but for now,
    // the last element in the array holds the number of entries in the
    // array.  I don't yet know how to allocate an STL container on the nodes.
    for (int mat = 0; mat < numMPMMatls; mat++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", mat);
      int matID             = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);
      constParticleVariable<Matrix3> pSize, pDefGrad_old;

      old_dw->get(pX[mat], d_mpmLabels->pXLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);

      // Initialize the gParticleList
      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        for (auto particle = 0; particle < 400; particle++) {
          gParticleList[mat][node][particle] = 0;
        }
      }

      // Loop over particles and find which multi-mat nodes they contribute to
      for (auto idx : *pset) {

        interpolator->findCellAndWeightsAndShapeDerivatives(pX[mat][idx],
                                                            ni,
                                                            S,
                                                            d_S,
                                                            pSize[idx],
                                                            pDefGrad_old[idx]);

        std::set<IntVector> nodeList;
        for (int k = 0; k < numInfluenceNodes; k++) {
          auto node = ni[k];
          if (patch->containsNode(node) && gNumMatlsOnNode[node] > 1 &&
              S[k] > 1.e-100) {
            nodeList.insert(node);
          } // conditional
        }   // loop over nodes returned by interpolator

        for (auto node : nodeList) {
          auto& particle                     = gParticleList[mat][node][399];
          gParticleList[mat][node][particle] = idx;
          particle++;
          gNumParticlesOnNode[node]++;
        }
        nodeList.clear();
      } // Loop over Particles
    }

    // This is the Logistic Regression code that finds the normal to the
    // plane that separates two materials.  This is as directly as possible
    // from Nairn & Hammerquist, 2019.
    using Vector4  = Eigen::Matrix<double, 4, 1>;
    using Matrix44 = Eigen::Matrix<double, 4, 4>;

    double lam     = 1.e-7 * dx.x() * dx.x();
    Vector4 lambda = { lam, lam, lam, 0 };
    double wp      = 1.0;
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      // Only work on multi-material nodes
      if (gAlphaMaterial[node] >= 0) {
        bool converged = false;
        int num_iters  = 0;
        double tol     = 1.e-5;
        Vector4 phi    = { 1., 0., 0., 0. };
        Vector nhat_k(phi[0], phi[1], phi[2]);
        Vector nhat_backup(0.);
        double error_min = 1.0;
        double error     = 1.0;
        while (!converged) {
          num_iters++;

          Vector4 g_phi          = -lambda.cwiseProduct(phi);
          Matrix44 g_prime_phi   = Matrix44::Zero();
          g_prime_phi.diagonal() = lambda;

          // std::cout << "phi_k = " << phi.transpose() << "\n";
          for (int mat = 0; mat < numMPMMatls; mat++) {
            double cp = 0.;
            if (gAlphaMaterial[node] == mat) {
              cp = -1.;
            } else {
              cp = 1.;
            }

            for (int part = 0; part < gParticleList[mat][node][399]; part++) {
              Point xp    = pX[mat][gParticleList[mat][node][part]];
              Vector4 xp4 = { xp.x(), xp.y(), xp.z(), 1.0 };
              Matrix44 xx = xp4 * xp4.transpose();

              double theta    = xp4.dot(phi);
              double exptheta = std::exp(-theta);
              if (!std::isfinite(exptheta)) {
                theta    = xp4.dot(phi / phi.norm());
                exptheta = std::exp(-theta);
                // std::cout << "**INTERNAL ERROR** MPM: In logistic regression:
                // xp . phi too large.\n"; std::cout << "xp = " <<
                // xp4.transpose() << " phi = " << phi.transpose() << "\n"
                //           << "theta = " << theta
                //           << "alpha = " << alpha << " psi = " << psi << " f =
                //           " << fEq20 << "\n";
              }

              double alpha    = 1.0 + exptheta;
              double psi      = 2.0 * exptheta / (alpha * alpha);
              double fEq20    = 2.0 / alpha - 1.0;
              double cp_fEq20 = cp - fEq20;

              g_phi += xp4 * (wp * cp_fEq20 * psi);
              g_prime_phi += xx * (psi * psi * wp);

              // double psi_deriv = psi * (2.0 / alpha * exptheta - 1.0);
              // g_prime_phi += xx * (psi * psi * wp - cp_fEq20 * psi_deriv);

            } // Loop over each material's particle list
          }   // Loop over materials

          Eigen::MatrixXd phi_inc =
            g_prime_phi.colPivHouseholderQr().solve(g_phi);
          phi += phi_inc;

          /*
          std::cout << "iter = " << num_iters << "\n";
          std::cout << "g_phi = " << g_phi.transpose() << "\n";
          std::cout << "g_prime_phi = " << g_prime_phi << "\n";
          std::cout << "phi_inc = " << phi_inc.transpose() <<"\n";
          std::cout << "phi_k+1 = " << phi.transpose() << "\n";
          */

          Vector nhat_kp1(phi[0], phi[1], phi[2]);
          nhat_kp1 /= (nhat_kp1.length() + 1.e-100);
          error = 1.0 - Dot(nhat_kp1, nhat_k);
          if (error < error_min) {
            error_min   = error;
            nhat_backup = nhat_kp1;
          }
          if (error < tol || num_iters > 15) {
            converged = true;
            if (num_iters > 15) {
              gNormAlphaToBeta[node] = nhat_backup;
            } else {
              gNormAlphaToBeta[node] = nhat_kp1;
            }
          } else {
            nhat_k = nhat_kp1;
          }
        } // while(!converged) loop

      } // If this node has more than one particle on it
    }   // Loop over nodes

    MPMBoundCond bc;
    bc.setBoundaryCondition(patch,
                            0,
                            "Symmetric",
                            gNormAlphaToBeta,
                            interp_type);

    // Renormalize normal std::vectors after setting BCs
    for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      gNormAlphaToBeta[node] /= (gNormAlphaToBeta[node].length() + 1.e-100);
      if (gAlphaMaterial[node] == -99) {
        gNormAlphaToBeta[node] = Vector(0.);
      }
      if (!(gNormAlphaToBeta[node].length() >= 0.0)) {
        std::cout << "Node  = " << node << "\n";
        std::cout << "gNormAlphaToBeta[node] = " << gNormAlphaToBeta[node]
                  << "\n";
        std::cout << "gAlphaMaterial[node] = " << gAlphaMaterial[node] << "\n";
        std::cout << "gNumMatlsOnNode[node] = " << gNumMatlsOnNode[node]
                  << "\n";
        std::cout << "gNumParticlesOnNode[node] = " << gNumParticlesOnNode[node]
                  << "\n";
      }
    } // Loop over nodes

    // Loop over all the particles, find the nodes they interact with
    // For the alpha material (gAlphaMaterial) find g.position as the
    // maximum of the dot product between each particle corner and the
    // gNormalAlphaToBeta std::vector.
    // For the beta materials (every other material) find g.position as the
    // minimum of the dot product between each particle corner and the
    // gNormalAlphaToBeta std::vector.
    // Compute "MatlProminence" as the min/max of the dot product between
    // the normal and the particle corners

    std::vector<NCVariable<double>> d_x_p_dot_n(numMPMMatls);

    for (int mat = 0; mat < numMPMMatls; mat++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", mat);
      int matID             = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      constParticleVariable<Matrix3> pSize, pDefGrad_old;
      constParticleVariable<double> pSurf;

      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);
      new_dw->get(pSurf, d_mpmLabels->pSurfLabel_preReloc, pset);

      new_dw->allocateAndPut(d_x_p_dot_n[mat],
                             d_mpmLabels->gMatlProminenceLabel,
                             matID,
                             patch);

      d_x_p_dot_n[mat].initialize(-99.);

      NCVariable<double> gProjMax, gProjMin;
      new_dw->allocateTemporary(gProjMax, patch, gnone);
      new_dw->allocateTemporary(gProjMin, patch, gnone);
      gProjMax.initialize(-9.e99);
      gProjMin.initialize(9.e99);

      for (auto idx : *pset) {

        if (pSurf[idx] > 0.9) {
          interpolator->findCellAndWeights(pX[mat][idx],
                                           ni,
                                           S,
                                           pSize[idx],
                                           pDefGrad_old[idx]);

          Matrix3 curSize = pDefGrad_old[idx] * pSize[idx];
          Matrix3 dsize =
            curSize * Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]);

#if 0
          // This version uses particle corners to compute prominence
          // Compute std::vectors from particle center to the corners
          Vector RNL[8];
          RNL[0] = Vector(-dsize(0,0)-dsize(0,1)+dsize(0,2),
                          -dsize(1,0)-dsize(1,1)+dsize(1,2),
                          -dsize(2,0)-dsize(2,1)+dsize(2,2))*0.5;
          RNL[1] = Vector( dsize(0,0)-dsize(0,1)+dsize(0,2),
                           dsize(1,0)-dsize(1,1)+dsize(1,2),
                           dsize(2,0)-dsize(2,1)+dsize(2,2))*0.5;
          RNL[2] = Vector( dsize(0,0)+dsize(0,1)+dsize(0,2),
                           dsize(1,0)+dsize(1,1)+dsize(1,2),
                           dsize(2,0)+dsize(2,1)+dsize(2,2))*0.5;
          RNL[3] = Vector(-dsize(0,0)+dsize(0,1)+dsize(0,2),
                          -dsize(1,0)+dsize(1,1)+dsize(1,2),
                          -dsize(2,0)+dsize(2,1)+dsize(2,2))*0.5;
          RNL[4] = Vector(-dsize(0,0)-dsize(0,1)-dsize(0,2),
                          -dsize(1,0)-dsize(1,1)-dsize(1,2),
                          -dsize(2,0)-dsize(2,1)-dsize(2,2))*0.5;
          RNL[5] = Vector( dsize(0,0)-dsize(0,1)-dsize(0,2),
                           dsize(1,0)-dsize(1,1)-dsize(1,2),
                           dsize(2,0)-dsize(2,1)-dsize(2,2))*0.5;
          RNL[6] = Vector( dsize(0,0)+dsize(0,1)-dsize(0,2),
                           dsize(1,0)+dsize(1,1)-dsize(1,2),
                           dsize(2,0)+dsize(2,1)-dsize(2,2))*0.5;
          RNL[7] = Vector(-dsize(0,0)+dsize(0,1)-dsize(0,2),
                          -dsize(1,0)+dsize(1,1)-dsize(1,2),
                          -dsize(2,0)+dsize(2,1)-dsize(2,2))*0.5;

          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              if (S[k] > 0. && gNumParticlesOnNode[node] > 1) {
                for (int ic = 0; ic < 8; ic++) {
                  Vector xp_xi = (pX[mat][idx].asVector()+RNL[ic]);
                  double proj = Dot(xp_xi, gNormAlphaToBeta[node]);
                  if (mat == gAlphaMaterial[node]) {
                    if (proj > gProjMax[node]) {
                       gProjMax[node] = proj;
                       d_x_p_dot_n[mat][node] = proj;
                    }
                  } else {
                    if (proj < gProjMin[node]) {
                       gProjMin[node] = proj;
                       d_x_p_dot_n[mat][node] = proj;
                    }
                  }
                } // Loop over all 8 particle corners
              }  // Only deal with nodes that this particle affects
            }  // If node is on the patch
          } // Loop over nodes near this particle
#endif
#if 0
          // This version uses constant particle radius, here assuming 2 PPC
          // in each direction
          double Rp = 0.25*dx.x();
          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              if (S[k] > 0. && gNumParticlesOnNode[node] > 1) {
                for (int ic = 0; ic < 8; ic++) {
                  Vector xp_xi = (pX[mat][idx].asVector());
                  double proj = Dot(xp_xi, gNormAlphaToBeta[node]);
                  if (mat == gAlphaMaterial[node]) {
                    proj += Rp;
                    if (proj > gProjMax[node]) {
                      gProjMax[node] = proj;
                      d_x_p_dot_n[mat][node] = proj;
                    }
                  } else {
                    proj -= Rp;
                    if (proj < gProjMin[node]) {
                      gProjMin[node] = proj;
                      d_x_p_dot_n[mat][node] = proj;
                    }
                  }
                } // Loop over all 8 particle corners
              }  // Only deal with nodes that this particle affects
            }  // If node is on the patch
          } // Loop over nodes near this particle
#endif

#if 1
          // This version uses particle faces to compute prominence.
          // Compute std::vectors from particle center to the faces
          Vector RFL[6];
          RFL[0] = Vector(-dsize(0, 0), -dsize(1, 0), -dsize(2, 0)) * 0.5;
          RFL[1] = Vector(dsize(0, 0), dsize(1, 0), dsize(2, 0)) * 0.5;
          RFL[2] = Vector(-dsize(0, 1), -dsize(1, 1), -dsize(2, 1)) * 0.5;
          RFL[3] = Vector(dsize(0, 1), dsize(1, 1), dsize(2, 1)) * 0.5;
          RFL[4] = Vector(-dsize(0, 2), -dsize(1, 2), -dsize(2, 2)) * 0.5;
          RFL[5] = Vector(dsize(0, 2), dsize(1, 2), dsize(2, 2)) * 0.5;

          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              if (S[k] > 0. && gNumParticlesOnNode[node] > 1) {
                for (int ic = 0; ic < 6; ic++) {
                  Vector xp_xi = pX[mat][idx].asVector() + RFL[ic];
                  double proj  = Dot(xp_xi, gNormAlphaToBeta[node]);
                  if (mat == gAlphaMaterial[node]) {
                    if (proj > gProjMax[node]) {
                      gProjMax[node]         = proj;
                      d_x_p_dot_n[mat][node] = proj;
                    }
                  } else {
                    if (proj < gProjMin[node]) {
                      gProjMin[node]         = proj;
                      d_x_p_dot_n[mat][node] = proj;
                    }
                  }
                } // Loop over all 8 particle corners
              }   // Only deal with nodes that this particle affects
            }     // If node is on the patch
          }       // Loop over nodes near this particle
#endif
        } // Is a surface particle
      }   // end Particle loop
    }     // loop over matls

  } // patches
}

/*!----------------------------------------------------------------------
 * scheduleExMomInterpolated
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleExMomInterpolated(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleExMomInterpolated");

  contactModel->addComputesAndRequires(sched,
                                       patches,
                                       matls,
                                       d_mpmLabels->gVelocityLabel);
}

/*!----------------------------------------------------------------------
 * scheduleUpdateCohesiveZones
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleUpdateCohesiveZones(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSubset* mpm_matls,
                                       const MaterialSubset* cz_matls,
                                       const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleUpdateCohesiveZones");

  Task* t = scinew Task("MPM::updateCohesiveZones",
                        this,
                        &SerialMPM::updateCohesiveZones);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityLabel,
              mpm_matls,
              gac,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gMassLabel,
              mpm_matls,
              gac,
              d_numGhostNodes);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czLengthLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czNormLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czTangLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czDispTopLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czDispBottomLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czSeparationLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czForceLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czTopMatLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czBotMatLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czFailedLabel, cz_matls, gnone);
  t->requires(Task::OldDW, d_mpmLabels->czIDLabel, cz_matls, gnone);

  t->computes(d_mpmLabels->pXLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czLengthLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czNormLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czTangLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czDispTopLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czDispBottomLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czSeparationLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czForceLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czTopMatLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czBotMatLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czFailedLabel_preReloc, cz_matls);
  t->computes(d_mpmLabels->czIDLabel_preReloc, cz_matls);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * updateCohesiveZones
 *   The following is adapted from "Simulation of dynamic crack growth
 *   using the generalized interpolation material point (GIMP) method"
 *   Daphalapurkar, N.P., et al., Int. J. Fracture, 143, 79-102, 2007.
 *-----------------------------------------------------------------------*/
void
SerialMPM::updateCohesiveZones(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing updateCohesiveZones");

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    std::vector<constNCVariable<Vector>> gVelocity(numMPMMatls);
    std::vector<constNCVariable<double>> gMass(numMPMMatls);
    // double rho_init[numMPMMatls];
    Vector dx = patch->dCell();
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // rho_init[m]=mpm_matl->getInitialDensity();
      Ghost::GhostType gac = Ghost::AroundCells;
      new_dw->get(gVelocity[m],
                  d_mpmLabels->gVelocityLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostNodes);
      new_dw->get(gMass[m],
                  d_mpmLabels->gMassLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostNodes);
    }

    /*
      double time = d_materialManager->getElapsedTime();
      std::string outfile_name = "force_sep.dat";
      ofstream dest;
      dest.open(outfile_name.c_str(),ios::app);
      if(!dest){
      std::cerr << "File " << outfile_name << " can't be opened." << "\n";
      }
    */

    int numCZMatls = d_materialManager->getNumMatls("CohesizeZone");
    for (int m = 0; m < numCZMatls; m++) {
      CZMaterial* cz_matl = d_materialManager->getMaterial("CohesiveZone", m);
      int matID           = cz_matl->getDWIndex();

      // Not populating the delset, but we need this to satisfy Relocate
      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);
      new_dw->deleteParticles(delset);

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> czx;
      ParticleVariable<Point> czx_new;
      constParticleVariable<double> czlength;
      ParticleVariable<double> czlength_new;
      constParticleVariable<long64> czids;
      ParticleVariable<long64> czids_new;
      constParticleVariable<Vector> cznorm, cztang, czDispTop;
      ParticleVariable<Vector> cznorm_new, cztang_new, czDispTop_new;
      constParticleVariable<Vector> czDispBot, czsep, czforce;
      ParticleVariable<Vector> czDispBot_new, czsep_new, czforce_new;
      constParticleVariable<int> czTopMat, czBotMat, czFailed;
      ParticleVariable<int> czTopMat_new, czBotMat_new, czFailed_new;

      old_dw->get(czx, d_mpmLabels->pXLabel, pset);
      old_dw->get(czlength, d_mpmLabels->czLengthLabel, pset);
      old_dw->get(cznorm, d_mpmLabels->czNormLabel, pset);
      old_dw->get(cztang, d_mpmLabels->czTangLabel, pset);
      old_dw->get(czDispTop, d_mpmLabels->czDispTopLabel, pset);
      old_dw->get(czDispBot, d_mpmLabels->czDispBottomLabel, pset);
      old_dw->get(czsep, d_mpmLabels->czSeparationLabel, pset);
      old_dw->get(czforce, d_mpmLabels->czForceLabel, pset);
      old_dw->get(czids, d_mpmLabels->czIDLabel, pset);
      old_dw->get(czTopMat, d_mpmLabels->czTopMatLabel, pset);
      old_dw->get(czBotMat, d_mpmLabels->czBotMatLabel, pset);
      old_dw->get(czFailed, d_mpmLabels->czFailedLabel, pset);

      new_dw->allocateAndPut(czx_new, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(czlength_new,
                             d_mpmLabels->czLengthLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(cznorm_new,
                             d_mpmLabels->czNormLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(cztang_new,
                             d_mpmLabels->czTangLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czDispTop_new,
                             d_mpmLabels->czDispTopLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czDispBot_new,
                             d_mpmLabels->czDispBottomLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czsep_new,
                             d_mpmLabels->czSeparationLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czforce_new,
                             d_mpmLabels->czForceLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czids_new, d_mpmLabels->czIDLabel_preReloc, pset);
      new_dw->allocateAndPut(czTopMat_new,
                             d_mpmLabels->czTopMatLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czBotMat_new,
                             d_mpmLabels->czBotMatLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(czFailed_new,
                             d_mpmLabels->czFailedLabel_preReloc,
                             pset);

      czlength_new.copyData(czlength);
      czids_new.copyData(czids);
      czTopMat_new.copyData(czTopMat);
      czBotMat_new.copyData(czBotMat);

      double sig_max  = cz_matl->getCohesiveNormalStrength();
      double delta_n  = cz_matl->getCharLengthNormal();
      double tau_max  = cz_matl->getCohesiveTangentialStrength();
      double delta_t  = cz_matl->getCharLengthTangential();
      double delta_s  = delta_t;
      bool rotate_CZs = cz_matl->getDoRotation();

      double phi_n = M_E * sig_max * delta_n;
      double phi_t = sqrt(M_E / 2) * tau_max * delta_t;
      double q     = phi_t / phi_n;
      // From the text following Eq. 15 in Nitin's paper it is a little hard
      // to tell what r should be, but zero seems like a reasonable value
      // based on the example problem in that paper
      double r = 0.;

      // Loop over particles
      for (auto idx : *pset) {

        //        double length = sqrt(czlength[idx]);
        //        Vector size(length,length,length);
        Matrix3 size(0.1, 0., 0., 0., 0.1, 0., 0., 0., 0.1);
        Matrix3 defgrad;
        defgrad.Identity();

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(czx[idx], ni, S, size, defgrad);

        Vector velTop(0.0, 0.0, 0.0);
        Vector velBot(0.0, 0.0, 0.0);
        double massTop     = 0.0;
        double massBot     = 0.0;
        double mass_ratio  = 0.0;
        int TopMat         = czTopMat[idx];
        int BotMat         = czBotMat[idx];
        double cell_volume = dx.x() * dx.y() * dx.z();
        // double denseTop = rho_init[TopMat];
        // double denseBot = rho_init[BotMat];
        double TOPMAX = 0.0;
        double BOTMAX = 0.0;

        //      if (denseBot != denseTop){
        //         throw ProblemSetupException("Different densities not allowed
        //         for Bottom and Top Material of Cohesive Zone",
        //                                 __FILE__, __LINE__);
        //      }

        // double density_ratio = denseTop/denseBot;
        //  Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          velTop += gVelocity[TopMat][node] * S[k];
          velBot += gVelocity[BotMat][node] * S[k];
          massTop += gMass[TopMat][node] * S[k];
          TOPMAX += cell_volume;
          massBot += gMass[BotMat][node] * S[k];
          BOTMAX += cell_volume;
        }
        massTop = massTop / TOPMAX;
        massBot = massBot / BOTMAX;
        if (massBot > 0.0) {
          mass_ratio = massTop / massBot;
          mass_ratio = std::min(mass_ratio, 1.0 / mass_ratio);
        } else {
          mass_ratio = 0.0;
        }
        double mass_correction_factor = mass_ratio;

        // Update the cohesive zone's position and displacements
        czx_new[idx]       = czx[idx] + .5 * (velTop + velBot) * delT;
        czDispTop_new[idx] = czDispTop[idx] + velTop * delT;
        czDispBot_new[idx] = czDispBot[idx] + velBot * delT;
        czsep_new[idx]     = czDispTop_new[idx] - czDispBot_new[idx];

        double disp = czsep_new[idx].length();
        if (disp > 0.0 && rotate_CZs) {
          Matrix3 Rotation;
          Matrix3 Rotation_tang;
          cz_matl->computeRotationMatrix(Rotation,
                                         Rotation_tang,
                                         cznorm[idx],
                                         czsep_new[idx]);

          cznorm_new[idx] = Rotation * cznorm[idx];
          cztang_new[idx] = Rotation_tang * cztang[idx];
        } else {
          cznorm_new[idx] = cznorm[idx];
          cztang_new[idx] = cztang[idx];
        }

        Vector cztang2 = Cross(cztang_new[idx], cznorm_new[idx]);

        double D_n  = Dot(czsep_new[idx], cznorm_new[idx]);
        double D_t1 = Dot(czsep_new[idx], cztang_new[idx]);
        double D_t2 = Dot(czsep_new[idx], cztang2);

        // Determine if a CZ has failed.  Currently harmatIDring failure
        // criteria to fail zone if normal sep is > 4*delta_n or 2*delta_t
        double czf = 0.0;
        if (czFailed[idx] > 0) {
          czFailed_new[idx] = czFailed[idx];
          czf               = 1.0;
        } else if (D_n > 4.0 * delta_n) {
          czFailed_new[idx] = 1;
          czf               = 1.0;
        } else if (fabs(D_t1) > 2.0 * delta_t) {
          czFailed_new[idx] = 2;
          czf               = 1.0;
        } else if (fabs(D_t2) > 2.0 * delta_s) {
          czFailed_new[idx] = 2;
          czf               = 1.0;
        } else {
          czFailed_new[idx] = 0;
        }

        double normal_stress =
          (phi_n / delta_n) * exp(-D_n / delta_n) *
          ((D_n / delta_n) * exp((-D_t1 * D_t1) / (delta_t * delta_t)) +
           ((1. - q) / (r - 1.)) *
             (1. - exp(-D_t1 * D_t1 / (delta_t * delta_t))) *
             (r - D_n / delta_n));

        double tang1_stress =
          (phi_n / delta_n) * (2. * delta_n / delta_t) * (D_t1 / delta_t) *
          (q + ((r - q) / (r - 1.)) * (D_n / delta_n)) * exp(-D_n / delta_n) *
          exp(-D_t1 * D_t1 / (delta_t * delta_t));

        double tang2_stress =
          (phi_n / delta_n) * (2. * delta_n / delta_s) * (D_t2 / delta_s) *
          (q + ((r - q) / (r - 1.)) * (D_n / delta_n)) * exp(-D_n / delta_n) *
          exp(-D_t2 * D_t2 / (delta_s * delta_s));

        czforce_new[idx] =
          mass_correction_factor *
          (normal_stress * cznorm_new[idx] * czlength_new[idx] +
           tang1_stress * cztang_new[idx] * czlength_new[idx] +
           tang2_stress * cztang2 * czlength_new[idx]) *
          (1.0 - czf);

        /*
          dest << time << " " << czsep_new[idx].x() << " " << czsep_new[idx].y()
          << " " << czforce_new[idx].x() << " " << czforce_new[idx].y() << "\n";
          if(fabs(normal_force) >= 0.0){
          std::cout << "czx_new " << czx_new[idx] << "\n";
          std::cout << "czforce_new " << czforce_new[idx] << "\n";
          std::cout << "czsep_new " << czsep_new[idx] << "\n";
          std::cout << "czDispTop_new " << czDispTop_new[idx] << "\n";
          std::cout << "czDispBot_new " << czDispBot_new[idx] << "\n";
          std::cout << "velTop " << velTop << "\n";
          std::cout << "velBot " << velBot << "\n";
          std::cout << "delT " << delT << "\n";
          }
        */
      }
    }

    // delete interpolator;
  }
}
/*!----------------------------------------------------------------------
 * scheduleAddCohesiveZoneForces
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleAddCohesiveZoneForces(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSubset* mpm_matls,
                                         const MaterialSubset* cz_matls,
                                         const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleAddCohesiveZoneForces");

  Task* t = scinew Task("MPM::addCohesiveZoneForces",
                        this,
                        &SerialMPM::addCohesiveZoneForces);

  Ghost::GhostType gan = Ghost::AroundNodes;
  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              cz_matls,
              gan,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->czLengthLabel_preReloc,
              cz_matls,
              gan,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->czForceLabel_preReloc,
              cz_matls,
              gan,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->czTopMatLabel_preReloc,
              cz_matls,
              gan,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->czBotMatLabel_preReloc,
              cz_matls,
              gan,
              d_numGhostParticles);

  t->modifies(d_mpmLabels->gExternalForceLabel, mpm_matls);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * addCohesiveZoneForces
 *-----------------------------------------------------------------------*/
void
SerialMPM::addCohesiveZoneForces(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing addCohesiveZoneForces");

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    std::vector<NCVariable<Vector>> gext_force(numMPMMatls);
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* cz_matl = d_materialManager->getMaterial("MPM", m);
      int matID            = cz_matl->getDWIndex();

      new_dw->getModifiable(gext_force[m],
                            d_mpmLabels->gExternalForceLabel,
                            matID,
                            patch);
    }

    Ghost::GhostType gan = Ghost::AroundNodes;
    int numCZMatls       = d_materialManager->getNumMatls("CohesizeZone");
    for (int m = 0; m < numCZMatls; m++) {
      CZMaterial* cz_matl = d_materialManager->getMaterial("CohesiveZone", m);
      int matID           = cz_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      // Get the arrays of particle values to be changed
      constParticleVariable<Point> czx;
      constParticleVariable<double> czlength;
      constParticleVariable<Vector> czforce;
      constParticleVariable<int> czTopMat, czBotMat;
      constParticleVariable<Matrix3> pDefGrad;

      old_dw->get(czx, d_mpmLabels->pXLabel, pset);
      new_dw->get(czlength, d_mpmLabels->czLengthLabel_preReloc, pset);
      new_dw->get(czforce, d_mpmLabels->czForceLabel_preReloc, pset);
      new_dw->get(czTopMat, d_mpmLabels->czTopMatLabel_preReloc, pset);
      new_dw->get(czBotMat, d_mpmLabels->czBotMatLabel_preReloc, pset);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        //        double length = sqrt(czlength[idx]);
        Matrix3 size(0.1, 0., 0., 0., 0.1, 0., 0., 0., 0.1);
        Matrix3 defgrad;
        defgrad.Identity();

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(czx[idx], ni, S, size, defgrad);

        int TopMat = czTopMat[idx];
        int BotMat = czBotMat[idx];

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          if (patch->containsNode(node)) {
            gext_force[BotMat][node] =
              gext_force[BotMat][node] + czforce[idx] * S[k];
            gext_force[TopMat][node] =
              gext_force[TopMat][node] - czforce[idx] * S[k];
          }
        }
      }
    }
    // delete interpolator;
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeContactArea
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeContactArea(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /** computeContactArea */
  if (d_boundaryTractionFaces.size() > 0) {

    printSchedule(patches, cout_doing, "MPM::scheduleComputeContactArea");
    Task* t = scinew Task("MPM::computeContactArea",
                          this,
                          &SerialMPM::computeContactArea);

    Ghost::GhostType gnone = Ghost::None;
    t->requires(Task::NewDW, d_mpmLabels->gVolumeLabel, gnone);
    for (auto face : d_boundaryTractionFaces) {
      int iface = (int)face;
      t->computes(d_mpmLabels->BndyContactCellAreaLabel[iface]);
    }
    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * computeContactArea
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeContactArea(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset*,
                              DataWarehouse* /*old_dw*/,
                              DataWarehouse* new_dw)
{
  // six indices for each of the faces
  double bndyCArea[6] = { 0, 0, 0, 0, 0, 0 };

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeContactArea");

    Vector dx = patch->dCell();

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      constNCVariable<double> gVolume;

      new_dw
        ->get(gVolume, d_mpmLabels->gVolumeLabel, matID, patch, Ghost::None, 0);

      for (auto face : d_boundaryTractionFaces) {
        int iface = (int)(face);

        // Check if the face is on an external boundary
        if (patch->getBCType(face) == Patch::Neighbor) {
          continue;
        }

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,

        // loop over face nodes to find boundary areas
        // Because this calculation uses gVolume, particle volumes interpolated
        // to the nodes, it will give 1/2 the expected value because the
        // particle values are distributed to all nodes, not just those on this
        // face.  It would require particles on the other side of the face to
        // "fill" the nodal volumes and give the correct area when divided by
        // the face normal cell dimension (celldepth). To correct for this,
        // nodearea incorporates a factor of two.

        IntVector projlow, projhigh;
        patch->getFaceNodes(face, 0, projlow, projhigh);
        const double celldepth = dx[iface / 2];

        for (int i = projlow.x(); i < projhigh.x(); i++) {
          for (int j = projlow.y(); j < projhigh.y(); j++) {
            for (int k = projlow.z(); k < projhigh.z(); k++) {
              IntVector ijk(i, j, k);
              double nodearea = 2.0 * gVolume[ijk] / celldepth; // node area
              bndyCArea[iface] += nodearea;
            }
          }
        }
      } // faces
    }   // materials
  }     // patches

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (auto face : d_boundaryTractionFaces) {
    int iface = (int)face;
    new_dw->put(sum_vartype(bndyCArea[iface]),
                d_mpmLabels->BndyContactCellAreaLabel[iface]);
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeInternalForce
 *
 * computeInternalForce
 *   in(P.CONMOD, P.NAT_X, P.VOLUME)
 *   operation(evaluate the divergence of the stress (stored in
 *   P.CONMOD) using P.NAT_X and the gradients of the
 *   shape functions)
 * out(G.F_INTERNAL)
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeInternalForce");

  Task* t = scinew Task("MPM::computeInternalForce",
                        this,
                        &SerialMPM::computeInternalForce);

  Ghost::GhostType gan   = Ghost::AroundNodes;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, d_mpmLabels->gVolumeLabel, gnone);
  t->requires(Task::NewDW,
              d_mpmLabels->gVolumeLabel,
              d_materialManager->getAllInOneMaterial(),
              Task::OutOfDomain,
              gnone);
  t->requires(Task::OldDW, d_mpmLabels->pStressLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pVolumeLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              gan,
              d_numGhostParticles);
#ifdef DEBUG_WITH_PARTICLE_ID
  t->requires(Task::OldDW,
              d_mpmLabels->pParticleIDLabel,
              gan,
              d_numGhostParticles);
#endif

  if (d_mpmFlags->d_withICE) {
    t->requires(Task::NewDW,
                d_mpmLabels->pPressureLabel,
                gan,
                d_numGhostParticles);
  }

  if (d_mpmFlags->d_artificialViscosity) {
    t->requires(Task::OldDW, d_mpmLabels->p_qLabel, gan, d_numGhostParticles);
  }

  t->computes(d_mpmLabels->gInternalForceLabel);

  for (auto face : d_boundaryTractionFaces) {
    int iface = (int)face;
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

/*!----------------------------------------------------------------------
 * computeInternalForce
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeInternalForce(const ProcessorGroup*,
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

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeInternalForce");

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();
    Matrix3 Id;
    Id.Identity();

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);
    std::vector<Vector> d_S(numInfluenceNodes);
    std::string interp_type = d_mpmFlags->d_interpolatorType;

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    NCVariable<Matrix3> gStressglobal;
    constNCVariable<double> gVolumeglobal;
    new_dw->get(gVolumeglobal,
                d_mpmLabels->gVolumeLabel,
                d_materialManager->getAllInOneMaterial()->get(0),
                patch,
                Ghost::None,
                0);
    new_dw->allocateAndPut(gStressglobal,
                           d_mpmLabels->gStressForSavingLabel,
                           d_materialManager->getAllInOneMaterial()->get(0),
                           patch);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // Create arrays for the particle position, volume
      // and the constitutive model
      constParticleVariable<Point> pX;
      constParticleVariable<double> pVol;
      constParticleVariable<double> p_pressure;
      constParticleVariable<double> p_q;
      constParticleVariable<Matrix3> pStress;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad_old;
      NCVariable<Vector> gInternalForce;
      NCVariable<Matrix3> gStress;
      constNCVariable<double> gVolume;

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pVol, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pStress, d_mpmLabels->pStressLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);

#ifdef DEBUG_WITH_PARTICLE_ID
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
#endif

      new_dw
        ->get(gVolume, d_mpmLabels->gVolumeLabel, matID, patch, Ghost::None, 0);

      new_dw->allocateAndPut(gStress,
                             d_mpmLabels->gStressForSavingLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gInternalForce,
                             d_mpmLabels->gInternalForceLabel,
                             matID,
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

      gInternalForce.initialize(Vector(0, 0, 0));

      Matrix3 stressvol;
      Matrix3 stresspress;

      // for the non axisymmetric case:
      if (!d_mpmFlags->d_axisymmetric) {
        for (auto idx : *pset) {

          // Get the node indices that surround the cell
          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[idx],
            ni,
            S,
            d_S,
            pSize[idx],
            pDefGrad_old[idx]);
          stressvol   = pStress[idx] * pVol[idx];
          stresspress = pStress[idx] + Id * (p_pressure[idx] - p_q[idx]);
          // std::cerr << " idx = " << idx << " pStress = " << pStress[idx] <<
          // "\n";

          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              Vector div(d_S[k].x() * oodx[0],
                         d_S[k].y() * oodx[1],
                         d_S[k].z() * oodx[2]);
              gInternalForce[node] -= (div * stresspress) * pVol[idx];
              gStress[node] += stressvol * S[k];

#ifdef DEBUG_WITH_PARTICLE_ID
              if (pParticleID[idx] == testParticleID) {
                if (node == IntVector(3, 38, 0)) {
                  proc0cout << "Particle ID = " << pParticleID[idx]
                            << " node = " << node << " dS = " << d_S[k]
                            << " div = " << div << " stress = " << pStress[idx]
                            << " damp = " << p_q[idx]
                            << " stresspress = " << stresspress
                            << " vol = " << pVol[idx]
                            << " fint_g = " << gInternalForce[node] << "\n";
                }
              }
#endif
#ifdef CHECK_ISFINITE
              if (!std::isfinite(gInternalForce[node].x()) ||
                  !std::isfinite(gInternalForce[node].y()) ||
                  !std::isfinite(gInternalForce[node].z())) {
                std::cout << "vol = " << pVol[idx] << " node = " << node
                          << " f_i = " << gInternalForce[node]
                          << " sig_g = " << gStress[node]
                          << " sig_p = " << stresspress << "\n";
              }
#endif
            }
          }
        }
      }

      // for the axisymmetric case
      if (d_mpmFlags->d_axisymmetric) {
        for (auto part : *pset) {

          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[part],
            ni,
            S,
            d_S,
            pSize[part],
            pDefGrad_old[part]);

          stressvol   = pStress[part] * pVol[part];
          stresspress = pStress[part] + Id * (p_pressure[part] - p_q[part]);

#ifdef CHECK_ISFINITE
          if (!std::isfinite(stresspress(0, 0)) ||
              !std::isfinite(stresspress(0, 1)) ||
              !std::isfinite(stresspress(0, 2)) ||
              !std::isfinite(stresspress(1, 1)) ||
              !std::isfinite(stresspress(1, 2)) ||
              !std::isfinite(stresspress(2, 2))) {
            std::cout << " p_pressure = " << p_pressure[part]
                      << " p_q = " << p_q[part] << "\n";
          }
#endif

          // r is the x direction, z (axial) is the y direction
          double IFr = 0., IFz = 0.;
          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              IFr = d_S[k].x() * oodx[0] * stresspress(0, 0) +
                    d_S[k].y() * oodx[1] * stresspress(0, 1) +
                    d_S[k].z() * stresspress(2, 2);
              IFz = d_S[k].x() * oodx[0] * stresspress(0, 1) +
                    d_S[k].y() * oodx[1] * stresspress(1, 1);
              gInternalForce[node] -= Vector(IFr, IFz, 0.0) * pVol[part];
              gStress[node] += stressvol * S[k];
#ifdef CHECK_ISFINITE
              if (!std::isfinite(gInternalForce[node].x()) ||
                  !std::isfinite(gInternalForce[node].y()) ||
                  !std::isfinite(gInternalForce[node].z())) {
                std::cout << "vol = " << pVol[part] << " node = " << ni[k]
                          << " f_i = " << gInternalForce[node]
                          << " sig_g = " << gStress[node]
                          << " sig_p = " << stresspress << " IFr = " << IFr
                          << " IFz = " << IFz << "\n";
              }
#endif
#ifdef DEBUG_WITH_PARTICLE_ID
              // if (pParticleID[part] == testParticleID) {
              if (node == IntVector(3, 38, 0)) {
                proc0cout << "Particle ID = " << pParticleID[part]
                          << " node = " << node << " dS = " << d_S[k]
                          << " IFr = " << IFr << " IFz = " << IFz
                          << " stress = " << pStress[part]
                          << " damp = " << p_q[part]
                          << " stresspress = " << stresspress
                          << " vol = " << pVol[part]
                          << " fint_g = " << gInternalForce[node] << "\n";
              }
              //}
#endif
            }
          }
        }
      }

      for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        gStressglobal[c] += gStress[c];
        gStress[c] /= gVolume[c];
      }

      // save boundary forces before apply symmetry boundary condition.
      for (auto face : d_boundaryTractionFaces) {

        // Check if the face is on an external boundary
        if (patch->getBCType(face) == Patch::Neighbor) {
          continue;
        }

        const int iface = (int)face;

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,

        IntVector projlow, projhigh;
        patch->getFaceNodes(face, 0, projlow, projhigh);
        Vector norm      = face_norm(face);
        double celldepth = dx[iface / 2]; // length in dir. perp. to boundary

        // loop over face nodes to find boundary forces, ave. stress (traction).
        // Note that nodearea incorporates a factor of two as described in the
        // bndyCellArea calculation in order to get node face areas.

        for (int i = projlow.x(); i < projhigh.x(); i++) {
          for (int j = projlow.y(); j < projhigh.y(); j++) {
            for (int k = projlow.z(); k < projhigh.z(); k++) {
              IntVector ijk(i, j, k);

              // flip sign so that pushing on boundary gives positive force
              bndyForce[iface] -= gInternalForce[ijk];

              double nodearea = 2.0 * gVolume[ijk] / celldepth; // node area
              for (int ic = 0; ic < 3; ic++) {
                for (int jc = 0; jc < 3; jc++) {
                  bndyTraction[iface][ic] +=
                    gStress[ijk](ic, jc) * norm[jc] * nodearea;
                }
              }
            }
          }
        }
      } // faces

#ifdef DEBUG_WITH_PARTICLE_ID
      IntVector node(3, 38, 0);
      if (patch->containsNode(node)) {
        proc0cout << "Before BC: Material = " << m << " Node = " << node
                  << " fint_g = " << gInternalForce[node] << "\n";
      }
#endif
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              matID,
                              "Symmetric",
                              gInternalForce,
                              interp_type);
#ifdef DEBUG_WITH_PARTICLE_ID
      // IntVector node(3,38,0);
      if (patch->containsNode(node)) {
        proc0cout << "After BC: Material = " << m << " Node = " << node
                  << " fint_g = " << gInternalForce[node] << "\n";
      }
#endif

      // for (NodeIterator iter = patch->getNodeIterator(); !iter.done();
      // iter++) {
      //   std::cout << "After internal force: node = " << *iter
      //             << " gInternalForce = " << gInternalForce[*iter] << "\n";
      // }
    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gStressglobal[c] /= gVolumeglobal[c];
    }
    // delete interpolator;
  }

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (auto face : d_boundaryTractionFaces) {
    int iface = (int)face;
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

    // Use the face force and traction calculations to provide a second estimate
    // of the contact area.
    double bndyContactArea_iface = bndyContactCellArea_iface;
    if (bndyTraction[iface][iface / 2] * bndyTraction[iface][iface / 2] >
        1.e-12) {
      bndyContactArea_iface =
        bndyForce[iface][iface / 2] / bndyTraction[iface][iface / 2];
    }

    new_dw->put(sum_vartype(bndyContactArea_iface),
                d_mpmLabels->BndyContactAreaLabel[iface]);
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeAndIntegrateacceleration
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
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
                        &SerialMPM::computeAndIntegrateAcceleration);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gBodyForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gExternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  t->computes(d_mpmLabels->gVelocityStarLabel);
  t->computes(d_mpmLabels->gAccelerationLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeAndIntegrateacceleration
 *-----------------------------------------------------------------------*/
void
SerialMPM::computeAndIntegrateAcceleration(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset*,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing computeAndIntegrateAcceleration");

    Ghost::GhostType gnone = Ghost::None;
    // Vector gravity = d_mpmFlags->d_gravity;
    for (int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<Vector> gInternalForce, gBodyForce, gExternalForce,
        gVelocity;
      constNCVariable<double> gMass;

      delt_vartype delT;
      old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

      new_dw->get(gInternalForce,
                  d_mpmLabels->gInternalForceLabel,
                  matID,
                  patch,
                  gnone,
                  0);
      new_dw
        ->get(gBodyForce, d_mpmLabels->gBodyForceLabel, matID, patch, gnone, 0);
      new_dw->get(gExternalForce,
                  d_mpmLabels->gExternalForceLabel,
                  matID,
                  patch,
                  gnone,
                  0);
      new_dw->get(gMass, d_mpmLabels->gMassLabel, matID, patch, gnone, 0);
      new_dw
        ->get(gVelocity, d_mpmLabels->gVelocityLabel, matID, patch, gnone, 0);

      // Create variables for the results
      NCVariable<Vector> gVelocity_star, gAcceleration;
      new_dw->allocateAndPut(gVelocity_star,
                             d_mpmLabels->gVelocityStarLabel,
                             matID,
                             patch);
      new_dw->allocateAndPut(gAcceleration,
                             d_mpmLabels->gAccelerationLabel,
                             matID,
                             patch);

      gAcceleration.initialize(Vector(0., 0., 0.));
      double damp_coef = d_mpmFlags->d_artificialDampCoeff;

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;

        Vector acc(0., 0., 0.);
        if (gMass[c] > d_mpmFlags->d_minMassForAcceleration) {
          acc =
            (gInternalForce[c] + gExternalForce[c] + gBodyForce[c]) / gMass[c];
          acc -= damp_coef * gVelocity[c];
        }
        // gAcceleration[c] = acc +  gravity;
        gAcceleration[c]  = acc;
        gVelocity_star[c] = gVelocity[c] + gAcceleration[c] * delT;
// std::cout << "After acceleration: material = " << m << " node = " << c
//           << " gMass = " << gMass[c]
//           << " gAcceleration = " << gAcceleration[c] << "\n";
#ifdef CHECK_ISFINITE
        if (!std::isfinite(gAcceleration[c].x()) ||
            !std::isfinite(gAcceleration[c].y()) ||
            !std::isfinite(gAcceleration[c].z())) {
          std::cout << " node = " << c << " f_i = " << gInternalForce[c]
                    << " f_e = " << gExternalForce[c]
                    << " f_b = " << gBodyForce[c] << " m = " << gMass[c]
                    << " v = " << gVelocity[c] << "\n";
        }
#endif
#ifdef DEBUG_WITH_PARTICLE_ID
        IntVector node(3, 38, 0);
        if (c == node) {
          proc0cout << "Node = " << node << " fint_g = " << gInternalForce[node]
                    << " fext_g = " << gExternalForce[node]
                    << " fbod_g = " << gBodyForce[node]
                    << " acc = " << gAcceleration[node] << "\n";
        }
#endif
      }
    } // matls
  }
}

/*!----------------------------------------------------------------------
 * scheduleExMomIntegrated
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleExMomIntegrated(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* exMomIntegrated
   *   in(G.MASS, G.VELOCITY_STAR, G.ACCELERATION)
   *   operation(peform operations which will cause each of
   *              velocity fields to feel the influence of the
   *              the others according to specific rules)
   *   out(G.VELOCITY_STAR, G.ACCELERATION) */
  printSchedule(patches, cout_doing, "MPM::scheduleExMomIntegrated");
  contactModel->addComputesAndRequires(sched,
                                       patches,
                                       matls,
                                       d_mpmLabels->gVelocityStarLabel);
}

/*!----------------------------------------------------------------------
 * scheduleSetGridBoundaryConditions
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)

{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleSetGridBoundaryConditions");
  Task* t = scinew Task("MPM::setGridBoundaryConditions",
                        this,
                        &SerialMPM::setGridBoundaryConditions);

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  t->modifies(d_mpmLabels->gAccelerationLabel, mss);
  t->modifies(d_mpmLabels->gVelocityStarLabel, mss);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  if (!d_mpmFlags->d_doGridReset) {
    t->requires(Task::OldDW, d_mpmLabels->gDisplacementLabel, Ghost::None);
    t->computes(d_mpmLabels->gDisplacementLabel);
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
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    std::string interp_type = d_mpmFlags->d_interpolatorType;
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gAcceleration;
      constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(gAcceleration,
                            d_mpmLabels->gAccelerationLabel,
                            matID,
                            patch);
      new_dw->getModifiable(gVelocity_star,
                            d_mpmLabels->gVelocityStarLabel,
                            matID,
                            patch);
      new_dw->get(gVelocity,
                  d_mpmLabels->gVelocityLabel,
                  matID,
                  patch,
                  Ghost::None,
                  0);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch,
                              matID,
                              "Velocity",
                              gVelocity_star,
                              interp_type);
      bc.setBoundaryCondition(patch,
                              matID,
                              "Symmetric",
                              gVelocity_star,
                              interp_type);

      // Now recompute acceleration as the difference between the velocity
      // interpolated to the grid (no bcs applied) and the new velocity_star
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c      = *iter;
        gAcceleration[c] = (gVelocity_star[c] - gVelocity[c]) / delT;
      }

      if (!d_mpmFlags->d_doGridReset) {
        NCVariable<Vector> displacement;
        constNCVariable<Vector> displacementOld;
        new_dw->allocateAndPut(displacement,
                               d_mpmLabels->gDisplacementLabel,
                               matID,
                               patch);
        old_dw->get(displacementOld,
                    d_mpmLabels->gDisplacementLabel,
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpmFlags->d_prescribeDeformation) {
    printSchedule(patches, cout_doing, "MPM::scheduleSetPrescribedMotion");

    Task* t = scinew Task("MPM::setPrescribedMotion",
                          this,
                          &SerialMPM::setPrescribedMotion);

    const MaterialSubset* mss = matls->getUnion();
    t->modifies(d_mpmLabels->gAccelerationLabel, mss);
    t->modifies(d_mpmLabels->gVelocityStarLabel, mss);
    t->requires(Task::OldDW, d_materialManager->get_delt_label());
    if (!d_mpmFlags->d_doGridReset) {
      t->requires(Task::OldDW, d_mpmLabels->gDisplacementLabel, Ghost::None);
      t->modifies(d_mpmLabels->gDisplacementLabel, mss);
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
    double time = d_materialManager->getElapsedTime();
    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matlID            = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity_star, gAcceleration;

      new_dw->getModifiable(gVelocity_star,
                            d_mpmLabels->gVelocityStarLabel,
                            matlID,
                            patch);
      new_dw->getModifiable(gAcceleration,
                            d_mpmLabels->gAccelerationLabel,
                            matlID,
                            patch);

      gAcceleration.initialize(Vector(0.0));

      // Get F and Q from file by interpolating between available times
      auto t_upper_iter = std::upper_bound(d_prescribedTimes.begin(),
                                           d_prescribedTimes.end(),
                                           time);

      auto t_upper_index = t_upper_iter - d_prescribedTimes.begin();
      if (t_upper_iter == d_prescribedTimes.end()) {
        t_upper_index = --t_upper_iter - d_prescribedTimes.begin();
      }

      auto t_lower = *(t_upper_iter - 1);
      auto t_upper = *t_upper_iter;
      auto ss      = (time - t_lower) / (t_upper - t_lower);

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
      if (d_mpmFlags->d_exactDeformation) {
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
                        d_materialManager->get_delt_label(),
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
            new_dw->put(delt_vartype(tst),
                        d_materialManager->get_delt_label(),
                        getLevel(patches));
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
        if (d_mpmFlags->d_exactDeformation) {
          gVelocity_star[n] = (F_upper * F_lower.Inverse() - Ident) *
                              R_previous_inv * QQ * position / delT;
        } else {
          gVelocity_star[n] =
            Fdotstar * Ft_inv * R_previous_inv * QQ * position;
        }

      } // Node Iterator

      if (!d_mpmFlags->d_doGridReset) {
        NCVariable<Vector> displacement;
        constNCVariable<Vector> displacementOld;
        new_dw->allocateAndPut(displacement,
                               d_mpmLabels->gDisplacementLabel,
                               matlID,
                               patch);
        old_dw->get(displacementOld,
                    d_mpmLabels->gDisplacementLabel,
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeXPICVelocities");

  // Particle velocities
  Task* t_part = scinew Task("MPM::computeParticleVelocityXPIC",
                             this,
                             &SerialMPM::computeParticleVelocityXPIC);

  t_part->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
  t_part->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t_part->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, Ghost::None);

  t_part->requires(Task::NewDW,
                   d_mpmLabels->gVelocityLabel,
                   Ghost::AroundCells,
                   d_numGhostNodes);

  t_part->computes(d_mpmLabels->pVelocityXPICLabel);

  sched->addTask(t_part, patches, matls);

  // Grid velocities
  Task* t_grid = scinew Task("MPM::computeGridVelocityXPIC",
                             this,
                             &SerialMPM::computeGridVelocityXPIC);

  t_grid->requires(Task::OldDW,
                   d_mpmLabels->pXLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->requires(Task::OldDW,
                   d_mpmLabels->pMassLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->requires(Task::OldDW,
                   d_mpmLabels->pSizeLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->requires(Task::OldDW,
                   d_mpmLabels->pDefGradLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);
  t_grid->requires(Task::NewDW,
                   d_mpmLabels->pVelocityXPICLabel,
                   Ghost::AroundNodes,
                   d_numGhostParticles);

  t_grid->requires(Task::NewDW,
                   d_mpmLabels->gMassLabel,
                   Ghost::AroundCells,
                   d_numGhostNodes);

  t_grid->computes(d_mpmLabels->gVelocityXPICLabel);

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

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    auto materials = d_materialManager->mpm_materials();
    for (const auto& material : materials) {

      auto mat_id    = material->getDWIndex();
      auto particles = old_dw->getParticleSubset(mat_id, patch);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpmLabels->pXLabel, particles);

      constParticleVariable<Matrix3> pSize, pDefGrad;
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, particles);
      old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, particles);

      constNCVariable<Vector> gVelocity;
      new_dw->get(gVelocity,
                  d_mpmLabels->gVelocityLabel,
                  mat_id,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      ParticleVariable<Vector> pVelocityXPIC;
      new_dw->allocateAndPut(pVelocityXPIC,
                             d_mpmLabels->pVelocityXPICLabel,
                             particles);

      for (auto particle : *particles) {

        interpolator->findCellAndWeights(pX[particle],
                                         ni,
                                         S,
                                         pSize[particle],
                                         pDefGrad[particle]);
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

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    auto materials = d_materialManager->mpm_materials();
    for (const auto& material : materials) {

      auto mat_id    = material->getDWIndex();
      auto particles = old_dw->getParticleSubset(mat_id, patch);

      constParticleVariable<double> pMass;
      old_dw->get(pMass, d_mpmLabels->pMassLabel, particles);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpmLabels->pXLabel, particles);

      constParticleVariable<Matrix3> pSize, pDefGrad;
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, particles);
      old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, particles);

      constParticleVariable<Vector> pVelocityXPIC;
      new_dw->get(pVelocityXPIC, d_mpmLabels->pVelocityXPICLabel, particles);

      constNCVariable<double> gMass;
      new_dw->get(gMass,
                  d_mpmLabels->gMassLabel,
                  mat_id,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);

      NCVariable<Vector> gVelocityXPIC;
      new_dw->allocateAndPut(gVelocityXPIC,
                             d_mpmLabels->gVelocityXPICLabel,
                             mat_id,
                             patch);
      gVelocityXPIC.initialize(Vector(0.0, 0.0, 0.0));

      IntVector node;
      for (auto particle : *particles) {

        interpolator->findCellAndWeights(pX[particle],
                                         ni,
                                         S,
                                         pSize[particle],
                                         pDefGrad[particle]);
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

/*!----------------------------------------------------------------------
 * scheduleComputeDeformationGradient
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeDeformationGradient(SchedulerP& sched,
                                              const PatchSet* patches,
                                              const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
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
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    // Add requires and computes for vel grad/def grad
    d_defGradComputer->addComputesAndRequires(t, mpm_matl, patches);
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
  printTask(patches,
            patches->get(0),
            cout_doing,
            "Doing computeDeformationGradient");

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
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    int matID = mpm_matl->getDWIndex();
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Matrix3> pVelGrad_mid, pDefGrad_mid;
      new_dw->get(pVelGrad_mid,  d_mpmLabels->pVelGradLabel_preReloc, pset);
      new_dw->get(pDefGrad_mid,  d_mpmLabels->pDefGradLabel_preReloc, pset);
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  scheduleUnrotateStressAndDeformationRate(sched, patches, matls);

  /* Create a task for computing the stress tensor */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeStressTensor");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeStressTensor",
                        this,
                        &SerialMPM::computeStressTensor);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    // Add requires and computes for constitutive model
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    t->computes(d_mpmLabels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_materialManager->get_delt_label(), getLevel(patches));

  if (d_mpmFlags->d_reductionVars->accStrainEnergy ||
      d_mpmFlags->d_reductionVars->strainEnergy) {
    t->computes(d_mpmLabels->StrainEnergyLabel);
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleUnrotateStressAndDeformationRate");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeUnrotatedStressAndDeformationRate",
                        this,
                        &SerialMPM::computeUnrotatedStressAndDeformationRate);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->requires(Task::OldDW,
                  d_mpmLabels->pParticleIDLabel,
                  matlset,
                  Ghost::None);
      t->requires(Task::OldDW,
                  d_mpmLabels->pPolarDecompRLabel,
                  matlset,
                  Ghost::None);
      t->requires(Task::NewDW,
                  d_mpmLabels->pPolarDecompRMidLabel,
                  matlset,
                  Ghost::None);
      t->requires(Task::OldDW, d_mpmLabels->pStressLabel, matlset, Ghost::None);
      t->requires(Task::NewDW,
                  d_mpmLabels->pVelGradLabel_preReloc,
                  matlset,
                  Ghost::None);

      t->computes(d_mpmLabels->pDeformRateMidLabel, matlset);
      t->computes(d_mpmLabels->pStressUnrotatedLabel, matlset);
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
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch   = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pStress_old, pR_old, pR_mid,
          pVelGrad_mid;
        old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
        old_dw->get(pR_old, d_mpmLabels->pPolarDecompRLabel, pset);
        new_dw->get(pR_mid, d_mpmLabels->pPolarDecompRMidLabel, pset);
        old_dw->get(pStress_old, d_mpmLabels->pStressLabel, pset);
        new_dw->get(pVelGrad_mid, d_mpmLabels->pVelGradLabel_preReloc, pset);

        ParticleVariable<Matrix3> pDeformRate_mid, pStress_old_unrotated;
        new_dw->allocateAndPut(pDeformRate_mid,
                               d_mpmLabels->pDeformRateMidLabel,
                               pset);
        new_dw->allocateAndPut(pStress_old_unrotated,
                               d_mpmLabels->pStressUnrotatedLabel,
                               pset);

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

  for (int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

    if (cout_dbg.active()) {
      cout_dbg << " Patch = " << (patches->get(0))->getID();
      cout_dbg << " Mat = " << m;
    }

    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleRotateStress");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeRotatedStress",
                        this,
                        &SerialMPM::computeRotatedStress);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->requires(Task::OldDW,
                  d_mpmLabels->pParticleIDLabel,
                  matlset,
                  Ghost::None);
      t->requires(Task::NewDW,
                  d_mpmLabels->pPolarDecompRLabel_preReloc,
                  matlset,
                  Ghost::None);

      t->modifies(d_mpmLabels->pStressLabel_preReloc, matlset);
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
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch   = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pR_new;
        ParticleVariable<Matrix3> pStress_new;
        old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
        new_dw->get(pR_new, d_mpmLabels->pPolarDecompRLabel_preReloc, pset);
        new_dw->getModifiable(pStress_new,
                              d_mpmLabels->pStressLabel_preReloc,
                              pset);

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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  /* Create a task for computing the damage variables */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeBasicDamage");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::computeBasicDamage",
                        this,
                        &SerialMPM::computeBasicDamage);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    // Add requires and computes for vel grad/def grad
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* d_basicDamageModel =
        mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addComputesAndRequires(t,
                                                 mpm_matl,
                                                 patches,
                                                 d_mpmLabels);
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

  for (int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

    if (cout_dbg.active()) {
      cout_dbg << " Patch = " << (patches->get(0))->getID() << " Mat = " << m;
    }

    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    if (cout_dbg.active()) {
      cout_dbg << " MPM_Mat = " << mpm_matl;
    }

    // Compute basic damage
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->computeBasicDamage(patches,
                                           mpm_matl,
                                           old_dw,
                                           new_dw,
                                           lb);
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleUpdateErosionParameter");

  Task* t      = scinew Task("MPM::updateErosionParameter",
                        this,
                        &SerialMPM::updateErosionParameter);
  int numMatls = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* d_basicDamageModel =
        mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addRequiresLocalizationParameter(t,
                                                           mpm_matl,
                                                           patches);
    }
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addRequiresDamageParameter(t, mpm_matl, patches);
  }
  /*
  t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, Ghost::None);
  */
  t->computes(d_mpmLabels->pLocalizedMPMLabel);

  if (d_mpmFlags->d_deleteRogueParticles) {
    t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
    t->computes(d_mpmLabels->numLocInCellLabel);
    t->computes(d_mpmLabels->numInCellLabel);
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

      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: mpm_matl* = " << mpm_matl
                 << " matID = " << matID << " pset* = " << pset << "\n";
      }

      // Get the localization info
      ParticleVariable<int> isLocalized;
      new_dw->allocateAndPut(isLocalized,
                             d_mpmLabels->pLocalizedMPMLabel,
                             pset);
      for (auto particle : *pset) {
        isLocalized[particle] = 0;
      }

      // Update the localization info from basic damage model
      if (mpm_matl->doBasicDamage()) {
        Vaango::BasicDamageModel* basicDamageModel =
          mpm_matl->getBasicDamageModel();
        basicDamageModel->getLocalizationParameter(patch,
                                                   isLocalized,
                                                   matID,
                                                   old_dw,
                                                   new_dw);
      }

      // Update the localization info from constitutive model
      mpm_matl->getConstitutiveModel()->getDamageParameter(patch,
                                                           isLocalized,
                                                           matID,
                                                           old_dw,
                                                           new_dw);

      /*
      constParticleVariable<long64> pParticleID;;
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
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

      if (d_mpmFlags->d_deleteRogueParticles) {
        // The following looks for localized particles that are isolated
        // either individually or in small groups
        // Ghost::GhostType  gac = Ghost::AroundCells;
        CCVariable<int> numLocInCell, numInCell;
        new_dw->allocateAndPut(numLocInCell,
                               d_mpmLabels->numLocInCellLabel,
                               matID,
                               patch);
        new_dw->allocateAndPut(numInCell,
                               d_mpmLabels->numInCellLabel,
                               matID,
                               patch);
        numLocInCell.initialize(0);
        numInCell.initialize(0);

        constParticleVariable<Point> pX;
        old_dw->get(pX, d_mpmLabels->pXLabel, pset);

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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpmFlags->d_deleteRogueParticles) {
    printSchedule(patches, cout_doing, "MPM::scheduleFindRogueParticles");

    Task* t              = scinew Task("MPM::findRogueParticles",
                          this,
                          &SerialMPM::findRogueParticles);
    Ghost::GhostType gac = Ghost::AroundCells;
    t->requires(Task::NewDW, d_mpmLabels->numLocInCellLabel, gac, 1);
    t->requires(Task::NewDW, d_mpmLabels->numInCellLabel, gac, 1);
    t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);
    t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, Ghost::None);
    t->modifies(d_mpmLabels->pLocalizedMPMLabel);

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

      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      // The following looks for localized particles that are isolated
      // either individually or in small groups
      Ghost::GhostType gac = Ghost::AroundCells;
      constCCVariable<int> numLocInCell, numInCell;
      constParticleVariable<Point> pX;

      ParticleVariable<int> isLocalized;
      constParticleVariable<long64> pParticleID;

      new_dw->get(numLocInCell,
                  d_mpmLabels->numLocInCellLabel,
                  matID,
                  patch,
                  gac,
                  1);
      new_dw->get(numInCell, d_mpmLabels->numInCellLabel, matID, patch, gac, 1);
      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
      new_dw->getModifiable(isLocalized, d_mpmLabels->pLocalizedMPMLabel, pset);

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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeAccStrainEnergy");

  Task* t = scinew Task("MPM::computeAccStrainEnergy",
                        this,
                        &SerialMPM::computeAccStrainEnergy);
  t->requires(Task::OldDW, d_mpmLabels->AccStrainEnergyLabel);
  t->requires(Task::NewDW, d_mpmLabels->StrainEnergyLabel);
  t->computes(d_mpmLabels->AccStrainEnergyLabel);
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
  old_dw->get(accStrainEnergy, d_mpmLabels->AccStrainEnergyLabel);

  // Get the incremental strain energy from the new datawarehouse
  sum_vartype incStrainEnergy;
  new_dw->get(incStrainEnergy, d_mpmLabels->StrainEnergyLabel);

  // Add the two a put into new dw
  double totalStrainEnergy = (double)accStrainEnergy + (double)incStrainEnergy;
  new_dw->put(max_vartype(totalStrainEnergy),
              d_mpmLabels->AccStrainEnergyLabel);
}

/*!----------------------------------------------------------------------
 * scheduleComputeHeatExchange
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeHeatExchange(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  /* computeHeatExchange
   *   in(G.MASS, G.TEMPERATURE, G.EXTERNAL_HEAT_RATE)
   *   operation(peform heat exchange which will cause each of
   *   velocity fields to exchange heat according to
   *   the temperature differences)
   *   out(G.EXTERNAL_HEAT_RATE) */

  printSchedule(patches, cout_doing, "MPM::scheduleComputeHeatExchange");

  Task* t = scinew Task("ThermalContact::computeHeatExchange",
                        thermalContactModel,
                        &ThermalContact::computeHeatExchange);

  thermalContactModel->addComputesAndRequires(t, patches, matls);
  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleComputeInternalHeatRate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeInternalHeatRate(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeInternalHeatRate");
  heatConductionModel->scheduleComputeInternalHeatRate(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleComputeNodalHeatFlux
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleComputeNodalHeatFlux(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleComputeNodalHeatFlux");
  heatConductionModel->scheduleComputeNodalHeatFlux(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleSolveHeatEquations
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleSolveHeatEquations(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleSolveHeatEquations");
  heatConductionModel->scheduleSolveHeatEquations(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleIntegrateTemperatureRate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleIntegrateTemperatureRate(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }
  printSchedule(patches, cout_doing, "MPM::scheduleIntegrateTemperatureRate");
  heatConductionModel->scheduleIntegrateTemperatureRate(sched, patches, matls);
}

/*!----------------------------------------------------------------------
 * scheduleAddNewParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleAddNewParticles(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  // if  manual_new_material==false, DON't do this task OR
  // if  create_new_particles==true, DON'T do this task
  if (!d_mpmFlags->d_addNewMaterial || d_mpmFlags->d_createNewParticles) {
    return;
  }

  // if  manual__new_material==true, DO this task OR
  // if  create_new_particles==false, DO this task

  printSchedule(patches, cout_doing, "MPM::scheduleAddNewParticles");
  Task* t =
    scinew Task("MPM::addNewParticles", this, &SerialMPM::addNewParticles);

  int numMatls = d_materialManager->getNumMaterials("MPM");

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    mpm_matl->getParticleCreator()->allocateVariablesAddRequires(t,
                                                                 mpm_matl,
                                                                 patches);

    // Deformation gradient related stuff
    d_defGradComputer->addRequiresForConvert(t, mpm_matl);

    // Constitutive model related stuff
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->allocateCMDataAddRequires(t, mpm_matl, patches, d_mpmLabels);

    // Basic damage model related stuff
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->allocateDamageDataAddRequires(t,
                                                      mpm_matl,
                                                      patches,
                                                      d_mpmLabels);
      basicDamageModel->addRequiresLocalizationParameter(t, mpm_matl, patches);
    } else {
      cm->addRequiresDamageParameter(t, mpm_matl, patches);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * addNewParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::addNewParticles(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing addNewParticles");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    // Find the mpm material that the void particles are going to change
    // into.
    MPMMaterial* null_matl = 0;
    int null_matID         = -1;
    for (int void_matl = 0; void_matl < numMPMMatls; void_matl++) {
      null_matID =
        d_materialManager->getMaterial("MPM", void_matl)->nullGeomObject();

      if (cout_dbg.active()) {
        cout_dbg << "Null DWI = " << null_matID << "\n";
      }

      if (null_matID != -1) {
        null_matl  = d_materialManager->getMaterial("MPM", void_matl);
        null_matID = null_matl->getDWIndex();
        break;
      }
    }
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      if (matID == null_matID) {
        continue;
      }

      ParticleVariable<int> damage;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      new_dw->allocateTemporary(damage, pset);
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        damage[*iter] = 0;
      }

      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

      if (mpm_matl->doBasicDamage()) {
        Vaango::BasicDamageModel* basicDamageModel =
          mpm_matl->getBasicDamageModel();
        basicDamageModel->getLocalizationParameter(patch,
                                                   damage,
                                                   matID,
                                                   old_dw,
                                                   new_dw);
      } else {
        mpm_matl->getConstitutiveModel()->getDamageParameter(patch,
                                                             damage,
                                                             matID,
                                                             old_dw,
                                                             new_dw);
      }

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        if (damage[*iter]) {

          if (cout_dbg.active()) {
            cout_dbg << "damage[" << *iter << "]=" << damage[*iter] << "\n";
          }

          delset->addParticle(*iter);
        }
      }

      // Find the mpm material that corresponds to the void particles.
      // Will probably be the same type as the deleted ones, but have
      // different parameters.

      int numparticles = delset->numParticles();

      if (cout_dbg.active()) {
        cout_dbg << "Num Failed Particles = " << numparticles << "\n";
      }

      if (numparticles != 0) {

        if (cout_dbg.active()) {
          cout_dbg << "Deleted " << numparticles << " particles"
                   << "\n";
        }

        ParticleCreator* particle_creator = null_matl->getParticleCreator();
        ParticleSubset* addset =
          scinew ParticleSubset(numparticles, null_matID, patch);

        if (cout_dbg.active()) {
          cout_dbg << "Address of delset = " << delset << "\n";
          cout_dbg << "Address of pset = " << pset << "\n";
          cout_dbg << "Address of addset = " << addset << "\n";
        }

        ParticleLabelVariableMap* newState = scinew ParticleLabelVariableMap;

        if (cout_dbg.active()) {
          cout_dbg << "Address of newState = " << newState << "\n";
          cout_dbg << "Null Material"
                   << "\n";
        }

        // std::vector<const VarLabel* > particle_labels =
        //   particle_creator->returnParticleState();

        // printParticleLabels(particle_labels, old_dw, null_matID, patch);

        if (cout_dbg.active()) {
          cout_dbg << "MPM Material"
                   << "\n";
        }

        // std::vector<const VarLabel* > mpm_particle_labels =
        //   mpm_matl->getParticleCreator()->returnParticleState();
        // printParticleLabels(mpm_particle_labels, old_dw, matID, patch);

        particle_creator->allocateVariablesAdd(new_dw,
                                               addset,
                                               newState,
                                               delset,
                                               old_dw);

        // Add null-matl deformation gradient etc.
        d_defGradComputer->copyAndDeleteForConvert(new_dw,
                                                   addset,
                                                   newState,
                                                   delset,
                                                   old_dw);

        // Need to do the constitutive models particle variables;
        null_matl->getConstitutiveModel()->allocateCMDataAdd(new_dw,
                                                             addset,
                                                             newState,
                                                             delset,
                                                             old_dw);

        // Add null material basic damage variables
        if (null_matl->doBasicDamage()) {
          null_matl->getBasicDamageModel()
            ->copyDamageDataFromDeletedToAddedParticle(new_dw,
                                                       addset,
                                                       newState,
                                                       delset,
                                                       old_dw);
        }

        // Need to carry forward the cellNAPID for each time step;
        // Move the particle variable declarations in ParticleCreator.h to one
        // of the functions to save on memory;

        if (cout_dbg.active()) {
          cout_dbg << "addset num particles = " << addset->numParticles()
                   << " for material " << addset->getMatlIndex() << "\n";
        }

        new_dw->addParticles(patch, null_matID, newState);

        if (cout_dbg.active()) {
          cout_dbg << "Calling deleteParticles for material: " << matID << "\n";
        }

        new_dw->deleteParticles(delset);

      } else {
        delete delset;
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleConvertLocalizedParticles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleConvertLocalizedParticles(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  // if  create_new_particles==false, DON't do this task OR
  // if  manual_create_new_matl==true, DON'T do this task
  if (!d_mpmFlags->d_createNewParticles || d_mpmFlags->d_addNewMaterial) {
    return;
  }

  // if  create_new_particles==true, DO this task OR
  // if  manual_create_new_matl==false, DO this task
  printSchedule(patches, cout_doing, "MPM::scheduleConvertLocalizedParticles");
  Task* t = scinew Task("MPM::convertLocalizedParticles",
                        this,
                        &SerialMPM::convertLocalizedParticles);

  int numMatls = d_materialManager->getNumMaterials("MPM");

  if (cout_convert.active()) {
    cout_convert << "MPM:scheduleConvertLocalizedParticles : numMatls = "
                 << numMatls << "\n";
  }

  for (int m = 0; m < numMatls; m += 2) {

    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    if (cout_convert.active()) {
      cout_convert << " Material = " << m << " mpm_matl = " << mpm_matl << "\n";
    }

    mpm_matl->getParticleCreator()->allocateVariablesAddRequires(t,
                                                                 mpm_matl,
                                                                 patches);
    if (cout_convert.active()) {
      cout_convert << "   Done ParticleCreator::allocateVariablesAddRequires\n";
    }

    // Deformation gradient related stuff
    d_defGradComputer->addRequiresForConvert(t, mpm_matl);

    // Constitutive model related stuff
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cout_convert.active()) {
      cout_convert << "   cm = " << cm << "\n";
    }

    cm->allocateCMDataAddRequires(t, mpm_matl, patches, d_mpmLabels);

    if (cout_convert.active()) {
      cout_convert << "   Done cm->allocateCMDataAddRequires = "
                   << "\n";
    }

    // Basic damage model related stuff
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->allocateDamageDataAddRequires(t,
                                                      mpm_matl,
                                                      patches,
                                                      d_mpmLabels);
      basicDamageModel->addRequiresLocalizationParameter(t, mpm_matl, patches);
    } else {
      cm->addRequiresDamageParameter(t, mpm_matl, patches);
      if (cout_convert.active()) {
        cout_convert << "   Done cm->addRequiresDamageParameter = "
                     << "\n";
      }
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * convertLocalizedParticles
 *   Convert the localized particles of material "i" into particles of
 *   material "i+1"
 *-----------------------------------------------------------------------*/
void
SerialMPM::convertLocalizedParticles(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  // This function is called only when the "createNewParticles" flag is on.
  // When this flag is on, every second material is a copy of the previous
  // material and is used the material into which particles of the previous
  // material are converted.
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing convertLocalizedParticles");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    if (cout_convert.active()) {
      cout_convert << "MPM::convertLocalizeParticles:: on patch"
                   << patch->getID() << " numMPMMaterials = " << numMPMMatls
                   << "\n";
    }

    for (int m = 0; m < numMPMMatls; m += 2) {

      if (cout_convert.active()) {
        cout_convert << " material # = " << m << "\n";
      }

      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      if (cout_convert.active()) {
        cout_convert << " mpm_matl* = " << mpm_matl << " matID = " << matID
                     << " pset* = " << pset << "\n";
      }

      ParticleVariable<int> isLocalized;

      // old_dw->allocateTemporary(isLocalized, pset);
      new_dw->allocateTemporary(isLocalized, pset);

      ParticleSubset::iterator iter = pset->begin();
      for (; iter != pset->end(); iter++) {
        isLocalized[*iter] = 0;
      }

      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

      if (mpm_matl->doBasicDamage()) {
        Vaango::BasicDamageModel* basicDamageModel =
          mpm_matl->getBasicDamageModel();
        basicDamageModel->getLocalizationParameter(patch,
                                                   isLocalized,
                                                   matID,
                                                   old_dw,
                                                   new_dw);
      } else {
        mpm_matl->getConstitutiveModel()->getDamageParameter(patch,
                                                             isLocalized,
                                                             matID,
                                                             old_dw,
                                                             new_dw);
      }

      if (cout_convert.active()) {
        cout_convert << " Got Damage Parameter"
                     << "\n";
      }

      iter = pset->begin();
      for (; iter != pset->end(); iter++) {
        if (isLocalized[*iter]) {

          if (cout_convert.active()) {
            cout_convert << "damage[" << *iter << "]=" << isLocalized[*iter]
                         << "\n";
          }
          delset->addParticle(*iter);
        }
      }

      if (cout_convert.active()) {
        cout_convert << " Created Delset ";
      }

      int numparticles = delset->numParticles();

      if (cout_convert.active()) {
        cout_convert << " numparticles = " << numparticles << "\n";
      }

      if (numparticles != 0) {

        if (cout_convert.active()) {
          cout_convert << " Converting " << numparticles
                       << " particles of material " << m
                       << " into particles of material " << (m + 1)
                       << " in patch " << p << "\n";
        }

        MPMMaterial* conv_matl = d_materialManager->getMaterial("MPM", m + 1);
        int conv_matID         = conv_matl->getDWIndex();

        ParticleCreator* particle_creator = conv_matl->getParticleCreator();
        ParticleSubset* addset =
          scinew ParticleSubset(numparticles, conv_matID, patch);

        ParticleLabelVariableMap* newState = scinew ParticleLabelVariableMap;

        if (cout_convert.active()) {
          cout_convert << "New Material"
                       << "\n";
        }

        // std::vector<const VarLabel* > particle_labels =
        //   particle_creator->returnParticleState();
        // printParticleLabels(particle_labels, old_dw, conv_matID, patch);

        if (cout_convert.active()) {
          cout_convert << "MPM Material"
                       << "\n";
        }

        // std::vector<const VarLabel* > mpm_particle_labels =
        //   mpm_matl->getParticleCreator()->returnParticleState();
        // printParticleLabels(mpm_particle_labels, old_dw, matID, patch);

        particle_creator->allocateVariablesAdd(new_dw,
                                               addset,
                                               newState,
                                               delset,
                                               old_dw);

        // Copy gradient data
        d_defGradComputer->copyAndDeleteForConvert(new_dw,
                                                   addset,
                                                   newState,
                                                   delset,
                                                   old_dw);

        // Copy constitutive model data
        conv_matl->getConstitutiveModel()->allocateCMDataAdd(new_dw,
                                                             addset,
                                                             newState,
                                                             delset,
                                                             old_dw);

        // Add conv material basic damage variables
        if (conv_matl->doBasicDamage()) {
          conv_matl->getBasicDamageModel()
            ->copyDamageDataFromDeletedToAddedParticle(new_dw,
                                                       addset,
                                                       newState,
                                                       delset,
                                                       old_dw);
        }

        if (cout_convert.active()) {
          cout_convert << "addset num particles = " << addset->numParticles()
                       << " for material " << addset->getMatlIndex() << "\n";
        }

        new_dw->addParticles(patch, conv_matID, newState);
        new_dw->deleteParticles(delset);

        // delete addset;
      } else {
        delete delset;
      }
    }

    if (cout_convert.active()) {
      cout_convert << "Done convertLocalizedParticles on patch "
                   << patch->getID() << "\t MPM"
                   << "\n";
    }
  }

  if (cout_convert.active()) {
    cout_convert << "Completed convertLocalizedParticles "
                 << "\n";
  }
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
  for (std::vector<const VarLabel*>::const_iterator it = labels.begin();
       it != labels.end();
       it++) {
    if (dw->exists(*it, matID, patch)) {
      std::cout << (*it)->getName() << " does exists"
                << "\n";
    } else {
      std::cout << (*it)->getName() << " does NOT exists"
                << "\n";
    }
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleInterpolateToParticlesAndUpdate");

  Task* t = scinew Task("MPM::interpolateToParticlesAndUpdate",
                        this,
                        &SerialMPM::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW,
              d_mpmLabels->gAccelerationLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityStarLabel,
              gac,
              d_numGhostNodes);
#ifdef XPIC2_UPDATE
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityXPICLabel,
              gac,
              d_numGhostNodes);
#endif
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureRateLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->frictionalWorkLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pTemperatureLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pVelocityLabel, gnone);
#ifdef XPIC2_UPDATE
  t->requires(Task::OldDW, d_mpmLabels->pVelocityXPICLabel, gnone);
#endif
  t->requires(Task::OldDW, d_mpmLabels->pDispLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pdTdtLabel_preReloc, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pLocalizedMPMLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel_preReloc, gnone);
  t->modifies(d_mpmLabels->pVolumeLabel_preReloc);

  if (d_mpmFlags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_mpmLabels->pLoadCurveIDLabel, Ghost::None);
  }

  if (d_mpmFlags->d_withICE) {
    t->requires(Task::NewDW, d_mpmLabels->dTdt_NCLabel, gac, d_numGhostNodes);
    t->requires(Task::NewDW,
                d_mpmLabels->massBurnFractionLabel,
                gac,
                d_numGhostNodes);
  }

  t->computes(d_mpmLabels->pDispLabel_preReloc);
  t->computes(d_mpmLabels->pVelocityLabel_preReloc);
  t->computes(d_mpmLabels->pAccelerationLabel_preReloc);
  t->computes(d_mpmLabels->pXLabel_preReloc);
  t->computes(d_mpmLabels->pParticleIDLabel_preReloc);
  t->computes(d_mpmLabels->pTemperatureLabel_preReloc);
  t->computes(d_mpmLabels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpmLabels->pMassLabel_preReloc);
  t->computes(d_mpmLabels->pSizeLabel_preReloc);
  t->computes(d_mpmLabels->pXXLabel);

  //__________________________________
  //  reduction variables
  if (d_mpmFlags->d_reductionVars->momentum) {
    t->computes(d_mpmLabels->TotalMomentumLabel);
  }
  if (d_mpmFlags->d_reductionVars->KE) {
    t->computes(d_mpmLabels->KineticEnergyLabel);
  }
  if (d_mpmFlags->d_reductionVars->thermalEnergy) {
    t->computes(d_mpmLabels->ThermalEnergyLabel);
  }
  if (d_mpmFlags->d_reductionVars->centerOfMass) {
    t->computes(d_mpmLabels->CenterOfMassPositionLabel);
  }
  if (d_mpmFlags->d_reductionVars->mass) {
    t->computes(d_mpmLabels->TotalMassLabel);
  }
  if (d_mpmFlags->d_reductionVars->volDeformed) {
    t->computes(d_mpmLabels->TotalVolumeDeformedLabel);
  }

  // debugging scalar
  if (d_mpmFlags->d_withColor) {
    t->requires(Task::OldDW, d_mpmLabels->pColorLabel, Ghost::None);
    t->computes(d_mpmLabels->pColorLabel_preReloc);
  }

  // Carry Forward particle refinement flag
  if (d_mpmFlags->d_refineParticles) {
    t->requires(Task::OldDW, d_mpmLabels->pRefinedLabel, Ghost::None);
    t->computes(d_mpmLabels->pRefinedLabel_preReloc);
  }

  // Carry forward external heat flux for switch from explicit to implicit
  t->requires(Task::OldDW, d_mpmLabels->pExternalHeatFluxLabel, Ghost::None);
  t->computes(d_mpmLabels->pExternalHeatFluxLabel_preReloc);

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->requires(Task::OldDW, d_mpmLabels->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(d_mpmLabels->NC_CCweightLabel, z_matl);

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
    printTask(patches,
              patch,
              cout_doing,
              "Doing interpolateToParticlesAndUpdate");

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
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
    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));
    // bool combustion_problem=false;

    // Material* reactant;
    // int RMI = -99;
    // reactant = d_materialManager->getMaterialByName("reactant");
    // if(reactant != 0){
    //   RMI = reactant->getDWIndex();
    //   combustion_problem=true;
    // }
    double move_particles = 1.;
    if (!d_mpmFlags->d_doGridReset) {
      move_particles = 0.;
    }

    // Copy NC_CCweight (only material 0)
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_new;
    old_dw->get(NC_CCweight, d_mpmLabels->NC_CCweightLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NC_CCweight_new,
                           d_mpmLabels->NC_CCweightLabel,
                           0,
                           patch);
    NC_CCweight_new.copyData(NC_CCweight);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      // Copy particle IDs and particle size
      constParticleVariable<long64> pParticleID;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<long64> pParticleID_new;
      ParticleVariable<Matrix3> pSize_new;
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pParticleID_new,
                             d_mpmLabels->pParticleIDLabel_preReloc,
                             pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      new_dw->allocateAndPut(pSize_new, d_mpmLabels->pSizeLabel_preReloc, pset);
      pParticleID_new.copyData(pParticleID);
      pSize_new.copyData(pSize);

      // Copy needed for switch from explicit to implicit MPM
      constParticleVariable<double> pExtHeatFlux;
      ParticleVariable<double> pExtHeatFlux_new;
      old_dw->get(pExtHeatFlux, d_mpmLabels->pExternalHeatFluxLabel, pset);
      new_dw->allocateAndPut(pExtHeatFlux_new,
                             d_mpmLabels->pExternalHeatFluxLabel_preReloc,
                             pset);
      pExtHeatFlux_new.copyData(pExtHeatFlux);

      // Get particle variables
      constParticleVariable<int> pLocalized_new;
      new_dw->get(pLocalized_new, d_mpmLabels->pLocalizedMPMLabel, pset);

      constParticleVariable<double> pMass, pTemperature, pdTdt_new;
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      new_dw->get(pdTdt_new, d_mpmLabels->pdTdtLabel_preReloc, pset);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_mpmLabels->pXLabel, pset);

      constParticleVariable<Vector> pVelocity, pDisp;
      old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);

#ifdef XPIC2_UPDATE
      constParticleVariable<Vector> pVelocityXPIC;
      new_dw->get(pVelocityXPIC, d_mpmLabels->pVelocityXPICLabel, pset);
#endif

      constParticleVariable<Matrix3> pDefGrad_new;
      new_dw->get(pDefGrad_new, d_mpmLabels->pDefGradLabel_preReloc, pset);

      // Allocate updated particle variables
      ParticleVariable<double> pMass_new, pVolume_new, pTemp_new, pTempPrev_new;
      new_dw->allocateAndPut(pMass_new, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->getModifiable(pVolume_new,
                            d_mpmLabels->pVolumeLabel_preReloc,
                            pset);
      new_dw->allocateAndPut(pTemp_new,
                             d_mpmLabels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempPrev_new,
                             d_mpmLabels->pTempPreviousLabel_preReloc,
                             pset);

      ParticleVariable<Point> pX_new, pXx;
      new_dw->allocateAndPut(pX_new, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pXx, d_mpmLabels->pXXLabel, pset);

      ParticleVariable<Vector> pVelocity_new, pDisp_new, pAcc_new;
      new_dw->allocateAndPut(pVelocity_new,
                             d_mpmLabels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pDisp_new, d_mpmLabels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(pAcc_new,
                             d_mpmLabels->pAccelerationLabel_preReloc,
                             pset);

      // Get grid variables
      constNCVariable<double> gTemperatureRate, frictionTempRate;
      new_dw->get(gTemperatureRate,
                  d_mpmLabels->gTemperatureRateLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpmLabels->frictionalWorkLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);

      constNCVariable<double> dTdt, massBurnFrac;
      if (d_mpmFlags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpmLabels->dTdt_NCLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpmLabels->massBurnFractionLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create, patch, gac, d_numGhostParticles);
        new_dw->allocateTemporary(massBurnFrac_create,
                                  patch,
                                  gac,
                                  d_numGhostParticles);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);

        dTdt         = dTdt_create;         // reference created data
        massBurnFrac = massBurnFrac_create; // reference created data
      }

      constNCVariable<Vector> gVelocityStar, gAcceleration;
      new_dw->get(gVelocityStar,
                  d_mpmLabels->gVelocityStarLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gAcceleration,
                  d_mpmLabels->gAccelerationLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);

#ifdef XPIC2_UPDATE
      constNCVariable<Vector> gVelocityXPIC;
      new_dw->get(gVelocityXPIC,
                  d_mpmLabels->gVelocityXPICLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
#endif

      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

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

        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad_new[idx]);

        Vector velocity(0.0, 0.0, 0.0);
        Vector acceleration(0.0, 0.0, 0.0);
#ifdef XPIC2_UPDATE
        Vector velocityXPIC(0.0, 0.0, 0.0);
#endif
        double fricTempRate = 0.0;
        double tempRate     = 0.0;
        double burnFraction = 0.0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          velocity += gVelocityStar[node] * S[k];
          acceleration += gAcceleration[node] * S[k];

#ifdef XPIC2_UPDATE
          velocityXPIC += gVelocityXPIC[node] * S[k];
#endif

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

          fricTempRate = frictionTempRate[node] * d_mpmFlags->d_addFrictionWork;
          tempRate +=
            (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
          burnFraction += massBurnFrac[node] * S[k];
        }

// Update the particle's position and velocity
#ifdef XPIC2_UPDATE
        pX_new[idx] = pX[idx] + velocity * delT -
                      0.5 *
                        (acceleration * delT + pVelocity[idx] -
                         2.0 * pVelocityXPIC[idx] + velocityXPIC) *
                        delT;
        pVelocity_new[idx] =
          2.0 * pVelocityXPIC[idx] - velocityXPIC + acceleration * delT;
        pDisp_new[idx] = pDisp[idx] + (pX_new[idx] - pX[idx]);
#else
        pX_new[idx]        = pX[idx] + velocity * delT * move_particles;
        pDisp_new[idx]     = pDisp[idx] + velocity * delT;
        pVelocity_new[idx] = pVelocity[idx] + acceleration * delT;
#endif

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

      // If load curves are being used with VelocityBC then apply
      // these BCs to the boundary particles
      if (d_mpmFlags->d_useLoadCurves) {

        std::vector<VelocityBC*> vbcP;
        bool do_VelocityBCs = false;
        for (auto bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          std::string bcType = bc->getType();
          if (bcType == "Velocity") {
            do_VelocityBCs  = true;
            VelocityBC* vbc = dynamic_cast<VelocityBC*>(bc.get());
            vbcP.push_back(vbc);
          }
        }

        // std::cout << "do_VelocityBCs = " << do_VelocityBCs << "\n";
        if (do_VelocityBCs) {

          // Get the current time
          double time = d_materialManager->getElapsedTime();

          // Get the load curve data
          constParticleVariable<int> pLoadCurveID;
          old_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);

          // Iterate over the particles
          for (auto iter = pset->begin(); iter != pset->end(); iter++) {
            particleIndex idx = *iter;
            int loadCurveID   = pLoadCurveID[idx] - 1;
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

        if ((pMass_new[idx] <= d_mpmFlags->d_minPartMass) ||
            (pTemp_new[idx] < 0.0) || (pLocalized_new[idx] == -999)) {
          if (d_mpmFlags->d_erosionAlgorithm != "none") {
            delset->addParticle(idx);
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

        if (pVelocity_new[idx].length() > d_mpmFlags->d_maxVel) {
          if (d_mpmFlags->d_deleteRogueParticles) {
            if (d_mpmFlags->d_erosionAlgorithm != "none") {
              delset->addParticle(idx);
            }
            std::cout << "\n Warning: particle " << pParticleID[idx]
                      << " hit speed ceiling #1. Deleting particle."
                      << "\n";
          } else {
            if (pVelocity_new[idx].length() >= pVelocity[idx].length()) {
              pVelocity_new[idx] =
                (pVelocity_new[idx] / pVelocity_new[idx].length()) *
                (d_mpmFlags->d_maxVel * .9);
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
    //__________________________________
    //  reduction variables
    if (d_mpmFlags->d_reductionVars->mass) {
      new_dw->put(sum_vartype(totalmass), d_mpmLabels->TotalMassLabel);
    }
    if (d_mpmFlags->d_reductionVars->volDeformed) {
      new_dw->put(sum_vartype(partvoldef),
                  d_mpmLabels->TotalVolumeDeformedLabel);
    }
    if (d_mpmFlags->d_reductionVars->momentum) {
      new_dw->put(sumvec_vartype(totalMom), d_mpmLabels->TotalMomentumLabel);
    }
    if (d_mpmFlags->d_reductionVars->KE) {
      new_dw->put(sum_vartype(ke), d_mpmLabels->KineticEnergyLabel);
    }
    if (d_mpmFlags->d_reductionVars->thermalEnergy) {
      new_dw->put(sum_vartype(thermal_energy), d_mpmLabels->ThermalEnergyLabel);
    }
    if (d_mpmFlags->d_reductionVars->centerOfMass) {
      new_dw->put(sumvec_vartype(CMX), d_mpmLabels->CenterOfMassPositionLabel);
    }

    // std::cout << "Solid mass lost this timestep = " << massLost << "\n";
    // std::cout << "Solid momentum after advection = " << totalMom << "\n";

    // std::cout << "THERMAL ENERGY " << thermal_energy << "\n";

    // delete interpolator;
  }
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateToParticlesAndUpdateMom1 - MOMENTUM FORM : STEP 1
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInterpolateToParticlesAndUpdateMom1(SchedulerP& sched,
                                                       const PatchSet* patches,
                                                       const MaterialSet* matls)

{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleInterpolateToParticlesAndUpdateMom1");

  Task* t = scinew Task("MPM::interpolateToParticlesAndUpdateMom1",
                        this,
                        &SerialMPM::interpolateToParticlesAndUpdateMom1);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW,
              d_mpmLabels->gAccelerationLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityStarLabel,
              gac,
              d_numGhostNodes);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pDispLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pVelocityLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, gnone);

  t->computes(d_mpmLabels->pVelocityLabel_preReloc);
  t->computes(d_mpmLabels->pXLabel_preReloc);
  t->computes(d_mpmLabels->pXXLabel);
  t->computes(d_mpmLabels->pDispLabel_preReloc);

  //__________________________________
  //  reduction variables
  if (d_mpmFlags->d_reductionVars->momentum) {
    t->computes(d_mpmLabels->TotalMomentumLabel);
  }
  if (d_mpmFlags->d_reductionVars->KE) {
    t->computes(d_mpmLabels->KineticEnergyLabel);
  }
  if (d_mpmFlags->d_reductionVars->centerOfMass) {
    t->computes(d_mpmLabels->CenterOfMassPositionLabel);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * interpolateToParticlesAndUpdateMom1 - MOMENTUM FORM : STEP 1
 *-----------------------------------------------------------------------*/
void
SerialMPM::interpolateToParticlesAndUpdateMom1(const ProcessorGroup*,
                                               const PatchSubset* patches,
                                               const MaterialSubset*,
                                               DataWarehouse* old_dw,
                                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing interpolateToParticlesAndUpdateMom1");

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively

    // DON'T MOVE THESE!!!
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke = 0;

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    double move_particles = 1.;
    if (!d_mpmFlags->d_doGridReset) {
      move_particles = 0.;
    }
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> pX;
      ParticleVariable<Point> pX_new, pXx;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Vector> pVelocity_new;
      constParticleVariable<Vector> pDisp;
      ParticleVariable<Vector> pDisp_new;
      constParticleVariable<double> pMass;
      constParticleVariable<Matrix3> pDefGrad_old;
      constParticleVariable<long64> pParticleID;

      // Get the arrays of grid data on which the new part. values depend
      constNCVariable<Vector> gVelocity_star, gAcceleration;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);

      new_dw->allocateAndPut(pVelocity_new,
                             d_mpmLabels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pX_new, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pXx, d_mpmLabels->pXXLabel, pset);
      new_dw->allocateAndPut(pDisp_new, d_mpmLabels->pDispLabel_preReloc, pset);

      Ghost::GhostType gac = Ghost::AroundCells;
      new_dw->get(gVelocity_star,
                  d_mpmLabels->gVelocityStarLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(gAcceleration,
                  d_mpmLabels->gAccelerationLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad_old[idx]);

        Vector vel(0.0, 0.0, 0.0);
        Vector acc(0.0, 0.0, 0.0);

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          vel += gVelocity_star[node] * S[k];
          acc += gAcceleration[node] * S[k];
        }

        // Update the particle's position and velocity
        pX_new[idx]        = pX[idx] + vel * delT * move_particles;
        pDisp_new[idx]     = pDisp[idx] + vel * delT;
        pVelocity_new[idx] = pVelocity[idx] + acc * delT;
        // pXx is only useful if we're not in normal grid resetting mode.
        pXx[idx] = pX[idx] + pDisp_new[idx];

        ke += .5 * pMass[idx] * pVelocity_new[idx].length2();
        CMX = CMX + (pX_new[idx] * pMass[idx]).asVector();
        totalMom += pVelocity_new[idx] * pMass[idx];
      }

      // Delete particles that have left the domain
      // This is only needed if extra cells are being used.
      // Also delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        if (pVelocity_new[idx].length() > d_mpmFlags->d_maxVel) {
          pVelocity_new[idx] =
            (pVelocity_new[idx] / pVelocity_new[idx].length()) *
            d_mpmFlags->d_maxVel;
          std::cout
            << "\n"
            << "Warning: particle " << pParticleID[idx]
            << " hit speed ceiling #2. Modifying particle velocity accordingly."
            << "\n";
          // pVelocity_new[idx]=pVelocity[idx];
        }
      }
    }

    // DON'T MOVE THESE!!!
    //__________________________________
    //  reduction variables
    if (d_mpmFlags->d_reductionVars->momentum) {
      new_dw->put(sumvec_vartype(totalMom), d_mpmLabels->TotalMomentumLabel);
    }
    if (d_mpmFlags->d_reductionVars->KE) {
      new_dw->put(sum_vartype(ke), d_mpmLabels->KineticEnergyLabel);
    }
    if (d_mpmFlags->d_reductionVars->centerOfMass) {
      new_dw->put(sumvec_vartype(CMX), d_mpmLabels->CenterOfMassPositionLabel);
    }

    // delete interpolator;
  }
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateParticleVelToGridMom
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInterpolateParticleVelToGridMom(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleInterpolateParticleVelToGridMom");

  Task* t              = scinew Task("MPM::interpolateParticleVelToGridMom",
                        this,
                        &SerialMPM::interpolateParticleVelToGridMom);
  Ghost::GhostType gan = Ghost::AroundNodes;
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gan, d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pVelocityLabel_preReloc,
              gan,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gan, d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              gan,
              d_numGhostParticles);

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
  t->modifies(d_mpmLabels->gVelocityStarLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * interpolateParticleVelToGridMom
 *-----------------------------------------------------------------------*/
void
SerialMPM::interpolateParticleVelToGridMom(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset*,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              cout_doing,
              "Doing interpolateParticleVelToGridMom");

    int numMatls           = d_materialManager->getNumMaterials("MPM");
    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    std::vector<IntVector> ni(numInfluenceNodes);
    std::vector<double> S(numInfluenceNodes);

    Ghost::GhostType gan   = Ghost::AroundNodes;
    Ghost::GhostType gnone = Ghost::None;
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad_old;

      ParticleSubset* pset = old_dw->getParticleSubset(matID,
                                                       patch,
                                                       gan,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      new_dw->get(pVelocity, d_mpmLabels->pVelocityLabel_preReloc, pset);
      old_dw->get(pDefGrad_old, d_mpmLabels->pDefGradLabel, pset);

      // Create arrays for the grid data
      constNCVariable<double> gMass;
      NCVariable<Vector> gVelocity_star;
      new_dw->get(gMass, d_mpmLabels->gMassLabel, matID, patch, gnone, 0);
      new_dw->getModifiable(gVelocity_star,
                            d_mpmLabels->gVelocityStarLabel,
                            matID,
                            patch);
      gVelocity_star.initialize(Vector(0, 0, 0));

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad_old[idx]);

        Vector pMom = pVelocity[idx] * pMass[idx];

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];
          if (patch->containsNode(node)) {
            gVelocity_star[node] += pMom * S[k];
          }
        }
      } // End of particle loop

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        gVelocity_star[c] /= gMass[c];
      }

      //    setGridBoundaryConditions handles the BCs for gVelocity_star
    } // end of materials loop

    // delete interpolator;

  } // End loop over patches
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateToParticlesAndUpdateMom2 - MOMENTUM FORM : STEP 2
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInterpolateToParticlesAndUpdateMom2(SchedulerP& sched,
                                                       const PatchSet* patches,
                                                       const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "MPM::scheduleInterpolateToParticlesAndUpdate2");

  Task* t = scinew Task("MPM::interpolateToParticlesAndUpdateMom2",
                        this,
                        &SerialMPM::interpolateToParticlesAndUpdateMom2);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureRateLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::NewDW,
              d_mpmLabels->frictionalWorkLabel,
              gac,
              d_numGhostNodes);
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pTemperatureLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel_preReloc, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pdTdtLabel_preReloc, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pLocalizedMPMLabel, gnone);

  if (d_mpmFlags->d_withICE) {
    t->requires(Task::NewDW, d_mpmLabels->dTdt_NCLabel, gac, d_numGhostNodes);
    t->requires(Task::NewDW,
                d_mpmLabels->massBurnFractionLabel,
                gac,
                d_numGhostNodes);
  }

  t->modifies(d_mpmLabels->pVolumeLabel_preReloc);

  t->computes(d_mpmLabels->pParticleIDLabel_preReloc);
  t->computes(d_mpmLabels->pTemperatureLabel_preReloc);
  t->computes(d_mpmLabels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpmLabels->pMassLabel_preReloc);
  t->computes(d_mpmLabels->pSizeLabel_preReloc);

  //__________________________________
  //  reduction variables
  if (d_mpmFlags->d_reductionVars->thermalEnergy) {
    t->computes(d_mpmLabels->ThermalEnergyLabel);
  }
  if (d_mpmFlags->d_reductionVars->mass) {
    t->computes(d_mpmLabels->TotalMassLabel);
  }
  if (d_mpmFlags->d_reductionVars->volDeformed) {
    t->computes(d_mpmLabels->TotalVolumeDeformedLabel);
  }

  // debugging scalar
  if (d_mpmFlags->d_withColor) {
    t->requires(Task::OldDW, d_mpmLabels->pColorLabel, Ghost::None);
    t->computes(d_mpmLabels->pColorLabel_preReloc);
  }

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->requires(Task::OldDW, d_mpmLabels->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(d_mpmLabels->NC_CCweightLabel, z_matl);

  sched->addTask(t, patches, matls);

  // The task will have a reference to z_matl
  if (z_matl->removeReference()) {
    delete z_matl; // shouln't happen, but...
  }
}

/*!----------------------------------------------------------------------
 * interpolateToParticlesAndUpdateMom2 - MOMENTUM FORM : STEP 2
 *-----------------------------------------------------------------------*/
void
SerialMPM::interpolateToParticlesAndUpdateMom2(const ProcessorGroup*,
                                               const PatchSubset* patches,
                                               const MaterialSubset*,
                                               DataWarehouse* old_dw,
                                               DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches,
              patch,
              cout_doing,
              "Doing interpolateToParticlesAndUpdateMom2");

    auto interpolator      = d_mpmFlags->d_interpolator->clone(patch);
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
    int numMPMMatls       = d_materialManager->getNumMaterials("MPM");
    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));
    // bool combustion_problem=false;

    // Material* reactant;
    // int RMI = -99;
    // reactant = d_materialManager->getMaterialByName("reactant");
    // if(reactant != 0){
    // RMI = reactant->getDWIndex();
    // combustion_problem=true;
    //}

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> pX;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pSize_new;
      constParticleVariable<double> pMass, pTemperature, pdTdt;
      ParticleVariable<double> pMass_new, pVolume, pTemp_new;
      constParticleVariable<long64> pParticleID;
      ParticleVariable<long64> pParticleID_new;
      constParticleVariable<int> pLocalized;
      constParticleVariable<Matrix3> pDefGrad_new, pDefGrad_old;

      // for thermal stress analysis
      ParticleVariable<double> pTempPreNew;

      // Get the arrays of grid data on which the new part. values depend
      constNCVariable<double> gTemperatureRate;
      constNCVariable<double> dTdt, massBurnFrac, frictionTempRate;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pParticleID, d_mpmLabels->pParticleIDLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      new_dw->get(pdTdt, d_mpmLabels->pdTdtLabel_preReloc, pset);
      new_dw->get(pDefGrad_new, d_mpmLabels->pDefGradLabel_preReloc, pset);
      new_dw->get(pLocalized, d_mpmLabels->pLocalizedMPMLabel, pset);

      new_dw->getModifiable(pVolume, d_mpmLabels->pVolumeLabel_preReloc, pset);

      new_dw->allocateAndPut(pMass_new, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pParticleID_new,
                             d_mpmLabels->pParticleIDLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTemp_new,
                             d_mpmLabels->pTemperatureLabel_preReloc,
                             pset);

      // for thermal stress analysis
      new_dw->allocateAndPut(pTempPreNew,
                             d_mpmLabels->pTempPreviousLabel_preReloc,
                             pset);

      // Carry forward NC_CCweight
      constNCVariable<double> NC_CCweight;
      NCVariable<double> NC_CCweight_new;
      Ghost::GhostType gnone = Ghost::None;
      old_dw
        ->get(NC_CCweight, d_mpmLabels->NC_CCweightLabel, 0, patch, gnone, 0);
      new_dw->allocateAndPut(NC_CCweight_new,
                             d_mpmLabels->NC_CCweightLabel,
                             0,
                             patch);
      NC_CCweight_new.copyData(NC_CCweight);

      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

      pParticleID_new.copyData(pParticleID);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      new_dw->allocateAndPut(pSize_new, d_mpmLabels->pSizeLabel_preReloc, pset);
      pSize_new.copyData(pSize);

      Ghost::GhostType gac = Ghost::AroundCells;
      new_dw->get(gTemperatureRate,
                  d_mpmLabels->gTemperatureRateLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      new_dw->get(frictionTempRate,
                  d_mpmLabels->frictionalWorkLabel,
                  matID,
                  patch,
                  gac,
                  d_numGhostParticles);
      if (d_mpmFlags->d_withICE) {
        new_dw->get(dTdt,
                    d_mpmLabels->dTdt_NCLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
        new_dw->get(massBurnFrac,
                    d_mpmLabels->massBurnFractionLabel,
                    matID,
                    patch,
                    gac,
                    d_numGhostParticles);
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create, patch, gac, d_numGhostParticles);
        new_dw->allocateTemporary(massBurnFrac_create,
                                  patch,
                                  gac,
                                  d_numGhostParticles);
        dTdt_create.initialize(0.);
        massBurnFrac_create.initialize(0.);
        dTdt         = dTdt_create;         // reference created data
        massBurnFrac = massBurnFrac_create; // reference created data
      }

      double Cp       = mpm_matl->getSpecificHeat();
      double rho_init = mpm_matl->getInitialDensity();

      double rho_frac_min = 0.;
      /*
      if(m == RMI){
        rho_frac_min = .1;
      }
      */

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad_new[idx]);

        Vector vel(0.0, 0.0, 0.0);
        Vector acc(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0.0;
        double burnFraction = 0.0;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < numInfluenceNodes; k++) {
          IntVector node = ni[k];

          fricTempRate = frictionTempRate[node] * d_mpmFlags->d_addFrictionWork;
          tempRate +=
            (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
          burnFraction += massBurnFrac[node] * S[k];
        }

        // Update the particle's position and velocity
        pTemp_new[idx]   = pTemperature[idx] + (tempRate + pdTdt[idx]) * delT;
        pTempPreNew[idx] = pTemperature[idx]; // for thermal stress

        if (cout_heat.active()) {
          cout_heat << "MPM::Particle = " << idx
                    << " T_old = " << pTemperature[idx]
                    << " Tdot = " << tempRate << " dT = " << (tempRate * delT)
                    << " T_new = " << pTemp_new[idx] << "\n";
        }

        double rho;
        if (pVolume[idx] > 0.) {
          rho = std::max(pMass[idx] / pVolume[idx], rho_frac_min * rho_init);
        } else {
          rho = rho_init;
        }
        pMass_new[idx] = Max(pMass[idx] * (1. - burnFraction), 0.);
        pVolume[idx]   = pMass_new[idx] / rho;

        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
      }

      // Delete particles that have left the domain
      // This is only needed if extra cells are being used.
      // Also delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for (auto idx : *pset) {
        if ((pMass_new[idx] <= d_mpmFlags->d_minPartMass) ||
            pTemp_new[idx] < 0. || (pLocalized[idx] == -999)) {
          if (d_mpmFlags->d_erosionAlgorithm != "none") {
            delset->addParticle(idx);
          }
          //        std::cout << "Material = " << m << " Deleted Particle = " <<
          //        idx
          //             << " xold = " << pX[idx] << " xnew = " << pX_new[idx]
          //             << " vold = " << pVelocity[idx] << " vnew = "<<
          //             pVelocity_new[idx]
          //             << " massold = " << pMass[idx] << " massnew = " <<
          //             pMass_new[idx]
          //             << " tempold = " << pTemperature[idx]
          //             << " tempnew = " << pTemp_new[idx]
          //             << " volnew = " << pVolume[idx] << "\n";
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
    //__________________________________
    //  reduction variables
    if (d_mpmFlags->d_reductionVars->mass) {
      new_dw->put(sum_vartype(totalmass), d_mpmLabels->TotalMassLabel);
    }
    if (d_mpmFlags->d_reductionVars->volDeformed) {
      new_dw->put(sum_vartype(partvoldef),
                  d_mpmLabels->TotalVolumeDeformedLabel);
    }
    if (d_mpmFlags->d_reductionVars->thermalEnergy) {
      new_dw->put(sum_vartype(thermal_energy), d_mpmLabels->ThermalEnergyLabel);
    }

    // delete interpolator;
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  if (d_mpmFlags->d_insertParticles) {
    printSchedule(patches, cout_doing, "MPM::scheduleInsertParticles");

    Task* t =
      scinew Task("MPM::insertParticles", this, &SerialMPM::insertParticles);

    t->requires(Task::OldDW, d_materialManager->get_delt_label());

    t->modifies(d_mpmLabels->pXLabel_preReloc);
    t->modifies(d_mpmLabels->pVelocityLabel_preReloc);
    t->requires(Task::OldDW, d_mpmLabels->pColorLabel, Ghost::None);

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
    double time = d_materialManager->getElapsedTime();
    delt_vartype delT;
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    int index = -999;

    for (auto i = 0u; i < d_IPTimes.size(); i++) {
      if (time + delT > d_IPTimes[i] && time <= d_IPTimes[i]) {
        index = i;
      }
    }

    if (index >= 0) {
      int numMPMMatls = d_materialManager->getNumMaterials("MPM");
      for (int m = 0; m < numMPMMatls; m++) {
        MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
        int matID             = mpm_matl->getDWIndex();
        ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

        // Get the arrays of particle values to be changed
        ParticleVariable<Point> pX;
        ParticleVariable<Vector> pVelocity;
        constParticleVariable<double> pcolor;

        old_dw->get(pcolor, d_mpmLabels->pColorLabel, pset);
        new_dw->getModifiable(pX, d_mpmLabels->pXLabel_preReloc, pset);
        new_dw->getModifiable(pVelocity,
                              d_mpmLabels->pVelocityLabel_preReloc,
                              pset);

        int numParticles = pset->end() - pset->begin();
        std::cout << "Insertion: Patch " << p << " now contains "
                  << numParticles << " particles\n";
        // Loop over particles here
        for (auto idx : *pset) {
          if (pcolor[idx] == d_IPColor[index]) {
            pVelocity[idx] = d_IPVelNew[index];
            pX[idx] += d_IPTranslate[index];
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleAddParticles");

  Task* t = scinew Task("MPM::addParticles", this, &SerialMPM::addParticles);

  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();

  t->modifies(d_mpmLabels->pParticleIDLabel_preReloc);
  t->modifies(d_mpmLabels->pXLabel_preReloc);
  t->modifies(d_mpmLabels->pVolumeLabel_preReloc);
  t->modifies(d_mpmLabels->pVelocityLabel_preReloc);
  t->modifies(d_mpmLabels->pMassLabel_preReloc);
  t->modifies(d_mpmLabels->pSizeLabel_preReloc);
  t->modifies(d_mpmLabels->pDispLabel_preReloc);
  t->modifies(d_mpmLabels->pStressLabel_preReloc);
  t->modifies(d_mpmLabels->pColorLabel_preReloc);
  t->modifies(d_mpmLabels->pLocalizedMPMLabel_preReloc);
  t->modifies(d_mpmLabels->pExtForceLabel_preReloc);
  t->modifies(d_mpmLabels->pTemperatureLabel_preReloc);
  t->modifies(d_mpmLabels->pTempPreviousLabel_preReloc);
  t->modifies(d_mpmLabels->pDefGradLabel_preReloc);
  t->modifies(d_mpmLabels->pRefinedLabel_preReloc);
  t->modifies(d_mpmLabels->pVelGradLabel_preReloc);

  t->requires(Task::OldDW,
              d_mpmLabels->pCellNAPIDLabel,
              zeroth_matl,
              Ghost::None);
  t->computes(d_mpmLabels->pCellNAPIDLabel, zeroth_matl);

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
    old_dw->get(NAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NAPID_new, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    NAPID_new.copyData(NAPID);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      ParticleVariable<Point> pX;
      ParticleVariable<Matrix3> pF, pSize, pStress, pvelgrad;
      ParticleVariable<long64> pParticleID;
      ParticleVariable<double> pVolume, pMass, ptemp, ptempP, pcolor;
      ParticleVariable<Vector> pVelocity, pextforce, pDisp;
      ParticleVariable<int> pref, ploc;
      new_dw->getModifiable(pX, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->getModifiable(pParticleID,
                            d_mpmLabels->pParticleIDLabel_preReloc,
                            pset);
      new_dw->getModifiable(pMass, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->getModifiable(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);
      new_dw->getModifiable(pDisp, d_mpmLabels->pDispLabel_preReloc, pset);
      new_dw->getModifiable(pStress, d_mpmLabels->pStressLabel_preReloc, pset);
      new_dw->getModifiable(pcolor, d_mpmLabels->pColorLabel_preReloc, pset);
      new_dw->getModifiable(pVolume, d_mpmLabels->pVolumeLabel_preReloc, pset);
      new_dw->getModifiable(pVelocity,
                            d_mpmLabels->pVelocityLabel_preReloc,
                            pset);
      new_dw->getModifiable(pextforce,
                            d_mpmLabels->pExtForceLabel_preReloc,
                            pset);
      new_dw->getModifiable(ptemp,
                            d_mpmLabels->pTemperatureLabel_preReloc,
                            pset);
      new_dw->getModifiable(ptempP,
                            d_mpmLabels->pTempPreviousLabel_preReloc,
                            pset);
      new_dw->getModifiable(pref, d_mpmLabels->pRefinedLabel_preReloc, pset);
      new_dw->getModifiable(ploc,
                            d_mpmLabels->pLocalizedMPMLabel_preReloc,
                            pset);
      new_dw->getModifiable(pvelgrad,
                            d_mpmLabels->pVelGradLabel_preReloc,
                            pset);
      new_dw->getModifiable(pF, d_mpmLabels->pDefGradLabel_preReloc, pset);

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
      new_dw->put(pParticleIDtmp, d_mpmLabels->pParticleIDLabel_preReloc, true);
      new_dw->put(pXtmp, d_mpmLabels->pXLabel_preReloc, true);
      new_dw->put(pVoltmp, d_mpmLabels->pVolumeLabel_preReloc, true);
      new_dw->put(pveltmp, d_mpmLabels->pVelocityLabel_preReloc, true);
      new_dw->put(pextFtmp, d_mpmLabels->pExtForceLabel_preReloc, true);
      new_dw->put(pMasstmp, d_mpmLabels->pMassLabel_preReloc, true);
      new_dw->put(ptemptmp, d_mpmLabels->pTemperatureLabel_preReloc, true);
      new_dw->put(ptempPtmp, d_mpmLabels->pTempPreviousLabel_preReloc, true);
      new_dw->put(pSizetmp, d_mpmLabels->pSizeLabel_preReloc, true);
      new_dw->put(pDisptmp, d_mpmLabels->pDispLabel_preReloc, true);
      new_dw->put(pstrstmp, d_mpmLabels->pStressLabel_preReloc, true);
      new_dw->put(pcolortmp, d_mpmLabels->pColorLabel_preReloc, true);
      new_dw->put(pFtmp, d_mpmLabels->pDefGradLabel_preReloc, true);
      new_dw->put(preftmp, d_mpmLabels->pRefinedLabel_preReloc, true);
      new_dw->put(ploctmp, d_mpmLabels->pLocalizedMPMLabel_preReloc, true);
      new_dw->put(pvgradtmp, d_mpmLabels->pVelGradLabel_preReloc, true);
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
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleScaleFactor");

  Task* t = scinew Task("MPM::computeParticleScaleFactor",
                        this,
                        &SerialMPM::computeParticleScaleFactor);

  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel_preReloc, Ghost::None);
  t->computes(d_mpmLabels->pScaleFactorLabel_preReloc);

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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      ParticleSubset* pset  = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pScaleFactor;
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      new_dw->allocateAndPut(pScaleFactor,
                             d_mpmLabels->pScaleFactorLabel_preReloc,
                             pset);

      if (dataArchiver->isOutputTimestep()) {
        Vector dx = patch->dCell();

        if (d_mpmFlags->d_interpolatorType != "cpdi" &&
            d_mpmFlags->d_interpolatorType != "cpti") {
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel_preReloc, pset);
          for (auto pidx : *pset) {
            pScaleFactor[pidx] =
              ((pDefGrad[pidx] * pSize[pidx]) *
               Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));
          }
        } else {
          for (ParticleSubset::iterator iter = pset->begin();
               iter != pset->end();
               iter++) {
            particleIndex idx = *iter;
            pScaleFactor[idx] =
              (pSize[idx] * Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));

          } // for particles
        }
      } // isOutputTimestep
    }   // matls
  }     // patches
}

/*!----------------------------------------------------------------------
 * scheduleCheckNeedAddMPMMaterial
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleCheckNeedAddMPMMaterial(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "MPM::scheduleCheckNeedAddMPMMaterial");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("MPM::checkNeedAddMPMMaterial",
                        this,
                        &SerialMPM::checkNeedAddMPMMaterial);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->scheduleCheckNeedAddMPMMaterial(t, mpm_matl, patches);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * checkNeedAddMPMMaterial
 *-----------------------------------------------------------------------*/
void
SerialMPM::checkNeedAddMPMMaterial(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  // printSchedule(patches, cout_doing, "MPM::checkNeedAddMPMMaterial");

  for (int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->checkNeedAddMPMMaterial(patches, mpm_matl, old_dw, new_dw);
  }
}

/*!----------------------------------------------------------------------
 * scheduleFinalParticleUpdate
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleFinalParticleUpdate(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)

{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleFinalParticleUpdate");

  Task* t = scinew Task("MPM::finalParticleUpdate",
                        this,
                        &SerialMPM::finalParticleUpdate);

  t->requires(Task::OldDW, d_materialManager->get_delt_label());

  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, d_mpmLabels->pdTdtLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pLocalizedMPMLabel_preReloc, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pMassLabel_preReloc, gnone);

  t->modifies(d_mpmLabels->pTemperatureLabel_preReloc);

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
    old_dw->get(delT, d_materialManager->get_delt_label(), getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<int> pLocalized;
      constParticleVariable<double> pdTdt, pMass_new;
      ParticleVariable<double> pTemp_new;

      ParticleSubset* pset   = old_dw->getParticleSubset(matID, patch);
      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

      new_dw->get(pdTdt, d_mpmLabels->pdTdtLabel, pset);
      new_dw->get(pMass_new, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->get(pLocalized, d_mpmLabels->pLocalizedMPMLabel_preReloc, pset);

      new_dw->getModifiable(pTemp_new,
                            d_mpmLabels->pTemperatureLabel_preReloc,
                            pset);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        pTemp_new[idx] += pdTdt[idx] * delT;

        // Delete particles whose mass is too small (due to combustion),
        // whose pLocalized flag has been set to -999 or who have a negative
        // temperature
        if ((pMass_new[idx] <= d_mpmFlags->d_minPartMass) ||
            pTemp_new[idx] < 0. || (pLocalized[idx] == -999)) {
          if (d_mpmFlags->d_erosionAlgorithm != "none") {
            delset->addParticle(idx);
          }
        }

      } // particles

      new_dw->deleteParticles(delset);

    } // materials
  }   // patches
}

/*!----------------------------------------------------------------------
 * scheduleRefine
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  printSchedule(patches, cout_doing, "MPM::scheduleRefine");
  Task* t = scinew Task("SerialMPM::refine", this, &SerialMPM::refine);

  t->computes(d_mpmLabels->pXLabel);
  t->computes(d_mpmLabels->p_qLabel);
  t->computes(d_mpmLabels->pDispLabel);
  t->computes(d_mpmLabels->pMassLabel);
  t->computes(d_mpmLabels->pVolumeLabel);
  t->computes(d_mpmLabels->pTemperatureLabel);
  t->computes(d_mpmLabels->pTempPreviousLabel); // for therma  stresm analysis
  t->computes(d_mpmLabels->pdTdtLabel);
  t->computes(d_mpmLabels->pVelocityLabel);
  t->computes(d_mpmLabels->pExternalForceLabel);
  t->computes(d_mpmLabels->pParticleIDLabel);
  t->computes(d_mpmLabels->pStressLabel);
  t->computes(d_mpmLabels->pSizeLabel);
  t->computes(d_mpmLabels->NC_CCweightLabel);
  t->computes(d_materialManager->get_delt_label(), getLevel(patches));

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

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);

    // For deformation gradient computer
    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    // For constitutive models
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Basic damage model related stuff
    if (mpm_matl->doBasicDamage()) {
      Vaango::BasicDamageModel* basicDamageModel =
        mpm_matl->getBasicDamageModel();
      basicDamageModel->addInitialComputesAndRequires(t,
                                                      mpm_matl,
                                                      patches,
                                                      d_mpmLabels);
    }
  }

  sched->addTask(t, patches, d_materialManager->allMPMMaterials());
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
    new_dw->allocateAndPut(NC_CCweight,
                           d_mpmLabels->NC_CCweightLabel,
                           0,
                           patch);
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
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

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

        new_dw->allocateAndPut(pX, d_mpmLabels->pXLabel, pset);
        new_dw->allocateAndPut(p_q, d_mpmLabels->p_qLabel, pset);
        new_dw->allocateAndPut(pMass, d_mpmLabels->pMassLabel, pset);
        new_dw->allocateAndPut(pVolume, d_mpmLabels->pVolumeLabel, pset);
        new_dw->allocateAndPut(pVelocity, d_mpmLabels->pVelocityLabel, pset);
        new_dw->allocateAndPut(pTemperature,
                               d_mpmLabels->pTemperatureLabel,
                               pset);
        new_dw->allocateAndPut(pTempPrev,
                               d_mpmLabels->pTempPreviousLabel,
                               pset);
        new_dw->allocateAndPut(pExternalForce,
                               d_mpmLabels->pExternalForceLabel,
                               pset);
        new_dw->allocateAndPut(pID, d_mpmLabels->pParticleIDLabel, pset);
        new_dw->allocateAndPut(pDisp, d_mpmLabels->pDispLabel, pset);
        if (d_mpmFlags->d_useLoadCurves) {
          new_dw->allocateAndPut(pLoadCurve,
                                 d_mpmLabels->pLoadCurveIDLabel,
                                 pset);
        }
        new_dw->allocateAndPut(pSize, d_mpmLabels->pSizeLabel, pset);

        // Init deformation gradient
        d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

        // Init constitutive model
        mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                           mpm_matl,
                                                           new_dw);

        // Initialize basic damage models
        if (mpm_matl->doBasicDamage()) {
          mpm_matl->getBasicDamageModel()->initializeDamageData(patch,
                                                                mpm_matl,
                                                                new_dw,
                                                                lb);
        }

#if 0
        if(d_mpmFlags->d_withColor) {
          ParticleVariable<double> pcolor;
          int index = mpm_matl->getDWIndex();
          ParticleSubset* pset = new_dw->getParticleSubset(index, patch);
          setParticleDefault<double>(pcolor, d_mpmLabels->pColorLabel, pset, new_dw, 0.0);
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
 *   Schedule to mark d_mpmFlags for AMR regridding
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
    task->requires(Task::NewDW, d_mpmLabels->pXLabel, Ghost::AroundCells, 0);
  } else {
    task->requires(Task::NewDW,
                   d_materialManager->get_refineFlag_label(),
                   0,
                   Task::FineLevel,
                   d_materialManager->refineFlagMaterials(),
                   Task::NormalDomain,
                   Ghost::None,
                   0);
  }
  task->modifies(d_materialManager->get_refineFlag_label(),
                 d_materialManager->refineFlagMaterials());
  task->modifies(d_materialManager->get_refinePatchFlag_label(),
                 d_materialManager->refineFlagMaterials());
  sched->addTask(task,
                 coarseLevel->eachPatch(),
                 d_materialManager->allMPMMaterials());
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

      new_dw->getModifiable(refineFlag,
                            d_materialManager->get_refineFlag_label(),
                            0,
                            coarsePatch);
      new_dw->get(refinePatchFlag,
                  d_materialManager->get_refinePatchFlag_label(),
                  0,
                  coarsePatch);

      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);

      // coarsen the fineLevel flag
      for (int i = 0; i < finePatches.size(); i++) {
        const Patch* finePatch = finePatches[i];

        IntVector cl, ch, fl, fh;
        getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);

        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }
        constCCVariable<int> fineErrorFlag;
        new_dw->getRegion(fineErrorFlag,
                          d_materialManager->get_refineFlag_label(),
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
 *   Schedule to mark initial d_mpmFlags for AMR regridding
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
    new_dw->getModifiable(refineFlag,
                          d_materialManager->get_refineFlag_label(),
                          0,
                          patch);
    new_dw->get(refinePatchFlag,
                d_materialManager->get_refinePatchFlag_label(),
                0,
                patch);

    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

    for (int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
      constParticleVariable<Point> pX;
      new_dw->get(pX, d_mpmLabels->pXLabel, pset);

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        refineFlag[patch->getLevel()->getCellIndex(pX[*iter])] = true;
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
 * scheduleTotalParticleCount
 *   Diagnostic task: compute the total number of particles
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleTotalParticleCount(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("SerialMPM::totalParticleCount",
                        this,
                        &SerialMPM::totalParticleCount);
  t->computes(d_mpmLabels->partCountLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * totalParticleCount
 *   Diagnostic task: compute the total number of particles
 *-----------------------------------------------------------------------*/
void
SerialMPM::totalParticleCount(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* matls,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch      = patches->get(p);
    long int totalParticles = 0;

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
      int matID             = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      int numParticles     = pset->end() - pset->begin();

      totalParticles += numParticles;
    }
    new_dw->put(sumlong_vartype(totalParticles), d_mpmLabels->partCountLabel);
  }
}

/*!----------------------------------------------------------------------
 * addMaterial
 *   For adding materials mid-Simulation
 *-----------------------------------------------------------------------*/
void
SerialMPM::addMaterial(const ProblemSpecP& prob_spec,
                       MaterialManagerP& matManager)
{
  d_recompile = true;
  ProblemSpecP mat_ps =
    prob_spec->findBlockWithAttribute("MaterialProperties", "add");

  std::string attr = "";
  mat_ps->getAttribute("add", attr);

  if (attr == "true") {
    ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");
    for (ProblemSpecP ps = mpm_mat_ps->findBlock("material"); ps != 0;
         ps              = ps->findNextBlock("material")) {
      // Create and register as an MPM material
      // MPMMaterial *mat = scinew MPMMaterial(ps, grid, d_materialManager,
      // d_mpmFlags);
      std::shared_ptr<MPMMaterial> mat =
        std::make_shared<MPMMaterial>(ps, matManager, d_mpmFlags, false);
      matManager->registerMaterial("MPM", std::move(mat));
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleInitializeAddedMaterial
 *-----------------------------------------------------------------------*/
void
SerialMPM::scheduleInitializeAddedMaterial(const LevelP& level,
                                           SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();
  printSchedule(patches, cout_doing, "MPM::scheduleInitializeAddedMaterial");

  Task* t = scinew Task("SerialMPM::actuallyInitializeAddedMaterial",
                        this,
                        &SerialMPM::actuallyInitializeAddedMaterial);

  int numALLMatls          = d_materialManager->getNumMaterials();
  int numMPMMatls          = d_materialManager->getNumMaterials("MPM");
  MaterialSubset* add_matl = scinew MaterialSubset();
  std::cout << "Added Material = " << numALLMatls - 1 << "\n";
  add_matl->add(numALLMatls - 1);
  add_matl->addReference();

  t->computes(d_materialManager->get_delt_label(), level.get_rep());

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
  t->computes(d_mpmLabels->pFiberDirLabel, add_matl);
  t->computes(d_mpmLabels->pRefinedLabel, add_matl);
  if (!d_mpmFlags->d_doGridReset) {
    t->computes(d_mpmLabels->gDisplacementLabel);
  }

  // Add initialization of body force and coriolis importance terms
  t->computes(d_mpmLabels->pCoriolisImportanceLabel);
  t->computes(d_mpmLabels->pBodyForceAccLabel);

  // Computes accumulated strain energy
  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    t->computes(d_mpmLabels->AccStrainEnergyLabel);
  }

  MPMMaterial* mpm_matl =
    d_materialManager->getMaterial("MPM", numMPMMatls - 1);

  // Add vel grad/def grad computes
  d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

  // Add cm computes
  ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
  cm->addInitialComputesAndRequires(t, mpm_matl, patches);

  // Add damage model computes
  if (mpm_matl->doBasicDamage()) {
    Vaango::BasicDamageModel* basicDamageModel =
      mpm_matl->getBasicDamageModel();
    basicDamageModel->addInitialComputesAndRequires(t,
                                                    mpm_matl,
                                                    patches,
                                                    d_mpmLabels);
  }

  sched->addTask(t, patches, d_materialManager->allMPMMaterials());

  // The task will have a reference to add_matl
  if (add_matl->removeReference()) {
    delete add_matl; // shouln't happen, but...
  }
}

/*!----------------------------------------------------------------------
 * actuallyInitializeAddedMaterial
 *-----------------------------------------------------------------------*/
void
SerialMPM::actuallyInitializeAddedMaterial(const ProcessorGroup*,
                                           const PatchSubset* patches,
                                           const MaterialSubset* /*matls*/,
                                           DataWarehouse*,
                                           DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              cout_doing,
              "Doing actuallyInitializeAddedMaterial");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    std::cout << "num MPM Matls = " << numMPMMatls << "\n";
    CCVariable<short int> cellNAPID;
    int m                 = numMPMMatls - 1;
    MPMMaterial* mpm_matl = d_materialManager->getMaterial("MPM", m);
    new_dw->unfinalize();
    // particleIndex numParticles =
    mpm_matl->createParticles(cellNAPID, patch, new_dw);

    // Initialize deformation gradient of added particles
    d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

    // Initialize constitutive models of added particles
    mpm_matl->getConstitutiveModel()->initializeCMData(patch, mpm_matl, new_dw);

    // Initialize basic damage models of added particles
    if (mpm_matl->doBasicDamage()) {
      mpm_matl->getBasicDamageModel()->initializeDamageData(patch,
                                                            mpm_matl,
                                                            new_dw,
                                                            lb);
    }

    new_dw->refinalize();
  }
}

/* _____________________________________________________________________
   Purpose:   Set variables that are normally set during the initialization
   phase, but get wiped clean when you restart
   _____________________________________________________________________*/
void
SerialMPM::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
}

/*!----------------------------------------------------------------------
 * restartInitialize
 *-----------------------------------------------------------------------*/
void
SerialMPM::restartInitialize()
{
  cout_doing << "Doing restartInitialize\t\t\t\t\t MPM"
             << "\n";

  if (d_analysisModules.size() != 0) {
    for (auto module : d_analysisModules) {
      module->restartInitialize();
    }
  }
}

/*!----------------------------------------------------------------------
 * needRecompile
 *-----------------------------------------------------------------------*/
bool
SerialMPM::needRecompile(double, double, const GridP&)
{
  if (d_recompile) {
    d_recompile = false;
    return true;
  } else {
    return false;
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
  for (auto particle : *pset) {
    pvar[particle] = val;
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
