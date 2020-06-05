/*
 * The MIT License
 *
 * Copyright (c) 2018-2020 Parresia Research Limited, New Zealand
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
#include <CCA/Components/MPM/UofU_MPM.h>
#include <CCA/Components/MPM/ConstitutiveModel/DamageModels/BasicDamageModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/Contact.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Components/MPM/MPMBoundCond.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/PhysicalBC/ForceBC.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/MomentBC.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/VelocityBC.h>
#include <CCA/Components/MPM/MMS/MMS.h>
#include <CCA/Components/Regridder/PerPatchVars.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/SimulationState.h>
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
#include <Core/Thread/Mutex.h>
#include <Core/Util/DebugStream.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

//#define TIME_COMPUTE_STRESS
//#define CHECK_ISFINITE
//#define DEBUG_WITH_PARTICLE_ID

using namespace Uintah;

//__________________________________
//  To turn on debug d_flags
//  csh/tcsh : setenv SCI_DEBUG "UofU_MPM_Doing:+,UofU_MPM_Debug:+"
//  bash     : export SCI_DEBUG="UofU_MPM_Doing:+,UofU_MPM_Debug:+"
//  default is OFF

static DebugStream cout_doing("UofU_MPM_Doing", false);
static DebugStream cout_dbg("UofU_MPM_Debug", false);

// From ThreadPool.cc:  Used for syncing cerr'ing so it is easier to read.
extern Mutex cerrLock;

// Constants and static functions
constexpr double SMALL_NUM_MPM = 1.0e-200;
const Uintah::Matrix3 UofU_MPM::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
const Uintah::Matrix3 UofU_MPM::Zero(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

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

UofU_MPM::UofU_MPM(const ProcessorGroup* myworld)
  : MPMCommon(myworld)
  , UintahParallelComponent(myworld)
{
  d_labels = scinew MPMLabel();
  d_flags = scinew MPMFlags(myworld);

  d_nextOutputTime = 0.;
  d_loadCurveIndex = 0;

  d_contactModel = 0;
  d_numGhostParticles = 1;
  d_numGhostNodes = 1;
  d_dataArchiver = 0;
}

UofU_MPM::~UofU_MPM()
{
  delete d_labels;
  delete d_flags;
  delete d_contactModel;
  MPMPhysicalBCFactory::clean();
}

/*!----------------------------------------------------------------------
 * problemSetup
 *-----------------------------------------------------------------------*/
void
UofU_MPM::problemSetup(const ProblemSpecP& prob_spec,
                       const ProblemSpecP& restart_prob_spec, GridP& grid,
                       SimulationStateP& sharedState)
{
  cout_doing << "Doing problemSetup\t\t\t\t\t MPM" << endl;
  d_sharedState = sharedState;
  dynamic_cast<Scheduler*>(getPort("scheduler"))
    ->setPositionVar(d_labels->pXLabel);

  d_dataArchiver = dynamic_cast<Output*>(getPort("output"));
  if (!d_dataArchiver) {
    throw InternalError("MPM:couldn't get output port", __FILE__, __LINE__);
  }

  ProblemSpecP mat_ps = prob_spec;
  if (restart_prob_spec) {
    mat_ps = restart_prob_spec;
  }

  ProblemSpecP mpm_soln_ps = mat_ps->findBlock("MPM");
  if (!mpm_soln_ps) {
    ostringstream warn;
    warn << "**ERROR**: MPM:\n missing MPM section in the input file\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // Read all MPM d_flags (look in MPMFlags.cc)
  d_flags->readMPMFlags(mat_ps, d_dataArchiver);
  if (d_flags->d_integratorType == "implicit") {
    throw ProblemSetupException("Can't use implicit integration with -uofu_mpm",
                                __FILE__, __LINE__);
  }

  // convert text representation of face into FaceType
  for (auto faceName : d_flags->d_boundaryTractionFaceStrings) {
    Patch::FaceType face = Patch::invalidFace;
    for (Patch::FaceType ft = Patch::startFace; ft <= Patch::endFace;
         ft = Patch::nextFace(ft)) {
      if (Patch::getFaceName(ft) == faceName)
        face = ft;
    }
    if (face != Patch::invalidFace) {
      d_boundaryTractionFaces.push_back(face);
    } else {
      std::cerr << "**WARNING**: Ignoring unknown face '" << faceName
                << "' in input file. \n";
    }
  }

  if (d_flags->d_8or27 == 8) {
    d_numGhostParticles = 1;
    d_numGhostNodes = 1;
  } else {
    d_numGhostParticles = 2;
    d_numGhostNodes = 2;
  }

  if (d_flags->d_prescribeDeformation) {
    readPrescribedDeformations(d_flags->d_prescribedDeformationFile);
  }

  d_sharedState->setParticleGhostLayer(Ghost::AroundNodes, d_numGhostParticles);

  MPMPhysicalBCFactory::create(mat_ps, grid, d_flags);

  d_contactModel =
    ContactFactory::create(UintahParallelComponent::d_myworld, mat_ps,
                           d_sharedState, d_labels, d_flags);

  // Creates MPM material w/ constitutive models and damage models
  materialProblemSetup(mat_ps, grid, d_sharedState, d_flags);

}

/*!----------------------------------------------------------------------
 * outputProblemSpec
 *-----------------------------------------------------------------------*/
void
UofU_MPM::outputProblemSpec(ProblemSpecP& root_ps)
{
  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP d_flags_ps = root->appendChild("MPM");
  d_flags->outputProblemSpec(d_flags_ps);

  ProblemSpecP mat_ps = 0;
  mat_ps = root->findBlockWithOutAttribute("MaterialProperties");

  if (mat_ps == 0)
    mat_ps = root->appendChild("MaterialProperties");

  ProblemSpecP mpm_ps = mat_ps->appendChild("MPM");
  for (int i = 0; i < d_sharedState->getNumMPMMatls(); i++) {
    MPMMaterial* mat = d_sharedState->getMPMMaterial(i);
    ProblemSpecP cm_ps = mat->outputProblemSpec(mpm_ps);
  }

  d_contactModel->outputProblemSpec(mpm_ps);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps = physical_bc_ps->appendChild("MPM");
  for (auto& physicalBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    physicalBC->outputProblemSpec(mpm_ph_bc_ps);
  }
}

/*!----------------------------------------------------------------------
 * readPrescribedDeformations
 *-----------------------------------------------------------------------*/
void
UofU_MPM::readPrescribedDeformations(string filename)
{
  if (filename != "") {
    std::ifstream is(filename.c_str());
    if (!is) {
      std::ostringstream msg;
      msg << "**ERROR** Could not open prescribed deformation file '" 
          << filename << "'\n";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    } else {
      std::cout << "**INFO** Reading prescribed deformations from file:"
                << filename << "\n";
    }
    double t0(-1.e9);
    while (is) {
      double t1, F11, F12, F13, F21, F22, F23, F31, F32, F33, Theta, a1, a2, a3;
      is >> t1 >> F11 >> F12 >> F13 >> F21 >> F22 >> F23 >> F31 >> F32 >> F33 >>
        Theta >> a1 >> a2 >> a3;
      if (is) {
        if (t1 <= t0) {
          std::ostringstream msg;
          msg << "**ERROR** Time in prescribed deformation file '" 
              << filename << "' is not monotonically increasing.\n";
          throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
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
        "**ERROR** Failed to generate valid precribed deformation profile",
        __FILE__, __LINE__);
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleInitialize
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!d_flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()))
    return;

  auto patches = level->eachPatch();
  printSchedule(patches, cout_doing, "MPM::scheduleInitialize");

  Task* t =
    scinew Task("UofU_MPM::actuallyInitialize", this, &UofU_MPM::actuallyInitialize);

  MaterialSubset* zeroth_matl = scinew MaterialSubset();
  zeroth_matl->add(0);
  zeroth_matl->addReference();
  t->computes(d_labels->pCellNAPIDLabel, zeroth_matl);
  t->computes(d_labels->NC_CCweightLabel, zeroth_matl);

  t->computes(d_sharedState->get_delt_label(), level.get_rep());

  t->computes(d_labels->partCountLabel);
  t->computes(d_labels->pXLabel);
  t->computes(d_labels->pDispLabel);
  t->computes(d_labels->pMassLabel);
  t->computes(d_labels->pVolumeLabel);
  t->computes(d_labels->pTemperatureLabel);
  t->computes(d_labels->pTempPreviousLabel);
  t->computes(d_labels->pVelocityLabel);
  t->computes(d_labels->pExternalForceLabel);
  t->computes(d_labels->pParticleIDLabel);
  t->computes(d_labels->pStressLabel);
  t->computes(d_labels->pSizeLabel);

  if (d_flags->d_artificialViscosity) {
    t->computes(d_labels->p_qLabel);
  }

  // Debugging Scalar
  if (d_flags->d_withColor) {
    t->computes(d_labels->pColorLabel);
  }

  // Computes the load curve ID associated with each particle
  if (d_flags->d_useLoadCurves) {
    t->computes(d_labels->pLoadCurveIDLabel);
  }

  // Computes accumulated strain energy
  if (d_flags->d_reductionVars->accStrainEnergy) {
    t->computes(d_labels->AccStrainEnergyLabel);
  }

  int numMPM = d_sharedState->getNumMPMMatls();
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    const MaterialSubset* matlset = mpm_matl->thisMaterial();

    // Add vel grad/def grad computes
    t->computes(d_labels->pVelGradLabel,      matlset);
    t->computes(d_labels->pDefGradLabel,      matlset);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
      t->computes(d_labels->pRemoveLabel,       matlset);
      t->computes(d_labels->pPolarDecompRLabel, matlset);
    }

    // Add constitutive model computes
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  // Add initialization of body force and coriolis importance terms
  // These are initialized to zero in ParticleCreator
  t->computes(d_labels->pCoriolisImportanceLabel);
  t->computes(d_labels->pBodyForceAccLabel);

  // Add task to scheduler
  sched->addTask(t, patches, d_sharedState->allMPMMaterials());

  // The task will have a reference to zeroth_matl
  if (zeroth_matl->removeReference())
    delete zeroth_matl; // shouldn't happen, but...

  // Print particle count
  schedulePrintParticleCount(level, sched);

  // Compute initial stresses due to body forces and recompute the initial
  // deformation
  // gradient
  if (d_flags->d_initializeStressFromBodyForce) {
    scheduleInitializeStressAndDefGradFromBodyForce(level, sched);
  }

  // Schedule the initialization of pressure BCs per particle
  if (d_flags->d_useLoadCurves) {
    if (MPMPhysicalBCFactory::mpmPhysicalBCs.size() > 0) {
      string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[0]->getType();
      if (bcs_type == "Pressure") {
        scheduleInitializePressureBCs(level, sched);
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * actuallyInitialize
 *-----------------------------------------------------------------------*/
void
UofU_MPM::actuallyInitialize(const ProcessorGroup*, const PatchSubset* patches,
                             const MaterialSubset* matls, DataWarehouse*,
                             DataWarehouse* new_dw)
{
  particleIndex totalParticles = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_labels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight, d_labels->NC_CCweightLabel, 0, patch);

    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NC_CCweight.initialize(0.125);
    for (auto face = Patch::startFace; face <= Patch::endFace;
         face = Patch::nextFace(face)) {
      int mat_id = 0;

      if (patch->haveBC(face, mat_id, "symmetry", "Symmetric")) {
        for (CellIterator iter = patch->getFaceIterator(face, Patch::FaceNodes);
             !iter.done(); iter++) {
          NC_CCweight[*iter] = 2.0 * NC_CCweight[*iter];
        }
      }
    }

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();


      particleIndex numParticles =
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      // Initialize deformation gradient and its polar decomposition
      ParticleSubset* pset = new_dw->getParticleSubset(mpm_matl->getDWIndex(), patch);
      ParticleVariable<Matrix3> pVelGrad, pDefGrad;

      ParticleVariable<int>     pRemove;
      ParticleVariable<Matrix3> pPolarDecompR;
      new_dw->allocateAndPut(pVelGrad,       d_labels->pVelGradLabel,      pset);
      new_dw->allocateAndPut(pDefGrad,       d_labels->pDefGradLabel,      pset);

      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        new_dw->allocateAndPut(pRemove,        d_labels->pRemoveLabel,       pset);
        new_dw->allocateAndPut(pPolarDecompR,  d_labels->pPolarDecompRLabel, pset);
        for (auto particle : *pset) {
          pRemove[particle] = 0;
          pVelGrad[particle] = Zero;
          pDefGrad[particle] = Identity;
          pPolarDecompR[particle] = Identity;
        }
      } else {
        for (auto particle : *pset) {
          pVelGrad[particle] = Zero;
          pDefGrad[particle] = Identity;
        }
      }

      // Initialize constitutive models
      mpm_matl->getConstitutiveModel()->initializeCMData(patch, mpm_matl,
                                                         new_dw);
    }

    IntVector num_extra_cells = patch->getExtraCells();
    IntVector periodic = patch->getLevel()->getPeriodicBoundaries();
    auto interp_type = d_flags->d_interpolatorType;
    if (interp_type == "linear" && num_extra_cells != IntVector(0, 0, 0)) {
      if (!d_flags->d_withICE) {
        std::ostringstream msg;
        msg << "\n ERROR: When using <interpolator>linear</interpolator> \n"
            << " you should also use <extraCells>[0,0,0]</extraCells> \n"
            << " unless you are running an MPMICE case.\n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
    } else if (((interp_type == "gimp" || interp_type == "3rdorderBS" ||
                 interp_type == "cpdi") &&
                ((num_extra_cells + periodic) != IntVector(1, 1, 1) &&
                 ((num_extra_cells + periodic) != IntVector(1, 1, 0) &&
                  d_flags->d_axisymmetric)))) {
      ostringstream msg;
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
    if (d_flags->d_axisymmetric) {
      IntVector patchLowNode = patch->getNodeLowIndex();
      IntVector patchHighNode = patch->getNodeHighIndex();
      int num_cells_in_theta = (patchHighNode.z() - patchLowNode.z()) - 1;
      if (num_cells_in_theta > 1) {
        ostringstream msg;
        msg << "\n ERROR: When using <axisymmetric>true</axisymmetric> \n"
            << "the grid can only have one cell in the circumferential "
               "direction.\n";
        throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
      }
    }
  }

  // Initialize the accumulated strain energy
  if (d_flags->d_reductionVars->accStrainEnergy) {
    new_dw->put(max_vartype(0.0), d_labels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_labels->partCountLabel);
}

/*!----------------------------------------------------------------------
 * schedulePrintParticleCount
 *-----------------------------------------------------------------------*/
void
UofU_MPM::schedulePrintParticleCount(const LevelP& level, SchedulerP& sched)
{
  Task* t =
    scinew Task("UofU_MPM::printParticleCount", this, &UofU_MPM::printParticleCount);
  t->requires(Task::NewDW, d_labels->partCountLabel);
  t->setType(Task::OncePerProc);
  sched->addTask(t, sched->getLoadBalancer()->getPerProcessorPatchSet(level),
                 d_sharedState->allMPMMaterials());
}

/*!----------------------------------------------------------------------
 * printParticleCount
 *-----------------------------------------------------------------------*/
void
UofU_MPM::printParticleCount(const ProcessorGroup* pg, const PatchSubset*,
                             const MaterialSubset*, DataWarehouse*,
                             DataWarehouse* new_dw)
{
  sumlong_vartype pcount;
  new_dw->get(pcount, d_labels->partCountLabel);

  if (pg->myrank() == 0) {
    std::cerr << "**INFO** Created " << (long)pcount << " total particles\n";
  }
}

/*!------------------------------------------------------------------------
 * Schedule the initialization of the stress and deformation gradient
 * based on the body forces (which also have to be computed)
 *------------------------------------------------------------------------*/
void
UofU_MPM::scheduleInitializeStressAndDefGradFromBodyForce(const LevelP& level,
                                                          SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();
  printSchedule(patches, cout_doing,
                "MPM::initializeStressAndDefGradFromBodyForce");

  // First compute the body force
  Task* t1 = scinew Task("UofU_MPM::initializeBodyForce", this,
                         &UofU_MPM::initializeBodyForce);
  t1->requires(Task::NewDW, d_labels->pXLabel, Ghost::None);
  t1->modifies(d_labels->pBodyForceAccLabel);
  sched->addTask(t1, patches, d_sharedState->allMPMMaterials());

  // Compute the stress and deformation gradient only for selected
  // constitutive models that have a "initializeWithBodyForce" flag as true.
  // This is because a more general implementation is quite involved and
  // not worth the effort at this time. BB
  Task* t2 = scinew Task("UofU_MPM::initializeStressAndDefGradFromBodyForce", this,
                         &UofU_MPM::initializeStressAndDefGradFromBodyForce);

  t2->requires(Task::NewDW, d_labels->pXLabel, Ghost::None);
  t2->requires(Task::NewDW, d_labels->pBodyForceAccLabel, Ghost::None);
  t2->modifies(d_labels->pStressLabel);
  t2->modifies(d_labels->pDefGradLabel);
  sched->addTask(t2, patches, d_sharedState->allMPMMaterials());
}

/*!------------------------------------------------------------------------
 * Actually initialize the body force acceleration
 *-------------------------------------------------------------------------*/
void
UofU_MPM::initializeBodyForce(const ProcessorGroup*, const PatchSubset* patches,
                              const MaterialSubset* matls, DataWarehouse*,
                              DataWarehouse* new_dw)
{
  // Get the MPM d_flags and make local copies
  Uintah::Point rotation_center = d_flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis = d_flags->d_coordRotationAxis;
  double rotation_speed = d_flags->d_coordRotationSpeed;

  // Compute angular velocity vector (omega)
  Uintah::Vector omega = rotation_axis * rotation_speed;

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleBodyForce");

    // Loop thru materials
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      // Get the material ID
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      new_dw->getModifiable(pBodyForceAcc, d_labels->pBodyForceAccLabel, pset);

      // Get the position data
      constParticleVariable<Point> pPosition;
      new_dw->get(pPosition, d_labels->pXLabel, pset);

      // Iterate over the particles
      for (auto iter = pset->begin(); iter != pset->end(); iter++) {
        particleIndex pidx = *iter;

        // Compute the body force acceleration (g)
        // Just use gravity if rotation is off
        pBodyForceAcc[pidx] = d_flags->d_gravity;

        // If rotating add centrifugal force
        if (d_flags->d_useCoordRotation) {
          // Compute the centrifugal term (omega x omega x r)
          // Simplified version where body ref point is not needed
          Vector rVec = pPosition[pidx] - rotation_center;
          Vector omega_x_r = Uintah::Cross(omega, rVec);
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
UofU_MPM::initializeStressAndDefGradFromBodyForce(const ProcessorGroup*,
                                                  const PatchSubset* patches,
                                                  const MaterialSubset* matls,
                                                  DataWarehouse*,
                                                  DataWarehouse* new_dw)
{
  // Loop over patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing,
              "Doing initializeStressAndDefGradFromBodyForce");

    // Loop over materials
    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

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
UofU_MPM::scheduleInitializePressureBCs(const LevelP& level, SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();

  d_loadCurveIndex = scinew MaterialSubset();
  d_loadCurveIndex->add(0);
  d_loadCurveIndex->addReference();

  int pressureBCId = 0;
  for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
       ii++) {
    string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();
    if (bcs_type == "Pressure") {
      d_loadCurveIndex->add(pressureBCId++);
    }
  }
  if (pressureBCId > 0) {
    printSchedule(patches, cout_doing, "MPM::countMaterialPointsPerLoadCurve");
    printSchedule(patches, cout_doing, "MPM::scheduleInitializePressureBCs");
    // Create a task that calculates the total number of particles
    // associated with each load curve.
    Task* t = scinew Task("UofU_MPM::countMaterialPointsPerLoadCurve", this,
                          &UofU_MPM::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, d_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_labels->materialPointsPerLoadCurveLabel, d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, patches, d_sharedState->allMPMMaterials());

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("UofU_MPM::initializePressureBC", this,
                    &UofU_MPM::initializePressureBC);
    t->requires(Task::NewDW, d_labels->pXLabel, Ghost::None);
    t->requires(Task::NewDW, d_labels->pSizeLabel, Ghost::None);
    t->requires(Task::NewDW, d_labels->pDispLabel, Ghost::None);
    t->requires(Task::NewDW, d_labels->pDefGradLabel, Ghost::None);
    t->requires(Task::NewDW, d_labels->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW, d_labels->materialPointsPerLoadCurveLabel,
                d_loadCurveIndex, Task::OutOfDomain, Ghost::None);
    t->modifies(d_labels->pExternalForceLabel);
    if (d_flags->d_useCBDI) {
      t->computes(d_labels->pExternalForceCorner1Label);
      t->computes(d_labels->pExternalForceCorner2Label);
      t->computes(d_labels->pExternalForceCorner3Label);
      t->computes(d_labels->pExternalForceCorner4Label);
    }
    sched->addTask(t, patches, d_sharedState->allMPMMaterials());
  }

  if (d_loadCurveIndex->removeReference())
    delete d_loadCurveIndex;
}

/*!----------------------------------------------------------------------
 * countMaterialPointsPerLoadCurve
 *   Calculate the number of material points per load curve
 *-----------------------------------------------------------------------*/
void
UofU_MPM::countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset*, DataWarehouse*,
                                          DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing,
            "countMaterialPointsPerLoadCurve");
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
        int numMPMMatls = d_sharedState->getNumMPMMatls();
        int numPts = 0;
        for (int m = 0; m < numMPMMatls; m++) {
          MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);

          for (auto particle : *pset) {
            if (pLoadCurveID[particle] == (nofPressureBCs))
              ++numPts;
          }
        } // matl loop
        new_dw->put(sumlong_vartype(numPts),
                    d_labels->materialPointsPerLoadCurveLabel, 0,
                    nofPressureBCs - 1);
      } // patch loop
    }
  }
}

/*!----------------------------------------------------------------------
 * initializePressureBC
 *-----------------------------------------------------------------------*/
void
UofU_MPM::initializePressureBC(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*, DataWarehouse*,
                               DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;
  printTask(patches, patches->get(0), cout_doing, "Doing initializePressureBC");
  if (cout_dbg.active())
    cout_dbg << "Current Time (Initialize Pressure BC) = " << time << endl;

  // Calculate the force vector at each particle
  int pressureBCId = 0;
  for (int ii = 0; ii < (int)MPMPhysicalBCFactory::mpmPhysicalBCs.size();
       ii++) {
    string bcs_type = MPMPhysicalBCFactory::mpmPhysicalBCs[ii]->getType();

    if (bcs_type != "Pressure")
      return;

    // Get the material points per load curve
    sumlong_vartype numPart = 0;
    new_dw->get(numPart, d_labels->materialPointsPerLoadCurveLabel, 0,
                pressureBCId++);

    // Save the material points per load curve in the PressureBC object
    PressureBC* pbc =
      dynamic_cast<PressureBC*>(MPMPhysicalBCFactory::mpmPhysicalBCs[ii].get());
    pbc->numMaterialPoints(numPart);

    if (cout_dbg.active())
      cout_dbg << "    Load Curve = " << pressureBCId
               << " Num Particles = " << numPart << endl;

    // Calculate the force per particle at t = 0.0
    double forcePerPart = pbc->forcePerParticle(time);

    // Loop through the patches and calculate the force vector
    // at each particle
    for (int p = 0; p < patches->size(); p++) {
      const Patch* patch = patches->get(p);
      int numMPMMatls = d_sharedState->getNumMPMMatls();
      for (int m = 0; m < numMPMMatls; m++) {
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
        int matID = mpm_matl->getDWIndex();

        ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

        constParticleVariable<Point> pX;
        constParticleVariable<Matrix3> pSize;
        constParticleVariable<Matrix3> pDefGrad;
        new_dw->get(pX, d_labels->pXLabel, pset);
        new_dw->get(pSize, d_labels->pSizeLabel, pset);
        new_dw->get(pDefGrad, d_labels->pDefGradLabel, pset);

        constParticleVariable<int> pLoadCurveID;
        new_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);
        ParticleVariable<Vector> pExternalForce;
        new_dw->getModifiable(pExternalForce, d_labels->pExternalForceLabel,
                              pset);

        constParticleVariable<Vector> pDisp;
        new_dw->get(pDisp, d_labels->pDispLabel, pset);

        ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
          pExternalForceCorner3, pExternalForceCorner4;
        if (d_flags->d_useCBDI) {
          if (ii == 0) {
            new_dw->allocateAndPut(pExternalForceCorner1,
                                   d_labels->pExternalForceCorner1Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner2,
                                   d_labels->pExternalForceCorner2Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner3,
                                   d_labels->pExternalForceCorner3Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner4,
                                   d_labels->pExternalForceCorner4Label, pset);
          } else {
            new_dw->getModifiable(pExternalForceCorner1,
                                  d_labels->pExternalForceCorner1Label, pset);
            new_dw->getModifiable(pExternalForceCorner2,
                                  d_labels->pExternalForceCorner2Label, pset);
            new_dw->getModifiable(pExternalForceCorner3,
                                  d_labels->pExternalForceCorner3Label, pset);
            new_dw->getModifiable(pExternalForceCorner4,
                                  d_labels->pExternalForceCorner4Label, pset);
          }
        }

        for (auto particle : *pset) {
          if (pLoadCurveID[particle] == pressureBCId) {
            if (d_flags->d_useCBDI) {
              Vector dxCell = patch->dCell();
              pExternalForce[particle] = pbc->getForceVectorCBDI(
                pX[particle], pDisp[particle], pSize[particle],
                pDefGrad[particle], forcePerPart, time,
                pExternalForceCorner1[particle],
                pExternalForceCorner2[particle],
                pExternalForceCorner3[particle],
                pExternalForceCorner4[particle], dxCell);
            } else {
              pExternalForce[particle] =
                pbc->getForceVector(pX[particle], pDisp[particle], forcePerPart,
                                    time, pDefGrad[particle]);
            }
          }
        }
      } // matl loop
    }   // patch loop
  }     // PhysicalBC loop
}

/*!----------------------------------------------------------------------
 * scheduleComputeStableTimsetp
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeStableTimestep(const LevelP& level, SchedulerP& sched)
{
  // Nothing to do here - delt is computed as a by-product of the
  // constitutive model
  // However, this task needs to do something in the case that MPM
  // is being run on more than one level.
  Task* t = 0;
  cout_doing << UintahParallelComponent::d_myworld->myrank()
             << " MPM::scheduleComputeStableTimestep \t\t\t\tL-"
             << level->getIndex() << endl;

  t = scinew Task("UofU_MPM::actuallyComputeStableTimestep", this,
                  &UofU_MPM::actuallyComputeStableTimestep);

  const MaterialSet* mpm_matls = d_sharedState->allMPMMaterials();

  t->computes(d_sharedState->get_delt_label(), level.get_rep());
  sched->addTask(t, level->eachPatch(), mpm_matls);
}

/*!----------------------------------------------------------------------
 * actuallyComputeStableTimestep
 *-----------------------------------------------------------------------*/
void
UofU_MPM::actuallyComputeStableTimestep(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  // Put something here to satisfy the need for a reduction operation in
  // the case that there are multiple levels present
  const Level* level = getLevel(patches);
  new_dw->put(delt_vartype(999.0), d_labels->delTLabel, level);
}

/*!----------------------------------------------------------------------
 * scheduleTimeAdvance
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  MALLOC_TRACE_TAG_SCOPE("UofU_MPM::scheduleTimeAdvance()");
  if (!d_flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()))
    return;

  const PatchSet* patches = level->eachPatch();
  const MaterialSet* matls = d_sharedState->allMPMMaterials();

  // Compute body forces first
  scheduleComputeParticleBodyForce(sched, patches, matls);

  scheduleApplyExternalLoads(sched, patches, matls);
  scheduleInterpolateParticlesToGrid(sched, patches, matls);

  scheduleContactMomentumExchange(sched, patches, matls, d_labels->gVelocityLabel);
  scheduleComputeContactArea(sched, patches, matls);

  scheduleComputeInternalForce(sched, patches, matls);
  scheduleComputeAcceleration(sched, patches, matls);

  scheduleSetGridBoundaryConditions(sched, patches, matls);
  scheduleSetPrescribedMotion(sched, patches, matls);

  //scheduleCheckGridVelocity(sched, patches, matls);

  scheduleComputeVelocityAndDeformationGradient(sched, patches, matls);
  scheduleUnrotateStressAndDeformationRate(sched, patches, matls);
  scheduleComputeStressTensor(sched, patches, matls);
  scheduleRotateStress(sched, patches, matls);

  scheduleComputeBasicDamage(sched, patches, matls);
  scheduleUpdateErosionParameter(sched, patches, matls);
  
  scheduleFindRogueParticles(sched, patches, matls);
  if (d_flags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }

  scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);

  if (d_flags->d_computeScaleFactor) {
    scheduleComputeParticleScaleFactor(sched, patches, matls);
  }

  sched->scheduleParticleRelocation(
    level, d_labels->pXLabel_preReloc, d_sharedState->d_particleState_preReloc,
    d_labels->pXLabel, d_sharedState->d_particleState,
    d_labels->pParticleIDLabel, matls, 1);
}

/*!====================================================================================
 * Method: scheduleComputeParticleBodyForce
 * Purpose: Schedule a task to compute particle body forces
 * Inputs:  p.x
 * Outputs: p.bodyForce
 *====================================================================================*/
void
UofU_MPM::scheduleComputeParticleBodyForce(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleBodyForce");

  Task* t = scinew Task("UofU_MPM::computeParticleBodyForce", this,
                        &UofU_MPM::computeParticleBodyForce);

  t->requires(Task::OldDW, d_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pVelocityLabel, Ghost::None);
  t->computes(d_labels->pBodyForceAccLabel_preReloc);
  t->computes(d_labels->pCoriolisImportanceLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*!====================================================================================
 * Method: computeParticleBodyForce
 * Purpose: Actually compute particle body forces
 * Inputs:  p.x
 * Outputs: p.bodyForce
 *====================================================================================*/
void
UofU_MPM::computeParticleBodyForce(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*, DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  // Get the MPM d_flags and make local copies
  Uintah::Point rotation_center = d_flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis = d_flags->d_coordRotationAxis;
  double rotation_speed = d_flags->d_coordRotationSpeed;
  //Uintah::Point body_ref_point = d_flags->d_coord_rotation_body_ref_point;

  // Compute angular velocity vector (omega)
  Uintah::Vector omega = rotation_axis * rotation_speed;

  // Loop thru patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleBodyForce");

    // Loop thru materials
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      // Get the material ID
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      //std::cout << "Material " << m << " Starting num particles = " << pset->numParticles() << "\n";

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      // new_dw->allocateAndPut(pBodyForceAcc, d_labels->pBodyForceAccLabel,
      // pset);
      new_dw->allocateAndPut(pBodyForceAcc,
                             d_labels->pBodyForceAccLabel_preReloc, pset);

      // Create space for particle coriolis importance
      ParticleVariable<double> pCoriolisImportance;
      // new_dw->allocateAndPut(pCoriolisImportance,
      // d_labels->pCoriolisImportanceLabel, pset);
      new_dw->allocateAndPut(pCoriolisImportance,
                             d_labels->pCoriolisImportanceLabel_preReloc, pset);

      // Don't do much if coord rotation is off
      if (!d_flags->d_useCoordRotation) {
        // Iterate over the particles
        for (auto particle : *pset) {
          // Compute the body force acceleration (g)
          pBodyForceAcc[particle] = d_flags->d_gravity;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[particle] = 0.0;
        } // particle loop

      } else { // Use coordinate rotation

        // Get the particle data
        constParticleVariable<Point> pPosition;
        old_dw->get(pPosition, d_labels->pXLabel, pset);

        constParticleVariable<Vector> pVelocity;
        old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);

        // Iterate over the particles
        // std::cout << "Mat id = " << matID << " patch = " << patch <<
        // std::endl;
        // std::cout << "Particle subset = " << *pset;
        // std::cout << "Num particles = " << pset->numParticles() << std::endl;
        for (auto particle : *pset) {
          // std::cout << " Particle # = " << pidx << std::endl;
          // Compute the local "x" vector wrt ref point in body
          // Vector xVec = pPosition[pidx].vector() - body_ref_point;

          // Compute reference vector R wrt rotation center
          // Uintah::Vector Rvec = body_ref_point - rotation_center;

          // Compute the local "r" vector with respect to rotation center
          // Vector rVec = Rvec + pPosition[pidx].vector();

          // Compute the Coriolis term (omega x v)
          Vector coriolis_accel =
            Uintah::Cross(omega, pVelocity[particle]) * 2.0;

          // Compute the centrifugal term (omega x omega x r)
          // Simplified version where body ref point is not needed
          Vector rVec = pPosition[particle] - rotation_center;
          Vector omega_x_r = Uintah::Cross(omega, rVec);
          Vector centrifugal_accel = Uintah::Cross(omega, omega_x_r);

          // Compute the body force acceleration (g - omega x omega x r - 2
          // omega x v)
          pBodyForceAcc[particle] =
            d_flags->d_gravity - centrifugal_accel - coriolis_accel;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[particle] =
            coriolis_accel.length() /
            (centrifugal_accel.length() + coriolis_accel.length());

          //std::cout << "After compute body force: material = " << m
          //          << " particle = " << particle
          //          << " pVel = " << pVelocity[particle] << "\n";

        } // particle loop
      }   // end if coordinate rotation
    }     // matl loop
  }       // patch loop
}

/*====================================================================================*/
// Apply external loads
//*
//* applyExternalLoads
//*   in(p.externalForce)
//*   out(p.externalForceNew) */
/*====================================================================================*/
void
UofU_MPM::scheduleApplyExternalLoads(SchedulerP& sched, const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "MPM::scheduleApplyExternalLoads");

  Task* t =
    scinew Task("UofU_MPM::applyExternalLoads", this, &UofU_MPM::applyExternalLoads);

  t->requires(Task::OldDW, d_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pMassLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pDispLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pDefGradLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pExternalForceLabel, Ghost::None);
  t->computes(d_labels->pExtForceLabel_preReloc);
  if (d_flags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_labels->pLoadCurveIDLabel, Ghost::None);
    t->computes(d_labels->pLoadCurveIDLabel_preReloc);
    if (d_flags->d_useCBDI) {
      t->computes(d_labels->pExternalForceCorner1Label);
      t->computes(d_labels->pExternalForceCorner2Label);
      t->computes(d_labels->pExternalForceCorner3Label);
      t->computes(d_labels->pExternalForceCorner4Label);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * addExternalLoads
 *-----------------------------------------------------------------------*/
void
UofU_MPM::applyExternalLoads(const ProcessorGroup*, const PatchSubset* patches,
                             const MaterialSubset*, DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  // Get the current time
  double time = d_sharedState->getElapsedTime();

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << endl;
  }

  // Calculate the force vector at each particle for each pressure bc
  std::vector<double> forcePerPart;
  std::vector<PressureBC*> pbcP;
  std::vector<MomentBC*> pbcM;
  if (d_flags->d_useLoadCurves) {
    for (auto physicalBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
      string bcs_type = physicalBC->getType();
      if (bcs_type != "Pressure")
        continue;
      PressureBC* pbc = dynamic_cast<PressureBC*>(physicalBC.get());
      pbcP.push_back(pbc);

      // Calculate the force per particle at current time
      forcePerPart.push_back(pbc->forcePerParticle(time));
    }
  }

  // Loop thru patches to update external force vector
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing applyExternalLoads");

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Get the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      ParticleVariable<Vector> pExternalForce_new;

      old_dw->get(pX, d_labels->pXLabel, pset);
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad, d_labels->pDefGradLabel, pset);
      new_dw->allocateAndPut(pExternalForce_new,
                             d_labels->pExtForceLabel_preReloc, pset);

      if (d_flags->d_useLoadCurves) {
        // Get the load curve data
        // Recycle the loadCurveIDs
        constParticleVariable<int> pLoadCurveID;
        old_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);

        ParticleVariable<int> pLoadCurveID_new;
        new_dw->allocateAndPut(pLoadCurveID_new,
                               d_labels->pLoadCurveIDLabel_preReloc, pset);
        pLoadCurveID_new.copyData(pLoadCurveID);
        // std::cout << " Recycled load curve ID" << std::endl;

        // Check whether it's a presure or moment bc
        bool do_PressureBCs = false;
        for (auto physicalBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          string bcs_type = physicalBC->getType();
          if (bcs_type != "Pressure")
            continue;
          do_PressureBCs = true;
        }

        if (!do_PressureBCs) {
          for (auto particle : *pset) {
            pExternalForce_new[particle] = 0.;
          }

        } else {
          // Get the external force data and allocate new space for
          // external force
          constParticleVariable<Vector> pDisp;
          old_dw->get(pDisp, d_labels->pDispLabel, pset);

          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_labels->pExternalForceLabel, pset);

          ParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
            pExternalForceCorner3, pExternalForceCorner4;
          if (d_flags->d_useCBDI) {
            new_dw->allocateAndPut(pExternalForceCorner1,
                                   d_labels->pExternalForceCorner1Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner2,
                                   d_labels->pExternalForceCorner2Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner3,
                                   d_labels->pExternalForceCorner3Label, pset);
            new_dw->allocateAndPut(pExternalForceCorner4,
                                   d_labels->pExternalForceCorner4Label, pset);
          }

          // Iterate over the particles
          for (auto particle : *pset) {
            int loadCurveID = pLoadCurveID[particle] - 1;
            if (loadCurveID < 0) {
              pExternalForce_new[particle] = Vector(0.0, 0.0, 0.0);
            } else {
              PressureBC* pbc = pbcP[loadCurveID];
              double force = forcePerPart[loadCurveID];

              if (d_flags->d_useCBDI) {
                Vector dxCell = patch->dCell();
                pExternalForce_new[particle] = pbc->getForceVectorCBDI(
                  pX[particle], pDisp[particle], pSize[particle],
                  pDefGrad[particle], force, time,
                  pExternalForceCorner1[particle],
                  pExternalForceCorner2[particle],
                  pExternalForceCorner3[particle],
                  pExternalForceCorner4[particle], dxCell);
              } else {
                pExternalForce_new[particle] =
                  pbc->getForceVector(pX[particle], pDisp[particle], force,
                                      time, pDefGrad[particle]);
              }
            } // end if (loadCurveID < 0)
          } // end particle loop
        } // end if (doPressureBCs)

        // MMS (compute body force)
        if (!d_flags->d_mmsType.empty()) {

          MMS MMSObject;
          MMSObject.computeBodyForceForMMS(old_dw, new_dw, time, pset, d_labels, d_flags,
                                           pExternalForce_new);
          /*
          for (auto particle : *pset) {
            if (!std::isfinite(pExternalForce_new[particle].length())) {
              std::cout << "MMS: pExt = " << pExternalForce_new[particle] << "\n";
            }
          }
          */
        } 

      } else { // if (!d_useLoadCurves)

        // MMS
        if (!d_flags->d_mmsType.empty()) {

          MMS MMSObject;
          MMSObject.computeExternalForceForMMS(old_dw, new_dw, time, pset, d_labels, d_flags,
                                               pExternalForce_new);
        } else {

          // Get the external force data and allocate new space for
          // external force and copy the data
          constParticleVariable<Vector> pExternalForce;
          old_dw->get(pExternalForce, d_labels->pExternalForceLabel, pset);

          for (auto particle : *pset) {
            pExternalForce_new[particle] = pExternalForce[particle] * d_flags->d_forceIncrementFactor;
          }
        }
      } // end if (d_useLoadCurves)

    } // end matl loop
  } // end patch loop
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateParticlesToGrid
 * interpolateParticlesToGrid
 *   in(P.MASS_ip_av, P.VELOCITY, P.NAT_X)
 *   operation(interpolate the P.MASS and P.VEL to the grid
 *             using P.NAT_X and some shape function evaluations)
 *   out(G.MASS_ip_av, G.VELOCITY)
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM::scheduleInterpolateParticlesToGrid");

  Task* t = scinew Task("UofU_MPM::interpolateParticlesToGrid", this,
                        &UofU_MPM::interpolateParticlesToGrid);
  t->requires(Task::OldDW, d_labels->pXLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pMassLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pVolumeLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pVelocityLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW, d_labels->pBodyForceAccLabel_preReloc,
              Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Task::NewDW, d_labels->pExtForceLabel_preReloc,
              Ghost::AroundNodes, d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pDefGradLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  if (d_flags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_labels->pLoadCurveIDLabel, Ghost::AroundNodes,
                d_numGhostParticles);
    if (d_flags->d_useCBDI) {
      t->requires(Task::NewDW, d_labels->pExternalForceCorner1Label,
                  Ghost::AroundNodes, d_numGhostParticles);
      t->requires(Task::NewDW, d_labels->pExternalForceCorner2Label,
                  Ghost::AroundNodes, d_numGhostParticles);
      t->requires(Task::NewDW, d_labels->pExternalForceCorner3Label,
                  Ghost::AroundNodes, d_numGhostParticles);
      t->requires(Task::NewDW, d_labels->pExternalForceCorner4Label,
                  Ghost::AroundNodes, d_numGhostParticles);
    }
  }

#ifdef DEBUG_WITH_PARTICLE_ID
  t->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::AroundNodes,
              d_numGhostParticles);
#endif

  t->computes(d_labels->gMassLabel, d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);
  t->computes(d_labels->gVolumeLabel, d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);
  t->computes(d_labels->gVelocityLabel, d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);

  t->computes(d_labels->gMassLabel);
  t->computes(d_labels->gVolumeLabel);
  t->computes(d_labels->gVelocityLabel);
  t->computes(d_labels->gBodyForceLabel);
  t->computes(d_labels->gExternalForceLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * interpolateParticlesToGrid
 *-----------------------------------------------------------------------*/
void
UofU_MPM::interpolateParticlesToGrid(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing interpolateParticlesToGrid");

    int numMatls = d_sharedState->getNumMPMMatls();
    auto interpolator = d_flags->d_interpolator->clone(patch);
    auto linear_interpolator = std::make_unique<LinearInterpolator>(patch);
    auto num_influence_nodes = interpolator->size();
    auto num_linear_influence_nodes = linear_interpolator->size();

    vector<IntVector> influenceNodes(num_influence_nodes);
    vector<double> S_ip_av(num_influence_nodes);
    string interp_type = d_flags->d_interpolatorType;

    NCVariable<double> gMassglobal, gVolumeglobal;
    NCVariable<Vector> gVelglobal;
    new_dw->allocateAndPut(gMassglobal, d_labels->gMassLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    new_dw->allocateAndPut(gVolumeglobal, d_labels->gVolumeLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    new_dw->allocateAndPut(gVelglobal, d_labels->gVelocityLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gMassglobal.initialize(SMALL_NUM_MPM);
    gVolumeglobal.initialize(SMALL_NUM_MPM);
    gVelglobal.initialize(Vector(0.0));

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume;
      constParticleVariable<Vector> pVelocity, pBodyForceAcc, pExternalForce;
      constParticleVariable<Point> pExternalForceCorner1, pExternalForceCorner2,
        pExternalForceCorner3, pExternalForceCorner4;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad_old;

      ParticleSubset* pset =
        old_dw->getParticleSubset(matID, patch, Ghost::AroundNodes,
                                  d_numGhostParticles, d_labels->pXLabel);

      old_dw->get(pX, d_labels->pXLabel, pset);
      old_dw->get(pMass, d_labels->pMassLabel, pset);
      old_dw->get(pVolume, d_labels->pVolumeLabel, pset);
      old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_labels->pDefGradLabel, pset);
      new_dw->get(pBodyForceAcc, d_labels->pBodyForceAccLabel_preReloc, pset);
      new_dw->get(pExternalForce, d_labels->pExtForceLabel_preReloc, pset);
      constParticleVariable<int> pLoadCurveID;
      if (d_flags->d_useLoadCurves) {
        old_dw->get(pLoadCurveID, d_labels->pLoadCurveIDLabel, pset);
        if (d_flags->d_useCBDI) {
          new_dw->get(pExternalForceCorner1,
                      d_labels->pExternalForceCorner1Label, pset);
          new_dw->get(pExternalForceCorner2,
                      d_labels->pExternalForceCorner2Label, pset);
          new_dw->get(pExternalForceCorner3,
                      d_labels->pExternalForceCorner3Label, pset);
          new_dw->get(pExternalForceCorner4,
                      d_labels->pExternalForceCorner4Label, pset);
        }
      }

#ifdef DEBUG_WITH_PARTICLE_ID
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);
#endif

      // Create arrays for the grid data
      NCVariable<double> gMass;
      NCVariable<double> gVolume;
      NCVariable<Vector> gVelocity;
      NCVariable<Vector> gBodyForce;
      NCVariable<Vector> gExternalForce;

      new_dw->allocateAndPut(gMass, d_labels->gMassLabel, matID, patch);
      new_dw->allocateAndPut(gVolume, d_labels->gVolumeLabel, matID, patch);
      new_dw->allocateAndPut(gVelocity, d_labels->gVelocityLabel, matID, patch);
      new_dw->allocateAndPut(gBodyForce, d_labels->gBodyForceLabel, matID,
                             patch);
      new_dw->allocateAndPut(gExternalForce, d_labels->gExternalForceLabel,
                             matID, patch);

      gMass.initialize(SMALL_NUM_MPM);
      gVolume.initialize(SMALL_NUM_MPM);
      gVelocity.initialize(Vector(0, 0, 0));
      gBodyForce.initialize(Vector(0, 0, 0));
      gExternalForce.initialize(Vector(0, 0, 0));

      // Interpolate particle data to grid
      Vector total_mom(0.0, 0.0, 0.0);
      for (auto particle : *pset) {
        interpolator->findCellAndWeights(pX[particle], influenceNodes, S_ip_av,
                                         pSize[particle],
                                         pDefGrad_old[particle]);
        auto pMom = pVelocity[particle] * pMass[particle];
        total_mom += pMom;

        //std::cout << "In interpolateToGrid:  material = " << m << " particle = " << particle
        //          << " pVel = " << pVelocity[particle]
        //          << " pMass = " << pMass[particle] << "\n";

        // Add each particles contribution to the local mass & velocity
        // Must use the node indices
        for (int k = 0; k < num_influence_nodes; k++) {
          auto node = influenceNodes[k];
          if (patch->containsNode(node)) {
            gMass[node] += pMass[particle] * S_ip_av[k];
            gVelocity[node] += pMom * S_ip_av[k];
            //std::cout << "particle = " << particle 
            //          << " node = " << node << " k = " << k
            //          << " S_ip = " << S_ip_av[k]
            //          << " pMom = " << pMom 
            //          << " gMass = " << gMass[node] 
            //          << " gVel = " << gVelocity[node] << "\n";
            gVolume[node] += pVolume[particle] * S_ip_av[k];
            if (!d_flags->d_useCBDI) {
              gExternalForce[node] += pExternalForce[particle] * S_ip_av[k];
            }
            gBodyForce[node] +=
              pBodyForceAcc[particle] * pMass[particle] * S_ip_av[k];
          }
        }
        if (d_flags->d_useLoadCurves && d_flags->d_useCBDI) {
          vector<IntVector> influenceNodesCorner1(num_linear_influence_nodes);
          vector<IntVector> influenceNodesCorner2(num_linear_influence_nodes);
          vector<IntVector> influenceNodesCorner3(num_linear_influence_nodes);
          vector<IntVector> influenceNodesCorner4(num_linear_influence_nodes);
          vector<double> S_ip_av_Corner1(num_linear_influence_nodes);
          vector<double> S_ip_av_Corner2(num_linear_influence_nodes);
          vector<double> S_ip_av_Corner3(num_linear_influence_nodes);
          vector<double> S_ip_av_Corner4(num_linear_influence_nodes);
          linear_interpolator->findCellAndWeights(
            pExternalForceCorner1[particle], influenceNodesCorner1,
            S_ip_av_Corner1, pSize[particle], pDefGrad_old[particle]);
          linear_interpolator->findCellAndWeights(
            pExternalForceCorner2[particle], influenceNodesCorner2,
            S_ip_av_Corner2, pSize[particle], pDefGrad_old[particle]);
          linear_interpolator->findCellAndWeights(
            pExternalForceCorner3[particle], influenceNodesCorner3,
            S_ip_av_Corner3, pSize[particle], pDefGrad_old[particle]);
          linear_interpolator->findCellAndWeights(
            pExternalForceCorner4[particle], influenceNodesCorner4,
            S_ip_av_Corner4, pSize[particle], pDefGrad_old[particle]);
          for (int k = 0; k < num_linear_influence_nodes; k++) {
            auto node = influenceNodesCorner1[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] +=
                pExternalForce[particle] * S_ip_av_Corner1[k];
            }
            node = influenceNodesCorner2[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] +=
                pExternalForce[particle] * S_ip_av_Corner2[k];
            }
            node = influenceNodesCorner3[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] +=
                pExternalForce[particle] * S_ip_av_Corner3[k];
            }
            node = influenceNodesCorner4[k];
            if (patch->containsNode(node)) {
              gExternalForce[node] +=
                pExternalForce[particle] * S_ip_av_Corner4[k];
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
        //std::cout << " patch = " << patch << " node = " << c 
        //          << " gVel = " << gVelocity[c] << "\n";
      }

      // Apply velocity boundary conditions (if symmetry)
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch, matID, "Velocity", gVelocity, interp_type);
      bc.setBoundaryCondition(patch, matID, "Symmetric", gVelocity, interp_type);

      //for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
      //  std::cout << " after BC: node = " << *iter 
      //            << " gVel = " << gVelocity[*iter] << "\n";
      //}

    } // End loop over materials

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gVelglobal[c] /= gMassglobal[c];
    }
  } // End loop over patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeContactArea
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeContactArea(SchedulerP& sched, const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  /** computeContactArea */
  if (d_boundaryTractionFaces.size() > 0) {
    printSchedule(patches, cout_doing, "MPM::scheduleComputeContactArea");
    Task* t = scinew Task("UofU_MPM::computeContactArea", this,
                          &UofU_MPM::computeContactArea);

    t->requires(Task::NewDW, d_labels->gVolumeLabel, Ghost::None);
    for (auto face : d_boundaryTractionFaces) {
      int iface = (int)(face);
      t->computes(d_labels->BndyContactCellAreaLabel[iface]);
    }
    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * computeContactArea
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeContactArea(const ProcessorGroup*, const PatchSubset* patches,
                             const MaterialSubset*, DataWarehouse* /*old_dw*/,
                             DataWarehouse* new_dw)
{
  // six indices for each of the faces
  double bndyCArea[6] = { 0, 0, 0, 0, 0, 0 };

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeContactArea");

    Vector dx = patch->dCell();

    int numMPMMatls = d_sharedState->getNumMPMMatls();

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      constNCVariable<double> gVolume;

      new_dw->get(gVolume, d_labels->gVolumeLabel, matID, patch, Ghost::None, 0);

      for (auto face : d_boundaryTractionFaces) {
        int iface = (int)(face);

        // Check if the face is on an external boundary
        if (patch->getBCType(face) == Patch::Neighbor)
          continue;

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,

        // loop over face nodes to find boundary areas
        // Because this calculation uses gVolume, particle volumes interpolated
        // to
        // the nodes, it will give 1/2 the expected value because the particle
        // values
        // are distributed to all nodes, not just those on this face.  It would
        // require
        // particles on the other side of the face to "fill" the nodal volumes
        // and give
        // the correct area when divided by the face normal cell dimension
        // (celldepth).
        // To correct for this, nodearea incorporates a factor of two.

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
    int iface = (int)(face);
    new_dw->put(sum_vartype(bndyCArea[iface]),
                d_labels->BndyContactCellAreaLabel[iface]);
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeInternalForce
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeInternalForce(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM::scheduleComputeInternalForce");

  Task* t = scinew Task("UofU_MPM::computeInternalForce", this,
                        &UofU_MPM::computeInternalForce);

  t->requires(Task::NewDW, d_labels->gVolumeLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gVolumeLabel,
              d_sharedState->getAllInOneMatl(), Task::OutOfDomain, Ghost::None);
  t->requires(Task::OldDW, d_labels->pStressLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pVolumeLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pXLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW, d_labels->pDefGradLabel, Ghost::AroundNodes,
              d_numGhostParticles);
#ifdef DEBUG_WITH_PARTICLE_ID
  t->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::AroundNodes,
              d_numGhostParticles);
#endif

  if (d_flags->d_artificialViscosity) {
    t->requires(Task::OldDW, d_labels->p_qLabel, Ghost::AroundNodes,
                d_numGhostParticles);
  }

  t->computes(d_labels->gInternalForceLabel);

  for (auto face : d_boundaryTractionFaces) {
    int iface = static_cast<int>(face);
    t->requires(Task::NewDW, d_labels->BndyContactCellAreaLabel[iface]);
    t->computes(d_labels->BndyForceLabel[iface]);
    t->computes(d_labels->BndyContactAreaLabel[iface]);
    t->computes(d_labels->BndyTractionLabel[iface]);
  }

  t->computes(d_labels->gStressForSavingLabel);
  t->computes(d_labels->gStressForSavingLabel, d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeInternalForce
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeInternalForce(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*, DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  // node based forces
  Vector bndyForce[6];
  Vector bndyTraction[6];
  for (int iface = 0; iface < 6; iface++) {
    bndyForce[iface] = Vector(0.);
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

    auto interpolator = d_flags->d_interpolator->clone(patch);
    auto numInfluenceNodes = interpolator->size();
    vector<IntVector> influenceNodes(numInfluenceNodes);
    vector<double> S_ip_av(numInfluenceNodes);
    vector<Vector> d_S_ip_av(numInfluenceNodes);
    string interp_type = d_flags->d_interpolatorType;

    int numMPMMatls = d_sharedState->getNumMPMMatls();

    NCVariable<Matrix3> gStressglobal;
    constNCVariable<double> gVolumeglobal;
    new_dw->get(gVolumeglobal, d_labels->gVolumeLabel,
                d_sharedState->getAllInOneMatl()->get(0), patch, Ghost::None,
                0);
    new_dw->allocateAndPut(gStressglobal, d_labels->gStressForSavingLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      constParticleVariable<Point> pX;
      constParticleVariable<double> pVol, p_q;
      constParticleVariable<Matrix3> pStress, pSize, pDefGrad_old;
      constNCVariable<double> gVolume;
      NCVariable<Vector> gInternalForce;
      NCVariable<Matrix3> gStress;

      ParticleSubset* pset = old_dw->getParticleSubset(
        matID, patch, Ghost::AroundNodes, d_numGhostParticles, d_labels->pXLabel);

      old_dw->get(pX, d_labels->pXLabel, pset);
      old_dw->get(pVol, d_labels->pVolumeLabel, pset);
      old_dw->get(pStress, d_labels->pStressLabel, pset);
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      old_dw->get(pDefGrad_old, d_labels->pDefGradLabel, pset);

#ifdef DEBUG_WITH_PARTICLE_ID
      constParticleVariable<long64> pParticleID;
      old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);
#endif

      new_dw->get(gVolume, d_labels->gVolumeLabel, matID, patch, Ghost::None, 0);

      new_dw->allocateAndPut(gStress, d_labels->gStressForSavingLabel, matID,
                             patch);
      new_dw->allocateAndPut(gInternalForce, d_labels->gInternalForceLabel, matID,
                             patch);

      if (d_flags->d_artificialViscosity) {
        old_dw->get(p_q, d_labels->p_qLabel, pset);
      } else {
        ParticleVariable<double> p_q_create;
        new_dw->allocateTemporary(p_q_create, pset);
        for (auto particle : *pset) {
          p_q_create[particle] = 0.0;
        }
        p_q = p_q_create; // reference created data
      }

      gInternalForce.initialize(Vector(0, 0, 0));

      Matrix3 stressvol;
      Matrix3 stresspress;

      // for the non axisymmetric case:
      if (!d_flags->d_axisymmetric) {

        for (auto particle : *pset) {
          // Get the node indices that surround the cell
          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[particle], influenceNodes, S_ip_av, d_S_ip_av, pSize[particle],
            pDefGrad_old[particle]);
          stressvol = pStress[particle] * pVol[particle];
          stresspress = pStress[particle] - Identity * p_q[particle];
          // cerr << " particle = " << particle << " pStress = " << pStress[particle] << endl;

          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = influenceNodes[k];
            if (patch->containsNode(node)) {
              Vector div(d_S_ip_av[k].x() * oodx[0], d_S_ip_av[k].y() * oodx[1],
                         d_S_ip_av[k].z() * oodx[2]);
              gInternalForce[node] -= (div * stresspress) * pVol[particle];
              gStress[node] += stressvol * S_ip_av[k];

#ifdef CHECK_ISFINITE
              if (!std::isfinite(gInternalForce[node].x()) ||
                  !std::isfinite(gInternalForce[node].y()) ||
                  !std::isfinite(gInternalForce[node].z())) {
                std::cout << "vol = " << pVol[particle] << " node = " << node
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
      if (d_flags->d_axisymmetric) {
        for (auto particle : *pset) {
          interpolator->findCellAndWeightsAndShapeDerivatives(
            pX[particle], influenceNodes, S_ip_av, d_S_ip_av, pSize[particle],
            pDefGrad_old[particle]);

          stressvol = pStress[particle] * pVol[particle];
          stresspress = pStress[particle] - Identity * p_q[particle];

#ifdef CHECK_ISFINITE
          if (!std::isfinite(stresspress(0, 0)) ||
              !std::isfinite(stresspress(0, 1)) ||
              !std::isfinite(stresspress(0, 2)) ||
              !std::isfinite(stresspress(1, 1)) ||
              !std::isfinite(stresspress(1, 2)) ||
              !std::isfinite(stresspress(2, 2))) {
            std::cout << " p_pressure = " << p_pressure[particle]
                      << " p_q = " << p_q[particle] << "\n";
          }
#endif

          // r is the x direction, z (axial) is the y direction
          double IFr = 0., IFz = 0.;
          for (int k = 0; k < numInfluenceNodes; k++) {
            auto node = influenceNodes[k];
            if (patch->containsNode(node)) {
              IFr = d_S_ip_av[k].x() * oodx[0] * stresspress(0, 0) +
                    d_S_ip_av[k].y() * oodx[1] * stresspress(0, 1) +
                    d_S_ip_av[k].z() * stresspress(2, 2);
              IFz = d_S_ip_av[k].x() * oodx[0] * stresspress(0, 1) +
                    d_S_ip_av[k].y() * oodx[1] * stresspress(1, 1);
              gInternalForce[node] -= Vector(IFr, IFz, 0.0) * pVol[particle];
              gStress[node] += stressvol * S_ip_av[k];
#ifdef CHECK_ISFINITE
              if (!std::isfinite(gInternalForce[node].x()) ||
                  !std::isfinite(gInternalForce[node].y()) ||
                  !std::isfinite(gInternalForce[node].z())) {
                std::cout << "vol = " << pVol[particle]
                          << " node = " << influenceNodes[k]
                          << " f_i = " << gInternalForce[node]
                          << " sig_g = " << gStress[node]
                          << " sig_p = " << stresspress << " IFr = " << IFr
                          << " IFz = " << IFz << "\n";
              }
#endif
#ifdef DEBUG_WITH_PARTICLE_ID
              // if (pParticleID[part] == 158913855488) {
              if (node == IntVector(3, 38, 0)) {
                proc0cout << "Particle ID = " << pParticleID[particle]
                          << " node = " << node << " dS = " << d_S_ip_av[k]
                          << " IFr = " << IFr << " IFz = " << IFz
                          << " stress = " << pStress[particle]
                          << " damp = " << p_q[particle]
                          << " stresspress = " << stresspress
                          << " vol = " << pVol[particle]
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
        if (patch->getBCType(face) == Patch::Neighbor)
          continue;

        const int iface = (int)face;

        // We are on the boundary, i.e. not on an interior patch
        // boundary, and also on the correct side,
        IntVector projlow, projhigh;
        patch->getFaceNodes(face, 0, projlow, projhigh);
        Vector norm = face_norm(face);
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
      bc.setBoundaryCondition(patch, matID, "Symmetric", gInternalForce,
                              interp_type);
#ifdef DEBUG_WITH_PARTICLE_ID
      // IntVector node(3,38,0);
      if (patch->containsNode(node)) {
        proc0cout << "After BC: Material = " << m << " Node = " << node
                  << " fint_g = " << gInternalForce[node] << "\n";
      }
#endif

      //for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      //  std::cout << "After internal force: node = " << *iter
      //            << " gInternalForce = " << gInternalForce[*iter] << "\n";
      //}
    }

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      gStressglobal[c] /= gVolumeglobal[c];
    }
  }

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (auto face : d_boundaryTractionFaces) {
    int iface = (int)(face);
    new_dw->put(sumvec_vartype(bndyForce[iface]),
                d_labels->BndyForceLabel[iface]);

    sum_vartype bndyContactCellArea_iface;
    new_dw->get(bndyContactCellArea_iface,
                d_labels->BndyContactCellAreaLabel[iface]);

    if (bndyContactCellArea_iface > 0)
      bndyTraction[iface] /= bndyContactCellArea_iface;

    new_dw->put(sumvec_vartype(bndyTraction[iface]),
                d_labels->BndyTractionLabel[iface]);

    // Use the face force and traction calculations to provide a second estimate
    // of the contact area.
    double bndyContactArea_iface = bndyContactCellArea_iface;
    if (bndyTraction[iface][iface / 2] * bndyTraction[iface][iface / 2] >
        1.e-12)
      bndyContactArea_iface =
        bndyForce[iface][iface / 2] / bndyTraction[iface][iface / 2];

    new_dw->put(sum_vartype(bndyContactArea_iface),
                d_labels->BndyContactAreaLabel[iface]);
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeAcceleration
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeAcceleration(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing,
                "UofU_MPM::scheduleComputeAcceleration");

  Task* t = scinew Task("UofU_MPM::computeAcceleration", this,
                        &UofU_MPM::computeAcceleration);

  t->requires(Task::OldDW, d_sharedState->get_delt_label());

  t->requires(Task::NewDW, d_labels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gBodyForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gExternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->gVelocityLabel, Ghost::None);

  t->computes(d_labels->gAccelerationLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeAcceleration
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeAcceleration(const ProcessorGroup* pg,
                              const PatchSubset* patches,
                              const MaterialSubset* ms,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing,
              "Doing computeAndIntegrateAcceleration");

    for (int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<Vector> gVelocity;
      constNCVariable<Vector> gInternalForce, gBodyForce, gExternalForce;
      constNCVariable<double> gMass;

      delt_vartype delT;
      old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches));

      new_dw->get(gVelocity, d_labels->gVelocityLabel, matID, patch,
                  Ghost::None, 0);
      new_dw->get(gInternalForce, d_labels->gInternalForceLabel, matID, patch,
                  Ghost::None, 0);
      new_dw->get(gBodyForce, d_labels->gBodyForceLabel, matID, patch, Ghost::None, 0);
      new_dw->get(gExternalForce, d_labels->gExternalForceLabel, matID, patch,
                  Ghost::None, 0);
      new_dw->get(gMass, d_labels->gMassLabel, matID, patch, Ghost::None, 0);

      NCVariable<Vector> gAcceleration;
      new_dw->allocateAndPut(gAcceleration, d_labels->gAccelerationLabel, matID,
                             patch);

      gAcceleration.initialize(Vector(0., 0., 0.));

      double damp_coef = d_flags->d_artificialDampCoeff;

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;

        Vector acc(0., 0., 0.);
        if (gMass[c] > d_flags->d_minMassForAcceleration) {
          acc =
            (gInternalForce[c] + gExternalForce[c] + gBodyForce[c]) / gMass[c];
          acc -= damp_coef * gVelocity[c];
        }

        gAcceleration[c] = acc;
        /*
        std::cout << "After acceleration: material = " << m << " node = " << c
                  << " gMass = " << gMass[c] << " minMass = " << d_flags->d_minMassForAcceleration << "\n"
                  << " gInt = " << gInternalForce[c] 
                  << " gExt = " << gExternalForce[c] << "\n"
                  << " gAcceleration = " << gAcceleration[c] << "\n";
        */
      }
    } // matls
  }
}

/*!----------------------------------------------------------------------
 * scheduleContactMomentumExchange
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleContactMomentumExchange(SchedulerP& sched, const PatchSet* patches,
                                          const MaterialSet* matls, 
                                          const VarLabel* gVelLabel)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM::scheduleContactMomentumExchange");
  d_contactModel->addComputesAndRequires(sched, patches, matls, gVelLabel);
}

/*!----------------------------------------------------------------------
 * scheduleSetGridBoundaryConditions
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                            const PatchSet* patches,
                                            const MaterialSet* matls)

{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;
  printSchedule(patches, cout_doing, "MPM::scheduleSetGridBoundaryConditions");
  Task* t = scinew Task("UofU_MPM::setGridBoundaryConditions", this,
                        &UofU_MPM::setGridBoundaryConditions);

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_sharedState->get_delt_label());

  t->modifies(d_labels->gVelocityLabel, mss);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * setGridBoundaryConditions
 *-----------------------------------------------------------------------*/
void
UofU_MPM::setGridBoundaryConditions(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset*,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing setGridBoundaryConditions");

    int numMPMMatls = d_sharedState->getNumMPMMatls();

    delt_vartype delT;
    old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches));

    string interp_type = d_flags->d_interpolatorType;
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity;

      new_dw->getModifiable(gVelocity, d_labels->gVelocityLabel, matID,
                            patch);

      // Apply grid boundary conditions to the velocity_star and
      // acceleration before interpolating back to the particles
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch, matID, "Velocity", gVelocity,
                              interp_type);
      bc.setBoundaryCondition(patch, matID, "Symmetric", gVelocity,
                              interp_type);

      //for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
      //     iter++) {
      //  std::cout << "After grid BC: node = " << *iter
      //            << " gVelocity = " << gVelocity[*iter] << "\n";
      //}
    }   // matl loop
  }     // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleSetPrescribedMotion
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleSetPrescribedMotion(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)

{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  if (d_flags->d_prescribeDeformation) {
    printSchedule(patches, cout_doing, "MPM::scheduleSetPrescribedMotion");

    Task* t = scinew Task("UofU_MPM::setPrescribedMotion", this,
                          &UofU_MPM::setPrescribedMotion);

    const MaterialSubset* mss = matls->getUnion();
    t->modifies(d_labels->gAccelerationLabel, mss);
    t->modifies(d_labels->gVelocityLabel, mss);
    t->requires(Task::OldDW, d_sharedState->get_delt_label());

    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * setPrescribedMotion
 *-----------------------------------------------------------------------*/
void
UofU_MPM::setPrescribedMotion(const ProcessorGroup*, const PatchSubset* patches,
                              const MaterialSubset*, DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing setPrescribedMotion");

    // Get the current time
    double time = d_sharedState->getElapsedTime();
    delt_vartype delT;
    old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches));

    int numMPMMatls = d_sharedState->getNumMPMMatls();

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matlID = mpm_matl->getDWIndex();
      NCVariable<Vector> gVelocity, gAcceleration;

      new_dw->getModifiable(gVelocity, d_labels->gVelocityLabel,
                            matlID, patch);
      new_dw->getModifiable(gAcceleration, d_labels->gAccelerationLabel, matlID,
                            patch);

      gAcceleration.initialize(Vector(0.0));

      // Get F and Q from file by interpolating between available times
      auto t_upper_iter = std::upper_bound(d_prescribedTimes.begin(),
                                           d_prescribedTimes.end(), time);

      auto t_upper_index = t_upper_iter - d_prescribedTimes.begin();
      if (t_upper_iter == d_prescribedTimes.end()) {
        t_upper_index = --t_upper_iter - d_prescribedTimes.begin();
      }

      auto t_lower = *(t_upper_iter - 1);
      auto t_upper = *t_upper_iter;
      auto ss = (time - t_lower) / (t_upper - t_lower);

      // Interpolate to get the deformation gradient at the current time:
      auto F_lower = d_prescribedF[t_upper_index - 1];
      auto F_upper = d_prescribedF[t_upper_index];
      auto Ft = (1 - ss) * F_lower + ss * F_upper;

      // Calculate the rate of the deformation gradient without the rotation:
      auto Fdot = (F_upper - F_lower) / (t_upper - t_lower);
      auto Ft_inv = Ft.Inverse();
      auto L = Fdot * Ft_inv;

      // Now we need to construct the rotation matrix and its time rate:
      // We are only interested in the rotation information at the next
      // specified time
      // since the rotations specified should be relative to the previously
      // specified time.
      // For example if I specify Theta=90 at time=1.0, and Theta = 91 and
      // time=2.0 the
      // total rotation at time=2.0 will be 181 degrees.
      const double pi = M_PI; // 3.1415926535897932384626433832795028841972;
      const double degtorad = pi / 180.0;

      auto theta_upper = d_prescribedAngle[t_upper_index];
      auto thetat = ss * theta_upper * degtorad;

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
      if (d_flags->d_exactDeformation) {
        // Check to see we do not exceed bounds
        int count = 0;
        auto t_upper_iter_copy = t_upper_iter;
        //auto t_upper_index_copy = t_upper_index;
        while (++t_upper_iter_copy != d_prescribedTimes.end()) {
          ++count;
          // t_upper_index_copy = t_upper_iter_copy - d_prescribedTimes.begin();
          // std::cout << "t_upper_index_copy = " << t_upper_index_copy << "
          // count = " << count << "\n";
          if (count > 1)
            break;
        }

        // If there are at least two extra data points
        if (count > 1) {
          double t3 = d_prescribedTimes[t_upper_index + 1];
          double t4 = d_prescribedTimes[t_upper_index + 2];
          if (time == 0 && t4 != 0) {
            new_dw->put(delt_vartype(t3 - t_upper),
                        d_sharedState->get_delt_label(), getLevel(patches));

          } else {
            F_lower = d_prescribedF[t_upper_index]; // last prescribed
            // deformation gradient
            F_upper = d_prescribedF[t_upper_index +
                                    1]; // next prescribed deformation gradient
            Ft = (1 - ss) * F_lower + ss * F_upper;
            Ft_inv = Ft.Inverse();
            Fdot = (F_upper - F_lower) / (t3 - t_upper);
            thetadot = theta_upper * degtorad / (t3 - t_upper);

            double tst = t4 - t3;
            new_dw->put(delt_vartype(tst), d_sharedState->get_delt_label(),
                        getLevel(patches));
          }
        }
      }

      // Construct Qdot:
      const double costhetat = cos(thetat);
      const double sinthetat = sin(thetat);
      Matrix3 aa(rot_axis_upper, rot_axis_upper);
      Matrix3 AA(rot_axis_upper);
      auto Qdot =
        (Identity - aa) * (-sinthetat * thetadot) + AA * costhetat * thetadot;

      // Now we need to compute the total previous rotation:
      Matrix3 R_previous;
      R_previous.Identity();
      for (auto ii = 0; ii < t_upper_index; ii++) {
        auto thetai = d_prescribedAngle[ii] * degtorad;
        auto ai = d_prescribedRotationAxis[ii];
        Matrix3 Qi(thetai, ai);
        R_previous = Qi * R_previous;
      }

      // Fstar is the deformation gradient with the superimposed rotations
      // included
      // Fdotstar is the rate of the deformation gradient with superimposed
      // rotations included
      auto Fstar = Qt * R_previous * Ft;
      auto Fdotstar = Qdot * R_previous * Ft + Qt * R_previous * Fdot;
      auto R_previous_inv = R_previous.Inverse();

      // Update grid velocities
      for (auto iter = patch->getExtraNodeIterator(); !iter.done(); iter++) {
        IntVector n = *iter;
        Vector position = patch->getNodePosition(n).asVector();

        // Exact Deformation Update
        if (d_flags->d_exactDeformation) {
          gVelocity[n] = (F_upper * F_lower.Inverse() - Identity) *
                              R_previous_inv * QQ * position / delT;
        } else {
          gVelocity[n] =
            Fdotstar * Ft_inv * R_previous_inv * QQ * position;
        }

      } // Node Iterator

    } // matl loop
  }   // patch loop
}

/*!----------------------------------------------------------------------
 * scheduleCheckGridVelocity
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleCheckGridVelocity(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, 
                "MPM::scheduleCheckGridVelocity");

  Task* t = scinew Task("UofU_MPM::checkGridVelocity", this,
                        &UofU_MPM::checkGridVelocity);

  t->requires(Task::OldDW, d_labels->pParticleIDLabel,   Ghost::None);
  t->requires(Task::OldDW, d_labels->pXLabel,            Ghost::None);
  t->requires(Task::OldDW, d_labels->pSizeLabel,         Ghost::None);
  t->requires(Task::OldDW, d_labels->pDefGradLabel,      Ghost::None);
  t->requires(Task::NewDW, d_labels->gVelocityLabel, Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::NewDW, d_labels->gAccelerationLabel, Ghost::AroundCells, d_numGhostNodes);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * checkGridVelocity
 *-----------------------------------------------------------------------*/
void
UofU_MPM::checkGridVelocity(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing,
            "Doing checkGridVelocity");

  for (auto patch : *patches) {

    //Vector dx = patch->dCell();
    //Vector oodx = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    auto interpolator = d_flags->d_interpolator->clone(patch);
    auto num_influence_nodes = interpolator->size();
    std::vector<IntVector> influenceNodes(num_influence_nodes);
    std::vector<double> S_ip_av(num_influence_nodes);
    std::vector<Vector> dS_ip_av(num_influence_nodes);

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){

      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );

      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<long64>  pParticleID;
      constParticleVariable<Point>   pX;
      constParticleVariable<Matrix3> pDefGrad_old, pSize;
      constNCVariable<Vector>        gVelocity, gAcceleration;

      old_dw->get(pParticleID,   d_labels->pParticleIDLabel,   pset);
      old_dw->get(pX,            d_labels->pXLabel,            pset);
      old_dw->get(pSize,         d_labels->pSizeLabel,         pset);
      old_dw->get(pDefGrad_old,  d_labels->pDefGradLabel,      pset);
      new_dw->get(gAcceleration, d_labels->gAccelerationLabel, matID, patch, 
                  Ghost::AroundCells, d_numGhostNodes);
      new_dw->get(gVelocity, d_labels->gVelocityLabel, matID, patch, 
                  Ghost::AroundCells, d_numGhostNodes);

      for (auto particle : *pset) {

        interpolator->findCellAndWeightsAndShapeDerivatives(pX[particle], 
                                                            influenceNodes,
                                                            S_ip_av, dS_ip_av,
                                                            pSize[particle],
                                                            pDefGrad_old[particle]);
        for (int k = 0; k < num_influence_nodes; k++) {
          IntVector node = influenceNodes[k];
          auto v_i_old = gVelocity[node];
          auto a_i_old = gAcceleration[node];
          std::cout << " patch = " << patch << " node = " << node 
                    << " gAcc = " << a_i_old 
                    << " gVel = " << v_i_old << "\n";
        }

      } // End of loop over particles

    } // end of materials loop
  } // end of patch loop
}

/*!----------------------------------------------------------------------
 * scheduleComputeVelocityAndDeformationGradient
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeVelocityAndDeformationGradient(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  /* Create a task for computing the deformation gradient */
  printSchedule(patches, cout_doing, 
                "MPM::scheduleComputeVelocityAndDeformationGradient");

  Task* t = scinew Task("UofU_MPM::computeVelocityAndDeformationGradient", this,
                        &UofU_MPM::computeVelocityAndDeformationGradient);

  t->requires(Task::OldDW, d_labels->delTLabel);
  t->requires(Task::NewDW, d_labels->gVelocityLabel, Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::NewDW, d_labels->gAccelerationLabel, Ghost::AroundCells, d_numGhostNodes);

  t->requires(Task::OldDW, d_labels->pParticleIDLabel,   Ghost::None);
  t->requires(Task::OldDW, d_labels->pXLabel,            Ghost::None);
  t->requires(Task::OldDW, d_labels->pMassLabel,         Ghost::None);
  t->requires(Task::OldDW, d_labels->pVolumeLabel,       Ghost::None);
  t->requires(Task::OldDW, d_labels->pSizeLabel,         Ghost::None);
  t->requires(Task::OldDW, d_labels->pDefGradLabel,      Ghost::None);

  t->computes(d_labels->pVelGradLabel_preReloc);
  t->computes(d_labels->pDefGradMidLabel);
  t->computes(d_labels->pDefGradLabel_preReloc);
  t->computes(d_labels->pVolumeMidLabel);
  t->computes(d_labels->pVolumeLabel_preReloc);

  int numMPMMatls = d_sharedState->getNumMPMMatls();
  for (int m = 0; m < numMPMMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    const MaterialSubset* matlset = mpm_matl->thisMaterial();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
      t->requires(Task::OldDW, d_labels->pRemoveLabel, matlset, Ghost::None);
      t->computes(d_labels->pRemoveLabel_preReloc, matlset);
      t->computes(d_labels->pPolarDecompRMidLabel, matlset);
      t->computes(d_labels->pPolarDecompRLabel_preReloc, matlset);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeVelocityAndDeformationGradient
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeVelocityAndDeformationGradient(const ProcessorGroup*,
                                                const PatchSubset* patches,
                                                const MaterialSubset*,
                                                DataWarehouse* old_dw,
                                                DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing,
            "Doing computeVelocityAndDeformationGradient");

  delt_vartype delT;
  old_dw->get(delT, d_labels->delTLabel, getLevel(patches));

  for (int pp = 0; pp < patches->size(); pp++) {
    const Patch* patch = patches->get(pp);
    Vector dx = patch->dCell();
    Vector oodx = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    auto interpolator = d_flags->d_interpolator->clone(patch);
    auto num_influence_nodes = interpolator->size();

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){

      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<double>  pVolume_old;
      constParticleVariable<long64>  pParticleID;
      constParticleVariable<Point>   pX;
      constParticleVariable<Matrix3> pDefGrad_old, pSize;
      constNCVariable<Vector>        gVelocity, gAcceleration;

      old_dw->get(pParticleID,   d_labels->pParticleIDLabel,   pset);
      old_dw->get(pX,            d_labels->pXLabel,            pset);
      old_dw->get(pVolume_old,   d_labels->pVolumeLabel,       pset);
      old_dw->get(pSize,         d_labels->pSizeLabel,         pset);
      old_dw->get(pDefGrad_old,  d_labels->pDefGradLabel,      pset);
      new_dw->get(gAcceleration, d_labels->gAccelerationLabel, matID, patch, 
                  Ghost::AroundCells, d_numGhostNodes);
      new_dw->get(gVelocity, d_labels->gVelocityLabel, matID, patch, 
                  Ghost::AroundCells, d_numGhostNodes);

      ParticleVariable<int>        pRemove_new;
      ParticleVariable<double>     pVolume_mid, pVolume_new;
      ParticleVariable<Matrix3>    pVelGrad_mid, pDefGrad_mid, pPolarDecompR_mid;
      ParticleVariable<Matrix3>    pDefGrad_new, pPolarDecompR_new;

      new_dw->allocateAndPut(pVolume_mid,       d_labels->pVolumeMidLabel,             pset);
      new_dw->allocateAndPut(pVolume_new,       d_labels->pVolumeLabel_preReloc,       pset);
      new_dw->allocateAndPut(pVelGrad_mid,      d_labels->pVelGradLabel_preReloc,      pset);
      new_dw->allocateAndPut(pDefGrad_mid,      d_labels->pDefGradMidLabel,            pset);
      new_dw->allocateAndPut(pDefGrad_new,      d_labels->pDefGradLabel_preReloc,      pset);
      if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
        new_dw->allocateAndPut(pRemove_new,       d_labels->pRemoveLabel_preReloc,       pset);
        new_dw->allocateAndPut(pPolarDecompR_new, d_labels->pPolarDecompRLabel_preReloc, pset);
        new_dw->allocateAndPut(pPolarDecompR_mid, d_labels->pPolarDecompRMidLabel,       pset);
      }
      
      for (auto particle : *pset) {

        std::vector<IntVector> influenceNodes(num_influence_nodes);
        std::vector<double> S_ip_av(num_influence_nodes);
        std::vector<Vector> dS_ip_av(num_influence_nodes);
        interpolator->findCellAndWeightsAndShapeDerivatives(pX[particle], 
                                                            influenceNodes,
                                                            S_ip_av, dS_ip_av,
                                                            pSize[particle],
                                                            pDefGrad_old[particle]);
        Matrix3 f_p_old(0.0);
        Matrix3 f_p_mid(0.0), f_p_new(0.0);
        Matrix3 fdot_p_mid(0.0), fdot_p_new(0.0);
        for (int k = 0; k < num_influence_nodes; k++) {
          IntVector node = influenceNodes[k];
          auto x_i_old = patch->getNodePosition(node);
          auto v_i_old = gVelocity[node];
          auto a_i_old = gAcceleration[node];
          //auto x_i_mid = x_i_old + v_i_old * (0.5 * delT) + a_i_old * (0.125 * delT * delT);
          auto v_i_quart = v_i_old + a_i_old * (0.25 * delT);
          auto x_i_mid = x_i_old + v_i_quart * (0.25 * delT);
          //auto x_i_new = x_i_old + v_i_old * delT + a_i_old * (0.5 * delT * delT);
          auto v_i_mid = v_i_old + a_i_old * (0.5 * delT);
          auto x_i_new = x_i_old + v_i_mid * delT;

          Vector G_ip_av = dS_ip_av[k] * oodx;
          Matrix3 xG_old(x_i_old.asVector(), G_ip_av);
          Matrix3 xG_mid(x_i_mid.asVector(), G_ip_av);
          Matrix3 xG_new(x_i_new.asVector(), G_ip_av);
          Matrix3 vG_mid(v_i_mid, G_ip_av);

          /*
          if (!std::isfinite(gAcceleration[node].length())) {
          }
          */
          /*
          if (f_p_new.Determinant() < 0.0) {
            std::cout << "k = " << k << " node = " << node 
                      << " v_i " << v_i_old
                      << " a_i " << a_i_old << "\n"
                      << " x_i = " << x_i_old << " G = " << G_ip_av << "\n"
                      << " xG_n " << xG_old << "\n"
                      << " xG_mid " << xG_mid << "\n"
                      << " xG_new " << xG_new << "\n";
          }
          */

          if (d_flags->d_axisymmetric) {
            // x -> r, y -> z, z -> theta
            xG_mid(2,2) = x_i_mid.x() * G_ip_av.z();
            xG_new(2,2) = x_i_new.x() * G_ip_av.z();
            vG_mid(2,2) = v_i_mid.x() * G_ip_av.z();
          }

          f_p_old += xG_old;
          f_p_mid += xG_mid;
          f_p_new += xG_new;
          fdot_p_mid += vG_mid;
        }

        //if (f_p_new.Determinant() < 0.0) {
        //}
        
        /*
        if (f_p_old.Determinant() != 1.0) {
          std::cout <<  "xG_n = " << f_p_old << "\n"
                    <<  "xG_mid = " << f_p_mid << "\n"
                    <<  "xG_new = " << f_p_new << "\n";
        }
        */

        pVelGrad_mid[particle] = fdot_p_mid * f_p_mid.Inverse();
        pDefGrad_mid[particle] = f_p_mid * pDefGrad_old[particle];
        pDefGrad_new[particle] = f_p_new * pDefGrad_old[particle];
        pVolume_mid[particle] = pVolume_old[particle] * f_p_mid.Determinant();
        pVolume_new[particle] = pVolume_old[particle] * f_p_new.Determinant();

        // Check 1: Look at Jacobian
        double J = pDefGrad_new[particle].Determinant();
        if (!(J > 0.0)) {
          std::ostringstream msg;
          msg << "**ERROR** Negative Jacobian of deformation gradient"
              << " J = " << J << " in material # ="
              << mpm_matl << " and particle " << pParticleID[particle]  << " which has volume "
              << pVolume_new[particle] << "\n";
          msg << "matl = "  << mpm_matl << " matID = " << matID << " particle = " << particle
              << " particleID = " << pParticleID[particle] << "\n";
          msg << "velGrad = " << pVelGrad_mid[particle] << "\n";
          msg << "F_old = " << pDefGrad_old[particle]     << "\n";
          msg << "F_inc = " << f_p_new       << "\n";
          msg << "F_new = " << pDefGrad_new[particle] << "\n";

          #ifdef IGNORE_NEGATIVE_JACOBIANS
            msg << "\t Ignoring new deformation gradient\n";
            J = pDefGrad_old[particle].Determinant();
            pDefGrad_new[particle] = pDefGrad_old[particle];
            std::cerr << msg.str();
          #else
            throw InvalidValue(msg.str(), __FILE__, __LINE__);
          #endif
        }

        if (!d_flags->d_doPressureStabilization) {
          if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
            Matrix3 RR, UU;
            Matrix3 FF_new = pDefGrad_new[particle];
            double Fmax_new = FF_new.MaxAbsElem();
            double JJ_new = FF_new.Determinant();
            // These checks prevent failure of the polar decomposition algorithm 
            // if [F_new] has some extreme values.
            if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
              pRemove_new[particle] = -999;
              proc0cout << "Deformation gradient component unphysical: [F] = " << FF_new
                        << std::endl;
              proc0cout << "Resetting [F]=[I] for this step and deleting particle"
                        << " idx = " << particle << " particleID = " << pParticleID[particle]
                        << std::endl;
              Identity.polarDecompositionRMB(UU, RR);
            } else {
              pRemove_new[particle] = 0;
              FF_new.polarDecompositionRMB(UU, RR);
            }
            pPolarDecompR_new[particle] = RR;

            pDefGrad_mid[particle].polarDecompositionRMB(UU, RR);
            pPolarDecompR_mid[particle] = RR;
          }
        }

        //std::cout << "After compute vel/def gradients: material = " << m
        //          << " particle = " << particle << "\n"
        //          << " L_mid = " << pVelGrad_mid[particle] << "\n"
        //          << " F_old = " << pDefGrad_old[particle] << "\n"
        //          << " F_mid = " << pDefGrad_mid[particle] << "\n"
        //          << " F_new = " << pDefGrad_new[particle] << "\n"
        //          << " R_mid = " << pPolarDecompR_mid[particle] << "\n"
        //          << " R_new = " << pPolarDecompR_new[particle] << "\n";

      } // End of loop over particles

      // The following is used only for pressure stabilization
      if (d_flags->d_doPressureStabilization) {

        double rho_orig = mpm_matl->getInitialDensity();

        constParticleVariable<double>  pMass;
        old_dw->get(pMass, d_labels->pMassLabel, pset);
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);

        CCVariable<double> vol_0_CC, vol_CC_mid, vol_CC_new;
        new_dw->allocateTemporary(vol_0_CC,   patch);
        new_dw->allocateTemporary(vol_CC_mid, patch);
        new_dw->allocateTemporary(vol_CC_new, patch);
        vol_0_CC.initialize(0.);
        vol_CC_mid.initialize(0.);
        vol_CC_new.initialize(0.);
      
        // Step 1: loop thru particles and cells to compute cell centered J
        for (auto particle : *pset) {

          IntVector cell_index;
          patch->findCell(pX[particle],cell_index);

          vol_CC_mid[cell_index] += pVolume_mid[particle];
          vol_CC_new[cell_index] += pVolume_new[particle];
          vol_0_CC[cell_index] += pMass[particle]/rho_orig;
        }

        // Compute cell centered J
        CCVariable<double> J_CC_mid, J_CC_new;
        new_dw->allocateTemporary(J_CC_mid, patch);
        new_dw->allocateTemporary(J_CC_new, patch);
        J_CC_mid.initialize(0.);
        J_CC_new.initialize(0.);
        for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
          IntVector c = *iter;
          J_CC_mid[c] = vol_CC_mid[c] / vol_0_CC[c];
          J_CC_new[c] = vol_CC_new[c] / vol_0_CC[c];
        }

        // Step 2: loop thru particles again to compute corrected def grad
        for (auto particle : *pset) {

          IntVector cell_index;
          patch->findCell(pX[particle], cell_index);

          // get the original volumetric part of the deformation
          double J_mid = pDefGrad_mid[particle].Determinant();
          double J_new = pDefGrad_new[particle].Determinant();

          // Change F such that the determinant is equal to the average for
          // the cell
          pDefGrad_mid[particle] *= cbrt(J_CC_mid[cell_index]/J_mid);
          pDefGrad_new[particle] *= cbrt(J_CC_new[cell_index]/J_new);

          // Update the deformed volume
          J_mid = pDefGrad_mid[particle].Determinant();
          J_new = pDefGrad_new[particle].Determinant();
          pVolume_mid[particle] = (pMass[particle]/rho_orig)*J_mid;
          pVolume_new[particle] = (pMass[particle]/rho_orig)*J_new;

          // Check 1: Look at Jacobian
          if (!(J_new > 0.0)) {

            std::ostringstream msg;
            msg << "**ERROR** Negative Jacobian of deformation gradient (J) = " << J_new 
                << " after pressure stablization in particle ID "
                << pParticleID[particle] << " matl = "  << mpm_matl << "\n";
            msg << "F_old = " << pDefGrad_old[particle] << "\n";
            msg << "F_new = " << pDefGrad_new[particle] << "\n";
            throw InvalidValue(msg.str(), __FILE__, __LINE__);
          }

          if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {
            Matrix3 RR, UU;
            Matrix3 FF_new = pDefGrad_new[particle];
            double Fmax_new = FF_new.MaxAbsElem();
            double JJ_new = FF_new.Determinant();

            // These checks prevent failure of the polar decomposition algorithm 
            // if [F_new] has some extreme values.
            if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
              pRemove_new[particle] = -999;
              proc0cout << "Deformation gradient component unphysical: [F] = " << FF_new
                        << std::endl;
              proc0cout << "Resetting [F]=[I] for this step and deleting particle"
                        << " idx = " << particle << " particleID = " << pParticleID[particle]
                        << std::endl;
              Identity.polarDecompositionRMB(UU, RR);
            } else {
              pRemove_new[particle] = 0;
              FF_new.polarDecompositionRMB(UU, RR);
            }
            pPolarDecompR_new[particle] = RR;

            pDefGrad_mid[particle].polarDecompositionRMB(UU, RR);
            pPolarDecompR_mid[particle] = RR;
          }

        }

      } //end of pressureStabilization loop
    } // end of materials loop
  } // end of patch loop
}

/*!----------------------------------------------------------------------
 * scheduleUnrotateStressAndDeformationRate
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleUnrotateStressAndDeformationRate(SchedulerP& sched,
                                                   const PatchSet* patches,
                                                   const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, 
                "MPM::scheduleUnrotateStressAndDeformationRate");

  int numMatls = d_sharedState->getNumMPMMatls();
  Task* t = scinew Task("UofU_MPM::computeUnrotatedStressAndDeformationRate", this,
                        &UofU_MPM::computeUnrotatedStressAndDeformationRate);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->requires(Task::OldDW, d_labels->pParticleIDLabel,      matlset, Ghost::None);
      t->requires(Task::OldDW, d_labels->pPolarDecompRLabel,    matlset, Ghost::None);
      t->requires(Task::OldDW, d_labels->pStressLabel,          matlset, Ghost::None);
      t->requires(Task::NewDW, d_labels->pVelGradLabel_preReloc, matlset, Ghost::None);
      t->requires(Task::NewDW, d_labels->pPolarDecompRMidLabel, matlset, Ghost::None);

      t->computes(d_labels->pDeformRateMidLabel,   matlset);
      t->computes(d_labels->pStressUnrotatedLabel, matlset);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeUnrotatedStressAndDeformationRate
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeUnrotatedStressAndDeformationRate(const ProcessorGroup*, 
                                                   const PatchSubset* patches,
                                                   const MaterialSubset*, 
                                                   DataWarehouse* old_dw,
                                                   DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing, "Doing computeUnrotated");

  int numMatls = d_sharedState->getNumMPMMatls();
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pStress_old, pR_old, pR_mid, pVelGrad_mid;
        old_dw->get(pParticleID,  d_labels->pParticleIDLabel,      pset);
        old_dw->get(pStress_old,  d_labels->pStressLabel,          pset);
        old_dw->get(pR_old,       d_labels->pPolarDecompRLabel,    pset);
        new_dw->get(pVelGrad_mid, d_labels->pVelGradLabel_preReloc, pset);
        new_dw->get(pR_mid,       d_labels->pPolarDecompRMidLabel, pset);

        ParticleVariable<Matrix3> pDeformRate_mid, pStress_old_unrotated;
        new_dw->allocateAndPut(pDeformRate_mid,       d_labels->pDeformRateMidLabel,   pset);
        new_dw->allocateAndPut(pStress_old_unrotated, d_labels->pStressUnrotatedLabel, pset);

        for (auto particle : *pset) {
          pStress_old_unrotated[particle] = (pR_old[particle].Transpose()) * (pStress_old[particle] *
                                             pR_old[particle]);
          Matrix3 DD = (pVelGrad_mid[particle] + pVelGrad_mid[particle].Transpose()) * 0.5;
          pDeformRate_mid[particle] = (pR_mid[particle].Transpose()) * (DD * pR_mid[particle]);
          //std::cout << "After unrotate: material = " << m
          //        << " particle = " << particle << "\n"
          //        << " sig_old = " << pStress_old_unrotated[particle] << "\n"
          //        << " D_mid = " << pDeformRate_mid[particle] << "\n";
        }
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * computeStressTensor
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeStressTensor(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  /* Create a task for computing the stress tensor */
  printSchedule(patches, cout_doing, "MPM::scheduleComputeStressTensor");

  int numMatls = d_sharedState->getNumMPMMatls();
  Task* t = scinew Task("UofU_MPM::computeStressTensor", this,
                        &UofU_MPM::computeStressTensor);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    const MaterialSubset* matlset = mpm_matl->thisMaterial();

    // Add requires and computes for constitutive model
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
    t->computes(d_labels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_sharedState->get_delt_label(), getLevel(patches));

  if (d_flags->d_reductionVars->accStrainEnergy || d_flags->d_reductionVars->strainEnergy) {
    t->computes(d_labels->StrainEnergyLabel);
  }

  sched->addTask(t, patches, matls);
}

void
UofU_MPM::computeStressTensor(const ProcessorGroup*, const PatchSubset* patches,
                              const MaterialSubset*, DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing, "Doing computeStressTensor");

  for (int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

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
  }
}

/*!----------------------------------------------------------------------
 * scheduleRotateStress
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleRotateStress(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, 
                "MPM::scheduleRotateStress");

  int numMatls = d_sharedState->getNumMPMMatls();
  Task* t = scinew Task("UofU_MPM::computeRotatedStress", this,
                        &UofU_MPM::computeRotatedStress);
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      const MaterialSubset* matlset = mpm_matl->thisMaterial();

      t->requires(Task::OldDW, d_labels->pParticleIDLabel,      matlset, Ghost::None);
      t->requires(Task::NewDW, d_labels->pPolarDecompRLabel_preReloc, matlset, Ghost::None);

      t->modifies(d_labels->pStressLabel_preReloc, matlset);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeRotatedStress
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeRotatedStress(const ProcessorGroup*, 
                               const PatchSubset* patches,
                               const MaterialSubset*, 
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing, "Doing computeRotatedStress");

  int numMatls = d_sharedState->getNumMPMMatls();
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    if (cm->modelType() == ConstitutiveModel::ModelType::INCREMENTAL) {

      int matID = mpm_matl->getDWIndex();
      for (int p = 0; p < patches->size(); p++) {
        const Patch* patch = patches->get(p);
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

        constParticleVariable<long64> pParticleID;
        constParticleVariable<Matrix3> pR_new;
        ParticleVariable<Matrix3> pStress_new;
        old_dw->get(pParticleID,           d_labels->pParticleIDLabel,            pset);
        new_dw->get(pR_new,                d_labels->pPolarDecompRLabel_preReloc, pset);
        new_dw->getModifiable(pStress_new, d_labels->pStressLabel_preReloc,       pset);

        for (auto particle : *pset) {
          pStress_new[particle] = (pR_new[particle] * pStress_new[particle]) *
                                             (pR_new[particle].Transpose());
          //std::cout << "After rotate stress: material = " << m
          //        << " particle = " << particle << "\n"
          //        << " pStress_new = " << pStress_new[particle] << "\n";
        }
      }
    }
  }
}

/*!----------------------------------------------------------------------
 * scheduleComputeBasicDamage
 *-----------------------------------------------------------------------*/
void 
UofU_MPM::scheduleComputeBasicDamage(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(), 
                           getLevel(patches)->getGrid()->numLevels()))
    return;
  
  /* Create a task for computing the damage variables */
  printSchedule(patches,cout_doing,"MPM::scheduleComputeBasicDamage");
  
  int numMatls = d_sharedState->getNumMPMMatls();
  Task* t = scinew Task("MPM::computeBasicDamage",
                        this, &UofU_MPM::computeBasicDamage);
  for(int m = 0; m < numMatls; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

    // Add requires and computes for vel grad/def grad
    if (mpm_matl->d_doBasicDamage) {    
      Vaango::BasicDamageModel* d_basicDamageModel = mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addComputesAndRequires(t, mpm_matl, patches, d_labels);
    }
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeBasicDamage
 *-----------------------------------------------------------------------*/
void 
UofU_MPM::computeBasicDamage(const ProcessorGroup*,
                              const PatchSubset* patches,
                              const MaterialSubset* ,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{

  printTask(patches, patches->get(0), cout_doing, "Doing computeBasicDamage");

  for(int m = 0; m < d_sharedState->getNumMPMMatls(); m++){

    if (cout_dbg.active()) cout_dbg << " Patch = " << (patches->get(0))->getID() << " Mat = " << m;

    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    if (cout_dbg.active()) cout_dbg << " MPM_Mat = " << mpm_matl;

    // Compute basic damage
    if (mpm_matl->d_doBasicDamage) { 
      Vaango::BasicDamageModel* basicDamageModel = mpm_matl->getBasicDamageModel();
      basicDamageModel->computeBasicDamage(patches, mpm_matl, old_dw, new_dw, d_labels);
      if (cout_dbg.active()) cout_dbg << " Damage model = " << basicDamageModel;
    }

    if (cout_dbg.active()) cout_dbg << " Exit\n" ;
  }
}

/*!----------------------------------------------------------------------
 * scheduleUpdateErosionParameter
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleUpdateErosionParameter(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM::scheduleUpdateErosionParameter");

  Task* t = scinew Task("UofU_MPM::updateErosionParameter", this,
                        &UofU_MPM::updateErosionParameter);
  int numMatls = d_sharedState->getNumMPMMatls();
  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

    if (mpm_matl->d_doBasicDamage) {
      Vaango::BasicDamageModel* d_basicDamageModel =
        mpm_matl->getBasicDamageModel();
      d_basicDamageModel->addRequiresLocalizationParameter(t, mpm_matl,
                                                           patches);
    } else {
      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
      cm->addRequiresDamageParameter(t, mpm_matl, patches);
    }
  }
  t->computes(d_labels->pLocalizedMPMLabel);

  if (d_flags->d_deleteRogueParticles) {
    t->requires(Task::OldDW, d_labels->pXLabel, Ghost::None);
    t->computes(d_labels->numLocInCellLabel);
    t->computes(d_labels->numInCellLabel);
  }

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * updateErosionParameter
 *-----------------------------------------------------------------------*/
void
UofU_MPM::updateErosionParameter(const ProcessorGroup*,
                                 const PatchSubset* patches,
                                 const MaterialSubset*, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing updateErosionParameter");

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      if (cout_dbg.active())
        cout_dbg << "updateErosionParameter:: material # = " << m << endl;

      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      if (cout_dbg.active()) {
        cout_dbg << "updateErosionParameter:: mpm_matl* = " << mpm_matl
                 << " matID = " << matID << " pset* = " << pset << endl;
      }

      // Get the localization info
      ParticleVariable<int> isLocalized;
      new_dw->allocateAndPut(isLocalized, d_labels->pLocalizedMPMLabel, pset);
      for (auto particle : *pset) {
        isLocalized[particle] = 0;
      }
      if (mpm_matl->d_doBasicDamage) {
        Vaango::BasicDamageModel* basicDamageModel =
          mpm_matl->getBasicDamageModel();
        basicDamageModel->getLocalizationParameter(patch, isLocalized, matID,
                                                   old_dw, new_dw);
      } else {
        mpm_matl->getConstitutiveModel()->getDamageParameter(
          patch, isLocalized, matID, old_dw, new_dw);
      }
      //for (auto particle : *pset) {
      //  std::cout << "After isLocalized: material = " << m << " particle = " << particle
      //            << " localized = " << isLocalized[particle] << "\n";
      //}

      if (cout_dbg.active())
        cout_dbg << "updateErosionParameter:: Got Damage Parameter" << endl;

      if (d_flags->d_deleteRogueParticles) {
        // The following looks for localized particles that are isolated
        // either individually or in small groups
        // Ghost::GhostType  Ghost::AroundCells = Ghost::AroundCells;
        CCVariable<int> numLocInCell, numInCell;
        new_dw->allocateAndPut(numLocInCell, d_labels->numLocInCellLabel, matID,
                               patch);
        new_dw->allocateAndPut(numInCell, d_labels->numInCellLabel, matID, patch);
        numLocInCell.initialize(0);
        numInCell.initialize(0);

        constParticleVariable<Point> pX;
        old_dw->get(pX, d_labels->pXLabel, pset);

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

      if (cout_dbg.active())
        cout_dbg << "updateErosionParameter:: Updated Erosion " << endl;
    }

    if (cout_dbg.active())
      cout_dbg << "Done updateErosionParamter on patch " << patch->getID()
               << "\t MPM" << endl;
  }
}

/*!----------------------------------------------------------------------
 * scheduleFindRogueParticles
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleFindRogueParticles(SchedulerP& sched, const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  if (d_flags->d_deleteRogueParticles) {
    printSchedule(patches, cout_doing, "MPM::scheduleFindRogueParticles");

    Task* t = scinew Task("UofU_MPM::findRogueParticles", this,
                          &UofU_MPM::findRogueParticles);
    t->requires(Task::NewDW, d_labels->numLocInCellLabel, Ghost::AroundCells, 1);
    t->requires(Task::NewDW, d_labels->numInCellLabel, Ghost::AroundCells, 1);
    t->requires(Task::OldDW, d_labels->pXLabel, Ghost::None);
    t->modifies(d_labels->pLocalizedMPMLabel);

    sched->addTask(t, patches, matls);
  }
}

/*!----------------------------------------------------------------------
 * findRogueParticles
 *-----------------------------------------------------------------------*/
void
UofU_MPM::findRogueParticles(const ProcessorGroup*, const PatchSubset* patches,
                             const MaterialSubset*, DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing findRogueParticles");

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // The following looks for localized particles that are isolated
      // either individually or in small groups
      constCCVariable<int> numLocInCell, numInCell;
      constParticleVariable<Point> pX;

      ParticleVariable<int> isLocalized;

      new_dw->get(numLocInCell, d_labels->numLocInCellLabel, matID, patch, Ghost::AroundCells,
                  1);
      new_dw->get(numInCell, d_labels->numInCellLabel, matID, patch, Ghost::AroundCells, 1);
      old_dw->get(pX, d_labels->pXLabel, pset);
      new_dw->getModifiable(isLocalized, d_labels->pLocalizedMPMLabel, pset);

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
            isLocalized[particle] = -999;
          }
        } // if localized
        //std::cout << "After findRogue: material = " << m << " particle = " << particle
        //          << " localized = " << isLocalized[particle] << "\n";
      }   // particles
    }     // matls
  }       // patches
}

/*!----------------------------------------------------------------------
 * scheduleComputeAccStrainEnergy
 *   Compute the accumulated strain energy
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeAccStrainEnergy(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;
  printSchedule(patches, cout_doing, "MPM::scheduleComputeAccStrainEnergy");

  Task* t = scinew Task("UofU_MPM::computeAccStrainEnergy", this,
                        &UofU_MPM::computeAccStrainEnergy);
  t->requires(Task::OldDW, d_labels->AccStrainEnergyLabel);
  t->requires(Task::NewDW, d_labels->StrainEnergyLabel);
  t->computes(d_labels->AccStrainEnergyLabel);
  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeAccStrainEnergy
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeAccStrainEnergy(const ProcessorGroup*, const PatchSubset*,
                                 const MaterialSubset*, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  // Get the totalStrainEnergy from the old datawarehouse
  max_vartype accStrainEnergy;
  old_dw->get(accStrainEnergy, d_labels->AccStrainEnergyLabel);

  // Get the incremental strain energy from the new datawarehouse
  sum_vartype incStrainEnergy;
  new_dw->get(incStrainEnergy, d_labels->StrainEnergyLabel);

  // Add the two a put into new dw
  double totalStrainEnergy = (double)accStrainEnergy + (double)incStrainEnergy;
  new_dw->put(max_vartype(totalStrainEnergy), d_labels->AccStrainEnergyLabel);
}

/*!----------------------------------------------------------------------
 * printParticleLabels
 *-----------------------------------------------------------------------*/
void
UofU_MPM::printParticleLabels(vector<const VarLabel*> labels, DataWarehouse* dw,
                              int matID, const Patch* patch)
{
  for (const auto& label : labels) {
    if (dw->exists(label, matID, patch))
      std::cout << label->getName() << " exists\n";
    else
      std::cout << label->getName() << " does NOT exist\n";
  }
}

/*!----------------------------------------------------------------------
 * scheduleInterpolateToParticlesAndUpdate
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* matls)

{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing,
                "MPM::scheduleInterpolateToParticlesAndUpdate");

  Task* t = scinew Task("UofU_MPM::interpolateToParticlesAndUpdate", this,
                        &UofU_MPM::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_sharedState->get_delt_label());

  t->requires(Task::NewDW, d_labels->gAccelerationLabel, Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::NewDW, d_labels->gVelocityLabel, Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::NewDW, d_labels->frictionalWorkLabel, Ghost::AroundCells, d_numGhostNodes);
  t->requires(Task::OldDW, d_labels->pXLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pMassLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pTemperatureLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pTempPreviousLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pParticleIDLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pVelocityLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pDispLabel, Ghost::None);
  t->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->pLocalizedMPMLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->pDefGradLabel_preReloc, Ghost::None);
  t->modifies(d_labels->pVolumeLabel_preReloc);

  if (d_flags->d_useLoadCurves) {
    t->requires(Task::OldDW, d_labels->pLoadCurveIDLabel, Ghost::None);
  }

  t->computes(d_labels->pXLabel_preReloc);
  t->computes(d_labels->pMassLabel_preReloc);
  t->computes(d_labels->pTemperatureLabel_preReloc);
  t->computes(d_labels->pTempPreviousLabel_preReloc);
  t->computes(d_labels->pDispLabel_preReloc);
  t->computes(d_labels->pVelocityLabel_preReloc);
  t->computes(d_labels->pParticleIDLabel_preReloc);
  t->computes(d_labels->pSizeLabel_preReloc);
  t->computes(d_labels->pXXLabel);

  //__________________________________
  //  reduction variables
  if (d_flags->d_reductionVars->momentum) {
    t->computes(d_labels->TotalMomentumLabel);
  }
  if (d_flags->d_reductionVars->KE) {
    t->computes(d_labels->KineticEnergyLabel);
  }
  if (d_flags->d_reductionVars->centerOfMass) {
    t->computes(d_labels->CenterOfMassPositionLabel);
  }
  if (d_flags->d_reductionVars->mass) {
    t->computes(d_labels->TotalMassLabel);
  }
  if (d_flags->d_reductionVars->volDeformed) {
    t->computes(d_labels->TotalVolumeDeformedLabel);
  }

  // debugging scalar
  if (d_flags->d_withColor) {
    t->requires(Task::OldDW, d_labels->pColorLabel, Ghost::None);
    t->computes(d_labels->pColorLabel_preReloc);
  }

  MaterialSubset* z_matl = scinew MaterialSubset();
  z_matl->add(0);
  z_matl->addReference();
  t->requires(Task::OldDW, d_labels->NC_CCweightLabel, z_matl, Ghost::None);
  t->computes(d_labels->NC_CCweightLabel, z_matl);

  sched->addTask(t, patches, matls);

  // The task will have a reference to z_matl
  if (z_matl->removeReference())
    delete z_matl; // shouln't happen, but...
}

/*!----------------------------------------------------------------------
 * interpolateToParticlesAndUpdate
 *-----------------------------------------------------------------------*/
void
UofU_MPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                          const PatchSubset* patches,
                                          const MaterialSubset*,
                                          DataWarehouse* old_dw,
                                          DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing,
              "Doing interpolateToParticlesAndUpdate");

    auto interpolator = d_flags->d_interpolator->clone(patch);
    auto num_influence_nodes = interpolator->size();
    vector<IntVector> influenceNodes(num_influence_nodes);
    vector<double> S_ip_av(num_influence_nodes);

    double totalmass = 0;
    double partvoldef = 0.;
    Vector centerOfMass(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke = 0;
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    delt_vartype delT;
    old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches));

    // Copy NC_CCweight (only material 0)
    constNCVariable<double> NC_CCweight;
    NCVariable<double> NC_CCweight_new;
    old_dw->get(NC_CCweight, d_labels->NC_CCweightLabel, 0, patch, Ghost::None, 0);
    new_dw->allocateAndPut(NC_CCweight_new, d_labels->NC_CCweightLabel, 0,
                           patch);
    NC_CCweight_new.copyData(NC_CCweight);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Copy data that doesn't change
      constParticleVariable<long64> pParticleID;
      ParticleVariable<long64> pParticleID_new;
      old_dw->get(pParticleID, d_labels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pParticleID_new,
                             d_labels->pParticleIDLabel_preReloc, pset);
      pParticleID_new.copyData(pParticleID);

      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pSize_new;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(pSize_new, d_labels->pSizeLabel_preReloc, pset);
      pSize_new.copyData(pSize);

      // Get particle variables
      constParticleVariable<int> pLocalized_new;
      new_dw->get(pLocalized_new, d_labels->pLocalizedMPMLabel, pset);

      constParticleVariable<double> pMass, pTemperature;
      old_dw->get(pMass, d_labels->pMassLabel, pset);
      old_dw->get(pTemperature, d_labels->pTemperatureLabel, pset);

      constParticleVariable<Point> pX;
      old_dw->get(pX, d_labels->pXLabel, pset);

      constParticleVariable<Vector> pVelocity, pDisp;
      old_dw->get(pVelocity, d_labels->pVelocityLabel, pset);
      old_dw->get(pDisp, d_labels->pDispLabel, pset);

      constParticleVariable<Matrix3> pDefGrad_new;
      new_dw->get(pDefGrad_new, d_labels->pDefGradLabel_preReloc, pset);

      // Allocate updated particle variables
      ParticleVariable<double> pMass_new, pVolume_new, pTemperature_new, pTempPrevious;
      new_dw->getModifiable(pVolume_new, d_labels->pVolumeLabel_preReloc, pset);
      new_dw->allocateAndPut(pMass_new, d_labels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pTemperature_new, d_labels->pTemperatureLabel_preReloc, pset);
      new_dw->allocateAndPut(pTempPrevious, d_labels->pTempPreviousLabel_preReloc, pset);

      ParticleVariable<Point> pX_new, pXX;
      new_dw->allocateAndPut(pX_new, d_labels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pXX, d_labels->pXXLabel, pset);

      ParticleVariable<Vector> pVelocity_new, pDisp_new;
      new_dw->allocateAndPut(pVelocity_new, d_labels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pDisp_new, d_labels->pDispLabel_preReloc, pset);

      constNCVariable<Vector> gVelocity, gAcceleration;
      new_dw->get(gVelocity, d_labels->gVelocityLabel, matID, patch, Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gAcceleration, d_labels->gAccelerationLabel, matID, patch, Ghost::AroundCells,
                  d_numGhostParticles);

      ParticleSubset* delset = scinew ParticleSubset(0, matID, patch);

      for (auto particle : *pset) {

        interpolator->findCellAndWeights(pX[particle], influenceNodes, S_ip_av,
                                         pSize[particle], pDefGrad_new[particle]);

        Vector pVelocity_old(0.0, 0.0, 0.0);
        Vector pVel_new(0.0, 0.0, 0.0);
        Vector pAcceleration_old(0.0, 0.0, 0.0);
        for (int k = 0; k < num_influence_nodes; k++) {
          IntVector node = influenceNodes[k];
          Vector gVelocity_new = gVelocity[node] + gAcceleration[node] * delT;
          pVelocity_old += gVelocity[node] * S_ip_av[k];
          pVel_new += gVelocity_new * S_ip_av[k];
          pAcceleration_old += gAcceleration[node] * S_ip_av[k];
        }

        // Copy temperature
        pTemperature_new[particle] = pTemperature[particle];
        pTempPrevious[particle] = pTemperature[particle];

        // Update the particle's position and velocity
        pMass_new[particle] = pMass[particle];
        //auto v_p_mid = pVelocity[particle] + pAcceleration_old * (0.5 * delT);
        //auto v_p_mid = pVelocity_old + pAcceleration_old * (0.5 * delT);
        pX_new[particle] = pX[particle] + pVel_new * delT;
        //pX_new[particle] = pX[particle] + v_p_mid * delT;
        pDisp_new[particle] = pDisp[particle] + (pX_new[particle] - pX[particle]);
        pXX[particle] = pX[particle] + pDisp_new[particle];
        pVelocity_new[particle] = pVelocity[particle] + pAcceleration_old * delT;
        //pVelocity_new[particle] = pVelocity_old + pAcceleration_old * delT;

        /*
        if (!std::isfinite(pVelocity_new[particle].length())) {
          std::cout << "After interpolateToParticles: " 
                    << " material = " << m << " particle = " << particle
                    << " pX_old = " << pX[particle]
                    << " pX_new = " << pX_new[particle]
                    << " pAcc_old = " << pAcceleration_old
                    << " pVel_new = " << pVelocity_new[particle] << "\n";
        }
        */

        if (pParticleID[particle] == 30066212865) {
          std::cout << "v_p = " << pVelocity[particle]
                    << " v_p (interp) = " << pVelocity_old
                    << " diff = " << pVelocity[particle] - pVelocity_old << "\n";
        }

        ke += .5 * pMass[particle] * pVelocity_new[particle].length2();
        centerOfMass = centerOfMass + (pX_new[particle] * pMass[particle]).asVector();
        totalMom += pVelocity_new[particle] * pMass[particle];
        totalmass += pMass_new[particle];
        partvoldef += pVolume_new[particle];

        // Flag particles for deletion
        if ((pMass_new[particle] <= d_flags->d_minPartMass) ||
            (pLocalized_new[particle] == -999)) {
          delset->addParticle(particle);
          std::cout << "\n Warning: Added to delset: material = " << m << " particle = " << particle
                    << " pMass_new = " << pMass_new[particle]
                    << " pLocalized_new = " << pLocalized_new[particle] << "\n";
        }

        double vel_new = pVelocity_new[particle].length();
        if (vel_new > d_flags->d_maxVel) {
          double vel_old = pVelocity[particle].length();
          if (d_flags->d_deleteRogueParticles) {
            delset->addParticle(particle);
            std::cout << "\n Warning: particle " << pParticleID[particle]
                 << " hit speed ceiling #1 (" << vel_new << ">" << d_flags->d_maxVel
                 << "). Deleting particle.\n";
          } else {
            if (vel_new >= vel_old) {
              double scale_factor =  (d_flags->d_maxVel * 0.9) / vel_new;
              pVelocity_new[particle] *= scale_factor;
              std::cout << "\n Warning: particle " << pParticleID[particle]
                   << " hit speed ceiling #1 (" << vel_new << ">" << d_flags->d_maxVel
                   << ") vel_old = " << vel_old << ". Modifying particle velocity "
                      "accordingly.\n";
            }
          }
        }

      } // End loop over particles

      new_dw->deleteParticles(delset);

      //__________________________________
      //  particle debugging label-- carry forward
      if (d_flags->d_withColor) {
        constParticleVariable<double> pColor;
        ParticleVariable<double> pColor_new;
        old_dw->get(pColor, d_labels->pColorLabel, pset);
        new_dw->allocateAndPut(pColor_new, d_labels->pColorLabel_preReloc,
                               pset);
        pColor_new.copyData(pColor);
      }
    }

    // DON'T MOVE THESE!!!
    //__________________________________
    //  reduction variables
    if (d_flags->d_reductionVars->mass) {
      new_dw->put(sum_vartype(totalmass), d_labels->TotalMassLabel);
    }
    if (d_flags->d_reductionVars->volDeformed) {
      new_dw->put(sum_vartype(partvoldef), d_labels->TotalVolumeDeformedLabel);
    }
    if (d_flags->d_reductionVars->momentum) {
      new_dw->put(sumvec_vartype(totalMom), d_labels->TotalMomentumLabel);
    }
    if (d_flags->d_reductionVars->KE) {
      new_dw->put(sum_vartype(ke), d_labels->KineticEnergyLabel);
    }
    if (d_flags->d_reductionVars->centerOfMass) {
      new_dw->put(sumvec_vartype(centerOfMass), d_labels->CenterOfMassPositionLabel);
    }

    // cout << "Solid momentum after advection = " << totalMom << endl;
  }
}


/*!----------------------------------------------------------------------
 * scheduleComputeParticleScaleFactor
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleComputeParticleScaleFactor(SchedulerP& sched,
                                             const PatchSet* patches,
                                             const MaterialSet* matls)

{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels()))
    return;

  printSchedule(patches, cout_doing, "MPM::scheduleComputeParticleScaleFactor");

  Task* t = scinew Task("UofU_MPM::computeParticleScaleFactor", this,
                        &UofU_MPM::computeParticleScaleFactor);

  t->requires(Task::OldDW, d_labels->pSizeLabel, Ghost::None);
  t->requires(Task::NewDW, d_labels->pDefGradLabel_preReloc, Ghost::None);
  t->computes(d_labels->pScaleFactorLabel_preReloc);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * computeParticleScaleFactor
 *   This task computes the particles initial physical size, to be used
 *   in scaling particles for the deformed particle vis feature
 *-----------------------------------------------------------------------*/
void
UofU_MPM::computeParticleScaleFactor(const ProcessorGroup*,
                                     const PatchSubset* patches,
                                     const MaterialSubset*,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleScaleFactor");

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Matrix3> pScaleFactor;
      old_dw->get(pSize, d_labels->pSizeLabel, pset);
      new_dw->allocateAndPut(pScaleFactor, d_labels->pScaleFactorLabel_preReloc,
                             pset);

      if (d_dataArchiver->isOutputTimestep()) {
        Vector dx = patch->dCell();

        if (d_flags->d_interpolatorType != "cpdi" &&
            d_flags->d_interpolatorType != "cpti") {
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pDefGrad, d_labels->pDefGradLabel_preReloc, pset);
          for (auto particle : *pset) {
            pScaleFactor[particle] =
              ((pDefGrad[particle] * pSize[particle]) *
               Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));
          }
        } else {
          for (auto particle : *pset) {
            pScaleFactor[particle] =
              (pSize[particle] * Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));

          } // for particles
        }
      } // isOutputTimestep
    }   // matls
  }     // patches
}

/*!----------------------------------------------------------------------
 * scheduleTotalParticleCount
 *   Diagnostic task: compute the total number of particles
 *-----------------------------------------------------------------------*/
void
UofU_MPM::scheduleTotalParticleCount(SchedulerP& sched, const PatchSet* patches,
                                     const MaterialSet* matls)
{
  if (!d_flags->doMPMOnLevel(getLevel(patches)->getIndex(),
                             getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("UofU_MPM::totalParticleCount", this,
                        &UofU_MPM::totalParticleCount);
  t->computes(d_labels->partCountLabel);

  sched->addTask(t, patches, matls);
}

/*!----------------------------------------------------------------------
 * totalParticleCount
 *   Diagnostic task: compute the total number of particles
 *-----------------------------------------------------------------------*/
void
UofU_MPM::totalParticleCount(const ProcessorGroup*, const PatchSubset* patches,
                             const MaterialSubset* matls, DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    long int totalParticles = 0;

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
      int matID = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      int numParticles = pset->end() - pset->begin();

      totalParticles += numParticles;
    }
    new_dw->put(sumlong_vartype(totalParticles), d_labels->partCountLabel);
  }
}

/* _____________________________________________________________________
   Purpose:   Set variables that are normally set during the initialization
   phase, but get wiped clean when you restart
   _____________________________________________________________________*/
void
UofU_MPM::scheduleRestartInitialize(const LevelP& level, SchedulerP& sched)
{
}

/*!----------------------------------------------------------------------
 * restartInitialize
 *-----------------------------------------------------------------------*/
void
UofU_MPM::restartInitialize()
{
}

/*!----------------------------------------------------------------------
 * Set particle default
 *-----------------------------------------------------------------------*/
template <typename T>
void
UofU_MPM::setParticleDefault(ParticleVariable<T>& pvar,
                             const VarLabel* label, ParticleSubset* pset,
                             DataWarehouse* new_dw, const T& val)
{
  new_dw->allocateAndPut(pvar, label, pset);
  for (auto particle : *pset) {
    pvar[particle] = val;
  }
}
