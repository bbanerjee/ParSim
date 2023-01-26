/*
 * The MIT License
 *
 * Copyright (c) 1997-2015 The University of Utah
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
#include <CCA/Components/MPM/AMRMPM.h>

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/Contact/ContactFactory.h>
#include <CCA/Components/MPM/Core/MPMBoundCond.h>
#include <CCA/Components/MPM/MMS/MMS.h>
#include <CCA/Components/MPM/ParticleCreator/ParticleCreator.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/ScalarFluxBC.h>

#include <CCA/Components/MPM/ReactionDiffusion/DiffusionModels/ScalarDiffusionModel.h>
#include <CCA/Components/MPM/ReactionDiffusion/ScalarDiffusionTasks.h>

#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Regridder.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Geometry/Point.h>
#include <Core/Geometry/Vector.h>
#include <Core/GeometryPiece/FileGeometryPiece.h>
#include <Core/GeometryPiece/GeometryObject.h>
#include <Core/GeometryPiece/GeometryPieceFactory.h>
#include <Core/GeometryPiece/UnionGeometryPiece.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/AMR_CoarsenRefine.h>
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
#include <Core/Math/MinMax.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <fstream>
#include <iostream>
#include <sstream>

using namespace Uintah;

static DebugStream cout_doing("AMRMPM_cout", false);
static DebugStream amr_doing("AMRMPM_amr", false);

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "AMRMPM_cout:+,AMRMPM_amr:+"
//  bash     : export SCI_DEBUG="AMRMPM_cout:+,AMRMPM_amr:+"
//  default is OFF
// #define USE_DEBUG_TASK
// #define DEBUG_VEL
// #define DEBUG_ACC
#undef CBDI_FLUXBCS

//__________________________________
//   TODO:
// - We only need to compute ZOI when the grid changes not every timestep
//
// - InterpolateParticlesToGrid_CFI()  Need to account for gimp when getting
// particles on coarse level.
//
// - scheduleTimeAdvance:  Do we need to schedule each task in a separate level
// loop?  I suspect that we only need
//                         to in the CFI tasks
//
//  What is going on in refine & coarsen
//  To Test:
//    Symetric BC
//    compilicated grids
//    multi processors
//
//  Need to Add gimp interpolation

AMRMPM::AMRMPM(const ProcessorGroup* myworld,
               const MaterialManagerP& mat_manager)
  : SerialMPM(myworld, mat_manager)
{
  d_mpmFlags->d_minGridLevel = 0;
  d_mpmFlags->d_maxGridLevel = 1000;

  pDbgLabel =
    VarLabel::create("p.dbg", ParticleVariable<double>::getTypeDescription());
  gSumSLabel =
    VarLabel::create("g.sum_S", NCVariable<double>::getTypeDescription());
  RefineFlagXMaxLabel =
    VarLabel::create("RefFlagXMax", max_vartype::getTypeDescription());
  RefineFlagXMinLabel =
    VarLabel::create("RefFlagXMin", min_vartype::getTypeDescription());
  RefineFlagYMaxLabel =
    VarLabel::create("RefFlagYMax", max_vartype::getTypeDescription());
  RefineFlagYMinLabel =
    VarLabel::create("RefFlagYMin", min_vartype::getTypeDescription());
  RefineFlagZMaxLabel =
    VarLabel::create("RefFlagZMax", max_vartype::getTypeDescription());
  RefineFlagZMinLabel =
    VarLabel::create("RefFlagZMin", min_vartype::getTypeDescription());

  d_amrmpmLabels = std::make_unique<AMRMPMLabel>();

  d_oneMaterial = std::make_unique<MaterialSubset>();
  d_oneMaterial->add(0);
  d_oneMaterial->addReference();
}

AMRMPM::~AMRMPM()
{
  VarLabel::destroy(pDbgLabel);
  VarLabel::destroy(gSumSLabel);
  VarLabel::destroy(RefineFlagXMaxLabel);
  VarLabel::destroy(RefineFlagYMaxLabel);
  VarLabel::destroy(RefineFlagZMaxLabel);
  VarLabel::destroy(RefineFlagXMinLabel);
  VarLabel::destroy(RefineFlagYMinLabel);
  VarLabel::destroy(RefineFlagZMinLabel);

  d_oneMaterial->removeReference();

  d_refineGeomObjs.clear();
}

void
AMRMPM::problemSetup(const ProblemSpecP& prob_spec,
                     const ProblemSpecP& restart_prob_spec,
                     GridP& grid)
{
  cout_doing << "Doing problemSetup\t\t\t\t\t AMRMPM" << std::endl;

  d_scheduler->setPositionVar(d_mpmLabels->pXLabel);

  ProblemSpecP mat_ps = 0;
  if (restart_prob_spec) {
    mat_ps      = restart_prob_spec;
    d_isRestart = true;
  } else {
    mat_ps = prob_spec;
  }

  ProblemSpecP mpm_soln_ps = mat_ps->findBlock("MPM");
  if (!mpm_soln_ps) {
    std::ostringstream warn;
    warn << "ERROR:MPM:\n missing MPM section in the AMRMPM input file\n";
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  // Read all MPM flags (look in MPMFlags.cc)
  d_mpmFlags->readMPMFlags(mat_ps, d_output);
  if (d_mpmFlags->d_integratorType == "implicit") {
    throw ProblemSetupException("Can't use implicit integration with AMRMPM",
                                __FILE__,
                                __LINE__);
  }

  //__________________________________
  //  Read in the AMR section
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  ProblemSpecP mpm_ps{ nullptr };
  if (amr_ps) {
    mpm_ps = amr_ps->findBlock("MPM");

    d_mpmFlags->d_AMR = true;
  } else {
    string warn;
    warn = "\n INPUT FILE ERROR:\n <AMR>  block not found in input file \n";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  if (!mpm_ps) {
    string warn;
    warn =
      "\n INPUT FILE ERROR:\n <MPM>  block not found inside of <AMR> block \n";
    throw ProblemSetupException(warn, __FILE__, __LINE__);
  }

  // override CFI_interpolator
  d_CFI_interpolator = d_mpmFlags->d_interpolatorType;
  mpm_ps->get("CFI_interpolator", d_CFI_interpolator);

  if (d_CFI_interpolator != d_mpmFlags->d_interpolatorType) {
    proc0cout << "______________________________________________________\n"
              << "          AMRMPM:  WARNING\n"
              << "  The particle to grid interpolator at the CFI is ("
              << d_CFI_interpolator
              << "), however the over rest of the domain it is: "
              << flags->d_interpolatorType
              << "\n___________________________________________________________"
                 "___________"
              << std::endl;
  }

  // bulletproofing
  int maxLevel = grid->numLevels();
  for (int L = 0; L < maxLevel; L++) {
    LevelP level = grid->getLevel(L);
    IntVector ec = level->getExtraCells();

    if (ec != IntVector(1, 1, 1) &&
        flags->d_interpolatorType == "gimp") { // This should be generalized
      std::ostringstream
        msg; // Each interpolator should know how many EC needed.
      msg << "\n AMRMPM ERROR:\n The number of extraCells on level ("
          << level->getIndex()
          << ") is not equal to [1,1,1] required for the GIMP particle "
             "interpolator";
      throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
    }
  }

  // Get refinement criteria and regions
  ProblemSpecP refine_ps = mpm_ps->findBlock("Refinement_Criteria_Thresholds");
  if (refine_ps) {
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
        std::cerr << "WARNING: AMRMPM.cc: stringstream failed...\n";
      }

      size_t numMatls = m_materialManager->getNumMatls();

      // if using "all" matls
      if (matl == "all") {
        for (size_t m = 0; m < numMatls; m++) {
          data.matl = m;
          d_thresholdVars.push_back(data);
        }
      } else {
        d_thresholdVars.push_back(data);
      }
    }
  } // end if refine_ps

#if 0
  ProblemSpecP refine_ps = mpm_ps->findBlock("Refine_Regions");
  if (refine_ps) {
    // Read in the refined regions geometry objects
    int piece_num = 0;
    std::list<GeometryObject::DataItem> geom_obj_data;
    geom_obj_data.push_back(
      GeometryObject::DataItem("level", GeometryObject::Integer));

    for (ProblemSpecP geom_obj_ps = refine_ps->findBlock("geom_object");
         geom_obj_ps != 0;
         geom_obj_ps = geom_obj_ps->findNextBlock("geom_object")) {

      std::vector<GeometryPieceP> pieces;
      GeometryPieceFactory::create(geom_obj_ps, pieces);

      GeometryPieceP mainpiece;
      if (pieces.size() == 0) {
        throw ParameterNotFound("No piece specified in geom_object",
                                __FILE__,
                                __LINE__);
      } else if (pieces.size() > 1) {
        mainpiece = scinew UnionGeometryPiece(pieces);
      } else {
        mainpiece = pieces[0];
      }
      piece_num++;
      d_refineGeomObjs.push_back(
        scinew GeometryObject(mainpiece, geom_obj_ps, geom_obj_data));
    }
  } // if(refine_ps)
#endif

  //  bulletproofing
  if (!isLockstepAMR()) {
    std::ostringstream msg;
    msg << "\n ERROR: You must add \n"
        << " <useLockStep> true </useLockStep> \n"
        << " inside of the <AMR> section. \n";
    throw ProblemSetupException(msg.str(), __FILE__, __LINE__);
  }

  if (d_mpmFlags->d_8or27 == 8) {
    d_numGhostParticles = 1;
    d_numGhostNodes     = 1;
  } else if (d_mpmFlags->d_8or27 == 27 || d_mpmFlags->d_8or27 == 64) {
    d_numGhostParticles = 2;
    d_numGhostNodes     = 2;
  }

  MPMPhysicalBCFactory::create(mat_ps, grid, d_mpmFlags.get());

  contactModel = ContactFactory::create(UintahParallelComponent::d_myworld,
                                        mat_ps,
                                        d_materialManager,
                                        d_mpmLabels.get(),
                                        d_mpmFlags.get());

  // Determine extents for coarser level particle data
  // Linear Interpolation:  1 layer of coarse level cells
  // Gimp Interpolation:    2 layers

  /*`==========TESTING==========*/
  d_nPaddingCells_Coarse = 1;
  //  d_numGhostParticles = 1;
  /*===========TESTING==========`*/

  setParticleGhostLayer(Ghost::AroundNodes, d_numGhostParticles);

  materialProblemSetup(mat_ps, d_mpmFlags.get(), d_isRestart);

  // Create deformation gradient computer
  d_defGradComputer =
    std::make_unique<DeformationGradientComputer>(d_materialManager,
                                                  d_mpmLabels.get(),
                                                  d_mpmFlags.get());

  // Scalar diffusion
  d_diffusionTasks = std::make_unique<ScalarDiffusionTasks>(mat_ps,
                                                            d_materialManager,
                                                            d_mpmLabels.get(),
                                                            d_mpmFlags.get());
}

void
AMRMPM::outputProblemSpec(ProblemSpecP& root_ps)
{
  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP flags_ps = root->appendChild("MPM");
  d_mpmFlags->outputProblemSpec(flags_ps);

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

  d_diffusionTasks->outputProblemSpec(mps_ps);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps   = physical_bc_ps->appendChild("MPM");
  for (auto& bc : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    bc->outputProblemSpec(mpm_ph_bc_ps);
  }
}

void
AMRMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{

  if (d_mpmFlags->doMPMOnLevel(level->getIndex(),
                               level->getGrid()->numLevels())) {
    proc0cout << "doMPMOnLevel = " << level->getIndex() << std::endl;
  } else {
    proc0cout << "DontDoMPMOnLevel = " << level->getIndex() << std::endl;
  }

  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  Task* t = scinew Task("AMRMPM::actuallyInitialize",
                        this,
                        &AMRMPM::actuallyInitialize);

  t->computes(d_mpmLabels->partCountLabel);
  t->computes(d_mpmLabels->pXLabel);
  t->computes(d_mpmLabels->pDispLabel);
  t->computes(d_mpmLabels->pFiberDirLabel);
  t->computes(d_mpmLabels->pMassLabel);
  t->computes(d_mpmLabels->pVolumeLabel);
  t->computes(d_mpmLabels->pTemperatureLabel);
  t->computes(d_mpmLabels->pTempPreviousLabel); // for therma  stress analysis
  t->computes(d_mpmLabels->pdTdtLabel);
  t->computes(d_mpmLabels->pVelocityLabel);
  t->computes(d_mpmLabels->pAccelerationLabel);
  t->computes(d_mpmLabels->pExternalForceLabel);
  t->computes(d_mpmLabels->pParticleIDLabel);
  t->computes(d_mpmLabels->pStressLabel);
  t->computes(d_mpmLabels->pTemperatureGradientLabel);
  t->computes(d_mpmLabels->pSizeLabel);
  t->computes(d_mpmLabels->pRefinedLabel);
  t->computes(d_mpmLabels->pLastLevelLabel);
  t->computes(d_mpmLabels->delTLabel, level.get_rep());
  t->computes(d_mpmLabels->pCellNAPIDLabel, d_oneMaterial.get());
  t->computes(d_mpmLabels->NC_CCweightLabel, d_oneMaterial.get());

  // Needed for switch from explicit to implicit MPM
  t->computes(d_mpmLabels->pExternalHeatFluxLabel);

  // For friction contact
  t->computes(d_mpmLabels->pSurfLabel);

  if (!d_mpmFlags->d_doGridReset) {
    t->computes(d_mpmLabels->gDisplacementLabel);
  }

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

  if (d_mpmFlags->d_artificialViscosity) {
    t->computes(d_mpmLabels->p_qLabel);
  }

  // Scalar diffusion
  d_diffusionTasks->scheduleInitialize(t);

  const PatchSet* patches = level->eachPatch();

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    // Add vel grad/def grad computes
    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    cm->addInitialComputesAndRequires(t, mpm_matl, patches);

    // Scalar diffusion
    d_diffusionTasks->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  // Add initialization of body force and coriolis importance terms
  // These are initialize to zero in ParticleCreator
  t->computes(d_mpmLabels->pCoriolisImportanceLabel);
  t->computes(d_mpmLabels->pBodyForceAccLabel);

  sched->addTask(t, level->eachPatch(), d_materialManager->allMaterials("MPM"));

  if (level->getIndex() == 0) {
    schedulePrintParticleCount(level, sched);
  }

  if (d_mpmFlags->d_useLoadCurves && !d_mpmFlags->d_doScalarDiffusion) {
    // Schedule the initialization of pressure BCs per particle
    std::cout << "Pressure load curves are untested for multiple levels"
              << std::endl;
    scheduleInitializePressureBCs(level, sched);
  }

  // Scalar diffusion flux particle BCs
  d_diffusionTasks->scheduleInitializeFluxBCs(level, sched);
}

void
AMRMPM::schedulePrintParticleCount(const LevelP& level, SchedulerP& sched)
{
  Task* t = scinew Task("AMRMPM::printParticleCount",
                        this,
                        &AMRMPM::printParticleCount);
  t->requires(Task::NewDW, d_mpmLabels->partCountLabel);
  t->setType(Task::OncePerProc);

  sched->addTask(t,
                 sched->getLoadBalancer()->getPerProcessorPatchSet(level),
                 d_materialManager->allMaterials("MPM"));
}

void
AMRMPM::scheduleComputeStableTimestep(const LevelP&, SchedulerP&)
{
  // Nothing to do here - delt is computed as a by-product of the
  // constitutive model
}

void
AMRMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  if (level->getIndex() > 0) { // only schedule once
    return;
  }

  const MaterialSet* matls = d_materialManager->allMaterials("MPM");
  GridP grid               = level->getGrid();
  int maxLevels            = grid->numLevels();

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    schedulePartitionOfUnity(sched, patches, matls);
    scheduleComputeZoneOfInfluence(sched, patches, matls);
    scheduleComputeParticleBodyForce(sched, patches, matls);
    scheduleComputeCurrentParticleSize(sched, patches, matls);
    scheduleApplyExternalLoads(sched, patches, matls);
    d_diffusionTasks->scheduleApplyExternalScalarFlux(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleInterpolateParticlesToGrid(sched, patches, matls);
    // Need to add a task to do the reductions on the max hydro stress - JG ???
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleInterpolateParticlesToGrid_CFI(sched, patches, matls);
  }

#ifdef USE_DEBUG_TASK
  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleDebug_CFI(sched, patches, matls);
  }
#endif

  for (int l = 0; l < maxLevels - 1; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleCoarsenNodalData_CFI(sched, patches, matls, coarsenData);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleNormalizeNodalVelTempConc(sched, patches, matls);
    scheduleComputeNormals(sched, patches, matls);
    scheduleFindSurfaceParticles(sched, patches, matls);
    scheduleComputeLogisticRegression(sched, patches, matls);
    scheduleMomentumExchangeInterpolated(sched, patches, matls);
    d_diffusionTasks->scheduleConcInterpolated(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleComputeInternalForce(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleComputeInternalForce_CFI(sched, patches, matls);
  }

  d_diffusionTasks->scheduleComputeAMR(sched, patches, matls);

  for (int l = 0; l < maxLevels - 1; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleCoarsenNodalData_CFI2(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleComputeAndIntegrateAcceleration(sched, patches, matls);
    scheduleMomentumExchangeIntegrated(sched, patches, matls);
    scheduleSetGridBoundaryConditions(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleComputeLAndF(sched, patches, matls);
  }

#if 0
  // zero the nodal data at the CFI on the coarse level 
  for (int l = 0; l < maxLevels-1; l++) {
    const LevelP& level = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleCoarsenNodalData_CFI( sched, patches, matls, zeroData);
  }
#endif

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleInterpolateToParticlesAndUpdate(sched, patches, matls);
  }

#if 0
  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleInterpolateToParticlesAndUpdate_CFI(sched, patches, matls);
  }
#endif

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleComputeStressTensor(sched, patches, matls);
  }

  if (d_mpmFlags->d_computeScaleFactor) {
    for (int l = 0; l < maxLevels; l++) {
      const LevelP& level     = grid->getLevel(l);
      const PatchSet* patches = level->eachPatch();
      scheduleComputeParticleScaleFactor(sched, patches, matls);
    }
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleFinalParticleUpdate(sched, patches, matls);
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    if (d_mpmFlags->d_refineParticles) {
      scheduleAddParticles(sched, patches, matls);
    }
  }

  for (int l = 0; l < maxLevels; l++) {
    const LevelP& level     = grid->getLevel(l);
    const PatchSet* patches = level->eachPatch();
    scheduleReduceFlagsExtents(sched, patches, matls);
  }
}

void
AMRMPM::scheduleFinalizeTimestep(const LevelP& level, SchedulerP& sched)
{

  const PatchSet* patches = level->eachPatch();

  if (level->getIndex() == 0) {
    const MaterialSet* matls = d_materialManager->allMaterials("MPM");
    sched->scheduleParticleRelocation(level,
                                      d_mpmLabels->pXLabel_preReloc,
                                      d_particleState_preReloc,
                                      d_mpmLabels->pXLabel,
                                      d_particleState,
                                      d_mpmLabels->pParticleIDLabel,
                                      matls);
  }
}

void
AMRMPM::scheduleAnalysis(const LevelP& level, SchedulerP& sched)
{
  const PatchSet* patches = level->eachPatch();

  scheduleCountParticles(patches, sched);
}

void
AMRMPM::schedulePartitionOfUnity(SchedulerP& sched,
                                 const PatchSet* patches,
                                 const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "AMRMPM::schedulePartitionOfUnity");
  Task* t =
    scinew Task("AMRMPM::partitionOfUnity", this, &AMRMPM::partitionOfUnity);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, Ghost::None);

  // Carry forward and update pSize if particles change levels
  t->requires(Task::OldDW, d_mpmLabels->pSizeLabel, Ghost::None);
  t->requires(Task::OldDW, d_mpmLabels->pLastLevelLabel, Ghost::None);

  t->computes(d_mpmLabels->pSizeLabel_preReloc);
  t->computes(d_mpmLabels->pLastLevelLabel_preReloc);
  t->computes(d_mpmLabels->pPartitionUnityLabel);
  t->computes(d_mpmLabels->MPMRefineCellLabel, d_oneMaterial.get());

  sched->addTask(t, patches, matls);
}

void
AMRMPM::scheduleComputeZoneOfInfluence(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  int L_indx         = level->getIndex();

  if (L_indx > 0) {

    printSchedule(patches,
                  cout_doing,
                  "AMRMPM::scheduleComputeZoneOfInfluence");
    Task* t = scinew Task("AMRMPM::computeZoneOfInfluence",
                          this,
                          &AMRMPM::computeZoneOfInfluence);

    t->computes(d_mpmLabels->gZOILabel, d_oneMaterial);

    sched->addTask(t, patches, matls);
  }
}

void
AMRMPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{

  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleInterpolateParticlesToGrid");

  Task* t = scinew Task("AMRMPM::interpolateParticlesToGrid",
                        this,
                        &AMRMPM::interpolateParticlesToGrid);

  Ghost::GhostType gan = Ghost::AroundNodes;
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
  if (d_mpmFlags->d_GEVelProj) {
    t->requires(Task::OldDW,
                d_mpmLabels->pVelGradLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }
  t->requires(Task::OldDW,
              d_mpmLabels->pXLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pExtForceLabel_preReloc,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pTemperatureLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::NewDW,
              d_mpmLabels->pSizeLabel_preReloc,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);

  t->computes(d_mpmLabels->gMassLabel);
  t->computes(d_mpmLabels->gVolumeLabel);
  t->computes(d_mpmLabels->gVelocityLabel);
  t->computes(d_mpmLabels->gTemperatureLabel);
  t->computes(d_mpmLabels->gTemperatureRateLabel);
  t->computes(d_mpmLabels->gExternalForceLabel);

  d_diffusionTasks->scheduleInterpolateParticlesToGrid(t, d_numGhostParticles);

  sched->addTask(t, patches, matls);
}

//  You need particle data from the coarse levels at the CFI on the fine level
void
AMRMPM::scheduleInterpolateParticlesToGrid_CFI(SchedulerP& sched,
                                               const PatchSet* patches,
                                               const MaterialSet* matls)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx             = fineLevel->getIndex();

  if (L_indx > 0) {
    printSchedule(patches,
                  cout_doing,
                  "AMRMPM::scheduleInterpolateParticlesToGrid_CFI");

    Task* t = nullptr;
    if (d_CFI_interpolator == "gimp") {

      t = scinew Task("AMRMPM::interpolateParticlesToGrid_CFI_GIMP",
                      this,
                      &AMRMPM::interpolateParticlesToGrid_CFI_GIMP);
    } else {
      t = scinew Task("AMRMPM::interpolateParticlesToGrid_CFI",
                      this,
                      &AMRMPM::interpolateParticlesToGrid_CFI);
    }

    /*`==========TESTING==========*/
    // Linear 1 coarse Level cells:
    // Gimp:  2 coarse level cells:
    int npc = d_nPaddingCells_Coarse;
    /*===========TESTING==========`*/

#define allPatches 0
#define allMatls 0
    //__________________________________
    //  Note: were using nPaddingCells to extract the region of coarse level
    // particles around every fine patch.   Technically, these are ghost
    // cells but somehow it works.
    t->requires(Task::NewDW,
                d_mpmLabels->gZOILabel,
                d_oneMaterial.get(),
                Ghost::None,
                0);
    t->requires(Task::OldDW,
                d_mpmLabels->pMassLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pVolumeLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pVelocityLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pXLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::NewDW,
                d_mpmLabels->pExtForceLabel_preReloc,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pTemperatureLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pDefGradLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);

    t->modifies(d_mpmLabels->gMassLabel);
    t->modifies(d_mpmLabels->gVolumeLabel);
    t->modifies(d_mpmLabels->gVelocityLabel);
    t->modifies(d_mpmLabels->gTemperatureLabel);
    t->modifies(d_mpmLabels->gExternalForceLabel);

    d_diffusionTasks->scheduleInterpolateParticlesToGrid_CFI(t, npc);

    sched->addTask(t, patches, matls);
  }
}

//______________________________________________________________________
//  This task does one of two operations on the coarse nodes along
//  the coarse fine interface.  The input parameter "flag" determines
//  which.
//  Coarsen:  copy fine patch node data to the coarse level at CFI
//  Zero:     zero the coarse level nodal data directly under the fine level
void
AMRMPM::scheduleCoarsenNodalData_CFI(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls,
                                     const coarsenFlag flag)
{
  string txt = "(zero)";
  if (flag == coarsenData) {
    txt = "(coarsen)";
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleCoarsenNodalData_CFI" + txt);

  Task* t = scinew Task("AMRMPM::coarsenNodalData_CFI",
                        this,
                        &AMRMPM::coarsenNodalData_CFI,
                        flag);

  Ghost::GhostType gn         = Ghost::None;
#define allPatches 0
#define allMatls 0

  t->requires(Task::NewDW,
              d_mpmLabels->gMassLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);
  t->requires(Task::NewDW,
              d_mpmLabels->gVolumeLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);
  t->requires(Task::NewDW,
              d_mpmLabels->gVelocityLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);
  t->requires(Task::NewDW,
              d_mpmLabels->gTemperatureLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);
  t->requires(Task::NewDW,
              d_mpmLabels->gExternalForceLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);

  t->modifies(d_mpmLabels->gMassLabel);
  t->modifies(d_mpmLabels->gVolumeLabel);
  t->modifies(d_mpmLabels->gVelocityLabel);
  t->modifies(d_mpmLabels->gTemperatureLabel);
  t->modifies(d_mpmLabels->gExternalForceLabel);

  if (d_mpmFlags->d_doScalarDiffusion) {
    t->requires(Task::NewDW,
                d_mpmLabels->gConcentrationLabel,
                allPatches,
                Task::FineLevel,
                allMatls,
                Task::NormalDomain,
                gn,
                0);
    t->modifies(d_mpmLabels->gConcentrationLabel);
    t->requires(Task::NewDW,
                d_mpmLabels->gExternalScalarFluxLabel,
                allPatches,
                Task::FineLevel,
                allMatls,
                Task::NormalDomain,
                gn,
                0);
    t->modifies(d_mpmLabels->gExternalScalarFluxLabel);
  }

  if (flag == zeroData) {
    t->modifies(d_mpmLabels->gAccelerationLabel);
    t->modifies(d_mpmLabels->gVelocityStarLabel);
  }

  sched->addTask(t, patches, matls);
}

//______________________________________________________________________
//  This task copies fine patch node data to the coarse level at CFI
void
AMRMPM::scheduleCoarsenNodalData_CFI2(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "AMRMPM::scheduleCoarsenNodalData_CFI2");

  Task* t = scinew Task("AMRMPM::coarsenNodalData_CFI2",
                        this,
                        &AMRMPM::coarsenNodalData_CFI2);

  Ghost::GhostType gn         = Ghost::None;
  Task::MaterialDomainSpec ND = Task::NormalDomain;
#define allPatches 0
#define allMatls 0

  t->requires(Task::NewDW,
              d_mpmLabels->gMassLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);
  t->requires(Task::NewDW,
              d_mpmLabels->gInternalForceLabel,
              allPatches,
              Task::FineLevel,
              allMatls,
              Task::NormalDomain,
              gn,
              0);

  t->modifies(d_mpmLabels->gInternalForceLabel);
  if (d_mpmFlags->d_doScalarDiffusion) {
    t->requires(Task::NewDW,
                d_mpmLabels->diffusion->gConcentrationRate,
                allPatches,
                Task::FineLevel,
                allMatls,
                Task::NormalDomain,
                gn,
                0);
    t->modifies(d_mpmLabels->diffusion->gConcentrationRate);
  }

  sched->addTask(t, patches, matls);
}

//______________________________________________________________________
//  compute the nodal velocity and temperature after coarsening the fine
//  nodal data
void
AMRMPM::scheduleNormalizeNodalVelTempConc(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleNormalizeNodalVelTempConc");

  Task* t = scinew Task("AMRMPM::normalizeNodalVelTempConc",
                        this,
                        &AMRMPM::normalizeNodalVelTempConc);

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);

  t->modifies(d_mpmLabels->gVelocityLabel);
  t->modifies(d_mpmLabels->gTemperatureLabel);
  if (d_mpmFlags->d_doScalarDiffusion) {
    t->modifies(d_mpmLabels->gConcentrationLabel);
    t->computes(d_mpmLabels->gConcentrationNoBCLabel);
    t->modifies(d_mpmLabels->diffusion->gHydrostaticStress);
  }

  sched->addTask(t, patches, matls);
}

//______________________________________________________________________
//
/////////////////////////////////////////////////////////////////////////
/*!  **WARNING** In addition to the stresses and deformations, the internal
 *               heat rate in the particles (pdTdtLabel)
 *               is computed here */
/////////////////////////////////////////////////////////////////////////
void
AMRMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  // Schedule compute of the deformation gradient
  scheduleComputeDeformationGradient(sched, patches, matls);

  printSchedule(patches, cout_doing, "AMRMPM::scheduleComputeStressTensor");

  int numMatls = d_materialManager->getNumMaterials("MPM");
  Task* t      = scinew Task("AMRMPM::computeStressTensor",
                        this,
                        &AMRMPM::computeStressTensor);

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    const MaterialSubset* matlset = mpm_matl->thisMaterial();

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);

    t->computes(d_mpmLabels->p_qLabel_preReloc, matlset);
  }

  t->computes(d_mpmLabels->delTLabel, getLevel(patches));
  t->computes(d_mpmLabels->StrainEnergyLabel);

  sched->addTask(t, patches, matls);

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    scheduleComputeAccStrainEnergy(sched, patches, matls);
  }
}

//______________________________________________________________________
//
void
AMRMPM::scheduleUpdateErosionParameter(SchedulerP& sched,
                                       const PatchSet* patches,
                                       const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "AMRMPM::scheduleUpdateErosionParameter");

  Task* t = scinew Task("AMRMPM::updateErosionParameter",
                        this,
                        &AMRMPM::updateErosionParameter);

  int numMatls = d_materialManager->getNumMaterials("MPM");

  for (int m = 0; m < numMatls; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addRequiresDamageParameter(t, mpm_matl, patches);
  }

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
//
void
AMRMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "AMRMPM::scheduleComputeInternalForce");

  Task* t = scinew Task("AMRMPM::computeInternalForce",
                        this,
                        &AMRMPM::computeInternalForce);

  Ghost::GhostType gan   = Ghost::AroundNodes;
  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, d_mpmLabels->gVolumeLabel, gnone);
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
  t->requires(Task::NewDW,
              d_mpmLabels->pSizeLabel_preReloc,
              Ghost::AroundNodes,
              d_numGhostParticles);
  t->requires(Task::OldDW,
              d_mpmLabels->pDefGradLabel,
              Ghost::AroundNodes,
              d_numGhostParticles);
  if (d_mpmFlags->d_artificialViscosity) {
    t->requires(Task::OldDW,
                d_mpmLabels->p_qLabel,
                Ghost::AroundNodes,
                d_numGhostParticles);
  }

  t->computes(gSumSLabel);
  t->computes(d_mpmLabels->gInternalForceLabel);
  t->computes(d_mpmLabels->TotalVolumeDeformedLabel);
  t->computes(d_mpmLabels->gStressForSavingLabel);

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
//
void
AMRMPM::scheduleComputeInternalForce_CFI(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  const Level* fineLevel = getLevel(patches);
  int L_indx             = fineLevel->getIndex();

  if (L_indx > 0) {
    printSchedule(patches,
                  cout_doing,
                  "AMRMPM::scheduleComputeInternalForce_CFI");

    Task* t = scinew Task("AMRMPM::computeInternalForce_CFI",
                          this,
                          &AMRMPM::computeInternalForce_CFI);

    Ghost::GhostType gac        = Ghost::AroundCells;
    Task::MaterialDomainSpec ND = Task::NormalDomain;

    /*`==========TESTING==========*/
    // Linear 1 coarse Level cells:
    // Gimp:  2 coarse level cells:
    int npc = d_nPaddingCells_Coarse;
    /*===========TESTING==========`*/

#define allPatches 0
#define allMatls 0
    //__________________________________
    //  Note: were using nPaddingCells to extract the region of coarse level
    // particles around every fine patch.   Technically, these are ghost
    // cells but somehow it works.
    t->requires(Task::NewDW,
                d_mpmLabels->gZOILabel,
                d_oneMaterial,
                Ghost::None,
                0);
    t->requires(Task::OldDW,
                d_mpmLabels->pXLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pStressLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);
    t->requires(Task::OldDW,
                d_mpmLabels->pVolumeLabel,
                allPatches,
                Task::CoarseLevel,
                allMatls,
                Task::NormalDomain,
                Ghost::AroundCells,
                npc);

    if (d_mpmFlags->d_artificialViscosity) {
      t->requires(Task::OldDW,
                  d_mpmLabels->p_qLabel,
                  allPatches,
                  Task::CoarseLevel,
                  allMatls,
                  Task::NormalDomain,
                  Ghost::AroundCells,
                  npc);
    }

    t->modifies(gSumSLabel);
    t->modifies(d_mpmLabels->gInternalForceLabel);

    sched->addTask(t, patches, matls);
  }
}

//______________________________________________________________________
//
void
AMRMPM::scheduleComputeAndIntegrateAcceleration(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleComputeAndIntegrateAcceleration");

  Task* t = scinew Task("AMRMPM::computeAndIntegrateAcceleration",
                        this,
                        &AMRMPM::computeAndIntegrateAcceleration);

  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  t->requires(Task::NewDW, d_mpmLabels->gMassLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gInternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gExternalForceLabel, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  t->computes(d_mpmLabels->gVelocityStarLabel);
  t->computes(d_mpmLabels->gAccelerationLabel);

  // This stuff should probably go in its own task, but for expediency...JG
  if (d_mpmFlags->d_doScalarDiffusion) {
    t->requires(Task::NewDW, d_mpmLabels->gConcentrationNoBCLabel, Ghost::None);
    t->requires(Task::NewDW, d_mpmLabels->gConcentrationLabel, Ghost::None);
    t->requires(Task::NewDW,
                d_mpmLabels->gExternalScalarFluxLabel,
                Ghost::None);
    t->modifies(d_mpmLabels->diffusion->gConcentrationRate);
    t->computes(d_mpmLabels->gConcentrationStarLabel);
  }

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
//
void
AMRMPM::scheduleSetGridBoundaryConditions(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)

{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleSetGridBoundaryConditions");

  Task* t = scinew Task("AMRMPM::setGridBoundaryConditions",
                        this,
                        &AMRMPM::setGridBoundaryConditions);

  const MaterialSubset* mss = matls->getUnion();
  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  t->modifies(d_mpmLabels->gAccelerationLabel, mss);
  t->modifies(d_mpmLabels->gVelocityStarLabel, mss);
  t->requires(Task::NewDW, d_mpmLabels->gVelocityLabel, Ghost::None);

  if (!d_mpmFlags->d_doGridReset) {
    t->requires(Task::OldDW, d_mpmLabels->gDisplacementLabel, Ghost::None);
    t->computes(d_mpmLabels->gDisplacementLabel);
  }

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
//
void
AMRMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)

{
  const Level* level = getLevel(patches);
  if (!d_mpmFlags->doMPMOnLevel(level->getIndex(),
                                level->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleInterpolateToParticlesAndUpdate");

  Task* t = scinew Task("AMRMPM::interpolateToParticlesAndUpdate",
                        this,
                        &AMRMPM::interpolateToParticlesAndUpdate);

  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;

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
              d_mpmLabels->frictionalWorkLabel,
              Ghost::AroundCells,
              d_numGhostNodes);

  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pMassLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pParticleIDLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pTemperatureLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pVelocityLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pDispLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pSizeLabel_preReloc, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pVolumeLabel, gnone);
  t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, gnone);

  t->computes(d_mpmLabels->pDispLabel_preReloc);
  t->computes(d_mpmLabels->pVelocityLabel_preReloc);
  t->computes(d_mpmLabels->pAccelerationLabel_preReloc);
  t->computes(d_mpmLabels->pXLabel_preReloc);
  t->computes(d_mpmLabels->pParticleIDLabel_preReloc);
  t->computes(d_mpmLabels->pTemperatureLabel_preReloc);
  t->computes(d_mpmLabels->pTempPreviousLabel_preReloc); // for thermal stress
  t->computes(d_mpmLabels->pMassLabel_preReloc);
  t->computes(d_mpmLabels->pXXLabel);

  // Carry forward external heat flux for switch from explicit to implicit
  t->requires(Task::OldDW, d_mpmLabels->pExternalHeatFluxLabel, Ghost::None);
  t->computes(d_mpmLabels->pExternalHeatFluxLabel_preReloc);

  // Carry Forward particle refinement flag
  if (d_mpmFlags->d_refineParticles) {
    t->requires(Task::OldDW, d_mpmLabels->pRefinedLabel, Ghost::None);
    t->computes(d_mpmLabels->pRefinedLabel_preReloc);
  }

  if (d_mpmFlags->d_doScalarDiffusion) {
    t->requires(Task::OldDW, d_mpmLabels->diffusion->pConcentration, gnone);
    t->requires(Task::NewDW,
                d_mpmLabels->diffusion->gConcentrationRate,
                Ghost::AroundCells,
                d_numGhostNodes);

    t->computes(d_mpmLabels->diffusion->pConcentration_preReloc);
    t->computes(d_mpmLabels->pConcPreviousLabel_preReloc);
  }

  t->computes(d_mpmLabels->TotalMassLabel);
  t->computes(d_mpmLabels->KineticEnergyLabel);
  t->computes(d_mpmLabels->ThermalEnergyLabel);
  t->computes(d_mpmLabels->CenterOfMassPositionLabel);
  t->computes(d_mpmLabels->TotalMomentumLabel);

#ifndef USE_DEBUG_TASK
  // debugging scalar
  if (d_mpmFlags->d_withColor) {
    t->requires(Task::OldDW, d_mpmLabels->pColorLabel, gnone);
    t->computes(d_mpmLabels->pColorLabel_preReloc);
  }
#endif
  sched->addTask(t, patches, matls);
}

void
AMRMPM::scheduleComputeParticleScaleFactor(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)

{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches,
                cout_doing,
                "AMRMPM::scheduleComputeParticleScaleFactor");

  Task* t = scinew Task("AMRMPM::computeParticleScaleFactor",
                        this,
                        &AMRMPM::computeParticleScaleFactor);

  t->requires(Task::NewDW, d_mpmLabels->pSizeLabel_preReloc, Ghost::None);
  t->requires(Task::NewDW, d_mpmLabels->pDefGradLabel_preReloc, Ghost::None);
  t->computes(d_mpmLabels->pScaleFactorLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
AMRMPM::scheduleFinalParticleUpdate(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "AMRMPM::scheduleFinalParticleUpdate");

  Task* t = scinew Task("AMRMPM::finalParticleUpdate",
                        this,
                        &AMRMPM::finalParticleUpdate);

  t->requires(Task::OldDW, d_mpmLabels->delTLabel);

  Ghost::GhostType gnone = Ghost::None;
  t->requires(Task::NewDW, d_mpmLabels->pdTdtLabel, gnone);
  t->requires(Task::NewDW, d_mpmLabels->pMassLabel_preReloc, gnone);

  t->modifies(d_mpmLabels->pTemperatureLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void
AMRMPM::scheduleAddParticles(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  if (!d_mpmFlags->doMPMOnLevel(getLevel(patches)->getIndex(),
                                getLevel(patches)->getGrid()->numLevels())) {
    return;
  }

  printSchedule(patches, cout_doing, "AMRMPM::scheduleAddParticles");

  Task* t = scinew Task("AMRMPM::addParticles", this, &AMRMPM::addParticles);

  t->modifies(d_mpmLabels->pParticleIDLabel_preReloc);
  t->modifies(d_mpmLabels->pXLabel_preReloc);
  t->modifies(d_mpmLabels->pVolumeLabel_preReloc);
  t->modifies(d_mpmLabels->pVelocityLabel_preReloc);
  t->modifies(d_mpmLabels->pAccelerationLabel_preReloc);
  t->modifies(d_mpmLabels->pMassLabel_preReloc);
  t->modifies(d_mpmLabels->pSizeLabel_preReloc);
  t->modifies(d_mpmLabels->pDispLabel_preReloc);
  t->modifies(d_mpmLabels->pStressLabel_preReloc);
  if (d_mpmFlags->d_withColor) {
    t->modifies(d_mpmLabels->pColorLabel_preReloc);
  }
  if (d_mpmFlags->d_doScalarDiffusion) {
    t->modifies(d_mpmLabels->diffusion->pConcentration_preReloc);
    t->modifies(d_mpmLabels->pConcPreviousLabel_preReloc);
    t->modifies(d_mpmLabels->pGradConcentrationLabel_preReloc);
    t->modifies(d_mpmLabels->pExternalScalarFluxLabel);
  }
  if (d_mpmFlags->d_useLoadCurves) {
    t->modifies(d_mpmLabels->pLoadCurveIDLabel_preReloc);
  }
  t->modifies(d_mpmLabels->pExtForceLabel_preReloc);
  t->modifies(d_mpmLabels->pTemperatureLabel_preReloc);
  t->modifies(d_mpmLabels->pTempPreviousLabel_preReloc);
  t->modifies(d_mpmLabels->pDefGradLabel_preReloc);
  t->modifies(d_mpmLabels->pRefinedLabel_preReloc);
  t->modifies(d_mpmLabels->pScaleFactorLabel_preReloc);
  t->modifies(d_mpmLabels->pLastLevelLabel_preReloc);
  t->modifies(d_mpmLabels->pVelGradLabel_preReloc);
  t->modifies(d_mpmLabels->MPMRefineCellLabel, d_oneMaterial);

  // For body dorce + coriolis importance
  t->modifies(d_mpmLabels->pBodyForceAccLabel_preReloc);
  t->modifies(d_mpmLabels->pCoriolisImportanceLabel_preReloc);

  t->requires(Task::OldDW,
              d_mpmLabels->pCellNAPIDLabel,
              d_oneMaterial,
              Ghost::None);
  t->computes(d_mpmLabels->pCellNAPIDLabel, d_oneMaterial);

  sched->addTask(t, patches, matls);
}

void
AMRMPM::scheduleReduceFlagsExtents(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
{
  const Level* level = getLevel(patches);

  //  if( !level->hasFinerLevel() ){
  if (level->getIndex() > 0) {
    printSchedule(patches, cout_doing, "AMRMPM::scheduleReduceFlagsExtents");

    Task* t = scinew Task("AMRMPM::reduceFlagsExtents",
                          this,
                          &AMRMPM::reduceFlagsExtents);

    t->requires(Task::NewDW,
                d_mpmLabels->MPMRefineCellLabel,
                d_oneMaterial,
                Ghost::None);

    t->computes(RefineFlagXMaxLabel);
    t->computes(RefineFlagXMinLabel);
    t->computes(RefineFlagYMaxLabel);
    t->computes(RefineFlagYMinLabel);
    t->computes(RefineFlagZMaxLabel);
    t->computes(RefineFlagZMinLabel);

    sched->addTask(t, patches, matls);
  }
}

//______________________________________________________________________
////
void
AMRMPM::scheduleRefine(const PatchSet* patches, SchedulerP& sched)
{
  printSchedule(patches, cout_doing, "AMRMPM::scheduleRefine");
  Task* t = scinew Task("AMRMPM::refineGrid", this, &AMRMPM::refineGrid);

  t->computes(d_mpmLabels->pXLabel);
  t->computes(d_mpmLabels->pDispLabel);
  t->computes(d_mpmLabels->pMassLabel);
  t->computes(d_mpmLabels->pTemperatureLabel);
  t->computes(d_mpmLabels->pTempPreviousLabel); // for thermal  stress analysis
  t->computes(d_mpmLabels->pdTdtLabel);
  t->computes(d_mpmLabels->pVelocityLabel);
  t->computes(d_mpmLabels->pExternalForceLabel);
  t->computes(d_mpmLabels->pExternalScalarFluxLabel);
  t->computes(d_mpmLabels->pParticleIDLabel);
  // t->computes(d_mpmLabels->pDefGradLabel);
  t->computes(d_mpmLabels->pVolumeLabel);
  t->computes(d_mpmLabels->pStressLabel);
  if (d_mpmFlags->d_doScalarDiffusion) {
    t->computes(d_mpmLabels->diffusion->pConcentration);
    t->computes(d_mpmLabels->pConcPreviousLabel);
    t->computes(d_mpmLabels->pGradConcentrationLabel);
  }
  t->computes(d_mpmLabels->pLastLevelLabel);
  t->computes(d_mpmLabels->pRefinedLabel);
  t->computes(d_mpmLabels->pSizeLabel);
  // t->computes(d_mpmLabels->pVelGradLabel);
  t->computes(d_mpmLabels->pCellNAPIDLabel, d_oneMaterial);

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

  if (d_mpmFlags->d_artificialViscosity) {
    t->computes(d_mpmLabels->p_qLabel);
  }

  int numMPM = d_materialManager->getNumMaterials("MPM");
  for (int m = 0; m < numMPM; m++) {
    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    // Add requires and computes for vel grad/def grad
    d_defGradComputer->addComputesOnly(t, mpm_matl, patches);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }
  t->computes(d_mpmLabels->delTLabel, getLevel(patches));

  // For body dorce + coriolis importance
  t->computes(d_mpmLabels->pBodyForceAccLabel);
  t->computes(d_mpmLabels->pCoriolisImportanceLabel);

  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
}
//______________________________________________________________________
//
void
AMRMPM::scheduleRefineInterface(const LevelP& /*fineLevel*/,
                                SchedulerP& /*scheduler*/,
                                bool,
                                bool)
{
  // do nothing for now
}
//______________________________________________________________________
//
void
AMRMPM::scheduleCoarsen(const LevelP& coarseLevel, SchedulerP& sched)
{
  // Coarsening the refineCell data so that errorEstimate will have it
  // on all levels
  Ghost::GhostType gn = Ghost::None;

  Task* task = scinew Task("AMRMPM::coarsen", this, &AMRMPM::coarsen);

  Task::MaterialDomainSpec oims = Task::OutOfDomain; // outside of ice matlSet.
  const MaterialSet* all_matls  = d_materialManager->allMaterials();
  const PatchSet* patch_set     = coarseLevel->eachPatch();

  bool fat = true; // possibly (F)rom (A)nother (T)askgraph

  task->requires(Task::NewDW,
                 d_mpmLabels->MPMRefineCellLabel,
                 0,
                 Task::FineLevel,
                 d_oneMaterial,
                 oims,
                 gn,
                 0,
                 fat);

  task->requires(Task::NewDW, RefineFlagXMaxLabel);
  task->requires(Task::NewDW, RefineFlagXMinLabel);
  task->requires(Task::NewDW, RefineFlagYMaxLabel);
  task->requires(Task::NewDW, RefineFlagYMinLabel);
  task->requires(Task::NewDW, RefineFlagZMaxLabel);
  task->requires(Task::NewDW, RefineFlagZMinLabel);

  task->modifies(d_mpmLabels->MPMRefineCellLabel, d_oneMaterial, oims, fat);

  sched->addTask(task, patch_set, all_matls);
}

void
AMRMPM::coarsen(const ProcessorGroup*,
                const PatchSubset* patches,
                const MaterialSubset* matls,
                DataWarehouse*,
                DataWarehouse* new_dw)
{
  const Level* coarseLevel = getLevel(patches);
  const Level* fineLevel   = coarseLevel->getFinerLevel().get_rep();
  GridP grid               = coarseLevel->getGrid();
  int numLevels            = grid->numLevels();
  IntVector RR             = fineLevel->getRefinementRatio();

  for (int p = 0; p < patches->size(); p++) {
    const Patch* coarsePatch = patches->get(p);
    cout_doing << "  patch " << coarsePatch->getID() << std::endl;

    CCVariable<double> refineCell;
    new_dw->getModifiable(refineCell,
                          d_mpmLabels->MPMRefineCellLabel,
                          0,
                          coarsePatch);
    bool computesAve = true;

    fineToCoarseOperator<double>(refineCell,
                                 computesAve,
                                 d_mpmLabels->MPMRefineCellLabel,
                                 0,
                                 new_dw,
                                 coarsePatch,
                                 coarseLevel,
                                 fineLevel);

    if (coarseLevel->getIndex() == numLevels - 2) {
      //    std::cout << "coarseLevelIndex = " << coarseLevel->getIndex() <<
      //    std::endl;
      max_vartype xmax, ymax, zmax;
      min_vartype xmin, ymin, zmin;
      new_dw->get(xmax, RefineFlagXMaxLabel);
      new_dw->get(ymax, RefineFlagYMaxLabel);
      new_dw->get(zmax, RefineFlagZMaxLabel);
      new_dw->get(xmin, RefineFlagXMinLabel);
      new_dw->get(ymin, RefineFlagYMinLabel);
      new_dw->get(zmin, RefineFlagZMinLabel);

      //    std::cout << "xmax = " << xmax << std::endl;
      //    std::cout << "ymax = " << ymax << std::endl;
      //    std::cout << "zmax = " << zmax << std::endl;
      //    std::cout << "xmin = " << xmin << std::endl;
      //    std::cout << "ymin = " << ymin << std::endl;
      //    std::cout << "zmin = " << zmin << std::endl;

      IntVector fineXYZMaxMin(xmax, ymax, zmax);
      IntVector fineXYZMinMax(xmin, ymin, zmin);
      IntVector fineXZMaxYMin(xmax, ymin, zmax);
      IntVector fineXZMinYMax(xmin, ymax, zmin);
      IntVector fineXYMaxZMin(xmax, ymax, zmin);
      IntVector fineXYMinZMax(xmin, ymin, zmax);
      IntVector fineXMinYZMax(xmin, ymax, zmax);
      IntVector fineXMaxYZMin(xmax, ymin, zmin);

      IntVector coarseMinMax[8];

      coarseMinMax[0] = fineXYZMaxMin / RR;
      coarseMinMax[1] = fineXYZMinMax / RR;

      // Set the refine flags to 1 in all cells in the interior of the minimum
      // and maximum to ensure a rectangular region is refined.
      int imax, jmax, kmax, imin, jmin, kmin;
      imax = coarseMinMax[0].x();
      jmax = coarseMinMax[0].y();
      kmax = coarseMinMax[0].z();
      imin = coarseMinMax[1].x();
      jmin = coarseMinMax[1].y();
      kmin = coarseMinMax[1].z();
      for (CellIterator iter = coarsePatch->getCellIterator(); !iter.done();
           iter++) {
        IntVector c = *iter;
        if (c.x() >= imin && c.x() <= imax && c.y() >= jmin && c.y() <= jmax &&
            c.z() >= kmin && c.z() <= kmax) {
          refineCell[c] = 1;
        }
      }

    } // end if level
  }   // end patch
}

//______________________________________________________________________
// Schedule to mark flags for AMR regridding
void
AMRMPM::scheduleErrorEstimate(const LevelP& coarseLevel, SchedulerP& sched)
{
  //  std::cout << "scheduleErrorEstimate" << std::endl;
  printSchedule(coarseLevel, cout_doing, "AMRMPM::scheduleErrorEstimate");

  Task* task =
    scinew Task("AMRMPM::errorEstimate", this, &AMRMPM::errorEstimate);

  task->modifies(d_regridder->getRefineFlagLabel(),
                 d_regridder->refineFlagMaterials());
  task->modifies(d_regridder->getRefinePatchFlagLabel(),
                 d_regridder->refineFlagMaterials());
  task->requires(Task::NewDW, d_mpmLabels->MPMRefineCellLabel, Ghost::None);

  sched->addTask(task,
                 coarseLevel->eachPatch(),
                 d_materialManager->allMaterials("MPM"));
}
//______________________________________________________________________
// Schedule to mark initial flags for AMR regridding
void
AMRMPM::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                     SchedulerP& sched)
{
  //  std::cout << "scheduleInitialErrorEstimate" << std::endl;
  //  std::cout << "Doing nothing for now" << std::endl;

  //  scheduleErrorEstimate(coarseLevel, sched);
}

//______________________________________________________________________
//
void
AMRMPM::printParticleCount(const ProcessorGroup* pg,
                           const PatchSubset*,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  sumlong_vartype pcount;
  new_dw->get(pcount, d_mpmLabels->partCountLabel);

  if (pg->myRank() == 0) {
    std::cout << "Created " << (long)pcount << " total particles" << std::endl;
  }
}
//______________________________________________________________________
//
void
AMRMPM::actuallyInitialize(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  const Level* level           = getLevel(patches);
  int levelIndex               = level->getIndex();
  particleIndex totalParticles = 0;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing AMRMPM::actuallyInitialize");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    for (int m = 0; m < matls->size(); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int indx = mpm_matl->getDWIndex();

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

      // Initialize deformation gradient
      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

      mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                         mpm_matl,
                                                         new_dw);

      //__________________________________
      // color particles according to what level they're on
      if (d_mpmFlags->d_withColor) {
        ParticleSubset* pset = new_dw->getParticleSubset(indx, patch);
        ParticleVariable<double> pColor;
        new_dw->getModifiable(pColor, d_mpmLabels->pColorLabel, pset);

        ParticleSubset::iterator iter = pset->begin();
        for (; iter != pset->end(); iter++) {
          particleIndex idx = *iter;
          pColor[idx]       = levelIndex;
        }
      }
    } // matl loop
  }

  if (d_mpmFlags->d_reductionVars->accStrainEnergy) {
    // Initialize the accumulated strain energy
    new_dw->put(max_vartype(0.0), d_mpmLabels->AccStrainEnergyLabel);
  }

  new_dw->put(sumlong_vartype(totalParticles), d_mpmLabels->partCountLabel);
}
//______________________________________________________________________
//
void
AMRMPM::actuallyComputeStableTimestep(const ProcessorGroup*,
                                      const PatchSubset*,
                                      const MaterialSubset*,
                                      DataWarehouse*,
                                      DataWarehouse*)
{
}
//______________________________________________________________________
//  This task computes the partition of unity for each particle
//
void
AMRMPM::partitionOfUnity(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset*,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw)
{
  const Level* curLevel = getLevel(patches);
  int curLevelIndex     = curLevel->getIndex();
  Vector dX             = curLevel->dCell();
  Vector dX_fine        = 0.5 * dX;
  Vector dX_coarse      = 2.0 * dX;
  Vector RRC            = dX / dX_coarse;
  Vector RRF            = dX / dX_fine;
  if (curLevel->hasFinerLevel()) {
    dX_fine = curLevel->getFinerLevel()->dCell();
    RRF     = dX / dX_fine;
    RRC     = Vector(1. / RRF.x(), 1. / RRF.y(), 1. / RRF.z());
  }
  if (curLevel->hasCoarserLevel()) {
    dX_coarse = curLevel->getCoarserLevel()->dCell();
    RRC       = dX / dX_coarse;
    RRF       = Vector(1. / RRC.x(), 1. / RRC.y(), 1. / RRC.z());
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing AMRMPM::partitionOfUnity");

    // Create and Initialize refine flags to be modified later
    CCVariable<double> refineCell;
    new_dw->allocateAndPut(refineCell,
                           d_mpmLabels->MPMRefineCellLabel,
                           0,
                           patch);
    refineCell.initialize(0.0);

    int numMatls      = d_materialManager->getNumMaterials("MPM");
    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    const Matrix3 notUsed;

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
      constParticleVariable<Point> pX;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<int> plastlevel;
      ParticleVariable<Matrix3> pSizenew;
      ParticleVariable<int> plastlevelnew;
      ParticleVariable<double> partitionUnity;

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
      old_dw->get(plastlevel, d_mpmLabels->pLastLevelLabel, pset);
      new_dw->allocateAndPut(pSizenew, d_mpmLabels->pSizeLabel_preReloc, pset);
      new_dw->allocateAndPut(plastlevelnew,
                             d_mpmLabels->pLastLevelLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(partitionUnity,
                             d_mpmLabels->pPartitionUnityLabel,
                             pset);

      int n8or27 = d_mpmFlags->d_8or27;

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        Matrix3 ps = pSize[idx];

        if (curLevelIndex < plastlevel[idx]) {
          pSizenew[idx] = Matrix3(ps(0, 0) * RRC.x(),
                                  ps(0, 1) * RRC.x(),
                                  ps(0, 2) * RRC.x(),
                                  ps(1, 0) * RRC.y(),
                                  ps(1, 1) * RRC.y(),
                                  ps(1, 2) * RRC.y(),
                                  ps(2, 0) * RRC.z(),
                                  ps(2, 1) * RRC.z(),
                                  ps(2, 2) * RRC.z());
        } else if (curLevelIndex > plastlevel[idx]) {
          pSizenew[idx] = Matrix3(ps(0, 0) * RRF.x(),
                                  ps(0, 1) * RRF.x(),
                                  ps(0, 2) * RRF.x(),
                                  ps(1, 0) * RRF.y(),
                                  ps(1, 1) * RRF.y(),
                                  ps(1, 2) * RRF.y(),
                                  ps(2, 0) * RRF.z(),
                                  ps(2, 1) * RRF.z(),
                                  ps(2, 2) * RRF.z());
        } else {
          pSizenew[idx] = pSize[idx];
        }

        plastlevelnew[idx] = curLevelIndex;

        partitionUnity[idx] = 0;

        interpolator->findCellAndWeights(pX[idx], ni, S, pSize[idx], notUsed);

        for (int k = 0; k < n8or27; k++) {
          partitionUnity[idx] += S[k];
        }
      }
    } // loop over materials
  }   // loop over patches
}
//______________________________________________________________________
//
void
AMRMPM::interpolateParticlesToGrid(const ProcessorGroup*,
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
              "Doing AMRMPM::interpolateParticlesToGrid");

    int numMatls      = d_materialManager->getNumMaterials("MPM");
    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());

#ifdef CBDI_FLUXBCS
    LinearInterpolator* LPI;
    LPI = scinew LinearInterpolator(patch);
    std::vector<IntVector> ni_LPI(LPI->size());
    std::vector<double> S_LPI(LPI->size());
#endif

    Ghost::GhostType gan = Ghost::AroundNodes;
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Create arrays for the particle data
      constParticleVariable<Point> pX;
      constParticleVariable<double> pMass, pVolume, pTemperature;
      constParticleVariable<double> pConcentration;
      constParticleVariable<double> pExternalScalarFlux;
      constParticleVariable<Vector> pVelocity, pExternalforce;
      constParticleVariable<Matrix3> pSize;
      constParticleVariable<Matrix3> pDefGrad;
      constParticleVariable<Matrix3> pStress;
      constParticleVariable<Matrix3> pVelGrad;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pVolume, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);

#ifdef CBDI_FLUXBCS
      constParticleVariable<int> pLoadCurveID;
      if (d_mpmFlags->d_useLoadCurves) {
        old_dw->get(pLoadCurveID, d_mpmLabels->pLoadCurveIDLabel, pset);
      }
#endif

      new_dw->get(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);
      old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);
      new_dw->get(pExternalforce, d_mpmLabels->pExtForceLabel_preReloc, pset);
      if (d_mpmFlags->d_GEVelProj) {
        old_dw->get(pVelGrad, d_mpmLabels->pVelGradLabel, pset);
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->get(pExternalScalarFlux,
                    d_mpmLabels->pExternalScalarFluxLabel,
                    pset);
        old_dw->get(pConcentration,
                    d_mpmLabels->diffusion->pConcentration,
                    pset);
        old_dw->get(pStress, d_mpmLabels->pStressLabel, pset);
      }

      // Create arrays for the grid data
      NCVariable<double> gMass;
      NCVariable<double> gVolume;
      NCVariable<Vector> gVelocity;
      NCVariable<Vector> ggExtForce;
      NCVariable<double> gTemperature;
      NCVariable<double> gTemperatureRate;
      NCVariable<double> gconcentration;
      NCVariable<double> gextscalarflux;
      NCVariable<double> ghydrostaticstress;

      new_dw->allocateAndPut(gMass, d_mpmLabels->gMassLabel, dwi, patch);
      new_dw->allocateAndPut(gVolume, d_mpmLabels->gVolumeLabel, dwi, patch);
      new_dw->allocateAndPut(gVelocity,
                             d_mpmLabels->gVelocityLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperature,
                             d_mpmLabels->gTemperatureLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gTemperatureRate,
                             d_mpmLabels->gTemperatureRateLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(ggExtForce,
                             d_mpmLabels->gExternalForceLabel,
                             dwi,
                             patch);

      gMass.initialize(d_SMALL_NUM_MPM);
      gVolume.initialize(d_SMALL_NUM_MPM);
      gVelocity.initialize(Vector(0, 0, 0));
      ggExtForce.initialize(Vector(0, 0, 0));
      gTemperature.initialize(0);
      gTemperatureRate.initialize(0);

      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->allocateAndPut(gconcentration,
                               d_mpmLabels->gConcentrationLabel,
                               dwi,
                               patch);
        new_dw->allocateAndPut(ghydrostaticstress,
                               d_mpmLabels->diffusion->gHydrostaticStress,
                               dwi,
                               patch);
        new_dw->allocateAndPut(gextscalarflux,
                               d_mpmLabels->gExternalScalarFluxLabel,
                               dwi,
                               patch);
        gconcentration.initialize(0);
        ghydrostaticstress.initialize(0);
        gextscalarflux.initialize(0);
      }

      Vector pmom;
      int n8or27 = d_mpmFlags->d_8or27;

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pDefGrad[idx]);

        pmom = pVelocity[idx] * pMass[idx];

        // Add each particles contribution to the local mass & velocity
        IntVector node;
        for (int k = 0; k < n8or27; k++) {
          node = ni[k];
          if (patch->containsNode(node)) {
            if (d_mpmFlags->d_GEVelProj) {
              Point gpos      = patch->getNodePosition(node);
              Vector distance = pX[idx] - gpos;
              Vector pVel_ext = pVelocity[idx] - pVelGrad[idx] * distance;
              pmom            = pVel_ext * pMass[idx];
            }
            gMass[node] += pMass[idx] * S[k];
            gVelocity[node] += pmom * S[k];
            gVolume[node] += pVolume[idx] * S[k];
            ggExtForce[node] += pExternalforce[idx] * S[k];
            gTemperature[node] += pTemperature[idx] * pMass[idx] * S[k];
          }
        }
        if (d_mpmFlags->d_doScalarDiffusion) {
          double one_third    = 1. / 3.;
          double phydrostress = one_third * pStress[idx].Trace();
          for (int k = 0; k < n8or27; k++) {
            node = ni[k];
            if (patch->containsNode(node)) {
              ghydrostaticstress[node] += phydrostress * pMass[idx] * S[k];
              gconcentration[node] += pConcentration[idx] * pMass[idx] * S[k];
#ifndef CBDI_FLUXBCS
              gextscalarflux[node] += pExternalScalarFlux[idx] * S[k];
#endif
            }
          }
        }
      } // End of particle loop

#ifdef CBDI_FLUXBCS
      Vector dx = patch->dCell();
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        Point flux_pos;
        if (pLoadCurveID[idx] == 1) {
          flux_pos = Point(pX[idx].x() - 0.5 * pSize[idx](0, 0) * dx.x(),
                           pX[idx].y(),
                           pX[idx].z());
        }
        if (pLoadCurveID[idx] == 2) {
          flux_pos = Point(pX[idx].x() + 0.5 * pSize[idx](0, 0) * dx.x(),
                           pX[idx].y(),
                           pX[idx].z());
        }
        if (pLoadCurveID[idx] == 3) {
          flux_pos = Point(pX[idx].x(),
                           pX[idx].y() + 0.5 * pSize[idx](1, 1) * dx.y(),
                           pX[idx].z());
        }
        LPI->findCellAndWeights(flux_pos,
                                ni_LPI,
                                S_LPI,
                                pSize[idx],
                                pDefGrad[idx]);
        for (int k = 0; k < (int)ni_LPI.size(); k++) {
          if (patch->containsNode(ni_LPI[k])) {
            gextscalarflux[ni_LPI[k]] += pExternalScalarFlux[idx] * S_LPI[k];
          }
        }
      }
#endif

      // gVelocity and gtemperature are divided by gMass in
      // AMRMPM::NormalizeNodalVelTempConc() task

    } // End loop over materials
  }   // End loop over patches
}

//______________________________________________________________________
//  At the CFI fine patch nodes add contributions from the coarse level
//  particles.
void
AMRMPM::interpolateParticlesToGrid_CFI(const ProcessorGroup*,
                                       const PatchSubset* finePatches,
                                       const MaterialSubset*,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  const Level* fineLevel   = getLevel(finePatches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
  IntVector refineRatio(fineLevel->getRefinementRatio());

  for (int fp = 0; fp < finePatches->size(); fp++) {
    const Patch* finePatch = finePatches->get(fp);
    printTask(finePatches,
              finePatch,
              cout_doing,
              "Doing AMRMPM::interpolateParticlesToGrid_CFI");

    int numMatls      = d_materialManager->getNumMaterials("MPM");
    auto interpolator = d_mpmFlags->d_interpolator->clone(finePatch);

    constNCVariable<Stencil7> zoi_fine;
    new_dw->get(zoi_fine, d_mpmLabels->gZOILabel, 0, finePatch, Ghost::None, 0);

    // Determine extents for coarser level particle data
    // Linear Interpolation:  1 layer of coarse level cells
    // Gimp Interpolation:    2 layers
    /*`==========TESTING==========*/
    IntVector nLayers(d_nPaddingCells_Coarse,
                      d_nPaddingCells_Coarse,
                      d_nPaddingCells_Coarse);
    IntVector nPaddingCells = nLayers * (fineLevel->getRefinementRatio());
    // std::cout << " nPaddingCells " << nPaddingCells << "nLayers " << nLayers
    // << endl;
    /*===========TESTING==========`*/

    int nGhostCells           = 0;
    bool returnExclusiveRange = false;
    IntVector cl_tmp, ch_tmp, fl, fh;

    getCoarseLevelRange(finePatch,
                        coarseLevel,
                        cl_tmp,
                        ch_tmp,
                        fl,
                        fh,
                        nPaddingCells,
                        nGhostCells,
                        returnExclusiveRange);

    //  expand cl_tmp when a neighor patch exists.
    //  This patch owns the low nodes.  You need particles
    //  from the neighbor patch.
    cl_tmp -= finePatch->neighborsLow() * nLayers;

    // find the coarse patches under the fine patch.
    // You must add a single layer of padding cells.
    int padding = 1;
    Level::selectType coarsePatches;
    finePatch->getOtherLevelPatches(-1, coarsePatches, padding);

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // get fine level nodal data
      NCVariable<double> gMass_fine;
      NCVariable<double> gVolume_fine;
      NCVariable<Vector> gVelocity_fine;
      NCVariable<Vector> gExternalforce_fine;
      NCVariable<double> gTemperature_fine;
      NCVariable<double> gConc_fine;
      NCVariable<double> gExtScalarFlux_fine;
      NCVariable<double> gHStress_fine;

      new_dw->getModifiable(gMass_fine,
                            d_mpmLabels->gMassLabel,
                            dwi,
                            finePatch);
      new_dw->getModifiable(gVolume_fine,
                            d_mpmLabels->gVolumeLabel,
                            dwi,
                            finePatch);
      new_dw->getModifiable(gVelocity_fine,
                            d_mpmLabels->gVelocityLabel,
                            dwi,
                            finePatch);
      new_dw->getModifiable(gTemperature_fine,
                            d_mpmLabels->gTemperatureLabel,
                            dwi,
                            finePatch);
      new_dw->getModifiable(gExternalforce_fine,
                            d_mpmLabels->gExternalForceLabel,
                            dwi,
                            finePatch);
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->getModifiable(gConc_fine,
                              d_mpmLabels->gConcentrationLabel,
                              dwi,
                              finePatch);
        new_dw->getModifiable(gExtScalarFlux_fine,
                              d_mpmLabels->gConcentrationLabel,
                              dwi,
                              finePatch);
        new_dw->getModifiable(gHStress_fine,
                              d_mpmLabels->diffusion->gHydrostaticStress,
                              dwi,
                              finePatch);
      }

      // loop over the coarse patches under the fine patches.
      for (size_t cp = 0; cp < coarsePatches.size(); cp++) {
        const Patch* coarsePatch = coarsePatches[cp];

        // get coarse level particle data
        constParticleVariable<Point> pX_coarse;
        constParticleVariable<double> pMass_coarse;
        constParticleVariable<double> pVolume_coarse;
        constParticleVariable<double> pTemperature_coarse;
        constParticleVariable<Vector> pVelocity_coarse;
        constParticleVariable<Vector> pExternalforce_coarse;
        constParticleVariable<double> pConc_coarse;
        constParticleVariable<double> pExtScalarFlux_c;
        constParticleVariable<Matrix3> pStress_coarse;

        // coarseLow and coarseHigh cannot lie outside of the coarse patch
        IntVector cl = Max(cl_tmp, coarsePatch->getCellLowIndex());
        IntVector ch = Min(ch_tmp, coarsePatch->getCellHighIndex());

        ParticleSubset* pset = 0;

        pset = old_dw->getParticleSubset(dwi,
                                         cl,
                                         ch,
                                         coarsePatch,
                                         d_mpmLabels->pXLabel);
#if 0
        std::cout << "  coarseLevel: " << coarseLevel->getIndex() << std::endl;
        std::cout << " cl_tmp: "<< cl_tmp << " ch_tmp: " << ch_tmp << std::endl;
        std::cout << " cl:     " << cl    << " ch:     " << ch<< " fl: " << fl << " fh " << fh << std::endl;
        std::cout << "  " << *pset << std::endl;
#endif
        old_dw->get(pX_coarse, d_mpmLabels->pXLabel, pset);
        old_dw->get(pMass_coarse, d_mpmLabels->pMassLabel, pset);
        old_dw->get(pVolume_coarse, d_mpmLabels->pVolumeLabel, pset);
        old_dw->get(pVelocity_coarse, d_mpmLabels->pVelocityLabel, pset);
        old_dw->get(pTemperature_coarse, d_mpmLabels->pTemperatureLabel, pset);
        new_dw->get(pExternalforce_coarse,
                    d_mpmLabels->pExtForceLabel_preReloc,
                    pset);
        if (d_mpmFlags->d_doScalarDiffusion) {
          old_dw->get(pConc_coarse,
                      d_mpmLabels->diffusion->pConcentration,
                      pset);
          new_dw->get(pExtScalarFlux_c,
                      d_mpmLabels->pExternalScalarFluxLabel,
                      pset);
          old_dw->get(pStress_coarse, d_mpmLabels->pStressLabel, pset);
        }

        for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
             iter++) {
          particleIndex idx = *iter;

          // Get the node indices that surround the fine patch cell
          std::vector<IntVector> ni;
          std::vector<double> S;

          interpolator->findCellAndWeights_CFI(pX_coarse[idx], ni, S, zoi_fine);

          Vector pmom = pVelocity_coarse[idx] * pMass_coarse[idx];

          // Add each particle's contribution to the local mass & velocity
          IntVector fineNode;

          for (int k = 0; k < (int)ni.size(); k++) {
            fineNode = ni[k];

            gMass_fine[fineNode] += pMass_coarse[idx] * S[k];
            gVelocity_fine[fineNode] += pmom * S[k];
            gVolume_fine[fineNode] += pVolume_coarse[idx] * S[k];
            gExternalforce_fine[fineNode] += pExternalforce_coarse[idx] * S[k];
            gTemperature_fine[fineNode] +=
              pTemperature_coarse[idx] * pMass_coarse[idx] * S[k];
          }
          if (d_mpmFlags->d_doScalarDiffusion) {
            double one_third = 1. / 3.;
            double ConcMass  = pConc_coarse[idx] * pMass_coarse[idx];
            double phydrostressmass =
              one_third * pStress_coarse[idx].Trace() * pMass_coarse[idx];
            double pESFlux_c = pExtScalarFlux_c[idx];

            for (int k = 0; k < (int)ni.size(); k++) {
              fineNode = ni[k];
              gConc_fine[fineNode] += ConcMass * S[k];
              gExtScalarFlux_fine[fineNode] += pESFlux_c * S[k];
              gHStress_fine[fineNode] += phydrostressmass * S[k];
            }
          }
        } // End of particle loop
      }   // loop over coarse patches
    }     // End loop over materials
  }       // End loop over fine patches
}

//______________________________________________________________________
//  copy the fine level nodal data to the underlying coarse nodes at the CFI.
void
AMRMPM::coarsenNodalData_CFI(const ProcessorGroup*,
                             const PatchSubset* coarsePatches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw,
                             const coarsenFlag flag)
{
  Level::selectType coarseCFI_Patches;

  coarseLevelCFI_Patches(coarsePatches, coarseCFI_Patches);

  //__________________________________
  // From the coarse patch look up to the fine patches that have
  // coarse fine interfaces.
  const Level* coarseLevel = getLevel(coarsePatches);

  for (size_t p = 0; p < coarseCFI_Patches.size(); p++) {
    const Patch* coarsePatch = coarseCFI_Patches[p];

    string txt = "(zero)";
    if (flag == coarsenData) {
      txt = "(coarsen)";
    }
    printTask(coarsePatch,
              cout_doing,
              "Doing AMRMPM::coarsenNodalData_CFI" + txt);

    int numMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // get coarse level data
      NCVariable<double> gMass_coarse;
      NCVariable<double> gVolume_coarse;
      NCVariable<Vector> gVelocity_coarse;
      NCVariable<Vector> gVelocityStar_coarse;
      NCVariable<Vector> gAcceleration_coarse;
      NCVariable<Vector> gExternalforce_coarse;
      NCVariable<double> gTemperature_coarse;
      NCVariable<double> gConcentration_coarse;
      NCVariable<double> gExtScalarFlux_coarse;

      new_dw->getModifiable(gMass_coarse,
                            d_mpmLabels->gMassLabel,
                            dwi,
                            coarsePatch);
      new_dw->getModifiable(gVolume_coarse,
                            d_mpmLabels->gVolumeLabel,
                            dwi,
                            coarsePatch);
      new_dw->getModifiable(gVelocity_coarse,
                            d_mpmLabels->gVelocityLabel,
                            dwi,
                            coarsePatch);
      new_dw->getModifiable(gTemperature_coarse,
                            d_mpmLabels->gTemperatureLabel,
                            dwi,
                            coarsePatch);
      new_dw->getModifiable(gExternalforce_coarse,
                            d_mpmLabels->gExternalForceLabel,
                            dwi,
                            coarsePatch);
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->getModifiable(gConcentration_coarse,
                              d_mpmLabels->gConcentrationLabel,
                              dwi,
                              coarsePatch);
        new_dw->getModifiable(gExtScalarFlux_coarse,
                              d_mpmLabels->gExternalScalarFluxLabel,
                              dwi,
                              coarsePatch);
      }

      if (flag == zeroData) {
        new_dw->getModifiable(gVelocityStar_coarse,
                              d_mpmLabels->gVelocityStarLabel,
                              dwi,
                              coarsePatch);
        new_dw->getModifiable(gAcceleration_coarse,
                              d_mpmLabels->gAccelerationLabel,
                              dwi,
                              coarsePatch);
      }

      //__________________________________
      // Iterate over coarse/fine interface faces
      ASSERT(coarseLevel->hasFinerLevel());
      const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);

      // loop over all the fine level patches
      for (size_t fp = 0; fp < finePatches.size(); fp++) {
        const Patch* finePatch = finePatches[fp];
        if (finePatch->hasCoarseFaces()) {

          // get fine level data
          constNCVariable<double> gMass_fine;
          constNCVariable<double> gVolume_fine;
          constNCVariable<Vector> gVelocity_fine;
          constNCVariable<double> gTemperature_fine;
          constNCVariable<Vector> gExternalforce_fine;
          constNCVariable<double> gConcentration_fine;
          constNCVariable<double> gExtScalarFlux_fine;

          if (flag == coarsenData) {
            // use getRegion() instead of get().  They should be equivalent but
            // get() throws assert on parallel runs.
            IntVector fl = finePatch->getNodeLowIndex();
            IntVector fh = finePatch->getNodeHighIndex();
            new_dw->getRegion(gMass_fine,
                              d_mpmLabels->gMassLabel,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
            new_dw->getRegion(gVolume_fine,
                              d_mpmLabels->gVolumeLabel,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
            new_dw->getRegion(gVelocity_fine,
                              d_mpmLabels->gVelocityLabel,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
            new_dw->getRegion(gTemperature_fine,
                              d_mpmLabels->gTemperatureLabel,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
            new_dw->getRegion(gExternalforce_fine,
                              d_mpmLabels->gExternalForceLabel,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
            if (d_mpmFlags->d_doScalarDiffusion) {
              new_dw->getRegion(gConcentration_fine,
                                d_mpmLabels->gConcentrationLabel,
                                dwi,
                                fineLevel,
                                fl,
                                fh);
              new_dw->getRegion(gExtScalarFlux_fine,
                                d_mpmLabels->gExternalScalarFluxLabel,
                                dwi,
                                fineLevel,
                                fl,
                                fh);
            }
          }

          std::vector<Patch::FaceType> cf;
          finePatch->getCoarseFaces(cf);

          // Iterate over coarse/fine interface faces
          std::vector<Patch::FaceType>::const_iterator iter;
          for (iter = cf.begin(); iter != cf.end(); ++iter) {
            Patch::FaceType patchFace = *iter;

            // determine the iterator on the coarse level.
            NodeIterator n_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
            bool isRight_CP_FP_pair;

            coarseLevel_CFI_NodeIterator(patchFace,
                                         coarsePatch,
                                         finePatch,
                                         fineLevel,
                                         n_iter,
                                         isRight_CP_FP_pair);

            // Is this the right coarse/fine patch pair
            if (isRight_CP_FP_pair) {
              switch (flag) {
                case coarsenData:
                  for (; !n_iter.done(); n_iter++) {
                    IntVector c_node = *n_iter;
                    IntVector f_node = coarseLevel->mapNodeToFiner(c_node);

                    // only overwrite coarse data if there is non-zero fine data
                    if (gMass_fine[f_node] > 2 * d_SMALL_NUM_MPM) {
                      gMass_coarse[c_node]        = gMass_fine[f_node];
                      gVolume_coarse[c_node]      = gVolume_fine[f_node];
                      gVelocity_coarse[c_node]    = gVelocity_fine[f_node];
                      gTemperature_coarse[c_node] = gTemperature_fine[f_node];
                      gExternalforce_coarse[c_node] =
                        gExternalforce_fine[f_node];
                      if (d_mpmFlags->d_doScalarDiffusion) {
                        gConcentration_coarse[c_node] =
                          gConcentration_fine[f_node];
                        gExtScalarFlux_coarse[c_node] =
                          gExtScalarFlux_fine[f_node];
                      }
                    } // if mass
                  }   // end node iterator loop
                  break;
                case zeroData:
                  for (; !n_iter.done(); n_iter++) {
                    IntVector c_node     = *n_iter;
                    IntVector f_node     = coarseLevel->mapNodeToFiner(c_node);
                    gMass_coarse[c_node] = 0;
                    gVolume_coarse[c_node]        = 0;
                    gVelocity_coarse[c_node]      = Vector(0, 0, 0);
                    gVelocityStar_coarse[c_node]  = Vector(0, 0, 0);
                    gAcceleration_coarse[c_node]  = Vector(0, 0, 0);
                    gTemperature_coarse[c_node]   = 0;
                    gExternalforce_coarse[c_node] = Vector(0, 0, 0);
                    if (d_mpmFlags->d_doScalarDiffusion) {
                      gConcentration_coarse[c_node] = 0;
                      gExtScalarFlux_coarse[c_node] = 0;
                    }
                  } // end node iterator loop
                  break;
              }
            } //  isRight_CP_FP_pair
          }   //  end CFI face loop
        }     //  if finepatch has coarse face
      }       //  end fine Patch loop
    }         //  end matl loop
  }           // end coarse patch loop
}

//______________________________________________________________________
//  copy the fine level nodal data to the underlying coarse nodes at the CFI.
void
AMRMPM::coarsenNodalData_CFI2(const ProcessorGroup*,
                              const PatchSubset* coarsePatches,
                              const MaterialSubset*,
                              DataWarehouse* old_dw,
                              DataWarehouse* new_dw)
{
  Level::selectType coarseCFI_Patches;
  coarseLevelCFI_Patches(coarsePatches, coarseCFI_Patches);

  //__________________________________
  // From the coarse patch look up to the fine patches that have
  // coarse fine interfaces.
  const Level* coarseLevel = getLevel(coarsePatches);

  for (size_t p = 0; p < coarseCFI_Patches.size(); p++) {
    const Patch* coarsePatch = coarseCFI_Patches[p];

    printTask(coarsePatch, cout_doing, "Doing AMRMPM::coarsenNodalData_CFI2");

    int numMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // get coarse level data
      NCVariable<Vector> internalForce_coarse;
      new_dw->getModifiable(internalForce_coarse,
                            d_mpmLabels->gInternalForceLabel,
                            dwi,
                            coarsePatch);
      NCVariable<double> gConcRate_coarse;
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->getModifiable(gConcRate_coarse,
                              d_mpmLabels->diffusion->gConcentrationRate,
                              dwi,
                              coarsePatch);
      }

      //__________________________________
      // Iterate over coarse/fine interface faces
      ASSERT(coarseLevel->hasFinerLevel());
      const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();

      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);

      // loop over all the coarse level patches
      for (int fp = 0; fp < finePatches.size(); fp++) {
        const Patch* finePatch = finePatches[fp];
        if (finePatch->hasCoarseFaces()) {

          // get fine level data
          constNCVariable<double> gMass_fine, gConcRate_fine;
          constNCVariable<Vector> internalForce_fine;

          // use getRegion() instead of get().  They should be equivalent but
          // get() throws assert on parallel runs.
          IntVector fl = finePatch->getNodeLowIndex();
          IntVector fh = finePatch->getNodeHighIndex();
          new_dw->getRegion(gMass_fine,
                            d_mpmLabels->gMassLabel,
                            dwi,
                            fineLevel,
                            fl,
                            fh);
          new_dw->getRegion(internalForce_fine,
                            d_mpmLabels->gInternalForceLabel,
                            dwi,
                            fineLevel,
                            fl,
                            fh);

          if (d_mpmFlags->d_doScalarDiffusion) {
            new_dw->getRegion(gConcRate_fine,
                              d_mpmLabels->diffusion->gConcentrationRate,
                              dwi,
                              fineLevel,
                              fl,
                              fh);
          }
          std::vector<Patch::FaceType> cf;
          finePatch->getCoarseFaces(cf);

          // Iterate over coarse/fine interface faces
          std::vector<Patch::FaceType>::const_iterator iter;
          for (iter = cf.begin(); iter != cf.end(); ++iter) {
            Patch::FaceType patchFace = *iter;

            // determine the iterator on the coarse level.
            NodeIterator n_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
            bool isRight_CP_FP_pair;

            coarseLevel_CFI_NodeIterator(patchFace,
                                         coarsePatch,
                                         finePatch,
                                         fineLevel,
                                         n_iter,
                                         isRight_CP_FP_pair);

            // Is this the right coarse/fine patch pair
            if (isRight_CP_FP_pair) {

              for (; !n_iter.done(); n_iter++) {
                IntVector c_node = *n_iter;

                IntVector f_node = coarseLevel->mapNodeToFiner(c_node);

                // only overwrite coarse data if there is non-zero fine data
                if (gMass_fine[f_node] > 2 * d_SMALL_NUM_MPM) {

                  internalForce_coarse[c_node] = internalForce_fine[f_node];

                  if (d_mpmFlags->d_doScalarDiffusion) {
                    gConcRate_coarse[c_node] = gConcRate_fine[f_node];
                  }

/*`==========TESTING==========*/
#if 0
                  if( internalForce_coarse[c_node].length()  >1e-8){
                     std::ostringstream warn;
                    warn << "Too Big: " << c_node << " f_node " << f_node 
                         << "    L-"<< fineLevel->getIndex()
                         <<" InternalForce_fine   " << internalForce_fine[f_node] 
                         <<" InternalForce_coarse " << internalForce_coarse[c_node] << std::endl;
                     
                    throw InternalError(warn.str(), __FILE__, __LINE__);
                  }
#endif
                  /*===========TESTING==========`*/
                }
              } //  node loop
            }   //  isRight_CP_FP_pair
          }     //  end CFI face loop
        }       //  if finepatch has coarse face
      }         //  end fine Patch loop
    }           //  end matl loop
  }             //  end coarse patch loop
}

//______________________________________________________________________
// Divide gVelocity and gTemperature by gMass
void
AMRMPM::normalizeNodalVelTempConc(const ProcessorGroup*,
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
              "Doing AMRMPM::normalizeNodalVelTempConc");

    int numMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // get  level nodal data
      constNCVariable<double> gMass;
      NCVariable<Vector> gVelocity;
      NCVariable<double> gTemperature;
      NCVariable<double> gConcentration;
      NCVariable<double> gConcentrationNoBC;
      NCVariable<double> gHydroStress;
      Ghost::GhostType gn = Ghost::None;

      new_dw->get(gMass, d_mpmLabels->gMassLabel, dwi, patch, gn, 0);
      new_dw->getModifiable(gVelocity,
                            d_mpmLabels->gVelocityLabel,
                            dwi,
                            patch,
                            gn,
                            0);
      new_dw->getModifiable(gTemperature,
                            d_mpmLabels->gTemperatureLabel,
                            dwi,
                            patch,
                            gn,
                            0);
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->getModifiable(gConcentration,
                              d_mpmLabels->gConcentrationLabel,
                              dwi,
                              patch,
                              gn,
                              0);
        new_dw->getModifiable(gHydroStress,
                              d_mpmLabels->diffusion->gHydrostaticStress,
                              dwi,
                              patch,
                              gn,
                              0);
        new_dw->allocateAndPut(gConcentrationNoBC,
                               d_mpmLabels->gConcentrationNoBCLabel,
                               dwi,
                               patch);
      }

      //__________________________________
      //  back out the nodal quantities
      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector n = *iter;
        gVelocity[n] /= gMass[n];
        gTemperature[n] /= gMass[n];
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector n = *iter;
          gConcentration[n] /= gMass[n];
          gHydroStress[n] /= gMass[n];
          gConcentrationNoBC[n] = gConcentration[n];
        }
      }

      // Apply boundary conditions to the temperature and velocity (if symmetry)
      MPMBoundCond bc;
      string interp_type = d_mpmFlags->d_interpolatorType;
      bc.setBoundaryCondition(patch,
                              dwi,
                              "Temperature",
                              gTemperature,
                              interp_type);
      bc.setBoundaryCondition(patch, dwi, "Symmetric", gVelocity, interp_type);
      if (d_mpmFlags->d_doScalarDiffusion) {
        bc.setBoundaryCondition(patch,
                                dwi,
                                "SD-Type",
                                gConcentration,
                                interp_type);
      }
    } // End loop over materials
  }   // End loop over fine patches
}

//______________________________________________________________________
//
void
AMRMPM::computeStressTensor(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  printTask(patches,
            patches->get(0),
            cout_doing,
            "Doing AMRMPM::computeStressTensor");

  for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {

    MPMMaterial* mpm_matl =
      static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();

    cm->setWorld(UintahParallelComponent::d_myworld);
    cm->computeStressTensor(patches, mpm_matl, old_dw, new_dw);
  }
}

//______________________________________________________________________
//
void
AMRMPM::updateErosionParameter(const ProcessorGroup*,
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
              "Doing AMRMPM::updateErosionParameter");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {

      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      // Get the localization info
      ParticleVariable<int> isLocalized;
      new_dw->allocateTemporary(isLocalized, pset);
      ParticleSubset::iterator iter = pset->begin();
      for (; iter != pset->end(); iter++) {
        isLocalized[*iter] = 0;
      }
      mpm_matl->getConstitutiveModel()->getDamageParameter(patch,
                                                           isLocalized,
                                                           dwi,
                                                           old_dw,
                                                           new_dw);
    }
  }
}
//______________________________________________________________________
//
void
AMRMPM::computeInternalForce(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing AMRMPM::computeInternalForce");

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0 / dx.x();
    oodx[1] = 1.0 / dx.y();
    oodx[2] = 1.0 / dx.z();
    Matrix3 Id;
    Id.Identity();

    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);

    string interp_type = d_mpmFlags->d_interpolatorType;

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      constParticleVariable<Point> pX;
      constParticleVariable<double> pVol;
      constParticleVariable<double> p_pressure;
      constParticleVariable<double> p_q;
      constParticleVariable<Matrix3> pStress;
      constParticleVariable<Matrix3> pSize;
      NCVariable<Vector> gIntForce;
      NCVariable<Matrix3> gStress;
      constNCVariable<double> gVolume;
      constParticleVariable<Matrix3> pDefGrad;

      ParticleSubset* pset = old_dw->getParticleSubset(dwi,
                                                       patch,
                                                       Ghost::AroundNodes,
                                                       d_numGhostParticles,
                                                       d_mpmLabels->pXLabel);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pVol, d_mpmLabels->pVolumeLabel, pset);
      old_dw->get(pStress, d_mpmLabels->pStressLabel, pset);
      new_dw->get(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);
      old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);

      new_dw
        ->get(gVolume, d_mpmLabels->gVolumeLabel, dwi, patch, Ghost::None, 0);
      new_dw->allocateAndPut(gStress,
                             d_mpmLabels->gStressForSavingLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gIntForce,
                             d_mpmLabels->gInternalForceLabel,
                             dwi,
                             patch);
      gStress.initialize(Matrix3(0));
      gIntForce.initialize(Vector(0, 0, 0));

      // load p_q
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

      /*`==========TESTING==========*/
      NCVariable<double> gSumS;
      new_dw->allocateAndPut(gSumS, gSumSLabel, dwi, patch);
      gSumS.initialize(0);
      /*===========TESTING==========`*/

      //__________________________________
      //  fine Patch
      gStress.initialize(Matrix3(0));

      Matrix3 stressvol;
      Matrix3 stresspress;
      int n8or27 = d_mpmFlags->d_8or27;
      std::vector<IntVector> ni(interpolator->size());
      std::vector<double> S(interpolator->size());
      std::vector<Vector> d_S(interpolator->size());

      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeightsAndShapeDerivatives(pX[idx],
                                                            ni,
                                                            S,
                                                            d_S,
                                                            pSize[idx],
                                                            pDefGrad[idx]);

        stresspress = pStress[idx] + Id * (/*p_pressure*/ -p_q[idx]);

        for (int k = 0; k < n8or27; k++) {

          if (patch->containsNode(ni[k])) {
            Vector div(d_S[k].x() * oodx[0],
                       d_S[k].y() * oodx[1],
                       d_S[k].z() * oodx[2]);

            gIntForce[ni[k]] -= (div * stresspress) * pVol[idx];

            // std::cout << " CIF: ni: " << ni[k] << " div " << div << "\t
            // internalForce " << gIntForce[ni[k]] << std::endl;
            // std::cout << " div " << div[k] << " stressPress: " << stresspress
            // << endl;

            if (std::isinf(gIntForce[ni[k]].length()) ||
                std::isnan(gIntForce[ni[k]].length())) {
              std::cout << "INF: " << ni[k] << " " << gIntForce[ni[k]]
                        << " div: " << div << " stressPress: " << stresspress
                        << " pVol " << pVol[idx] << std::endl;
            }
            /*`==========TESTING==========*/
            gSumS[ni[k]] += S[k];
            /*===========TESTING==========`*/
          }
        }
      }

      string interp_type = d_mpmFlags->d_interpolatorType;
      MPMBoundCond bc;
      bc.setBoundaryCondition(patch, dwi, "Symmetric", gIntForce, interp_type);
    } // End matl loop
  }   // End patch loop
}

//______________________________________________________________________
//
void
AMRMPM::computeInternalForce_CFI(const ProcessorGroup*,
                                 const PatchSubset* finePatches,
                                 const MaterialSubset*,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  const Level* fineLevel   = getLevel(finePatches);
  const Level* coarseLevel = fineLevel->getCoarserLevel().get_rep();
  IntVector refineRatio(fineLevel->getRefinementRatio());

  for (int p = 0; p < finePatches->size(); p++) {
    const Patch* finePatch = finePatches->get(p);
    printTask(finePatches,
              finePatch,
              cout_doing,
              "Doing AMRMPM::computeInternalForce_CFI");

    auto interpolator = d_mpmFlags->d_interpolator->clone(finePatch);

    //__________________________________
    //          AT CFI
    if (fineLevel->hasCoarserLevel() && finePatch->hasCoarseFaces()) {

      // Determine extents for coarser level particle data
      // Linear Interpolation:  1 layer of coarse level cells
      // Gimp Interpolation:    2 layers
      /*`==========TESTING==========*/
      IntVector nLayers(d_nPaddingCells_Coarse,
                        d_nPaddingCells_Coarse,
                        d_nPaddingCells_Coarse);
      IntVector nPaddingCells = nLayers * (fineLevel->getRefinementRatio());
      // std::cout << " nPaddingCells " << nPaddingCells << "nLayers " <<
      // nLayers << endl;
      /*===========TESTING==========`*/

      int nGhostCells           = 0;
      bool returnExclusiveRange = false;
      IntVector cl_tmp, ch_tmp, fl, fh;

      getCoarseLevelRange(finePatch,
                          coarseLevel,
                          cl_tmp,
                          ch_tmp,
                          fl,
                          fh,
                          nPaddingCells,
                          nGhostCells,
                          returnExclusiveRange);

      //  expand cl_tmp when a neighor patch exists.
      //  This patch owns the low nodes.  You need particles
      //  from the neighbor patch.
      cl_tmp -= finePatch->neighborsLow() * nLayers;

      // find the coarse patches under the fine patch.
      // You must add a single layer of padding cells.
      int padding = 1;
      Level::selectType coarsePatches;
      finePatch->getOtherLevelPatches(-1, coarsePatches, padding);

      Matrix3 Id;
      Id.Identity();

      constNCVariable<Stencil7> zoi_fine;
      new_dw
        ->get(zoi_fine, d_mpmLabels->gZOILabel, 0, finePatch, Ghost::None, 0);

      int numMPMMatls = d_materialManager->getNumMaterials("MPM");

      for (int m = 0; m < numMPMMatls; m++) {
        MPMMaterial* mpm_matl =
          static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
        int dwi = mpm_matl->getDWIndex();

        NCVariable<Vector> gIntForce;
        new_dw->getModifiable(gIntForce,
                              d_mpmLabels->gInternalForceLabel,
                              dwi,
                              finePatch);

        /*`==========TESTING==========*/
        NCVariable<double> gSumS;
        new_dw->getModifiable(gSumS, gSumSLabel, dwi, finePatch);
        /*===========TESTING==========`*/

        // loop over the coarse patches under the fine patches.
        for (size_t cp = 0; cp < coarsePatches.size(); cp++) {
          const Patch* coarsePatch = coarsePatches[cp];

          // get coarse level particle data
          ParticleSubset* pset_coarse;
          constParticleVariable<Point> pX_coarse;
          constParticleVariable<Matrix3> pStress_coarse;
          constParticleVariable<double> pVol_coarse;
          constParticleVariable<double> p_q_coarse;

          // coarseLow and coarseHigh cannot lie outside of the coarse patch
          IntVector cl = Max(cl_tmp, coarsePatch->getCellLowIndex());
          IntVector ch = Min(ch_tmp, coarsePatch->getCellHighIndex());

          pset_coarse = old_dw->getParticleSubset(dwi,
                                                  cl,
                                                  ch,
                                                  coarsePatch,
                                                  d_mpmLabels->pXLabel);

#if 0
          std::cout << " fine patch : " << finePatch->getGridIndex() << std::endl;
          std::cout << " cl_tmp: "<< cl_tmp << " ch_tmp: " << ch_tmp << std::endl;
          std::cout << " cl:     " << cl    << " ch:     " << ch<< " fl: " << fl << " fh " << fh << std::endl;                                                     
          std::cout << "  " << *pset_coarse << std::endl;
#endif

          // coarse level data
          old_dw->get(pX_coarse, d_mpmLabels->pXLabel, pset_coarse);
          old_dw->get(pVol_coarse, d_mpmLabels->pVolumeLabel, pset_coarse);
          old_dw->get(pStress_coarse, d_mpmLabels->pStressLabel, pset_coarse);

          // Artificial Viscosity
          if (d_mpmFlags->d_artificialViscosity) {
            old_dw->get(p_q_coarse, d_mpmLabels->p_qLabel, pset_coarse);
          } else {
            ParticleVariable<double> p_q_create;
            new_dw->allocateTemporary(p_q_create, pset_coarse);
            for (ParticleSubset::iterator it = pset_coarse->begin();
                 it != pset_coarse->end();
                 it++) {
              p_q_create[*it] = 0.0;
            }
            p_q_coarse = p_q_create; // reference created data
          }

          //__________________________________
          //  Iterate over the coarse level particles and
          // add their contribution to the internal stress on the fine patch
          for (ParticleSubset::iterator iter = pset_coarse->begin();
               iter != pset_coarse->end();
               iter++) {
            particleIndex idx = *iter;

            std::vector<IntVector> ni;
            std::vector<double> S;
            std::vector<Vector> div;
            interpolator->findCellAndWeightsAndShapeDerivatives_CFI(
              pX_coarse[idx],
              ni,
              S,
              div,
              zoi_fine);

            Matrix3 stresspress = pStress_coarse[idx] + Id * (-p_q_coarse[idx]);

            IntVector fineNode;
            for (int k = 0; k < (int)ni.size(); k++) {
              fineNode = ni[k];

              if (finePatch->containsNode(fineNode)) {
                gSumS[fineNode] += S[k];

                Vector Increment((div[k] * stresspress) * pVol_coarse[idx]);
                // Vector Before = gIntForce[fineNode];
                // Vector After  = Before - Increment;

                gIntForce[fineNode] -= Increment;

                //  std::cout << " CIF_CFI: ni: " << ni[k] << " div " << div[k]
                //  <<
                //  "\t internalForce " << gIntForce[fineNode] << std::endl;
                //  std::cout << "    before " << Before << " After " << After
                //  << " Increment " << Increment << std::endl; std::cout << "
                //  div " << div[k] << " stressPress: " << stresspress << "
                //  pVol_coarse " << pVol_coarse[idx] << std::endl;

                /*`==========TESTING==========*/
                if (std::isinf(gIntForce[fineNode].length()) ||
                    std::isnan(gIntForce[fineNode].length())) {
                  std::cout << "INF: " << fineNode << " " << gIntForce[fineNode]
                            << " div[k]:" << div[k]
                            << " stressPress: " << stresspress << " pVol "
                            << pVol_coarse[idx] << std::endl;
                }
#if 0             
                if( gIntForce[fineNode].length()  >1e-10){
                  std::cout << "CIF_CFI: " << fineNode
                       << "    L-"<< getLevel(finePatches)->getIndex()
                       <<" InternalForce " << gIntForce[fineNode] << " div[k]: " << div[k] << " stressPress: " << stresspress 
                       << " pVol " << pVol_coarse[idx] << std::endl;
                  std::cout << "          Before: " << Before << " Increment " << Increment << std::endl;
                }
#endif
                /*===========TESTING==========`*/
              } // contains node
            }   // node loop
          }     // pset loop
        }       // coarse Patch loop

        //__________________________________
        //  Set boundary conditions
        string interp_type = d_mpmFlags->d_interpolatorType;
        MPMBoundCond bc;
        bc.setBoundaryCondition(finePatch,
                                dwi,
                                "Symmetric",
                                gIntForce,
                                interp_type);

      } // End matl loop
    }   // patch has CFI faces
  }     // End fine patch loop
}

//______________________________________________________________________
//
void
AMRMPM::computeAndIntegrateAcceleration(const ProcessorGroup*,
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
              "Doing AMRMPM::computeAndIntegrateAcceleration");

    Ghost::GhostType gnone = Ghost::None;
    Vector gravity         = d_mpmFlags->d_gravity;

    for (size_t m = 0; m < d_materialManager->getNumMaterials("MPM"); m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      // Get required variables for this patch
      constNCVariable<Vector> gIntForce;
      constNCVariable<Vector> gExtForce;
      constNCVariable<Vector> gVelocity;
      constNCVariable<double> gMass;
      constNCVariable<double> gConcentration, gConcNoBC, gExtScalarFlux;

      delt_vartype delT;
      old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

      new_dw->get(gIntForce,
                  d_mpmLabels->gInternalForceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(gExtForce,
                  d_mpmLabels->gExternalForceLabel,
                  dwi,
                  patch,
                  gnone,
                  0);
      new_dw->get(gMass, d_mpmLabels->gMassLabel, dwi, patch, gnone, 0);
      new_dw->get(gVelocity, d_mpmLabels->gVelocityLabel, dwi, patch, gnone, 0);

      // Create variables for the results
      NCVariable<Vector> gVelocity_star;
      NCVariable<Vector> gAcceleration;
      NCVariable<double> gConcStar, gConcRate;
      new_dw->allocateAndPut(gVelocity_star,
                             d_mpmLabels->gVelocityStarLabel,
                             dwi,
                             patch);
      new_dw->allocateAndPut(gAcceleration,
                             d_mpmLabels->gAccelerationLabel,
                             dwi,
                             patch);

      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->get(gConcentration,
                    d_mpmLabels->gConcentrationLabel,
                    dwi,
                    patch,
                    gnone,
                    0);
        new_dw->get(gConcNoBC,
                    d_mpmLabels->gConcentrationNoBCLabel,
                    dwi,
                    patch,
                    gnone,
                    0);
        new_dw->get(gExtScalarFlux,
                    d_mpmLabels->gExternalScalarFluxLabel,
                    dwi,
                    patch,
                    gnone,
                    0);

        new_dw->getModifiable(gConcRate,
                              d_mpmLabels->diffusion->gConcentrationRate,
                              dwi,
                              patch);
        new_dw->allocateAndPut(gConcStar,
                               d_mpmLabels->gConcentrationStarLabel,
                               dwi,
                               patch);
      }

      gAcceleration.initialize(Vector(0., 0., 0.));
      double damp_coef = d_mpmFlags->d_artificialDampCoeff;
      gVelocity_star.initialize(Vector(0., 0., 0.));

      for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
           iter++) {
        IntVector n = *iter;

        Vector acc(0, 0, 0);
        if (gMass[n] > d_mpmFlags->d_minMassForAcceleration) {
          acc = (gIntForce[n] + gExtForce[n]) / gMass[n];
          acc -= damp_coef * gVelocity[n];
        }
        gAcceleration[n]  = acc + gravity;
        gVelocity_star[n] = gVelocity[n] + gAcceleration[n] * delT;

/*`==========TESTING==========*/
#ifdef DEBUG_ACC
        if (abs(gAcceleration[n].length() - d_acc_ans.length()) > d_acc_tol) {
          Vector diff = gAcceleration[n] - d_acc_ans;
          std::cout << "    L-" << getLevel(patches)->getIndex()
                    << " node: " << n << " gAcceleration: " << gAcceleration[n]
                    << " gExtForce: " << gExtForce[n]
                    << " gIntForce: " << gIntForce[n] << " diff: " << diff
                    << " gMass: " << gMass[n] << " gravity: " << gravity
                    << std::endl;
        }
#endif
        /*===========TESTING==========`*/
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c = *iter;
          //          gConcRate[c] /= gMass[c];
          gConcStar[c] = gConcentration[c] + gConcRate[c] * delT;
        }

        MPMBoundCond bc;
        bc.setBoundaryCondition(patch,
                                dwi,
                                "SD-Type",
                                gConcStar,
                                d_mpmFlags->d_interpolatorType);

        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c = *iter;
          gConcRate[c] =
            (gConcStar[c] - gConcNoBC[c]) / delT + gExtScalarFlux[c];
        }
      } // if doScalarDiffusion
    }   // matls
  }     // patches
}
//______________________________________________________________________
//
void
AMRMPM::setGridBoundaryConditions(const ProcessorGroup*,
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
              "Doing AMRMPM::setGridBoundaryConditions");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    delt_vartype delT;
    old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

    string interp_type = d_mpmFlags->d_interpolatorType;

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      NCVariable<Vector> gVelocity_star;
      NCVariable<Vector> gAcceleration;
      constNCVariable<Vector> gVelocity;

      new_dw->getModifiable(gAcceleration,
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

      //__________________________________
      // Apply grid boundary conditions to velocity_star and acceleration
      if (patch->hasBoundaryFaces()) {
        IntVector node(0, 4, 4);

        MPMBoundCond bc;
        bc.setBoundaryCondition(patch,
                                dwi,
                                "Velocity",
                                gVelocity_star,
                                interp_type);
        bc.setBoundaryCondition(patch,
                                dwi,
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
      }

      //__________________________________
      //
      if (!d_mpmFlags->d_doGridReset) {
        NCVariable<Vector> displacement;
        constNCVariable<Vector> displacementOld;
        new_dw->allocateAndPut(displacement,
                               d_mpmLabels->gDisplacementLabel,
                               dwi,
                               patch);
        old_dw->get(displacementOld,
                    d_mpmLabels->gDisplacementLabel,
                    dwi,
                    patch,
                    Ghost::None,
                    0);
        for (NodeIterator iter = patch->getExtraNodeIterator(); !iter.done();
             iter++) {
          IntVector c     = *iter;
          displacement[c] = displacementOld[c] + gVelocity_star[c] * delT;
        }
      } // d_doGridReset

    } // matl loop
  }   // patch loop
}
//______________________________________________________________________

void
AMRMPM::computeZoneOfInfluence(const ProcessorGroup*,
                               const PatchSubset* patches,
                               const MaterialSubset*,
                               DataWarehouse* old_dw,
                               DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  ASSERT(level->hasCoarserLevel());

  //__________________________________
  //  Initialize the interior nodes
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();

    printTask(patches,
              patch,
              cout_doing,
              "Doing AMRMPM::computeZoneOfInfluence");
    NCVariable<Stencil7> zoi;
    new_dw->allocateAndPut(zoi, d_mpmLabels->gZOILabel, 0, patch);

    for (NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      zoi[c].p    = -9876543210e99;
      zoi[c].w    = dx.x();
      zoi[c].e    = dx.x();
      zoi[c].s    = dx.y();
      zoi[c].n    = dx.y();
      zoi[c].b    = dx.z();
      zoi[c].t    = dx.z();
    }
  }

  //__________________________________
  // Set the ZOI on the current level.
  // Look up at the finer level patches
  // for coarse-fine interfaces
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    NCVariable<Stencil7> zoi;
    new_dw->getModifiable(zoi, d_mpmLabels->gZOILabel, 0, patch);

    if (level->hasFinerLevel()) {
      const Level* fineLevel = level->getFinerLevel().get_rep();

      Level::selectType finePatches;
      patch->getFineLevelPatches(finePatches);

      for (size_t p = 0; p < finePatches.size(); p++) {
        const Patch* finePatch = finePatches[p];

        Vector fine_dx = finePatch->dCell();

        //__________________________________
        // Iterate over coarsefine interface faces
        if (finePatch->hasCoarseFaces()) {
          std::vector<Patch::FaceType> cf;
          finePatch->getCoarseFaces(cf);

          std::vector<Patch::FaceType>::const_iterator iter;
          for (iter = cf.begin(); iter != cf.end(); ++iter) {
            Patch::FaceType patchFace = *iter;

            // determine the iterator on the coarse level.
            NodeIterator n_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
            bool isRight_CP_FP_pair;

            coarseLevel_CFI_NodeIterator(patchFace,
                                         patch,
                                         finePatch,
                                         fineLevel,
                                         n_iter,
                                         isRight_CP_FP_pair);

            // The ZOI element is opposite
            // of the patch face
            int element = patchFace;
            if (patchFace == Patch::xminus || patchFace == Patch::yminus ||
                patchFace == Patch::zminus) {
              element += 1; // e, n, t
            }
            if (patchFace == Patch::xplus || patchFace == Patch::yplus ||
                patchFace == Patch::zplus) {
              element -= 1; // w, s, b
            }
            IntVector dir = patch->getFaceAxes(patchFace); // face axes
            int p_dir     = dir[0];                        // normal direction

            // eject if this is not the right coarse/fine patch pair
            if (isRight_CP_FP_pair) {

              //              std::cout << "  A) Setting ZOI  "
              //                   << " \t On L-" << level->getIndex() << "
              //                   patch  " << patch->getID()
              //                   << ", beneath patch " << finePatch->getID()
              //                   << ", face: "  <<
              //                   finePatch->getFaceName(patchFace)
              //                   << ", isRight_CP_FP_pair: " <<
              //                   isRight_CP_FP_pair  << " n_iter: " << n_iter
              //                   << std::endl;

              for (; !n_iter.done(); n_iter++) {
                IntVector c     = *n_iter;
                zoi[c][element] = fine_dx[p_dir];
              }
            }

          } // patch face loop
        }   // hasCoarseFaces
      }     // finePatches loop
    }       // has finer level
  }         // patches loop

  //__________________________________
  // set the ZOI in cells in which there are overlaping coarse level nodes
  // look down for coarse level patches
  for (int p = 0; p < patches->size(); p++) {
    const Patch* finePatch = patches->get(p);
    NCVariable<Stencil7> zoi_fine;
    new_dw->getModifiable(zoi_fine, d_mpmLabels->gZOILabel, 0, finePatch);

    // underlying coarse level
    if (level->hasCoarserLevel()) {
      Level::selectType coarsePatches;
      finePatch->getOtherLevelPatchesNB(-1, coarsePatches, 0);
      //__________________________________
      // Iterate over coarsefine interface faces
      if (finePatch->hasCoarseFaces()) {
        std::vector<Patch::FaceType> cf;
        finePatch->getCoarseFaces(cf);

        std::vector<Patch::FaceType>::const_iterator iter;
        for (iter = cf.begin(); iter != cf.end(); ++iter) {
          Patch::FaceType patchFace = *iter;
          bool setFace              = false;

          for (size_t p = 0; p < coarsePatches.size(); p++) {
            const Patch* coarsePatch = coarsePatches[p];
            Vector coarse_dx         = coarsePatch->dCell();

            // determine the iterator on the coarse level.
            NodeIterator n_iter(IntVector(-8, -8, -8), IntVector(-9, -9, -9));
            bool isRight_CP_FP_pair;

            fineLevel_CFI_NodeIterator(patchFace,
                                       coarsePatch,
                                       finePatch,
                                       n_iter,
                                       isRight_CP_FP_pair);

            int element   = patchFace;
            IntVector dir = finePatch->getFaceAxes(patchFace); // face axes
            int p_dir     = dir[0];                            // normal dir

            // Is this the right coarse/fine patch pair
            if (isRight_CP_FP_pair) {
              setFace = true;

              //              std::cout << "  B) Setting ZOI  "
              //                   << " \t On L-" << level->getIndex() << "
              //                   patch  " << finePatch->getID()
              //                   << "   CFI face: "  <<
              //                   finePatch->getFaceName(patchFace)
              //                   << " isRight_CP_FP_pair: " <<
              //                   isRight_CP_FP_pair  << " n_iter: " << n_iter
              //                   << std::endl;

              for (; !n_iter.done(); n_iter++) {
                IntVector c          = *n_iter;
                zoi_fine[c][element] = coarse_dx[p_dir];
              }
            }
          } // coarsePatches loop

          // bulletproofing
          if (!setFace) {
            std::ostringstream warn;
            warn << "\n ERROR: computeZoneOfInfluence:Fine Level: Did not find "
                    "node iterator! "
                 << "\n coarse: L-" << level->getIndex()
                 << "\n coarePatches size: " << coarsePatches.size()
                 << "\n fine patch:   " << *finePatch
                 << "\n fine patch face: " << finePatch->getFaceName(patchFace);
            throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
          }
        } // face interator
      }   // patch has coarse face
    }     // has finer level
  }       // patch loop
}

//______________________________________________________________________
//
void
AMRMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
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
              "Doing AMRMPM::interpolateToParticlesAndUpdate");

    auto interpolator = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
    // Vector dx = patch->dCell();

    // Performs the interpolation from the cell vertices of the grid
    // acceleration and velocity to the particles to update their
    // velocity and position respectively

    // DON'T MOVE THESE!!!
    double thermal_energy = 0.0;
    double totalmass      = 0;
    Vector CMX(0.0, 0.0, 0.0);
    Vector totalMom(0.0, 0.0, 0.0);
    double ke = 0;

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    delt_vartype delT;
    old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

    double move_particles = 1.;
    if (!d_mpmFlags->d_doGridReset) {
      move_particles = 0.;
    }

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<Point> pX;
      ParticleVariable<Point> pXnew, pXx;
      constParticleVariable<Vector> pVelocity;
      constParticleVariable<Matrix3> pSize;
      ParticleVariable<Vector> pVelocitynew, pAcceleration_new;
      ParticleVariable<Matrix3> pSizeNew;
      constParticleVariable<double> pMass, pTemperature;
      ParticleVariable<double> pMassNew, pVolume, pTempNew;
      constParticleVariable<long64> pIDs;
      ParticleVariable<long64> pIDs_new;
      constParticleVariable<Vector> pDisp;
      ParticleVariable<Vector> pDispnew;
      ParticleVariable<double> pTempPreNew;
      constParticleVariable<Matrix3> pFOld;

      constParticleVariable<double> pConcentration;
      ParticleVariable<double> pConcentrationNew;
      ParticleVariable<double> pConcPreviousNew;
      constNCVariable<double> gConcentrationRate;

      // Get the arrays of grid data on which the new particle values depend
      constNCVariable<Vector> gVelocity_star, gAcceleration;
      constNCVariable<double> gTemperatureRate;
      constNCVariable<double> dTdt, frictionTempRate;
      double Cp = mpm_matl->getSpecificHeat();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      old_dw->get(pX, d_mpmLabels->pXLabel, pset);
      old_dw->get(pDisp, d_mpmLabels->pDispLabel, pset);
      old_dw->get(pMass, d_mpmLabels->pMassLabel, pset);
      old_dw->get(pVelocity, d_mpmLabels->pVelocityLabel, pset);
      old_dw->get(pTemperature, d_mpmLabels->pTemperatureLabel, pset);
      old_dw->get(pFOld, d_mpmLabels->pDefGradLabel, pset);
      new_dw->get(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);

      new_dw->allocateAndPut(pVelocitynew,
                             d_mpmLabels->pVelocityLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pAcceleration_new,
                             d_mpmLabels->pAccelerationLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pXnew, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->allocateAndPut(pXx, d_mpmLabels->pXXLabel, pset);
      new_dw->allocateAndPut(pDispnew, d_mpmLabels->pDispLabel_preReloc, pset);
      new_dw->allocateAndPut(pMassNew, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->allocateAndPut(pTempNew,
                             d_mpmLabels->pTemperatureLabel_preReloc,
                             pset);
      new_dw->allocateAndPut(pTempPreNew,
                             d_mpmLabels->pTempPreviousLabel_preReloc,
                             pset);

      // Copy needed for switch from explicit to implicit MPM
      constParticleVariable<double> pExtHeatFlux;
      ParticleVariable<double> pExtHeatFlux_new;
      old_dw->get(pExtHeatFlux, d_mpmLabels->pExternalHeatFluxLabel, pset);
      new_dw->allocateAndPut(pExtHeatFlux_new,
                             d_mpmLabels->pExternalHeatFluxLabel_preReloc,
                             pset);
      pExtHeatFlux_new.copyData(pExtHeatFlux);

      Ghost::GhostType gac = Ghost::AroundCells;
      if (d_mpmFlags->d_doScalarDiffusion) {
        old_dw->get(pConcentration,
                    d_mpmLabels->diffusion->pConcentration,
                    pset);
        new_dw->get(gConcentrationRate,
                    d_mpmLabels->diffusion->gConcentrationRate,
                    dwi,
                    patch,
                    Ghost::AroundCells,
                    d_numGhostParticles);
        new_dw->allocateAndPut(pConcentrationNew,
                               d_mpmLabels->diffusion->pConcentration_preReloc,
                               pset);
        new_dw->allocateAndPut(pConcPreviousNew,
                               d_mpmLabels->pConcPreviousLabel_preReloc,
                               pset);
      }

      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      // Carry forward ParticleID and pSize
      old_dw->get(pIDs, d_mpmLabels->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pIDs_new,
                             d_mpmLabels->pParticleIDLabel_preReloc,
                             pset);
      pIDs_new.copyData(pIDs);

      new_dw->get(gVelocity_star,
                  d_mpmLabels->gVelocityStarLabel,
                  dwi,
                  patch,
                  Ghost::AroundCells,
                  d_numGhostParticles);
      new_dw->get(gAcceleration,
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
      new_dw->get(frictionTempRate,
                  d_mpmLabels->frictionalWorkLabel,
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
      } else {
        NCVariable<double> dTdt_create, massBurnFrac_create;
        new_dw->allocateTemporary(dTdt_create, patch, Ghost::AroundCells, d_numGhostParticles);
        dTdt_create.initialize(0.);
        dTdt = dTdt_create; // reference created data
      }

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;

        // Get the node indices that surround the cell
        interpolator->findCellAndWeights(pX[idx],
                                         ni,
                                         S,
                                         pSize[idx],
                                         pFOld[idx]);

        Vector vel(0.0, 0.0, 0.0);
        Vector acc(0.0, 0.0, 0.0);
        double fricTempRate = 0.0;
        double tempRate     = 0.0;
        double concRate     = 0.0;

        // Accumulate the contribution from vertices on this level
        for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
          IntVector node = ni[k];
          vel += gVelocity_star[node] * S[k];
          acc += gAcceleration[node] * S[k];

          fricTempRate = frictionTempRate[node] * d_mpmFlags->d_addFrictionWork;
          tempRate +=
            (gTemperatureRate[node] + dTdt[node] + fricTempRate) * S[k];
        }

        // Update the particle's position and velocity
        pXnew[idx]             = pX[idx] + vel * delT * move_particles;
        pDispnew[idx]          = pDisp[idx] + vel * delT;
        pVelocitynew[idx]      = pVelocity[idx] + acc * delT;
        pAcceleration_new[idx] = acc;

        // pXx is only useful if we're not in normal grid resetting mode.
        pXx[idx]         = pX[idx] + pDispnew[idx];
        pTempNew[idx]    = pTemperature[idx] + tempRate * delT;
        pTempPreNew[idx] = pTemperature[idx]; // for thermal stress
        pMassNew[idx]    = pMass[idx];

        if (d_mpmFlags->d_doScalarDiffusion) {
          for (int k = 0; k < d_mpmFlags->d_8or27; k++) {
            IntVector node = ni[k];
            concRate += gConcentrationRate[node] * S[k];
          }

          pConcentrationNew[idx] = pConcentration[idx] + concRate * delT;
          pConcPreviousNew[idx]  = pConcentration[idx];
        }
/*`==========TESTING==========*/
#ifdef DEBUG_VEL
        Vector diff = (pVelocitynew[idx] - d_vel_ans);
        if (abs(diff.length()) > d_vel_tol) {
          std::cout << "    L-" << getLevel(patches)->getIndex()
                    << " pX: " << pXnew[idx]
                    << " pVelocitynew: " << pVelocitynew[idx] << " pVelocity "
                    << pVelocity[idx] << " diff " << diff << std::endl;
        }
#endif
#ifdef DEBUG_ACC
#endif
        /*===========TESTING==========`*/

        totalmass += pMass[idx];
        thermal_energy += pTemperature[idx] * pMass[idx] * Cp;
        ke += .5 * pMass[idx] * pVelocitynew[idx].length2();
        CMX = CMX + (pXnew[idx] * pMass[idx]).asVector();
        totalMom += pVelocitynew[idx] * pMass[idx];
      }

#if 0 // Until Todd is ready for this, leave inactive
      // Delete particles that have left the domain
      // This is only needed if extra cells are being used.
      // Also delete particles whose mass is too small (due to combustion)
      // For particles whose new velocity exceeds a maximum set in the input
      // file, set their velocity back to the velocity that it came into
      // this step with
      for(ParticleSubset::iterator iter  = pset->begin();
          iter != pset->end(); iter++){
        particleIndex idx = *iter;
        if ((pMassNew[idx] <= d_mpmFlags->d_minPartMass) || pTempNew[idx] < 0. ||
            (pLocalized[idx]==-999)){
          delset->addParticle(idx);
//        std::cout << "Material = " << m << " Deleted Particle = " << pIDs_new[idx] 
//             << " xold = " << pX[idx] << " xnew = " << pXnew[idx]
//             << " vold = " << pVelocity[idx] << " vnew = "<< pVelocitynew[idx]
//             << " massold = " << pMass[idx] << " massnew = " << pMassNew[idx]
//             << " tempold = " << pTemperature[idx] 
//             << " tempnew = " << pTempNew[idx]
//             << " pLocalized = " << pLocalized[idx]
//             << " volnew = " << pVolume[idx] << std::endl;
        }
        if(pVelocitynew[idx].length() > d_mpmFlags->d_max_vel){
          if(pVelocitynew[idx].length() >= pVelocity[idx].length()){
            pVelocitynew[idx]=(pVelocitynew[idx]/pVelocitynew[idx].length())*(d_mpmFlags->d_max_vel*.9);  
            cout<<endl<<"Warning: particle "<<pIDs[idx]<<" hit speed ceiling #1. Modifying particle velocity accordingly."<<endl;
            //pVelocitynew[idx]=pVelocity[idx];
          }
        }
      }
#endif

      new_dw->deleteParticles(delset);

      new_dw->put(sum_vartype(totalmass), d_mpmLabels->TotalMassLabel);
      new_dw->put(sum_vartype(ke), d_mpmLabels->KineticEnergyLabel);
      new_dw->put(sum_vartype(thermal_energy), d_mpmLabels->ThermalEnergyLabel);
      new_dw->put(sumvec_vartype(CMX), d_mpmLabels->CenterOfMassPositionLabel);
      new_dw->put(sumvec_vartype(totalMom), d_mpmLabels->TotalMomentumLabel);
#ifndef USE_DEBUG_TASK
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
#endif
      if (d_mpmFlags->d_refineParticles) {
        constParticleVariable<int> pRefinedOld;
        ParticleVariable<int> pRefinedNew;
        old_dw->get(pRefinedOld, d_mpmLabels->pRefinedLabel, pset);
        new_dw->allocateAndPut(pRefinedNew,
                               d_mpmLabels->pRefinedLabel_preReloc,
                               pset);
        pRefinedNew.copyData(pRefinedOld);
      }
    }
  }
}

void
AMRMPM::finalParticleUpdate(const ProcessorGroup*,
                            const PatchSubset* patches,
                            const MaterialSubset*,
                            DataWarehouse* old_dw,
                            DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing finalParticleUpdate");

    delt_vartype delT;
    old_dw->get(delT, d_mpmLabels->delTLabel, getLevel(patches));

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();
      // Get the arrays of particle values to be changed
      constParticleVariable<double> pdTdt, pMassNew;
      ParticleVariable<double> pTempNew;

      ParticleSubset* pset   = old_dw->getParticleSubset(dwi, patch);
      ParticleSubset* delset = scinew ParticleSubset(0, dwi, patch);

      new_dw->get(pdTdt, d_mpmLabels->pdTdtLabel, pset);
      new_dw->get(pMassNew, d_mpmLabels->pMassLabel_preReloc, pset);

      new_dw->getModifiable(pTempNew,
                            d_mpmLabels->pTemperatureLabel_preReloc,
                            pset);

      // Loop over particles
      for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
           iter++) {
        particleIndex idx = *iter;
        pTempNew[idx] += pdTdt[idx] * delT;

        // Delete particles whose mass is too small (due to combustion),
        // whose pLocalized flag has been set to -999 or who have a negative
        // temperature
        if ((pMassNew[idx] <= d_mpmFlags->d_minPartMass) ||
            pTempNew[idx] < 0.) {
          delset->addParticle(idx);
        }

      } // particles
      new_dw->deleteParticles(delset);
    } // materials
  }   // patches
}

void
AMRMPM::addParticles(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset*,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing addParticles");
    int numMPMMatls = d_materialManager->getNumMaterials("MPM");

    const Level* level = getLevel(patches);
    int levelIndex     = level->getIndex();
    bool hasCoarser    = false;
    if (level->hasCoarserLevel()) {
      hasCoarser = true;
    }

    // Carry forward CellNAPID
    constCCVariable<short int> NAPID;
    CCVariable<short int> NAPID_new;
    Ghost::GhostType gnone = Ghost::None;
    old_dw->get(NAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch, gnone, 0);
    new_dw->allocateAndPut(NAPID_new, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    NAPID_new.copyData(NAPID);

    // Mark cells where particles are refined for grid refinement
    CCVariable<double> refineCell;
    new_dw->getModifiable(refineCell,
                          d_mpmLabels->MPMRefineCellLabel,
                          0,
                          patch);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      ParticleVariable<Point> pX;
      ParticleVariable<Matrix3> pF, pSize, pStress, pVelgrad, pScalefac;
      ParticleVariable<long64> pIDs;
      ParticleVariable<double> pVolume, pMass, pTemp, pTempP, pColor, pConc,
        pConcpre;
      ParticleVariable<double> pESF;
      ParticleVariable<Vector> pVelocity, pAcc, pExtforce, pDisp, pConcgrad;
      ParticleVariable<int> pRef, ploc, pLaL, pRefOld, pLoadCID;
      new_dw->getModifiable(pX, d_mpmLabels->pXLabel_preReloc, pset);
      new_dw->getModifiable(pIDs, d_mpmLabels->pParticleIDLabel_preReloc, pset);
      new_dw->getModifiable(pMass, d_mpmLabels->pMassLabel_preReloc, pset);
      new_dw->getModifiable(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);
      new_dw->getModifiable(pDisp, d_mpmLabels->pDispLabel_preReloc, pset);
      new_dw->getModifiable(pStress, d_mpmLabels->pStressLabel_preReloc, pset);
      new_dw->getModifiable(pVolume, d_mpmLabels->pVolumeLabel_preReloc, pset);
      new_dw->getModifiable(pVelocity,
                            d_mpmLabels->pVelocityLabel_preReloc,
                            pset);
      new_dw->getModifiable(pAcc,
                            d_mpmLabels->pAccelerationLabel_preReloc,
                            pset);
      new_dw->getModifiable(pScalefac,
                            d_mpmLabels->pScaleFactorLabel_preReloc,
                            pset);
      new_dw->getModifiable(pExtforce,
                            d_mpmLabels->pExtForceLabel_preReloc,
                            pset);
      new_dw->getModifiable(pTemp,
                            d_mpmLabels->pTemperatureLabel_preReloc,
                            pset);
      new_dw->getModifiable(pTempP,
                            d_mpmLabels->pTempPreviousLabel_preReloc,
                            pset);
      new_dw->getModifiable(pRef, d_mpmLabels->pRefinedLabel_preReloc, pset);
      new_dw->getModifiable(pLaL, d_mpmLabels->pLastLevelLabel_preReloc, pset);
      new_dw->getModifiable(pVelgrad,
                            d_mpmLabels->pVelGradLabel_preReloc,
                            pset);
      new_dw->getModifiable(pF, d_mpmLabels->pDefGradLabel_preReloc, pset);
      if (d_mpmFlags->d_withColor) {
        new_dw->getModifiable(pColor, d_mpmLabels->pColorLabel_preReloc, pset);
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->getModifiable(pConc,
                              d_mpmLabels->diffusion->pConcentration_preReloc,
                              pset);
        new_dw->getModifiable(pConcpre,
                              d_mpmLabels->pConcPreviousLabel_preReloc,
                              pset);
        new_dw->getModifiable(pConcgrad,
                              d_mpmLabels->pGradConcentrationLabel_preReloc,
                              pset);
        new_dw->getModifiable(pESF,
                              d_mpmLabels->pExternalScalarFluxLabel,
                              pset);
      }
      if (d_mpmFlags->d_useLoadCurves) {
        new_dw->getModifiable(pLoadCID,
                              d_mpmLabels->pLoadCurveIDLabel_preReloc,
                              pset);
      }

      // Body force acceleration
      ParticleVariable<double> pCoriolis;
      ParticleVariable<Vector> pBodyFAcc;
      new_dw->getModifiable(pCoriolis,
                            d_mpmLabels->pCoriolisImportanceLabel_preReloc,
                            pset);
      new_dw->getModifiable(pBodyFAcc,
                            d_mpmLabels->pBodyForceAccLabel_preReloc,
                            pset);

      new_dw->allocateTemporary(pRefOld, pset);

      int numNewPartNeeded = 0;
      // Put refinement criteria here
      const unsigned int origNParticles = pset->addParticles(0);
      for (unsigned int pp = 0; pp < origNParticles; ++pp) {
        pRefOld[pp] = pRef[pp];
        // Conditions to refine particle based on physical state
        // TODO:  Check below, should be < or <= in first conditional
        if (pRef[pp] < levelIndex && pStress[pp].Norm() > 1) {
          pRef[pp]++;
          numNewPartNeeded++;
        }
        if (pRef[pp] > pRefOld[pp] || pStress[pp].Norm() > 1) {
          IntVector c = level->getCellIndex(pX[pp]);
          if (patch->containsCell(c)) {
            refineCell[c] = 3.0; // Why did I use 3 here?  JG
          }
        } else {
          if (hasCoarser) { /* see comment below */
            IntVector c = level->getCellIndex(pX[pp]);
            if (patch->containsCell(c)) {
              refineCell[c] = -100.;
            }
          }
        }
        // Refine particle if it is too big relative to the cell size
        // of the level it is on.  Don't refine the grid.
        if (pRef[pp] < levelIndex) {
          pRef[pp]++;
          numNewPartNeeded++;
        }
      }
      numNewPartNeeded *= 8;

      /*  This tomfoolery is in place to keep refined regions that contain
          particles refined.  If a patch with particles coarsens, the particles
          on that patch disappear when the fine patch is deleted.  This
          prevents the deletion of those patches.  Ideally, we'd allow
          coarsening and relocate the orphan particles, but I don't know how to
          do that.  JG */
      bool keep_patch_refined = false;
      IntVector low           = patch->getCellLowIndex();
      IntVector high          = patch->getCellHighIndex();
      IntVector middle        = (low + high) / IntVector(2, 2, 2);

      for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        if (refineCell[c] < 0.0) {
          keep_patch_refined = true;
          refineCell[c]      = 0.0;
        }
      }
      if (keep_patch_refined == true) {
        refineCell[middle] = -100.0;
      }
      /*  End tomfoolery */

      const unsigned int oldNumPar = pset->addParticles(numNewPartNeeded);

      ParticleVariable<Point> pXtmp;
      ParticleVariable<Matrix3> pFtmp, pSizetmp, pStrstmp, pVgradtmp, pSFtmp;
      ParticleVariable<long64> pIDstmp;
      ParticleVariable<double> pVoltmp, pMasstmp, pTemptmp, pTempPtmp, pESFtmp;
      ParticleVariable<double> pColortmp, pConctmp, pConcpretmp;
      ParticleVariable<Vector> pVeltmp, pAcctmp, pExtFtmp, pDisptmp,
        pConcgradtmp;
      ParticleVariable<int> pReftmp, ploctmp, pLaLtmp, pLoadCIDtmp;
      new_dw->allocateTemporary(pIDstmp, pset);
      new_dw->allocateTemporary(pXtmp, pset);
      new_dw->allocateTemporary(pVoltmp, pset);
      new_dw->allocateTemporary(pVeltmp, pset);
      new_dw->allocateTemporary(pAcctmp, pset);
      new_dw->allocateTemporary(pSFtmp, pset);
      new_dw->allocateTemporary(pExtFtmp, pset);
      new_dw->allocateTemporary(pTemptmp, pset);
      new_dw->allocateTemporary(pTempPtmp, pset);
      new_dw->allocateTemporary(pFtmp, pset);
      new_dw->allocateTemporary(pSizetmp, pset);
      new_dw->allocateTemporary(pDisptmp, pset);
      new_dw->allocateTemporary(pStrstmp, pset);
      if (d_mpmFlags->d_withColor) {
        new_dw->allocateTemporary(pColortmp, pset);
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->allocateTemporary(pConctmp, pset);
        new_dw->allocateTemporary(pConcpretmp, pset);
        new_dw->allocateTemporary(pConcgradtmp, pset);
        new_dw->allocateTemporary(pESFtmp, pset);
      }
      if (d_mpmFlags->d_useLoadCurves) {
        new_dw->allocateTemporary(pLoadCIDtmp, pset);
      }
      new_dw->allocateTemporary(pMasstmp, pset);
      new_dw->allocateTemporary(pReftmp, pset);
      new_dw->allocateTemporary(pLaLtmp, pset);
      // new_dw->allocateTemporary(ploctmp,  pset);
      new_dw->allocateTemporary(pVgradtmp, pset);

      // Body force acceleration
      ParticleVariable<double> pCoriolis_tmp;
      ParticleVariable<Vector> pBodyFAcc_tmp;
      new_dw->allocateTemporary(pCoriolis_tmp, pset);
      new_dw->allocateTemporary(pBodyFAcc_tmp, pset);

      // copy data from old variables for particle IDs and the position vector
      for (unsigned int pp = 0; pp < oldNumPar; ++pp) {
        pIDstmp[pp]   = pIDs[pp];
        pXtmp[pp]     = pX[pp];
        pVoltmp[pp]   = pVolume[pp];
        pVeltmp[pp]   = pVelocity[pp];
        pAcctmp[pp]   = pAcc[pp];
        pSFtmp[pp]    = pScalefac[pp];
        pExtFtmp[pp]  = pExtforce[pp];
        pTemptmp[pp]  = pTemp[pp];
        pTempPtmp[pp] = pTempP[pp];
        pFtmp[pp]     = pF[pp];
        pSizetmp[pp]  = pSize[pp];
        pDisptmp[pp]  = pDisp[pp];
        pStrstmp[pp]  = pStress[pp];
        if (d_mpmFlags->d_withColor) {
          pColortmp[pp] = pColor[pp];
        }
        if (d_mpmFlags->d_useLoadCurves) {
          pLoadCIDtmp[pp] = pLoadCID[pp];
        }
        pMasstmp[pp] = pMass[pp];
        pReftmp[pp]  = pRef[pp];
        pLaLtmp[pp]  = pLaL[pp];
        // ploctmp[pp]  = ploc[pp];
        pVgradtmp[pp] = pVelgrad[pp];

        // Body force quantities
        pCoriolis_tmp[pp] = pCoriolis[pp];
        pBodyFAcc_tmp[pp] = pBodyFAcc[pp];
      }

      if (d_mpmFlags->d_doScalarDiffusion) {
        for (unsigned int pp = 0; pp < oldNumPar; ++pp) {
          pConctmp[pp]     = pConc[pp];
          pConcpretmp[pp]  = pConcpre[pp];
          pConcgradtmp[pp] = pConcgrad[pp];
          pESFtmp[pp]      = pESF[pp];
        }
      }

      Vector dx     = patch->dCell();
      int numRefPar = 0;
      for (unsigned int idx = 0; idx < oldNumPar; ++idx) {
        if (pRef[idx] != pRefOld[idx]) {
          IntVector c_orig;
          patch->findCell(pX[idx], c_orig);
          std::vector<Point> new_part_pos;

          Matrix3 dsize = (pF[idx] * pSize[idx] *
                           Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));

          // Find vectors to new particle locations, based on particle size and
          // deformation (patterned after CPDI interpolator code)
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

          std::cout << "OPP = " << pX[idx] << std::endl;
          for (int i = 0; i < 8; i++) {
            //        std::cout << "NPP = " << new_part_pos[i] << std::endl;
            if (!level->containsPoint(new_part_pos[i])) {
              Point anchor    = level->getAnchor();
              Point orig      = new_part_pos[i];
              new_part_pos[i] = Point(std::max(orig.x(), anchor.x()),
                                      std::max(orig.y(), anchor.y()),
                                      std::max(orig.z(), anchor.z()));
            }

            long64 cellID = ((long64)c_orig.x() << 16) |
                            ((long64)c_orig.y() << 32) |
                            ((long64)c_orig.z() << 48);

            short int& myCellNAPID = NAPID_new[c_orig];
            int new_index;
            if (i == 0) {
              new_index = idx;
            } else {
              new_index = oldNumPar + 7 * numRefPar + i;
            }
            //          std::cout << "new_index = " << new_index << std::endl;
            pIDstmp[new_index]  = (cellID | (long64)myCellNAPID);
            pXtmp[new_index]    = new_part_pos[i];
            pVoltmp[new_index]  = .125 * pVolume[idx];
            pMasstmp[new_index] = .125 * pMass[idx];
            pVeltmp[new_index]  = pVelocity[idx];
            pAcctmp[new_index]  = pAcc[idx];
            pSFtmp[new_index]   = 0.5 * pScalefac[idx];
            pExtFtmp[new_index] = pExtforce[idx];
            pFtmp[new_index]    = pF[idx];
            pSizetmp[new_index] = 0.5 * pSize[idx];
            pDisptmp[new_index] = pDisp[idx];
            pStrstmp[new_index] = pStress[idx];
            if (d_mpmFlags->d_withColor) {
              pColortmp[new_index] = pColor[idx];
            }
            if (d_mpmFlags->d_doScalarDiffusion) {
              pConctmp[new_index]     = pConc[idx];
              pConcpretmp[new_index]  = pConcpre[idx];
              pConcgradtmp[new_index] = pConcgrad[idx];
              pESFtmp[new_index]      = pESF[idx];
            }
            if (d_mpmFlags->d_useLoadCurves) {
              pLoadCIDtmp[new_index] = pLoadCID[idx];
            }
            pTemptmp[new_index]  = pTemp[idx];
            pTempPtmp[new_index] = pTempP[idx];
            pReftmp[new_index]   = pRef[idx];
            pLaLtmp[new_index]   = pLaL[idx];
            // ploctmp[new_index]    = ploc[idx];
            pVgradtmp[new_index] = pVelgrad[idx];

            // Body force quantities
            pCoriolis_tmp[new_index] = pCoriolis[idx];
            pBodyFAcc_tmp[new_index] = pBodyFAcc[idx];

            NAPID_new[c_orig]++;
          }
          numRefPar++;
        } // if particle flagged for refinement
      }   // for particles

      // put back temporary data
      new_dw->put(pIDstmp, d_mpmLabels->pParticleIDLabel_preReloc, true);
      new_dw->put(pXtmp, d_mpmLabels->pXLabel_preReloc, true);
      new_dw->put(pVoltmp, d_mpmLabels->pVolumeLabel_preReloc, true);
      new_dw->put(pVeltmp, d_mpmLabels->pVelocityLabel_preReloc, true);
      new_dw->put(pAcctmp, d_mpmLabels->pAccelerationLabel_preReloc, true);
      new_dw->put(pSFtmp, d_mpmLabels->pScaleFactorLabel_preReloc, true);
      new_dw->put(pExtFtmp, d_mpmLabels->pExtForceLabel_preReloc, true);
      new_dw->put(pMasstmp, d_mpmLabels->pMassLabel_preReloc, true);
      new_dw->put(pTemptmp, d_mpmLabels->pTemperatureLabel_preReloc, true);
      new_dw->put(pTempPtmp, d_mpmLabels->pTempPreviousLabel_preReloc, true);
      new_dw->put(pSizetmp, d_mpmLabels->pSizeLabel_preReloc, true);
      new_dw->put(pDisptmp, d_mpmLabels->pDispLabel_preReloc, true);
      new_dw->put(pStrstmp, d_mpmLabels->pStressLabel_preReloc, true);
      if (d_mpmFlags->d_withColor) {
        new_dw->put(pColortmp, d_mpmLabels->pColorLabel_preReloc, true);
      }
      if (d_mpmFlags->d_doScalarDiffusion) {
        new_dw->put(pConctmp,
                    d_mpmLabels->diffusion->pConcentration_preReloc,
                    true);
        new_dw->put(pConcpretmp,
                    d_mpmLabels->pConcPreviousLabel_preReloc,
                    true);
        new_dw->put(pConcgradtmp,
                    d_mpmLabels->pGradConcentrationLabel_preReloc,
                    true);
        new_dw->put(pESFtmp, d_mpmLabels->pExternalScalarFluxLabel, true);
      }
      if (d_mpmFlags->d_useLoadCurves) {
        new_dw->put(pLoadCIDtmp, d_mpmLabels->pLoadCurveIDLabel_preReloc, true);
      }
      new_dw->put(pFtmp, d_mpmLabels->pDefGradLabel_preReloc, true);
      new_dw->put(pReftmp, d_mpmLabels->pRefinedLabel_preReloc, true);
      new_dw->put(pLaLtmp, d_mpmLabels->pLastLevelLabel_preReloc, true);
      new_dw->put(pVgradtmp, d_mpmLabels->pVelGradLabel_preReloc, true);

      // Body force terms
      new_dw->put(pBodyFAcc_tmp,
                  d_mpmLabels->pBodyForceAccLabel_preReloc,
                  true);
      new_dw->put(pCoriolis_tmp,
                  d_mpmLabels->pCoriolisImportanceLabel_preReloc,
                  true);

      // put back temporary data
    } // for matls
  }   // for patches
}

void
AMRMPM::reduceFlagsExtents(const ProcessorGroup*,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{

  // Currently doing for levels > 0
  const Level* level     = getLevel(patches);
  int levelIndex         = level->getIndex();
  IntVector RR_thisLevel = level->getRefinementRatio();
  int numLevels          = level->getGrid()->numLevels();

  IntVector RR_RelToFinest = IntVector(1, 1, 1);
  if (level->hasFinerLevel()) {
    RR_RelToFinest = RR_thisLevel * (numLevels - levelIndex - 1);
  }

  //  std::cout << "rFE levelIndex = " << levelIndex << std::endl;
  //  std::cout << "RR_RelToFinest = " << RR_RelToFinest << std::endl;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing reduceFlagsExtents");

    // Mark cells where particles are refined for grid refinement
    Ghost::GhostType gnone = Ghost::None;
    constCCVariable<double> refineCell;
    new_dw
      ->get(refineCell, d_mpmLabels->MPMRefineCellLabel, 0, patch, gnone, 0);

    int xmax, xmin, ymax, ymin, zmax, zmin;
    xmax = -999;
    ymax = -999;
    zmax = -999;
    xmin = 999999;
    ymin = 999999;
    zmin = 999999;
    //    int print = 0;
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      if (refineCell[c] > 0) {
        xmax = std::max(xmax, c.x());
        ymax = std::max(ymax, c.y());
        zmax = std::max(zmax, c.z());
        xmin = std::min(xmin, c.x());
        ymin = std::min(ymin, c.y());
        zmin = std::min(zmin, c.z());
        //        print = 1;
      }
    }

    xmax = xmax * RR_RelToFinest.x();
    ymax = ymax * RR_RelToFinest.y();
    zmax = zmax * RR_RelToFinest.z();
    xmin = xmin * RR_RelToFinest.x();
    ymin = ymin * RR_RelToFinest.y();
    zmin = zmin * RR_RelToFinest.z();

    /*
      if (print==1){
      std::cout << "Xmax = " << xmax << std::endl;
      std::cout << "Ymax = " << ymax << std::endl;
      std::cout << "Zmax = " << zmax << std::endl;
      std::cout << "Xmin = " << xmin << std::endl;
      std::cout << "Ymin = " << ymin << std::endl;
      std::cout << "Zmin = " << zmin << std::endl;
      }
    */

    new_dw->put(max_vartype(xmax), RefineFlagXMaxLabel);
    new_dw->put(max_vartype(ymax), RefineFlagYMaxLabel);
    new_dw->put(max_vartype(zmax), RefineFlagZMaxLabel);
    new_dw->put(min_vartype(xmin), RefineFlagXMinLabel);
    new_dw->put(min_vartype(ymin), RefineFlagYMinLabel);
    new_dw->put(min_vartype(zmin), RefineFlagZMinLabel);
  } // for patches
}

void
AMRMPM::computeParticleScaleFactor(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset*,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  // This task computes the particles initial physical size, to be used
  // in scaling particles for the deformed particle vis feature

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing computeParticleScaleFactor");

    int numMPMMatls = d_materialManager->getNumMaterials("MPM");
    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi              = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

      constParticleVariable<Matrix3> pSize, pF;
      ParticleVariable<Matrix3> pScaleFactor;
      new_dw->get(pSize, d_mpmLabels->pSizeLabel_preReloc, pset);
      new_dw->get(pF, d_mpmLabels->pDefGradLabel_preReloc, pset);
      new_dw->allocateAndPut(pScaleFactor,
                             d_mpmLabels->pScaleFactorLabel_preReloc,
                             pset);

      if (d_dataArchiver->isOutputTimestep()) {
        Vector dx = patch->dCell();
        for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
             iter++) {
          particleIndex idx = *iter;
          pScaleFactor[idx] = (pF[idx] * pSize[idx] *
                               Matrix3(dx[0], 0, 0, 0, dx[1], 0, 0, 0, dx[2]));
        } // for particles
      }   // isOutputTimestep
    }     // matls
  }       // patches
}

//______________________________________________________________________
//
void
AMRMPM::errorEstimate(const ProcessorGroup*,
                      const PatchSubset* patches,
                      const MaterialSubset* /*matls*/,
                      DataWarehouse*,
                      DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing AMRMPM::errorEstimate");

    Ghost::GhostType gnone = Ghost::None;
    constCCVariable<double> refineCell;
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
    new_dw
      ->get(refineCell, d_mpmLabels->MPMRefineCellLabel, 0, patch, gnone, 0);

    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();

    IntVector low    = patch->getCellLowIndex();
    IntVector high   = patch->getCellHighIndex();
    IntVector middle = (low + high) / IntVector(2, 2, 2);

    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      if (refineCell[c] > 0.0 || refineFlag[c] == true) {
        refineFlag[c] = 1;
        refinePatch->set();
      } else if (refineCell[c] < 0.0 && ((int)refineCell[c]) % 100 != 0) {
        refineFlag[c] = 1;
        refinePatch->set();
      } else {
        refineFlag[c] = 0;
      }
    }

#if 0
    // loop over all the geometry objects
    for(int obj=0; obj<(int)d_refineGeomObjs.size(); obj++){
      GeometryPieceP piece = d_refineGeomObjs[obj]->getPiece();
      Vector dx = patch->dCell();
      
      int geom_level =  d_refineGeomObjs[obj]->getInitialData_int("level");
     
      // don't add refinement flags if the current level is greater than
      // the geometry level specification
      if(geom_level!=-1 && level->getIndex()>=geom_level)
        continue;

      for(CellIterator iter=patch->getExtraCellIterator(); !iter.done();iter++){
        IntVector c = *iter;
        Point  lower  = patch->nodePosition(c);
        Vector upperV = lower.asVector() + dx; 
        Point  upper  = upperV.asPoint();
        
        if(piece->inside(upper) && piece->inside(lower))
          refineFlag[c] = true;
        refinePatch->set();
      }
    }
#endif

#if 0
    for(int m = 0; m < d_materialManager->getNumMaterials("MPM"); m++){
      MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM",  m ));
      int dwi = mpm_matl->getDWIndex();
      
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
      
      constParticleVariable<Point> pX;
      constParticleVariable<int> pRefined;
      new_dw->get(pX,       d_mpmLabels->pXLabel,       pset);
      new_dw->get(pRefined, d_mpmLabels->pRefinedLabel, pset);
#if 0
      for(CellIterator iter=patch->getExtraCellIterator(); !iter.done();iter++){
        IntVector c = *iter;
        
        if(level->getIndex()==0 &&(c==IntVector(26,1,0) && step > 48)){
          refineFlag[c] = true;
          refinePatch->set();
        } else{
          refineFlag[c] = false;
        }
      }
#endif

#if 1

      for(ParticleSubset::iterator iter = pset->begin();
          iter!= pset->end();  iter++){
        if(pRefined[*iter]==1){
          IntVector c = level->getCellIndex(pX[*iter]);
          refineFlag[c] = true;
          std::cout << "refineFlag Cell = " << c << std::endl;
          refinePatch->set();
        }
      }
#endif
    }
#endif
  }
}
//______________________________________________________________________
//
void
AMRMPM::refineGrid(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* /*matls*/,
                   DataWarehouse*,
                   DataWarehouse* new_dw)
{
  // just create a particle subset if one doesn't exist
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing AMRMPM::refineGrid");

    CCVariable<short int> cellNAPID;
    new_dw->allocateAndPut(cellNAPID, d_mpmLabels->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

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
        ParticleVariable<Vector> pVelocity, pExternalforce, pDisp, pConcGrad;
        ParticleVariable<Matrix3> pSize;
        ParticleVariable<double> pTempPrev, pColor, pConc, pConcPrev,
          pExtScalFlux;
        ParticleVariable<int> pLoadCurve, pLastLevel, pLocalized, pRefined;
        ParticleVariable<long64> pID;
        ParticleVariable<Matrix3> pdeform, pStress;
        // ParticleVariable<Matrix3> pVelGrad;

        new_dw->allocateAndPut(pX, d_mpmLabels->pXLabel, pset);
        new_dw->allocateAndPut(pMass, d_mpmLabels->pMassLabel, pset);
        new_dw->allocateAndPut(pVelocity, d_mpmLabels->pVelocityLabel, pset);
        new_dw->allocateAndPut(pTemperature,
                               d_mpmLabels->pTemperatureLabel,
                               pset);
        new_dw->allocateAndPut(pTempPrev,
                               d_mpmLabels->pTempPreviousLabel,
                               pset);
        new_dw->allocateAndPut(pExternalforce,
                               d_mpmLabels->pExternalForceLabel,
                               pset);
        new_dw->allocateAndPut(pExtScalFlux,
                               d_mpmLabels->pExternalScalarFluxLabel,
                               pset);
        new_dw->allocateAndPut(pID, d_mpmLabels->pParticleIDLabel, pset);
        new_dw->allocateAndPut(pDisp, d_mpmLabels->pDispLabel, pset);
        new_dw->allocateAndPut(pLastLevel, d_mpmLabels->pLastLevelLabel, pset);
        new_dw->allocateAndPut(pRefined, d_mpmLabels->pRefinedLabel, pset);
        // new_dw->allocateAndPut(pVelGrad,       d_mpmLabels->pVelGradLabel,
        // pset);
        new_dw->allocateAndPut(pVolume, d_mpmLabels->pVolumeLabel, pset);
        if (d_mpmFlags->d_useLoadCurves) {
          new_dw->allocateAndPut(pLoadCurve,
                                 d_mpmLabels->pLoadCurveIDLabel,
                                 pset);
        }
        if (d_mpmFlags->d_withColor) {
          new_dw->allocateAndPut(pColor, d_mpmLabels->pColorLabel, pset);
        }
        if (d_mpmFlags->d_doScalarDiffusion) {
          new_dw->allocateAndPut(pConc,
                                 d_mpmLabels->diffusion->pConcentration,
                                 pset);
          new_dw->allocateAndPut(pConcPrev,
                                 d_mpmLabels->pConcPreviousLabel,
                                 pset);
          new_dw->allocateAndPut(pConcGrad,
                                 d_mpmLabels->pGradConcentrationLabel,
                                 pset);
        }
        new_dw->allocateAndPut(pSize, d_mpmLabels->pSizeLabel, pset);

        // Init deformation gradient
        d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);

        mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                           mpm_matl,
                                                           new_dw);

        // Body force quantities
        ParticleVariable<Vector> pBodyFAcc;
        ParticleVariable<double> pCoriolis;
        new_dw->allocateAndPut(pBodyFAcc,
                               d_mpmLabels->pBodyForceAccLabel,
                               pset);
        new_dw->allocateAndPut(pCoriolis,
                               d_mpmLabels->pCoriolisImportanceLabel,
                               pset);
      }
    }
  }
} // end refine()

//______________________________________________________________________
// Debugging Task that counts the number of particles in the domain.
void
AMRMPM::scheduleCountParticles(const PatchSet* patches, SchedulerP& sched)
{
  printSchedule(patches, cout_doing, "AMRMPM::scheduleCountParticles");
  Task* t =
    scinew Task("AMRMPM::countParticles", this, &AMRMPM::countParticles);
  t->computes(d_mpmLabels->partCountLabel);
  sched->addTask(t, patches, d_materialManager->allMaterials("MPM"));
}
//______________________________________________________________________
//
void
AMRMPM::countParticles(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset*,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw)
{
  long int totalParticles = 0;
  int numMPMMatls         = d_materialManager->getNumMaterials("MPM");

  //  const Level* level = getLevel(patches);
  //  std::cout << "Level " << level->getIndex() << " has " <<
  //  level->numPatches() << " patches" << std::endl;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing AMRMPM::countParticles");

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl =
        static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM", m));
      int dwi = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
      totalParticles += pset->end() - pset->begin();
    }
    //    std::cout << "patch = " << patch->getID()
    //         << ", numParticles = " << totalParticles << std::endl;
  }
  new_dw->put(sumlong_vartype(totalParticles), d_mpmLabels->partCountLabel);
}

//______________________________________________________________________
//

//______________________________________________________________________
// This task colors the particles that are retrieved from the coarse level and
// used on the CFI.  This task mimics interpolateParticlesToGrid_CFI
void
AMRMPM::scheduleDebug_CFI(SchedulerP& sched,
                          const PatchSet* patches,
                          const MaterialSet* matls)
{
  const Level* level = getLevel(patches);

  printSchedule(patches, cout_doing, "AMRMPM::scheduleDebug_CFI");

  Task* t = scinew Task("AMRMPM::debug_CFI", this, &AMRMPM::debug_CFI);

  Ghost::GhostType gn = Ghost::None;
  t->requires(Task::OldDW, d_mpmLabels->pXLabel, gn, 0);
  t->requires(Task::NewDW, d_mpmLabels->pSizeLabel, gn, 0);
  t->requires(Task::OldDW, d_mpmLabels->pDefGradLabel, gn, 0);

  if (level->hasFinerLevel()) {
#define allPatches 0
    Task::MaterialDomainSpec ND = Task::NormalDomain;
    t->requires(Task::NewDW,
                d_mpmLabels->gZOILabel,
                allPatches,
                Task::FineLevel,
                d_oneMaterial,
                Task::NormalDomain,
                gn,
                0);
  }

  t->computes(d_mpmLabels->pColorLabel_preReloc);

  sched->addTask(t, patches, matls);
}
//______________________________________________________________________
void
AMRMPM::debug_CFI(const ProcessorGroup*,
                  const PatchSubset* patches,
                  const MaterialSubset*,
                  DataWarehouse* old_dw,
                  DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  for (int cp = 0; cp < patches->size(); cp++) {
    const Patch* patch = patches->get(cp);

    printTask(patches, patch, cout_doing, "Doing AMRMPM::debug_CFI");

    //__________________________________
    //  Write p.color all levels all patches
    ParticleSubset* pset = 0;
    int dwi              = 0;
    pset                 = old_dw->getParticleSubset(dwi, patch);

    constParticleVariable<Point> pX;
    constParticleVariable<Matrix3> pSize;
    constParticleVariable<Matrix3> pDefGrad;
    ParticleVariable<double> pColor;

    old_dw->get(pX, d_mpmLabels->pXLabel, pset);
    new_dw->get(pSize, d_mpmLabels->pSizeLabel, pset);
    old_dw->get(pDefGrad, d_mpmLabels->pDefGradLabel, pset);
    new_dw->allocateAndPut(pColor, d_mpmLabels->pColorLabel_preReloc, pset);

    auto interpolatorCoarse = d_mpmFlags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolatorCoarse->size());
    std::vector<double> S(interpolatorCoarse->size());

    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
         iter++) {
      particleIndex idx = *iter;
      pColor[idx]       = 0;

      interpolatorCoarse->findCellAndWeights(pX[idx],
                                             ni,
                                             S,
                                             pSize[idx],
                                             pDefGrad[idx]);

      for (int k = 0; k < (int)ni.size(); k++) {
        pColor[idx] += S[k];
      }
    }

    //__________________________________
    //  Mark the particles that are accessed at the CFI.
    if (level->hasFinerLevel()) {

      // find the fine patches over the coarse patch.  Add a single layer of
      // cells
      // so you will get the required patches when coarse patch and fine patch
      // boundaries coincide.
      Level::selectType finePatches;
      patch->getOtherLevelPatches(1, finePatches, 1);

      const Level* fineLevel = level->getFinerLevel().get_rep();
      IntVector refineRatio(fineLevel->getRefinementRatio());

      for (size_t fp = 0; fp < finePatches.size(); fp++) {
        const Patch* finePatch = finePatches[fp];

        // Determine extents for coarser level particle data
        // Linear Interpolation:  1 layer of coarse level cells
        // Gimp Interpolation:    2 layers
        /*`==========TESTING==========*/
        IntVector nLayers(d_nPaddingCells_Coarse,
                          d_nPaddingCells_Coarse,
                          d_nPaddingCells_Coarse);
        IntVector nPaddingCells = nLayers * (refineRatio);
        // std::cout << " nPaddingCells " << nPaddingCells << "nLayers " <<
        // nLayers
        // << std::endl;
        /*===========TESTING==========`*/

        int nGhostCells           = 0;
        bool returnExclusiveRange = false;
        IntVector cl_tmp, ch_tmp, fl, fh;

        getCoarseLevelRange(finePatch,
                            level,
                            cl_tmp,
                            ch_tmp,
                            fl,
                            fh,
                            nPaddingCells,
                            nGhostCells,
                            returnExclusiveRange);

        cl_tmp -= finePatch->neighborsLow() *
                  nLayers; //  expand cl_tmp when a neighor patch exists.
        //  This patch owns the low nodes.  You need particles
        //  from the neighbor patch.

        // coarseLow and coarseHigh cannot lie outside of the coarse patch
        IntVector cl = Max(cl_tmp, patch->getCellLowIndex());
        IntVector ch = Min(ch_tmp, patch->getCellHighIndex());

        ParticleSubset* pset2 = 0;
        pset2 =
          old_dw->getParticleSubset(dwi, cl, ch, patch, d_mpmLabels->pXLabel);

        constParticleVariable<Point> pX_CFI;
        constNCVariable<Stencil7> zoi;
        old_dw->get(pX_CFI, d_mpmLabels->pXLabel, pset2);
        new_dw->get(zoi, d_mpmLabels->gZOILabel, 0, finePatch, Ghost::None, 0);

        auto interpolatorFine = d_mpmFlags->d_interpolator->clone(finePatch);

        for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
             iter++) {
          particleIndex idx = *iter;

          for (ParticleSubset::iterator iter2 = pset2->begin();
               iter2 != pset2->end();
               iter2++) {
            particleIndex idx2 = *iter2;

            if (pX[idx] == pX_CFI[idx2]) {
              pColor[idx] = 0;
              std::vector<IntVector> ni;
              std::vector<double> S;
              interpolatorFine->findCellAndWeights_CFI(pX[idx], ni, S, zoi);
              for (int k = 0; k < (int)ni.size(); k++) {
                pColor[idx] += S[k];
              }
            }
          } // pset2 loop
        }   // pset loop
      }     // loop over fine patches
    }       //// hasFinerLevel
  }         // End loop over coarse patches
}

//______________________________________________________________________
//  Returns the patches that have coarse fine interfaces
void
AMRMPM::coarseLevelCFI_Patches(const PatchSubset* coarsePatches,
                               Level::selectType& CFI_patches)
{
  for (int p = 0; p < coarsePatches->size(); p++) {
    const Patch* coarsePatch = coarsePatches->get(p);

    Level::selectType finePatches;
    coarsePatch->getFineLevelPatches(finePatches);
    // loop over all the coarse level patches

    for (size_t fp = 0; fp < finePatches.size(); fp++) {
      const Patch* finePatch = finePatches[fp];

      if (finePatch->hasCoarseFaces()) {
        CFI_patches.push_back(coarsePatch);
      }
    }
  }
}

#if 0
//______________________________________________________________________
//
void AMRMPM::scheduleInterpolateToParticlesAndUpdate_CFI(SchedulerP& sched,
                                                         const PatchSet* patches,
                                                         const MaterialSet* matls)

{
  const Level* level = getLevel(patches);
  
  if(level->hasFinerLevel()){

    printSchedule(patches,cout_doing,"AMRMPM::scheduleInterpolateToParticlesAndUpdate_CFI");

    Task* t=scinew Task("AMRMPM::interpolateToParticlesAndUpdate_CFI",
                        this, &AMRMPM::interpolateToParticlesAndUpdate_CFI);

    Ghost::GhostType  gn  = Ghost::None;
    Task::MaterialDomainSpec  ND  = Task::NormalDomain;
#define allPatches 0
#define allMatls 0
    t->requires(Task::OldDW, d_mpmLabels->delTLabel );
    
    t->requires(Task::OldDW, d_mpmLabels->pXLabel, gn);
    t->requires(Task::NewDW, d_mpmLabels->gVelocityStarLabel, allPatches, Task::FineLevel,allMatls,   ND, gn,0);
    t->requires(Task::NewDW, d_mpmLabels->gAccelerationLabel, allPatches, Task::FineLevel,allMatls,   ND, gn,0);
    t->requires(Task::NewDW, d_mpmLabels->gZOILabel,          allPatches, Task::FineLevel,d_oneMaterial, ND, gn,0);
    
    t->modifies(d_mpmLabels->pXLabel_preReloc);
    t->modifies(d_mpmLabels->pDispLabel_preReloc);
    t->modifies(d_mpmLabels->pVelocityLabel_preReloc);

    sched->addTask(t, patches, matls);
  }
}
#endif

#if 0
//______________________________________________________________________
//
void AMRMPM::interpolateToParticlesAndUpdate_CFI(const ProcessorGroup*,
                                                 const PatchSubset* coarsePatches,
                                                 const MaterialSubset* ,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw)
{
  const Level* coarseLevel = getLevel(coarsePatches);
  const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
  
  delt_vartype delT;
  old_dw->get(delT, d_mpmLabels->delTLabel, coarseLevel );
  
  double move_particles=1.;
  if(!d_mpmFlags->d_doGridReset){
    move_particles=0.;
  }
  
  //__________________________________
  //Loop over the coarse level patches
  for(int p=0;p<coarsePatches->size();p++){
    const Patch* coarsePatch = coarsePatches->get(p);
    printTask(coarsePatches,coarsePatch,cout_doing,"AMRMPM::interpolateToParticlesAndUpdate_CFI");

    int numMatls = d_materialManager->getNumMaterials("MPM");
        
    Level::selectType finePatches;
    coarsePatch->getFineLevelPatches(finePatches);
    
    
    //__________________________________
    //  Fine patch loop
    for(int i=0;i<finePatches.size();i++){
      const Patch* finePatch = finePatches[i]; 
      
      if(finePatch->hasCoarseFaces()){

        auto interpolator = d_mpmFlags->d_interpolator->clone(finePatch);
        
        constNCVariable<Stencil7> zoi_fine;
        new_dw->get(zoi_fine, d_mpmLabels->gZOILabel, 0, finePatch, Ghost::None, 0 );

        for(int m = 0; m < numMatls; m++){
          MPMMaterial* mpm_matl = static_cast<MPMMaterial*>(d_materialManager->getMaterial("MPM",  m ));
          int dwi = mpm_matl->getDWIndex();

          // get fine level grid data
          constNCVariable<double> gMass_fine;
          constNCVariable<Vector> gVelocity_star_fine;
          constNCVariable<Vector> gAcceleration_fine;

          // use getRegion() instead of get().  They should be equivalent but 
          // get() throws assert on parallel runs.
          IntVector fl = finePatch->getNodeLowIndex();
          IntVector fh = finePatch->getNodeHighIndex();
          new_dw->getRegion(gVelocity_star_fine,  d_mpmLabels->gVelocityStarLabel, dwi, fineLevel,fl, fh);   
          new_dw->getRegion(gAcceleration_fine,   d_mpmLabels->gAccelerationLabel, dwi, fineLevel,fl, fh); 
            
          
          // get coarse level particle data
          ParticleVariable<Point>  pXnew_coarse;
          ParticleVariable<Vector> pDispnew_coarse;
          ParticleVariable<Vector> pVelocitynew_coarse;
          constParticleVariable<Point>  pXold_coarse;
          
          ParticleSubset* pset=0;
          
/*`==========TESTING==========*/
          // get the particles for the entire coarse patch
          // Ideally you only need the particle subset in the cells
          // that surround the fine patch.  Currently, getModifiable doesn't
          // allow you to get a pset with a high/low index that does not match
          // the patch low high index  
          pset = old_dw->getParticleSubset(dwi, coarsePatch);
          //std::cout << *pset << std::endl; 
/*===========TESTING==========`*/
          old_dw->get(pXold_coarse,                  d_mpmLabels->pXLabel,                 pset);
          new_dw->getModifiable(pXnew_coarse,        d_mpmLabels->pXLabel_preReloc,        pset);
          new_dw->getModifiable(pDispnew_coarse,     d_mpmLabels->pDispLabel_preReloc,     pset);
          new_dw->getModifiable(pVelocitynew_coarse, d_mpmLabels->pVelocityLabel_preReloc, pset);


          for (ParticleSubset::iterator iter = pset->begin();
               iter != pset->end(); iter++){
            particleIndex idx = *iter;

            // Get the node indices that surround the fine patch cell
            std::vector<IntVector> ni;
            std::vector<double> S;
            
            interpolator->findCellAndWeights_CFI(pXold_coarse[idx],ni,S,zoi_fine);

            Vector acc(0.0, 0.0, 0.0); 
            Vector vel(0.0, 0.0, 0.0);

            // Add each nodes contribution to the particle's velocity & acceleration 
            IntVector fineNode;
            for(int k = 0; k < (int)ni.size(); k++) {
              
              fineNode = ni[k];
            
              vel  += gVelocity_star_fine[fineNode] * S[k];
              acc  += gAcceleration_fine[fineNode]  * S[k];
/*`==========TESTING==========*/
#ifdef DEBUG_ACC 
              Vector diff = acc - d_acc_ans;
              if( abs(acc.length() - d_acc_ans.length() > d_acc_tol ) ) {
                const Level* fineLevel = coarseLevel->getFinerLevel().get_rep();
                std::cout << "    L-"<< fineLevel->getIndex() << " node: "<< fineNode << " gAcceleration: " << gAcceleration_fine[fineNode] 
                     << "  diff " << diff << std::endl;
              }
#endif 
/*===========TESTING==========`*/
            }
            
//            std::cout << " pVelocitynew_coarse  "<< idx << " "  << pVelocitynew_coarse[idx] << " p.x " << pXnew_coarse[idx] ;
            
            // Update the particle's position and velocity
            pXnew_coarse[idx]         += vel*delT*move_particles;  
            pDispnew_coarse[idx]      += vel*delT;                 
            pVelocitynew_coarse[idx]  += acc*delT; 
            
          } // End of particle loop
        } // End loop over materials 
      
      }  // if has coarse face
    }  // End loop over fine patches 
  }  // End loop over patches
}
#endif
