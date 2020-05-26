/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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
#include <CCA/Components/MPM/ImpMPMFlags.h> 
#include <CCA/Components/MPM/MPMUtils.h> 
#include <Core/Math/Matrix3.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ImplicitCM.h>
#include <CCA/Components/MPM/PhysicalBC/MPMPhysicalBCFactory.h>
#include <CCA/Components/MPM/PhysicalBC/PressureBC.h>
#include <CCA/Components/MPM/PhysicalBC/HeatFluxBC.h>
#include <CCA/Components/MPM/HeatConduction/ImplicitHeatConduction.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContact.h>
#include <CCA/Components/MPM/ThermalContact/ThermalContactFactory.h>
#include <CCA/Components/MPM/MPMBoundCond.h>
#include <CCA/Ports/DataWarehouse.h>
#include <CCA/Ports/Scheduler.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/AMR.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Variables/SoleVariable.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/BoundaryConditions/BoundCond.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Parallel/ProcessorGroup.h>
#include <Core/Parallel/Parallel.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/Point.h>
#include <Core/Math/MinMax.h>
#include <CCA/Ports/LoadBalancer.h>
#include <CCA/Ports/Output.h>
#include <Core/Util/DebugStream.h>
#include <CCA/Components/MPM/PetscSolver.h>
#include <CCA/Components/MPM/SimpleSolver.h>
#include <CCA/Components/Regridder/PerPatchVars.h>
#include <Core/Grid/BoundaryConditions/BCDataArray.h>
#include <Core/Math/FastMatrix.h>
#include <set>
#include <map>
#include <numeric>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>
#include <cmath>


using namespace Uintah;

using NCdoubleArray  = std::vector<NCVariable<double> >;
using NCVectorArray  = std::vector<NCVariable<Vector> >;
using NCMatrix3Array = std::vector<NCVariable<Matrix3> >;
using CellParticleTempMap = std::multimap<IntVector, ParticleTempShape>;
using CellParticleTempPair = std::pair<IntVector, ParticleTempShape>;
using CellParticleTempMapArray = std::vector<CellParticleTempMap>;


//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "IMPM:+,IMPM_debug:+".....
//  bash     : export SCI_DEBUG="IMPM:+,IMPM_debug:+" 
//  default is OFF
static DebugStream cout_doing("IMPM", false);
static DebugStream cout_dbg("IMPM_Debug", false);


ImpMPM::ImpMPM(const ProcessorGroup* myworld) :
  MPMCommon(myworld), UintahParallelComponent(myworld)
{
  lb = scinew MPMLabel();
  flags = scinew ImpMPMFlags(myworld);
  d_nextOutputTime = 0.;
  d_SMALL_NUM_MPM = 1e-200;
  d_rigid_body = false;
  d_numIterations = 0;

  heatConductionModel = nullptr;
  thermalContactModel = nullptr;
  d_perproc_patches = nullptr;
  d_switchCriteria = nullptr;
  d_defGradComputer = nullptr;

  one_matl = scinew MaterialSubset();
  one_matl->add(0);
  one_matl->addReference();
  NGP = 1;
  NGN = 1;
  d_loadCurveIndex=0;
}

bool 
ImpMPM::restartableTimesteps()
{
  return true;
}

ImpMPM::~ImpMPM()
{
  delete lb;
  delete flags;

  if(d_perproc_patches && d_perproc_patches->removeReference()) { 
    delete d_perproc_patches;
    std::cout << "Freeing patches!!\n";
  }

  if(one_matl->removeReference()) {
    delete one_matl;
  }

  delete d_solver;
  delete heatConductionModel;
  delete thermalContactModel;
  MPMPhysicalBCFactory::clean();

  if(d_switchCriteria) {
    delete d_switchCriteria;
  }

  if (d_defGradComputer) {
    delete d_defGradComputer;
  }

}

void 
ImpMPM::problemSetup(const ProblemSpecP& prob_spec, 
                     const ProblemSpecP& restart_prob_spec,
                     GridP& grid,
                     SimulationStateP& sharedState)
{
  cout_doing << " Doing ImpMPM::problemSetup " << endl;
  d_sharedState = sharedState;
  dynamic_cast<Scheduler*>(getPort("scheduler"))->setPositionVar(lb->pXLabel);
  
  Output* dataArchiver = dynamic_cast<Output*>(getPort("output"));
  if(!dataArchiver){
    throw InternalError("ImpMPM:couldn't get output port", __FILE__, __LINE__);
  }

  ProblemSpecP mpm_ps = 0;
  ProblemSpecP restart_mat_ps = 0;

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

  std::string integrator_type;
  if (mpm_soln_ps) {

    // Read all MPM flags (look in MPMFlags.cc)
    flags->readMPMFlags(restart_mat_ps, dataArchiver);
     
    if (flags->d_integratorType != "implicit") {
      throw ProblemSetupException("Can't use explicit integration with -impm", __FILE__, __LINE__);
    }

    // convert text representation of face into FaceType
    for (auto face : flags->d_bndyFaceTxtList) {
      auto faceType = Patch::invalidFace;
      for (auto ft = Patch::startFace; ft <= Patch::endFace; ft = Patch::nextFace(ft)) {
        if (Patch::getFaceName(ft) == face) {
          faceType = ft;
          break;
        }
      }
      if (faceType != Patch::invalidFace) {
        d_bndy_traction_faces.push_back(faceType);
      } else {
        std::cerr << "warning: ignoring unknown face '" << face << "'\n";
      }
    }
  }
   
  // read in AMR flags from the main ups file
  ProblemSpecP amr_ps = prob_spec->findBlock("AMR");
  if (amr_ps) {
    ProblemSpecP mpm_amr_ps = amr_ps->findBlock("MPM");
    mpm_amr_ps->getWithDefault("min_grid_level", flags->d_minGridLevel, 0);
    mpm_amr_ps->getWithDefault("max_grid_level", flags->d_maxGridLevel, 1000);
  }

  if (flags->d_8or27==8) {
    NGP=1;
    NGN=1;
  } else if (flags->d_8or27==27 || flags->d_8or27==64) {
    NGP=2;
    NGN=2;
  }

  //Search for the MaterialProperties block and then get the MPM section
  ProblemSpecP mat_ps =  
     restart_mat_ps->findBlockWithOutAttribute("MaterialProperties");
  ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");

  d_con_type = "null";
  ProblemSpecP contact_ps = mpm_mat_ps->findBlock("contact");
  if (contact_ps){
    contact_ps->getWithDefault("type", d_con_type, "null");
  }
  d_rigid_body = false;
  if (d_con_type == "rigid") {
    d_rigid_body = true;
    Vector defaultDir(1,1,1);
    contact_ps->getWithDefault("direction", d_contact_dirs, defaultDir);
    contact_ps->getWithDefault("stop_time", d_stop_time, 
                                            std::numeric_limits<double>::max());
    contact_ps->getWithDefault("velocity_after_stop", d_vel_after_stop, 
                                                      Vector(0,0,0));
  }

  d_sharedState->setParticleGhostLayer(Ghost::AroundNodes, 1);

  MPMPhysicalBCFactory::create(restart_mat_ps, grid, flags);
  if (MPMPhysicalBCFactory::mpmPhysicalBCs.size() == 0 && flags->d_useLoadCurves) {
    throw ProblemSetupException("No load curve in ups, d_useLoadCurve==true?", __FILE__, __LINE__);
  }

  materialProblemSetup(restart_mat_ps, grid, d_sharedState,flags);
   
  if (flags->d_solverType == "petsc") {
    d_solver = scinew MPMPetscSolver();
    d_solver->initialize();
  } else {
    d_solver = scinew SimpleSolver();
    d_solver->initialize();
  }

  d_defGradComputer = scinew DeformationGradientComputer(flags, d_sharedState);

  // setup sub scheduler
  Scheduler* sched = dynamic_cast<Scheduler*>(getPort("scheduler"));
  sched->setRestartable(true);

  d_recompileSubsched = true;
  d_subsched = sched->createSubScheduler();
  std::cout << "Creating sub scheduler: " << d_subsched << "\n";
  d_subsched->initialize(3, 1);
  d_subsched->clearMappings();
  d_subsched->mapDataWarehouse(Task::ParentOldDW, 0);
  d_subsched->mapDataWarehouse(Task::ParentNewDW, 1);
  d_subsched->mapDataWarehouse(Task::OldDW, 2);
  d_subsched->mapDataWarehouse(Task::NewDW, 3);
   
  heatConductionModel = scinew ImplicitHeatConduction(sharedState, lb, flags);
  heatConductionModel->problemSetup(flags->d_solverType);

  thermalContactModel =
     ThermalContactFactory::create(restart_mat_ps, sharedState, lb,flags);

  d_switchCriteria = dynamic_cast<SwitchingCriteria*>(getPort("switch_criteria"));
  if (d_switchCriteria) {
    d_switchCriteria->problemSetup(restart_mat_ps,restart_prob_spec,d_sharedState);
  }

  // Pull out from Time section
  d_initialDt = 10000.0;
  ProblemSpecP time_ps = restart_mat_ps->findBlock("Time");
  if ( time_ps ) {
    time_ps->get("delt_init", d_initialDt);
  }
}

void 
ImpMPM::outputProblemSpec(ProblemSpecP& root_ps)
{
  ProblemSpecP root = root_ps->getRootNode();

  ProblemSpecP flags_ps = root->appendChild("MPM");
  flags->outputProblemSpec(flags_ps);

  ProblemSpecP mat_ps = nullptr;
  mat_ps = root->findBlockWithOutAttribute("MaterialProperties");

  if( mat_ps == nullptr ) {
    mat_ps = root->appendChild("MaterialProperties");
  }
    
  ProblemSpecP mpm_ps = mat_ps->appendChild("MPM");
  for (int i = 0; i < d_sharedState->getNumMPMMatls(); i++) {
    MPMMaterial* mat = d_sharedState->getMPMMaterial(i);
    ProblemSpecP cm_ps = mat->outputProblemSpec(mpm_ps);
  }

  ProblemSpecP contact_ps = mpm_ps->appendChild("contact");
  contact_ps->appendElement("type",d_con_type);
  contact_ps->appendElement("direction",d_contact_dirs);
  contact_ps->appendElement("stop_time",d_stop_time);
  contact_ps->appendElement("velocity_after_stop",d_vel_after_stop);

  ProblemSpecP physical_bc_ps = root->appendChild("PhysicalBC");
  ProblemSpecP mpm_ph_bc_ps = physical_bc_ps->appendChild("MPM");
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    particleBC->outputProblemSpec(mpm_ph_bc_ps);
  }
}

void 
ImpMPM::scheduleInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels())) {
    return;
  }

  LoadBalancer* loadbal = sched->getLoadBalancer();
  d_perproc_patches = loadbal->getPerProcessorPatchSet(level);

  Task* t = scinew Task("ImpMPM::actuallyInitialize",
                        this, &ImpMPM::actuallyInitialize);

  const PatchSet* patches = level->eachPatch();

  t->computes(lb->partCountLabel);
  t->computes(lb->pXLabel);
  t->computes(lb->pDispLabel);
  t->computes(lb->pMassLabel);
  t->computes(lb->pVolumeLabel);
  t->computes(lb->pFiberDirLabel);
  t->computes(lb->pVelocityLabel);
  t->computes(lb->pAccelerationLabel);
  t->computes(lb->pExternalForceLabel);
  t->computes(lb->pTemperatureLabel);
  t->computes(lb->pTempPreviousLabel);
  t->computes(lb->pSizeLabel);
  t->computes(lb->pParticleIDLabel);
  //t->computes(lb->pDefGradLabel);
  //t->computes(lb->pVelGradLabel);
  //t->computes(lb->pDispGradLabel);
  t->computes(lb->pStressLabel);
  t->computes(lb->pRefinedLabel);
  t->computes(lb->pCellNAPIDLabel);

  t->computes(lb->pCoriolisImportanceLabel);
  t->computes(lb->pBodyForceAccLabel);

  t->computes(lb->pExternalHeatFluxLabel);
  t->computes(lb->heatRate_CCLabel);

  t->computes(d_sharedState->get_delt_label(), level.get_rep());

  int numMPM = d_sharedState->getNumMPMMatls();
  for(int m = 0; m < numMPM; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

    d_defGradComputer->addInitialComputesAndRequires(t, mpm_matl, patches);

    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addInitialComputesAndRequires(t, mpm_matl, patches);
  }

  if (flags->d_artificialViscosity) {
    t->computes(lb->p_qLabel);        //  only used for imp -> exp transition
  }

  if (flags->d_useLoadCurves) {
    t->computes(lb->pLoadCurveIDLabel);
  }

  // For friction contact
  t->computes(lb->pSurfLabel);

  if (!flags->d_doGridReset) {
    t->computes(lb->gDisplacementLabel);
  }

  t->computes(lb->NC_CCweightLabel, one_matl);
  if (!flags->d_tempSolve) {
    t->computes(lb->gTemperatureLabel, one_matl);
  }

  sched->addTask(t, d_perproc_patches, d_sharedState->allMPMMaterials());

  t = scinew Task("ImpMPM::printParticleCount",
                  this, &ImpMPM::printParticleCount);
  t->requires(Task::NewDW, lb->partCountLabel);
  t->setType(Task::OncePerProc);
  sched->addTask(t, d_perproc_patches, d_sharedState->allMPMMaterials());

  if (flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    scheduleInitializeHeatFluxBCs(level, sched);
    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
  }

  if (d_switchCriteria) {
    d_switchCriteria->scheduleInitialize(level, sched);
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
    new_dw->allocateAndPut(cellNAPID, lb->pCellNAPIDLabel, 0, patch);
    cellNAPID.initialize(0);

    CCVariable<double> heatFlux;
    new_dw->allocateAndPut(heatFlux,lb->heatRate_CCLabel,0,patch);
    heatFlux.initialize(1.0);

    for (int m = 0; m < matls->size(); m++) {
      int matl = matls->get(m);
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( matl );
      particleIndex numParticles = 
        mpm_matl->createParticles(cellNAPID, patch, new_dw);
      totalParticles += numParticles;

      d_defGradComputer->initializeGradient(patch, mpm_matl, new_dw);
      mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                         mpm_matl, new_dw);
      if (!flags->d_doGridReset) {
        int matID = mpm_matl->getDWIndex();
        NCVariable<Vector> gDisplacement;
        new_dw->allocateAndPut(gDisplacement, lb->gDisplacementLabel, matID, patch);        
        gDisplacement.initialize(Vector(0.));
      }
    }
    
    std::string interp_type = flags->d_interpolatorType;
    if ((interp_type == "gimp" || 
         interp_type == "3rdorderBS" || 
         interp_type == "cpdi" || 
         interp_type == "cpti" || 
         interp_type == "cpgimp")) {
      proc0cout << "__________________________________\n"
                << "WARNING: Use of GIMP/3rdorderBS/cpdi/cpgimp with Implicit "
                   "MPM is untested and may not work at this time.\n\n";
    }
    
    IntVector num_extra_cells = patch->getExtraCells();
    IntVector periodic = patch->getLevel()->getPeriodicBoundaries();
    
    if (interp_type=="linear" && num_extra_cells != IntVector(0,0,0)) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>linear</interpolator> \n"
          << " you should also use <extraCells>[0,0,0]</extraCells> \n";
      throw ProblemSetupException(msg.str(),__FILE__, __LINE__);
    } else if ((interp_type == "gimp" || 
                interp_type =="3rdorderBS" || 
                interp_type =="cpdi" || 
                interp_type =="cpti" || 
                interp_type =="cpgimp") && 
                (num_extra_cells + periodic) != IntVector(1,1,1)) {
      std::ostringstream msg;
      msg << "\n ERROR: When using <interpolator>gimp</interpolator> \n"
          << " or <interpolator>3rdorderBS</interpolator> \n"
          << " or <interpolator>cpdi</interpolator> \n"
          << " or <interpolator>cpti</interpolator> \n"
          << " or <interpolator>cpgimp</interpolator> \n"
          << " you must also use extraCells and/or periodicBCs such that\n"
          << " the sum of the two is [1,1,1].\n";
      throw ProblemSetupException(msg.str(),__FILE__, __LINE__);
    }

    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight, lb->NC_CCweightLabel,    0, patch);
    NC_CCweight.initialize(0.125);
    for (auto face = Patch::startFace; face <= Patch::endFace;
              face = Patch::nextFace(face)) {
      int mat_id = 0;
      if (patch->haveBC(face,mat_id,"symmetry","Symmetric")) {
        for(auto iter = patch->getFaceIterator(face, Patch::FaceNodes);
                 !iter.done(); iter++) {
          NC_CCweight[*iter] = 2.0*NC_CCweight[*iter];
        }
      }
    }
    if (!flags->d_tempSolve) {
      NCVariable<double> gTemperature;
      new_dw->allocateAndPut(gTemperature, lb->gTemperatureLabel,    0, patch);
      gTemperature.initialize(0.);
    }
  }
  new_dw->put(sumlong_vartype(totalParticles), lb->partCountLabel);
}

void 
ImpMPM::printParticleCount(const ProcessorGroup* pg,
                           const PatchSubset*,
                           const MaterialSubset*,
                           DataWarehouse*,
                           DataWarehouse* new_dw)
{
  if (pg->myrank() == 0) {
    sumlong_vartype pcount;
    new_dw->get(pcount, lb->partCountLabel);
    std::cerr << "Created " << (long) pcount << " total particles\n";
  }
}

void 
ImpMPM::switchInitialize(const LevelP& level, SchedulerP& sched)
{
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels())) {
    return;
  }

  if (flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    if(UintahParallelComponent::d_myworld->myrank() == 0){
      std::cout << " \n--------------------------------------------------------------"<< endl;
      std::cout << " ImpMPM: the heat flux BC cannot be applied on the timestep" << endl; 
      std::cout << " immediately after a component switch.  The computes/requires " << endl;
      std::cout << " cannot be met and one pseudo timestep must take place" << endl;
      std::cout << " ---------------------------------------------------------------\n"<< endl;
    }
    scheduleInitializeHeatFluxBCs(level, sched);

    // Schedule the initialization of pressure BCs per particle
    scheduleInitializePressureBCs(level, sched);
  }
}

void 
ImpMPM::scheduleInitializeHeatFluxBCs(const LevelP& level,
                                      SchedulerP& sched)
{
  MaterialSubset* loadCurveIndex = scinew MaterialSubset();
  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {
      loadCurveIndex->add(nofHeatFluxBCs++);
    }
  }

  if (nofHeatFluxBCs > 0) {

    // Create a task that calculates the total number of particles
    // associated with each load curve.  
    Task* t = scinew Task("ImpMPM::countMaterialPointsPerLoadCurve",
                          this, &ImpMPM::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, lb->pLoadCurveIDLabel, Ghost::None);
    t->computes(lb->materialPointsPerLoadCurveLabel, loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, level->eachPatch(), d_sharedState->allMPMMaterials());

    // Create a task that calculates the heatflux to be associated with
    // each particle based on the HeatFluxBCs
    t = scinew Task("ImpMPM::initializeHeatFluxBC",
                    this, &ImpMPM::initializeHeatFluxBC);
    t->requires(Task::NewDW, lb->pXLabel, Ghost::None);
    t->requires(Task::NewDW, lb->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW, lb->materialPointsPerLoadCurveLabel, loadCurveIndex, 
                Task::OutOfDomain, Ghost::None);
    t->modifies(lb->pExternalHeatFluxLabel);
    sched->addTask(t, level->eachPatch(), d_sharedState->allMPMMaterials());

  } else {
    delete loadCurveIndex;
  }
}

void 
ImpMPM::scheduleInitializePressureBCs(const LevelP& level,
                                      SchedulerP& sched)
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
                          this, &ImpMPM::countMaterialPointsPerLoadCurve);
    t->requires(Task::NewDW, lb->pLoadCurveIDLabel, Ghost::None);
    t->computes(lb->materialPointsPerLoadCurveLabel, d_loadCurveIndex,
                Task::OutOfDomain);
    sched->addTask(t, level->eachPatch(), d_sharedState->allMPMMaterials());

    // Create a task that calculates the force to be associated with
    // each particle based on the pressure BCs
    t = scinew Task("ImpMPM::initializePressureBC",
                    this, &ImpMPM::initializePressureBC);
    t->requires(Task::NewDW, lb->pXLabel, Ghost::None);
    t->requires(Task::NewDW, lb->pLoadCurveIDLabel, Ghost::None);
    t->requires(Task::NewDW, lb->materialPointsPerLoadCurveLabel, d_loadCurveIndex, 
                Task::OutOfDomain, Ghost::None);
    t->modifies(lb->pExternalForceLabel);
    sched->addTask(t, level->eachPatch(), d_sharedState->allMPMMaterials());
  }

  if(d_loadCurveIndex->removeReference()) {
    delete d_loadCurveIndex;
  }
}

void 
ImpMPM::countMaterialPointsPerLoadCurve(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset*,
                                        DataWarehouse* ,
                                        DataWarehouse* new_dw)
{
  // Find the number of pressure BCs in the problem
  int nofPressureBCs = 0;
  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {
      nofHeatFluxBCs++;
      //std::cout << "nofHeatFluxBCs = " << nofHeatFluxBCs << endl;

      for(int p=0;p<patches->size();p++){
        const Patch* patch = patches->get(p);
        int numMPMMatls=d_sharedState->getNumMPMMatls();
        int numPts = 0;
        for(int m = 0; m < numMPMMatls; m++){
          MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);

          for(auto idx : *pset) {
            if (pLoadCurveID[idx] == (nofHeatFluxBCs)) {
              ++numPts;
            }
          }
        } // matl loop
        //std::cout << "numPts found = " << numPts << endl;
        new_dw->put(sumlong_vartype(numPts), 
                    lb->materialPointsPerLoadCurveLabel, 0, nofHeatFluxBCs-1);
      }  // patch loop
    } else if (bcType == "Pressure") {
      nofPressureBCs++;

      for(int p=0;p<patches->size();p++){
        const Patch* patch = patches->get(p);
        int numMPMMatls=d_sharedState->getNumMPMMatls();
        int numPts = 0;
        for(int m = 0; m < numMPMMatls; m++){
          MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<int> pLoadCurveID;
          new_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);

          for(auto idx : *pset) {
            if (pLoadCurveID[idx] == (nofPressureBCs)) ++numPts;
          }
        } // matl loop
        new_dw->put(sumlong_vartype(numPts),
                    lb->materialPointsPerLoadCurveLabel, 0, nofPressureBCs-1);
      }  // patch loop
    }
  }
}

void 
ImpMPM::initializeHeatFluxBC(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* ,
                             DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;

  // Calculate the heat flux at each particle
  int nofHeatFluxBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "HeatFlux") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart, lb->materialPointsPerLoadCurveLabel,
                  0, nofHeatFluxBCs++);

      double fluxPerPart = 0.;
      HeatFluxBC* phf = nullptr;
      if (bcType == "HeatFlux") {
        phf = dynamic_cast<HeatFluxBC*>(particleBC.get());
        phf->numMaterialPoints(numPart);
        fluxPerPart = phf->fluxPerParticle(time);
        //std::cout << "numPart = " << numPart << endl;
        //std::cout << "fluxPerPart = " << fluxPerPart << endl;
      }
      
      // Loop through the patches and calculate the force vector
      // at each particle
      for(int p=0;p<patches->size();p++){
        const Patch* patch = patches->get(p);
        int numMPMMatls=d_sharedState->getNumMPMMatls();
        for(int m = 0; m < numMPMMatls; m++){
          MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
          int matID = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
          constParticleVariable<Point>  pX;
          constParticleVariable<int>    pLoadCurveID;
          ParticleVariable<double>      pExternalHeatFlux;
          new_dw->get(pX,           lb->pXLabel,           pset);
          new_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);
          new_dw->getModifiable(pExternalHeatFlux, lb->pExternalHeatFluxLabel, 
                                pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == nofHeatFluxBCs) {
              if (bcType == "HeatFlux") {
                pExternalHeatFlux[idx] = phf->getFlux(pX[idx], fluxPerPart);
              }
              // std::cout << "pExternalHeatFlux[idx] = " << pExternalHeatFlux[idx] << endl;
            }
          }
        } // matl loop
      }  // patch loop
    }
  }
}

void 
ImpMPM::initializePressureBC(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse* ,
                             DataWarehouse* new_dw)
{
  // Get the current time
  double time = 0.0;

  if (cout_dbg.active()) {
    cout_dbg << "Current Time (Initialize Pressure BC) = " << time << endl;
  }

  // Calculate the force vector at each particle
  int nofPressureBCs = 0;
  for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
    auto bcType = particleBC->getType();
    if (bcType == "Pressure") {

      // Get the material points per load curve
      sumlong_vartype numPart = 0;
      new_dw->get(numPart, lb->materialPointsPerLoadCurveLabel,
                  0, nofPressureBCs++);

      // Save the material points per load curve in the PressureBC object
      PressureBC* pbc = dynamic_cast<PressureBC*>(particleBC.get());
      pbc->numMaterialPoints(numPart);

      if (cout_dbg.active()) {
        cout_dbg << "    Load Curve = " << nofPressureBCs << " Num Particles = " << numPart << endl;
      }

      // Calculate the force per particle at t = 0.0
      double forcePerPart = pbc->forcePerParticle(time);

      // Loop through the patches and calculate the force vector
      // at each particle
      for(int p=0;p<patches->size();p++){
        const Patch* patch = patches->get(p);
        int numMPMMatls=d_sharedState->getNumMPMMatls();
        for(int m = 0; m < numMPMMatls; m++){
          MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
          int dwi = mpm_matl->getDWIndex();

          ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
          constParticleVariable<Point>   pX;
          constParticleVariable<int>     pLoadCurveID;
          constParticleVariable<Matrix3> pDefGrad;
          new_dw->get(pX,           lb->pXLabel,           pset);
          new_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);
          new_dw->get(pDefGrad,     lb->pDefGradLabel,     pset);

          constParticleVariable<Vector> pDisp;
          new_dw->get(pDisp,        lb->pDispLabel,        pset);

          ParticleVariable<Vector> pExternalForce;
          new_dw->getModifiable(pExternalForce, lb->pExternalForceLabel, pset);

          for (auto idx : *pset) {
            if (pLoadCurveID[idx] == nofPressureBCs) {
              pExternalForce[idx] = pbc->getForceVector(pX[idx], pDisp[idx], forcePerPart,
                                                        time, pDefGrad[idx]);
            }
          }

        } // matl loop
      }  // patch loop
    }
  }
}

void 
ImpMPM::scheduleComputeStableTimestep(const LevelP& level,
                                      SchedulerP& sched)
{
  printSchedule(level,cout_doing,"IMPM::scheduleComputeStableTimestep");

  Task* t = scinew Task("ImpMPM::actuallyComputeStableTimestep",
                        this, &ImpMPM::actuallyComputeStableTimestep);

  const MaterialSet* matls = d_sharedState->allMPMMaterials();
  
  if (flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()) ) {
    t->requires(Task::OldDW, d_sharedState->get_delt_label());
    t->requires(Task::NewDW, lb->pVelocityLabel,   Ghost::None);
  }

  // compute a delT on all levels even levels where impm is not running
  t->computes( d_sharedState->get_delt_label(),level.get_rep() );

  sched->addTask(t,level->eachPatch(), matls);
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
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels()) ) {
    new_dw->put(delt_vartype(999), lb->delTLabel, level);
    return;
  }  
  
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
   
    printTask(patches, patch,cout_doing,"Doing ImpMPM::actuallyComputeStableTimestep");

    if (d_numIterations==0) {
      new_dw->put(delt_vartype(d_initialDt), lb->delTLabel, patch->getLevel());
    } else {
      Vector dx = patch->dCell();
      delt_vartype old_delT;
      old_dw->get(old_delT, d_sharedState->get_delt_label(), patch->getLevel());

      int numMPMMatls=d_sharedState->getNumMPMMatls();
      for(int m = 0; m < numMPMMatls; m++){
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
        int dwindex = mpm_matl->getDWIndex();

        ParticleSubset* pset = new_dw->getParticleSubset(dwindex, patch);

        constParticleVariable<Vector> pVelocity;
        new_dw->get(pVelocity, lb->pVelocityLabel, pset);

        Vector ParticleSpeed(1.e-12,1.e-12,1.e-12);

        for (auto idx : *pset) {
          ParticleSpeed=Vector(Max(std::abs(pVelocity[idx].x()),ParticleSpeed.x()),
                               Max(std::abs(pVelocity[idx].y()),ParticleSpeed.y()),
                               Max(std::abs(pVelocity[idx].z()),ParticleSpeed.z()));
        }
        ParticleSpeed = dx/ParticleSpeed;
        double delT_new = .8*ParticleSpeed.minComponent();

        double old_dt=old_delT;
        if (d_numIterations <= flags->d_numItersToIncreaseDelT) {
          old_dt = flags->d_delTIncreaseFactor*old_delT;
        }
        if (d_numIterations >= flags->d_numItersToDecreaseDelT) {
          old_dt = flags->d_delTDecreaseFactor*old_delT;
        }
        delT_new = std::min(delT_new, old_dt);

        new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
      }
    }
  }
}

void
ImpMPM::scheduleTimeAdvance(const LevelP& level, SchedulerP& sched)
{
  if (!flags->doMPMOnLevel(level->getIndex(), level->getGrid()->numLevels())) {
    return;
  }

  const MaterialSet* matls = d_sharedState->allMPMMaterials();
  LoadBalancer* loadbal = sched->getLoadBalancer();
  d_perproc_patches = loadbal->getPerProcessorPatchSet(level);
  d_perproc_patches->addReference();

  scheduleComputeParticleBodyForce(    sched, d_perproc_patches,           matls);
  scheduleApplyExternalLoads(          sched, d_perproc_patches,           matls);
  scheduleInterpolateParticlesToGrid(  sched, d_perproc_patches, one_matl, matls);

  scheduleFindSurfaceParticles(           sched, d_perproc_patches,        matls);

  scheduleDestroyMatrix(                  sched, d_perproc_patches, matls, false);
  scheduleCreateMatrix(                   sched, d_perproc_patches, matls);
  scheduleDestroyHCMatrix(                sched, d_perproc_patches, matls);
  scheduleCreateHCMatrix(                 sched, d_perproc_patches, matls);
  scheduleApplyBoundaryConditions(        sched, d_perproc_patches, matls);
  scheduleApplyHCBoundaryConditions(      sched, d_perproc_patches, matls);
  scheduleComputeContact(                 sched, d_perproc_patches, matls);
  scheduleFindFixedDOF(                   sched, d_perproc_patches, matls);
  scheduleFindFixedHCDOF(                 sched, d_perproc_patches, matls);
  scheduleFormHCStiffnessMatrix(          sched, d_perproc_patches, matls);
  scheduleFormHCQ(                        sched, d_perproc_patches, matls);
  scheduleAdjustHCQAndHCKForBCs(          sched, d_perproc_patches, matls);
  scheduleSolveForTemp(                   sched, d_perproc_patches, matls);
  scheduleGetTemperatureIncrement(        sched, d_perproc_patches, matls);

  scheduleIterate(                        sched, level, d_perproc_patches, matls);

  if (!flags->d_doGridReset) {
    scheduleUpdateTotalDisplacement(      sched, d_perproc_patches, matls);
  }

  scheduleComputeDeformationGradient(     sched, d_perproc_patches, matls);
  scheduleComputeStressTensor(            sched, d_perproc_patches, matls);
  scheduleComputeAcceleration(            sched, d_perproc_patches, matls);
  scheduleInterpolateToParticlesAndUpdate(sched, d_perproc_patches, matls);
  scheduleInterpolateStressToGrid(        sched, d_perproc_patches, matls);

  sched->scheduleParticleRelocation(level, 
                                    lb->pXLabel_preReloc, 
                                    d_sharedState->d_particleState_preReloc,
                                    lb->pXLabel, 
                                    d_sharedState->d_particleState,
                                    lb->pParticleIDLabel, matls);
}

void
ImpMPM::scheduleComputeParticleBodyForce(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeParticleBodyForce");
  Task* t=scinew Task("IMPM::computeParticleBodyForce",
                      this, &ImpMPM::computeParticleBodyForce);

  t->requires(Task::OldDW, lb->pXLabel,        Ghost::None);
  t->requires(Task::OldDW, lb->pVelocityLabel, Ghost::None);
  t->computes(lb->pBodyForceAccLabel_preReloc);
  t->computes(lb->pCoriolisImportanceLabel_preReloc);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeParticleBodyForce(const ProcessorGroup* ,
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  // Get the MPM flags and make local copies
  Uintah::Point rotation_center = flags->d_coordRotationCenter;
  Uintah::Vector rotation_axis = flags->d_coordRotationAxis;
  double rotation_speed = flags->d_coordRotationSpeed;
  //Uintah::Point body_ref_point = flags->d_coord_rotation_body_ref_point;

  // Compute angular velocity vector (omega)
  Uintah::Vector omega = rotation_axis*rotation_speed;

  // Loop thru patches 
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing,"Doing computeParticleBodyForce");

    // Loop thru materials
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMPMMatls; m++) {

      // Get the material ID
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      // Get the particle subset
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      // Create space for particle body force
      ParticleVariable<Vector> pBodyForceAcc;
      new_dw->allocateAndPut(pBodyForceAcc, lb->pBodyForceAccLabel_preReloc, pset);

      // Create space for particle coriolis importance
      ParticleVariable<double> pCoriolisImportance;
      new_dw->allocateAndPut(pCoriolisImportance, lb->pCoriolisImportanceLabel_preReloc, pset);

      // Don't do much if coord rotation is off
      if (!flags->d_useCoordRotation) {

        // Iterate over the particles
        for (auto pidx : *pset) {

          // Compute the body force acceleration (g)
          pBodyForceAcc[pidx] = flags->d_gravity;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[pidx] = 0.0;
        } // particle loop

      } else { // Use coordinate rotation

        // Get the particle data
        constParticleVariable<Point> pPosition;
        old_dw->get(pPosition, lb->pXLabel, pset);

        constParticleVariable<Vector> pVelocity;
        old_dw->get(pVelocity, lb->pVelocityLabel, pset);

        for (auto pidx : *pset) {

          // Compute the Coriolis term (omega x v)
          Vector coriolis_accel = Uintah::Cross(omega, pVelocity[pidx])*2.0;

          // Compute the centrifugal term (omega x omega x r)
          // Simplified version where body ref point is not needed
          Vector rVec = pPosition[pidx] - rotation_center;
          Vector omega_x_r = Uintah::Cross(omega, rVec);
          Vector centrifugal_accel = Uintah::Cross(omega, omega_x_r);

          // Compute the body force acceleration (g - omega x omega x r - 2 omega x v)
          pBodyForceAcc[pidx] = flags->d_gravity - centrifugal_accel - coriolis_accel;

          // Compute relative importance of Coriolis term
          pCoriolisImportance[pidx] =
            coriolis_accel.length()/(centrifugal_accel.length() + coriolis_accel.length());

        } // particle loop
      } // end if coordinate rotation
    } // matl loop
  }  // patch loop
}

void 
ImpMPM::scheduleApplyExternalLoads(SchedulerP& sched,
                                   const PatchSet* patches,
                                   const MaterialSet* matls)
                                                                                
{
  printSchedule(patches,cout_doing,"IMPM::scheduleApplyExternalLoads");
  Task* t = scinew Task("IMPM::applyExternalLoads",
                        this, &ImpMPM::applyExternalLoads);
                                                                                
  t->requires(Task::OldDW, lb->pExternalForceLabel,     Ghost::None);
  t->requires(Task::OldDW, lb->pExternalHeatFluxLabel,  Ghost::None);
  t->requires(Task::OldDW, lb->pXLabel,                 Ghost::None);
  t->requires(Task::OldDW, lb->pDefGradLabel,           Ghost::None);
  t->requires(Task::OldDW, lb->pDispLabel,              Ghost::None);
  t->computes(             lb->pExtForceLabel_preReloc);
  t->computes(             lb->pExternalHeatRateLabel);
  t->computes(             lb->pExternalHeatFluxLabel_preReloc);
  if (flags->d_useLoadCurves) {
    t->requires(Task::OldDW, lb->pLoadCurveIDLabel,    Ghost::None);
    t->computes(             lb->pLoadCurveIDLabel_preReloc);
  }

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::applyExternalLoads(const ProcessorGroup* ,
                           const PatchSubset* patches,
                           const MaterialSubset*,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  // Get the current time
  double time = d_sharedState->getElapsedTime();

  if (cout_doing.active()) {
    cout_doing << "Current Time (applyExternalLoads) = " << time << endl;
  }
                                                                                
  // Calculate the force vector at each particle for each bc
  std::vector<double> forceMagPerPart;
  std::vector<PressureBC*> pbcP;
  std::vector<double> heatFluxMagPerPart;
  std::vector<HeatFluxBC*> hfbcP;
  if (flags->d_useLoadCurves) {

    // Currently, only one load curve at a time is supported, but
    // I've left the infrastructure in place to go to multiple
    for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
      auto bcType = particleBC->getType();
      if (bcType == "Pressure") {

        // std::cerr << "Pressure BCs is being supported in ImpMPM" << endl;
        PressureBC* pbc = dynamic_cast<PressureBC*>(particleBC.get());
        pbcP.push_back(pbc);
        forceMagPerPart.push_back(pbc->forcePerParticle(time));

      } else if (bcType == "HeatFlux") {

        HeatFluxBC* hfbc = dynamic_cast<HeatFluxBC*>(particleBC.get());
        #if 0
        std::cout << *hfbc << endl;
        std::cout << "hfbc type = " << hfbc->getType() << endl;
        std::cout << "surface area = " << hfbc->getSurfaceArea() << endl;
        std::cout << "heat flux = " << hfbc->heatflux(time) << endl;
        std::cout << "flux per particle = " << hfbc->fluxPerParticle(time) << endl;
        #endif
        hfbcP.push_back(hfbc);
        heatFluxMagPerPart.push_back(hfbc->fluxPerParticle(time));
      }
    }
  }
                                                                                
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing applyExternalLoads");
    
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      
      constParticleVariable<Point>    pX;
      constParticleVariable<double>   pExternalHeatFlux;
      constParticleVariable<Vector>   pDisp, pExternalForce;
      constParticleVariable<Matrix3>  pDefGrad;
      old_dw->get(pX,                lb->pXLabel,                pset);
      old_dw->get(pExternalHeatFlux, lb->pExternalHeatFluxLabel, pset);
      old_dw->get(pDisp,             lb->pDispLabel,             pset);
      old_dw->get(pExternalForce,    lb->pExternalForceLabel,    pset);
      old_dw->get(pDefGrad,          lb->pDefGradLabel,          pset);

      ParticleVariable<double> pExternalHeatFlux_new, pExtHeatRate;
      ParticleVariable<Vector> pExternalForce_new;
      new_dw->allocateAndPut(pExternalHeatFlux_new,
                             lb->pExternalHeatFluxLabel_preReloc, pset);
      new_dw->allocateAndPut(pExtHeatRate, 
                             lb->pExternalHeatRateLabel,          pset);
      new_dw->allocateAndPut(pExternalForce_new,
                             lb->pExtForceLabel_preReloc,         pset);
      
      for (auto idx : *pset) {
        pExternalHeatFlux_new[idx] = 0.;
        pExtHeatRate[idx] = 0.0;
        #if 0
        // Prescribe an external heat rate to some particles
        if(pX[idx].x()*pX[idx].x() + pX[idx].y()*pX[idx].y() > 0.0562*0.0562 ||
           pX[idx].z()>.0562 || pX[idx].z()<-.0562){
          pExtHeatRate[idx]=0.001;
        }
        #endif
      }

      if (flags->d_useLoadCurves) {
        bool do_PressureBCs=false;
        for (auto particleBC : MPMPhysicalBCFactory::mpmPhysicalBCs) {
          auto bcType = particleBC->getType();
          if (bcType == "Pressure") {
            do_PressureBCs = true;
          }
        }
        
        // Get the load curve data
        constParticleVariable<int> pLoadCurveID;
        old_dw->get(pLoadCurveID, lb->pLoadCurveIDLabel, pset);
        
        if (do_PressureBCs) {
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx]-1;
            if (loadCurveID < 0) {
              pExternalForce_new[idx] = pExternalForce[idx];
            } else {
              PressureBC* pbc = pbcP[loadCurveID];
              double force = forceMagPerPart[loadCurveID];
              pExternalForce_new[idx] = pbc->getForceVector(pX[idx], pDisp[idx],
                                                            force,time, pDefGrad[idx]);
            }
          }
        } else {
          for (auto idx : *pset) {
            pExternalForce_new[idx] = pExternalForce[idx]
              *flags->d_forceIncrementFactor;
          }
        } //end do_PressureBCs
        
        if (!heatFluxMagPerPart.empty()) {
          //double mag = heatFluxMagPerPart[0];
          //std::cout << "heat flux mag = " << mag << endl;
          for (auto idx : *pset) {
            int loadCurveID = pLoadCurveID[idx]-1;
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
                               lb->pLoadCurveIDLabel_preReloc, pset);
        pLoadCurveID_new.copyData(pLoadCurveID);

      } else { //not use pLoadCurve
        
        for (auto idx : *pset) {
          pExternalForce_new[idx] = pExternalForce[idx]
            *flags->d_forceIncrementFactor;
          pExternalHeatFlux_new[idx] = pExternalHeatFlux[idx];
        }
      }
    } // matl loop

    printTask(patches, patch, cout_doing, "Completed applyExternalLoads");

  }  // patch loop
}

void 
ImpMPM::scheduleInterpolateParticlesToGrid(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSubset* one_matl,
                                           const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleInterpolateParticlesToGrid");
  Task* t = scinew Task("ImpMPM::interpolateParticlesToGrid",
                        this, &ImpMPM::interpolateParticlesToGrid);

  t->requires(Task::OldDW, lb->pMassLabel,                  Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pVolumeLabel,                Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pTemperatureLabel,           Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pXLabel,                     Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pVelocityLabel,              Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pAccelerationLabel,          Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pSizeLabel,                  Ghost::AroundNodes, 1);
  t->requires(Task::OldDW, lb->pDefGradLabel,               Ghost::AroundNodes, 1);

  t->requires(Task::NewDW, lb->pBodyForceAccLabel_preReloc, Ghost::AroundNodes, 1);
  t->requires(Task::NewDW, lb->pExtForceLabel_preReloc,     Ghost::AroundNodes, 1);
  t->requires(Task::NewDW, lb->pExternalHeatRateLabel,      Ghost::AroundNodes, 1);
  t->requires(Task::NewDW, lb->pExternalHeatFluxLabel_preReloc, 
                                                            Ghost::AroundNodes, 1);

  if (!flags->d_doGridReset) {
    t->requires(Task::OldDW, lb->gDisplacementLabel,        Ghost::None);
    t->computes(lb->gDisplacementLabel);
  }

  t->requires(Task::OldDW,lb->NC_CCweightLabel, one_matl,Ghost::AroundCells,1);
  
  if (!flags->d_tempSolve) {
    t->requires(Task::OldDW, lb->gTemperatureLabel, one_matl, Ghost::None, 0);
  }

  t->computes(lb->gMassLabel,        d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);
  t->computes(lb->gVolumeLabel,      d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain);

  t->computes(lb->gMassLabel);
  t->computes(lb->gMassAllLabel);
  t->computes(lb->gVolumeLabel);
  t->computes(lb->gVelocityOldLabel);
  t->computes(lb->gVelocityLabel);
  t->computes(lb->dispNewLabel);
  t->computes(lb->gAccelerationLabel);
  t->computes(lb->gBodyForceLabel);
  t->computes(lb->gExternalForceLabel);
  t->computes(lb->gInternalForceLabel);
  t->computes(lb->TotalMassLabel);
  t->computes(lb->gTemperatureLabel,one_matl);
  t->computes(lb->gExternalHeatRateLabel);
  t->computes(lb->gExternalHeatFluxLabel);
  t->computes(lb->NC_CCweightLabel, one_matl);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);

  if (flags->d_projectHeatSource) {
    scheduleComputeCCVolume(           sched, patches, one_matl, matls);
    scheduleProjectCCHeatSourceToNodes(sched, patches, one_matl, matls);
  }
}

void 
ImpMPM::interpolateParticlesToGrid(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  Ghost::GhostType gac   = Ghost::AroundCells;
  Ghost::GhostType gnone = Ghost::None;

  static int timestep=0;

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing interpolateParticlesToGrid");

    auto interpolator = flags->d_interpolator->clone(patch);
    int i_size = interpolator->size();
    std::vector<IntVector> ni(i_size);
    std::vector<double> S(i_size);

    NCVariable<double> gMass_global, gVolume_global;
    new_dw->allocateAndPut(gMass_global, lb->gMassLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    new_dw->allocateAndPut(gVolume_global, lb->gVolumeLabel,
                           d_sharedState->getAllInOneMatl()->get(0), patch);
    gMass_global.initialize(d_SMALL_NUM_MPM);
    gVolume_global.initialize(d_SMALL_NUM_MPM);

    NCVariable<double> gMass_sum, gVolume_sum, gSpecificHeat;
    NCVariable<Vector> gVelocity_sum_old, gAcceleration_sum, 
                       gBodyForce_sum, gExternalForce_sum;
    new_dw->allocateTemporary(gMass_sum,          patch, gnone, 0);
    new_dw->allocateTemporary(gVolume_sum,        patch, gnone, 0);
    new_dw->allocateTemporary(gSpecificHeat,      patch, gnone, 0);
    new_dw->allocateTemporary(gVelocity_sum_old,  patch, gnone, 0);
    new_dw->allocateTemporary(gAcceleration_sum,  patch, gnone, 0);
    new_dw->allocateTemporary(gBodyForce_sum,     patch, gnone, 0);
    new_dw->allocateTemporary(gExternalForce_sum, patch, gnone, 0);
    gMass_sum.initialize(0.);
    gVolume_sum.initialize(0.);
    gSpecificHeat.initialize(0.);
    gVelocity_sum_old.initialize(Vector(0.,0.,0.));
    gAcceleration_sum.initialize(Vector(0.,0.,0.));
    gBodyForce_sum.initialize(Vector(0.,0.,0.));
    gExternalForce_sum.initialize(Vector(0.,0.,0.));

    constNCVariable<double> gTemperature_old;
    NCVariable<double>      gTemperature;
    bool switching_to_implicit_from_explicit = false;

    if (!flags->d_tempSolve) {
      if (old_dw->exists(lb->gTemperatureLabel, 0, patch)) {
        old_dw->get(gTemperature_old, lb->gTemperatureLabel, 0, patch, gnone, 0);
        switching_to_implicit_from_explicit = true;
      }
    }
    new_dw->allocateAndPut(gTemperature, lb->gTemperatureLabel, 0, patch);
    gTemperature.initialize(0.0);
    
    // carry forward interpolation weight
    IntVector low = patch->getExtraNodeLowIndex();
    IntVector hi  = patch->getExtraNodeHighIndex();
    constNCVariable<double> NC_CCweight;
    NCVariable<double>      NC_CCweight_copy;
    old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gac, 1);
    new_dw->allocateAndPut(NC_CCweight_copy, lb->NC_CCweightLabel, 0, patch);
    NC_CCweight_copy.copyPatch(NC_CCweight, low, hi);

    int numMatls = d_sharedState->getNumMPMMatls();
    NCdoubleArray gMass(numMatls), gVolume(numMatls), gExternalHeatRate(numMatls),
                  gExternalHeatFlux(numMatls), gMass_all(numMatls);
    NCVectorArray gVelocity_old(numMatls), gVelocity(numMatls), gAcceleration(numMatls),
                  dispNew(numMatls), gBodyForce(numMatls), gExternalForce(numMatls),
                  gInternalForce(numMatls);

    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      double Cp = mpm_matl->getSpecificHeat();
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch,
                                                       Ghost::AroundNodes, 1,
                                                       lb->pXLabel);

      constParticleVariable<Point>   pX;
      constParticleVariable<double>  pMass, pVolume, pTemperature,
                                     pExternalHeatRate, pExternalHeatFlux;
      constParticleVariable<Vector>  pVelocity, pAcceleration, pBodyForceAcc, 
                                     pExternalForce;
      constParticleVariable<Matrix3> pSize, pDefGrad;

      old_dw->get(pX,                lb->pXLabel,                         pset);
      old_dw->get(pMass,             lb->pMassLabel,                      pset);
      old_dw->get(pSize,             lb->pSizeLabel,                      pset);
      old_dw->get(pVolume,           lb->pVolumeLabel,                    pset);
      old_dw->get(pVelocity,         lb->pVelocityLabel,                  pset);
      old_dw->get(pTemperature,      lb->pTemperatureLabel,               pset);
      old_dw->get(pAcceleration,     lb->pAccelerationLabel,              pset);
      old_dw->get(pDefGrad,          lb->pDefGradLabel,                   pset);
      new_dw->get(pBodyForceAcc,     lb->pBodyForceAccLabel_preReloc,     pset);
      new_dw->get(pExternalForce,    lb->pExtForceLabel_preReloc,         pset);
      new_dw->get(pExternalHeatRate, lb->pExternalHeatRateLabel,          pset);
      new_dw->get(pExternalHeatFlux, lb->pExternalHeatFluxLabel_preReloc, pset);

      new_dw->allocateAndPut(gMass[m],             lb->gMassLabel,             matID, patch);
      new_dw->allocateAndPut(gMass_all[m],         lb->gMassAllLabel,          matID, patch);
      new_dw->allocateAndPut(gVolume[m],           lb->gVolumeLabel,           matID, patch);
      new_dw->allocateAndPut(gVelocity_old[m],     lb->gVelocityOldLabel,      matID, patch);
      new_dw->allocateAndPut(gVelocity[m],         lb->gVelocityLabel,         matID, patch);
      new_dw->allocateAndPut(dispNew[m],           lb->dispNewLabel,           matID, patch);
      new_dw->allocateAndPut(gAcceleration[m],     lb->gAccelerationLabel,     matID, patch);
      new_dw->allocateAndPut(gBodyForce[m],        lb->gBodyForceLabel,        matID, patch);
      new_dw->allocateAndPut(gExternalForce[m],    lb->gExternalForceLabel,    matID, patch);
      new_dw->allocateAndPut(gInternalForce[m],    lb->gInternalForceLabel,    matID, patch);
      new_dw->allocateAndPut(gExternalHeatRate[m], lb->gExternalHeatRateLabel, matID, patch);
      new_dw->allocateAndPut(gExternalHeatFlux[m], lb->gExternalHeatFluxLabel, matID, patch);

      if (!flags->d_doGridReset) {
        constNCVariable<Vector> gDisplacement;
        NCVariable<Vector>      gDisplacement_new;
        old_dw->get(gDisplacement, lb->gDisplacementLabel, matID, patch, gnone, 0);
        new_dw->allocateAndPut(gDisplacement_new, lb->gDisplacementLabel, matID, patch);
        gDisplacement_new.copyData(gDisplacement);
      }

      gMass[m].initialize(d_SMALL_NUM_MPM);
      gVolume[m].initialize(0);
      gExternalHeatRate[m].initialize(0.0);
      gExternalHeatFlux[m].initialize(0.0);
      gVelocity_old[m].initialize(Vector(0,0,0));
      gVelocity[m].initialize(Vector(0,0,0));
      dispNew[m].initialize(Vector(0,0,0));
      gAcceleration[m].initialize(Vector(0,0,0));
      gBodyForce[m].initialize(Vector(0,0,0));
      gExternalForce[m].initialize(Vector(0,0,0));
      gInternalForce[m].initialize(Vector(0,0,0));

      double totalMass = 0;
      for (auto idx : *pset) {

        interpolator->findCellAndWeights(pX[idx], ni, S, pSize[idx], pDefGrad[idx]);

        auto pMassAcc    = pAcceleration[idx]*pMass[idx];
        auto pMom        = pVelocity[idx]*pMass[idx];
        totalMass += pMass[idx];

        // Add each particles contribution to the local mass & velocity 
        // Must use the node indices
        for (int k = 0; k < i_size; k++) {
          auto node = ni[k];
          if (patch->containsNode(node)) {
            gMass[m][node]             += pMass[idx]                      * S[k];
            gMass_global[node]         += pMass[idx]                      * S[k];
            gVolume[m][node]           += pVolume[idx]                    * S[k];
            gVolume_global[node]       += pVolume[idx]                    * S[k];
            gSpecificHeat[node]        += Cp                              * S[k];
            gExternalHeatRate[m][node] += pExternalHeatRate[idx]          * S[k];
            gExternalHeatFlux[m][node] += pExternalHeatFlux[idx]          * S[k];
            if (!flags->d_tempSolve) {
              gTemperature[node]       += pTemperature[idx]  * pMass[idx] * S[k];
            }
            gBodyForce[m][node]        += pBodyForceAcc[idx] * pMass[idx] * S[k];
            gExternalForce[m][node]    += pExternalForce[idx]             * S[k];
            gVelocity_old[m][node]     += pMom                            * S[k];
            gAcceleration[m][node]     += pMassAcc                        * S[k];
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
          gMass_sum[c]          += gMass[m][c];
          gVolume_sum[c]        += gVolume[m][c];
          gBodyForce_sum[c]     += gBodyForce[m][c];
          gExternalForce_sum[c] += gExternalForce[m][c];
          gVelocity_sum_old[c]  += gVelocity_old[m][c];
          gAcceleration_sum[c]  += gAcceleration[m][c];
        }
      }

      new_dw->put(sum_vartype(totalMass), lb->TotalMassLabel);
    }  // End loop over materials

    // Single material temperature field
    if (!flags->d_tempSolve) {

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

        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
        int matID = mpm_matl->getDWIndex();
        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch,
                                                         Ghost::AroundNodes, 1,
                                                         lb->pXLabel);
        
        constParticleVariable<double>  pTemperature;
        constParticleVariable<Point>   pX;
        constParticleVariable<Matrix3> pSize, pDefGrad;
       
        old_dw->get(pX,             lb->pXLabel,                 pset);
        old_dw->get(pTemperature,   lb->pTemperatureLabel,       pset);
        old_dw->get(pSize,          lb->pSizeLabel,              pset);
        old_dw->get(pDefGrad,       lb->pDefGradLabel,           pset);
 
        for (auto idx : *pset) {

          interpolator->findCellAndWeights(pX[idx], ni_cell, S, pSize[idx], pDefGrad[idx]);
        
          ParticleTempShape ptshape;
          ptshape.pTemperature = pTemperature[idx];
          ptshape.cellNodes = ni_cell;
          ptshape.shapeFnValues = S;

          IntVector cellID = ni_cell[0];
          cell_map.insert(CellParticleTempPair(cellID,ptshape));
        }
      }

      #ifdef debug
      std::cout << "size of cell_map before = " << cell_map.size() << endl;
      #endif
      for (auto iter = cell_map.begin(); iter != cell_map.end(); 
                iter = cell_map.upper_bound(iter->first)) {
        #ifdef debug
        std::cout << "cell = " << iter->first << " temp = " 
                  << iter->second.particleTemps << " count = " 
                  << cell_map.count(iter->first) << endl;
        #endif

        if (cell_map.count(iter->first) < 8 ) {
          #ifdef debug
          std::cout << "Inserting cell " << iter->first << " into sparse_cell_map\n" ;
          #endif
          auto& smap = sparse_cell_map[cell_map.count(iter->first)-1];
          auto eq_range = cell_map.equal_range(iter->first);
          IntVector cellID = iter->first;
          smap.insert(eq_range.first, eq_range.second);
          cell_map.erase(eq_range.first, eq_range.second);
        }  
      }
      #ifdef debug
      std::cout << "size of cell_map after = " << cell_map.size() << endl;
      for (int i = 0; i < 7; i++) {
        std::cout << "size of sparse_cell_map[" << i << "] after = " 
             << sparse_cell_map[i].size() << endl;
      }
      #endif

      // Process all of the cells with 8 particles in them
      FastMatrix A(8,8);
      double B[8];
      #ifdef debug    
      std::cout << "Working on cells with 8 particles" << endl;
      #endif
      for (auto iter = cell_map.begin(); iter != cell_map.end(); 
                iter = cell_map.upper_bound(iter->first)) {
        #ifdef debug        
        std::cout << "working on cell " << iter->first << endl;
        #endif

        auto eq_range = cell_map.equal_range(iter->first);
        int count = 0;

        ParticleTempShape ptshape;
        for (auto it = eq_range.first; it != eq_range.second; ++it) {
          ptshape = it->second;
          B[count] = ptshape.pTemperature;
          for (int j = 0; j < 8; j++) {
            A(count,j) = ptshape.shapeFnValues[j];
          }
          count++;
        }
        
        A.destructiveSolve(B);
        A.zero();
        for (int j = 0; j < 8; j++) {
          if (patch->containsNode(ptshape.cellNodes[j])) {
            gTemperature[ptshape.cellNodes[j]] = B[j];
            #ifdef debug
            std::cout << "gTemperature[" << ptshape.cellNodes[j] << "] = " 
                 << gTemperature[ptshape.cellNodes[j]] << endl;
            #endif
          }
        }
      }

      // Work on the cells that have fewer than 8 particles in them
      for (int i = 6; i >= 0; i--) {
        #ifdef debug
        std::cout << "Working on cells with " << i + 1 << " particles" << endl;
        #endif
        auto& smap = sparse_cell_map[i];

        FastMatrix A(8,8);
        double B[8];
        for (auto it = smap.begin(); it != smap.end();
                  it = smap.upper_bound(it->first)) {
          #ifdef debug
          std::cout << "working on cell " << it->first << endl;
          #endif
          
          auto eq_range = smap.equal_range(it->first);
          int count = 0;
          A.zero();
          for (int i = 0; i < 8; i++) {
            B[i] = 0.;
          }
          ParticleTempShape ptshape;
          for (auto it = eq_range.first; it != eq_range.second; ++it) {
            ptshape = it->second;
            B[count] = ptshape.pTemperature;
            for (int j = 0; j < 8; j++) {
              A(count, j) = ptshape.shapeFnValues[j];
            }
            count++;
          }

          FastMatrix A_t(8,8);
          A_t.transpose(A);
          double A_tB[8];
          A_t.multiply(B,A_tB);
          FastMatrix A_tA(8,8);
          A_tA.multiply(A_t,A);
          
          for (int i = 0; i < 8; i++) {
            if (patch->containsNode(ptshape.cellNodes[i])) {
              if (gTemperature[ptshape.cellNodes[i]] != 0.0) {
                #ifdef debug
                std::cout << "i = " << i << " setting gTemperature[" 
                     << ptshape.cellNodes[i] << "]=" 
                     << gTemperature[ptshape.cellNodes[i]] << endl;
                #endif
                for (int j = 0; j < 8; j++) {
                  A_tA(i,j) = 0.;
                }
                A_tA(i,i) = 1.0;
                A_tB[i] = gTemperature[ptshape.cellNodes[i]];
              }
            }
          }
          
          A_tA.destructiveSolve(A_tB);
          for (int j = 0; j < 8; j++) {
            if (patch->containsNode(ptshape.cellNodes[j])) {
              gTemperature[ptshape.cellNodes[j]] = A_tB[j];
              #ifdef debug
              std::cout << "gTemperature[" << ptshape.cellNodes[j] << "] = " 
                   << gTemperature[ptshape.cellNodes[j]] << endl;
              #endif
            }
          }
        }
      }
    }
    if (!flags->d_interpolateParticleTempToGridEveryStep) {
      if (timestep > 0) {
        if (!switching_to_implicit_from_explicit) {
          for(auto iter = patch->getNodeIterator(); !iter.done(); iter++){
            IntVector c = *iter;
            gTemperature[c] = gTemperature_old[c];
          }
        }
      }
    }

    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if(!mpm_matl->getIsRigid()){
        for(auto iter = patch->getNodeIterator(); !iter.done(); iter++){
          IntVector c = *iter;
          gMass_all[m][c]      = gMass_sum[c];
          gVolume[m][c]        = gVolume_sum[c];
          gBodyForce[m][c]     = gBodyForce_sum[c];
          gExternalForce[m][c] = gExternalForce_sum[c];
          gVelocity_old[m][c]  = gVelocity_sum_old[c]/(gMass_sum[c] + 1.e-200);
          gAcceleration[m][c]  = gAcceleration_sum[c]/(gMass_sum[c] + 1.e-200);
        }
      }
    }  // End loop over materials
  }  // End loop over patches
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

  Task* t = scinew Task("ImpMPM::findSurfaceParticles", this,
                        &ImpMPM::findSurfaceParticles);

  t->requires(Task::OldDW, lb->pSurfLabel, Ghost::AroundNodes, NGP);
  t->computes(lb->pSurfLabel_preReloc);

  sched->addTask(t, patches, matls);
}

void 
ImpMPM::findSurfaceParticles(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* ,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  auto numMPMMatls = d_sharedState->getNumMPMMatls();

  for (int p = 0; p<patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches, patch, cout_doing, "Doing findSurfaceParticles");

    for(int mat = 0; mat < numMPMMatls; mat++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(mat);
      int matID = mpm_matl->getDWIndex();

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

      constParticleVariable<double> pSurf_old;
      ParticleVariable<double> pSurf;

      old_dw->get(pSurf_old,        lb->pSurfLabel,          pset);
      new_dw->allocateAndPut(pSurf, lb->pSurfLabel_preReloc, pset);

      // For now carry forward the particle surface data
      for (auto particle : *pset) {
        pSurf[particle] = pSurf_old[particle];
      }
    }   // matl loop
  }    // patches
}

void 
ImpMPM::scheduleComputeCCVolume(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSubset* one_matl,
                                const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleComputeCCVolume");
  Task* t = scinew Task("ImpMPM::computeCCVolume",
                        this, &ImpMPM::computeCCVolume);
                                                                                
  t->requires(Task::OldDW,lb->NC_CCweightLabel, one_matl, Ghost::AroundCells, 1);
  t->requires(Task::NewDW,lb->gVolumeLabel,               Ghost::AroundCells, 1);

  t->computes(lb->cVolumeLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::computeCCVolume(const ProcessorGroup*,
                        const PatchSubset* patches,
                        const MaterialSubset* ,
                        DataWarehouse* old_dw,
                        DataWarehouse* new_dw)
{
  Ghost::GhostType gac = Ghost::AroundCells;

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing computeCCVolume");
      
    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gac, 1);

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      constNCVariable<double> gVolume;
      CCVariable<double> cVolume;

      new_dw->get(gVolume,            lb->gVolumeLabel, matID, patch, gac, 1);
      new_dw->allocateAndPut(cVolume, lb->cVolumeLabel, matID, patch);
      cVolume.initialize(1.e-20);

      for(auto iter = patch->getExtraCellIterator(); !iter.done(); iter++) {
        IntVector c = *iter;
        IntVector nodeIdx[8];
        patch->findNodesFromCell(c, nodeIdx);
        for (int in = 0; in < 8; in++) {
          cVolume[c] += NC_CCweight[nodeIdx[in]] * gVolume[nodeIdx[in]];
        }
      }
    }
  }
}

void 
ImpMPM::scheduleProjectCCHeatSourceToNodes(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSubset* one_matl,
                                           const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleProjectCCHeatSourceToNodes");
  Task* t = scinew Task("ImpMPM::projectCCHeatSourceToNodes",
                        this, &ImpMPM::projectCCHeatSourceToNodes);

  t->requires(Task::OldDW, lb->NC_CCweightLabel, one_matl, Ghost::AroundCells, 1);
  t->requires(Task::OldDW, lb->heatRate_CCLabel,           Ghost::AroundCells, 1);
  t->requires(Task::NewDW, lb->gVolumeLabel,               Ghost::AroundCells, 1);
  t->requires(Task::NewDW, lb->cVolumeLabel,               Ghost::AroundCells, 1);

  t->computes(lb->heatRate_CCLabel);
  t->modifies(lb->gExternalHeatRateLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::projectCCHeatSourceToNodes(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  Ghost::GhostType  gac = Ghost::AroundCells;

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches,patch,cout_doing,
              "Doing projectCCHeatSourceToNodes on patch\t\t");

    constNCVariable<double> NC_CCweight;
    old_dw->get(NC_CCweight, lb->NC_CCweightLabel, 0, patch, gac, 1);

    int numMPMMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      constNCVariable<double> gVolume;
      constCCVariable<double> cHeatRate, cVolume;

      new_dw->get(gVolume,   lb->gVolumeLabel,     matID, patch, gac, 1);
      old_dw->get(cHeatRate, lb->heatRate_CCLabel, matID, patch, gac, 1);
      new_dw->get(cVolume,   lb->cVolumeLabel,     matID, patch, gac, 1);

      NCVariable<double> gExternalHeatRate;
      CCVariable<double> cHeatRate_copy;

      new_dw->getModifiable(gExternalHeatRate, lb->gExternalHeatRateLabel, matID, patch);
      new_dw->allocateAndPut(cHeatRate_copy,   lb->heatRate_CCLabel,       matID, patch);
      cHeatRate_copy.copyData(cHeatRate);

      // Project  CC heat rate to nodes
      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        IntVector cIdx[8];
        patch->findCellsFromNode(node, cIdx);
        for (int ic = 0; ic < 8; ic++) {
          double solid_volume = cVolume[cIdx[ic]];
          gExternalHeatRate[node] += cHeatRate[cIdx[ic]]*(NC_CCweight[node]*gVolume[node])
                       /solid_volume;
        }
      }
    }
  }
}

void 
ImpMPM::scheduleDestroyMatrix(SchedulerP& sched,
                              const PatchSet* patches,
                              const MaterialSet* matls,
                              bool recursion)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleDestroyMatrix");
  Task* t = scinew Task("ImpMPM::destroyMatrix", this, &ImpMPM::destroyMatrix,
                         recursion);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::destroyMatrix(const ProcessorGroup*,
                      const PatchSubset* /*patches*/,
                      const MaterialSubset* ,
                      DataWarehouse* /* old_dw */,
                      DataWarehouse* /* new_dw */,
                      bool recursion)
{
  if (cout_doing.active()) {
    cout_doing <<"Doing destroyMatrix \t\t\t\t\t IMPM"  << "\n";
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

  t->requires(Task::OldDW, lb->pXLabel, Ghost::AroundNodes, 1);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::createMatrix(const ProcessorGroup*,
                     const PatchSubset* patches,
                     const MaterialSubset* ,
                     DataWarehouse* old_dw,
                     DataWarehouse* new_dw)
{
  std::map<int,int> dof_diag;
  d_solver->createLocalToGlobalMapping(UintahParallelComponent::d_myworld,
                                       d_perproc_patches,
                                       patches, 3, flags->d_8or27);
  int global_offset = 0;
  int numMatls = d_sharedState->getNumMPMMatls();
  int n8or27 = flags->d_8or27;

  for(int pp = 0; pp < patches->size(); pp++){
    const Patch* patch = patches->get(pp);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::createMatrix");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, n8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;

    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g, patch);
    
    //set the global offset if this is the first patch
    if (pp == 0) {
      global_offset = l2g[lowIndex];
    }

    CCVariable<int> visited;
    new_dw->allocateTemporary(visited, patch, Ghost::AroundCells, 1);
    visited.initialize(0);

    for (int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();    
      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch, 
                                                       Ghost::AroundNodes, 1,
                                                       lb->pXLabel);

      constParticleVariable<Point> pX;
      old_dw->get(pX, lb->pXLabel, pset);

      for (auto idx : *pset) {
        IntVector cell, ni[27];
        patch->findCell(pX[idx], cell);
        if (visited[cell] == 0) {
          visited[cell] = 1;
          patch->findNodesFromCell(cell, ni);
          std::vector<int> dof(0);
          int l2g_node_num;
          for (int k = 0; k < n8or27; k++) {
            if (patch->containsNode(ni[k]) ) {
              //subtract global offset in order to map into array correctly
              l2g_node_num = l2g[ni[k]] - global_offset; 
              dof.push_back(l2g_node_num);
              dof.push_back(l2g_node_num+1);
              dof.push_back(l2g_node_num+2);
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
ImpMPM::scheduleDestroyHCMatrix(SchedulerP& sched,
                                const PatchSet* patches,
                                const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleDestroyHCMatrix");
  heatConductionModel->scheduleDestroyHCMatrix(sched,patches,matls);
}

void 
ImpMPM::scheduleCreateHCMatrix(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleCreateHCMatrix");
  heatConductionModel->scheduleCreateHCMatrix(sched, patches, matls);
}

void 
ImpMPM::scheduleApplyBoundaryConditions(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleApplyBoundaryConditions");
  Task* t = scinew Task("ImpMPM::applyBoundaryCondition",
                        this, &ImpMPM::applyBoundaryConditions);

  t->modifies(lb->gVelocityOldLabel);
  t->modifies(lb->gAccelerationLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::applyBoundaryConditions(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* ,
                                DataWarehouse* /*old_dw*/,
                                DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing applyBoundaryConditions");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;

    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g,patch);

    // Apply grid boundary conditions to the velocity before storing the data
    IntVector offset =  IntVector(0,0,0);
    for (int m = 0; m < d_sharedState->getNumMPMMatls(); m++ ) {

      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();
      
      NCVariable<Vector> gAcceleration, gVelocity_old;
      new_dw->getModifiable(gVelocity_old, lb->gVelocityOldLabel,  matID, patch);
      new_dw->getModifiable(gAcceleration, lb->gAccelerationLabel, matID, patch);

      for(auto face = Patch::startFace; face <= Patch::endFace; 
               face=Patch::nextFace(face)) {

        if (patch->getBCType(face) == Patch::None) {

          int numChildren = patch->getBCDataArray(face)->getNumberChildren(matID);
          for (int child = 0; child < numChildren; child++) {

            Iterator nbound_ptr;
            Iterator nu;     // not used;
            
            BoundCondBaseP vel_bcs = patch->getArrayBCValues(face, matID, "Velocity", nu,
                                                             nbound_ptr, child);

            auto bc = std::dynamic_pointer_cast<BoundCond<Vector> >(vel_bcs);

            if (bc != 0) {
              if (bc->getBCType__NEW() == "Dirichlet") {
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] = bc->getValue();
                  gAcceleration[*nbound_ptr] = bc->getValue();
                }
                IntVector l,h;
                patch->getFaceNodes(face,0,l,h);
                for (NodeIterator it(l,h); !it.done(); it++) {
                  IntVector n = *it;
                  int l2g_node_num = l2g[n];
                  d_solver->d_DOF.insert(l2g_node_num);
                  d_solver->d_DOF.insert(l2g_node_num+1);
                  d_solver->d_DOF.insert(l2g_node_num+2);
                }
              }
            } 

            BoundCondBaseP sym_bcs  = patch->getArrayBCValues(face, matID, "Symmetric", nu,
                                                              nbound_ptr, child);
            auto sbc = std::dynamic_pointer_cast<BoundCond<NoValue> >(sym_bcs);

            if (sbc != 0) {
              if (face == Patch::xplus || face == Patch::xminus)
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] = 
                    Vector(0.,gVelocity_old[*nbound_ptr].y(),
                           gVelocity_old[*nbound_ptr].z());
                  gAcceleration[*nbound_ptr] = 
                    Vector(0.,gAcceleration[*nbound_ptr].y(),
                           gAcceleration[*nbound_ptr].z());
                }
              if (face == Patch::yplus || face == Patch::yminus)
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] = 
                    Vector(gVelocity_old[*nbound_ptr].x(),0.,
                           gVelocity_old[*nbound_ptr].z());
                  gAcceleration[*nbound_ptr] = 
                    Vector(gAcceleration[*nbound_ptr].x(),0.,
                           gAcceleration[*nbound_ptr].z());
                }
              if (face == Patch::zplus || face == Patch::zminus)
                for (nbound_ptr.reset(); !nbound_ptr.done(); nbound_ptr++) {
                  gVelocity_old[*nbound_ptr] = 
                    Vector(gVelocity_old[*nbound_ptr].x(),
                           gVelocity_old[*nbound_ptr].y(),0.);
                  gAcceleration[*nbound_ptr] = 
                    Vector(gAcceleration[*nbound_ptr].x(),
                           gAcceleration[*nbound_ptr].y(),0.);
                }
              IntVector l,h;
              patch->getFaceNodes(face,0,l,h);
              for(NodeIterator it(l,h); !it.done(); it++) {
                IntVector n = *it;
                // The DOF is an IntVector which is initially (0,0,0).
                // Inserting a 1 into any of the components indicates that 
                // the component should be inserted into the DOF array.
                IntVector DOF(0,0,0);
                if (face == Patch::xminus || face == Patch::xplus)
                  DOF=IntVector(std::max(DOF.x(),1),std::max(DOF.y(),0),std::max(DOF.z(),0));
                if (face == Patch::yminus || face == Patch::yplus)
                  DOF=IntVector(std::max(DOF.x(),0),std::max(DOF.y(),1),std::max(DOF.z(),0));
                if (face == Patch::zminus || face == Patch::zplus)
                  DOF=IntVector(std::max(DOF.x(),0),std::max(DOF.y(),0),std::max(DOF.z(),1));
                
                int l2g_node_num = l2g[n];
                if (DOF.x())
                  d_solver->d_DOF.insert(l2g_node_num);
                if (DOF.y())
                  d_solver->d_DOF.insert(l2g_node_num+1);
                if (DOF.z())
                  d_solver->d_DOF.insert(l2g_node_num+2);
              }
            } // endif (sbc)
          } // endfor child 
        } // endif patchBC
      } // endfor face
    } // endfor mpmmat
  } // endfor patch
}

void 
ImpMPM::scheduleApplyHCBoundaryConditions(SchedulerP& sched,
                                          const PatchSet* patches,
                                          const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleApplyHCBoundaryConditions");
  heatConductionModel->scheduleApplyHCBoundaryConditions(sched, patches, matls);
}

void 
ImpMPM::scheduleComputeContact(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleComputeContact");
  Task* t = scinew Task("ImpMPM::computeContact",
                         this, &ImpMPM::computeContact);

  t->requires(Task::OldDW, d_sharedState->get_delt_label());
  if (d_rigid_body) {
    t->requires(Task::NewDW, lb->gMassLabel,        Ghost::None);
    t->requires(Task::NewDW, lb->gVelocityOldLabel, Ghost::None);
    t->modifies(lb->dispNewLabel);
    if (!flags->d_doGridReset) {
      t->requires(Task::OldDW,lb->gDisplacementLabel, Ghost::None);
      t->modifies(lb->gDisplacementLabel);
    }
  }
  t->computes(lb->gContactLabel);  

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::computeContact(const ProcessorGroup*,
                       const PatchSubset* patches,
                       const MaterialSubset* ,
                       DataWarehouse* old_dw,
                       DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing <<"Doing computeContact on patch " << patch->getID()
                 <<"\t\t\t\t IMPM"<< "\n";
    }

    delt_vartype dt;

    int numMatls = d_sharedState->getNumMPMMatls();
    std::vector<NCVariable<int> >  contact(numMatls);
    for(int n = 0; n < numMatls; n++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( n );
      int matID = mpm_matl->getDWIndex();
      new_dw->allocateAndPut(contact[n], lb->gContactLabel, matID, patch);
      contact[n].initialize(0);
    }

    if (d_rigid_body) {
      constNCVariable<Vector> vel_rigid;
      constNCVariable<double> mass_rigid;
      int numMatls = d_sharedState->getNumMPMMatls();
      for(int n = 0; n < numMatls; n++){
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( n );
        if(mpm_matl->getIsRigid()){
          int matID = mpm_matl->getDWIndex();
          new_dw->get(vel_rigid, lb->gVelocityOldLabel, matID, patch, Ghost::None, 0);
          new_dw->get(mass_rigid,lb->gMassLabel,        matID, patch, Ghost::None, 0);
        }
      }

      // Get and modify non-rigid data
      for(int m = 0; m < numMatls; m++){
        MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
        int matID = mpm_matl->getDWIndex();
        NCVariable<Vector> dispNew;                     
        new_dw->getModifiable(dispNew, lb->dispNewLabel, matID, patch);

        delt_vartype dt;
        old_dw->get(dt, d_sharedState->get_delt_label(), patch->getLevel() );

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++){
          IntVector node = *iter;
          if (!compare(mass_rigid[node],0.0)) {
            dispNew[node] = Vector(vel_rigid[node].x()*dt*d_contact_dirs.x(),
                                   vel_rigid[node].y()*dt*d_contact_dirs.y(),
                                   vel_rigid[node].z()*dt*d_contact_dirs.z());
            contact[m][node] = 2;
          }
        } 

        if (!flags->d_doGridReset) {
          NCVariable<Vector> gDisplacement_new;
          constNCVariable<Vector> gDisplacement_old;
          new_dw->get(gDisplacement_old, lb->gDisplacementLabel, matID, patch, Ghost::None, 0);
          new_dw->getModifiable(gDisplacement_new, lb->gDisplacementLabel, matID, patch);
          for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
            IntVector node = *iter;
            gDisplacement_new[node] = gDisplacement_old[node] + dispNew[node];
          }
        }
      } // endfor numMatls
    } // endif rigid_body
  } // endfor patches
}

void 
ImpMPM::scheduleFindFixedDOF(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFindFixedDOF");
  Task* t = scinew Task("ImpMPM::findFixedDOF", this, 
                        &ImpMPM::findFixedDOF);

  t->requires(Task::NewDW, lb->gMassAllLabel, Ghost::None, 0);
  t->requires(Task::NewDW, lb->gContactLabel, Ghost::None, 0);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::findFixedDOF(const ProcessorGroup*, 
                     const PatchSubset* patches,
                     const MaterialSubset*, 
                     DataWarehouse* ,
                     DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    
    printTask(patches, patch,cout_doing,"Doing ImpMPM::findFixedDOF");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g,patch);

    bool firstTimeThrough=true;
    int numMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!mpm_matl->getIsRigid() && firstTimeThrough) { 
        firstTimeThrough=false;
        int matID = mpm_matl->getDWIndex();
        constNCVariable<double> mass;
        constNCVariable<int> contact;
        new_dw->get(mass,    lb->gMassAllLabel, matID, patch, Ghost::None, 0);
        new_dw->get(contact, lb->gContactLabel, matID, patch, Ghost::None, 0);

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          int l2g_node_num = l2g[node];

          // Just look on the grid to see if the gMass is 0 and then remove that
          if (compare(mass[node], 0.)) {
            d_solver->d_DOF.insert(l2g_node_num);
            d_solver->d_DOF.insert(l2g_node_num+1);
            d_solver->d_DOF.insert(l2g_node_num+2);
          }
          if (contact[node] == 2) {  // Rigid Contact imposed on these nodes
            for (int i=0; i<3; i++) {
              if (d_contact_dirs[i]==1) {
                d_solver->d_DOF.insert(l2g_node_num+i);  // specifically, these DOFs
              }
            }
          } // contact ==2
        } // node iterator
      } // if not rigid
    } // loop over matls
  } // patches
}

void 
ImpMPM::scheduleFindFixedHCDOF(SchedulerP& sched,
                               const PatchSet* patches,
                               const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFindFixedHCDOF");
  heatConductionModel->scheduleFindFixedHCDOF(sched,patches,matls);
}

void 
ImpMPM::scheduleFormHCStiffnessMatrix(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFormHCStiffnessMatrix");
  heatConductionModel->scheduleFormHCStiffnessMatrix(sched,patches,matls);
}

void 
ImpMPM::scheduleFormHCQ(SchedulerP& sched,
                        const PatchSet* patches,
                        const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFormHCQ");
  heatConductionModel->scheduleFormHCQ(sched,patches,matls);
}

void 
ImpMPM::scheduleAdjustHCQAndHCKForBCs(SchedulerP& sched,
                                      const PatchSet* patches,
                                      const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFormHCQAndHCKForBCs");
  heatConductionModel->scheduleAdjustHCQAndHCKForBCs(sched,patches,matls);
}

void 
ImpMPM::scheduleSolveForTemp(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleSolveForTemp");
  heatConductionModel->scheduleSolveForTemp(sched,patches,matls);
}

void 
ImpMPM::scheduleGetTemperatureIncrement(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleGetTemperatureIncrement");
  heatConductionModel->scheduleGetTemperatureIncrement(sched,patches,matls);
}

void 
ImpMPM::scheduleIterate(SchedulerP& sched,
                        const LevelP& level,
                        const PatchSet* patches, 
                        const MaterialSet*)
{
  d_recompileSubsched = true;
  printSchedule(patches,cout_doing,"IMPM::scheduleIterate");
  Task* task = scinew Task("ImpMPM::iterate", this, &ImpMPM::iterate,level,
                           sched.get_rep());

  task->hasSubScheduler();

  task->requires(Task::OldDW,lb->pXLabel,       Ghost::None, 0);
  task->requires(Task::OldDW,lb->pMassLabel,    Ghost::None, 0);
  task->requires(Task::OldDW,lb->pSizeLabel,    Ghost::None, 0);
  task->requires(Task::OldDW,lb->pVolumeLabel,  Ghost::None, 0);
  task->requires(Task::OldDW,lb->pDefGradLabel, Ghost::None, 0);

  task->modifies(lb->dispNewLabel);
  task->modifies(lb->gVelocityLabel);
  task->modifies(lb->gInternalForceLabel);
  if (!flags->d_doGridReset) {
    task->requires(Task::OldDW,lb->gDisplacementLabel, Ghost::None, 0);
    task->modifies(lb->gDisplacementLabel);
  }

  task->requires(Task::NewDW,lb->gVelocityOldLabel,    Ghost::None, 0);
  task->requires(Task::NewDW,lb->gMassAllLabel,        Ghost::None, 0);
  task->requires(Task::NewDW,lb->gBodyForceLabel,      Ghost::None, 0);
  task->requires(Task::NewDW,lb->gExternalForceLabel,  Ghost::None, 0);
  task->requires(Task::NewDW,lb->gAccelerationLabel,   Ghost::None, 0);
  task->requires(Task::NewDW,lb->gContactLabel,        Ghost::None, 0);

  if (flags->d_doMechanics) {
    int numMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);

      d_defGradComputer->addComputesAndRequires(task, mpm_matl, patches, true, false);

      ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
      cm->addComputesAndRequires(task, mpm_matl, patches, true, false);
    }
  }

  task->requires(Task::OldDW, d_sharedState->get_delt_label());

  task->setType(Task::OncePerProc);
  sched->addTask(task, patches, d_sharedState->allMaterials());
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
  d_subsched->setSimulationState(d_sharedState);
  d_subsched->advanceDataWarehouse(grid);
  d_subsched->setInitTimestep(true);
  const MaterialSet* matls = d_sharedState->allMPMMaterials();

  if (d_recompileSubsched) {
    int numOldDW = 3, numNewDW = 1;
    d_subsched->initialize(numOldDW, numNewDW);
    
    // This task only zeros out the stiffness matrix it doesn't free any memory.
    scheduleDestroyMatrix(           d_subsched, d_perproc_patches, matls, true);
    
    if (flags->d_doMechanics) {
      scheduleComputeDeformationGradient(d_subsched, d_perproc_patches, matls, true);
      scheduleComputeStressTensor(       d_subsched, d_perproc_patches, matls, true);
      scheduleFormStiffnessMatrix(       d_subsched, d_perproc_patches, matls);
      scheduleComputeInternalForce(      d_subsched, d_perproc_patches, matls);
      scheduleFormQ(                     d_subsched, d_perproc_patches, matls);
      scheduleSolveForDuCG(              d_subsched, d_perproc_patches, matls);
    }
    
    scheduleGetDisplacementIncrement(d_subsched,        d_perproc_patches, matls);
    scheduleUpdateGridKinematics(    d_subsched,        d_perproc_patches, matls);
    scheduleCheckConvergence(        d_subsched, level, d_perproc_patches, matls);
    
    d_subsched->compile();
    d_recompileSubsched = false;  
  }

  int count = 0;
  bool dispInc = false;
  bool dispIncQ = false;
  sum_vartype dispIncQNorm, dispIncNorm, dispIncQNorm0, dispIncNormMax;

  // Get all of the required particle data that is in the old_dw and put it 
  // in the subscheduler's  new_dw.  Then once dw is advanced, subscheduler
  // will be pulling data out of the old_dw.
  auto subsched_parent_old_dw = d_subsched->get_dw(0);
  auto subsched_new_dw = d_subsched->get_dw(3);
  //std::cout << "subsched: parent_old_dw = " << subsched_parent_old_dw << "\n";
  //std::cout << "subsched: new_dw = " << subsched_new_dw << "\n";

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::iterate-----------------------");

    for (auto m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();
      ParticleSubset* pset = subsched_parent_old_dw->getParticleSubset(matID, patch);

      delt_vartype dt;
      old_dw->get(dt, d_sharedState->get_delt_label(), patch->getLevel());

      subsched_new_dw->put(sum_vartype(0.0), lb->dispIncQNorm0);
      subsched_new_dw->put(sum_vartype(0.0), lb->dispIncNormMax);

      // New data to be stored in the subscheduler
      NCVariable<Vector> dispNew, dispNew_sub;
      new_dw->getModifiable(dispNew, lb->dispNewLabel, matID, patch);
      subsched_new_dw->allocateAndPut(dispNew_sub, lb->dispNewLabel, matID, patch);
      dispNew_sub.copyData(dispNew);

      if (!flags->d_doGridReset) {
        NCVariable<Vector> gDisplacement, gDisplacement_new;
        new_dw->getModifiable(gDisplacement, lb->gDisplacementLabel, matID, patch);
        subsched_new_dw->allocateAndPut(gDisplacement_new, lb->gDisplacementLabel,
                                        matID, patch);
        gDisplacement_new.copyData(gDisplacement);
      }

      subsched_new_dw->saveParticleSubset(pset, matID, patch);

      // These variables are ultimately retrieved from the subschedulers
      // old datawarehouse after the advancement of the data warehouse.
      double new_dt;
      new_dt = dt;
      subsched_new_dw->put(delt_vartype(new_dt), d_sharedState->get_delt_label());
    }
  }

  subsched_new_dw->finalize();
  d_subsched->advanceDataWarehouse(grid);
  d_subsched->setInitTimestep(false);

  d_numIterations = 0;
  while (!(dispInc && dispIncQ)) {
    proc0cout << "    Beginning Iteration = " << count << "\n";
    
    auto subsched_old_dw = d_subsched->get_dw(2);
    subsched_new_dw = d_subsched->get_dw(3);
    
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

    d_subsched->execute();  // THIS ACTUALLY GETS THE WORK DONE
    subsched_new_dw->get(dispIncNorm,    lb->dispIncNorm);
    subsched_new_dw->get(dispIncQNorm,   lb->dispIncQNorm); 
    subsched_new_dw->get(dispIncNormMax, lb->dispIncNormMax);
    subsched_new_dw->get(dispIncQNorm0,  lb->dispIncQNorm0);

    double frac_Norm  = dispIncNorm/(dispIncNormMax + 1.e-100);
    double frac_QNorm = dispIncQNorm/(dispIncQNorm0 + 1.e-100);

    if (UintahParallelComponent::d_myworld->myrank() == 0) {
      std::cerr << "  dispIncNorm/dispIncNormMax = " << frac_Norm << "\n";
      std::cerr << "  dispIncQNorm/dispIncQNorm0 = "<< frac_QNorm << "\n";
    }
    if ( (frac_Norm  <= flags->d_convCritDisp) || 
         (dispIncNormMax <= flags->d_convCritDisp) ) {
      dispInc = true;
    }  
    if ( (frac_QNorm <= flags->d_convCritEnergy) || 
         (dispIncQNorm0 <= flags->d_convCritEnergy) ) {
      dispIncQ = true;
    }
    
    // Check to see if the residual is likely a nan, if so, we'll restart.
    bool restart_nan = false;
    if ((std::isnan(dispIncQNorm/dispIncQNorm0) ||
         std::isnan(dispIncNorm/dispIncNormMax))
         && dispIncQNorm0!=0.) {
      restart_nan = true;
      if (UintahParallelComponent::d_myworld->myrank() == 0) {
        std::cerr << "Restarting due to a nan residual\n";
      }
    }

    bool restart_neg_residual = false;
    if (dispIncQNorm/(dispIncQNorm0 + 1e-100) < 0. ||
        dispIncNorm/(dispIncNormMax+1e-100) < 0.) {
      restart_neg_residual = true;
      if (UintahParallelComponent::d_myworld->myrank() == 0) {
        std::cerr << "Restarting due to a negative residual\n";
      }
    }

    bool restart_num_iters = false;
    if (count > flags->d_maxNumIterations) {
      restart_num_iters = true;
      if (UintahParallelComponent::d_myworld->myrank() == 0) {
        std::cerr << "Restarting due to exceeding max number of iterations\n";
      }
    }

    if (restart_nan || restart_neg_residual || restart_num_iters) {
      new_dw->abortTimestep();
      new_dw->restartTimestep();
      return;
    }

    d_subsched->advanceDataWarehouse(grid);

  } // endwhile

  d_numIterations = count;

  // Move the particle data from subscheduler to scheduler.
  auto subsched_old_dw = d_subsched->get_dw(2);
  for (int p = 0; p < patches->size();p++) {
    const Patch* patch = patches->get(p);
    if (cout_doing.active()) {
      cout_doing <<"  Getting the recursive data on patch " << patch->getID()
                 <<"\t\t\t IMPM"<< "\n" << "\n";
    }

    for(int m = 0; m < d_sharedState->getNumMPMMatls(); m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      // Needed in computeAcceleration 
      constNCVariable<Vector> gVelocity, dispNew, gInternalForce;
      subsched_old_dw->get(gVelocity, lb->gVelocityLabel, matID, patch, gnone, 0);
      subsched_old_dw->get(dispNew,   lb->dispNewLabel,   matID, patch, gnone, 0);
      if (flags->d_doMechanics) {
        subsched_old_dw->get(gInternalForce, lb->gInternalForceLabel, 
                             matID, patch, gnone, 0);
      }

      NCVariable<Vector> gVelocity_new, dispNew_new, gInternalForce_new;
      new_dw->getModifiable(gVelocity_new,      lb->gVelocityLabel,      matID, patch);
      new_dw->getModifiable(dispNew_new,        lb->dispNewLabel,        matID, patch);
      new_dw->getModifiable(gInternalForce_new, lb->gInternalForceLabel, matID, patch);
      gVelocity_new.copyData(gVelocity);
      dispNew_new.copyData(dispNew);
      if (flags->d_doMechanics) {
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
  printSchedule(patches, cout_doing, "ImpMPM::scheduleComputeDeformationGradient");
  Task* t = scinew Task("ImpMPM::computeDeformationGradient",
                        this, &ImpMPM::computeDeformationGradient, recursion);

  //std::cout << "OldDW = " << Task::OldDW << __FILE__ << __LINE__ << "\n";
  //std::cout << "ParentOldDW = " << Task::ParentOldDW << __FILE__ << __LINE__ << "\n";
  //t->requires(Task::ParentOldDW, d_sharedState->get_delt_label());
  int numMatls = d_sharedState->getNumMPMMatls();
  for(int m = 0; m < numMatls; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    d_defGradComputer->addComputesAndRequires(t, mpm_matl, patches, recursion, true);
  }

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeDeformationGradient(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw,
                                   bool recursion)
{
  printTask(patches, patches->get(0), cout_doing,
            "Doing ImpMPM::computeDeformationGradient with recursion");
  d_defGradComputer->computeDeformationGradient(patches, old_dw, new_dw, recursion);
}

void 
ImpMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls,
                                    bool recursion)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleComputeStressTensor");
  Task* t = scinew Task("ImpMPM::computeStressTensor",
                    this, &ImpMPM::computeStressTensorImplicit, recursion);

  t->requires(Task::ParentOldDW, d_sharedState->get_delt_label());
  int numMatls = d_sharedState->getNumMPMMatls();
  for(int m = 0; m < numMatls; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches, recursion, true);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::computeStressTensorImplicit(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset* ,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw,
                                    bool recursion)
{
  if (cout_doing.active()) {
    cout_doing <<"Doing computeStressTensor (wrapper) " <<"\t\t\t IMPM"<< "\n";
  }

  for(int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    ImplicitCM* cmi = dynamic_cast<ImplicitCM*>(cm);
    if (cmi) {
      cmi->computeStressTensorImplicit(patches, mpm_matl, old_dw, new_dw,
                               d_solver, recursion);
    }
  }
}

void
ImpMPM::scheduleComputeDeformationGradient(SchedulerP& sched,
                                           const PatchSet* patches,
                                           const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeDeformationGradientImplicit");
  Task* t = scinew Task("ImpMPM::computeDeformationGradientImplicit",
                        this, &ImpMPM::computeDeformationGradient);
  int numMatls = d_sharedState->getNumMPMMatls();
  for(int m = 0; m < numMatls; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    d_defGradComputer->addComputesAndRequires(t, mpm_matl, patches);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void
ImpMPM::computeDeformationGradient(const ProcessorGroup*,
                                   const PatchSubset* patches,
                                   const MaterialSubset* ,
                                   DataWarehouse* old_dw,
                                   DataWarehouse* new_dw)
{
  printTask(patches, patches->get(0), cout_doing,
            "Doing IMPM::computeDeformationGradient");
  d_defGradComputer->computeDeformationGradient(patches, old_dw, new_dw);
}

void 
ImpMPM::scheduleComputeStressTensor(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  int numMatls = d_sharedState->getNumMPMMatls();
  printSchedule(patches, cout_doing,"IMPM::scheduleComputeStressTensorImplicit");
  Task* t = scinew Task("ImpMPM::computeStressTensorImplicit",
                        this, &ImpMPM::computeStressTensorImplicit);

  for(int m = 0; m < numMatls; m++){
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->addComputesAndRequires(t, mpm_matl, patches);
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::computeStressTensorImplicit(const ProcessorGroup*,
                                    const PatchSubset* patches,
                                    const MaterialSubset* ,
                                    DataWarehouse* old_dw,
                                    DataWarehouse* new_dw)
{
  if (cout_doing.active()) {
    cout_doing <<"Doing computeStressTensorImplicit (wrapper)" <<"\t\t IMPM"<< "\n";
  }

  for(int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
    MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial(m);
    ConstitutiveModel* cm = mpm_matl->getConstitutiveModel();
    cm->computeStressTensorImplicit(patches, mpm_matl, old_dw, new_dw);
  }
}

void 
ImpMPM::scheduleFormStiffnessMatrix(SchedulerP& sched,
                                    const PatchSet* patches,
                                    const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFormStiffnessMatrix");
  Task* t = scinew Task("ImpMPM::formStiffnessMatrix",
                        this, &ImpMPM::formStiffnessMatrix);

  t->requires(Task::ParentOldDW, d_sharedState->get_delt_label());
  t->requires(Task::ParentNewDW, lb->gMassAllLabel, Ghost::None);

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
  if (!flags->d_dynamic) {
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::formStiffnessMatrix");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex, highIndex);

    bool firstTimeThrough=true;
    int numMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMatls; m++) {

      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!mpm_matl->getIsRigid() && firstTimeThrough) { 

        firstTimeThrough=false;
        int matID = mpm_matl->getDWIndex();
        d_solver->copyL2G(l2g, patch);
   
        DataWarehouse* parent_old_dw = new_dw->getOtherDataWarehouse(Task::ParentOldDW);
        DataWarehouse* parent_new_dw = new_dw->getOtherDataWarehouse(Task::ParentNewDW);

        delt_vartype dt;
        constNCVariable<double> gMass;

        parent_old_dw->get(dt,    d_sharedState->get_delt_label(), patch->getLevel());
        parent_new_dw->get(gMass, lb->gMassAllLabel, matID, patch, Ghost::None, 0);

        double v[1];
        int dof[3];
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          int l2g_node_num = l2g[node];
          v[0] = gMass[node]*(4./(dt*dt));
          dof[0] = l2g_node_num;
          dof[1] = l2g_node_num+1;
          dof[2] = l2g_node_num+2;

          d_solver->fillMatrix(1, &dof[0], 1, &dof[0], v);
          d_solver->fillMatrix(1, &dof[1], 1, &dof[1], v);
          d_solver->fillMatrix(1, &dof[2], 1, &dof[2], v);
        } // end node iterator

      } // endif
    } // endfor matls
  }
}

void 
ImpMPM::scheduleComputeInternalForce(SchedulerP& sched,
                                     const PatchSet* patches,
                                     const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleComputeInternalForce");
  Task* t = scinew Task("ImpMPM::computeInternalForce",
                         this, &ImpMPM::computeInternalForce);

  t->requires(Task::ParentOldDW, lb->pXLabel,                Ghost::AroundNodes, 1);
  t->requires(Task::ParentOldDW, lb->pSizeLabel,             Ghost::AroundNodes, 1);
  t->requires(Task::NewDW,       lb->pDefGradLabel_preReloc, Ghost::AroundNodes, 1);
  t->requires(Task::NewDW,       lb->pStressLabel_preReloc,  Ghost::AroundNodes, 1);
  t->requires(Task::NewDW,       lb->pVolumeLabel_preReloc,  Ghost::AroundNodes, 1);

  t->computes(lb->gInternalForceLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::computeInternalForce(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* ,
                             DataWarehouse* /*old_dw*/,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing ImpMPM::computeInternalForce");
    
    auto interpolator = flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0/dx.x();
    oodx[1] = 1.0/dx.y();
    oodx[2] = 1.0/dx.z();
    
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    int n8or27 = flags->d_8or27;

    NCVectorArray       gInternalForce(numMPMMatls);
    NCVariable<Vector>  gInternalForce_sum;
    new_dw->allocateTemporary(gInternalForce_sum, patch, Ghost::None, 0);
    gInternalForce_sum.initialize(Vector(0,0,0));

    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      new_dw->allocateAndPut(gInternalForce[m], lb->gInternalForceLabel, matID, patch);
      gInternalForce[m].initialize(Vector(0,0,0));

      if (!mpm_matl->getIsRigid()) {

        DataWarehouse* parent_old_dw = new_dw->getOtherDataWarehouse(Task::ParentOldDW);
        ParticleSubset* pset = parent_old_dw->getParticleSubset(matID, patch,
                                                                Ghost::AroundNodes, 1,
                                                                lb->pXLabel);

        constParticleVariable<Point>   pX;
        constParticleVariable<double>  pVolume;
        constParticleVariable<Matrix3> pStress;
        constParticleVariable<Matrix3> pSize;
        constParticleVariable<Matrix3> pDefGrad;

        parent_old_dw->get(pX,    lb->pXLabel,                pset);
        parent_old_dw->get(pSize, lb->pSizeLabel,             pset);
        new_dw->get(pVolume,      lb->pVolumeLabel_preReloc,  pset);
        new_dw->get(pDefGrad,     lb->pDefGradLabel_preReloc, pset);
        new_dw->get(pStress,      lb->pStressLabel_preReloc,  pset);

        Matrix3 pStressVol;

        for (auto idx : *pset) {

          interpolator->findCellAndShapeDerivatives(pX[idx], ni, d_S, pSize[idx],
                                                    pDefGrad[idx]);
          pStressVol  = pStress[idx]*pVolume[idx];
          for (int k = 0; k < n8or27; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
              Vector div(d_S[k].x()*oodx[0], d_S[k].y()*oodx[1],
                                             d_S[k].z()*oodx[2]);
              gInternalForce[m][node] -= (div * pStress[idx])  * pVolume[idx];
            }
          }
        }

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gInternalForce_sum[node] += gInternalForce[m][node];
        }
      }  // if matl isn't rigid
    }  // matls

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gInternalForce[m][node] = gInternalForce_sum[node];
        }
      }
    }  // matls
  } // patches
}

void 
ImpMPM::scheduleFormQ(SchedulerP& sched,
                      const PatchSet* patches,
                      const MaterialSet* matls)
{
  printSchedule(patches,cout_doing,"IMPM::scheduleFormQ");
  Task* t = scinew Task("ImpMPM::formQ", this, 
                        &ImpMPM::formQ);

  Ghost::GhostType  gnone = Ghost::None;

  t->requires(Task::ParentOldDW, d_sharedState->get_delt_label());
  t->requires(Task::ParentNewDW, lb->gMassAllLabel,       gnone, 0);
  t->requires(Task::ParentNewDW, lb->gBodyForceLabel,     gnone, 0);
  t->requires(Task::ParentNewDW, lb->gExternalForceLabel, gnone, 0);
  t->requires(Task::ParentNewDW, lb->gVelocityOldLabel,   gnone, 0);
  t->requires(Task::ParentNewDW, lb->gAccelerationLabel,  gnone, 0);
  t->requires(Task::OldDW,       lb->dispNewLabel,        gnone, 0);
  t->requires(Task::NewDW,       lb->gInternalForceLabel, gnone, 0);

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

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g,patch);

    bool firstTimeThrough = true;
    int numMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!mpm_matl->getIsRigid() && firstTimeThrough) {
        firstTimeThrough=false;
        int matID = mpm_matl->getDWIndex();

        DataWarehouse* parent_new_dw = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
        DataWarehouse* parent_old_dw = new_dw->getOtherDataWarehouse(Task::ParentOldDW);

        delt_vartype dt;
        constNCVariable<double> gMass;
        constNCVariable<Vector> gBodyForce, gExternalForce, gInternalForce;
        constNCVariable<Vector> dispNew, gVelocity, gAcceleration;

        parent_old_dw->get(dt, d_sharedState->get_delt_label(), patch->getLevel());
        parent_new_dw->get(gBodyForce,     lb->gBodyForceLabel,     matID, patch, gnone, 0);
        parent_new_dw->get(gExternalForce, lb->gExternalForceLabel, matID, patch, gnone, 0);
        parent_new_dw->get(gVelocity,      lb->gVelocityOldLabel,   matID, patch, gnone, 0);
        parent_new_dw->get(gAcceleration,  lb->gAccelerationLabel,  matID, patch, gnone, 0);
        parent_new_dw->get(gMass,          lb->gMassAllLabel,       matID, patch, gnone, 0);

        old_dw->get(       dispNew,        lb->dispNewLabel,        matID, patch, gnone, 0);
        new_dw->get(       gInternalForce, lb->gInternalForceLabel, matID, patch, gnone, 0);

        double fodts = 4./(dt*dt);
        double fodt = 4./dt;

        double Q = 0.;

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          int l2g_node_num = l2g[node];

          Vector force = gExternalForce[node] + gInternalForce[node] + gBodyForce[node];
          double v[3];
          v[0] = force.x();
          v[1] = force.y();
          v[2] = force.z();

          // temp2 = M*a^(k-1)(t+dt)
          if (flags->d_dynamic) {
            Vector damping = (dispNew[node]*fodts - gVelocity[node]*fodt -
                              gAcceleration[node])*gMass[node];
            v[0] -= damping.x();
            v[1] -= damping.y();
            v[2] -= damping.z();
          }
          d_solver->fillVector(l2g_node_num,   double(v[0]));
          d_solver->fillVector(l2g_node_num+1, double(v[1]));
          d_solver->fillVector(l2g_node_num+2, double(v[2]));
          Q += v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        }
        if (std::isnan(Q)) {
          std::cout << "RHS contains a nan, restarting timestep" << endl;
          new_dw->abortTimestep();
          new_dw->restartTimestep();
          return;
        }
      } // first time through non-rigid
    }  // matls
  }    // patches
}

void 
ImpMPM::scheduleSolveForDuCG(SchedulerP& sched,
                             const PatchSet* patches,
                             const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleSolveForDuCG");
  Task* t = scinew Task("ImpMPM::solveForDuCG", this, 
                        &ImpMPM::solveForDuCG);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::solveForDuCG(const ProcessorGroup* /*pg*/,
                     const PatchSubset* patches,
                     const MaterialSubset* ,
                     DataWarehouse*,
                     DataWarehouse* new_dw)

{
  if (cout_doing.active()) {
    for(int p = 0; p<patches->size(); p++) {
      const Patch* patch = patches->get(p);
      cout_doing <<"Doing solveForDuCG on patch " << patch->getID()
                  <<"\t\t\t\t IMPM"<< "\n";
    }
  }

  DataWarehouse* parent_new_dw = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
  bool tsr = parent_new_dw->timestepRestarted();

  if (!tsr) {  // if a tsr has already been called for don't do the solve
    d_solver->assembleVector();
    d_solver->removeFixedDOF();
    std::vector<double> guess;
    d_solver->solve(guess);   
  } else {
    std::cout << "skipping solve, timestep has already called for a restart" << endl;
  }
}

void 
ImpMPM::scheduleGetDisplacementIncrement(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleGetDisplacementIncrement");
  Task* t = scinew Task("ImpMPM::getDisplacementIncrement", this, 
                        &ImpMPM::getDisplacementIncrement);

  t->computes(lb->dispIncLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::getDisplacementIncrement(const ProcessorGroup* /*pg*/,
                                 const PatchSubset* patches,
                                 const MaterialSubset* ,
                                 DataWarehouse*,
                                 DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::getDisplacementIncrement");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g,patch);
 
    std::vector<double> x;
    int begin = d_solver->getSolution(x);
  
    int numMatls = d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      NCVariable<Vector> dispInc;
      new_dw->allocateAndPut(dispInc, lb->dispIncLabel, matID, patch);
      dispInc.initialize(Vector(0.));

      if (flags->d_doMechanics) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++){
          IntVector node = *iter;
          int l2g_node_num = l2g[node] - begin;
          dispInc[node] = Vector(x[l2g_node_num], x[l2g_node_num+1], x[l2g_node_num+2]);
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
  printSchedule(patches, cout_doing,"IMPM::scheduleUpdateGridKinematics");
  Task* t = scinew Task("ImpMPM::updateGridKinematics", this, 
                        &ImpMPM::updateGridKinematics);

  t->requires(Task::ParentOldDW,   d_sharedState->get_delt_label() );
  t->requires(Task::ParentNewDW,   lb->gVelocityOldLabel,    Ghost::None, 0);
  t->requires(Task::ParentNewDW,   lb->gContactLabel,        Ghost::None, 0);
  t->requires(Task::OldDW,         lb->dispNewLabel,         Ghost::None, 0);
  t->requires(Task::NewDW,         lb->dispIncLabel,         Ghost::None, 0);
  if (!flags->d_doGridReset) {
    t->requires(Task::ParentOldDW, lb->gDisplacementLabel,   Ghost::None, 0);
    t->computes(lb->gDisplacementLabel);
  }
  t->computes(lb->dispNewLabel);
  t->computes(lb->gVelocityLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::updateGridKinematics(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset* ,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{
  Ghost::GhostType gnone = Ghost::None;

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::updateGridKinematics");

    int matID_rigid =-99;
    int numMatls = d_sharedState->getNumMPMMatls();
    for (int m = 0; m < numMatls; m++) {
       MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
       if (mpm_matl->getIsRigid()) {
         matID_rigid = mpm_matl->getDWIndex();
       }
    }

    constNCVariable<Vector> gVelocity_rigid;
    if (d_rigid_body) {
      DataWarehouse* parent_new_dw = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      parent_new_dw->get(gVelocity_rigid, lb->gVelocityOldLabel, matID_rigid, patch,
                         gnone, 0);
    }

    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      delt_vartype            dt;
      constNCVariable<int>    gContact;
      constNCVariable<Vector> dispInc, dispNew_old, gVelocity_old;
      NCVariable<Vector>      dispNew, gVelocity;

      DataWarehouse* parent_new_dw = new_dw->getOtherDataWarehouse(Task::ParentNewDW);
      DataWarehouse* parent_old_dw = new_dw->getOtherDataWarehouse(Task::ParentOldDW);

      parent_old_dw->get(dt,            d_sharedState->get_delt_label(), patch->getLevel());
      parent_new_dw->get(gVelocity_old, lb->gVelocityOldLabel, matID, patch, gnone, 0);
      parent_new_dw->get(gContact,      lb->gContactLabel,     matID, patch, gnone, 0);
      old_dw->get(dispNew_old,          lb->dispNewLabel,      matID, patch, gnone, 0);
      new_dw->get(dispInc,              lb->dispIncLabel,      matID, patch, gnone, 0);

      new_dw->allocateAndPut(dispNew,   lb->dispNewLabel,      matID, patch);
      new_dw->allocateAndPut(gVelocity, lb->gVelocityLabel,    matID, patch);
      dispNew.copyData(dispNew_old);

      double oneifdyn = 0.;
      if (flags->d_dynamic) {
        oneifdyn = 1.;
      }

      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          dispNew[node] += dispInc[node];
          gVelocity[node] = dispNew[node]*(2./dt) - oneifdyn*gVelocity_old[node];
        }
      }

      if (d_rigid_body) {  // overwrite some of the values computed above
        for (auto iter = patch->getNodeIterator();!iter.done();iter++){
          IntVector node = *iter;
          if (gContact[node]==2) {
            dispNew[node] = Vector((1.-d_contact_dirs.x())*dispNew[node].x() +
                                   d_contact_dirs.x()*gVelocity_rigid[node].x()*dt,
                                   (1.-d_contact_dirs.y())*dispNew[node].y() +
                                   d_contact_dirs.y()*gVelocity_rigid[node].y()*dt,
                                   (1.-d_contact_dirs.z())*dispNew[node].z() +
                                   d_contact_dirs.z()*gVelocity_rigid[node].z()*dt);

            gVelocity[node] = dispNew[node]*(2./dt) - oneifdyn*gVelocity_old[node];
          } // if contact == 2
        } // for
      } // if d_rigid_body

      if (!flags->d_doGridReset) {
        constNCVariable<Vector> gDisplacement_old;
        NCVariable<Vector>      gDisplacement;
        parent_old_dw->get(gDisplacement_old, lb->gDisplacementLabel, matID, patch, gnone, 0);
        new_dw->allocateAndPut(gDisplacement, lb->gDisplacementLabel, matID, patch);
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gDisplacement[node] = gDisplacement_old[node] + dispNew[node];
        }
      }
    }   // matls
  }
}

void 
ImpMPM::scheduleCheckConvergence(SchedulerP& sched, 
                                 const LevelP& /* level */,
                                 const PatchSet* patches,
                                 const MaterialSet* matls)
{
  printSchedule(patches, cout_doing,"IMPM::scheduleCheckConvergence");
  Task* t = scinew Task("ImpMPM::checkConvergence", this,
                        &ImpMPM::checkConvergence);

  t->requires(Task::OldDW, lb->dispIncQNorm0);
  t->requires(Task::OldDW, lb->dispIncNormMax);
  t->requires(Task::NewDW, lb->dispIncLabel, Ghost::None, 0);

  t->computes(lb->dispIncNormMax);
  t->computes(lb->dispIncQNorm0);
  t->computes(lb->dispIncNorm);
  t->computes(lb->dispIncQNorm);

  t->setType(Task::OncePerProc);
  sched->addTask(t,patches,matls);
}

void 
ImpMPM::checkConvergence(const ProcessorGroup*,
                         const PatchSubset* patches,
                         const MaterialSubset* ,
                         DataWarehouse* old_dw,
                         DataWarehouse* new_dw)
{
  int global_offset = 0;
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::checkConvergence");

    auto lohiNodes = Util::getPatchLoHiNodes(patch, flags->d_8or27);
    auto lowIndex = lohiNodes.first;
    auto highIndex = lohiNodes.second;
    Array3<int> l2g(lowIndex,highIndex);
    d_solver->copyL2G(l2g,patch);

    int matID = 0;
    constNCVariable<Vector> dispInc;
    new_dw->get(dispInc, lb->dispIncLabel, matID, patch, Ghost::None, 0);
    
    double dispIncNorm  = 0.;
    double dispIncQNorm = 0.;
    std::vector<double> getQ;

    d_solver->getRHS(getQ);
    if (p == 0) {
      global_offset=l2g[lowIndex]; 
    }

    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      int l2g_node_num = l2g[node] - global_offset;
      dispIncNorm += Dot(dispInc[node], dispInc[node]);
      dispIncQNorm += dispInc[node].x()*getQ[l2g_node_num] +
                      dispInc[node].y()*getQ[l2g_node_num+1] + 
                      dispInc[node].z()*getQ[l2g_node_num+2];
    }

    // We are computing both dispIncQNorm0 and dispIncNormMax (max residuals)
    // We are computing both dispIncQNorm and dispIncNorm (current residuals)
    double dispIncQNorm0, dispIncNormMax;
    sum_vartype dispIncQNorm0_var, dispIncNormMax_var;
    old_dw->get(dispIncQNorm0_var,  lb->dispIncQNorm0);
    old_dw->get(dispIncNormMax_var, lb->dispIncNormMax);

    dispIncQNorm0 = dispIncQNorm0_var;
    dispIncNormMax = dispIncNormMax_var;

    bool first_iteration=false;
    if (compare(dispIncQNorm0,0.)) {
      first_iteration = true;
      dispIncQNorm0 = dispIncQNorm;
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
      dispIncQNorm0 /= ((double) numPatches);
      if (dispIncNormMax != dispIncNorm){
        dispIncNormMax /= ((double) numPatches);
      }
    }

    new_dw->put(sum_vartype(dispIncNorm),   lb->dispIncNorm);
    new_dw->put(sum_vartype(dispIncQNorm),  lb->dispIncQNorm);
    new_dw->put(sum_vartype(dispIncNormMax),lb->dispIncNormMax);
    new_dw->put(sum_vartype(dispIncQNorm0), lb->dispIncQNorm0);
  }  // End of loop over patches
}

void 
ImpMPM::scheduleUpdateTotalDisplacement(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  if (!flags->d_doGridReset) {
    printSchedule(patches, cout_doing, "IMPM::scheduleUpdateTotalDisplacement");
    Task* t = scinew Task("ImpMPM::updateTotalDisplacement",
                          this, &ImpMPM::updateTotalDisplacement);

    t->requires(Task::OldDW, lb->gDisplacementLabel, Ghost::None);
    t->requires(Task::NewDW, lb->dispNewLabel,       Ghost::None);
    t->modifies(lb->gDisplacementLabel);

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
    printTask(patches, patch, cout_doing, "Doing ImpMPM::updateTotalDisplacement");

    for (int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      constNCVariable<Vector> dispNew, gDisplacement_old;
      NCVariable<Vector>      gDisplacement;
      old_dw->get(gDisplacement_old,       lb->gDisplacementLabel, matID, patch, gnone, 0);
      new_dw->getModifiable(gDisplacement, lb->gDisplacementLabel, matID, patch);
      new_dw->get(dispNew,                 lb->dispNewLabel,       matID, patch, gnone, 0);

      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
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
  printSchedule(patches, cout_doing,"IMPM::scheduleComputeAcceleration");
  Task* t = scinew Task("ImpMPM::computeAcceleration",
                        this, &ImpMPM::computeAcceleration);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );
  t->requires(Task::NewDW, lb->gVelocityOldLabel,          Ghost::None);
  t->requires(Task::NewDW, lb->dispNewLabel,               Ghost::None);
  if (!flags->d_doGridReset) {
    t->requires(Task::OldDW, lb->gDisplacementLabel,       Ghost::None);
    t->modifies(lb->gDisplacementLabel);
  }
  t->modifies(lb->gAccelerationLabel);

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
  if (!flags->d_dynamic) {
    return;
  }
  Ghost::GhostType gnone = Ghost::None;

  delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::computeAcceleration");

    for (int m = 0; m < d_sharedState->getNumMPMMatls(); m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      constNCVariable<Vector> gVelocity, dispNew;
      NCVariable<Vector>      gAcceleration;

      new_dw->get(dispNew,                 lb->dispNewLabel,       matID, patch, gnone, 0);
      new_dw->get(gVelocity,               lb->gVelocityOldLabel,  matID, patch, gnone, 0);
      new_dw->getModifiable(gAcceleration, lb->gAccelerationLabel, matID, patch);

      double fodts = 4./(delT*delT);
      double fodt = 4./(delT);

      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        gAcceleration[node] = dispNew[node]*fodts - gVelocity[node]*fodt - gAcceleration[node];
      }
    }
  }
}

void 
ImpMPM::scheduleInterpolateToParticlesAndUpdate(SchedulerP& sched,
                                                const PatchSet* patches,
                                                const MaterialSet* matls)
{
  printSchedule(patches, cout_doing, "IMPM::scheduleInterpolateToParticlesAndUpdate");
  Task* t = scinew Task("ImpMPM::interpolateToParticlesAndUpdate",
                        this, &ImpMPM::interpolateToParticlesAndUpdate);

  t->requires(Task::OldDW, d_sharedState->get_delt_label() );

  t->requires(Task::OldDW, lb->pXLabel,                Ghost::None);
  t->requires(Task::OldDW, lb->pMassLabel,             Ghost::None);
  t->requires(Task::OldDW, lb->pParticleIDLabel,       Ghost::None);
  t->requires(Task::OldDW, lb->pVelocityLabel,         Ghost::None);
  t->requires(Task::OldDW, lb->pAccelerationLabel,     Ghost::None);
  t->requires(Task::OldDW, lb->pTemperatureLabel,      Ghost::None);
  t->requires(Task::OldDW, lb->pTempPreviousLabel,     Ghost::None);
  t->requires(Task::OldDW, lb->pDispLabel,             Ghost::None);
  t->requires(Task::OldDW, lb->pSizeLabel,             Ghost::None);

  t->requires(Task::NewDW, lb->pVolumeLabel_preReloc,  Ghost::None);
  t->requires(Task::NewDW, lb->pDefGradLabel_preReloc, Ghost::None);

  t->requires(Task::NewDW, lb->gAccelerationLabel,     Ghost::AroundCells,1);
  t->requires(Task::NewDW, lb->dispNewLabel,           Ghost::AroundCells,1);
  t->requires(Task::NewDW, lb->gTemperatureRateLabel, one_matl, Ghost::AroundCells, 1);

  t->computes(lb->pVelocityLabel_preReloc);
  t->computes(lb->pAccelerationLabel_preReloc);
  t->computes(lb->pXLabel_preReloc);
  t->computes(lb->pXXLabel);
  t->computes(lb->pParticleIDLabel_preReloc);
  t->computes(lb->pMassLabel_preReloc);
  t->computes(lb->pTemperatureLabel_preReloc);
  t->computes(lb->pDispLabel_preReloc);
  t->computes(lb->pSizeLabel_preReloc);
  t->computes(lb->pTempPreviousLabel_preReloc);

  if (flags->d_artificialViscosity) {
    t->requires(Task::OldDW, lb->p_qLabel,               Ghost::None);
    t->computes(lb->p_qLabel_preReloc);
  }

  t->computes(lb->KineticEnergyLabel);
  t->computes(lb->CenterOfMassPositionLabel);
  t->computes(lb->TotalMomentumLabel);
  t->computes(lb->ThermalEnergyLabel);

  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::interpolateToParticlesAndUpdate(const ProcessorGroup*,
                                        const PatchSubset* patches,
                                        const MaterialSubset* ,
                                        DataWarehouse* old_dw,
                                        DataWarehouse* new_dw)
{
  Ghost::GhostType gac = Ghost::AroundCells;

  delt_vartype delT;
  old_dw->get(delT, d_sharedState->get_delt_label(), getLevel(patches) );

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::interpolateToParticlesAndUpdate");

    auto interpolator = flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());
 
    // Performs the interpolation from the cell vertices of the grid
    // acceleration and displacement to the particles to update their
    // velocity and position respectively
    Vector disp(0.0,0.0,0.0);
    Vector acc(0.0,0.0,0.0);
  
    // DON'T MOVE THESE!!!
    Vector CMX(0.0,0.0,0.0);
    Vector totalMom(0.0,0.0,0.0);
    double ke = 0;
    double thermal_energy = 0.0;
    int numMPMMatls = d_sharedState->getNumMPMMatls();
    int n8or27 = flags->d_8or27;

    double move_particles = 1.;
    if (!flags->d_doGridReset) {
      move_particles = 0.;
    }

    constNCVariable<double> gTemperatureRate;
    new_dw->get(gTemperatureRate, lb->gTemperatureRateLabel, 0, patch, gac, 1);

    for (int m = 0; m < numMPMMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();
      double Cp = mpm_matl->getSpecificHeat();

      constParticleVariable<Point>   pX;
      constParticleVariable<double>  pMass, pVolume_new, pTemp, pq;
      constParticleVariable<Vector>  pVelocity, pAcceleration, pDisp;
      constParticleVariable<Matrix3> pSize, pDefGrad;

      ParticleVariable<double>  pMass_new, pTemp_new, pq_new;
      ParticleVariable<double>  pTempPre_new;
      ParticleVariable<Point>   pX_new, pXx;
      ParticleVariable<Vector>  pVelocity_new, pAcceleration_new, pDisp_new;
      ParticleVariable<Matrix3> pSize_new;
 
      constNCVariable<Vector> dispNew, gAcceleration;

      ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);
      ParticleSubset* delete_particles = scinew ParticleSubset(0, matID, patch);
    
      old_dw->get(pX,                    lb->pXLabel,                    pset);
      old_dw->get(pMass,                 lb->pMassLabel,                 pset);
      old_dw->get(pVelocity,             lb->pVelocityLabel,             pset);
      old_dw->get(pAcceleration,         lb->pAccelerationLabel,         pset);
      old_dw->get(pTemp,                 lb->pTemperatureLabel,          pset);
      old_dw->get(pDisp,                 lb->pDispLabel,                 pset);
      old_dw->get(pSize,                 lb->pSizeLabel,                 pset);
      new_dw->get(pVolume_new,           lb->pVolumeLabel_preReloc,      pset);
      new_dw->get(pDefGrad,              lb->pDefGradLabel_preReloc,     pset);

      new_dw->allocateAndPut(pMass_new,         lb->pMassLabel_preReloc,         pset);
      new_dw->allocateAndPut(pTemp_new,         lb->pTemperatureLabel_preReloc,  pset);
      new_dw->allocateAndPut(pTempPre_new,      lb->pTempPreviousLabel_preReloc, pset);
      new_dw->allocateAndPut(pVelocity_new,     lb->pVelocityLabel_preReloc,     pset);
      new_dw->allocateAndPut(pAcceleration_new, lb->pAccelerationLabel_preReloc, pset);
      new_dw->allocateAndPut(pX_new,            lb->pXLabel_preReloc,            pset);
      new_dw->allocateAndPut(pXx,               lb->pXXLabel,                    pset);
      new_dw->allocateAndPut(pDisp_new,         lb->pDispLabel_preReloc,         pset);
      new_dw->allocateAndPut(pSize_new,   lb->pSizeLabel_preReloc,        pset);

      new_dw->get(dispNew,        lb->dispNewLabel,       matID, patch, gac, 1);
      new_dw->get(gAcceleration,  lb->gAccelerationLabel, matID, patch, gac, 1);

      if (flags->d_artificialViscosity) {
        old_dw->get(pq,                lb->p_qLabel,          pset);
        new_dw->allocateAndPut(pq_new, lb->p_qLabel_preReloc, pset);
        pq_new.copyData(pq);
      }
      pSize_new.copyData(pSize);
      pTempPre_new.copyData(pTemp);

      for (auto idx : *pset) {

        interpolator->findCellAndWeightsAndShapeDerivatives(pX[idx], ni, S, d_S,
                                                            pSize[idx], pDefGrad[idx]);

        disp = Vector(0.0,0.0,0.0);
        acc = Vector(0.0,0.0,0.0);
        double tempRate = 0.;

        // Accumulate the contribution from each surrounding vertex
        for (int k = 0; k < n8or27; k++) {
          auto node = ni[k];
          disp      += dispNew[node]          * S[k];
          acc       += gAcceleration[node]    * S[k];
          tempRate  += gTemperatureRate[node] * S[k];
        }

        // Update the particle's position and velocity
        pX_new[idx]        = pX[idx] + disp*move_particles;
        pDisp_new[idx]     = pDisp[idx] + disp;
        pVelocity_new[idx] = pVelocity[idx] + (pAcceleration[idx]+acc)*(.5* delT);

        // pXx is only useful if we're not in normal grid resetting mode.
        pXx[idx]           = pX[idx] + pDisp_new[idx];

        pAcceleration_new[idx] = acc;
        pMass_new[idx]         = pMass[idx];
        pTemp_new[idx]         = pTemp[idx] + tempRate*delT;

        if (pMass_new[idx] <= 0.0) {
          delete_particles->addParticle(idx);
          pVelocity_new[idx] = Vector(0.,0.,0);
          pX_new[idx] = pX[idx];
        }

        thermal_energy += pTemp_new[idx] * pMass[idx] * Cp;
        // Thermal energy due to temperature flux (spatially varying part).
        //thermal_energy2 += potential_energy* pVolume_new[idx];
        ke += .5*pMass[idx] * pVelocity_new[idx].length2();
        CMX = CMX + (pX_new[idx] * pMass[idx]).asVector();
        totalMom += pVelocity_new[idx] * pMass[idx];
      }
 
      if (mpm_matl->getIsRigid()) {
        const double tcurr = d_sharedState->getElapsedTime(); 
        if (tcurr >= d_stop_time) {
          for (auto idx : *pset) { 
            pVelocity_new[idx] = d_vel_after_stop;
          }
        }
      }

      new_dw->deleteParticles(delete_particles);
      
      constParticleVariable<long64> pids;
      ParticleVariable<long64> pids_new;
      old_dw->get(pids, lb->pParticleIDLabel, pset);
      new_dw->allocateAndPut(pids_new, lb->pParticleIDLabel_preReloc, pset);
      pids_new.copyData(pids);
    }

    // DON'T MOVE THESE!!!
    new_dw->put(sum_vartype(ke),             lb->KineticEnergyLabel);
    new_dw->put(sumvec_vartype(CMX),         lb->CenterOfMassPositionLabel);
    new_dw->put(sumvec_vartype(totalMom),    lb->TotalMomentumLabel);
    new_dw->put(sum_vartype(thermal_energy), lb->ThermalEnergyLabel);
  }
}

void 
ImpMPM::scheduleInterpolateStressToGrid(SchedulerP& sched,
                                        const PatchSet* patches,
                                        const MaterialSet* matls)
{
  // This task is done for visualization only
  Task* t = scinew Task("ImpMPM::interpolateStressToGrid",
                        this, &ImpMPM::interpolateStressToGrid);

  t->requires(Task::OldDW, lb->pXLabel,               Ghost::AroundNodes,1);
  t->requires(Task::OldDW, lb->pSizeLabel,            Ghost::AroundNodes,1);
  t->requires(Task::NewDW, lb->pVolumeLabel_preReloc, Ghost::AroundNodes,1);
  t->requires(Task::NewDW, lb->pStressLabel_preReloc, Ghost::AroundNodes,1);
  t->requires(Task::NewDW, lb->gVolumeLabel,          Ghost::None);
  t->requires(Task::NewDW, lb->gVolumeLabel, d_sharedState->getAllInOneMatl(),
              Task::OutOfDomain, Ghost::None);
  t->requires(Task::OldDW, lb->pDefGradLabel,         Ghost::AroundNodes,1);

  t->modifies(lb->gInternalForceLabel);
  t->computes(lb->gStressForSavingLabel);

  if (!d_bndy_traction_faces.empty()) {
    for (auto face : d_bndy_traction_faces) {
      t->computes(lb->BndyForceLabel[face]);       // node based
      t->computes(lb->BndyContactAreaLabel[face]); // node based
    }
  }
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}

void 
ImpMPM::interpolateStressToGrid(const ProcessorGroup*,
                                const PatchSubset* patches,
                                const MaterialSubset* ,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  // node based forces
  Vector bndyForce[6];
  double bndyArea[6];
  for (int iface = 0; iface < 6; iface++) {
    bndyForce[iface]  = Vector(0.);
    bndyArea [iface]  = 0.;
  }
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch, cout_doing, "Doing ImpMPM::interpolateStressToGrid");

    auto interpolator = flags->d_interpolator->clone(patch);
    std::vector<IntVector> ni(interpolator->size());
    std::vector<double> S(interpolator->size());
    std::vector<Vector> d_S(interpolator->size());

    int numMatls = d_sharedState->getNumMPMMatls();
    int n8or27 = flags->d_8or27;

    constNCVariable<double>   gVolume_sum;
    new_dw->get(gVolume_sum, lb->gVolumeLabel,
                d_sharedState->getAllInOneMatl()->get(0), patch, Ghost::None, 0);

    NCVariable<Vector>  gInternalForce_sum;
    NCVariable<Matrix3> gStress_sum;
    new_dw->allocateTemporary(gInternalForce_sum, patch, Ghost::None, 0);
    new_dw->allocateTemporary(gStress_sum,        patch, Ghost::None, 0);
    gInternalForce_sum.initialize(Vector(0,0,0));
    gStress_sum.initialize(Matrix3(0.));

    std::vector<constNCVariable<double> > gVolume(numMatls);
    NCMatrix3Array                        gStress(numMatls);
    NCVectorArray                         gInternalForce(numMatls);

    Vector dx = patch->dCell();
    double oodx[3];
    oodx[0] = 1.0/dx.x();
    oodx[1] = 1.0/dx.y();
    oodx[2] = 1.0/dx.z();

    for(int m = 0; m < numMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int matID = mpm_matl->getDWIndex();

      new_dw->get(gVolume[m],                  lb->gVolumeLabel, matID, patch, Ghost::None, 0);
      new_dw->getModifiable(gInternalForce[m], lb->gInternalForceLabel, matID, patch);
      new_dw->allocateAndPut(gStress[m], lb->gStressForSavingLabel, matID, patch);
      gInternalForce[m].initialize(Vector(0.));
      gStress[m].initialize(Matrix3(0.));

      if (!mpm_matl->getIsRigid()) {

        constParticleVariable<Point>   pX;
        constParticleVariable<double>  pVolume;
        constParticleVariable<Matrix3> pSize, pStress, pDefGrad;

        ParticleSubset* pset = old_dw->getParticleSubset(matID, patch,
                                                         Ghost::AroundNodes, 1,
                                                         lb->pXLabel);
        old_dw->get(pX,       lb->pXLabel,               pset);
        new_dw->get(pVolume,  lb->pVolumeLabel_preReloc, pset);
        old_dw->get(pSize,    lb->pSizeLabel,            pset);
        old_dw->get(pDefGrad, lb->pDefGradLabel,         pset);
        new_dw->get(pStress,  lb->pStressLabel_preReloc, pset);

        Matrix3 pStressVol;
        for (auto idx : *pset) {

          interpolator->findCellAndWeightsAndShapeDerivatives(pX[idx], ni, S, d_S,
                                                              pSize[idx], pDefGrad[idx]);
          pStressVol = pStress[idx] * pVolume[idx];
          for (int k = 0; k < n8or27; k++) {
            auto node = ni[k];
            if (patch->containsNode(node)) {
               gStress[m][node] += pStressVol * S[k];
               Vector div(d_S[k].x()*oodx[0], d_S[k].y()*oodx[1], d_S[k].z()*oodx[2]);           
               gInternalForce[m][node] -= div * pStressVol;
            }
          }
        }

        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gStress_sum[node] += (gStress[m][node]);
          gInternalForce_sum[node]+=gInternalForce[m][node];
        }
      }  // if matl isn't rigid
    }  // Loop over matls

    // gStress will be normalized by gVolume (same for all matls)
    for (int m = 0; m < numMatls; m++) {
      for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
        IntVector node = *iter;
        gStress[m][node] = gStress_sum[node]/(gVolume[m][node]+1.e-200);
      }
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!mpm_matl->getIsRigid()) {
        for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
          IntVector node = *iter;
          gInternalForce[m][node] = gInternalForce_sum[node];
        }
      }
    }  // Loop over matls

    // Fill in the value for the all in one material
    // gStress_sum will be normalized by gVolume_global
    for (auto iter = patch->getNodeIterator(); !iter.done(); iter++) {
      IntVector node = *iter;
      gStress_sum[node] /= (gVolume_sum[node]+1.e-200);
    }

    // save boundary forces before apply symmetry boundary condition.
    bool did_it_already = false;
    for (int m = 0; m < numMatls; m++) {
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      if (!did_it_already && !mpm_matl->getIsRigid()) {
        did_it_already = true;
        for (auto face : d_bndy_traction_faces) {
                                                                                
          // Check if the face is on an external boundary
          if (patch->getBCType(face)==Patch::Neighbor) {
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
                IntVector node(i,j,k);
                // flip sign so that pushing on boundary gives positive force
                bndyForce[face] -= gInternalForce[m][node];
                                                                                
                double celldepth  = dx[face/2];
                bndyArea [face] += gVolume[m][node]/celldepth;
              }
            }
          }
        } // faces
      } // if
    }  // matls
  }

  // be careful only to put the fields that we have built
  // that way if the user asks to output a field that has not been built
  // it will fail early rather than just giving zeros.
  for (auto face : d_bndy_traction_faces) {
    new_dw->put(sumvec_vartype(bndyForce[face]), lb->BndyForceLabel[face]);
    new_dw->put(sum_vartype(bndyArea[face]), lb->BndyContactAreaLabel[face]);
  }
}

void ImpMPM::scheduleComputeHeatExchange(SchedulerP& sched,
                                         const PatchSet* patches,
                                         const MaterialSet* matls)
{
  /* computeHeatExchange
   *   in(G.MASS, G.TEMPERATURE, G.EXTERNAL_HEAT_RATE)
   *   operation(peform heat exchange which will cause each of
   *   velocity fields to exchange heat according to 
   *   the temperature differences)
   *   out(G.EXTERNAL_HEAT_RATE) */

  printSchedule(patches,cout_doing,"IMPM::scheduleComputeHeatExchange");
  Task* t = scinew Task("ThermalContact::computeHeatExchange",
                        thermalContactModel,
                        &ThermalContact::computeHeatExchange);

  thermalContactModel->addComputesAndRequires(t, patches, matls);
  t->setType(Task::OncePerProc);
  sched->addTask(t, patches, matls);
}



void ImpMPM::scheduleRefine(const PatchSet* patches,
                            SchedulerP& sched)
{
  printSchedule(patches,cout_doing,"ImpMPM::scheduleRefine");
  Task* t = scinew Task("ImpMPM::refine", this, &ImpMPM::refine);
                                                                                
  t->computes(lb->partCountLabel);
  t->computes(lb->pXLabel);
  t->computes(lb->pMassLabel);
  t->computes(lb->pVolumeLabel);
  t->computes(lb->pDispLabel);
  t->computes(lb->pVelocityLabel);
  t->computes(lb->pAccelerationLabel);
  t->computes(lb->pExternalForceLabel);
  t->computes(lb->pTemperatureLabel);
  t->computes(lb->pTempPreviousLabel);
  t->computes(lb->pSizeLabel);
  t->computes(lb->pParticleIDLabel);
  t->computes(lb->pDefGradLabel);
  t->computes(lb->pStressLabel);
  t->computes(lb->pCellNAPIDLabel);
  t->computes(d_sharedState->get_delt_label(),getLevel(patches));

  t->computes(lb->pExternalHeatFluxLabel);

  t->computes(lb->heatRate_CCLabel);
  if(!flags->d_doGridReset){
    t->computes(lb->gDisplacementLabel);
  }

  if (flags->d_useLoadCurves) {
    // Computes the load curve ID associated with each particle
    t->computes(lb->pLoadCurveIDLabel);
  }

  t->computes(lb->NC_CCweightLabel, one_matl);

  sched->addTask(t, patches, d_sharedState->allMPMMaterials());

  Level* level = const_cast<Level*>(getLevel(patches));
#if 1
  if (flags->d_useLoadCurves) {
    // Schedule the initialization of HeatFlux BCs per particle
    scheduleInitializeHeatFluxBCs(level, sched);
  }
#endif
}
                                                                                
void ImpMPM::scheduleRefineInterface(const LevelP& /*fineLevel*/,
                                        SchedulerP& /*scheduler*/,
                                        bool, bool)
{
  // do nothing for now
}
                                                                                
void ImpMPM::scheduleCoarsen(const LevelP& /*coarseLevel*/,
                                SchedulerP& /*sched*/)
{
  // do nothing for now
}
//______________________________________________________________________
// Schedule to mark flags for AMR regridding
void ImpMPM::scheduleErrorEstimate(const LevelP& coarseLevel,
                                      SchedulerP& sched)
{
  // main way is to count particles, but for now we only want particles on
  // the finest level.  Thus to schedule cells for regridding during the
  // execution, we'll coarsen the flagged cells (see coarsen).
  printSchedule(coarseLevel,cout_doing,"ImpMPM::scheduleErrorEstimate");     
                                                              
  if (cout_doing.active())
    cout_doing << "ImpMPM::scheduleErrorEstimate on level " << coarseLevel->getIndex() << '\n';
              
  // The simulation controller should not schedule it every time step
  Task* task = scinew Task("errorEstimate", this, &ImpMPM::errorEstimate);
                                                                                
  // if the finest level, compute flagged cells
  if (coarseLevel->getIndex() == coarseLevel->getGrid()->numLevels()-1) {
    task->requires(Task::NewDW, lb->pXLabel, Ghost::AroundCells, 0);
  }
  else {
    task->requires(Task::NewDW, d_sharedState->get_refineFlag_label(),
                   0, Task::FineLevel, d_sharedState->refineFlagMaterials(),
                   Task::NormalDomain, Ghost::None, 0);
  }
  task->modifies(d_sharedState->get_refineFlag_label(),      d_sharedState->refineFlagMaterials());
  task->modifies(d_sharedState->get_refinePatchFlag_label(), d_sharedState->refineFlagMaterials());
  sched->addTask(task, coarseLevel->eachPatch(), d_sharedState->allMPMMaterials());
                                                                                
}
//______________________________________________________________________
// Schedule to mark initial flags for AMR regridding
void ImpMPM::scheduleInitialErrorEstimate(const LevelP& coarseLevel,
                                             SchedulerP& sched)
{
  scheduleErrorEstimate(coarseLevel, sched);
}

void ImpMPM::scheduleSwitchTest(const LevelP& level, SchedulerP& sched)
{
  if (d_switchCriteria) {
    d_switchCriteria->scheduleSwitchTest(level,sched);
  }
}























void ImpMPM::setSharedState(SimulationStateP& ssp)
{
  d_sharedState = ssp;
}



double ImpMPM::recomputeTimestep(double current_dt)
{
  return current_dt*flags->d_delTDecreaseFactor;
}

void ImpMPM::initialErrorEstimate(const ProcessorGroup*,
                                  const PatchSubset* patches,
                                  const MaterialSubset* /*matls*/,
                                  DataWarehouse*,
                                  DataWarehouse* new_dw)
{
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing initialErrorEstimate");
                                                                                
    CCVariable<int> refineFlag;
    PerPatch<PatchFlagP> refinePatchFlag;
    new_dw->getModifiable(refineFlag, d_sharedState->get_refineFlag_label(),
                          0, patch);
    new_dw->get(refinePatchFlag, d_sharedState->get_refinePatchFlag_label(),
                0, patch);
                                                                                
    PatchFlag* refinePatch = refinePatchFlag.get().get_rep();
                                                                                
                                                                                
    for(int m = 0; m < d_sharedState->getNumMPMMatls(); m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int dwi = mpm_matl->getDWIndex();
      // Loop over particles
      ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
      constParticleVariable<Point> pX;
      new_dw->get(pX, lb->pXLabel, pset);
                                                                                
      for(ParticleSubset::iterator iter = pset->begin();
          iter != pset->end(); iter++){
        refineFlag[patch->getLevel()->getCellIndex(pX[*iter])] = true;
        refinePatch->set();
      }
    }
  }
}

void ImpMPM::errorEstimate(const ProcessorGroup* group,
                           const PatchSubset* patches,
                           const MaterialSubset* matls,
                           DataWarehouse* old_dw,
                           DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);
  if (level->getIndex() == level->getGrid()->numLevels()-1) {
    // on finest level, we do the same thing as initialErrorEstimate, so call it
    initialErrorEstimate(group, patches, matls, old_dw, new_dw);
  }
  else {
    // coarsen the errorflag.
    const Level* fineLevel = level->getFinerLevel().get_rep();
                                                                                
    for(int p=0;p<patches->size();p++){
      const Patch* coarsePatch = patches->get(p);
      
      printTask(patches, coarsePatch,cout_doing,
                "Doing IMPM::errorEstimate\t\t\t\t\t");
                                                                                
      CCVariable<int> refineFlag;
      PerPatch<PatchFlagP> refinePatchFlag;
                                                                                
      new_dw->getModifiable(refineFlag, d_sharedState->get_refineFlag_label(),
                            0, coarsePatch);
      new_dw->get(refinePatchFlag, d_sharedState->get_refinePatchFlag_label(),
                  0, coarsePatch);
                                                                                
      PatchFlag* refinePatch = refinePatchFlag.get().get_rep();
                                                                                
      Level::selectType finePatches;
      coarsePatch->getFineLevelPatches(finePatches);
                                                                                
      // coarsen the fineLevel flag
      for(int i=0;i<finePatches.size();i++){
        const Patch* finePatch = finePatches[i];
                                                                                
        IntVector cl, ch, fl, fh;
        getFineLevelRange(coarsePatch, finePatch, cl, ch, fl, fh);
        if (fh.x() <= fl.x() || fh.y() <= fl.y() || fh.z() <= fl.z()) {
          continue;
        }
        constCCVariable<int> fineErrorFlag;
        new_dw->getRegion(fineErrorFlag,
                          d_sharedState->get_refineFlag_label(), 0,
                          fineLevel,fl, fh, false);
                                                                                
        //__________________________________
        //if the fine level flag has been set
        // then set the corrsponding coarse level flag
        for(CellIterator iter(cl, ch); !iter.done(); iter++){
          IntVector fineStart(level->mapCellToFiner(*iter));
                                                                                
          for(CellIterator inside(IntVector(0,0,0),
               fineLevel->getRefinementRatio()); !inside.done(); inside++){
                                                                                
            if (fineErrorFlag[fineStart+*inside]) {
              refineFlag[*iter] = 1;
              refinePatch->set();
            }
          }
        }  // coarse patch iterator
      }  // fine patch loop
    } // coarse patch loop
  }
}

void ImpMPM::refine(const ProcessorGroup*,
                   const PatchSubset* patches,
                   const MaterialSubset* /*matls*/,
                   DataWarehouse*,
                   DataWarehouse* new_dw)
{
  // just create a particle subset if one doesn't exist
  for (int p = 0; p<patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch,cout_doing,"Doing refine");
                                                                                
    int numMPMMatls=d_sharedState->getNumMPMMatls();
    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = d_sharedState->getMPMMaterial( m );
      int dwi = mpm_matl->getDWIndex();
                                                                                
      if (cout_doing.active()) {
        cout_doing <<"Doing refine on patch "
                   << patch->getID() << " material # = " << dwi << endl;
      }
                                                                                
      // this is a new patch, so create empty particle variables.
      if (!new_dw->haveParticleSubset(dwi, patch)) {
        ParticleSubset* pset = new_dw->createParticleSubset(0, dwi, patch);
                                                                                
        // Create arrays for the particle data
        ParticleVariable<Point>  pX;
        ParticleVariable<double> pMass, pVolume, pTemperature;
        ParticleVariable<Vector> pVelocity, pExternalForce, pDisp;
        ParticleVariable<Matrix3> pSize;
        ParticleVariable<double> pTempPrev;
        ParticleVariable<int>    pLoadCurve;
        ParticleVariable<long64> pID;
        ParticleVariable<Matrix3> pDefGrad, pStress;
                                                                                
        new_dw->allocateAndPut(pX,             lb->pXLabel,             pset);
        new_dw->allocateAndPut(pMass,          lb->pMassLabel,          pset);
        new_dw->allocateAndPut(pVolume,        lb->pVolumeLabel,        pset);
        new_dw->allocateAndPut(pVelocity,      lb->pVelocityLabel,      pset);
        new_dw->allocateAndPut(pTemperature,   lb->pTemperatureLabel,   pset);
        new_dw->allocateAndPut(pTempPrev,      lb->pTempPreviousLabel,  pset);
        new_dw->allocateAndPut(pExternalForce, lb->pExternalForceLabel, pset);
        new_dw->allocateAndPut(pID,            lb->pParticleIDLabel,    pset);
        new_dw->allocateAndPut(pDisp,          lb->pDispLabel,          pset);
        if (flags->d_useLoadCurves){
          new_dw->allocateAndPut(pLoadCurve,   lb->pLoadCurveIDLabel,   pset);
        }
        new_dw->allocateAndPut(pSize,          lb->pSizeLabel,          pset);
                                                                                
        mpm_matl->getConstitutiveModel()->initializeCMData(patch,
                                                           mpm_matl,new_dw);
      }
    }
    //__________________________________
    // - Initialize NC_CCweight = 0.125
    // - Find the walls with symmetry BC and double NC_CCweight
    NCVariable<double> NC_CCweight;
    new_dw->allocateAndPut(NC_CCweight, lb->NC_CCweightLabel,    0, patch);
    NC_CCweight.initialize(0.125);
    for(Patch::FaceType face = Patch::startFace; face <= Patch::endFace;
        face=Patch::nextFace(face)){
      int mat_id = 0;
      if (patch->haveBC(face,mat_id,"symmetry","Symmetric")) {
        
        for(CellIterator iter = patch->getFaceIterator(face,Patch::FaceNodes);
                                                  !iter.done(); iter++) {
          NC_CCweight[*iter] = 2.0*NC_CCweight[*iter];
        }
      }
    }
  }
} // end refine()

