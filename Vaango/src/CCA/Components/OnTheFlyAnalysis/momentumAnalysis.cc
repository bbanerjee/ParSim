/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
 * Copyright (c) 2022-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/OnTheFlyAnalysis/FileInfoVar.h>
#include <CCA/Components/OnTheFlyAnalysis/momentumAnalysis.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/Variables/PerPatch.h>

#include <Core/Util/DOUT.hpp>
#include <Core/Util/DebugStream.h>
#include <Core/Util/FileUtils.h>

#include <cstdio>
#include <dirent.h>
#include <iostream>

using namespace Uintah;
using namespace std;
//______________________________________________________________________
//     T O D O
//
//  Create a vector of control volumes.  Remove the assumption that
//  the entire domain is the CV.
//  This assumes that the control volume is cubic and aligned with the grid.
//  The face centered velocities are used to compute the fluxes through
//  the control surface.
//  This assumes that the variables all come from the new_dw
//  This assumes the entire computational domain is being used as the control
//  volume!!!         <<<<<<<<<<<<,

Dout dout_OTF_MA("momentumAnalysis",
                 "OnTheFlyAnalysis",
                 "Task scheduling and execution.",
                 false);
Dout dbg_OTF_MA("momentumAnalysis_dbg",
                "OnTheFlyAnalysis",
                "Displays detailed debugging info.",
                false);
//______________________________________________________________________
//
momentumAnalysis::momentumAnalysis(const ProcessorGroup* myworld,
                                   const MaterialManagerP& materialManager,
                                   const ProblemSpecP& module_spec)
  : AnalysisModule(myworld, materialManager, module_spec)
{
  d_zeroPatch = 0;
  d_matlIndx  = -9;
  d_pressIndx = 0; // For ICE/MPMICE it's always 0.
  d_pressMatl = m_zeroMatl;

  labels = scinew MA_Labels();

  labels->lastCompTime =
    VarLabel::create("lastCompTime", max_vartype::getTypeDescription());
  labels->fileVarsStruct =
    VarLabel::create("FileInfo_MA", PerPatch<FileInfoP>::getTypeDescription());
  labels->totalCVMomentum =
    VarLabel::create("totalCVMomentum", sumvec_vartype::getTypeDescription());
  labels->convectMom_fluxes =
    VarLabel::create("convectMom_fluxes", sumvec_vartype::getTypeDescription());
  labels->viscousMom_fluxes =
    VarLabel::create("viscousMom_fluxes", sumvec_vartype::getTypeDescription());
  labels->pressForces =
    VarLabel::create("pressForces", sumvec_vartype::getTypeDescription());
}

//______________________________________________________________________
//
momentumAnalysis::~momentumAnalysis()
{
  DOUTR(dout_OTF_MA, " Doing: Destructor momentumAnalysis ");

  if (d_zeroPatch && d_zeroPatch->removeReference()) {
    delete d_zeroPatch;
  }
  if (d_matl_set && d_matl_set->removeReference()) {
    delete d_matl_set;
  }

  VarLabel::destroy(labels->lastCompTime);
  VarLabel::destroy(labels->fileVarsStruct);
  VarLabel::destroy(labels->totalCVMomentum);
  VarLabel::destroy(labels->convectMom_fluxes);
  VarLabel::destroy(labels->viscousMom_fluxes);
  VarLabel::destroy(labels->pressForces);

  delete labels;
}

//______________________________________________________________________
//
void
momentumAnalysis::problemSetup(
  const ProblemSpecP&,
  const ProblemSpecP&,
  GridP& grid,
  [[maybe_unused]] std::vector<std::vector<const VarLabel*>>& PState,
  [[maybe_unused]] std::vector<std::vector<const VarLabel*>>& PState_preReloc)
{
  DOUTR(dout_OTF_MA, "Doing momentumAnalysis::problemSetup");
  //__________________________________
  //  Read in timing information
  m_module_spec->require("samplingFrequency", m_analysisFreq);
  m_module_spec->require("timeStart", d_startTime);
  m_module_spec->require("timeStop", d_stopTime);

  // one patch
  const Patch* p = grid->getPatchByID(0, 0);
  d_zeroPatch    = scinew PatchSet();
  d_zeroPatch->add(p);
  d_zeroPatch->addReference();

  //__________________________________
  // find the material .  Default is matl 0.
  Material* matl{ nullptr };

  // find the material to extract data from.
  if (m_module_spec->findBlock("material")) {
    matl = d_materialManager->parseAndLookupMaterial(m_module_spec, "material");
  } else {
    throw ProblemSetupException(
      "ERROR:AnalysisModule:momentumAnalysis: Missing <material> tag. \n",
      __FILE__,
      __LINE__);
  }

  d_matlIndx = matl->getDWIndex();

  vector<int> m;
  m.push_back(0); // matl index for FileInfo label
  m.push_back(d_matlIndx);

  d_matl_set = scinew MaterialSet();
  d_matl_set->addAll_unique(m);
  d_matl_set->addReference();

  // HARDWIRED FOR ICE/MPMICE
  labels->vel_CC =
    VarLabel::find("vel_CC", "ERROR momentumAnalysis::problemSetup");
  labels->rho_CC =
    VarLabel::find("rho_CC", "ERROR momentumAnalysis::problemSetup");

  labels->uvel_FC =
    VarLabel::find("uvel_FCME", "ERROR momentumAnalysis::problemSetup");
  labels->vvel_FC =
    VarLabel::find("vvel_FCME", "ERROR momentumAnalysis::problemSetup");
  labels->wvel_FC =
    VarLabel::find("wvel_FCME", "ERROR momentumAnalysis::problemSetup");

  labels->pressX_FC =
    VarLabel::find("pressX_FC", "ERROR momentumAnalysis::problemSetup");
  labels->pressY_FC =
    VarLabel::find("pressY_FC", "ERROR momentumAnalysis::problemSetup");
  labels->pressZ_FC =
    VarLabel::find("pressZ_FC", "ERROR momentumAnalysis::problemSetup");

  labels->tau_X_FC =
    VarLabel::find("tau_X_FC", "ERROR momentumAnalysis::problemSetup");
  labels->tau_Y_FC =
    VarLabel::find("tau_Y_FC", "ERROR momentumAnalysis::problemSetup");
  labels->tau_Z_FC =
    VarLabel::find("tau_Z_FC", "ERROR momentumAnalysis::problemSetup");

  //__________________________________
  // Loop over each face and find the extents
  ProblemSpecP ma_ps = m_module_spec->findBlock("controlVolume");
  if (!ma_ps) {
    throw ProblemSetupException(
      "ERROR Radiometer: Couldn't find <controlVolume> xml node",
      __FILE__,
      __LINE__);
  }

  for (ProblemSpecP face_ps = ma_ps->findBlock("Face"); face_ps != nullptr;
       face_ps              = face_ps->findNextBlock("Face")) {

    map<string, string> faceMap;
    face_ps->getAttributes(faceMap);

    string side = faceMap["side"];
    int p_dir;
    int index;
    Vector norm;
    Point start(-9, -9, -9);
    Point end(-9, -9, -9);
    FaceType type = none;

    faceInfo(side, norm, p_dir, index);

    if (faceMap["extents"] == "partialFace") {

      face_ps->get("startPt", start);
      face_ps->get("endPt", end);
      type = partialFace;

      bulletProofing(grid, side, start, end);
    } else {
      type = entireFace;
    }

    // put the input variables into the global struct
    cv_face* cvFace   = scinew cv_face;
    cvFace->p_dir     = p_dir;
    cvFace->normalDir = norm;
    cvFace->face      = type;
    cvFace->startPt   = start;
    cvFace->endPt     = end;
    d_cv_faces[index] = cvFace;
  }
}

//______________________________________________________________________
//
void
momentumAnalysis::scheduleInitialize(SchedulerP& sched, const LevelP& level)
{
  printSchedule(level, dout_OTF_MA, "momentumAnalysis::scheduleInitialize");

  Task* t = scinew Task(
    "momentumAnalysis::initialize", this, &momentumAnalysis::initialize);

  t->computes(labels->lastCompTime);
  t->computes(labels->fileVarsStruct, m_zeroMatl);
  sched->addTask(t, d_zeroPatch, m_zeroMatlSet);
}

//______________________________________________________________________
//
void
momentumAnalysis::initialize(const ProcessorGroup*,
                             const PatchSubset* patches,
                             const MaterialSubset*,
                             DataWarehouse*,
                             DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(
      patches, patch, dout_OTF_MA, "Doing momentumAnalysis::initialize");

    double tminus = d_startTime - 1.0 / m_analysisFreq;
    new_dw->put(max_vartype(tminus), labels->lastCompTime);

    //__________________________________
    //  initialize fileInfo struct
    PerPatch<FileInfoP> fileInfo;
    FileInfo* myFileInfo = scinew FileInfo();
    fileInfo.get()       = myFileInfo;

    new_dw->put(fileInfo, labels->fileVarsStruct, 0, patch);

    if (patch->getGridIndex() == 0) { // only need to do this once
      string udaDir = d_output->getOutputLocation();

      //  Bulletproofing
      DIR* check = opendir(udaDir.c_str());
      if (check == nullptr) {
        ostringstream warn;
        warn
          << "ERROR:momentumAnalysis  The main uda directory does not exist. ";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      closedir(check);
    }
  }
}

//______________________________________________________________________
//
void
momentumAnalysis::scheduleRestartInitialize(SchedulerP& sched,
                                            const LevelP& level)
{
  scheduleInitialize(sched, level);
}

//______________________________________________________________________
//
void
momentumAnalysis::scheduleDoAnalysis(SchedulerP& sched, const LevelP& level)
{

  // Tell the scheduler to not copy this variable to a new AMR grid and
  // do not checkpoint it.
  sched->overrideVariableBehavior(
    "FileInfo_MA", false, false, false, true, true);

  //__________________________________
  //  compute the total momentum and fluxes
  printSchedule(level, dout_OTF_MA, "momentumAnalysis::scheduleDoAnalysis");

  Task* t0 = scinew Task("momentumAnalysis::integrateMomentumField",
                         this,
                         &momentumAnalysis::integrateMomentumField);

  Ghost::GhostType gn = Ghost::None;

  MaterialSubset* matl_SS = scinew MaterialSubset();
  matl_SS->add(d_matlIndx);
  matl_SS->addReference();

  sched_TimeVars(t0, level, labels->lastCompTime, false);

  t0->needs(Task::NewDW, labels->vel_CC, matl_SS, gn);
  t0->needs(Task::NewDW, labels->rho_CC, matl_SS, gn);

  t0->needs(Task::NewDW, labels->uvel_FC, matl_SS, gn);
  t0->needs(Task::NewDW, labels->vvel_FC, matl_SS, gn);
  t0->needs(Task::NewDW, labels->wvel_FC, matl_SS, gn);

  t0->needs(Task::NewDW, labels->pressX_FC, d_pressMatl, gn);
  t0->needs(Task::NewDW, labels->pressY_FC, d_pressMatl, gn);
  t0->needs(Task::NewDW, labels->pressZ_FC, d_pressMatl, gn);

  t0->needs(Task::NewDW, labels->tau_X_FC, matl_SS, gn);
  t0->needs(Task::NewDW, labels->tau_Y_FC, matl_SS, gn);
  t0->needs(Task::NewDW, labels->tau_Z_FC, matl_SS, gn);

  t0->computes(labels->totalCVMomentum);
  t0->computes(labels->convectMom_fluxes);
  t0->computes(labels->viscousMom_fluxes);
  t0->computes(labels->pressForces);

  sched->addTask(t0, level->eachPatch(), d_matl_set);

  //__________________________________
  //  Task that outputs the contributions
  Task* t1 = scinew Task(
    "momentumAnalysis::doAnalysis", this, &momentumAnalysis::doAnalysis);

  sched_TimeVars(t1, level, labels->lastCompTime, true);
  t1->needs(Task::OldDW, labels->fileVarsStruct, m_zeroMatl, gn, 0);

  t1->needs(Task::NewDW, labels->totalCVMomentum);
  t1->needs(Task::NewDW, labels->convectMom_fluxes);
  t1->needs(Task::NewDW, labels->viscousMom_fluxes);
  t1->needs(Task::NewDW, labels->pressForces);

  t1->computes(labels->fileVarsStruct, m_zeroMatl);
  sched->addTask(
    t1, d_zeroPatch, m_zeroMatlSet); // you only need to schedule patch 0 since
                                     // all you're doing is writing out data
}

//______________________________________________________________________
//  Compute the total momentum of the control volume, the fluxes passing
//  through the control surfaces and the pressure forces
//______________________________________________________________________
//
void
momentumAnalysis::integrateMomentumField([[maybe_unused]] const ProcessorGroup* pg,
                                         const PatchSubset* patches,
                                         [[maybe_unused]] const MaterialSubset* matl_sub,
                                         DataWarehouse* old_dw,
                                         DataWarehouse* new_dw)
{

  const Level* level = getLevel(patches);

  // Ignore the task if a recompute time step has been requested upstream
  bool rts      = new_dw->recomputeTimestep();
  bool itItTime = isItTime(old_dw, level, labels->lastCompTime);

  if (itItTime == false || rts) {
    return;
  }

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);

    printTask(patches,
              patch,
              dout_OTF_MA,
              "Doing momentumAnalysis::integrateMomentumField");

    Vector totalCVMomentum = Vector(0., 0., 0.);

    faceQuantities* faceQ = scinew faceQuantities;

    initializeVars(faceQ);

    Vector dx  = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();
    constCCVariable<double> rho_CC;
    constCCVariable<Vector> vel_CC;

    Ghost::GhostType gn = Ghost::None;
    new_dw->get(rho_CC, labels->rho_CC, d_matlIndx, patch, gn, 0);
    new_dw->get(vel_CC, labels->vel_CC, d_matlIndx, patch, gn, 0);

    //__________________________________
    //  Sum the total momentum over the patch
    // This assumes the entire computational domain is being used as the control
    // volume!!!         <<<<<<<<<<<<,
    for (CellIterator iter = patch->getCellIterator(); !iter.done(); iter++) {
      IntVector c = *iter;
      totalCVMomentum += rho_CC[c] * vol * vel_CC[c];
    }

    //__________________________________
    //
    if (patch->hasBoundaryFaces()) {

      constSFCXVariable<double> uvel_FC, pressX_FC;
      constSFCYVariable<double> vvel_FC, pressY_FC;
      constSFCZVariable<double> wvel_FC, pressZ_FC;

      constSFCXVariable<Vector> tau_X_FC;
      constSFCYVariable<Vector> tau_Y_FC;
      constSFCZVariable<Vector> tau_Z_FC;

      new_dw->get(uvel_FC, labels->uvel_FC, d_matlIndx, patch, gn, 0);
      new_dw->get(vvel_FC, labels->vvel_FC, d_matlIndx, patch, gn, 0);
      new_dw->get(wvel_FC, labels->wvel_FC, d_matlIndx, patch, gn, 0);

      new_dw->get(pressX_FC, labels->pressX_FC, d_pressIndx, patch, gn, 0);
      new_dw->get(pressY_FC, labels->pressY_FC, d_pressIndx, patch, gn, 0);
      new_dw->get(pressZ_FC, labels->pressZ_FC, d_pressIndx, patch, gn, 0);

      new_dw->get(tau_X_FC, labels->tau_X_FC, d_matlIndx, patch, gn, 0);
      new_dw->get(tau_Y_FC, labels->tau_Y_FC, d_matlIndx, patch, gn, 0);
      new_dw->get(tau_Z_FC, labels->tau_Z_FC, d_matlIndx, patch, gn, 0);

      //__________________________________
      // Sum the fluxes passing through control volume surface
      // and the pressure forces
      vector<Patch::FaceType> bf;
      patch->getBoundaryFaces(bf);

      for (vector<Patch::FaceType>::const_iterator itr = bf.begin();
           itr != bf.end();
           ++itr) {

        Patch::FaceType face = *itr;
        string faceName      = patch->getFaceName(face);
        cv_face* cvFace      = d_cv_faces[face];

        DOUTR(dbg_OTF_MA,
              "cvFace: " << faceName << " faceType " << cvFace->face
                         << " startPt: " << cvFace->startPt
                         << " endPt: " << cvFace->endPt);
        DOUTR(dbg_OTF_MA,
              "          norm: " << cvFace->normalDir
                                 << " p_dir: " << cvFace->p_dir);

        // define the iterator on this face  The default is the entire face
        Patch::FaceIteratorType SFC = Patch::SFCVars;
        CellIterator iterLimits     = patch->getFaceIterator(face, SFC);

        //__________________________________
        //  partial face iterator
        if (cvFace->face == partialFace) {

          IntVector lo  = level->getCellIndex(cvFace->startPt);
          IntVector hi  = level->getCellIndex(cvFace->endPt);
          IntVector pLo = patch->getCellLowIndex();
          IntVector pHi = patch->getCellHighIndex();

          IntVector low  = Max(lo, pLo); // find the intersection
          IntVector high = Min(hi, pHi);

          //__________________________________
          // enlarge the iterator by oneCell
          // x-           x+        y-       y+       z-        z+
          // (-1,0,0)  (1,0,0)  (0,-1,0)  (0,1,0)  (0,0,-1)  (0,0,1)
          IntVector oneCell = patch->faceDirection(face);
          if (face == Patch::xminus || face == Patch::yminus ||
              face == Patch::zminus) {
            low += oneCell;
          }
          if (face == Patch::xplus || face == Patch::yplus ||
              face == Patch::zplus) {
            high += oneCell;
          }

          iterLimits = CellIterator(low, high);
        }

        //__________________________________
        //           X faces
        if (face == Patch::xminus || face == Patch::xplus) {
          double area = dx.y() * dx.z();
          integrateOverFace(faceName,
                            area,
                            iterLimits,
                            faceQ,
                            uvel_FC,
                            pressX_FC,
                            tau_X_FC,
                            rho_CC,
                            vel_CC);
        }

        //__________________________________
        //        Y faces
        if (face == Patch::yminus || face == Patch::yplus) {
          double area = dx.x() * dx.z();
          integrateOverFace(faceName,
                            area,
                            iterLimits,
                            faceQ,
                            vvel_FC,
                            pressY_FC,
                            tau_Y_FC,
                            rho_CC,
                            vel_CC);
        }
        //__________________________________
        //        Z faces
        if (face == Patch::zminus || face == Patch::zplus) {
          double area = dx.x() * dx.y();
          integrateOverFace(faceName,
                            area,
                            iterLimits,
                            faceQ,
                            wvel_FC,
                            pressZ_FC,
                            tau_Z_FC,
                            rho_CC,
                            vel_CC);
        }
      } // boundary faces
    }   // patch has faces

    //__________________________________
    //  Now compute the net fluxes from the face quantites
    Vector net_convect_flux = L_minus_R(faceQ->convect_faceFlux);

    Vector net_viscous_flux = L_minus_R(faceQ->viscous_faceFlux);

    // net force on control volume due to pressure forces
    map<int, double> pressForce = faceQ->pressForce_face; // for readability
    double pressForceX          = pressForce[0] - pressForce[1];
    double pressForceY          = pressForce[2] - pressForce[3];
    double pressForceZ          = pressForce[4] - pressForce[5];

    Vector net_press_forces(pressForceX, pressForceY, pressForceZ);

    //__________________________________
    // put in the dw
    // cout.precision(15);
    // cout <<  " Total CV momentum: " << totalCVMomentum << " Net
    // convectiveFlux: " << net_convect_flux << " viscousFlux: " <<
    // net_viscous_flux << " pressForce " << net_press_forces << std::endl;

    new_dw->put(sumvec_vartype(totalCVMomentum), labels->totalCVMomentum);
    new_dw->put(sumvec_vartype(net_convect_flux), labels->convectMom_fluxes);
    new_dw->put(sumvec_vartype(net_viscous_flux), labels->viscousMom_fluxes);
    new_dw->put(sumvec_vartype(net_press_forces), labels->pressForces);
  } // patch loop
}

//______________________________________________________________________
//
void
momentumAnalysis::doAnalysis([[maybe_unused]] const ProcessorGroup* pg,
                             const PatchSubset* patches,
                             [[maybe_unused]] const MaterialSubset* matls,
                             DataWarehouse* old_dw,
                             DataWarehouse* new_dw)
{

  // Ignore the task if a recompute time step has been requested upstream
  bool rts = new_dw->recomputeTimestep();

  const Level* level = getLevel(patches);
  timeVars tv;

  getTimeVars(old_dw, level, labels->lastCompTime, tv);
  putTimeVars(new_dw, labels->lastCompTime, tv);

  if (rts || tv.isItTime == false) {
    return;
  }

  //__________________________________
  //
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    printTask(
      patches, patch, dout_OTF_MA, "Doing momentumAnalysis::doAnalysis");

    //__________________________________
    // open the struct that contains the file pointer map.  We use FileInfoP
    // types and store them in the DW to avoid doing system calls (SLOW). Note:
    // after regridding this may not exist for this patch in the old_dw
    PerPatch<FileInfoP> fileInfo;

    if (old_dw->exists(labels->fileVarsStruct, 0, patch)) {
      old_dw->get(fileInfo, labels->fileVarsStruct, 0, patch);
    } else {
      FileInfo* myFileInfo = scinew FileInfo();
      fileInfo.get()       = myFileInfo;
    }

    std::map<string, FILE*> myFiles;

    if (fileInfo.get().get_rep()) {
      myFiles = fileInfo.get().get_rep()->files;
    }

    string udaDir   = d_output->getOutputLocation();
    string filename = udaDir + "/" + "momentumAnalysis.dat";
    FILE* fp        = nullptr;

    if (myFiles.count(filename) == 0) {
      createFile(filename, fp);
      myFiles[filename] = fp;
    } else {
      fp = myFiles[filename];
    }

    if (!fp) {
      throw InternalError(
        "\nERROR:dataAnalysisModule:momentumAnalysis:  failed opening file" +
          filename,
        __FILE__,
        __LINE__);
    }
    //__________________________________
    //
    sumvec_vartype totalCVMomentum, convectFlux, viscousFlux, pressForce;
    new_dw->get(totalCVMomentum, labels->totalCVMomentum);
    new_dw->get(convectFlux, labels->convectMom_fluxes);
    new_dw->get(viscousFlux, labels->viscousMom_fluxes);
    new_dw->get(pressForce, labels->pressForces);

    // so fprintf can deal with it
    Vector momentum = totalCVMomentum;
    Vector conFlux  = convectFlux;
    Vector visFlux  = viscousFlux;
    Vector pForce   = pressForce;

    fprintf(fp,
            "%16.15E,   %16.15E,   %16.15E,   %16.15E,   %16.15E,   %16.15E,   "
            "%16.15E,   %16.15E,   %16.15E,   %16.15E,   %16.15E,   %16.15E,   "
            "%16.15E\n",
            tv.now,
            (double)momentum.x(),
            (double)momentum.y(),
            (double)momentum.z(),
            (double)conFlux.x(),
            (double)conFlux.y(),
            (double)conFlux.z(),
            (double)visFlux.x(),
            (double)visFlux.y(),
            (double)visFlux.z(),
            (double)pForce.x(),
            (double)pForce.y(),
            (double)pForce.z());

    //      fflush(fp);   If you want to write the data right now, no buffering.

    //__________________________________
    // put the file pointers into the DataWarehouse
    // these could have been altered. You must
    // reuse the Handle fileInfo and just replace the contents
    fileInfo.get().get_rep()->files = myFiles;
    new_dw->put(fileInfo, labels->fileVarsStruct, 0, patch);
  }
}

//______________________________________________________________________
//  Integrate fluxes/forces over the control volume face
//______________________________________________________________________
//
template<class SFC_D, class SFC_V>
void
momentumAnalysis::integrateOverFace(const std::string faceName,
                                    const double faceArea,
                                    CellIterator iter,
                                    faceQuantities* faceQ,
                                    SFC_D& vel_FC,
                                    SFC_D& press_FC,
                                    SFC_V& tau_FC,
                                    constCCVariable<double>& rho_CC,
                                    constCCVariable<Vector>& vel_CC)
{
  int dir; // x, y or z
  int f;   // face index
  Vector norm;
  faceInfo(faceName, norm, dir, f);

  Vector convect_flux(0.);
  Vector viscous_flux(0.);
  double pressForce(0.);

  //__________________________________
  //  Loop over a face
  for (; !iter.done(); iter++) {
    IntVector c = *iter;
    double vel  = vel_FC[c];

    // find upwind cell
    IntVector uw = c;
    if (vel > 0) {
      uw[dir] = c[dir] - 1;
    }

    // One way to define m dot through face
    double mdot = vel * faceArea * rho_CC[uw];

    // Another way
    // Vector mdotV  = faceArea * vel_CC[uw] * rho_CC[uw];
    // double mdot = mdotV[dir];

    convect_flux += mdot * vel_CC[uw];
    viscous_flux += (faceArea * tau_FC[c]);
    pressForce += (faceArea * press_FC[c]);
  }
  faceQ->convect_faceFlux[f] = convect_flux;
  faceQ->viscous_faceFlux[f] = viscous_flux;
  faceQ->pressForce_face[f]  = pressForce;

  //  cout << "face: " << faceName << "\t dir: " << dir << " convect_Flux = " <<
  //  convect_flux << " ViscousFlux " << viscous_flux << " pressForce " <<
  //  pressForce << std::endl;
}

//______________________________________________________________________
//  Open the file if it doesn't exist and write the file header
//______________________________________________________________________
//
void
momentumAnalysis::createFile(string& filename, FILE*& fp)
{
  // if the file already exists then exit.  The file could exist but not be
  // owned by this processor
  ifstream doExists(filename.c_str());
  if (doExists) {
    fp = fopen(filename.c_str(), "a");
    return;
  }

  fp = fopen(filename.c_str(), "w");
  fprintf(fp,
          "#                                                 total momentum in "
          "the control volume                                          Net "
          "convective momentum flux                                            "
          "   net viscous flux                                                 "
          "            pressure force on control vol.\n");
  fprintf(fp,
          "#Time                    CV_mom.x                 CV_mom.y          "
          "        CV_mom.z                  momFlux.x               momFlux.y "
          "               momFlux.z                 visFlux.x                 "
          "visFlux.y                visFlux.z                 pressForce.x     "
          "         pressForce.y             pressForce.z\n");

  proc0cout << Parallel::getMPIRank() << " momentumAnalysis:Created file "
            << filename << std::endl;
}

//______________________________________________________________________
//   This is a rip off of what's done in the boundary condition code
//______________________________________________________________________
//
void
momentumAnalysis::faceInfo(const std::string fc,
                           Vector& norm,
                           int& p_dir,
                           int& index)
{
  if (fc == "x-" || fc == "xminus") {
    norm  = Vector(-1, 0, 0);
    p_dir = 0;
    index = 0;
    return;
  }
  if (fc == "x+" || fc == "xplus") {
    norm  = Vector(1, 0, 0);
    p_dir = 0;
    index = 1;
    return;
  }
  if (fc == "y-" || fc == "yminus") {
    norm  = Vector(0, -1, 0);
    p_dir = 1;
    index = 2;
    return;
  }
  if (fc == "y+" || fc == "yplus") {
    norm  = Vector(0, 1, 0);
    p_dir = 1;
    index = 3;
    return;
  }
  if (fc == "z-" || fc == "zminus") {
    norm  = Vector(0, 0, -1);
    p_dir = 2;
    index = 4;
    return;
  }
  if (fc == "z+" || fc == "zplus") {
    norm  = Vector(0, 0, 1);
    p_dir = 2;
    index = 5;
    return;
  }

  ostringstream warn;
  warn << " ERROR:MomentumAnalysis face name (" << fc << ") unknown. ";

  throw InternalError(warn.str(), __FILE__, __LINE__);
}
//______________________________________________________________________
//  bulletProofing on the user inputs
//______________________________________________________________________
//
void
momentumAnalysis::bulletProofing(GridP& grid,
                                 const string& side,
                                 const Point& start,
                                 const Point& end)
{
  //__________________________________
  // plane must be parallel to the coordinate system
  // plane must be contained in the domain

  bulletProofing_LinesPlanes(
    objectType::plane, grid, "1stLawThermo", start, end);

  //__________________________________
  //  plane must be on the edge of the domain
  bool validPlane = true;
  BBox compDomain;
  grid->getInteriorSpatialRange(compDomain);
  Point min = compDomain.min();
  Point max = compDomain.max();

  Point me = min;
  if (side == "x+" || side == "y+" || side == "z+") {
    me = max;
  }

  if (side == "x+" || side == "x-") {
    if (start.x() != me.x()) {
      validPlane = false;
    }
  }
  if (side == "y+" || side == "y-") {
    if (start.y() != me.y()) {
      validPlane = false;
    }
  }
  if (side == "z+" || side == "z-") {
    if (start.z() != me.z()) {
      validPlane = false;
    }
  }
  if (validPlane == false) {
    ostringstream warn;
    warn << "\n ERROR:1stLawThermo: the plane on face (" << side
         << ") that you've specified " << start << " to " << end
         << " is not at the edge of the computational domain. \n"
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

//______________________________________________________________________
//    Initialize the face quantities
//______________________________________________________________________
//
void
momentumAnalysis::initializeVars(faceQuantities* faceQ)
{
  for (int f = 0; f < 6; f++) {
    faceQ->convect_faceFlux[f] = Vector(0.);
    faceQ->viscous_faceFlux[f] = Vector(0.);
    faceQ->pressForce_face[f]  = 0.;
  }
}

//______________________________________________________________________
//
Vector
momentumAnalysis::L_minus_R(std::map<int, Vector>& faceFlux)
{
  double X_flux = (faceFlux[0].x() - faceFlux[1].x()) +
                  (faceFlux[2].x() - faceFlux[3].x()) +
                  (faceFlux[4].x() - faceFlux[5].x());

  double Y_flux = (faceFlux[0].y() - faceFlux[1].y()) +
                  (faceFlux[2].y() - faceFlux[3].y()) +
                  (faceFlux[4].y() - faceFlux[5].y());

  double Z_flux = (faceFlux[0].z() - faceFlux[1].z()) +
                  (faceFlux[2].z() - faceFlux[3].z()) +
                  (faceFlux[4].z() - faceFlux[5].z());

  Vector net_flux(X_flux, Y_flux, Z_flux);
  return net_flux;
}
