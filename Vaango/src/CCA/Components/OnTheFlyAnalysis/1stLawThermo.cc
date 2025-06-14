/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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

#include <CCA/Components/OnTheFlyAnalysis/1stLawThermo.h>
#include <CCA/Components/OnTheFlyAnalysis/FileInfoVar.h>
#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/DbgOutput.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Variables/SFCXVariable.h>
#include <Core/Grid/Variables/SFCYVariable.h>
#include <Core/Grid/Variables/SFCZVariable.h>

#include <Core/OS/Dir.h> // for MKDIR
#include <Core/Util/FileUtils.h>
#include <Core/Util/DebugStream.h>

#include <dirent.h>
#include <iostream>
#include <cstdio>

using namespace Uintah;
using namespace std;

Dout dout_OTF_FLT("FirstLawThermo",     "OnTheFlyAnalysis", "Task scheduling and execution.", false);
Dout dbg_OTF_FLT("FirstLawThermo_dbg",  "OnTheFlyAnalysis", "Displays detailed debugging info.", false);

//______________________________________________________________________
FirstLawThermo::FirstLawThermo(const ProcessorGroup* myworld,
                               const MaterialManagerP& materialManager,
                               const ProblemSpecP& module_spec)

  : AnalysisModule(myworld, materialManager, module_spec)
{

  d_zeroPatch    = 0;
  d_conversion   = 1.0/1000.;     // in SI units this is J/KJ

  FL_lb = scinew FL_Labels();
  I_lb  = scinew ICELabel();
  M_lb  = scinew MPMLabel();

  FL_lb->lastCompTimeLabel    = VarLabel::create( "lastCompTime",     max_vartype::getTypeDescription() );
  FL_lb->fileVarsStructLabel  = VarLabel::create( "FileInfo_1stLaw",  PerPatch<FileInfoP>::getTypeDescription() );
  FL_lb->ICE_totalIntEngLabel = VarLabel::create( "ICE_totalIntEng",  sum_vartype::getTypeDescription() );
  FL_lb->MPM_totalIntEngLabel = VarLabel::create( "MPM_totalIntEng",  sum_vartype::getTypeDescription() );
  FL_lb->totalFluxesLabel     = VarLabel::create( "totalFluxes",      sum_vartype::getTypeDescription() );
}

//______________________________________________________________________
//
FirstLawThermo::~FirstLawThermo()
{
  DOUTR(dout_OTF_FLT, " Doing: destructor FirstLawThermo " );

  if( d_zeroPatch && d_zeroPatch->removeReference() ) {
    delete d_zeroPatch;
  }

  VarLabel::destroy( FL_lb->lastCompTimeLabel );
  VarLabel::destroy( FL_lb->fileVarsStructLabel );
  VarLabel::destroy( FL_lb->ICE_totalIntEngLabel );
  VarLabel::destroy( FL_lb->MPM_totalIntEngLabel );
  VarLabel::destroy( FL_lb->totalFluxesLabel );

  // delete each plane
  for( auto iter  = d_cv_faces.begin();iter != d_cv_faces.end(); iter++){
    delete iter->second;
  }

  delete FL_lb;
  delete M_lb;
  delete I_lb;
}

//______________________________________________________________________
//     P R O B L E M   S E T U P
void FirstLawThermo::problemSetup(const ProblemSpecP &,
                                  const ProblemSpecP & restart_prob_spec,
                                  GridP              & grid,
                                  [[maybe_unused]] std::vector<std::vector<const VarLabel* > > &PState,
                                  [[maybe_unused]] std::vector<std::vector<const VarLabel* > > &PState_preReloc)
{
  DOUTR(dout_OTF_FLT, "Doing FirstLawThermo::problemSetup");

  //__________________________________
  //  Read in timing information
  m_module_spec->require( "samplingFrequency", m_analysisFreq );
  m_module_spec->require( "timeStart",         d_StartTime );
  m_module_spec->require( "timeStop",          d_StopTime );
  m_module_spec->get(     "engy_convt_factor", d_conversion );   // energy conversion factor in SI it KJ->J

  // one patch
  const Patch* p = grid->getPatchByID(0,0);
  d_zeroPatch = scinew PatchSet();
  d_zeroPatch->add(p);
  d_zeroPatch->addReference();

  //__________________________________
  // Loop over each face and find the extents
  ProblemSpecP cv_ps = m_module_spec->findBlock("controlVolume");

  for( ProblemSpecP face_ps = cv_ps->findBlock( "Face" ); face_ps != nullptr; face_ps=face_ps->findNextBlock( "Face" ) ) {

    map<string,string> faceMap;
    face_ps->getAttributes(faceMap);

    string side = faceMap["side"];
    int p_dir;
    Vector norm;
    Patch::FaceType f;
    Point start(-9,-9,-9);
    Point end(-9,-9,-9);
    FaceType type=none;

    faceInfo(side, f, norm, p_dir);

    if (faceMap["extents"] == "partialFace"){

      face_ps->get( "startPt", start );
      face_ps->get( "endPt",   end );
      type = partialFace;

      bulletProofing(grid, side, start, end);
    }else{
      type = entireFace;
    }
    // put the input variables into the global struct
    cv_face* cvFace   = scinew cv_face;
    cvFace->p_dir     = p_dir;
    cvFace->normalDir = norm;
    cvFace->face      = type;
    cvFace->startPt   = start;
    cvFace->endPt     = end;

    d_cv_faces[f] = cvFace;
  }

  //__________________________________
  //  Loop over all the MPM Matls and pull out the specific heat
  //  This assumes that MPM matls are listed from 0 to N
  //  It also assumes that cp is constant
  // If we are doing a restart, then use the "timestep.xml"

  // first try prob_spec
  ProblemSpecP root_ps = m_module_spec->getRootNode();
  ProblemSpecP mat_ps = root_ps->findBlockWithOutAttribute( "MaterialProperties" );

  if ( !mat_ps && restart_prob_spec ){   // read in from checkpoint/timestep.xml on a restart
    ProblemSpecP silly = restart_prob_spec;
    root_ps = silly->getRootNode();
    mat_ps  = root_ps->findBlockWithOutAttribute( "MaterialProperties" );
  }

  int matl = 0;

  ProblemSpecP mpm_mat_ps = mat_ps->findBlock("MPM");
  if( mpm_mat_ps ){
    for( ProblemSpecP ps = mpm_mat_ps->findBlock( "material" ); ps != nullptr; ps = ps->findNextBlock( "material" ) ) {
      double cp;
      ps->require("specific_heat",cp);
      d_mpm_specificHeat[matl] = cp;
      matl +=1;
    }
  }
}

//______________________________________________________________________
//
void
FirstLawThermo::scheduleInitialize( SchedulerP   & sched,
                                    const LevelP & level )
{
  printSchedule(level,dout_OTF_FLT,"FirstLawThermo::scheduleInitialize");

  Task* t = scinew Task( "FirstLawThermo::initialize", this, &FirstLawThermo::initialize );

  t->computes(FL_lb->lastCompTimeLabel);
  t->computes(FL_lb->fileVarsStructLabel, m_zeroMatl);

  sched->addTask(t, level->eachPatch(), m_zeroMatlSet);
}

//______________________________________________________________________
//
void
FirstLawThermo::initialize( const ProcessorGroup *,
                            const PatchSubset    * patches,
                            const MaterialSubset *,
                                  DataWarehouse  *,
                                  DataWarehouse  * new_dw)
{
  double tminus = d_startTime - 1.0/m_analysisFreq;
  new_dw->put(max_vartype(tminus), FL_lb->lastCompTimeLabel);

  for( int p = 0; p < patches->size(); p++ ) {
    const Patch* patch = patches->get(p);
    printTask(patches, patch,dout_OTF_FLT,"Doing FirstLawThermo::initialize");

    //__________________________________
    //  initialize fileInfo struct
    PerPatch<FileInfoP> fileInfo;
    FileInfo* myFileInfo = scinew FileInfo();
    fileInfo.get() = myFileInfo;

    new_dw->put(fileInfo,    FL_lb->fileVarsStructLabel, 0, patch);

    if(patch->getGridIndex() == 0){   // only need to do this once
      string udaDir = d_output->getOutputLocation();

      //  Bulletproofing
      DIR *check = opendir(udaDir.c_str());
      if ( check == nullptr){
        ostringstream warn;
        warn << "ERROR:FirstLawThermo  The main uda directory does not exist. ";
        throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
      }
      closedir(check);
    }
  }
}


//______________________________________________________________________
//
void
FirstLawThermo::scheduleRestartInitialize(SchedulerP   & sched,
                                          const LevelP & level)
{
  scheduleInitialize( sched, level);
}

//______________________________________________________________________
//
void FirstLawThermo::scheduleDoAnalysis(SchedulerP   & sched,
                                        const LevelP & level)
{

  // Tell the scheduler to not copy this variable to a new AMR grid and
  // do not checkpoint it.
  sched->overrideVariableBehavior("FileInfo_1stLaw", false, false, false, true, true);

  //__________________________________
  //  compute the ICE contributions
  printSchedule( level, dout_OTF_FLT,"FirstLawThermo::scheduleDoAnalysis" );

  Task* t0 = scinew Task( "FirstLawThermo::compute_ICE_Contributions",
                     this,&FirstLawThermo::compute_ICE_Contributions );

  Ghost::GhostType  gn  = Ghost::None;
  const MaterialSet* all_matls = d_materialManager->allMaterials();
  const MaterialSet* ice_matls = d_materialManager->allMaterials( "ICE" );
  const MaterialSubset* ice_ss = ice_matls->getUnion();

  sched_TimeVars( t0, level, FL_lb->lastCompTimeLabel, false );

  t0->needs( Task::NewDW, I_lb->rho_CCLabel,        ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->temperature_CCLabel,       ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->specific_heatLabel, ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->gammaLabel,         ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->uvel_FCMELabel,     ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->vvel_FCMELabel,     ice_ss, gn );
  t0->needs( Task::NewDW, I_lb->wvel_FCMELabel,     ice_ss, gn );

  t0->computes( FL_lb->ICE_totalIntEngLabel );
  t0->computes( FL_lb->totalFluxesLabel );

  sched->addTask( t0, level->eachPatch(), all_matls );

  //__________________________________
  //  compute the MPM contributions
  Task* t1 = scinew Task( "FirstLawThermo::compute_MPM_Contributions",
                     this,&FirstLawThermo::compute_MPM_Contributions );

  const MaterialSet* mpm_matls = d_materialManager->allMaterials( "MPM" );
  const MaterialSubset* mpm_ss = mpm_matls->getUnion();

  sched_TimeVars( t1, level, FL_lb->lastCompTimeLabel, false );

  t1->needs( Task::NewDW, M_lb->pMassLabel_preReloc,        mpm_ss, gn );
  t1->needs( Task::NewDW, M_lb->pTemperatureLabel_preReloc, mpm_ss, gn );
  t1->computes( FL_lb->MPM_totalIntEngLabel );

  sched->addTask( t1, level->eachPatch(), all_matls );


  //__________________________________
  //  output the contributions
  Task* t2 = scinew Task("FirstLawThermo::doAnalysis",
                    this,&FirstLawThermo::doAnalysis );

  sched_TimeVars( t2, level, FL_lb->lastCompTimeLabel, true );

  t2->needs( Task::OldDW, FL_lb->fileVarsStructLabel, m_zeroMatl, gn, 0 );
  t2->needs( Task::NewDW, FL_lb->ICE_totalIntEngLabel );
  t2->needs( Task::NewDW, FL_lb->MPM_totalIntEngLabel );
  t2->needs( Task::NewDW, FL_lb->totalFluxesLabel );

  t2->computes( FL_lb->fileVarsStructLabel, m_zeroMatl );

  // you only need to schedule this  patch 0 since you're writing out data
  sched->addTask( t2, d_zeroPatch, m_zeroMatlSet);

}

//______________________________________________________________________
//        ICE Contributions to the energy
void FirstLawThermo::compute_ICE_Contributions([[maybe_unused]] const ProcessorGroup * pg,
                                               const PatchSubset    * patches,
                                               [[maybe_unused]] const MaterialSubset * matl_sub ,
                                               DataWarehouse        * old_dw,
                                               DataWarehouse        * new_dw)
{
  const Level* level = getLevel(patches);

  if( isItTime(old_dw, level, FL_lb->lastCompTimeLabel) == false){
    return;
  }

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch, dout_OTF_FLT,"Doing FirstLawThermo::compute_ICE_Contributions");

    CCVariable<double> int_eng;
    constCCVariable<double> temp_CC;
    constCCVariable<double> rho_CC;
    constCCVariable<double> cv;
    constCCVariable<double> gamma;

    constSFCXVariable<double> uvel_FC;
    constSFCYVariable<double> vvel_FC;
    constSFCZVariable<double> wvel_FC;

    new_dw->allocateTemporary(int_eng,patch);

    Vector dx = patch->dCell();
    double vol = dx.x() * dx.y() * dx.z();

    int numICEmatls = d_materialManager->getNumMaterials( "ICE" );

    double ICE_totalIntEng = 0.0;
    double total_flux      = 0.0;

    for (int m = 0; m < numICEmatls; m++ ) {

      ICEMaterial* ice_matl = (ICEMaterial*) d_materialManager->getMaterial( "ICE", m);
      int indx = ice_matl->getDWIndex();
      new_dw->get(rho_CC,  I_lb->rho_CCLabel,       indx, patch, m_gn,0);
      new_dw->get(temp_CC, I_lb->temperature_CCLabel,      indx, patch, m_gn,0);
      new_dw->get(cv,      I_lb->specific_heatLabel,indx, patch, m_gn,0);
      new_dw->get(uvel_FC, I_lb->uvel_FCMELabel,    indx, patch, m_gn,0);
      new_dw->get(vvel_FC, I_lb->vvel_FCMELabel,    indx, patch, m_gn,0);
      new_dw->get(wvel_FC, I_lb->wvel_FCMELabel,    indx, patch, m_gn,0);
      new_dw->get(gamma,   I_lb->gammaLabel,        indx, patch, m_gn,0);

      double mat_int_eng = 0.0;

      //__________________________________
      //  Sum contributions over patch
      for (CellIterator iter=patch->getCellIterator();!iter.done();iter++){
        IntVector c = *iter;
        int_eng[c] = rho_CC[c] * vol * cv[c] * temp_CC[c];
        mat_int_eng += int_eng[c];
      }

      ICE_totalIntEng += mat_int_eng;
      //__________________________________
      // Sum the fluxes passing through the boundaries
      vector<Patch::FaceType> bf;
      patch->getBoundaryFaces(bf);
      double mat_fluxes = 0.0;

      for( vector<Patch::FaceType>::const_iterator itr = bf.begin(); itr != bf.end(); ++itr ){
        Patch::FaceType face = *itr;
        string faceName = patch->getFaceName(face );
        cv_face* cvFace = d_cv_faces[face];

        DOUTR(dbg_OTF_FLT, "cvFace: " <<  faceName << " faceType " << cvFace->face
                        << " startPt: " << cvFace->startPt << " endPt: " << cvFace->endPt);
        DOUTR(dbg_OTF_FLT, "          norm: " << cvFace->normalDir << " p_dir: " << cvFace->p_dir );


        // define the iterator on this face  The defauls is the entire face
        Patch::FaceIteratorType SFC = Patch::SFCVars;
        CellIterator iterLimits=patch->getFaceIterator(face, SFC);

        //__________________________________
        //
        if( cvFace->face == partialFace ){

          IntVector lo  = level->getCellIndex( cvFace->startPt );
          IntVector hi  = level->getCellIndex( cvFace->endPt );
          IntVector pLo = patch->getExtraCellLowIndex();
          IntVector pHi = patch->getExtraCellHighIndex();

          IntVector low  = Max(lo, pLo);    // find the intersection
          IntVector high = Min(hi, pHi);

          //__________________________________
          // enlarge the iterator by oneCell
          // x-           x+        y-       y+       z-        z+
          // (-1,0,0)  (1,0,0)  (0,-1,0)  (0,1,0)  (0,0,-1)  (0,0,1)
          IntVector oneCell = patch->faceDirection( face );
          if( face == Patch::xminus || face == Patch::yminus || face == Patch::zminus) {
            low += oneCell;
          }
          if( face == Patch::xplus || face == Patch::yplus || face == Patch::zplus) {
            high += oneCell;
          }

          iterLimits = CellIterator(low,high);
        }

        IntVector axes = patch->getFaceAxes(face);
        int P_dir = axes[0];  // principal direction
        double plus_minus_one = (double) patch->faceDirection(face)[P_dir];

        DOUTR(dbg_OTF_FLT, "    face Direction " << patch->faceDirection(face) );

        //__________________________________
        //           X faces
        if (face == Patch::xminus || face == Patch::xplus) {
          double area = dx.y() * dx.z();
          double sumKE   = 0;
          double sumH    = 0;
          double sumMdot = 0;
          DOUTR(dbg_OTF_FLT, "    iterLimits: " << iterLimits );

          for(CellIterator iter = iterLimits; !iter.done();iter++) {
            IntVector c = *iter;
            double vel = uvel_FC[c];

 #if 0
            // upwinding
            IntVector offset(0,0,0);
            if (vel > 0 ){
              offset = IntVector (-1,0,0);
            }
            c + offset;
 #endif
            // compute the average values
            IntVector cc = c - IntVector(1,0,0);

            double mdot   = plus_minus_one * vel * area * (rho_CC[c] + rho_CC[cc])/2.0;
            double KE     = 0.5 * vel * vel;
            double enthpy = ( temp_CC[c]  * gamma[c]  * cv[c] +
                              temp_CC[cc] * gamma[cc] * cv[cc] )/2.0;

            sumKE   += mdot * KE;
            sumH    += mdot * enthpy;
            sumMdot += mdot;

            mat_fluxes +=  mdot * (enthpy + KE * d_conversion);
            //cout << "face: " << faceName << " c: " << c << " offset: " << offset << " vel = " << vel << " mdot = " << mdot << std::endl;
          }
          DOUTR(dbg_OTF_FLT, "    face: " << faceName << " mdot = " << sumMdot << "      sum of KE = " << sumKE << "     sum H = " << sumH <<  "      sum mat_fluxes = " << mat_fluxes );
        }

        //__________________________________
        //        Y faces
        if (face == Patch::yminus || face == Patch::yplus) {
          double area = dx.x() * dx.z();
          double sumKE   = 0;
          double sumH    = 0;
          double sumMdot = 0;
          DOUTR(dbg_OTF_FLT, "    iterLimits: " << iterLimits );


          for(CellIterator iter = iterLimits; !iter.done();iter++) {
            IntVector c = *iter;
            double vel = vvel_FC[c];

 #if 0
            // upwinding
            IntVector offset(0,0,0);

            if (vel > 0 ){
              offset = IntVector (0,-1,0);
            }
            c + offset;
#endif

            // compute the average values
            IntVector cc = c - IntVector(0,1,0);

            double mdot   = plus_minus_one * vel * area * (rho_CC[c] + rho_CC[cc])/2.0;
            double KE     = 0.5 * vel * vel;
            double enthpy = ( temp_CC[c]  * gamma[c]  * cv[c] +
                              temp_CC[cc] * gamma[cc] * cv[cc] )/2.0;

            sumH    += mdot * enthpy;
            sumKE   += mdot * KE;
            sumMdot += mdot;

            mat_fluxes +=  mdot * (enthpy + KE * d_conversion);
            //cout << "face: " << faceName << " c: " << c << " offset: " << offset << " vel = " << vel << " mdot = " << mdot << std::endl;
          }
           DOUTR(dbg_OTF_FLT, "    face: " << faceName << " mdot = "<< sumMdot << "     sum of KE = " << sumKE << "     sum H = " << sumH << "      sum mat_fluxes = " << mat_fluxes );
        }

        //__________________________________
        //        Z faces
        if (face == Patch::zminus || face == Patch::zplus) {
          double area = dx.x() * dx.y();
          double sumKE   = 0;
          double sumH    = 0;
          double sumMdot = 0;
          DOUTR(dbg_OTF_FLT, "    iterLimits: " << iterLimits );

          for(CellIterator iter = iterLimits; !iter.done();iter++) {
            IntVector c = *iter;
            double vel = wvel_FC[c];
#if 0
            // upwinding
            IntVector offset(0,0,0);
            if (vel > 0 ){
              offset = IntVector (0,0,-1);
            }
            c + offset;
#endif

            // compute the average values
            IntVector cc = c - IntVector(0,0,1);

            double mdot   = plus_minus_one * vel * area * (rho_CC[c] + rho_CC[cc])/2.0;
            double KE     = 0.5 * vel * vel;
            double enthpy = ( temp_CC[c]  * gamma[c]  * cv[c] +
                              temp_CC[cc] * gamma[cc] * cv[cc] )/2.0;

            sumH  += mdot * enthpy;
            sumKE += mdot * KE;
            sumMdot += mdot;
            mat_fluxes +=  mdot * (enthpy + KE * d_conversion);
          }
           DOUTR(dbg_OTF_FLT, "    face: " << faceName << " mdot = "<< sumMdot << "     sum of KE = " << sumKE << "     sum H = " << sumH << "      sum mat_fluxes = " << mat_fluxes );
        }
      }  // boundary faces

      total_flux += mat_fluxes;
    }  // ICE Matls loop

    DOUTR(dbg_OTF_FLT, "Patch: " << patch->getID() << " totalFlux: " << total_flux );

    new_dw->put( sum_vartype(ICE_totalIntEng), FL_lb->ICE_totalIntEngLabel );
    new_dw->put( sum_vartype(total_flux),      FL_lb->totalFluxesLabel );
  }  // patch loop
}


//______________________________________________________________________
//        MPM Contributions to the energy
void FirstLawThermo::compute_MPM_Contributions([[maybe_unused]] const ProcessorGroup * pg,
                                               const PatchSubset    * patches,
                                               [[maybe_unused]] const MaterialSubset * matl_sub ,
                                               DataWarehouse        * old_dw,
                                               DataWarehouse        * new_dw)
{
  const Level* level = getLevel(patches);
  if( isItTime(old_dw, level, FL_lb->lastCompTimeLabel) == false){
    return;
  }

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch, dout_OTF_FLT, "Doing FirstLawThermo::compute_MPM_Contributions");

    //__________________________________
    //  compute the thermal energy of the solids
    int    numMPMMatls     = d_materialManager->getNumMaterials( "MPM" );
    double MPM_totalIntEng = 0.0;

    for(int m = 0; m < numMPMMatls; m++){
      MPMMaterial* mpm_matl = (MPMMaterial*) d_materialManager->getMaterial( "MPM",  m );
      int dwi = mpm_matl->getDWIndex();

      constParticleVariable<double> pMassNew;
      constParticleVariable<double> pTempNew;

      ParticleSubset* pset = old_dw->getParticleSubset( dwi, patch );
      new_dw->get( pTempNew, M_lb->pTemperatureLabel_preReloc, pset );
      new_dw->get( pMassNew, M_lb->pMassLabel_preReloc,        pset );

      double Cp = d_mpm_specificHeat[dwi];

      for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
        particleIndex idx = *iter;

        MPM_totalIntEng += pTempNew[idx] * pMassNew[idx] * Cp;
      }
    }  // matl loop

    new_dw->put( sum_vartype(MPM_totalIntEng), FL_lb->MPM_totalIntEngLabel );
  }
}

//______________________________________________________________________
//
void FirstLawThermo::doAnalysis([[maybe_unused]] const ProcessorGroup * pg,
                                const PatchSubset    * patches,
                                [[maybe_unused]] const MaterialSubset * matls ,
                                DataWarehouse        * old_dw,
                                DataWarehouse        * new_dw)
{
  const Level* level = getLevel(patches);
  timeVars tv;

  getTimeVars( old_dw, level, FL_lb->lastCompTimeLabel, tv );
  putTimeVars( new_dw, FL_lb->lastCompTimeLabel, tv );

  if( tv.isItTime == false){
    return;
  }

  //__________________________________
  //
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    printTask(patches, patch, dout_OTF_FLT, "Doing FirstLawThermo::doAnalysis");

    //__________________________________
    // open the struct that contains the file pointer map.  We use FileInfoP types
    // and store them in the DW to avoid doing system calls (SLOW).
    // Note: after regridding this may not exist for this patch in the old_dw
    PerPatch<FileInfoP> fileInfo;

    if( old_dw->exists( FL_lb->fileVarsStructLabel, 0, patch ) ){
      old_dw->get(fileInfo, FL_lb->fileVarsStructLabel, 0, patch);
    }
    else{
      FileInfo* myFileInfo = scinew FileInfo();
      fileInfo.get() = myFileInfo;
    }

    std::map<string, FILE *> myFiles;

    if( fileInfo.get().get_rep() ){
      myFiles = fileInfo.get().get_rep()->files;
    }

    string udaDir = d_output->getOutputLocation();
    string filename = udaDir + "/" + "1stLawThermo.dat";
    FILE *fp=nullptr;

    if( myFiles.count(filename) == 0 ){
      createFile(filename, fp);
      myFiles[filename] = fp;

    }
    else {
      fp = myFiles[filename];
    }

    if (!fp){
      throw InternalError("\nERROR:dataAnalysisModule:1stLawThermo:  failed opening file"+filename,__FILE__, __LINE__);
    }

    //__________________________________
    //
    sum_vartype ICE_totalIntEng, MPM_totalIntEng, total_flux;
    new_dw->get( ICE_totalIntEng, FL_lb->ICE_totalIntEngLabel );
    new_dw->get( MPM_totalIntEng, FL_lb->MPM_totalIntEngLabel );
    new_dw->get( total_flux,      FL_lb->totalFluxesLabel );

    double totalIntEng = (double)ICE_totalIntEng + (double)MPM_totalIntEng;

    fprintf(fp, "%16.15E      %16.15E      %16.15E       %16.15E       %16.15E\n", tv.now,
                (double)ICE_totalIntEng,
                (double)MPM_totalIntEng,
                totalIntEng,
                (double)total_flux );

//      fflush(fp);   If you want to write the data right now, no buffering.

    //__________________________________
    // put the file pointers into the DataWarehouse
    // these could have been altered. You must
    // reuse the Handle fileInfo and just replace the contents
    fileInfo.get().get_rep()->files = myFiles;

    new_dw->put(fileInfo, FL_lb->fileVarsStructLabel, 0, patch);
  }
}


//______________________________________________________________________
//  Open the file if it doesn't exist and write the file header
void FirstLawThermo::createFile(string& filename,  FILE*& fp)
{
  // if the file already exists then exit.  The file could exist but not be owned by this processor
  ifstream doExists( filename.c_str() );
  if(doExists){
    fp = fopen(filename.c_str(), "a");
    return;
  }

  fp = fopen(filename.c_str(), "w");
  fprintf(fp,"# This assumes:\n");
  fprintf(fp,"#    - mpm matls have constant specific heat\n");
  fprintf(fp,"#    - mpm matls are listed in order 0, 1, 2, 3\n");
  fprintf(fp,"#    - Energy conversion factor, in SI units KJ ->J %E\n",d_conversion);
  fprintf(fp,"#Time                      ICE_totalIntEng            MPM_totalIntEng             totalIntEng                 total_ICE_Flux\n");
  cout << Parallel::getMPIRank() << " FirstLawThermo:Created file " << filename << std::endl;
}


//______________________________________________________________________
//   This is a rip off of what's done int the boundary condition code
void FirstLawThermo::faceInfo(const std::string fc,
                              Patch::FaceType& face_side,
                              Vector& norm,
                              int& p_dir)
{
  if (fc ==  "x-"){
    norm = Vector(-1, 0, 0);
    p_dir = 0;
    face_side = Patch::xminus;
  }
  if (fc == "x+"){
    norm = Vector(1, 0, 0);
    p_dir = 0;
    face_side = Patch::xplus;
  }
  if (fc == "y-"){
    norm = Vector(0, -1, 0);
    p_dir = 1;
    face_side = Patch::yminus;
  }
  if (fc == "y+"){
    norm = Vector(0, 1, 0);
    p_dir = 1;
    face_side = Patch::yplus;
  }
  if (fc == "z-"){
    norm = Vector(0, 0, -1);
    p_dir = 2;
    face_side = Patch::zminus;
  }
  if (fc == "z+"){
    norm = Vector(0, 0, 1);
    p_dir = 2;
    face_side = Patch::zplus;
  }
}
//______________________________________________________________________
//  bulletProofing on the user inputs
void FirstLawThermo::bulletProofing(GridP& grid,
                                    const string& side,
                                    const Point& start,
                                    const Point& end)
{
  bulletProofing_LinesPlanes( objectType::plane, grid, "1stLawThermo", start,end );

  //__________________________________
  //  plane must be on the edge of the domain
  bool validPlane = true;
  BBox compDomain;
  grid->getInteriorSpatialRange(compDomain);
  Point min = compDomain.min();
  Point max = compDomain.max();

  Point me = min;
  if (side == "x+" || side == "y+" || side == "z+" ){
    me = max;
  }

  if(side == "x+" || side == "x-" ){
    if(start.x() != me.x() ){
      validPlane = false;
    }
  }
  if(side == "y+" || side == "y-" ){
    if(start.y() != me.y() ){
      validPlane = false;
    }
  }
  if(side == "z+" || side == "z-" ){
    if(start.z() != me.z() ){
      validPlane = false;
    }
  }
  if( validPlane == false ){
    ostringstream warn;
    warn << "\n ERROR:1stLawThermo: the plane on face ("<< side
         << ")  specified " << start << " to " << end
         << " is not at the edge of the computational domain. \n" << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

}
