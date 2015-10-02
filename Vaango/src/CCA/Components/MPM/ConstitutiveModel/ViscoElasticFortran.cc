/*
 * The MIT License
 *
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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


#include <CCA/Components/MPM/ConstitutiveModel/ViscoElasticFortran.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h> 
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Math/MinMax.h>

#include <sci_defs/uintah_defs.h>

#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
// The following functions are found in fortran/visco.F90

extern "C"
{

#  define PROPCHECK propcheck_
#  define VISCOINI viscoini_
#  define VISCORELAX viscorelax_

  void PROPCHECK(int* nprop, double* props);
  void VISCOINI(int* nprop, double* props, int* nstatev, double* statev);
  void VISCORELAX(double* dtime, double* time, double* tempold, double* dtemp,
                  int* nprop, double* props, double* F, int* nstatev, double* statev,
                  double* sigo, double* sig, double* cfac);
}

// End fortran functions.
////////////////////////////////////////////////////////////////////////////////
  

using std::cerr;
using namespace Vaango;
using namespace Uintah;

ViscoElasticFortran::ViscoElasticFortran( ProblemSpecP& ps,
                                          MPMFlags* Mflag ) 
  : ConstitutiveModel(Mflag)
{
  ps->require("K", d_param.K);
  ps->require("G00", d_param.G00);
  ps->require("G01", d_param.G01);
  ps->require("G02", d_param.G02);
  ps->require("G03", d_param.G03);
  ps->require("G04", d_param.G04);
  ps->getWithDefault("G05", d_param.G05, 0.0);
  ps->getWithDefault("G06", d_param.G06, 0.0);
  ps->getWithDefault("G07", d_param.G07, 0.0);
  ps->getWithDefault("G08", d_param.G08, 0.0);
  ps->getWithDefault("G09", d_param.G09, 0.0);
  ps->getWithDefault("G10", d_param.G10, 0.0);
  ps->require("Tau00", d_param.Tau00);
  ps->require("Tau01", d_param.Tau01);
  ps->require("Tau02", d_param.Tau02);
  ps->require("Tau03", d_param.Tau03);
  ps->require("Tau04", d_param.Tau04);
  ps->getWithDefault("Tau05", d_param.Tau05, 0.0);
  ps->getWithDefault("Tau06", d_param.Tau06, 0.0);
  ps->getWithDefault("Tau07", d_param.Tau07, 0.0);
  ps->getWithDefault("Tau08", d_param.Tau08, 0.0);
  ps->getWithDefault("Tau09", d_param.Tau09, 0.0);
  ps->getWithDefault("Tau10", d_param.Tau10, 0.0);
  ps->getWithDefault("C1_WLF", d_param.C1_WLF, 0.0);
  ps->getWithDefault("C2_WLF", d_param.C2_WLF, 0.0);
  ps->getWithDefault("Tref_WLF", d_param.Tref_WLF, 300.0);

  /* Check inputs */
  d_nProp = 24;
  std::vector<double> d_props(d_nProp);
  d_props.push_back(d_param.C1_WLF);
  d_props.push_back(d_param.C2_WLF);
  d_props.push_back(d_param.Tref_WLF);
  d_props.push_back(d_param.G00);
  d_props.push_back(d_param.G01);
  d_props.push_back(d_param.G02);
  d_props.push_back(d_param.G03);
  d_props.push_back(d_param.G04);
  d_props.push_back(d_param.G05);
  d_props.push_back(d_param.G06);
  d_props.push_back(d_param.G07);
  d_props.push_back(d_param.G08);
  d_props.push_back(d_param.G09);
  d_props.push_back(d_param.G10);
  d_props.push_back(d_param.Tau01);
  d_props.push_back(d_param.Tau02);
  d_props.push_back(d_param.Tau03);
  d_props.push_back(d_param.Tau04);
  d_props.push_back(d_param.Tau05);
  d_props.push_back(d_param.Tau06);
  d_props.push_back(d_param.Tau07);
  d_props.push_back(d_param.Tau08);
  d_props.push_back(d_param.Tau09);
  d_props.push_back(d_param.Tau10);

  PROPCHECK(&d_nProp, &d_props[0]);

  d_nStateV = 68;
  initializeLocalMPMLabels();
}

ViscoElasticFortran::ViscoElasticFortran(const ViscoElasticFortran* cm) 
  : ConstitutiveModel(cm)
{
  d_nProp = cm->d_nProp;
  d_props = cm->d_props;

  d_param.K = cm->d_param.K;
  d_param.G00 = cm->d_param.G00;

  d_param.G01 = cm->d_param.G01;
  d_param.G02 = cm->d_param.G02;
  d_param.G03 = cm->d_param.G03;
  d_param.G04 = cm->d_param.G04;
  d_param.G05 = cm->d_param.G05;
  d_param.G06 = cm->d_param.G06;
  d_param.G07 = cm->d_param.G07;
  d_param.G08 = cm->d_param.G08;
  d_param.G09 = cm->d_param.G09;
  d_param.G10 = cm->d_param.G10;

  d_param.Tau01 = cm->d_param.Tau01;
  d_param.Tau02 = cm->d_param.Tau02;
  d_param.Tau03 = cm->d_param.Tau03;
  d_param.Tau04 = cm->d_param.Tau04;
  d_param.Tau05 = cm->d_param.Tau05;
  d_param.Tau06 = cm->d_param.Tau06;
  d_param.Tau07 = cm->d_param.Tau07;
  d_param.Tau08 = cm->d_param.Tau08;
  d_param.Tau09 = cm->d_param.Tau09;
  d_param.Tau10 = cm->d_param.Tau10;

  d_param.C1_WLF = cm->d_param.C1_WLF;
  d_param.C2_WLF = cm->d_param.C2_WLF;
  d_param.Tref_WLF = cm->d_param.Tref_WLF;

  d_nStateV = cm->d_nStateV;
  initializeLocalMPMLabels();
}

ViscoElasticFortran::~ViscoElasticFortran()
{
  for (unsigned int i = 0; i < stateVLabels.size(); i++) {
     VarLabel::destroy(stateVLabels[i]);
  }
}

void
ViscoElasticFortran::outputProblemSpec( ProblemSpecP& ps,bool output_cm_tag )
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type","viscoelastic_fortran");
  }

  cm_ps->appendElement("K", d_param.K);
  cm_ps->appendElement("G00", d_param.G00);
  cm_ps->appendElement("G01", d_param.G01);
  cm_ps->appendElement("G02", d_param.G02);
  cm_ps->appendElement("G03", d_param.G03);
  cm_ps->appendElement("G04", d_param.G04);
  cm_ps->appendElement("G05", d_param.G05);
  cm_ps->appendElement("G06", d_param.G06);
  cm_ps->appendElement("G07", d_param.G07);
  cm_ps->appendElement("G08", d_param.G08);
  cm_ps->appendElement("G09", d_param.G09);
  cm_ps->appendElement("G10", d_param.G10);
  cm_ps->appendElement("Tau01", d_param.Tau01);
  cm_ps->appendElement("Tau02", d_param.Tau02);
  cm_ps->appendElement("Tau03", d_param.Tau03);
  cm_ps->appendElement("Tau04", d_param.Tau04);
  cm_ps->appendElement("Tau05", d_param.Tau05);
  cm_ps->appendElement("Tau06", d_param.Tau06);
  cm_ps->appendElement("Tau07", d_param.Tau07);
  cm_ps->appendElement("Tau08", d_param.Tau08);
  cm_ps->appendElement("Tau09", d_param.Tau09);
  cm_ps->appendElement("Tau10", d_param.Tau10);
  cm_ps->appendElement("C1_WLF", d_param.C1_WLF);
  cm_ps->appendElement("C2_WLF", d_param.C2_WLF);
  cm_ps->appendElement("Tref_WLF", d_param.Tref_WLF);
}

ViscoElasticFortran*
ViscoElasticFortran::clone()
{
  return scinew ViscoElasticFortran(*this);
}

void
ViscoElasticFortran::initializeLocalMPMLabels()
{
  std::vector<std::string> stateVNames(d_nStateV);
  stateVNames.push_back("WLF_AV_SHIFT");
  stateVNames.push_back("NUM_AV_SHIFT");

  stateVNames.push_back("DEVPK2_0_XX");
  stateVNames.push_back("DEVPK2_0_YY");
  stateVNames.push_back("DEVPK2_0_ZZ");
  stateVNames.push_back("DEVPK2_0_XY");
  stateVNames.push_back("DEVPK2_0_YZ");
  stateVNames.push_back("DEVPK2_0_ZX");

  stateVNames.push_back("DEVPK2_XX_1");
  stateVNames.push_back("DEVPK2_YY_1");
  stateVNames.push_back("DEVPK2_ZZ_1");
  stateVNames.push_back("DEVPK2_XY_1");
  stateVNames.push_back("DEVPK2_YZ_1");
  stateVNames.push_back("DEVPK2_ZX_1");

  stateVNames.push_back("DEVPK2_XX_2");
  stateVNames.push_back("DEVPK2_YY_2");
  stateVNames.push_back("DEVPK2_ZZ_2");
  stateVNames.push_back("DEVPK2_XY_2");
  stateVNames.push_back("DEVPK2_YZ_2");
  stateVNames.push_back("DEVPK2_ZX_2");

  stateVNames.push_back("DEVPK2_XX_3");
  stateVNames.push_back("DEVPK2_YY_3");
  stateVNames.push_back("DEVPK2_ZZ_3");
  stateVNames.push_back("DEVPK2_XY_3");
  stateVNames.push_back("DEVPK2_YZ_3");
  stateVNames.push_back("DEVPK2_ZX_3");

  stateVNames.push_back("DEVPK2_XX_4");
  stateVNames.push_back("DEVPK2_YY_4");
  stateVNames.push_back("DEVPK2_ZZ_4");
  stateVNames.push_back("DEVPK2_XY_4");
  stateVNames.push_back("DEVPK2_YZ_4");
  stateVNames.push_back("DEVPK2_ZX_4");

  stateVNames.push_back("DEVPK2_XX_5");
  stateVNames.push_back("DEVPK2_YY_5");
  stateVNames.push_back("DEVPK2_ZZ_5");
  stateVNames.push_back("DEVPK2_XY_5");
  stateVNames.push_back("DEVPK2_YZ_5");
  stateVNames.push_back("DEVPK2_ZX_5");

  stateVNames.push_back("DEVPK2_XX_6");
  stateVNames.push_back("DEVPK2_YY_6");
  stateVNames.push_back("DEVPK2_ZZ_6");
  stateVNames.push_back("DEVPK2_XY_6");
  stateVNames.push_back("DEVPK2_YZ_6");
  stateVNames.push_back("DEVPK2_ZX_6");

  stateVNames.push_back("DEVPK2_XX_7");
  stateVNames.push_back("DEVPK2_YY_7");
  stateVNames.push_back("DEVPK2_ZZ_7");
  stateVNames.push_back("DEVPK2_XY_7");
  stateVNames.push_back("DEVPK2_YZ_7");
  stateVNames.push_back("DEVPK2_ZX_7");

  stateVNames.push_back("DEVPK2_XX_8");
  stateVNames.push_back("DEVPK2_YY_8");
  stateVNames.push_back("DEVPK2_ZZ_8");
  stateVNames.push_back("DEVPK2_XY_8");
  stateVNames.push_back("DEVPK2_YZ_8");
  stateVNames.push_back("DEVPK2_ZX_8");

  stateVNames.push_back("DEVPK2_XX_9");
  stateVNames.push_back("DEVPK2_YY_9");
  stateVNames.push_back("DEVPK2_ZZ_9");
  stateVNames.push_back("DEVPK2_XY_9");
  stateVNames.push_back("DEVPK2_YZ_9");
  stateVNames.push_back("DEVPK2_ZX_9");

  stateVNames.push_back("DEVPK2_XX_10");
  stateVNames.push_back("DEVPK2_YY_10");
  stateVNames.push_back("DEVPK2_ZZ_10");
  stateVNames.push_back("DEVPK2_XY_10");
  stateVNames.push_back("DEVPK2_YZ_10");
  stateVNames.push_back("DEVPK2_ZX_10");

  for (int i = 0; i < d_nStateV; i++) {
    stateVLabels.push_back(VarLabel::create(stateVNames[i],
                           ParticleVariable<double>::getTypeDescription()));
    stateVLabels_preReloc.push_back(VarLabel::create(stateVNames[i]+"+",
                                    ParticleVariable<double>::getTypeDescription()));
  }
}

void
ViscoElasticFortran::initializeCMData( const Patch* patch,
				       const MPMMaterial* matl,
				       DataWarehouse* new_dw )
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  double statev[d_nStateV];
  VISCOINI(d_nProp, &d_props[0], &d_nStateV, &stateV);

  StaticArray<ParticleVariable<double> > stateV(d_nStateV+1);
  for (int i=0; i < d_nStateV; i++) {
    new_dw->allocateAndPut(stateV[i], stateVLabels[i], pset);
    for(auto iter = pset->begin(); iter != pset->end(); iter++){
      stateV[i][*iter] = statev[i];
    }
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
ViscoElasticFortran::allocateCMDataAddRequires( Task* task,
						const MPMMaterial* matl,
						const PatchSet* patches,
						MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();

  // Allocate the variables shared by all constitutive models
  // for the particle convert operation
  // This method is defined in the ConstitutiveModel base class.
  addSharedRForConvertExplicit(task, matlset, patches);

  // Add requires local to this model
  for (int i=0; i < d_nStateV; i++) {
    task->requires(Task::NewDW, stateVLabels_preReloc[i], matlset, Ghost::None);
  }
}

void
ViscoElasticFortran::allocateCMDataAdd( DataWarehouse* new_dw,
					ParticleSubset* addset,
					ParticleLabelVariableMap* newState,
					ParticleSubset* delset,
					DataWarehouse* )
{
  // Copy the data common to all constitutive models from the particle to be 
  // deleted to the particle to be added. 
  // This method is defined in the ConstitutiveModel base class.
  copyDelToAddSetForConvertExplicit(new_dw, delset, addset, newState);

  StaticArray<ParticleVariable<double> > stateV(d_nStateV+1);
  StaticArray<constParticleVariable<double> > o_stateV(d_nStateV+1);
  for(int i=0; i < d_nStateV; i++){
    new_dw->allocateTemporary(stateV[i], addset);
    new_dw->get(o_stateV[i],stateVLabels_preReloc[i], delset);

    ParticleSubset::iterator o,n = addset->begin();
    for (o=delset->begin(); o != delset->end(); o++, n++) {
      stateV[i][*n] = o_stateV[i][*n];
    }
    (*newState)[stateVLabels[i]]=stateV[i].clone();
  }
}

void
ViscoElasticFortran::addParticleState(std::vector<const VarLabel*>& from,
				      std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  for(int i=0;i<d_nStateV;i++){
    from.push_back(stateVLabels[i]);
    to.push_back(stateVLabels_preReloc[i]);
  }
}

void
ViscoElasticFortran::computeStableTimestep( const Patch* patch,
					    const MPMMaterial* matl,
					    DataWarehouse* new_dw )
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int matID = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);
  constParticleVariable<double> pmass, pvolume;
  constParticleVariable<Vector> pvelocity;

  new_dw->get(pmass,     lb->pMassLabel,     pset);
  new_dw->get(pvolume,   lb->pVolumeLabel,   pset);
  new_dw->get(pvelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector waveSpeed(1.e-12,1.e-12,1.e-12);

  double G = d_param.G;
  double bulk = d_param.K;
  for(ParticleSubset::iterator iter = pset->begin();iter != pset->end();iter++){
    particleIndex idx = *iter;

    // Compute wave speed at each particle, store the maximum
    c_dil = sqrt((bulk + 4.*G/3.)*pvolume[idx]/pmass[idx]);
    waveSpeed=Vector(Max(c_dil+fabs(pvelocity[idx].x()),waveSpeed.x()),
		     Max(c_dil+fabs(pvelocity[idx].y()),waveSpeed.y()),
		     Max(c_dil+fabs(pvelocity[idx].z()),waveSpeed.z()));
  }
  waveSpeed = dx/waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
ViscoElasticFortran::computeStressTensor( const PatchSubset* patches,
					  const MPMMaterial* matl,
					  DataWarehouse* old_dw,
					  DataWarehouse* new_dw )
{
  double time = d_sharedState->getElapsedTime();
  double rho_orig = matl->getInitialDensity();
  Matrix3 Identity, zero(0.), One(1.);
  Identity.Identity();
  // double onethird = (1.0/3.0);

  std::vector<double> d_props(d_nProp);
  d_props.push_back(d_param.C1_WLF);
  d_props.push_back(d_param.C2_WLF);
  d_props.push_back(d_param.Tref_WLF);
  d_props.push_back(d_param.G00);
  d_props.push_back(d_param.G01);
  d_props.push_back(d_param.G02);
  d_props.push_back(d_param.G03);
  d_props.push_back(d_param.G04);
  d_props.push_back(d_param.G05);
  d_props.push_back(d_param.G06);
  d_props.push_back(d_param.G07);
  d_props.push_back(d_param.G08);
  d_props.push_back(d_param.G09);
  d_props.push_back(d_param.G10);
  d_props.push_back(d_param.Tau01);
  d_props.push_back(d_param.Tau02);
  d_props.push_back(d_param.Tau03);
  d_props.push_back(d_param.Tau04);
  d_props.push_back(d_param.Tau05);
  d_props.push_back(d_param.Tau06);
  d_props.push_back(d_param.Tau07);
  d_props.push_back(d_param.Tau08);
  d_props.push_back(d_param.Tau09);
  d_props.push_back(d_param.Tau10);

  double se = 0.0;
  Vector waveSpeed(1.e-12,1.e-12,1.e-12);
  for(int p=0;p<patches->size();p++){

    const Patch* patch = patches->get(p);

    Vector dx = patch->dCell();
    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Get the deformation gradient (F) and velocity gradient (L)
    constParticleVariable<Matrix3>  pDefGrad_old, pDefGrad_new;
    constParticleVariable<Matrix3>  pVelGrad_old, pVelGrad_new;
    old_dw->get(pDefGrad_old, lb->pDefGradLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);
    old_dw->get(pVelGrad_old, lb->pVelGradLabel, pset);
    new_dw->get(pVelGrad_new, lb->pVelGradLabel_preReloc, pset);

    // Get the particle location, particle size, particle mass, particle volume
    constParticleVariable<Point>  px;
    constParticleVariable<Matrix3> psize;
    constParticleVariable<double> pMass, pVol_new, pTemp;
    old_dw->get(px,       lb->pXLabel,           pset);
    old_dw->get(psize,    lb->pSizeLabel,        pset);
    old_dw->get(pMass,    lb->pMassLabel,        pset);
    old_dw->get(pTemp,    lb->pTemperatureLabel, pset);
    new_dw->get(pVol_new, lb->pVolumeLabel_preReloc, pset);

    // Get the particle velocity
    constParticleVariable<Vector> pVelocity;
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);

    // Get the particle stress 
    constParticleVariable<Matrix3> pStress_old;
    old_dw->get(pStress_old, lb->pStressLabel,       pset);

    // Get the time increment (delT)
    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    // Get the state variables
    StaticArray<constParticleVariable<double> > stateV(d_nStateV+1);
    for (int i=0; i < d_nStateV; i++) {
      old_dw->get(stateV[i], stateVLabels[i], pset);
    }

    // Create and allocate arrays for storing the updated information
    ParticleVariable<Matrix3> pStress_new;
    ParticleVariable<double>  pdTdt, p_q;
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt,       lb->pdTdtLabel_preReloc,   pset);
    new_dw->allocateAndPut(p_q,         lb->p_qLabel_preReloc,     pset);

    StaticArray<ParticleVariable<double> > stateV_new(d_nStateV+1);
    for(int i=0; i < d_nStateV; i++){
      new_dw->allocateAndPut(stateV_new[i], stateVLabels_preReloc[i], pset);
    }

    for (auto iter = pset->begin(); iter != pset->end(); iter++) {
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      //-----------------------------------------------------------------------
      // Stage 1: Compute rate of deformation
      //-----------------------------------------------------------------------
      Matrix3 defGrad_new = pDefGrad_new[idx];
      double J_new = defGrad_new.Determinant();
      Matrix3 velGrad_new = pVelGrad_new[idx];
      Matrix3 rateOfDef_new = (velGrad_new + velGrad_new.Transpose())*0.5;

      // Calculate the current mass density and deformed volume
      double rho_cur = rho_0/J_new;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      Matrix3 D = (velGrad + velGrad.Transpose())*.5;

      // Relax the stress (use visco.F90)
      double dt = delT;
      double temp_new = pTemp[idx];
      double dtemp = 0.0;

      std::vector<double> F(3*3);
      F.push_back(defGrad_new(0,0));
      F.push_back(defGrad_new(1,0));
      F.push_back(defGrad_new(2,0));
      F.push_back(defGrad_new(0,1));
      F.push_back(defGrad_new(1,1));
      F.push_back(defGrad_new(2,1));
      F.push_back(defGrad_new(0,2));
      F.push_back(defGrad_new(1,2));
      F.push_back(defGrad_new(2,2));

      double statev[d_nStateV];
      for (int i = 0; i < d_nStateV; i++) {
        statev[i] = stateV[i][idx];
      }

      double sig_old[6] = {0}; 
      sig_old[0]=pStress_old[idx](0,0);
      sig_old[1]=pStress_old[idx](1,1);
      sig_old[2]=pStress_old[idx](2,2);
      sig_old[3]=pStress_old[idx](0,1);
      sig_old[4]=pStress_old[idx](1,2);
      sig_old[5]=pStress_old[idx](2,0);

      double sig_new[6] = {0};
      VISCORELAX(&dt, &time, &temp_new, &dtemp,
                  &d_nProp, &d_props[0], &F[0], &d_nStateV, &statev[0],
                  &sig_old[0], &sig_new[0], double* cfac);

      pStress_new[idx](0,0) = sigarg[0];
      pStress_new[idx](1,1) = sigarg[1];
      pStress_new[idx](2,2) = sigarg[2];
      pStress_new[idx](0,1) = sigarg[3];
      pStress_new[idx](1,0) = sigarg[3];
      pStress_new[idx](2,1) = sigarg[4];
      pStress_new[idx](1,2) = sigarg[4];
      pStress_new[idx](2,0) = sigarg[5];
      pStress_new[idx](0,2) = sigarg[5];

#if 0
      cout << pStress_new[idx] << endl;
#endif

      c_dil = sqrt(USM/rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pStress_new[idx] + pStress[idx])*.5;

      double e = (D(0,0)*AvgStress(0,0) +
                  D(1,1)*AvgStress(1,1) +
                  D(2,2)*AvgStress(2,2) +
		  2.*(D(0,1)*AvgStress(0,1) +
		      D(0,2)*AvgStress(0,2) +
		      D(1,2)*AvgStress(1,2))) * pvolume_new[idx]*delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pvelocity_idx = pvelocity[idx];
      waveSpeed=Vector(Max(c_dil+fabs(pvelocity_idx.x()),waveSpeed.x()),
                       Max(c_dil+fabs(pvelocity_idx.y()),waveSpeed.y()),
                       Max(c_dil+fabs(pvelocity_idx.z()),waveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z())/3.0;
        double c_bulk = sqrt(UI[0]/rho_cur);
        Matrix3 D=(velGrad + velGrad.Transpose())*0.5;
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    }  // end loop over particles

    waveSpeed = dx/waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se),     lb->StrainEnergyLabel);
    }
  }
}

void
ViscoElasticFortran::carryForward( const PatchSubset* patches,
				   const MPMMaterial* matl,
				   DataWarehouse* old_dw,
				   DataWarehouse* new_dw )
{
  for(int p=0; p<patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Carry forward the data common to all constitutive models 
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    StaticArray<constParticleVariable<double> > stateV(d_nStateV+1);
    StaticArray<ParticleVariable<double> > stateV_new(d_nStateV+1);
    for(int i=0;i<d_nStateV;i++){
      old_dw->get(stateV[i],stateVLabels[i], pset);
      new_dw->allocateAndPut(stateV_new[i],stateVLabels_preReloc[i], pset);
      stateV_new[i].copyData(stateV[i]);
    }

    // Carry forward the data local to this constitutive model 
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());
    
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.),     lb->StrainEnergyLabel);
    }
  }
}


void 
ViscoElasticFortran::addInitialComputesAndRequires(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet* ) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  cout << "In add InitialComputesAnd" << endl;

  // Other constitutive model and input dependent computes and requires
  for(int i=0;i<d_nStateV;i++){
    task->computes(stateVLabels[i], matlset);
  }
}

void
ViscoElasticFortran::addComputesAndRequires( Task* task,
					     const MPMMaterial* matl,
					     const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit 
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Computes and requires for internal state data
  for(int i=0;i<d_nStateV;i++){
    task->requires(Task::OldDW, stateVLabels[i],          matlset, Ghost::None);
    task->computes(             stateVLabels_preReloc[i], matlset);
  }
}

void
ViscoElasticFortran::addComputesAndRequires( Task*,
					     const MPMMaterial*,
					     const PatchSet*,
					     const bool, 
					     const bool ) const
{
}

double
ViscoElasticFortran::computeRhoMicroCM( double pressure,
					const double p_ref,
					const MPMMaterial* matl, 
					double temperature,
					double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge = pressure - p_ref;
  double rho_cur;
  double bulk = d_param.K;

  rho_cur = rho_orig/(1-p_gauge/bulk);

  return rho_cur;

#if 0
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR ViscoElasticFortran"
       << endl;
#endif
}

void
ViscoElasticFortran::computePressEOSCM( double rho_cur, double& pressure,
					double p_ref,
					double& dp_drho,      double& tmp,
					const MPMMaterial* matl, 
					double temperature )
{
  double bulk = d_param.K;
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk*(1.0 - rho_orig/rho_cur);
  pressure = p_ref + p_g;
  dp_drho  = bulk*rho_orig/(rho_cur*rho_cur);
  tmp = bulk/rho_cur;  // speed of sound squared

#if 0
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR ViscoElasticFortran"
       << endl;
#endif
}

double
ViscoElasticFortran::getCompressibility()
{
  return 1.0/d_param.K;
}

