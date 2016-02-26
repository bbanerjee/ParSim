/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2016 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVar_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_MasonSand.h>
#include <Core/Labels/MPMLabel.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/InternalError.h>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;

/*!-----------------------------------------------------*/
InternalVar_MasonSand::InternalVar_MasonSand(ProblemSpecP& ps,
                                             ElasticModuliModel* elastic)
{
  d_elastic = elastic;
  d_shear = 0;

  ps->require("p0", d_crushParam.p0);  // Crush Curve Parameter
  ps->require("p1", d_crushParam.p1);  // Crush Curve Parameter
  ps->require("p2", d_crushParam.p2);  // Crush Curve Parameter 
  ps->require("p3", d_crushParam.p3);  // Crush Curve Parameter

  ps->getWithDefault("use_disaggregation_algorithm", d_use_disaggregation_algorithm, false);

  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}
         
/*!-----------------------------------------------------*/
InternalVar_MasonSand::InternalVar_MasonSand(const InternalVar_MasonSand* cm)
{
  d_elastic = cm->d_elastic;
  d_shear = cm->d_shear;

  d_crushParam = cm->d_crushParam;
  d_use_disaggregation_algorithm = cm->d_use_disaggregation_algorithm;

  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}
         
/*!-----------------------------------------------------*/
InternalVar_MasonSand::~InternalVar_MasonSand()
{
  VarLabel::destroy(pKappaLabel);
  VarLabel::destroy(pKappaLabel_preReloc);

  VarLabel::destroy(pCapXLabel);
  VarLabel::destroy(pCapXLabel_preReloc);

  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);

  VarLabel::destroy(pPlasticVolStrainLabel);
  VarLabel::destroy(pPlasticVolStrainLabel_preReloc);

  VarLabel::destroy(pP3Label);
  VarLabel::destroy(pP3Label_preReloc);
}

/*!-----------------------------------------------------*/
void InternalVar_MasonSand::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type","mason_sand");

  int_var_ps->appendElement("p0", d_crushParam.p0);
  int_var_ps->appendElement("p1", d_crushParam.p1);
  int_var_ps->appendElement("p2", d_crushParam.p2);
  int_var_ps->appendElement("p3", d_crushParam.p3);

  int_var_ps->appendElement("use_disaggregation_algorithm", d_use_disaggregation_algorithm);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::addInitialComputesAndRequires(Task* task,
                                                     const MPMMaterial* matl ,
                                                     const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pKappaLabel,            matlset);
  task->computes(pCapXLabel,             matlset);
  task->computes(pPlasticStrainLabel,    matlset);
  task->computes(pPlasticVolStrainLabel, matlset);
  task->computes(pP3Label,               matlset);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::initializeInternalVariable(const Patch* patch,
                                                  const MPMMaterial* matl,
                                                  ParticleSubset* pset,
                                                  DataWarehouse* new_dw,
                                                  MPMLabel* lb,
                                                  ParameterDict& params)
{
  Uintah::constParticleVariable<double> pMass, pVolume;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pMass,   lb->pMassLabel,   pset);

  Uintah::ParticleVariable<double>  pKappa, pCapX;
  Uintah::ParticleVariable<double>  pPlasticVolStrain, pP3;
  Uintah::ParticleVariable<Matrix3> pPlasticStrain;
  new_dw->allocateAndPut(pKappa,            pKappaLabel,            pset);
  new_dw->allocateAndPut(pCapX,             pCapXLabel,             pset);
  new_dw->allocateAndPut(pPlasticStrain,    pPlasticStrainLabel,    pset);
  new_dw->allocateAndPut(pPlasticVolStrain, pPlasticVolStrainLabel, pset);
  new_dw->allocateAndPut(pP3,               pP3Label,               pset);

  double PEAKI1 = params.at("PEAKI1");
  double CR = params.at("CR");
  double phi0 = params.at("phi0");
  double Sw0 = params.at("Sw0");

  for(auto iter = pset->begin();iter != pset->end(); iter++) {
    if (d_use_disaggregation_algorithm) {
      pP3[*iter] = log(pVolume[*iter]*(matl->getInitialDensity())/pMass[*iter]);
    } else {
      pP3[*iter] = d_crushParam.p3;
    }
    pCapX[*iter] = computeX(0.0, 0.0, phi0, Sw0, params, pP3[*iter]);
    pKappa[*iter] = PEAKI1 - CR*(PEAKI1 - pCapX[*iter]); // Branch Point
    pPlasticStrain[*iter].set(0.0);
    pPlasticVolStrain[*iter] = 0.0;
  }
}

/*!-----------------------------------------------------*/
double 
InternalVar_MasonSand::computeX(const double& ev_p, 
                                const double& I1,
                                const double& phi,
                                const double& Sw,
                                ParameterDict& params,
                                const double& pP3)
{
  // Convert to bar quantities
  double ev_p_bar = -ev_p;
  double I1_bar = -I1;

  //------------Plastic strain exceeds allowable limit--------------------------
  // The plastic strain for this iteration has exceed the allowable
  // value.  X is not defined in this region, so we set it to a large
  // negative number.
  //
  // The code should never have ev_p_bar > p3, but will have ev_p_bar = p3 if the
  // porosity approaches zero (within the specified tolerance).  By setting
  // X=1e12, the material will respond as though there is no porosity.
  if (ev_p_bar > d_crushParam.p3) {
    double XX = params["Kmax"];  // Set X to equal the maximum allowable bulk modulus
    return XX;
  }

  // ------------------Plastic strain is within allowable domain------------------------
  
  // Compute elastic volumetric strain at yield
  double ev_e_yield = elasticVolStrainYield(ev_p_bar, params);

  // Compute partially saturated bulk modulus
  //double K_part_sat = d_elastic->bulkModulusParSatSand(I1_bar, ev_p_bar, phi, Sw, params);
  double K_part_sat = 0.0;

  // Compute hydrostatic strength
  double X_bar_part_sat = 3.0*K_part_sat*ev_e_yield;
  
  return -X_bar_part_sat;
}

/*!
 * --------------------------------------------------------------------------
 * Compute elastic volume strain at yield
 * --------------------------------------------------------------------------
 */
double 
InternalVar_MasonSand::elasticVolStrainYield(const double& ev_p_bar,
                                             ParameterDict& params)
{
  // Compute X(ev_p) using crush curve model for dry sand
  double X_bar_ev_p = crushCurveDrainedSandX(ev_p_bar);

  // Compute K(ev_p, I1) using bulk modulus model for dry sand
  double I1_bar = 0.5*X_bar_ev_p;
  double bulkModulus = 0.0, shearModulus = 0.0;
  //d_elastic->computeDrainedModuli(I1_bar, ev_p_bar, bulkModulus, shearModulus);

  // Compute elastic vol strain at yield
  double ev_e_yield = X_bar_ev_p/(3.0*bulkModulus);

  return ev_e_yield;
}

/*!
 * -------------------------------------------------------------
 *  Crush curve model for the drained sand
 * -------------------------------------------------------------
 */
double
InternalVar_MasonSand::crushCurveDrainedSandX(const double& ev_p_bar)
{
  double p0 = d_crushParam.p0;
  double p1 = d_crushParam.p1;
  double p2 = d_crushParam.p2;
  double p3 = d_crushParam.p3;

  double X_bar_drained = p0;
  if (ev_p_bar > 0.0) {
    double Phi0 = 1.0 - std::exp(-p3);
    double Phi = 1.0 - std::exp(-p3 + ev_p_bar);
    double term1 = (Phi0/Phi - 1.0)/p1;
    double xi_bar = std::pow(term1, 1.0/p2);
    X_bar_drained += xi_bar;
  }

  return X_bar_drained;
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl ,
                                              const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pKappaLabel,            matlset, Ghost::None);
  task->requires(Task::OldDW, pCapXLabel,             matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticStrainLabel,    matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticVolStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pP3Label,               matlset, Ghost::None);
  task->computes(pKappaLabel_preReloc,            matlset);
  task->computes(pCapXLabel_preReloc,             matlset);
  task->computes(pPlasticStrainLabel_preReloc,    matlset);
  task->computes(pPlasticVolStrainLabel_preReloc, matlset);
  task->computes(pP3Label_preReloc,               matlset);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::addParticleState(std::vector<const VarLabel*>& from,
                                        std::vector<const VarLabel*>& to)
{
  from.push_back(pKappaLabel);
  to.push_back(pKappaLabel_preReloc);

  from.push_back(pCapXLabel);
  to.push_back(pCapXLabel_preReloc);

  from.push_back(pPlasticStrainLabel);
  to.push_back(pPlasticStrainLabel_preReloc);

  from.push_back(pPlasticVolStrainLabel);
  to.push_back(pPlasticVolStrainLabel_preReloc);

  from.push_back(pP3Label);
  to.push_back(pP3Label_preReloc);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::allocateCMDataAddRequires(Task* task,
                                                 const MPMMaterial* matl ,
                                                 const PatchSet* ,
                                                 MPMLabel* )
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pKappaLabel_preReloc,            matlset, Ghost::None);
  task->requires(Task::NewDW, pCapXLabel_preReloc,             matlset, Ghost::None);
  task->requires(Task::NewDW, pPlasticStrainLabel_preReloc,    matlset, Ghost::None);
  task->requires(Task::NewDW, pPlasticVolStrainLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pP3Label_preReloc,               matlset, Ghost::None);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::allocateCMDataAdd(DataWarehouse* old_dw,
                                         ParticleSubset* addset,
                                         ParticleLabelVariableMap* newState,
                                         ParticleSubset* delset,
                                         DataWarehouse* new_dw )
{
  ParticleVariable<double> pKappa, pCapX, pEpv, pP3;
  ParticleVariable<Matrix3> pEp;
  constParticleVariable<double> o_kappa, o_capX, o_Epv, o_P3;
  constParticleVariable<Matrix3> o_Ep;

  new_dw->allocateTemporary(pKappa, addset);
  new_dw->allocateTemporary(pCapX,  addset);
  new_dw->allocateTemporary(pEp,    addset);
  new_dw->allocateTemporary(pEpv,   addset);
  new_dw->allocateTemporary(pP3,    addset);

  new_dw->get(o_kappa, pKappaLabel_preReloc,            delset);
  new_dw->get(o_capX,  pCapXLabel_preReloc,             delset);
  new_dw->get(o_Ep,    pPlasticStrainLabel_preReloc,    delset);
  new_dw->get(o_Epv,   pPlasticVolStrainLabel_preReloc, delset);
  new_dw->get(o_P3,    pP3Label_preReloc,               delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for(o = delset->begin(); o != delset->end(); o++, n++) {
    pKappa[*n] = o_kappa[*o];
    pCapX[*n]  = o_capX[*o];
    pEp[*n]    = o_Ep[*o];
    pEpv[*n]   = o_Epv[*o];
    pP3[*n]    = o_P3[*o];
  }

  (*newState)[pKappaLabel]            = pKappa.clone();
  (*newState)[pCapXLabel]             = pCapX.clone();
  (*newState)[pPlasticStrainLabel]    = pEp.clone();
  (*newState)[pPlasticVolStrainLabel] = pEpv.clone();
  (*newState)[pP3Label]               = pP3.clone();
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::getInternalVariable(ParticleSubset* pset ,
                                           DataWarehouse* old_dw,
                                           constParticleLabelVariableMap& var)
{
  constParticleVariable<double>  pKappa, pCapX, pEpv, pP3;
  constParticleVariable<Matrix3> pEp;
  old_dw->get(pKappa, pKappaLabel,            pset);
  old_dw->get(pCapX,  pCapXLabel,             pset);
  old_dw->get(pEp,    pPlasticStrainLabel,    pset);
  old_dw->get(pEpv,   pPlasticVolStrainLabel, pset);
  old_dw->get(pP3,    pP3Label,               pset);

  var[pKappaLabel]            = &pKappa;
  var[pCapXLabel]             = &pCapX;
  var[pPlasticStrainLabel]    = &pEp;
  var[pPlasticVolStrainLabel] = &pEpv;
  var[pP3Label]               = &pP3;
}

void 
InternalVar_MasonSand::allocateAndPutInternalVariable(ParticleSubset* pset,
                                                      DataWarehouse* new_dw,
                                                      ParticleLabelVariableMap& var_new) 
{
  new_dw->allocateAndPut(*var_new[pKappaLabel],            pKappaLabel_preReloc,            pset);
  new_dw->allocateAndPut(*var_new[pCapXLabel],             pCapXLabel_preReloc,             pset);
  new_dw->allocateAndPut(*var_new[pPlasticStrainLabel],    pPlasticStrainLabel_preReloc,    pset);
  new_dw->allocateAndPut(*var_new[pPlasticVolStrainLabel], pPlasticVolStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(*var_new[pP3Label],               pP3Label_preReloc,               pset);
}

void
InternalVar_MasonSand::allocateAndPutRigid(ParticleSubset* pset,
                                           DataWarehouse* new_dw,
                                           constParticleLabelVariableMap& var)
{
  ParticleVariable<double> pKappa_new, pCapX_new, pEpv_new, pP3_new;
  ParticleVariable<Matrix3> pEp_new;
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc,            pset);
  new_dw->allocateAndPut(pCapX_new,  pCapXLabel_preReloc,             pset);
  new_dw->allocateAndPut(pEp_new,    pPlasticStrainLabel_preReloc,    pset);
  new_dw->allocateAndPut(pEpv_new,   pPlasticVolStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(pP3_new,    pP3Label_preReloc,               pset);
  for(auto iter = pset->begin(); iter != pset->end(); iter++){
     pKappa_new[*iter] = 
       dynamic_cast<constParticleVariable<double>& >(*var[pKappaLabel])[*iter];
     pCapX_new[*iter]  = 
       dynamic_cast<constParticleVariable<double>& >(*var[pCapXLabel])[*iter];
     pEp_new[*iter]    = 
       dynamic_cast<constParticleVariable<Matrix3>& >(*var[pPlasticStrainLabel])[*iter];
     pEpv_new[*iter]    = 
       dynamic_cast<constParticleVariable<double>& >(*var[pPlasticVolStrainLabel])[*iter];
     pP3_new[*iter]    = 
       dynamic_cast<constParticleVariable<double>& >(*var[pP3Label])[*iter];
  }
}

//--------------------------------------------------------------------------------------
// Compute kappa_new using Newton's method
//--------------------------------------------------------------------------------------
double 
InternalVar_MasonSand::computeInternalVariable(const ModelStateBase* state_input) const
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  // return the new kappa
  double kappa_new = 0.0;
  return kappa_new;
}

