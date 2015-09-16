/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVar_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
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
InternalVar_MasonSand::InternalVar_MasonSand(ProblemSpecP& ps)
{
  ps->require("p0", d_crushParam.p0);  // Crush Curve Parameter
  ps->require("p1", d_crushParam.p1);  // Crush Curve Parameter
  ps->require("p2", d_crushParam.p2);  // Crush Curve Parameter 
  ps->require("p3", d_crushParam.p3);  // Crush Curve Parameter

  ps->require("initial_porosity",   d_fluidParam.phi0);  // Initial porosity
  ps->require("initial_saturation", d_fluidParam.S0);    // Initial water saturation

  // Initialize internal variable labels for evolution
  pKappaLabel = VarLabel::create("p.kappa",
        ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc = VarLabel::create("p.kappa+",
        ParticleVariable<double>::getTypeDescription());

  pCapXLabel = VarLabel::create("p.CapX",
        ParticleVariable<double>::getTypeDescription());
  pCapXLabel_preReloc = VarLabel::create("p.CapX+",
        ParticleVariable<double>::getTypeDescription());

  pPorosityLabel = VarLabel::create("p.phi",
        ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = VarLabel::create("p.phi+",
        ParticleVariable<double>::getTypeDescription());

  pSaturationLabel = VarLabel::create("p.Sw",
        ParticleVariable<double>::getTypeDescription());
  pSaturationLabel_preReloc = VarLabel::create("p.Sw+",
        ParticleVariable<double>::getTypeDescription());
}
         
/*!-----------------------------------------------------*/
InternalVar_MasonSand::InternalVar_MasonSand(const InternalVar_MasonSand* cm)
{
  d_crushParam = cm->d_crushParam;

  // Initialize internal variable labels for evolution
  pKappaLabel = VarLabel::create("p.kappa",
        ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc = VarLabel::create("p.kappa+",
        ParticleVariable<double>::getTypeDescription());

  pCapXLabel = VarLabel::create("p.CapX",
        ParticleVariable<double>::getTypeDescription());
  pCapXLabel_preReloc = VarLabel::create("p.CapX+",
        ParticleVariable<double>::getTypeDescription());

  pPorosityLabel = VarLabel::create("p.phi",
        ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = VarLabel::create("p.phi+",
        ParticleVariable<double>::getTypeDescription());

  pSaturationLabel = VarLabel::create("p.Sw",
        ParticleVariable<double>::getTypeDescription());
  pSaturationLabel_preReloc = VarLabel::create("p.Sw+",
        ParticleVariable<double>::getTypeDescription());
}
         
/*!-----------------------------------------------------*/
InternalVar_MasonSand::~InternalVar_MasonSand()
{
  VarLabel::destroy(pKappaLabel);
  VarLabel::destroy(pKappaLabel_preReloc);

  VarLabel::destroy(pCapXLabel);
  VarLabel::destroy(pCapXLabel_preReloc);

  VarLabel::destroy(pPorosityLabel);
  VarLabel::destroy(pPorosityLabel_preReloc);

  VarLabel::destroy(pSaturationLabel);
  VarLabel::destroy(pSaturationLabel_preReloc);
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

  int_var_ps->appendElement("initial_porosity", d_fluidParam.phi0);
  int_var_ps->appendElement("initial_saturation", d_fluidParam.S0);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::addInitialComputesAndRequires(Task* task,
                                                     const MPMMaterial* matl ,
                                                     const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pKappaLabel, matlset);
  task->computes(pCapXLabel, matlset);
  task->computes(pPorosityLabel, matlset);
  task->computes(pSaturationLabel, matlset);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::initializeInternalVariable(ParticleSubset* pset,
                                                  DataWarehouse* new_dw,
                                                  ParameterDict& params)
{
  Uintah::ParticleVariable<double> pKappa, pCapX, pPorosity, pSaturation;
  new_dw->allocateAndPut(pKappa,      pKappaLabel,      pset);
  new_dw->allocateAndPut(pCapX,       pCapXLabel,       pset);
  new_dw->allocateAndPut(pPorosity,   pPorosityLabel,   pset);
  new_dw->allocateAndPut(pSaturation, pSaturationLabel, pset);

  double PEAKI1 = params["PEAKI1"];
  double CR = params["CR"];

  for(auto iter = pset->begin();iter != pset->end(); iter++) {
    pCapX[*iter] = computeX(0.0, 0.0, d_fluidParam.phi0, d_fluidParam.S0, params);
    pKappa[*iter] = PEAKI1 - CR*(PEAKI1 - pCapX[*iter]); // Branch Point
    pPorosity[*iter] = d_fluidParam.phi0;
    pSaturation[*iter] = d_fluidParam.S0;
  }
}

/*!-----------------------------------------------------*/
double 
InternalVar_MasonSand::computeX(const double& ev_p, 
                                const double& I1,
                                const double& phi,
                                const double& Sw,
                                ParameterDict& params)
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
  double K_part_sat = bulkModulusParSatSand(I1_bar, ev_p_bar, phi, Sw, params);

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
  double K_I1_ev_p = bulkModulusDrainedSand(I1_bar, ev_p_bar, params);

  // Compute elastic vol strain at yield
  double ev_e_yield = X_bar_ev_p/(3.0*K_I1_ev_p);

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

/*!
 * -------------------------------------------------------------
 *  Efective bulk modulus model for partially saturated soil
 * -------------------------------------------------------------
 */
double
InternalVar_MasonSand::bulkModulusParSatSand(const double& I1_bar, 
                                             const double& ev_p_bar,
                                             const double& phi,
                                             const double& S_w,
                                             ParameterDict& params) 
{
  double K_a = bulkModulusAir(I1_bar, params);
  double K_w = bulkModulusWater(I1_bar, params);
  double K_d = bulkModulusDrainedSand(I1_bar, ev_p_bar, params);
  double K_s = 1.0e10; // K_max

  double K_f = 1.0/(S_w/K_w + (1.0-S_w)/K_a);

  double numer = (1.0 - K_d/K_s)*(1.0 - K_d/K_s);
  double denom = 1.0/K_s*(1.0 - K_d/K_s) + phi*(1.0/K_f - 1.0/K_s);
  double K_eff = K_d + numer/denom;

  return K_eff;
}

/*!
 * -------------------------------------------------------------
 *  Bulk modulus model for air
 * -------------------------------------------------------------
 */
double
InternalVar_MasonSand::bulkModulusAir(const double& I1_bar,
                                      ParameterDict& params) 
{
  double gamma = 1.4;
  double ps = I1_bar/3.0;
  double pr = 101325.0;
  double eps_v_a = 1.0/gamma*log(1.0 + ps/pr);
  double K_a = gamma*pr*exp(gamma*eps_v_a);
  return K_a;
}

/*!
 * -------------------------------------------------------------
 *  Bulk modulus model for water
 * -------------------------------------------------------------
 */
double
InternalVar_MasonSand::bulkModulusWater(const double& I1_bar,
                                        ParameterDict& params) 
{
  double K_w = params["Kw"];
  return K_w;
}

/*!
 * -------------------------------------------------------------
 *  Bulk modulus model for the drained sand
 * -------------------------------------------------------------
 */
double
InternalVar_MasonSand::bulkModulusDrainedSand(const double& I1_bar, 
                                              const double& ev_p_bar,
                                              ParameterDict& params) 
{
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl ,
                                              const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pKappaLabel,      matlset, Ghost::None);
  task->requires(Task::OldDW, pCapXLabel,       matlset, Ghost::None);
  task->requires(Task::OldDW, pPorosityLabel,   matlset, Ghost::None);
  task->requires(Task::OldDW, pSaturationLabel, matlset, Ghost::None);
  task->computes(pKappaLabel_preReloc,      matlset);
  task->computes(pCapXLabel_preReloc,       matlset);
  task->computes(pPorosityLabel_preReloc,   matlset);
  task->computes(pSaturationLabel_preReloc, matlset);
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

  from.push_back(pPorosityLabel);
  to.push_back(pPorosityLabel_preReloc);

  from.push_back(pSaturationLabel);
  to.push_back(pSaturationLabel_preReloc);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::allocateCMDataAddRequires(Task* task,
                                                 const MPMMaterial* matl ,
                                                 const PatchSet* ,
                                                 MPMLabel* )
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pKappaLabel_preReloc,      matlset, Ghost::None);
  task->requires(Task::NewDW, pCapXLabel_preReloc,       matlset, Ghost::None);
  task->requires(Task::NewDW, pPorosityLabel_preReloc,   matlset, Ghost::None);
  task->requires(Task::NewDW, pSaturationLabel_preReloc, matlset, Ghost::None);
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::allocateCMDataAdd(DataWarehouse* old_dw,
                                         ParticleSubset* addset,
                                         ParticleLabelVariableMap* newState,
                                         ParticleSubset* delset,
                                         DataWarehouse* new_dw )
{
  ParticleVariable<double> pKappa, pCapX, pPhi, pSw;
  constParticleVariable<double> o_kappa, o_capX, o_p3, o_phi, o_Sw;

  new_dw->allocateTemporary(pKappa,addset);
  new_dw->allocateTemporary(pCapX, addset);
  new_dw->allocateTemporary(pPhi,  addset);
  new_dw->allocateTemporary(pSw,   addset);

  new_dw->get(o_kappa,pKappaLabel_preReloc,      delset);
  new_dw->get(o_capX, pCapXLabel_preReloc,       delset);
  new_dw->get(o_phi,  pPorosityLabel_preReloc,   delset);
  new_dw->get(o_Sw,   pSaturationLabel_preReloc, delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for(o = delset->begin(); o != delset->end(); o++, n++) {
    pKappa[*n] = o_kappa[*o];
    pCapX[*n]  = o_capX[*o];
    pPhi[*n]    = o_phi[*o];
    pSw[*n]    = o_Sw[*o];
  }

  (*newState)[pKappaLabel]      = pKappa.clone();
  (*newState)[pCapXLabel]       = pCapX.clone();
  (*newState)[pPorosityLabel]   = pPhi.clone();
  (*newState)[pSaturationLabel] = pSw.clone();
}

/*!-----------------------------------------------------*/
void 
InternalVar_MasonSand::getInternalVariable(ParticleSubset* pset ,
                                           DataWarehouse* old_dw,
                                           constParticleLabelVariableMap& var)
{
  constParticleVariable<double> pKappa, pCapX, pPhi, pSw;
  old_dw->get(pKappa, pKappaLabel,      pset);
  old_dw->get(pCapX,  pCapXLabel,       pset);
  old_dw->get(pPhi,   pPorosityLabel,   pset);
  old_dw->get(pSw,    pSaturationLabel, pset);

  var[pKappaLabel]      = &pKappa;
  var[pCapXLabel]       = &pCapX;
  var[pPorosityLabel]   = &pPhi;
  var[pSaturationLabel] = &pSw;
}

void 
InternalVar_MasonSand::allocateAndPutInternalVariable(ParticleSubset* pset,
                                                      DataWarehouse* new_dw,
                                                      ParticleLabelVariableMap& var_new) 
{
  new_dw->allocateAndPut(*var_new[pKappaLabel],      pKappaLabel_preReloc,      pset);
  new_dw->allocateAndPut(*var_new[pCapXLabel],       pCapXLabel_preReloc,       pset);
  new_dw->allocateAndPut(*var_new[pPorosityLabel],   pPorosityLabel_preReloc,   pset);
  new_dw->allocateAndPut(*var_new[pSaturationLabel], pSaturationLabel_preReloc, pset);
}

void
InternalVar_MasonSand::allocateAndPutRigid(ParticleSubset* pset,
                                               DataWarehouse* new_dw,
                                               constParticleLabelVariableMap& var)
{
  ParticleVariable<double> pKappa_new, pCapX_new, pPhi_new, pSw_new;
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc,      pset);
  new_dw->allocateAndPut(pCapX_new,  pCapXLabel_preReloc,       pset);
  new_dw->allocateAndPut(pPhi_new,   pPorosityLabel_preReloc,   pset);
  new_dw->allocateAndPut(pSw_new,    pSaturationLabel_preReloc, pset);
  for(auto iter = pset->begin(); iter != pset->end(); iter++){
     pKappa_new[*iter] = 
       dynamic_cast<constParticleVariable<double>& >(*var[pKappaLabel])[*iter];
     pCapX_new[*iter]  = 
       dynamic_cast<constParticleVariable<double>& >(*var[pCapXLabel])[*iter];
     pPhi_new[*iter]    = 
       dynamic_cast<constParticleVariable<double>& >(*var[pPorosityLabel])[*iter];
     pSw_new[*iter]    = 
       dynamic_cast<constParticleVariable<double>& >(*var[pSaturationLabel])[*iter];
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

