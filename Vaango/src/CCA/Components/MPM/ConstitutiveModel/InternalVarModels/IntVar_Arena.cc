/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/IntVar_Arena.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include<CCA/Components/MPM/Core/MPMLabel.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;

/*!-----------------------------------------------------*/
IntVar_Arena::IntVar_Arena(ProblemSpecP& ps, ElasticModuliModel* elastic)
{
  d_elastic = elastic;
  d_shear   = nullptr;

  ps->require("p0", d_crushParam.p0); // Crush Curve Parameter
  ps->require("p1", d_crushParam.p1); // Crush Curve Parameter
  ps->require("p2", d_crushParam.p2); // Crush Curve Parameter
  ps->require("p3", d_crushParam.p3); // Crush Curve Parameter

  ps->getWithDefault(
    "use_disaggregation_algorithm", d_use_disaggregation_algorithm, false);

  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}

/*!-----------------------------------------------------*/
IntVar_Arena::IntVar_Arena(const IntVar_Arena* cm)
{
  d_elastic = cm->d_elastic;
  d_shear   = cm->d_shear;

  d_crushParam                   = cm->d_crushParam;
  d_use_disaggregation_algorithm = cm->d_use_disaggregation_algorithm;

  // Initialize internal variable labels for evolution
  initializeLocalMPMLabels();
}

void
IntVar_Arena::initializeLocalMPMLabels()
{
  pKappaLabel = Uintah::VarLabel::create(
    "p.kappa", Uintah::ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc = Uintah::VarLabel::create(
    "p.kappa+", Uintah::ParticleVariable<double>::getTypeDescription());

  pCapXLabel = Uintah::VarLabel::create(
    "p.capX", Uintah::ParticleVariable<double>::getTypeDescription());
  pCapXLabel_preReloc = Uintah::VarLabel::create(
    "p.capX+", Uintah::ParticleVariable<double>::getTypeDescription());

  pPlasticStrainLabel = Uintah::VarLabel::create(
    "p.plasticStrain",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());
  pPlasticStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.plasticStrain+",
    Uintah::ParticleVariable<Uintah::Matrix3>::getTypeDescription());

  pPlasticVolStrainLabel = Uintah::VarLabel::create(
    "p.plasticVolStrain",
    Uintah::ParticleVariable<double>::getTypeDescription());
  pPlasticVolStrainLabel_preReloc = Uintah::VarLabel::create(
    "p.plasticVolStrain+",
    Uintah::ParticleVariable<double>::getTypeDescription());

  pP3Label = Uintah::VarLabel::create(
    "p.p3", Uintah::ParticleVariable<double>::getTypeDescription());
  pP3Label_preReloc = Uintah::VarLabel::create(
    "p.p3+", Uintah::ParticleVariable<double>::getTypeDescription());
}

/*!-----------------------------------------------------*/
IntVar_Arena::~IntVar_Arena()
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
void
IntVar_Arena::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "arena");

  int_var_ps->appendElement("p0", d_crushParam.p0);
  int_var_ps->appendElement("p1", d_crushParam.p1);
  int_var_ps->appendElement("p2", d_crushParam.p2);
  int_var_ps->appendElement("p3", d_crushParam.p3);

  int_var_ps->appendElement("use_disaggregation_algorithm",
                            d_use_disaggregation_algorithm);
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::addInitialComputesAndRequires(Task* task,
                                            const MPMMaterial* matl,
                                            const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pKappaLabel, matlset);
  task->computes(pCapXLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(pPlasticVolStrainLabel, matlset);
  task->computes(pP3Label, matlset);
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::initializeInternalVariable(const Patch* patch,
                                         const MPMMaterial* matl,
                                         ParticleSubset* pset,
                                         DataWarehouse* new_dw,
                                         MPMLabel* lb,
                                         ParameterDict& params)
{
  Uintah::constParticleVariable<double> pMass, pVolume;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pMass, lb->pMassLabel, pset);

  Uintah::ParticleVariable<double> pKappa, pCapX;
  Uintah::ParticleVariable<double> pPlasticVolStrain, pP3;
  Uintah::ParticleVariable<Matrix3> pPlasticStrain;
  new_dw->allocateAndPut(pKappa, pKappaLabel, pset);
  new_dw->allocateAndPut(pCapX, pCapXLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticVolStrain, pPlasticVolStrainLabel, pset);
  new_dw->allocateAndPut(pP3, pP3Label, pset);

  double PEAKI1 = 0.0;
  double CR     = 0.0;
  double phi0   = 0.0;
  double Sw0    = 0.0;
  double pf0    = 0.0;
  try {
    PEAKI1 = params.at("PEAKI1");
    CR     = params.at("CR");
    phi0   = params.at("phi0");
    Sw0    = params.at("Sw0");
    pf0    = params.at("Pf0");
  } catch (const std::out_of_range& e) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameters PEAKI1, CR, phi0, Sw0"
        << std::endl;
    err << "\t Available parameters are:" << std::endl;
    for (auto param : params) {
      err << "\t \t" << param.first << " " << param.second << std::endl;
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
  }

  for (auto iter : *pset) {
    if (d_use_disaggregation_algorithm) {
      pP3[iter] =
        log(pVolume[iter] * (matl->getInitialDensity()) / pMass[iter]);
    } else {
      pP3[iter] = d_crushParam.p3;
    }
    double I1_eff_bar = 0.0;
    double ep_v_bar   = 0.0;
    double pw_bar     = pf0;
    if (Sw0 > 0.0) {
      double phi   = computePorosity(ep_v_bar, pP3[iter]);
      double X_bar = computePartSatHydrostaticStrength(
        I1_eff_bar, pw_bar, ep_v_bar, phi0, Sw0, phi);
      pCapX[iter] = -X_bar;
    } else {
      double X_bar = computeDrainedHydrostaticStrength(ep_v_bar, phi0);
      pCapX[iter]  = -X_bar;
    }
    // std::cout << "pCapX = " << pCapX[*iter] << std::endl;
    pKappa[iter] = PEAKI1 - CR * (PEAKI1 - pCapX[iter]); // Branch Point
    pPlasticStrain[iter].set(0.0);
    pPlasticVolStrain[iter] = 0.0;
  }
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::addComputesAndRequires(Task* task,
                                     const MPMMaterial* matl,
                                     const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->needs(Task::OldDW, pKappaLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pCapXLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pPlasticStrainLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pPlasticVolStrainLabel, matlset, Ghost::None);
  task->needs(Task::OldDW, pP3Label, matlset, Ghost::None);
  task->computes(pKappaLabel_preReloc, matlset);
  task->computes(pCapXLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(pPlasticVolStrainLabel_preReloc, matlset);
  task->computes(pP3Label_preReloc, matlset);
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::addParticleState(std::vector<const VarLabel*>& from,
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

/* Get one (possibly composite) internal variable */
template <>
void
IntVar_Arena::getInternalVariable(ParticleSubset* pset,
                                  DataWarehouse* old_dw,
                                  constParticleVariable<ArenaIntVar>& intvar)
{
  std::ostringstream err;
  err << "**ERROR** getInternalVariable with ArenaIntVar not implemented yet.";
  throw InternalError(err.str(), __FILE__, __LINE__);
}

/* Get multiple local <int/double/Vector/Matrix3> internal variables */
template <>
std::vector<Uintah::constParticleVariable<double>>
IntVar_Arena::getInternalVariables(ParticleSubset* pset, DataWarehouse* old_dw)
{
  constParticleVariable<double> pKappa, pCapX, pPlasticVolStrain, pP3;
  old_dw->get(pKappa, pKappaLabel, pset);
  old_dw->get(pCapX, pCapXLabel, pset);
  old_dw->get(pPlasticVolStrain, pPlasticVolStrainLabel, pset);
  old_dw->get(pP3, pP3Label, pset);

  std::vector<constParticleVariable<double>> pIntVars;
  pIntVars.emplace_back(pKappa);
  pIntVars.emplace_back(pCapX);
  pIntVars.emplace_back(pPlasticVolStrain);
  pIntVars.emplace_back(pP3);

  return pIntVars;
}

template <>
std::vector<Uintah::constParticleVariable<Matrix3>>
IntVar_Arena::getInternalVariables(ParticleSubset* pset, DataWarehouse* old_dw)
{
  constParticleVariable<Uintah::Matrix3> pPlasticStrain;
  old_dw->get(pPlasticStrain, pPlasticStrainLabel, pset);

  std::vector<Uintah::constParticleVariable<Uintah::Matrix3>> pIntVars;
  pIntVars.emplace_back(pPlasticStrain);

  return pIntVars;
}

/* Allocate multiple local <int/double/Vector/Matrix3> internal variables */
template <>
void
IntVar_Arena::allocateAndPutInternalVariable(
  ParticleSubset* pset,
  DataWarehouse* new_dw,
  std::vector<ParticleVariable<double>>& pVars)
{
  new_dw->allocateAndPut(pVars[0], pKappaLabel_preReloc, pset);
  new_dw->allocateAndPut(pVars[1], pCapXLabel_preReloc, pset);
  new_dw->allocateAndPut(pVars[2], pPlasticVolStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(pVars[3], pP3Label_preReloc, pset);
}

template <>
void
IntVar_Arena::allocateAndPutInternalVariable(
  ParticleSubset* pset,
  DataWarehouse* new_dw,
  std::vector<::ParticleVariable<Matrix3>>& pVars)
{
  new_dw->allocateAndPut(pVars[0], pPlasticStrainLabel_preReloc, pset);
}

template <>
void
IntVar_Arena::evolveInternalVariable(
  Uintah::particleIndex pidx,
  const ModelStateBase* state,
  Uintah::constParticleVariable<ArenaIntVar>& var_old,
  Uintah::ParticleVariable<ArenaIntVar>& var)
{
}

//--------------------------------------------------------------------------------------
// Compute hydrostatic strength
//--------------------------------------------------------------------------------------
double
IntVar_Arena::computeInternalVariable(const std::string& label,
                                      const ModelStateBase* state_old,
                                      const ModelStateBase* state_cur) const
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_cur);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the stress, backstress, plastic strain, and initial porosity
  double I1_eff = state->I1_eff;
  // double pbar_w = state->pbar_w;
  double ep_v = state->ep_v;
  double phi0 = computePorosity(0.0, state->p3);
  // std::cout << "IntVar_Arena:: I1_eff = " << I1_eff << " pbar_w = " << pbar_w
  //          << " ep_v = " << ep_v
  //          << " phi0 = " << phi0 << std::endl;

  if (ep_v > 0.0) { // tension
    double X_bar = computeDrainedHydrostaticStrength(0.0, phi0);
    return -X_bar;
  }

  // Convert to bar quantities
  double ep_v_bar   = -ep_v;
  double I1_eff_bar = -I1_eff; // This means that the bulk modulus of the
                               // partially saturated material is computed
                               // at zero pore pressure
  double pw_bar = state->pbar_w;

  // Compute the hydrostatic compressive strength
  double X_bar = 0.0;
  if (state->saturation > 0.0) {
    X_bar = computePartSatHydrostaticStrength(
      I1_eff_bar, pw_bar, ep_v_bar, state->porosity, state->saturation, phi0);
  } else {
    X_bar = computeDrainedHydrostaticStrength(ep_v_bar, phi0);
  }

  // Return the new capX
  double capX_new = -X_bar;
  return capX_new;
}

/**
 *  Compute drained hydrostatic strength
 */
double
IntVar_Arena::computeDrainedHydrostaticStrength(const double& ep_v_bar,
                                                const double& phi0) const
{
  double p0 = d_crushParam.p0;
  double p1 = d_crushParam.p1;
  double p2 = d_crushParam.p2;
  double p3 = computeP3(phi0);

  double X_bar_drained = std::max(p0, 1000.0);
  // std::cout << "\t\t ep_v_bar = " << ep_v_bar << std::endl;
  if (ep_v_bar > 0.0) {
    // double phi0 = 1.0 - std::exp(-p3);
    double phi    = 1.0 - std::exp(-p3 + ep_v_bar);
    double term1  = (phi0 / phi - 1.0);
    double xi_bar = p1 * std::pow(term1, 1.0 / p2);
    X_bar_drained += xi_bar;
    // std::cout << "\t\t phi = " << phi << " xi_bar = " << xi_bar
    //          << " X_bar_drained = " << X_bar_drained << std::endl;
  }

  return X_bar_drained;
}

/**
 * Compute elastic volume strain at yield
 */
double
IntVar_Arena::computeElasticVolStrainAtYield(const double& ep_v_bar,
                                             const double& phi0) const
{
  // Compute X(ep_v) using crush curve model for dry sand
  double X_bar_ep_v = computeDrainedHydrostaticStrength(ep_v_bar, phi0);

  // Compute K(ep_v, I1) using bulk modulus model for dry sand
  double I1_eff_bar  = 0.5 * X_bar_ep_v;
  double bulkModulus = 0.0, shearModulus = 0.0;
  d_elastic->computeDrainedModuli(
    I1_eff_bar, ep_v_bar, bulkModulus, shearModulus);

  // Compute elastic vol strain at yield
  double ev_e_yield = X_bar_ep_v / (3.0 * bulkModulus);

  return ev_e_yield;
}

/**
 * Compute partially saturated hydrostatic strength
 */
double
IntVar_Arena::computePartSatHydrostaticStrength(const double& I1_eff_bar,
                                                const double& pw_bar,
                                                const double& ep_v_bar,
                                                const double& phi,
                                                const double& Sw,
                                                const double& phi0) const
{
  //------------Plastic strain exceeds allowable limit--------------------------
  // The plastic strain for this iteration has exceed the allowable
  // value.  X is not defined in this region, so we set it to a large
  // negative number.
  //
  // The code should never have ep_v_bar > p3, but will have ep_v_bar = p3 if
  // the
  // porosity approaches zero (within the specified tolerance).  By setting
  // X=1e12, the material will respond as though there is no porosity.
  double p3 = computeP3(phi0);
  if (ep_v_bar > p3) {
    double X_bar_limit = 1.0e12;
    return X_bar_limit;
  }

  // ------------------Plastic strain is within allowable
  // domain------------------------
  // Compute elastic volumetric strain at yield
  double ev_e_yield = computeElasticVolStrainAtYield(ep_v_bar, phi0);

  // Compute partially saturated bulk modulus
  double K_part_sat = 0.0, G_part_sat = 0.0;
  d_elastic->computePartialSaturatedModuli(
    I1_eff_bar, pw_bar, ep_v_bar, phi, Sw, K_part_sat, G_part_sat);

  // Compute hydrostatic strength
  double X_bar_part_sat = 3.0 * K_part_sat * ev_e_yield;

  return X_bar_part_sat;
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::allocateCMDataAddRequires(Task* task,
                                        const MPMMaterial* matl,
                                        const PatchSet*,
                                        MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->needs(Task::NewDW, pKappaLabel_preReloc, matlset, Ghost::None);
  task->needs(Task::NewDW, pCapXLabel_preReloc, matlset, Ghost::None);
  task->needs(
    Task::NewDW, pPlasticStrainLabel_preReloc, matlset, Ghost::None);
  task->needs(
    Task::NewDW, pPlasticVolStrainLabel_preReloc, matlset, Ghost::None);
  task->needs(Task::NewDW, pP3Label_preReloc, matlset, Ghost::None);
}

/*!-----------------------------------------------------*/
void
IntVar_Arena::allocateCMDataAdd(DataWarehouse* old_dw,
                                ParticleSubset* addset,
                                ParticleLabelVariableMap* newState,
                                ParticleSubset* delset,
                                DataWarehouse* new_dw)
{
  ParticleVariable<double> pKappa, pCapX, pEpv, pP3;
  ParticleVariable<Matrix3> pEp;
  constParticleVariable<double> o_kappa, o_capX, o_Epv, o_P3;
  constParticleVariable<Matrix3> o_Ep;

  new_dw->allocateTemporary(pKappa, addset);
  new_dw->allocateTemporary(pCapX, addset);
  new_dw->allocateTemporary(pEp, addset);
  new_dw->allocateTemporary(pEpv, addset);
  new_dw->allocateTemporary(pP3, addset);

  new_dw->get(o_kappa, pKappaLabel_preReloc, delset);
  new_dw->get(o_capX, pCapXLabel_preReloc, delset);
  new_dw->get(o_Ep, pPlasticStrainLabel_preReloc, delset);
  new_dw->get(o_Epv, pPlasticVolStrainLabel_preReloc, delset);
  new_dw->get(o_P3, pP3Label_preReloc, delset);

  auto o = addset->begin();
  auto n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
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
IntVar_Arena::allocateAndPutRigid(ParticleSubset* pset,
                                  DataWarehouse* new_dw,
                                  constParticleLabelVariableMap& var)
{
  ParticleVariable<double> pKappa_new, pCapX_new, pEpv_new, pP3_new;
  ParticleVariable<Matrix3> pEp_new;
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc, pset);
  new_dw->allocateAndPut(pCapX_new, pCapXLabel_preReloc, pset);
  new_dw->allocateAndPut(pEp_new, pPlasticStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(pEpv_new, pPlasticVolStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(pP3_new, pP3Label_preReloc, pset);
  for (int& iter : *pset) {
    pKappa_new[iter] =
      dynamic_cast<constParticleVariable<double>&>(*var[pKappaLabel])[iter];
    pCapX_new[iter] =
      dynamic_cast<constParticleVariable<double>&>(*var[pCapXLabel])[iter];
    pEp_new[iter] = dynamic_cast<constParticleVariable<Matrix3>&>(
      *var[pPlasticStrainLabel])[iter];
    pEpv_new[iter] = dynamic_cast<constParticleVariable<double>&>(
      *var[pPlasticVolStrainLabel])[iter];
    pP3_new[iter] =
      dynamic_cast<constParticleVariable<double>&>(*var[pP3Label])[iter];
  }
}
