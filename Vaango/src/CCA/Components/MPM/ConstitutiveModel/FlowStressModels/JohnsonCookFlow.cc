/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/FlowStressModels/JohnsonCookFlow.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <cmath>

using namespace Uintah;
using Vaango::ModelStateBase;

JohnsonCookFlow::JohnsonCookFlow(ProblemSpecP& ps)
{
  ps->require("A", d_CM.A);
  ps->require("B", d_CM.B);
  ps->require("C", d_CM.C);
  ps->require("n", d_CM.n);
  ps->require("m", d_CM.m);
  d_CM.epdot_0 = 1.0;
  ps->get("epdot_0", d_CM.epdot_0);
  d_CM.TRoom = 298;
  ps->get("T_r", d_CM.TRoom);
  d_CM.TMelt = 1793;
  ps->get("T_m", d_CM.TMelt);
}

JohnsonCookFlow::JohnsonCookFlow(const JohnsonCookFlow* cm)
{
  d_CM.A = cm->d_CM.A;
  d_CM.B = cm->d_CM.B;
  d_CM.C = cm->d_CM.C;
  d_CM.n = cm->d_CM.n;
  d_CM.m = cm->d_CM.m;
  d_CM.epdot_0 = cm->d_CM.epdot_0;
  d_CM.TRoom = cm->d_CM.TRoom;
  d_CM.TMelt = cm->d_CM.TMelt;
}

JohnsonCookFlow::~JohnsonCookFlow() = default;

void
JohnsonCookFlow::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP flow_ps = ps->appendChild("flow_model");
  flow_ps->setAttribute("type", "johnson_cook");

  flow_ps->appendElement("A", d_CM.A);
  flow_ps->appendElement("B", d_CM.B);
  flow_ps->appendElement("C", d_CM.C);
  flow_ps->appendElement("n", d_CM.n);
  flow_ps->appendElement("m", d_CM.m);
  flow_ps->appendElement("epdot_0", d_CM.epdot_0);
  flow_ps->appendElement("T_r", d_CM.TRoom);
  flow_ps->appendElement("T_m", d_CM.TMelt);
}

void
JohnsonCookFlow::addInitialComputesAndRequires(Task*, const MPMMaterial*,
                                               const PatchSet*)
{
}

void
JohnsonCookFlow::addComputesAndRequires(Task*, const MPMMaterial*,
                                        const PatchSet*)
{
}

void
JohnsonCookFlow::addComputesAndRequires([[maybe_unused]] Task* task, [[maybe_unused]] const MPMMaterial* matl,
                                        const PatchSet*, bool /*recurse*/,
                                        bool /*SchedParent*/)
{
}

void
JohnsonCookFlow::addParticleState(std::vector<const VarLabel*>&,
                                  std::vector<const VarLabel*>&)
{
}

void
JohnsonCookFlow::allocateCMDataAddRequires(Task*, const MPMMaterial*,
                                           const PatchSet*, MPMLabel*)
{
}

void
JohnsonCookFlow::allocateCMDataAdd(DataWarehouse*, ParticleSubset*,
                                   ParticleLabelVariableMap*,
                                   ParticleSubset*, DataWarehouse*)
{
}

void
JohnsonCookFlow::initializeInternalVars(ParticleSubset*, DataWarehouse*)
{
}

void
JohnsonCookFlow::getInternalVars(ParticleSubset*, DataWarehouse*)
{
}

void
JohnsonCookFlow::allocateAndPutInternalVars(ParticleSubset*, DataWarehouse*)
{
}

void
JohnsonCookFlow::allocateAndPutRigid(ParticleSubset*, DataWarehouse*)
{
}

void
JohnsonCookFlow::updateElastic(const particleIndex)
{
}

void
JohnsonCookFlow::updatePlastic(const particleIndex, const double&)
{
}

double
JohnsonCookFlow::computeFlowStress(const ModelStateBase* state, const double&,
                                   const double&, [[maybe_unused]] const MPMMaterial* matl,
                                   const particleIndex idx)
{
  //double epdot = state->eqPlasticStrainRate/d_CM.epdot_0;
  double epdot = state->eqStrainRate / d_CM.epdot_0;
  double ep = state->eqPlasticStrain;
  double T = state->temperature;
  // double Tr = matl->getRoomTemperature();
  // double Tm = state->meltingTemp;
  double Tr = d_CM.TRoom;
  double Tm = d_CM.TMelt;

  double strainPart = d_CM.A + d_CM.B * pow(ep, d_CM.n);
  double strainRatePart = 1.0;
  if (epdot < 1.0) {
    strainRatePart = pow((1.0 + epdot), d_CM.C);
  } else {
    strainRatePart = 1.0 + d_CM.C * log(epdot);
  }
  double m = d_CM.m;
  double Tstar = (T > Tm) ? 1.0 : ((T - Tr) / (Tm - Tr));
  double tempPart = (Tstar < 0.0) ? (1.0 - Tstar) : (1.0 - pow(Tstar, m));
  double sigy = strainPart * strainRatePart * tempPart;
  if (std::isnan(sigy)) {
    std::cout << "**ERROR** JohnsonCook: sig_y == nan "
         << " Particle = " << idx << " epdot = " << epdot << " ep = " << ep
         << " T = " << T << " strainPart = " << strainPart
         << " strainRatePart = " << strainRatePart << " Tstar = " << Tstar
         << " tempPart = " << tempPart << "\n";
  }
  return sigy;
}

double
JohnsonCookFlow::computeEpdot(const ModelStateBase* state, const double&,
                              const double&, [[maybe_unused]] const MPMMaterial* matl,
                              const particleIndex)
{
  // All quantities should be at the beginning of the
  // time step
  double tau = state->yieldStress;
  double ep = state->eqPlasticStrain;
  double T = state->temperature;
  // double Tr = matl->getRoomTemperature();
  // double Tm = state->meltingTemp;
  double Tr = d_CM.TRoom;
  double Tm = d_CM.TMelt;

  double strainPart = d_CM.A + d_CM.B * pow(ep, d_CM.n);
  d_CM.TRoom = Tr;
  d_CM.TMelt = Tm;
  double m = d_CM.m;
  double Tstar = (T > Tm) ? 1.0 : ((T - Tr) / (Tm - Tr));
  double tempPart = (Tstar < 0.0) ? (1.0 - Tstar) : (1.0 - pow(Tstar, m));

  double fac1 = tau / (strainPart * tempPart);
  double fac2 = (1.0 / d_CM.C) * (fac1 - 1.0);
  double epdot = exp(fac2) * d_CM.epdot_0;
  if (std::isnan(epdot)) {
    std::cout << "**ERROR** JohnsonCook: epdot == nan " << "\n";
  }
  return epdot;
}

void
JohnsonCookFlow::computeTangentModulus([[maybe_unused]] const Matrix3& stress,
                                       const ModelStateBase*, const double&,
                                       const MPMMaterial*, const particleIndex,
                                       TangentModulusTensor&,
                                       TangentModulusTensor&)
{
  throw InternalError("Empty Function: JohnsonCookFlow::computeTangentModulus",
                      __FILE__, __LINE__);
}

void
JohnsonCookFlow::evalDerivativeWRTScalarVars(const ModelStateBase* state,
                                             const particleIndex idx,
                                             Vector& derivs) const
{
  derivs[0] = evalDerivativeWRTStrainRate(state, idx);
  derivs[1] = evalDerivativeWRTTemperature(state, idx);
  derivs[2] = evalDerivativeWRTPlasticStrain(state, idx);
}

double
JohnsonCookFlow::evalDerivativeWRTPlasticStrain(const ModelStateBase* state,
                                                const particleIndex) const
{
  // Get the state data
  //double epdot = state->eqPlasticStrainRate/d_CM.epdot_0;
  double epdot = state->eqStrainRate / d_CM.epdot_0;
  double ep = state->eqPlasticStrain;
  double T = state->temperature;
  // double Tm = state->meltingTemp;
  double Tr = d_CM.TRoom;
  double Tm = d_CM.TMelt;

  // Calculate strain rate part
  double strainRatePart =
    (epdot < 1.0) ? (pow((1.0 + epdot), d_CM.C)) : (1.0 + d_CM.C * log(epdot));

  // Calculate temperature part
  double m = d_CM.m;
  double Tstar = (T > Tm) ? 1.0 : (T - Tr) / (Tm - Tr);
  double tempPart = (Tstar < 0.0) ? (1.0 - Tstar) : (1.0 - pow(Tstar, m));

  double D = strainRatePart * tempPart;

  double deriv = (ep > 0.0) ? (d_CM.B * d_CM.n * D * pow(ep, d_CM.n - 1)) : 0.0;
  if (std::isnan(deriv)) {
    std::cout << "**ERROR** JohnsonCook: dsig/dep == nan " << "\n";
  }
  return deriv;
}

///////////////////////////////////////////////////////////////////////////
/*  Compute the shear modulus. */
///////////////////////////////////////////////////////////////////////////
double
JohnsonCookFlow::computeShearModulus(const ModelStateBase* state)
{
  return state->shearModulus;
}

///////////////////////////////////////////////////////////////////////////
/* Compute the melting temperature */
///////////////////////////////////////////////////////////////////////////
double
JohnsonCookFlow::computeMeltingTemp(const ModelStateBase* state)
{
  return state->meltingTemp;
}

double
JohnsonCookFlow::evalDerivativeWRTTemperature(const ModelStateBase* state,
                                              const particleIndex) const
{
  // Get the state data
  //double epdot = state->eqPlasticStrainRate/d_CM.epdot_0;
  double epdot = state->eqStrainRate / d_CM.epdot_0;
  double ep = state->eqPlasticStrain;
  double T = state->temperature;
  // double Tm = state->meltingTemp;
  double Tr = d_CM.TRoom;
  double Tm = d_CM.TMelt;

  // Calculate strain part
  double strainPart = d_CM.A + d_CM.B * pow(ep, d_CM.n);

  // Calculate strain rate part
  double strainRatePart = 1.0;
  if (epdot < 1.0)
    strainRatePart = pow((1.0 + epdot), d_CM.C);
  else
    strainRatePart = 1.0 + d_CM.C * log(epdot);

  // Calculate temperature part
  double m = d_CM.m;
  double Tstar = (T > Tm) ? 1.0 : (T - Tr) / (Tm - Tr);

  double F = strainPart * strainRatePart;
  double deriv =
    (T < Tr) ? -F / (Tm - Tr) : -m * F * pow(Tstar, m - 1) / (Tm - Tr);
  if (std::isnan(deriv)) {
    std::cout << "**ERROR** JohnsonCook: dsig/dT == nan " << "\n";
  }
  return deriv;
}

double
JohnsonCookFlow::evalDerivativeWRTStrainRate(const ModelStateBase* state,
                                             const particleIndex) const
{
  // Get the state data
  //double epdot = state->eqPlasticStrainRate/d_CM.epdot_0;
  double epdot = state->eqStrainRate / d_CM.epdot_0;
  double ep = state->eqPlasticStrain;
  double T = state->temperature;
  // double Tm = state->meltingTemp;
  double Tr = d_CM.TRoom;
  double Tm = d_CM.TMelt;

  // Calculate strain part
  double strainPart = d_CM.A + d_CM.B * pow(ep, d_CM.n);

  // Calculate temperature part
  double m = d_CM.m;
  double Tstar = (T > Tm) ? 1.0 : (T - Tr) / (Tm - Tr);
  double tempPart = (Tstar < 0.0) ? (1.0 - Tstar) : (1.0 - pow(Tstar, m));

  double E = strainPart * tempPart;

  double deriv = 0.0;
  if (epdot < 1.0)
    deriv = E * d_CM.C * pow((1.0 + epdot), (d_CM.C - 1.0));
  else
    deriv = E * d_CM.C / epdot;
  if (std::isnan(deriv)) {
    std::cout << "**ERROR** JohnsonCook: dsig/depdot == nan " << "\n";
  }
  return deriv;
}
