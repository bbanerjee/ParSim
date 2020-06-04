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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModels/HypoElasticFracture.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/Constants.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <fstream>
#include <iostream>

using namespace Uintah;

HypoElasticFracture::HypoElasticFracture(ProblemSpecP& ps, MPMFlags* Mflag)
  : HypoElastic(ps, Mflag)
{
  // Read in fracture criterion and the toughness curve
  ProblemSpecP curve_ps = ps->findBlock("fracture_toughness_curve");
  if (curve_ps != 0) {
    d_crackPropagationCriterion = "max_hoop_stress";
    curve_ps->get("crack_propagation_criterion", d_crackPropagationCriterion);

    if (d_crackPropagationCriterion != "max_hoop_stress" &&
        d_crackPropagationCriterion != "max_principal_stress" &&
        d_crackPropagationCriterion != "max_energy_release_rate" &&
        d_crackPropagationCriterion != "strain_energy_density" &&
        d_crackPropagationCriterion != "empirical_criterion") {
      std::ostringstream err;
      err << "**Error** Undefinded crack propagation criterion: "
          << d_crackPropagationCriterion
          << " for hypoelastic materials. " << "\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }

    d_fracParam.p = 2.0; // Default elliptical fracture locus
    d_fracParam.q = 2.0; 
    d_fracParam.r = -1.;
    if (d_crackPropagationCriterion == "empirical_criterion") {
      // Get parameters p & q in the fracture locus equation
      // (KI/Ic)^p+(KII/KIIc)^q=1 and KIIc=r*KIc
      curve_ps->get("p", d_fracParam.p);
      curve_ps->get("q", d_fracParam.q);
      curve_ps->get("r", d_fracParam.r);
    }

    for (ProblemSpecP child_ps = curve_ps->findBlock("point"); child_ps != 0;
         child_ps = child_ps->findNextBlock("point")) {
      double Vc, KIc, KIIc;
      child_ps->require("Vc", Vc);
      child_ps->require("KIc", KIc);
      if (d_fracParam.r < 0.) { // Input KIIc manually
        child_ps->get("KIIc", KIIc);
      } else { // The ratio of KIIc to KIc is a constant (r)
        KIIc = d_fracParam.r * KIc;
      }
      d_fracParam.Kc.push_back(Toughness(Vc, KIc, KIIc));
    }
  }
}

HypoElasticFracture::HypoElasticFracture(const HypoElasticFracture* cm)
  : HypoElastic(cm)
{
  d_crackPropagationCriterion = cm->d_crackPropagationCriterion;
  d_fracParam.p = cm->d_fracParam.p;
  d_fracParam.q = cm->d_fracParam.q;
  d_fracParam.r = cm->d_fracParam.r;
  d_fracParam.Kc = cm->d_fracParam.Kc;
}

void
HypoElasticFracture::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "hypo_elastic_fracture");
  }

  ProblemSpecP curve_ps = cm_ps->appendChild("fracture_toughness_curve");
  curve_ps->appendElement("crack_propagation_criterion", d_crackPropagationCriterion);
  curve_ps->appendElement("p", d_fracParam.p);
  curve_ps->appendElement("q", d_fracParam.q);
  curve_ps->appendElement("r", d_fracParam.r);
  for (auto Kc : d_fracParam.Kc) {
    ProblemSpecP Kc_ps = cm_ps->appendChild("point");
    Kc_ps->appendElement("Vc", Kc.Vc);
    Kc_ps->appendElement("KIc", Kc.KIc);
    Kc_ps->appendElement("KIIc", Kc.KIIc);
  }
}

HypoElasticFracture*
HypoElasticFracture::clone()
{
  return scinew HypoElasticFracture(*this);
}

void
HypoElasticFracture::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                              DataWarehouse* new_dw)
{
  HypoElastic::initializeCMData(patch, matl, new_dw);

  // for J-Integral
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);
  ParticleVariable<Matrix3> pDispGrads;
  ParticleVariable<double> pStrainEnergyDensity;
  new_dw->allocateAndPut(pDispGrads, lb->pDispGradsLabel, pset);
  new_dw->allocateAndPut(pStrainEnergyDensity, lb->pStrainEnergyDensityLabel,
                         pset);

  for (auto particle : *pset) {
    pDispGrads[particle] = Vaango::Util::Zero;
    pStrainEnergyDensity[particle] = 0.0;
  }
}

void
HypoElasticFracture::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  from.push_back(lb->pDispGradsLabel);
  from.push_back(lb->pStrainEnergyDensityLabel);
  to.push_back(lb->pDispGradsLabel_preReloc);
  to.push_back(lb->pStrainEnergyDensityLabel_preReloc);
}

void
HypoElasticFracture::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                    const PatchSet* patches) const
{
  HypoElastic::addComputesAndRequires(task, matl, patches);

  Ghost::GhostType gnone = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, lb->pDispGradsLabel, matlset, gnone);
  task->requires(Task::OldDW, lb->pStrainEnergyDensityLabel, matlset, gnone);

  task->computes(lb->pDispGradsLabel_preReloc, matlset);
  task->computes(lb->pVelGradsLabel, matlset);
  task->computes(lb->pStrainEnergyDensityLabel_preReloc, matlset);
}

void
HypoElasticFracture::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  double rho_orig = matl->getInitialDensity();
  int matID = matl->getDWIndex();

  double K = d_modelParam.K;
  double mu = d_modelParam.G;
  double lambda = d_modelParam.K - 2.0 / 3.0 * d_modelParam.G;
  double M_dil = lambda + 2.0 * mu;

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    Vector dx = patch->dCell();

    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<double> pVolume_new, pStrainEnergyDensity;
    constParticleVariable<Vector> pVelocity;
    constParticleVariable<Matrix3> pDefRate_mid, pVelGrad_mid, pDispGrads,
                                   pDefGrad_new, pStress_old;

    old_dw->get(pStrainEnergyDensity, lb->pStrainEnergyDensityLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);
    old_dw->get(pDispGrads, lb->pDispGradsLabel, pset);
    old_dw->get(pStress_old, lb->pStressLabel, pset);

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pVelGrad_mid, lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q, pStrainEnergyDensity_new;
    ParticleVariable<Matrix3> pDispGrads_new, pVelGrads, pStress_new;

    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStrainEnergyDensity_new,
                           lb->pStrainEnergyDensityLabel_preReloc, pset);
    new_dw->allocateAndPut(pVelGrads, lb->pVelGradsLabel, pset);
    new_dw->allocateAndPut(pDispGrads_new, lb->pDispGradsLabel_preReloc,
                           pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
    for (int idx : *pset) {

      pdTdt[idx] = 0.0;

      pVelGrads[idx] = pVelGrad_mid[idx];
      pDispGrads_new[idx] = pDispGrads[idx] + pVelGrad_mid[idx] * delT;

      Matrix3 D = pDefRate_mid[idx];
      double trD = D.Trace();
      pStress_new[idx] = pStress_old[idx] +
        (Vaango::Util::Identity * (lambda * trD) + D * (2.0 * mu)) * delT;

      // Compute the strain energy for all the particles
      Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * .5;
      double rateOfWork = computeRateOfWork(avgStress, pDefRate_mid[idx]);
      double energy = (rateOfWork * pVolume_new[idx] * delT);
      strainEnergy += energy;

      pStrainEnergyDensity_new[idx] =
          pStrainEnergyDensity[idx] + rateOfWork * delT;

      // Compute wave speed at each particle, store the maximum
      double J = pDefGrad_new[idx].Determinant();
      double rho_cur = rho_orig / J;
      double c_dil = std::sqrt(M_dil / rho_cur);
      Vector velMax = pVelocity[idx].cwiseAbs() + c_dil;
      waveSpeed     = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificialViscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = std::sqrt(K / rho_cur);
        p_q[idx] = artificialBulkViscosity(trD, c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

void
HypoElasticFracture::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                                    const bool, const bool) const
{
}

void
HypoElasticFracture::allocateCMDataAddRequires(Task* task, const MPMMaterial* matl,
                                       const PatchSet* patches,
                                       MPMLabel* lb) const
{
  HypoElastic::allocateCMDataAddRequires(task, matl, patches, lb);

  Ghost::GhostType gnone = Ghost::None;
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, lb->pDispGradsLabel_preReloc, matlset, gnone);
  task->requires(Task::NewDW, lb->pStrainEnergyDensityLabel_preReloc, matlset,
                 gnone);
}

void
HypoElasticFracture::allocateCMDataAdd(DataWarehouse* new_dw, ParticleSubset* addset,
                               ParticleLabelVariableMap* newState,
                               ParticleSubset* delset, DataWarehouse* old_dw)
{
  HypoElastic::allocateCMDataAdd(new_dw, addset, newState, delset, old_dw);

  constParticleVariable<double> o_StrainEnergyDensity;
  constParticleVariable<Matrix3> o_DispGrads;
  new_dw->get(o_StrainEnergyDensity, lb->pStrainEnergyDensityLabel_preReloc,
              delset);
  new_dw->get(o_DispGrads, lb->pDispGradsLabel_preReloc, delset);

  ParticleVariable<double> pStrainEnergyDensity;
  ParticleVariable<Matrix3> pDispGrads;
  new_dw->allocateTemporary(pDispGrads, addset);
  new_dw->allocateTemporary(pStrainEnergyDensity, addset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pDispGrads[*n] = o_DispGrads[*o];
    pStrainEnergyDensity[*n] = o_StrainEnergyDensity[*o];
  }

  (*newState)[lb->pDispGradsLabel] = pDispGrads.clone();
  (*newState)[lb->pStrainEnergyDensityLabel] = pStrainEnergyDensity.clone();
}

// Convert J-integral into stress intensity factors
// for Fracture
void
HypoElasticFracture::convertJToK(const MPMMaterial* matl, const string& stressState,
                         const Vector& J, const double& C, const Vector& D,
                         Vector& SIF)
{
  // J--J-integral vector,
  // C--Crack velocity,
  // D--COD near crack tip in local coordinates.

  double GT = std::abs(J.x()); // total energy release rate
  double CC = C * C;       // square of crack propagating velocity
  double D1 = D.y();
  double D2 = D.x();
  double D3 = D.z(); // D1,D2,D3: opening, sliding, tearing COD

  // Material properties
  double rho = matl->getInitialDensity();           // mass density
  double G = d_modelParam.G;                       // shear modulus
  double K = d_modelParam.K;                       // bulk modulus
  double v = 0.5 * (3. * K - 2. * G) / (3 * K + G); // Poisson ratio
  double k = (stressState == "planeStress") ? (3. - v) / (1. + v) : (3. - 4. * v);

  // Calculate stress intensity
  double KI, KII, KIII, K1 = 0, K2 = 0, K3 = 0;
  if (D1 == 0. && D2 == 0. && D3 == 0.) { // COD is zero
    KI = KII = KIII = 0.;
  } else { // COD is not zero
    // Parameters (A1,A2,A3) related to crack velocity
    double A1, A2, A3;
    if (std::sqrt(CC) < 1.e-16) { // for stationary crack
      A1 = (k + 1.) / 4.;
      A2 = (k + 1.) / 4.;
      A3 = 1.;
    } else { // for dynamic crack
      double Cs2 = G / rho;
      double Cd2 = (k + 1.) / (k - 1.) * Cs2;
      if (CC > Cs2)
        CC = Cs2;

      double B1 = std::sqrt(1. - CC / Cd2);
      double B2 = std::sqrt(1. - CC / Cs2);
      double DC = 4. * B1 * B2 - (1. + B2 * B2) * (1. + B2 * B2);
      A1 = B1 * (1. - B2 * B2) / DC;
      A2 = B2 * (1. - B2 * B2) / DC;
      A3 = 1. / B2;
    }

    // Solve stress intensity factors (absolute values)
    short CASE = 1;
    if (std::abs(D2) > std::abs(D1) && std::abs(D2) > std::abs(D3))
      CASE = 2;
    if (std::abs(D3) > std::abs(D1) && std::abs(D3) > std::abs(D2))
      CASE = 3;

    if (CASE == 1) { // Mode I COD is dominated
      double g21 = D2 / D1;
      double g31 = (1. - v) * D3 / D1;
      K1 = std::sqrt(2. * G * GT / (A1 + A2 * g21 * g21 + A3 * g31 * g31));
      K2 = std::abs(g21 * K1);
      K3 = std::abs(g31 * K1);
    }

    if (CASE == 2) { // Mode II COD is dominated
      double g12 = D1 / D2;
      double g32 = (1. - v) * D3 / D2;
      K2 = std::sqrt(2. * G * GT / (A1 * g12 * g12 + A2 + A3 * g32 * g32));
      K1 = std::abs(g12 * K2);
      K3 = std::abs(g32 * K2);
    }

    if (CASE == 3) { // Mode III COD is dominated
      double g13 = D1 / D3 / (1. - v);
      double g23 = D2 / D3 / (1. - v);
      K3 = std::sqrt(2. * G * GT / (A1 * g13 * g13 + A2 * g23 * g23 + A3));
      K1 = std::abs(g13 * K3);
      K2 = std::abs(g23 * K3);
    }

    // The signs of stress intensity are determined by the signs of the CODs
    double sign1 = D1 > 0. ? 1. : -1.;
    double sign2 = D2 > 0. ? 1. : -1.;
    double sign3 = D3 > 0. ? 1. : -1.;
    KI = D1 == 0. ? 0. : sign1 * K1;
    KII = D2 == 0. ? 0. : sign2 * K2;
    KIII = D3 == 0. ? 0. : sign3 * K3;
  }

  // Stress intensity vector
  SIF = Vector(KI, KII, KIII);
}

// Detect if crack propagates and the propagation direction
// for FRACTURE
short
HypoElasticFracture::crackPropagates(const double& Vc, const double& KI,
                             const double& KII, double& theta)
{
  /* Task 1: Determine fracture toughness KIc at velocity Vc
  */
  // Dynamic fracture toughness Kc(Vc,KIc,KIIC)
  std::vector<Toughness> Kc = d_fracParam.Kc;
  int num = (int)Kc.size();
  double KIc = -1., KIIc = -1.0;
  if (Vc <= Kc[0].Vc) { // Beyond the left bound
    KIc = Kc[0].KIc;
    KIIc = Kc[0].KIIc;
  } else if (Vc >= Kc[num - 1].Vc) { // Beyond the right bound
    KIc = Kc[num - 1].KIc;
    KIIc = Kc[num - 1].KIIc;
  } else { // In between
    for (int i = 0; i < num - 1; i++) {
      double Vi = Kc[i].Vc;
      double Vj = Kc[i + 1].Vc;
      if (Vc >= Vi && Vc < Vj) {
        double KIi = Kc[i].KIc;
        double KIj = Kc[i + 1].KIc;
        double KIIi = Kc[i].KIIc;
        double KIIj = Kc[i + 1].KIIc;
        KIc = KIi + (KIj - KIi) * (Vc - Vi) / (Vj - Vi);
        KIIc = KIIi + (KIIj - KIIi) * (Vc - Vi) / (Vj - Vi);
        break;
      }
    } // End of loop over i
  }

  /* Task 2: Determine crack propagation direction (theta) and
             the equivalent stress intensity factor (Kq)
  */
  double Kq = -9e32;
  if (d_crackPropagationCriterion == "max_hoop_stress" ||
      d_crackPropagationCriterion == "max_energy_release_rate") {
    // Crack propagation direction
    double sinTheta, cosTheta, value;
    if (KI == 0.0 || (KI != 0. && std::abs(KII / KI) > 1000.)) { // Pure mode II
      cosTheta = 1. / 3.;
      sinTheta = (KII >= 0.) ? -std::sqrt(8. / 9.) : sqrt(8. / 9.);
    } else { // Mixed mode or pure mode I
      double R = KII / KI;
      cosTheta = (3. * R * R + std::sqrt(1. + 8. * R * R)) / (1. + 9. * R * R);
      value = std::abs(R * (3. * cosTheta - 1.));
      sinTheta = (KII >= 0.) ? -value : value;
    }
    theta = std::asin(sinTheta);

    // Equivalent stress intensity
    double ct = std::cos(theta / 2.);
    double st = std::sin(theta / 2.);
    Kq = KI * std::pow(ct, 3) - 3 * KII * ct * ct * st;
  } // End of max_hoop_stress criterion

  if (d_crackPropagationCriterion == "max_principal_stress") {
    if (KII == 0. || (KII != 0. && std::abs(KI / KII) > 1000.)) { // Pure mode I
      theta = 0.;
      Kq = KI;
    } else { // Mixed mode or pure mode II
      double R = KI / KII;
      int sign = (KII > 0.) ? -1 : 1;
      double tanTheta2 = (R + sign * std::sqrt(R * R + 8.)) / 4.;
      double theta2 = std::atan(tanTheta2);
      double ct = std::cos(theta2);
      double st = std::sin(theta2);
      Kq = KI * std::pow(ct, 3) - 3 * KII * ct * ct * st;
      theta = 2 * theta2;
    }
  } // End of max_principal_stress criterion

  if (d_crackPropagationCriterion == "strain_energy_density") {
    // Calculate parameter k
    double k;
    std::string stressState = "planeStress";             // Plane stress
    double G = d_modelParam.G;                     // Shear modulus
    double K = d_modelParam.K;                     // Bulk modulus
    double v = 0.5 * (3 * K - 2 * G) / (3 * K + G); // Poisson ratio
    k = (stressState == "planeStress") ? (3. - v) / (1. + v) : (3. - 4. * v);

    // Crack propagation direction
    if (KII == 0.0 || (KII != 0. && std::abs(KI / KII) > 1000.)) { // Pure mode I
      theta = 0.0;
      Kq = KI;
    } else { // Mixed mode or pure mode II
      theta = crackPropagationAngleFromStrainEnergyDensityCriterion(k, KI, KII);
      // Equivalent stress intensity
      double ct = std::cos(theta), st = std::sin(theta);
      double a11 = (1 + ct) * (k - ct);
      double a12 = st * (2 * ct - k + 1);
      double a22 = (k + 1) * (1 - ct) + (1 + ct) * (3 * ct - 1);
      Kq = std::sqrt((a11 * KI * KI + 2 * a12 * KI * KII + a22 * KII * KII) / 2 /
                (k - 1));
    }
  } // End of strain_energy_density criterion

  if (d_crackPropagationCriterion == "empirical_criterion") {
    if (KII == 0. || (KII != 0. && std::abs(KI / KII) > 1000.)) { // Pure mode I
      theta = 0.;
      Kq = KI;
    } else { // For mixed mode or pure mode II, use maximum pricipal criterion
             // to
             // determine the crack propagation direction and the emprical
             // criterion to determine if crack will propagate.
      double R = KI / KII;
      int sign = (KII > 0.) ? -1 : 1;
      double tanTheta2 = (R + sign * sqrt(R * R + 8.)) / 4.;
      theta = 2 * std::atan(tanTheta2);
      Kq = (std::pow(KI / KIc, d_fracParam.p) + std::pow(KII / KIIc, d_fracParam.q)) * KIc;
    }
  } // End of empirical_criterion

  if (Kq >= KIc)
    return 1;
  else
    return 0;
}

// Obtain crack propagation angle numerically from strain energy density
// criterion
// for FRACTURE
double
HypoElasticFracture::crackPropagationAngleFromStrainEnergyDensityCriterion(
  const double& k, const double& KI, const double& KII)
{
  double errF = 1.e-6, errV = 1.e-2, PI = 3.141592654;
  double a, b, c, fa, fb, fc;

  double A = -PI, B = PI; // The region of the roots
  int n = 36;             // Divide [A,B] into n intervals
  double h = (B - A) / n; // Subinterval length

  double theta = 0.0;
  double theta0 = std::atan(KI / KII);
  std::vector<double> root; // Store the solutions of the equation
  // Solve the equation numerically
  for (int i = 0; i < n; i++) { // Loop over the whole interval [A,B]
    a = A + i * h;
    b = A + (i + 1) * h;
    fa = (k - 1) * std::sin(a - 2 * theta0) - 2 * std::sin(2 * (a - theta0)) - std::sin(2 * a);
    fb = (k - 1) * std::sin(b - 2 * theta0) - 2 * std::sin(2 * (b - theta0)) - std::sin(2 * b);

    // Find the root in [a,b)
    if (std::abs(fa) < errF) { // Where f(a)=0
      root.push_back(a);
    } else if (fa * fb < 0.) {          // There is a root in (a,b)
      double cp = 2 * B;                // Set the value beyond [A,B]
      for (int j = 0; j < 32768; j++) { // 32768=2^15 (a big int)
        c = b - (a - b) * fb / (fa - fb);
        fc = (k - 1) * std::sin(c - 2 * theta0) - 2 * std::sin(2 * (c - theta0)) -
             std::sin(2 * c);
        if (std::abs(fc) < errF || std::abs(c - cp) < errV) { // c is the root
          root.push_back(c);
          break;
        } else { // Record the cross point with axis x
          cp = c;
        }

        // Narrow the region of the root
        if (fc * fa < 0.) { // The root is in (a,c)
          fb = fc;
          b = c;
        } else if (fc * fb < 0.) { // The root is in (c,b)
          fa = fc;
          a = c;
        }
      } // End of loop over j
    }   // End of if(fa*fb<0.)
  }     // End of loop over i

  // Select the direction from the solutions
  // along which there exists the minimum strain energy density
  int count = 0;
  double S0 = 0.0;
  for (double r : root) {
    // The signs of propagation angle and KII must be opposite
    if (KII * r > 0.)
      continue;

    // Calculate the second derivative of the strain energy density
    double sr = std::sin(r), cr = std::cos(r), sr2 = std::sin(2 * r), cr2 = std::cos(2 * r);
    double dsdr2 = KI * KI * ((1 - k) * cr + 2 * cr2) -
                   2 * KI * KII * (4 * sr2 + (1 - k) * sr) +
                   KII * KII * ((k - 1) * cr - 6 * cr2);
    if (dsdr2 > 0.) {
      // Determine propagation angle by comparison of strain energy density.
      // Along the angle there exists the minimum strain energy density.
      double S = (1 + cr) * (k - cr) * KI * KI +
                 2 * sr * (2 * cr - k + 1) * KI * KII +
                 ((k + 1) * (1 - cr) + (1 + cr) * (3 * cr - 1)) * KII * KII;
      if (count == 0 || (count > 0 && S < S0)) {
        theta = r;
        S0 = S;
        count++;
      }
    }
  } // Enf of loop over i
  root.clear();

  return theta;
}

