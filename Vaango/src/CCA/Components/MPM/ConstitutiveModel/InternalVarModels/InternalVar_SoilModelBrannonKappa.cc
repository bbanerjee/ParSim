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

#include <CCA/Components/MPM/ConstitutiveModel/InternalVarModels/InternalVar_SoilModelBrannonKappa.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_SoilModelBrannon.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;

InternalVar_SoilModelBrannonKappa::InternalVar_SoilModelBrannonKappa(
  ProblemSpecP& ps)
{
  d_elastic = nullptr;
  d_shear = nullptr;

  ps->require("soil_model_brannon_p0", d_p0);
  ps->require("soil_model_brannon_p1", d_p1);
  ps->require("soil_model_brannon_p3", d_p3);
  ps->require("soil_model_brannon_p4", d_p4);
  ps->require("soil_model_brannon_B0", d_B0);
  ps->require("soil_model_brannon_Cr", d_Cr);
  ps->require("soil_model_brannon_fSlope", d_fSlope);
  ps->require("soil_model_brannon_peakI1", d_peakI1);

  // Initialize internal variable labels for evolution
  pKappaLabel =
    VarLabel::create("p.soil_model_brannon_kappa",
                     ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc =
    VarLabel::create("p.soil_model_brannon_kappa+",
                     ParticleVariable<double>::getTypeDescription());
}

InternalVar_SoilModelBrannonKappa::InternalVar_SoilModelBrannonKappa(
  const InternalVar_SoilModelBrannonKappa* cm)
{
  d_elastic = cm->d_elastic;
  d_shear = cm->d_shear;

  d_p0 = cm->d_p0;
  d_p1 = cm->d_p1;
  d_p3 = cm->d_p3;
  d_p4 = cm->d_p4;
  d_B0 = cm->d_B0;
  d_Cr = cm->d_Cr;
  d_fSlope = cm->d_fSlope;
  d_peakI1 = cm->d_peakI1;

  // Initialize internal variable labels for evolution
  pKappaLabel =
    VarLabel::create("p.soil_model_brannon_kappa",
                     ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc =
    VarLabel::create("p.soil_model_brannon_kappa+",
                     ParticleVariable<double>::getTypeDescription());
}

InternalVar_SoilModelBrannonKappa::~InternalVar_SoilModelBrannonKappa()
{
  VarLabel::destroy(pKappaLabel);
  VarLabel::destroy(pKappaLabel_preReloc);
}

void
InternalVar_SoilModelBrannonKappa::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type", "soil_model_brannon_kappa");

  int_var_ps->appendElement("soil_model_brannon_p0", d_p0);
  int_var_ps->appendElement("soil_model_brannon_p1", d_p1);
  int_var_ps->appendElement("soil_model_brannon_p3", d_p3);
  int_var_ps->appendElement("soil_model_brannon_p4", d_p4);
  int_var_ps->appendElement("soil_model_brannon_B0", d_B0);
  int_var_ps->appendElement("soil_model_brannon_Cr", d_Cr);
  int_var_ps->appendElement("soil_model_brannon_fSlope", d_fSlope);
  int_var_ps->appendElement("soil_model_brannon_peakI1", d_peakI1);
}

void
InternalVar_SoilModelBrannonKappa::addInitialComputesAndRequires(
  Task* task, const MPMMaterial* matl, const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pKappaLabel, matlset);
}

void
InternalVar_SoilModelBrannonKappa::addComputesAndRequires(
  Task* task, const MPMMaterial* matl, const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pKappaLabel, matlset, Ghost::None);
  task->computes(pKappaLabel_preReloc, matlset);
}

void
InternalVar_SoilModelBrannonKappa::addParticleState(
  std::vector<const VarLabel*>& from, std::vector<const VarLabel*>& to)
{
  from.push_back(pKappaLabel);
  to.push_back(pKappaLabel_preReloc);
}

void
InternalVar_SoilModelBrannonKappa::allocateCMDataAddRequires(
  Task* task, const MPMMaterial* matl, const PatchSet*, MPMLabel*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pKappaLabel_preReloc, matlset, Ghost::None);
}

void
InternalVar_SoilModelBrannonKappa::allocateCMDataAdd(
  DataWarehouse* old_dw, ParticleSubset* addset,
  ParticleLabelVariableMap* newState, ParticleSubset* delset,
  DataWarehouse* new_dw)
{
  ParticleVariable<double> pKappa;
  constParticleVariable<double> o_kappa;

  new_dw->allocateTemporary(pKappa, addset);

  new_dw->get(o_kappa, pKappaLabel_preReloc, delset);

  ParticleSubset::iterator o, n = addset->begin();
  for (o = delset->begin(); o != delset->end(); o++, n++) {
    pKappa[*n] = o_kappa[*o];
  }

  (*newState)[pKappaLabel] = pKappa.clone();
}

void
InternalVar_SoilModelBrannonKappa::initializeInternalVariable(
  ParticleSubset* pset, DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<double> pKappa;
  new_dw->allocateAndPut(pKappa, pKappaLabel, pset);
  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    pKappa[*iter] =
      (d_p0 + d_Cr * d_fSlope * d_peakI1) / (d_Cr * d_fSlope + 1.0);
  }
}

void
InternalVar_SoilModelBrannonKappa::getInternalVariable(
  ParticleSubset* pset, DataWarehouse* old_dw,
  constParticleVariableBase& pKappa)
{
  old_dw->get(pKappa, pKappaLabel, pset);
}

void
InternalVar_SoilModelBrannonKappa::allocateAndPutInternalVariable(
  ParticleSubset* pset, DataWarehouse* new_dw, ParticleVariableBase& pKappa_new)
{
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc, pset);
}

void
InternalVar_SoilModelBrannonKappa::allocateAndPutRigid(
  ParticleSubset* pset, DataWarehouse* new_dw,
  constParticleVariableBase& pKappa)
{
  ParticleVariable<double> pKappa_new;
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc, pset);
  ParticleSubset::iterator iter = pset->begin();
  for (; iter != pset->end(); iter++) {
    pKappa_new[*iter] =
      dynamic_cast<constParticleVariable<double>&>(pKappa)[*iter];
  }
}

//--------------------------------------------------------------------------------------
// Compute kappa_new using Newton's method
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeInternalVariable(
  const ModelStateBase* state_input) const
{
  const ModelState_SoilModelBrannon* state =
    static_cast<const ModelState_SoilModelBrannon*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_SoilModelBrannon.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Get the local variables needed
  double kappa_old = state->kappa;         // old value of kappa may have
                                           // been modified before
  double cap_radius = state->CR;           // C_radius
  double max_kappa = state->maxX;          // X_max
  double eps_v = state->eps_v;             // eps_v_p
  double delta_eps_v = state->delta_eps_v; // Delta eps_v_p
  double scale_fac = state->scale_eps_v;   // factor = 1+f_s c_r or 1

  // Scale the volumetric plastic strain
  delta_eps_v /= scale_fac;

  // Subtract cap_radius from kappa and init new kappa
  kappa_old -= cap_radius;
  double kappa_new = kappa_old;

  // Calculate new kappa
  double tolerance = 1.0e-3;
  int maxiter = 10;

  if (kappa_old - d_p0 < 0.0) {
    // Update \kappa in the cae of X < p0
    // (see "fig:Arenisca_YieldSurface" in the Arenisca manual)
    // (see "eq:evolutionOfKappaFluidEffect" in the Arenisca manual)
    kappa_new =
      computeKappaFromX1(kappa_old, eps_v, delta_eps_v, tolerance, maxiter);
    // if (std::isnan(kappa_new) || kappa_new > 0) {
    if (std::isnan(kappa_new)) {
      std::cerr << " kappa_new = " << kappa_new << " kappa_old = " << kappa_old
           << " eps_v = " << eps_v << " delta_eps_v = " << delta_eps_v << "\n";
      throw InvalidValue("**ERROR**: Nan in kappa_new - case 1", __FILE__,
                         __LINE__);
    }
  } else if (kappa_old < max_kappa) {
    // Update \kappa in the cae of p0 <= X < max_X
    // (see "fig:Arenisca_YieldSurface" in the Arenisca manual)
    // (see "eq:evolutionOfKappaFluidEffect1" in the Arenisca manual)
    // (for the limitation of max_X see "eq:limitationForX" in the Arenisca
    // manual)
    kappa_new =
      computeKappaFromX2(kappa_old, eps_v, delta_eps_v, tolerance, maxiter);
    // if (std::isnan(kappa_new) || kappa_new > 0) {
    if (std::isnan(kappa_new)) {
      std::cerr << " kappa_new = " << kappa_new << " kappa_old = " << kappa_old
           << " eps_v = " << eps_v << " delta_eps_v = " << delta_eps_v << "\n";
      throw InvalidValue("**ERROR**: Nan in kappa_new - case 2", __FILE__,
                         __LINE__);
    }
  } else {
    // Update \kappa in the cae of X >= max_X
    // (see "fig:Arenisca_YieldSurface" in the Arenisca manual)
    // In this case it is assumed that X=max_X
    // (for the limitation of max_X see "eq:limitationForX" in the Arenisca
    // manual)
    // pKappaState is a particle variable variable which defines if the particle
    // meet any of the limitation for \kappa and X or not?
    // pKappaState=1: means that the particle met the max_X limitation
    kappa_new = max_kappa;
  }

  // Add back the cap_radius
  kappa_new += cap_radius;

  // return the new kappa
  return kappa_new;
}

//--------------------------------------------------------------------------------------
// Compute kappa_new from the function X1(kappa_{n+1})
//  where
//       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//
// ** NOTE** (should be replaced with function pointers)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeKappaFromX1(const double& kappa_old,
                                                      const double& epsv,
                                                      const double& deltaEpsv,
                                                      const double& tolerance,
                                                      const int& maxiter) const
{
  // Compute B, G and H
  double B = computeB(); // This should be computed while creating the oebject
  double G = computeG(epsv, B);
  double H = computeH(epsv, B);

  // Calculate new kappa
  double kappa_new_iter = kappa_old;
  double kappa_new_iter_old = kappa_new_iter;
  double X1 = 0.0;
  double dX1dkappa = 0.0;
  int iter = 0;
  do {
    X1 = computeX1(kappa_old, kappa_new_iter, G, H, deltaEpsv);
    dX1dkappa = computeDerivX1dkappa(kappa_old, kappa_new_iter, deltaEpsv);
    kappa_new_iter_old = kappa_new_iter;
    kappa_new_iter -= X1 / dX1dkappa;
    ++iter;
  } while ((fabs(kappa_new_iter - kappa_new_iter_old) > tolerance) &&
           (iter < maxiter));
  kappa_new_iter_old = kappa_new_iter;
  if (!(iter < maxiter) || std::isnan(kappa_new_iter)) {
    kappa_new_iter = computeKappaAtX1Min(deltaEpsv);
    if (std::isnan(kappa_new_iter) & !(std::isnan(kappa_new_iter_old)))
      kappa_new_iter = kappa_new_iter_old;
    // std::cerr << "Func 1: kappa_new = " << kappa_new_iter << "\n";
    // std::cerr << "    epsv = " << epsv << " deltaEpsv = " << deltaEpsv << "
    // kappa_old = " << kappa_old << "\n";
  }
  // if (std::isnan(kappa_new_iter)) {
  //   std::cerr << " iter = " << iter << " maxiter = " << maxiter << " kappa[k+1] =
  //   " << kappa_new_iter
  //        << " kappa[k] = " << kappa_new_iter_old << " X1 = " << X1
  //        << " dX1dkappa = " << dX1dkappa << "\n";
  //   std::cerr << "    epsv = " << epsv << " deltaEpsv = " << deltaEpsv << "
  //   kappa_old = " << kappa_old
  //        << " kappa[k+1] - kappa[k] = " << fabs(kappa_new_iter -
  //        kappa_new_iter_old)
  //        << " tol = " << tolerance << "\n";
  // }

  return kappa_new_iter;
}

//--------------------------------------------------------------------------------------
// Compute kappa_new from the function X2(kappa_{n+1})
//  where
//       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//
// ** NOTE** (should be replaced with function pointers)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeKappaFromX2(const double& kappa_old,
                                                      const double& epsv,
                                                      const double& deltaEpsv,
                                                      const double& tolerance,
                                                      const int& maxiter) const
{
  // Compute B, G and H
  double B = computeB(); // This should be computed while creating the oebject
  double G = computeG(epsv, B);
  double H = computeH(epsv, B);

  // Calculate new kappa
  double kappa_new_iter = kappa_old;
  double kappa_new_iter_old = kappa_new_iter;
  int iter = 0;
  do {
    double X2 = computeX2(kappa_old, kappa_new_iter, G, H, deltaEpsv);
    double dX2dkappa =
      computeDerivX2dkappa(kappa_old, kappa_new_iter, deltaEpsv);
    kappa_new_iter_old = kappa_new_iter;
    kappa_new_iter -= X2 / dX2dkappa;
    ++iter;
    if (std::isnan(kappa_new_iter) || kappa_new_iter > 0) {
      std::cerr << " iter = " << iter << " maxiter = " << maxiter
           << " kappa[k+1] = " << kappa_new_iter
           << " kappa[k] = " << kappa_new_iter_old << " X2 = " << X2
           << " dX2dkappa = " << dX2dkappa << "\n";
      std::cerr << "    epsv = " << epsv << " deltaEpsv = " << deltaEpsv
           << " kappa[k+1] - kappa[k] "
           << fabs(kappa_new_iter - kappa_new_iter_old)
           << " tol = " << tolerance << "\n";
    }
  } while ((fabs(kappa_new_iter - kappa_new_iter_old) > tolerance) &&
           (iter < maxiter));

  return kappa_new_iter;
}

//--------------------------------------------------------------------------------------
// Compute the function X1(kappa_{n+1})
//  where
//       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeX1(const double& kappa_old,
                                             const double& kappa_new,
                                             const double& G, const double& H,
                                             const double& delEpsv) const
{
  double F1 = computeF1(kappa_new, G, H);
  double X1 = kappa_new - kappa_old - F1 * delEpsv;
  return X1;
}

//--------------------------------------------------------------------------------------
// Compute the function dX1/dkappa(kappa_{n+1})
//  where
//       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeDerivX1dkappa(
  const double& kappa_old, const double& kappa_new, const double& delEpsv) const
{
  double dF1dkappa = computeDerivF1dkappa(kappa_new);
  double dX1dkappa = 1.0 - dF1dkappa * delEpsv;
  return dX1dkappa;
}

//--------------------------------------------------------------------------------------
// Compute the value of kappa at which function X1 is a minimum
//  where
//       X1(kappa_{n+1}) = kappa_{n+1} - kappa_n - F1(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeKappaAtX1Min(
  const double& delEpsv) const
{
  double kappa = d_p0 - 1.0 / d_p1 * log(-d_p3 / delEpsv);
  return kappa;
}

//--------------------------------------------------------------------------------------
// Compute the function X2(kappa_{n+1})
//  where
//       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeX2(const double& kappa_old,
                                             const double& kappa_new,
                                             const double& G, const double& H,
                                             const double& delEpsv) const
{
  double F2 = computeF2(kappa_new, G, H);
  double X2 = kappa_new - kappa_old - F2 * delEpsv;
  return X2;
}

//--------------------------------------------------------------------------------------
// Compute the function dX2/dkappa(kappa_{n+1})
//  where
//       X2(kappa_{n+1}) = kappa_{n+1} - kappa_n - F2(kappa_{n+1},epsv_{n+1})
//       Delta epsv = 0
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeDerivX2dkappa(
  const double& kappa_old, const double& kappa_new, const double& delEpsv) const
{
  double dF2dkappa = computeDerivF2dkappa(kappa_new);
  double dX2dkappa = 1.0 - dF2dkappa * delEpsv;
  return dX2dkappa;
}

//--------------------------------------------------------------------------------------
// Compute the constant B
//  where
//        B = 3 B0 [exp(p3+p4) - 1]
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeB() const
{
  double B = 3.0 * d_B0 * (exp(d_p3 + d_p4) - 1.0);
  return B;
}

//--------------------------------------------------------------------------------------
// Compute the function G(epsv)
//  where
//        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
//        B = 3 B0 [exp(p3+p4) - 1]
//        g34 = exp(p3+p4+epsv)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeG(const double& epsv,
                                            const double& B) const
{
  double g34 = exp(d_p3 + d_p4 + epsv);
  double G = B * g34 / ((g34 - 1.0) * (g34 - 1.0));
  return G;
}

//--------------------------------------------------------------------------------------
// Compute the function H(epsv)
//  where
//        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
//        B = 3 B0 [exp(p3+p4) - 1]
//        h3 = exp(p3+epsv)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeH(const double& epsv,
                                            const double& B) const
{
  double h3 = exp(d_p3 + epsv);
  double H = B * h3 / ((h3 - 1.0) * (h3 - 1.0));
  return H;
}

//--------------------------------------------------------------------------------------
// Compute the function F1(kappa, epsv) = f1(kappa) - G(epsv) + H(epsv)
//  where f1(kappa) = 1/(p1 p3) exp(-p1 kappa - p0)
//        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
//        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
//        B = 3 B0 [exp(p3+p4) - 1]
//        g34 = exp(p3+p4+epsv)
//        h3 = exp(p3+epsv)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeF1(const double& kappa,
                                             const double& G,
                                             const double& H) const
{
  double f1 = 1.0 / (d_p1 * d_p3) * exp(-d_p1 * (kappa - d_p0));
  return (f1 - G + H);
}

//--------------------------------------------------------------------------------------
// Compute the function dF1/dkappa(kappa, epsv) = df1/dkappa(kappa)
//  where f1(kappa) = 1/(p1 p3) exp(-p1 kappa - p0)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeDerivF1dkappa(
  const double& kappa) const
{
  double df1dkappa = -1.0 / d_p3 * exp(-d_p1 * (kappa - d_p0));
  return df1dkappa;
}

//--------------------------------------------------------------------------------------
// Compute the function F2(kappa, epsv) = f2(kappa) - G(epsv) + H(epsv)
//  where f2(kappa) = 1/(p1 p3) [kappa/p0]^(1-p0p1p3)
//        G(epsv) = B g34(epsv)/[g34(epsv) - 1]^2
//        H(epsv) = B h3(epsv)/[h3(epsv) - 1]^2
//        B = 3 B0 [exp(p3+p4) - 1]
//        g34 = exp(p3+p4+epsv)
//        h3 = exp(p3+epsv)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeF2(const double& kappa,
                                             const double& G,
                                             const double& H) const
{
  double p1p3 = d_p1 * d_p3;
  double p0p1p3 = d_p0 * p1p3;
  double kappa_p0 = kappa / d_p0;
  if (fabs(kappa_p0 - 1.0) < 1.0e-10)
    kappa_p0 = 1.0;

  // feclearexcept(FE_ALL_EXCEPT);
  // double pow_kappa_p0 = pow(kappa_p0, p0p1p3);
  double pow_kappa_p0 = kappa_p0 / pow(kappa_p0, p0p1p3);
  // std::cerr << "Location 3: Floating point exception in particle ? " << hex <<
  // fetestexcept(FE_ALL_EXCEPT) << "\n";
  // std::cerr << setiosflags(ios::fixed) << setprecision(20) << scientific << "\n";
  // std::cerr << " kappa0 - p0 " << kappa - d_p0 << "\n";
  // std::cerr << " kappa/p0 = " << kappa_p0 << " p0p1p3 = " << p0p1p3 << " power = "
  // << pow_kappa_p0 << "\n";

  // double lnb = log(kappa_p0);
  // double clnb = p0p1p3*lnb;
  // pow_kappa_p0 = exp(clnb);
  // std::cerr << "Location 3: Floating point exception in particle ? " << hex <<
  // fetestexcept(FE_ALL_EXCEPT) << "\n";
  // std::cerr << setiosflags(ios::fixed) << setprecision(10) << " kappa/p0 = " <<
  // kappa_p0 << " p0p1p3 = " << p0p1p3 << " lnb = " << lnb
  //     << " c lnb = " << clnb << " power = " << pow_kappa_p0 << "\n";

  double f2 = (1.0 / p1p3) * (double)pow_kappa_p0;
  return (f2 - G + H);
}

//--------------------------------------------------------------------------------------
// Compute the function dF2/dkappa(kappa, epsv) = df2/dkappa(kappa)
//  where f2(kappa) = 1/(p1 p3) [kappa/p0]^(1-p0p1p3)
//--------------------------------------------------------------------------------------
double
InternalVar_SoilModelBrannonKappa::computeDerivF2dkappa(
  const double& kappa) const
{
  double p0p1p3 = d_p0 * d_p1 * d_p3;
  double df2dkappa = (1.0 / p0p1p3 - 1.0) * pow((kappa / d_p0), (-p0p1p3));
  return df2dkappa;
}
