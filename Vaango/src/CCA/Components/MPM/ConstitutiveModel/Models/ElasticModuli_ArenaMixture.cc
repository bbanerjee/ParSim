/*
 * The MIT License
 *
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_ArenaMixture.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Arena.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Uintah;
using namespace Vaango;

// Construct a default elasticity model.
/* From Arenisca3:
 * If the user has specified a nonzero nu1 and nu2, these are used to define a
 * pressure
 * dependent Poisson ratio, which is used to adjust the shear modulus along with
 * the
 * bulk modulus.  The high pressure limit has nu=nu1+nu2;
 * if ((d_cm.nu1!=0.0)&&(d_cm.nu2!=0.0)){
 *     // High pressure bulk modulus:
 *     double nu = d_cm.nu1+d_cm.nu2;
 *     shear = 1.5*bulk*(1.0-2.0*nu)/(1.0+nu);
 * }
 */

ElasticModuli_ArenaMixture::ElasticModuli_ArenaMixture(Uintah::ProblemSpecP& ps)
{
  ps->require("vol_frac.phase1", d_volfrac[0]); // Volume fractions
  d_volfrac[1] = 1.0 - d_volfrac[0];

  ps->require("b0.phase1",
              d_bulk[0].b0); // Tangent Elastic Bulk Modulus Parameters
  ps->require("b1.phase1", d_bulk[0].b1); // Phase1
  ps->require("b2.phase1", d_bulk[0].b2);
  ps->require("b3.phase1", d_bulk[0].b3);
  ps->require("b4.phase1", d_bulk[0].b4);

  ps->require("G0.phase1",
              d_shear[0].G0); // Tangent Elastic Shear Modulus Parameters
  ps->require("nu1.phase1",
              d_shear[0].nu1); // Phase 1// Low pressure Poisson ratio
  ps->require("nu2.phase1",
              d_shear[0].nu2); // Pressure-dependent Poisson ratio term

  ps->require("b0.phase2",
              d_bulk[1].b0); // Tangent Elastic Bulk Modulus Parameters
  ps->require("b1.phase2", d_bulk[1].b1); // Phase 2
  ps->require("b2.phase2", d_bulk[1].b2);
  ps->require("b3.phase2", d_bulk[1].b3);
  ps->require("b4.phase2", d_bulk[1].b4);

  ps->require("G0.phase2",
              d_shear[1].G0); // Tangent Elastic Shear Modulus Parameters
  ps->require("nu1.phase2", d_shear[1].nu1); // Phase 2
  ps->require("nu2.phase2", d_shear[1].nu2);

  checkInputParameters();
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
ElasticModuli_ArenaMixture::checkInputParameters()
{
  std::ostringstream warn;

  if (d_bulk[0].b0 <= 0.0) {
    warn << "Phase 1: b0 must be positive. b0 = " << d_bulk[0].b0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[0].b1 < 0.0) {
    warn << "Phase 1: b1 must be nonnegative. b1 = " << d_bulk[0].b1
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[0].b2 < 0.0) {
    warn << "Phase 1: b2 must be nonnegative. b2 = " << d_bulk[0].b2
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[0].b3 < 0.0) {
    warn << "Phase 1: b3 must be nonnegative. b3 = " << d_bulk[0].b3
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[0].b4 < 1.0) {
    warn << "Phase 1: b4 must be >= 1. b4 = " << d_bulk[0].b4 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_shear[0].G0 <= 0.0) {
    warn << "Phase 1: G0 must be positive. G0 = " << d_shear[0].G0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[1].b0 <= 0.0) {
    warn << "Phase 2: b0 must be positive. b0 = " << d_bulk[1].b0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[1].b1 < 0.0) {
    warn << "Phase 2: b1 must be nonnegative. b1 = " << d_bulk[1].b1
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[1].b2 < 0.0) {
    warn << "Phase 2: b2 must be nonnegative. b2 = " << d_bulk[1].b2
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[1].b3 < 0.0) {
    warn << "Phase 2: b3 must be nonnegative. b3 = " << d_bulk[1].b3
         << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_bulk[1].b4 < 1.0) {
    warn << "Phase 2: b4 must be >= 1. b4 = " << d_bulk[1].b4 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_shear[1].G0 <= 0.0) {
    warn << "Phase 2: G0 must be positive. G0 = " << d_shear[1].G0 << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
  if (d_volfrac[0] < 0.0 || d_volfrac[0] > 1.0) {
    warn << "Phase 1: Volume fraction must be between 0 and 1.  vf_phase1 = "
         << d_volfrac[0] << std::endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_ArenaMixture::ElasticModuli_ArenaMixture(
  const ElasticModuli_ArenaMixture* model)
{
  for (int ii = 0; ii < 2; ii++) {
    d_volfrac[ii] = model->d_volfrac[ii];
    d_bulk[ii] = model->d_bulk[ii];
    d_shear[ii] = model->d_shear[ii];
  }
}

// Destructor of elasticity model.
ElasticModuli_ArenaMixture::~ElasticModuli_ArenaMixture() = default;

void
ElasticModuli_ArenaMixture::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "arena_mixture");

  elasticModuli_ps->appendElement("vol_frac.phase1", d_volfrac[0]);

  elasticModuli_ps->appendElement("b0.phase1", d_bulk[0].b0);
  elasticModuli_ps->appendElement("b1.phase1", d_bulk[0].b1);
  elasticModuli_ps->appendElement("b2.phase1", d_bulk[0].b2);
  elasticModuli_ps->appendElement("b3.phase1", d_bulk[0].b3);
  elasticModuli_ps->appendElement("b4.phase1", d_bulk[0].b4);

  elasticModuli_ps->appendElement("G0.phase1", d_shear[0].G0);
  elasticModuli_ps->appendElement("nu1.phase1", d_shear[0].nu1);
  elasticModuli_ps->appendElement("nu2.phase1", d_shear[0].nu2);

  elasticModuli_ps->appendElement("b0.phase2", d_bulk[1].b0);
  elasticModuli_ps->appendElement("b1.phase2", d_bulk[1].b1);
  elasticModuli_ps->appendElement("b2.phase2", d_bulk[1].b2);
  elasticModuli_ps->appendElement("b3.phase2", d_bulk[1].b3);
  elasticModuli_ps->appendElement("b4.phase2", d_bulk[1].b4);

  elasticModuli_ps->appendElement("G0.phase2", d_shear[1].G0);
  elasticModuli_ps->appendElement("nu1.phase2", d_shear[1].nu1);
  elasticModuli_ps->appendElement("nu2.phase2", d_shear[1].nu2);
}

// Compute the elastic moduli
ElasticModuli
ElasticModuli_ArenaMixture::getInitialElasticModuli() const
{
  double Ks = d_granite.getBulkModulus();
  double KK0 = Ks * d_bulk[0].b0;
  double KK1 = Ks * d_bulk[1].b0;
  double KRatio0 = KK0 / Ks;
  double KRatio1 = KK1 / Ks;
  double nu0 = d_shear[0].nu1 + d_shear[0].nu2 * exp(-KRatio0);
  double nu1 = d_shear[1].nu1 + d_shear[1].nu2 * exp(-KRatio1);
  double GG0 =
    (nu0 > 0.0) ? 1.5 * KK0 * (1.0 - 2.0 * nu0) / (1.0 + nu0) : d_shear[0].G0;
  double GG1 =
    (nu1 > 0.0) ? 1.5 * KK1 * (1.0 - 2.0 * nu1) / (1.0 + nu1) : d_shear[1].G0;
  double K_mix =
    Ks / (d_volfrac[0] / d_bulk[0].b0 + d_volfrac[1] / d_bulk[1].b0);
  double G_mix = 1.0 / (d_volfrac[0] / GG0 + d_volfrac[1] / GG1);
  return ElasticModuli(K_mix, G_mix);
}

ElasticModuli
ElasticModuli_ArenaMixture::getCurrentElasticModuli(
  const ModelStateBase* state_input)
{
  const ModelState_Arena* state =
    static_cast<const ModelState_Arena*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Arena.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

  // Make sure the quantities are positive in compression
  double I1_eff_bar = -state->I1_eff;
  double ev_p_bar = -(state->plasticStrainTensor).Trace();
  double pw_bar = state->pbar_w;

  // Compute the elastic moduli
  double K_mix = 0.0, G_mix = 0.0;

  // *Note*
  // 1) If the mean stress is tensile (I1_eff_bar <= 0.0), the soil has the
  // properties of
  //    drained soil
  // 2) If the saturation is > 0.0 and the mean stress is compressive, the soil
  // has
  //    properties of the partially saturated soil
  if (state->saturation > 0.0) {
    // Partially saturated material
    double phi = state->porosity;
    double Sw = state->saturation;
    computePartialSaturatedModuli(I1_eff_bar, pw_bar, ev_p_bar, phi, Sw, K_mix,
                                  G_mix);
#ifdef DEBUG_BULK_MODULUS
    std::cout << "Computing bulk modulus for saturated material:" << std::endl;
    std::cout << "  phi = " << phi << " Sw = " << Sw
              << " I1bar = " << I1_eff_bar << " K = " << K_mix
              << " G = " << G_mix << std::endl;
#endif
  } else {
    // Drained material
    computeDrainedModuli(I1_eff_bar, ev_p_bar, K_mix, G_mix);
#ifdef DEBUG_BULK_MODULUS
    std::cout << " I1bar = " << I1_eff_bar << " K = " << K_mix
              << " G = " << G_mix << std::endl;
#endif
  }

  return ElasticModuli(K_mix, G_mix);
}

void
ElasticModuli_ArenaMixture::computeDrainedModuli(const double& I1_eff_bar,
                                                 const double& ev_p_bar,
                                                 double& K_mix, double& G_mix)
{
  double KK[2] = { 0.0, 0.0 };
  double GG[2] = { 0.0, 0.0 };
  computeDrainedModuli(0, I1_eff_bar, ev_p_bar, KK[0], GG[0]);
  computeDrainedModuli(1, I1_eff_bar, ev_p_bar, KK[1], GG[1]);
  K_mix = 1.0 / (d_volfrac[0] / KK[0] + d_volfrac[1] / KK[1]);
  G_mix = 1.0 / (d_volfrac[0] / GG[0] + d_volfrac[1] / GG[1]);
  return;
}

void
ElasticModuli_ArenaMixture::computePartialSaturatedModuli(
  const double& I1_eff_bar, const double& pw_bar, const double& ev_p_bar,
  const double& phi, const double& S_w, double& K_mix, double& G_mix)
{
  double KK[2] = { 0.0, 0.0 };
  double GG[2] = { 0.0, 0.0 };
  computePartialSaturatedModuli(0, I1_eff_bar, pw_bar, ev_p_bar, phi, S_w,
                                KK[0], GG[0]);
  computePartialSaturatedModuli(1, I1_eff_bar, pw_bar, ev_p_bar, phi, S_w,
                                KK[1], GG[1]);
  K_mix = 1.0 / (d_volfrac[0] / KK[0] + d_volfrac[1] / KK[1]);
  G_mix = 1.0 / (d_volfrac[0] / GG[0] + d_volfrac[1] / GG[1]);
  return;
}

void
ElasticModuli_ArenaMixture::computeDrainedModuli(int phase,
                                                 const double& I1_eff_bar,
                                                 const double& ev_p_bar,
                                                 double& KK, double& GG)
{
  if (I1_eff_bar > 0.0) { // Compressive mean stress

    double pressure = I1_eff_bar / 3.0;

    // Compute solid matrix bulk modulus
    double K_s = d_granite.computeBulkModulus(pressure);
    double ns = d_granite.computeDerivBulkModulusPressure(pressure);
    double KsRatio = K_s / (1.0 - ns * pressure / K_s);

    // Compute Ev_e
    double ev_e =
      std::pow((d_bulk[phase].b3 * pressure) /
                 (d_bulk[phase].b1 * K_s - d_bulk[phase].b2 * pressure),
               (1.0 / d_bulk[phase].b4));

    // Compute y, z
    double y = std::pow(ev_e, d_bulk[phase].b4);
    double z = d_bulk[phase].b2 * y + d_bulk[phase].b3;

    // Compute compressive bulk modulus
    KK = KsRatio * (d_bulk[phase].b0 +
                    (1 / ev_e) * d_bulk[phase].b1 * d_bulk[phase].b3 *
                      d_bulk[phase].b4 * y / (z * z));

    // Update the shear modulus (if needed, i.e., nu1 & nu2 > 0)
    double KRatio = KK / K_s;
    double nu = d_shear[phase].nu1 + d_shear[phase].nu2 * exp(-KRatio);
    GG =
      (nu > 0.0) ? 1.5 * KK * (1.0 - 2.0 * nu) / (1.0 + nu) : d_shear[phase].G0;

  } else {

    // Tensile bulk modulus = Bulk modulus at p = 0
    double K_s0 = d_granite.computeBulkModulus(0.0);
    KK = d_bulk[phase].b0 * K_s0;
    double KRatio = KK / K_s0;

    // Tensile shear modulus
    double nu = d_shear[phase].nu1 + d_shear[phase].nu2 * exp(-KRatio);
    GG =
      (nu > 0.0) ? 1.5 * KK * (1.0 - 2.0 * nu) / (1.0 + nu) : d_shear[phase].G0;
  }

  return;
}

void
ElasticModuli_ArenaMixture::computePartialSaturatedModuli(
  int phase, const double& I1_eff_bar, const double& pw_bar,
  const double& ev_p_bar, const double& phi, const double& S_w, double& KK,
  double& GG)
{
  if (I1_eff_bar > 0.0) { // Compressive mean stress

    // Bulk modulus of grains
    double pressure = I1_eff_bar / 3.0;
    double K_s = d_granite.computeBulkModulus(pressure);

    // Bulk modulus of air
    double K_a = d_air.computeBulkModulus(pw_bar);

    // Bulk modulus of water
    double K_w = d_water.computeBulkModulus(pw_bar);

    // Bulk modulus of drained material
    double K_d = 0.0;
    GG = 0.0;
    computeDrainedModuli(phase, I1_eff_bar, ev_p_bar, K_d, GG);

    // Bulk modulus of air + water mixture
    double K_f = 1.0 / (S_w / K_w + (1.0 - S_w) / K_a);

    // Bulk modulus of partially saturated material (Biot-Grassman model)
    double numer = (1.0 - K_d / K_s) * (1.0 - K_d / K_s);
    double denom =
      1.0 / K_s * (1.0 - K_d / K_s) + phi * (1.0 / K_f - 1.0 / K_s);
    KK = K_d + numer / denom;

  } else { // Tensile mean stress

    computeDrainedModuli(phase, I1_eff_bar, ev_p_bar, KK, GG);
  }

  return;
}
