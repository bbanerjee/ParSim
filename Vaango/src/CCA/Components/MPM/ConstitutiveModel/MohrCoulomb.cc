/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 Center for the Simulation of Accidental Fires and
 * Explosions (CSAFE), and  Scientific Computing and Imaging Institute (SCI),
 * University of Utah.
 * Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
 *
 * License for the specific language governing rights and limitations under
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombClassic.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/MohrCoulombSheng.h>
#include <CCA/Components/MPM/ConstitutiveModel/MohrCoulomb.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Exceptions/InvalidValue.h>
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
#include <Core/Math/Matrix3.h>
#include <Core/Math/Short27.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <Core/Malloc/Allocator.h>
#include <Core/Math/MinMax.h>
#include <Core/Util/DebugStream.h>

#include <sci_defs/uintah_defs.h>

#include <iostream>
#include <string>

using namespace Uintah;

static DebugStream dbg_doing("MohrCoulomb_doing", false);
static DebugStream dbg("MohrCoulomb", false);

/**
 * Constructors
 */
MohrCoulomb::MohrCoulomb(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  // Read model parameters from the input file
  getInputParameters(ps);

  // Create the model
  if (d_modelType == "classic" || d_modelType == "classic_semiimplicit") {

    d_modelP = std::make_unique<MohrCoulombClassic>(d_params.G,
                                                    d_params.K,
                                                    d_params.c,
                                                    d_params.phi,
                                                    d_params.psi,
                                                    d_params.pMin);
  } else if (d_modelType == "sheng") {

    d_modelP = std::make_unique<MohrCoulombSheng>(d_params.G,
                                                  d_params.K,
                                                  d_params.c,
                                                  d_params.phi,
                                                  d_params.psi,
                                                  d_params.pMin);
  } else {
    std::ostringstream err;
    err << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  // Set the integration parameters
  setIntegrationParameters();

  // Set up internal variables that can vary with each particle
  initializeLocalMPMLabels();
}

/**
 * Copy constructor
 */
MohrCoulomb::MohrCoulomb(const MohrCoulomb* cm)
  : ConstitutiveModel(cm)
{
  *this = cm;
}

MohrCoulomb::MohrCoulomb(const MohrCoulomb& cm)
  : ConstitutiveModel(cm)
{
  *this = cm;
}

MohrCoulomb&
MohrCoulomb::operator=(const MohrCoulomb& cm)
{
  dbg_doing << "Doing MohrCoulomb::operator=\n";

  if (this != &cm) {
    d_modelType = cm.d_modelType;
    if (d_modelType == "classic" || d_modelType == "classic_semiimplicit") {

      d_modelP = std::make_unique<MohrCoulombClassic>(
        static_cast<MohrCoulombClassic*>(cm.d_modelP.get()));

    } else if (d_modelType == "sheng") {

      d_modelP = std::make_unique<MohrCoulombSheng>(
        static_cast<MohrCoulombSheng*>(cm.d_modelP.get()));

    } else {
      std::ostringstream err;
      err
        << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
      throw ProblemSetupException(err.str(), __FILE__, __LINE__);
    }

    d_flags  = cm.d_flags;
    d_params = cm.d_params;
    d_int    = cm.d_int;

    initializeLocalMPMLabels();
  }
  return *this;
}

/*
 * Create clone
 */
MohrCoulomb*
MohrCoulomb::clone()
{
  return scinew MohrCoulomb(*this);
}

/**
 * Get parameters
 */
void
MohrCoulomb::getInputParameters(ProblemSpecP& ps)
{
  dbg_doing << "Doing MohrCoulomb::getInputParameters\n";

  d_modelType = "classic"; // options: classic, classic_semiimplicit, sheng
  ps->require("model_type", d_modelType);

  getModelParameters(ps);
  getIntegrationParameters(ps);
}

void
MohrCoulomb::getModelParameters(ProblemSpecP& cm_ps)
{
  ProblemSpecP ps = cm_ps->findBlock("model_parameters");
  if (!ps) {
    std::ostringstream out;
    out << "**Error** No  model parameters provided for Mohr-Coulomb model."
        << " Include the <model_parameters> tag in the input file.\n";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  d_params.G = 10000.0;
  d_params.K = 20000.0;
  ps->require("shear_modulus", d_params.G);
  ps->require("bulk_modulus", d_params.K);

  d_params.c    = 0.0;
  d_params.phi  = 30.0;
  d_params.psi  = 30.0;
  d_params.pMin = -1;
  ps->require("cohesion", d_params.c);
  ps->require("angle_internal_friction", d_params.phi);
  ps->require("angle_dilation", d_params.psi);
  ps->getWithDefault("max_hydrostatic_tension", d_params.pMin, 0.0);

  d_params.initialSuction = 0.0;
  d_params.phi_b          = 0.0;
  ps->getWithDefault("initial_suction", d_params.initialSuction, 0.0);
  ps->getWithDefault("phi_b", d_params.phi_b, 0.0);

  // whether water retention curve should be used
  d_flags.useWaterRetention = false;
  ps->getWithDefault("use_water_retention", d_flags.useWaterRetention, false);

  ps->getWithDefault(
    "water_retention_param_1", d_params.waterRetentionParams[0], 0.0);
  ps->getWithDefault(
    "water_retention_param_2", d_params.waterRetentionParams[1], 0.0);
  ps->getWithDefault(
    "water_retention_param_3", d_params.waterRetentionParams[2], 0.0);
  ps->getWithDefault(
    "water_retention_param_4", d_params.waterRetentionParams[3], 0.0);

  // whether undrained shear strength transition should be used
  d_flags.useUndrainedShearTransition = false;
  ps->getWithDefault("use_undrained_shear_transition",
                     d_flags.useUndrainedShearTransition,
                     false);

  // water influence parameters
  d_params.waterInfluenceA1 = 0.0;
  d_params.waterInfluenceB1 = 0.0;
  d_params.waterInfluenceW  = 0.0;
  ps->getWithDefault("water_influence_A1",
                     d_params.waterInfluenceA1,
                     0.0); // water influence parameter
  ps->getWithDefault("water_influence_B1",
                     d_params.waterInfluenceB1,
                     0.0); // water influence parameter
  ps->getWithDefault("water_influence_W",
                     d_params.waterInfluenceW,
                     0.0); // water content

  // strain rate influence parameters
  d_params.betaStrainRate  = 0.0;
  d_params.refStrainRate   = 0.0;
  ps->getWithDefault("beta_strain_rate", d_params.betaStrainRate, 0.0);
  ps->getWithDefault("ref_strain_rate", d_params.refStrainRate, 0.0);

  // use variable elastic modulus
  d_flags.useVariableElasticModulus = false;
  ps->getWithDefault(
    "use_variable_elastic_modulus", d_flags.useVariableElasticModulus, false);

  d_params.variableModulusM           = 0.0;
  d_params.variableModulusNuY         = 0.0;
  ps->getWithDefault("variable_modulus_m", d_params.variableModulusM, 0.0);
  ps->getWithDefault("variable_modulus_nu_y", d_params.variableModulusNuY, 0.0);

  // use linearly varying cohesion with depth
  d_flags.useLinearlyVaryingCohesion = false;
  ps->getWithDefault(
    "use_linearly_varying_cohesion", d_flags.useLinearlyVaryingCohesion, false);

  d_params.linearCohesionA = 0.0, d_params.linearCohesionYRef = 0.0;
  d_params.depthDirection = "y-"; // "x+", "x-", "y+", "y-", "z+", "z-"
  ps->getWithDefault("linear_cohesion_a", d_params.linearCohesionA, 0.0);
  ps->getWithDefault("linear_cohesion_y_ref", d_params.linearCohesionYRef, 0.0);
  ps->getWithDefault(
    "linear_cohesion_depth_direction", d_params.depthDirection, "y-");

  // use softening model
  d_flags.useSoftening = false;
  ps->getWithDefault("use_softening", d_flags.useSoftening, false);

  d_params.softeningSt = 0.0, d_params.softeningStrain95 = 0.0;
  ps->getWithDefault("softening_St", d_params.softeningSt, 0.0);
  ps->getWithDefault("softening_strain_95", d_params.softeningStrain95, 0.0);

  // use regularized nonlocal softening
  d_flags.useRegularizedNonlocalSoftening = false;
  ps->getWithDefault("use_regularized_nonlocal_softening",
                     d_flags.useRegularizedNonlocalSoftening,
                     false);

  d_params.regularizationTFE = 0.0, d_params.regularizationTShear = 0.0;
  ps->getWithDefault("regularization_t_FE", d_params.regularizationTFE, 0.0);
  ps->getWithDefault(
    "regularization_t_shear", d_params.regularizationTShear, 0.0);

  // use nonlocal correction
  d_flags.useNonlocalCorrection = false;
  ps->getWithDefault(
    "use_nonlocal_correction", d_flags.useNonlocalCorrection, false);

  d_params.nonlocalN = 0.0, d_params.nonlocalL = 0.0;
  ps->getWithDefault("nonlocal_n", d_params.nonlocalN, 0.0);
  ps->getWithDefault("nonlocal_l", d_params.nonlocalL, 0.0);

  // retention model
  std::string retentionModel =
    "gallipoli"; // "state_surface", "van_genuchten", "gallipoli"
  ps->getWithDefault("retention_model", retentionModel, "gallipoli");
  if (retentionModel == "state_surface") {
    d_params.retentionModel = RetentionModel::STATE_SURFACE;
  } else if (retentionModel == "van_genuchten") {
    d_params.retentionModel = RetentionModel::VAN_GENUCHTEN;
  } else {
    d_params.retentionModel = RetentionModel::GALLIPOLI;
  }

  // Check that model parameters are valid
  checkModelParameters();
}

void
MohrCoulomb::checkModelParameters() const
{
  if (d_params.G <= 0.0) {
    std::ostringstream err;
    err << "Shear Modulus in the Mohr-Coulomb Model equal to: " << d_params.G
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  if (d_params.c < 0.0) {
    std::ostringstream err;
    err << "Cohesion in the Mohr-Coulomb Model is set to negative value: "
        << d_params.c << " This will cause the code to malfuncion. Any results "
                         "obtained are invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  if (d_params.phi < 0.0 || d_params.phi > 90.0) {
    std::ostringstream err;
    err << "Friction angle in the Mohr-Coulomb Model is set to: "
        << d_params.phi
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
  if (d_params.psi < 0.0 || d_params.psi > 90.0) {
    std::ostringstream err;
    err << "Dilation angle in the Mohr-Coulomb Model is set to: "
        << d_params.psi
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }
}

void
MohrCoulomb::getIntegrationParameters(ProblemSpecP& cm_ps)
{
  dbg_doing << "Doing MohrCoulomb::getIntegrationParameters\n";

  ProblemSpecP ps = cm_ps->findBlock("integration_parameters");
  if (!ps) {
    std::ostringstream out;
    out << "**Error** No integration parameters provided for Mohr-Coulomb "
           "model."
        << " Include the <integration_parameters> tag in the input file.\n";
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  d_int.maxIter    = 200;
  d_int.alfaCheck  = 1;
  d_int.alfaChange = 0.05;
  d_int.alfaRatio  = 10;
  ps->getWithDefault("max_iterations_pegasus", d_int.maxIter, 200);
  ps->getWithDefault("alpha_check_pegasus", d_int.alfaCheck, 1);
  ps->getWithDefault("alpha_change_pegasus", d_int.alfaChange, 0.05);
  ps->getWithDefault("alpha_ratio_pegasus", d_int.alfaRatio, 10);

  d_int.yieldTol       = 1.0e-6;
  d_int.integrationTol = 0.01;
  d_int.betaFactor     = 0.9;
  d_int.minMeanStress  = -1.0e8;
  d_int.suctionTol     = 1.0e-8;
  ps->getWithDefault("yield_tolerance", d_int.yieldTol, 1.0e-6);
  ps->getWithDefault("integration_tolerance", d_int.integrationTol, 0.01);
  ps->getWithDefault("beta_safety_factor", d_int.betaFactor, 0.9);
  ps->getWithDefault("minimum_mean_stress", d_int.minMeanStress, -1.0e8);
  ps->getWithDefault("suction_tolerance", d_int.suctionTol, 1.0e-8);

  d_int.driftCorrection = "at_end"; // "none", "at_begin", "at_end"
  d_int.tolMethod       = "sloan";  // "relative", "sloan";
  d_int.solutionAlgorithm =
    "modified_Euler"; // "RK3", "RK3_Bogacki", "RK4", "RK5_England",
                      // "RK5_Cash", "RK5_Dormand", "RK5_Bogacki",
                      // "extrapolation"
  ps->getWithDefault(
    "drift_correction_algorithm", d_int.driftCorrection, "at_end");
  ps->getWithDefault("tolerance_algorithm", d_int.tolMethod, "sloan");
  ps->getWithDefault(
    "solution_algorithm", d_int.solutionAlgorithm, "modified_euler");
}

/**
 * Set parameters
 */
void
MohrCoulomb::setIntegrationParameters()
{
  dbg_doing << "Doing MohrCoulomb::setIntegrationParameters\n";

  DriftCorrection drift;
  if (d_int.driftCorrection == "none") {
    drift = DriftCorrection::NO_CORRECTION;
  } else if (d_int.driftCorrection == "at_begin") {
    drift = DriftCorrection::CORRECTION_AT_BEGIN;
  } else {
    drift = DriftCorrection::CORRECTION_AT_BEGIN;
  }

  ToleranceMethod tol;
  if (d_int.tolMethod == "relative") {
    tol = ToleranceMethod::EPUS_RELATIVE_ERROR;
  } else {
    tol = ToleranceMethod::SLOAN;
  }

  SolutionAlgorithm sol;
  if (d_int.solutionAlgorithm == "modified_Euler") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER;
  } else if (d_int.solutionAlgorithm == "RK3") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_NYSTROM;
  } else if (d_int.solutionAlgorithm == "RK3_Bogacki") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_BOGACKI;
  } else if (d_int.solutionAlgorithm == "RK4") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FOURTH_ORDER;
  } else if (d_int.solutionAlgorithm == "RK5_England") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_ENGLAND;
  } else if (d_int.solutionAlgorithm == "RK5_Cash") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_CASH;
  } else if (d_int.solutionAlgorithm == "RK5_Dormand") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_DORMAND;
  } else if (d_int.solutionAlgorithm == "RK5_Bogacki") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_BOGACKI;
  } else if (d_int.solutionAlgorithm == "extrapolation") {
    sol = SolutionAlgorithm::EXTRAPOLATION_BULIRSCH;
  } else {
    sol = SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER;
  }

  d_modelP->setIntegrationParameters(d_int.maxIter,
                                     d_int.alfaCheck,
                                     d_int.alfaChange,
                                     d_int.alfaRatio,
                                     d_int.yieldTol,
                                     d_int.integrationTol,
                                     d_int.betaFactor,
                                     d_int.minMeanStress,
                                     d_int.suctionTol,
                                     sol,
                                     tol,
                                     drift);
}

/**
 * Set up internal state variables
 */
void
MohrCoulomb::initializeLocalMPMLabels()
{
  dbg_doing << "Doing MohrCoulomb::initializeLocalMPMLabels\n";

  pStrainLabel = VarLabel::create(
    "p.strainMC", ParticleVariable<Matrix3>::getTypeDescription());
  pStrainLabel_preReloc = VarLabel::create(
    "p.strainMC+", ParticleVariable<Matrix3>::getTypeDescription());

  pPlasticStrainLabel = VarLabel::create(
    "p.plasticstrainMC", ParticleVariable<Matrix3>::getTypeDescription());
  pPlasticStrainLabel_preReloc = VarLabel::create(
    "p.plasticstrainMC+", ParticleVariable<Matrix3>::getTypeDescription());

  pShearModulusLabel = VarLabel::create(
    "p.shearModulusMC", ParticleVariable<double>::getTypeDescription());
  pShearModulusLabel_preReloc = VarLabel::create(
    "p.shearModulusMC+", ParticleVariable<double>::getTypeDescription());

  pBulkModulusLabel = VarLabel::create(
    "p.bulkModulusMC", ParticleVariable<double>::getTypeDescription());
  pBulkModulusLabel_preReloc = VarLabel::create(
    "p.bulkModulusMC+", ParticleVariable<double>::getTypeDescription());

  pCohesionLabel = VarLabel::create(
    "p.cohesionMC", ParticleVariable<double>::getTypeDescription());
  pCohesionLabel_preReloc = VarLabel::create(
    "p.cohesionMC+", ParticleVariable<double>::getTypeDescription());

  pSuctionLabel = VarLabel::create(
    "p.suctionMC", ParticleVariable<double>::getTypeDescription());
  pSuctionLabel_preReloc = VarLabel::create(
    "p.suctionMC+", ParticleVariable<double>::getTypeDescription());

  pSpVolLabel = VarLabel::create(
    "p.spvolMC", ParticleVariable<double>::getTypeDescription());
  pSpVolLabel_preReloc = VarLabel::create(
    "p.spvolMC+", ParticleVariable<double>::getTypeDescription());

  pShearStrainLabel = VarLabel::create(
    "p.shearstrainMC", ParticleVariable<double>::getTypeDescription());
  pShearStrainLabel_preReloc = VarLabel::create(
    "p.shearstrainMC+", ParticleVariable<double>::getTypeDescription());

  pShearStrainRateLabel = VarLabel::create(
    "p.shearstrainrateMC", ParticleVariable<double>::getTypeDescription());
  pShearStrainRateLabel_preReloc = VarLabel::create(
    "p.shearstrainrateMC+", ParticleVariable<double>::getTypeDescription());
}

/**
 * Destructor
 */
MohrCoulomb::~MohrCoulomb()
{
  dbg_doing << "Doing MohrCoulomb::destruct\n";

  VarLabel::destroy(pStrainLabel);
  VarLabel::destroy(pStrainLabel_preReloc);

  VarLabel::destroy(pPlasticStrainLabel);
  VarLabel::destroy(pPlasticStrainLabel_preReloc);

  VarLabel::destroy(pShearModulusLabel);
  VarLabel::destroy(pShearModulusLabel_preReloc);

  VarLabel::destroy(pBulkModulusLabel);
  VarLabel::destroy(pBulkModulusLabel_preReloc);

  VarLabel::destroy(pCohesionLabel);
  VarLabel::destroy(pCohesionLabel_preReloc);

  VarLabel::destroy(pSuctionLabel);
  VarLabel::destroy(pSuctionLabel_preReloc);

  VarLabel::destroy(pSpVolLabel);
  VarLabel::destroy(pSpVolLabel_preReloc);

  VarLabel::destroy(pShearStrainLabel);
  VarLabel::destroy(pShearStrainLabel_preReloc);

  VarLabel::destroy(pShearStrainRateLabel);
  VarLabel::destroy(pShearStrainRateLabel_preReloc);
}

/**
 * For restart files
 */
void
MohrCoulomb::outputProblemSpec(ProblemSpecP& ps, bool output_cm_tag)
{
  dbg_doing << "Doing MohrCoulomb::outputProblemSpec\n";

  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "mohr_coulomb");
  }

  cm_ps->appendElement("model_type", d_modelType);

  outputModelProblemSpec(ps);
  outputIntegrationProblemSpec(ps);
}

void
MohrCoulomb::outputModelProblemSpec(ProblemSpecP& ps) const
{
  dbg_doing << "Doing MohrCoulomb::outputModelProblemSpec\n";

  ProblemSpecP cm_ps = ps->appendChild("model_parameters");

  cm_ps->appendElement("shear_modulus", d_params.G);
  cm_ps->appendElement("bulk_modulus", d_params.K);

  cm_ps->appendElement("cohesion", d_params.c);
  cm_ps->appendElement("angle_internal_friction", d_params.phi);
  cm_ps->appendElement("angle_dilation", d_params.psi);
  cm_ps->appendElement("max_hydrostatic_tension", d_params.pMin);

  cm_ps->appendElement("initial_suction", d_params.initialSuction);
  cm_ps->appendElement("phi_b", d_params.phi_b);

  cm_ps->appendElement("use_water_retention", d_flags.useWaterRetention);

  cm_ps->appendElement("water_retention_param_1",
                       d_params.waterRetentionParams[0]);
  cm_ps->appendElement("water_retention_param_2",
                       d_params.waterRetentionParams[1]);
  cm_ps->appendElement("water_retention_param_3",
                       d_params.waterRetentionParams[2]);
  cm_ps->appendElement("water_retention_param_4",
                       d_params.waterRetentionParams[3]);

  cm_ps->appendElement("use_undrained_shear_transition",
                       d_flags.useUndrainedShearTransition);

  cm_ps->appendElement("water_influence_A1", d_params.waterInfluenceA1);
  cm_ps->appendElement("water_influence_B1", d_params.waterInfluenceB1);
  cm_ps->appendElement("water_influence_W", d_params.waterInfluenceW);

  cm_ps->appendElement("beta_strain_rate", d_params.betaStrainRate);
  cm_ps->appendElement("ref_strain_rate", d_params.refStrainRate);

  cm_ps->appendElement("use_variable_elastic_modulus",
                       d_flags.useVariableElasticModulus);

  cm_ps->appendElement("variable_modulus_m", d_params.variableModulusM);
  cm_ps->appendElement("variable_modulus_nu_y", d_params.variableModulusNuY);

  cm_ps->appendElement("use_linearly_varying_cohesion",
                       d_flags.useLinearlyVaryingCohesion);

  cm_ps->appendElement("linear_cohesion_a", d_params.linearCohesionA);
  cm_ps->appendElement("linear_cohesion_y_ref", d_params.linearCohesionYRef);
  cm_ps->appendElement("linear_cohesion_depth_direction",
                       d_params.depthDirection);

  cm_ps->appendElement("use_softening", d_flags.useSoftening);

  cm_ps->appendElement("softening_St", d_params.softeningSt);
  cm_ps->appendElement("softening_strain_95", d_params.softeningStrain95);

  cm_ps->appendElement("use_regularized_nonlocal_softening",
                       d_flags.useRegularizedNonlocalSoftening);

  cm_ps->appendElement("regularization_t_FE", d_params.regularizationTFE);
  cm_ps->appendElement("regularization_t_shear", d_params.regularizationTShear);

  cm_ps->appendElement("use_nonlocal_correction",
                       d_flags.useNonlocalCorrection);

  cm_ps->appendElement("nonlocal_n", d_params.nonlocalN);
  cm_ps->appendElement("nonlocal_l", d_params.nonlocalL);

  switch (d_params.retentionModel) {
    case RetentionModel::STATE_SURFACE:
      cm_ps->appendElement("retention_model", "state_surface");
      break;
    case RetentionModel::VAN_GENUCHTEN:
      cm_ps->appendElement("retention_model", "van_genuchten");
      break;
    case RetentionModel::GALLIPOLI:
      cm_ps->appendElement("retention_model", "gallipoli");
      break;
    default:
      break;
  }
}

void
MohrCoulomb::outputIntegrationProblemSpec(ProblemSpecP& ps) const
{
  dbg_doing << "Doing MohrCoulomb::outputIntegrationSpec\n";

  ProblemSpecP cm_ps = ps->appendChild("integration_parameters");

  cm_ps->appendElement("max_iterations_pegasus", d_int.maxIter);
  cm_ps->appendElement("alpha_check_pegasus", d_int.alfaCheck);
  cm_ps->appendElement("alpha_change_pegasus", d_int.alfaChange);
  cm_ps->appendElement("alpha_ratio_pegasus", d_int.alfaRatio);

  cm_ps->appendElement("yield_tolerance", d_int.yieldTol);
  cm_ps->appendElement("integration_tolerance", d_int.integrationTol);
  cm_ps->appendElement("beta_safety_factor", d_int.betaFactor);
  cm_ps->appendElement("minimum_mean_stress", d_int.minMeanStress);
  cm_ps->appendElement("suction_tolerance", d_int.suctionTol);

  cm_ps->appendElement("drift_correction_algorithm", d_int.driftCorrection);
  cm_ps->appendElement("tolerance_algorithm", d_int.tolMethod);
  cm_ps->appendElement("solution_algorithm", d_int.solutionAlgorithm);
}

/**
 * Add the internal state variables to the particle state
 * for relocation
 */
void
MohrCoulomb::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  dbg_doing << "Doing MohrCoulomb::addParticleState\n";

  // This is an INCREMENTAL model. Needs polar decomp R to be saved.
  from.push_back(lb->pPolarDecompRLabel);
  to.push_back(lb->pPolarDecompRLabel_preReloc);

  from.push_back(pStrainLabel);
  to.push_back(pStrainLabel_preReloc);

  from.push_back(pPlasticStrainLabel);
  to.push_back(pPlasticStrainLabel_preReloc);

  from.push_back(pShearModulusLabel);
  to.push_back(pShearModulusLabel_preReloc);

  from.push_back(pBulkModulusLabel);
  to.push_back(pBulkModulusLabel_preReloc);

  from.push_back(pCohesionLabel);
  to.push_back(pCohesionLabel_preReloc);

  from.push_back(pSuctionLabel);
  to.push_back(pSuctionLabel_preReloc);

  from.push_back(pSpVolLabel);
  to.push_back(pSpVolLabel_preReloc);

  from.push_back(pShearStrainLabel);
  to.push_back(pShearStrainLabel_preReloc);

  from.push_back(pShearStrainRateLabel);
  to.push_back(pShearStrainRateLabel_preReloc);
}

/**
 * Initialization of explicit time step
 */
void
MohrCoulomb::addInitialComputesAndRequires(Task* task,
                                           const MPMMaterial* matl,
                                           const PatchSet*) const
{
  dbg_doing << "Doing MohrCoulomb::addInitialComputesAndRequires\n";

  const MaterialSubset* matlset = matl->thisMaterial();

  task->computes(pStrainLabel, matlset);
  task->computes(pPlasticStrainLabel, matlset);
  task->computes(pShearModulusLabel, matlset);
  task->computes(pBulkModulusLabel, matlset);
  task->computes(pCohesionLabel, matlset);
  task->computes(pSuctionLabel, matlset);
  task->computes(pSpVolLabel, matlset);
  task->computes(pShearStrainLabel, matlset);
  task->computes(pShearStrainRateLabel, matlset);
}

void
MohrCoulomb::initializeCMData(const Patch* patch,
                              const MPMMaterial* matl,
                              DataWarehouse* new_dw)
{
  dbg_doing << "Doing MohrCoulomb::initializeCMData\n";

  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  constParticleVariable<double> pMass, pVolume;
  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);

  ParticleVariable<Matrix3> pStrain, pPlasticStrain;
  new_dw->allocateAndPut(pStrain, pStrainLabel, pset);
  new_dw->allocateAndPut(pPlasticStrain, pPlasticStrainLabel, pset);

  ParticleVariable<double> pShearModulus, pBulkModulus, pCohesion, pSuction,
    pSpVol, pShearStrain, pShearStrainRate;
  new_dw->allocateAndPut(pShearModulus, pShearModulusLabel, pset);
  new_dw->allocateAndPut(pBulkModulus, pBulkModulusLabel, pset);
  new_dw->allocateAndPut(pCohesion, pCohesionLabel, pset);
  new_dw->allocateAndPut(pSuction, pSuctionLabel, pset);
  new_dw->allocateAndPut(pSpVol, pSpVolLabel, pset);
  new_dw->allocateAndPut(pShearStrain, pShearStrainLabel, pset);
  new_dw->allocateAndPut(pShearStrainRate, pShearStrainRateLabel, pset);

  Matrix3 Zero(0.0);
  for (auto idx : *pset) {
    pStrain[idx]          = Zero;
    pPlasticStrain[idx]   = Zero;
    pShearModulus[idx]    = d_params.G;
    pBulkModulus[idx]     = d_params.K;
    pCohesion[idx]        = d_params.c;
    pSuction[idx]         = 0.0;
    pSpVol[idx]           = pVolume[idx] / pMass[idx];
    pShearStrain[idx]     = 0.0;
    pShearStrainRate[idx] = 0.0;
  }

  computeStableTimestep(patch, matl, new_dw);
}

/**
 * This is only called for the initial timestep - all other timesteps
 * are computed as a side-effect of computeStressTensor
 */
void
MohrCoulomb::computeStableTimestep(const Patch* patch,
                                   const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  dbg_doing << "Doing MohrCoulomb::computeStableTimestep\n";

  int matID            = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pVelocity, lb->pVelocityLabel, pset);

  double K = d_params.K;
  double G = d_params.G;

  // Compute bulk wave speed at each particle, store the maximum
  Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);
  for (auto idx : *pset) {

    double p_wave_speed =
      std::sqrt((K + 4.0 * G / 3.0) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + p_wave_speed;
    waveSpeed     = Max(velMax, waveSpeed);
  }

  Vector dx       = patch->dCell();
  waveSpeed       = dx / waveSpeed;
  double delT_new = waveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

/**
 * Computing the stress tensor
 */
void
MohrCoulomb::addComputesAndRequires(Task* task,
                                    const MPMMaterial* matl,
                                    const PatchSet* patches) const
{
  dbg_doing << "Doing MohrCoulomb::addComputesAndRequires\n";

  // Add the computes and requires that are common to all rotated explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);

  task->computes(lb->pdTdtLabel_preReloc, matlset);
  task->computes(lb->p_qLabel_preReloc, matlset);

  // Computes and requires for internal state data
  task->requires(Task::OldDW, pStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pPlasticStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearModulusLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pBulkModulusLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pCohesionLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pSuctionLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pSpVolLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearStrainRateLabel, matlset, Ghost::None);

  task->computes(pStrainLabel_preReloc, matlset);
  task->computes(pPlasticStrainLabel_preReloc, matlset);
  task->computes(pShearModulusLabel_preReloc, matlset);
  task->computes(pBulkModulusLabel_preReloc, matlset);
  task->computes(pCohesionLabel_preReloc, matlset);
  task->computes(pSuctionLabel_preReloc, matlset);
  task->computes(pSpVolLabel_preReloc, matlset);
  task->computes(pShearStrainLabel_preReloc, matlset);
  task->computes(pShearStrainRateLabel_preReloc, matlset);
}

void
MohrCoulomb::addComputesAndRequires(Task*,
                                    const MPMMaterial*,
                                    const PatchSet*,
                                    const bool,
                                    const bool) const
{
  std::ostringstream err;
  err << "This Mohr-Coulomb model does not support implicit time"
         "stepping and caanot me used with Immplicit MPM.\n";
  throw ProblemSetupException(err.str(), __FILE__, __LINE__);
}

void
MohrCoulomb::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  dbg_doing << "Doing MohrCoulomb::comuteStressTensor\n";

  Matrix3 Identity;
  Identity.Identity();

  double rho_orig = matl->getInitialDensity();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();

    double strainEnergy = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    int matID            = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<Point> pX;
    constParticleVariable<long64> pParticleID;
    constParticleVariable<double> pMass, pVolume, pTemperature;
    constParticleVariable<Vector> pVelocity;

    old_dw->get(pX, lb->pXLabel, pset);
    old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pStrain_old, pPlasticStrain_old;
    constParticleVariable<double> pShearModulus_old, pBulkModulus_old,
      pCohesion_old, pSuction_old, pSpVol_old, pShearStrain_old,
      pShearStrainRate_old;

    old_dw->get(pStrain_old, pStrainLabel, pset);
    old_dw->get(pPlasticStrain_old, pPlasticStrainLabel, pset);
    old_dw->get(pShearModulus_old, pShearModulusLabel, pset);
    old_dw->get(pBulkModulus_old, pBulkModulusLabel, pset);
    old_dw->get(pCohesion_old, pCohesionLabel, pset);
    old_dw->get(pSuction_old, pSuctionLabel, pset);
    old_dw->get(pSpVol_old, pSpVolLabel, pset);
    old_dw->get(pShearStrain_old, pShearStrainLabel, pset);
    old_dw->get(pShearStrainRate_old, pShearStrainRateLabel, pset);

    constParticleVariable<double> pVolume_new;
    constParticleVariable<Matrix3> pDefRate_mid, pDefGrad_new, pStress_old;

    new_dw->get(pVolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pDefRate_mid, lb->pDeformRateMidLabel, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;

    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel_preReloc, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStrain_new, pPlasticStrain_new;
    ParticleVariable<double> pShearModulus_new, pBulkModulus_new, pCohesion_new,
      pSuction_new, pSpVol_new, pShearStrain_new, pShearStrainRate_new;

    new_dw->allocateAndPut(pStrain_new, pStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pPlasticStrain_new, pPlasticStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pShearModulus_new, pShearModulusLabel_preReloc, pset);
    new_dw->allocateAndPut(pBulkModulus_new, pBulkModulusLabel_preReloc, pset);
    new_dw->allocateAndPut(pCohesion_new, pCohesionLabel_preReloc, pset);
    new_dw->allocateAndPut(pSuction_new, pSuctionLabel_preReloc, pset);
    new_dw->allocateAndPut(pSpVol_new, pSpVolLabel_preReloc, pset);
    new_dw->allocateAndPut(pShearStrain_new, pShearStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(
      pShearStrainRate_new, pShearStrainRateLabel_preReloc, pset);

    for (auto idx : *pset) {

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Shear strain and shear strain rate
      Vector6 D_vec      = Uintah::toVector6(pDefRate_mid[idx]);
      Vector6 strain_vec = Uintah::toVector6(pStrain_old[idx]);
      strain_vec += D_vec * delT;

      double shearStrain     = computeShearStrain(strain_vec);
      double shearStrainRate = computeShearStrain(D_vec);

      // Non local correction
      bool nonlocal = d_flags.useNonlocalCorrection;

      ParticleSubset* pset_neighbor = nullptr;

      constParticleVariable<Point> pX_neighbor;
      constParticleVariable<double> pVolume_neighbor;

      if (nonlocal) {

        Vector cellSize        = patch->dCell();
        double nonlocal_length = d_params.nonlocalL;
        Vector nonlocal_length_vec(
          nonlocal_length, nonlocal_length, nonlocal_length);
        nonlocal_length_vec += (0.5 * cellSize);

        auto cellIndexLow =
          patch->getLevel()->getCellIndex(pX[idx] - nonlocal_length_vec);
        auto cellIndexHigh =
          patch->getLevel()->getCellIndex(pX[idx] + nonlocal_length_vec);

        pset_neighbor =
          old_dw->getParticleSubset(matID, patch, cellIndexLow, cellIndexHigh);

        old_dw->get(pX_neighbor, lb->pXLabel, pset_neighbor);
        new_dw->get(pVolume_neighbor, lb->pVolumeLabel_preReloc, pset_neighbor);

        std::tie(shearStrain, shearStrainRate) =
          computeNonlocalShearStrains(pX[idx],
                                      pShearModulus_old[idx],
                                      pCohesion_old[idx],
                                      shearStrain,
                                      shearStrainRate,
                                      pset_neighbor,
                                      pX_neighbor,
                                      pVolume_neighbor);
      }

      pShearStrain_new[idx]     = shearStrain;
      pShearStrainRate_new[idx] = shearStrainRate;

      // Update internal variables and do integration
      Vector6 strainInc         = D_vec * delT;
      Vector6 stress_vec        = Uintah::toVector6(pStress_old[idx]);
      Vector6 plasticStrain_vec = Uintah::toVector6(pPlasticStrain_old[idx]);
      double cohesion           = pCohesion_old[idx];
      double shearModulus       = pShearModulus_old[idx];
      double bulkModulus        = pBulkModulus_old[idx];
      double suction            = pSuction_old[idx];
      double specificVol        = pSpVol_old[idx];

      dbg << __LINE__ << "MohrCoulomb::Particle ID = " << pParticleID[idx] << "\n";

      calculateStress(pParticleID[idx],
                      pX[idx],
                      strainInc,
                      shearStrain,
                      shearStrainRate,
                      stress_vec,
                      strain_vec,
                      plasticStrain_vec,
                      cohesion,
                      shearModulus,
                      bulkModulus,
                      suction,
                      specificVol);

      pStress_new[idx]        = Uintah::toMatrix3(stress_vec);
      pStrain_new[idx]        = Uintah::toMatrix3(strain_vec);
      pPlasticStrain_new[idx] = Uintah::toMatrix3(plasticStrain_vec);
      pCohesion_new[idx]      = cohesion;
      pShearModulus_new[idx]  = shearModulus;
      pBulkModulus_new[idx]   = bulkModulus;
      pSuction_new[idx]       = suction;
      pSpVol_new[idx]         = specificVol;

      // Nonlocal correction to stress
      if (nonlocal && pset_neighbor) {
        Vector6 stress_nonlocal = computeNonlocalStress(
          pX[idx], stress_vec, pset_neighbor, pX_neighbor, pVolume_neighbor);
        pStress_new[idx] = Uintah::toMatrix3(stress_nonlocal);
      }

      // Compute the local p-wave speed
      double J       = pDefGrad_new[idx].Determinant();
      double rho_cur = rho_orig / J;
      double p_wave_speed =
        std::sqrt((bulkModulus + 4.0 / 3.0 * shearModulus) / rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 avgStress = (pStress_new[idx] + pStress_old[idx]) * 0.5;

      double rateOfWork = computeRateOfWork(avgStress, pDefRate_mid[idx]);
      strainEnergy += (rateOfWork * pVolume_new[idx] * delT);

      // Compute wave speed at each particle, store the maximum
      Vector velMax = pVelocity[idx].cwiseAbs() + p_wave_speed;
      waveSpeed     = Max(velMax, waveSpeed);

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave          = (dx.x() + dx.y() + dx.z()) / 3.0;
        double bulk_wave_speed = std::sqrt(bulkModulus / rho_cur);
        p_q[idx]               = artificialBulkViscosity(
          pDefRate_mid[idx].Trace(), bulk_wave_speed, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.0;
      }

    } // end loop over particles

    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(strainEnergy), lb->StrainEnergyLabel);
    }
  }
}

double
MohrCoulomb::computeShearStrain(const Vector6& strain) const
{
  dbg_doing << "Doing MohrCoulomb::computeShearStrain\n";

  double diff12   = strain(0) - strain(1);
  double diff13   = strain(0) - strain(2);
  double diff23   = strain(1) - strain(2);
  double strain12 = strain(3);
  double strain13 = strain(4);
  double strain23 = strain(5);

  double shearStrain =
    1.0 / 2.0 *
    std::sqrt(
      2.0 * (diff12 * diff12 + diff13 * diff13 + diff23 * diff23) +
      3.0 * (strain12 * strain12 + strain13 * strain13 + strain23 * strain23));
  return shearStrain;
}

std::tuple<double, double>
MohrCoulomb::computeNonlocalShearStrains(
  const Point& pX,
  double pShearModulus,
  double pCohesion,
  double shearStrain,
  double shearStrainRate,
  ParticleSubset* pset_neighbor,
  constParticleVariable<Point>& pX_neighbor,
  constParticleVariable<double>& pVolume_neighbor)
{
  dbg_doing << "Doing MohrCoulomb::computeNonlocalSHearStrains\n";

  double nonlocal_length = d_params.nonlocalL;
  double domain_nonlocal = nonlocal_length * nonlocal_length;
  double up              = 0;
  double up1             = 0;
  double down            = 0;

  for (auto idx : *pset_neighbor) {

    double distSq = (pX_neighbor[idx] - pX).length2();
    double dist   = std::sqrt(distSq);
    double weight =
      dist * std::exp(-distSq / domain_nonlocal) / domain_nonlocal;
    double weightedVolume = weight * pVolume_neighbor[idx];
    up += (shearStrain * weightedVolume);
    up1 += (shearStrainRate * weightedVolume);
    down += weightedVolume;
  }

  double n_nonlocal = d_params.nonlocalN;
  double shearStrain_nonlocal =
    (1 - n_nonlocal) * shearStrain + (n_nonlocal * up / down);
  double shearStrainRate_nonlocal =
    (1 - n_nonlocal) * shearStrainRate + (n_nonlocal * up1 / down);

  double elastic_strain = pCohesion / pShearModulus;
  bool regularize       = d_flags.useRegularizedNonlocalSoftening;
  if (regularize) {
    double tFE    = d_params.regularizationTFE;
    double tShear = d_params.regularizationTShear;
    if (shearStrain_nonlocal > elastic_strain) {
      shearStrain_nonlocal     = shearStrain_nonlocal * tFE / tShear;
      shearStrainRate_nonlocal = shearStrainRate_nonlocal * tFE / tShear;
    }
  }

  return std::make_tuple(shearStrain_nonlocal, shearStrainRate_nonlocal);
}

Vector6
MohrCoulomb::computeNonlocalStress(
  const Point& pX,
  const Vector6& pStress_vec,
  ParticleSubset* pset_neighbor,
  constParticleVariable<Point>& pX_neighbor,
  constParticleVariable<double>& pVolume_neighbor)
{
  dbg_doing << "Doing MohrCoulomb::computeNonlocalStress\n";

  double nonlocal_length = d_params.nonlocalL;
  double domain_nonlocal = nonlocal_length * nonlocal_length;

  Vector6 up_vec = Vector6::Zero();
  double down    = 0;

  for (auto idx : *pset_neighbor) {

    double distSq = (pX_neighbor[idx] - pX).length2();
    double dist   = std::sqrt(distSq);
    double weight =
      dist * std::exp(-distSq / domain_nonlocal) / domain_nonlocal;
    double weightedVolume = weight * pVolume_neighbor[idx];

    up_vec += (pStress_vec * weightedVolume);
    down += weightedVolume;
  }

  double n_nonlocal = d_params.nonlocalN;
  Vector6 stress_nonlocal =
    pStress_vec * (1.0 - n_nonlocal) + up_vec * (n_nonlocal / down);

  return stress_nonlocal;
}

/**
 ***********************************************************************
 * The procedure calculates stress for the Mohr-Coulomb model, with several
 * additional options. Those include the variation of the flow rule,
 * several choices of Mohr-Coulomb like yield surfaces and the dependency
 * on the strain rate. Most of the above-described features are 'work in
 * progress'
 *
 *    input arguments
 *    ===============
 *     NBLK       int                   Number of blocks to be processed
 *     NINSV      int                   Number of internal state vars
 *     DTARG      dp                    Current time increment
 *     UI       dp,ar(nprop)            User inputs
 *     D          dp,ar(6)              Strain increment
 *
 *    input output arguments
 *    ======================
 *     STRESS   dp,ar(6)                stress
 *     SVARG    dp,ar(ninsv)            state variables
 *
 *    output arguments
 *    ================
 *     USM      dp                      uniaxial strain modulus
 *
 ***********************************************************************
 *
 *     stresss and strains, plastic strain tensors
 *         11, 22, 33, 12, 23, 13
 *
 ***********************************************************************
 * Flavour
 * 1- classic Mohr - Coulomb,
 * 2 - classic Mohr-Coulomb with tension cut-off (not implemented yet)
 * 3 - ShengMohrCoulomb (rounded MC surface see eq 13 in Sheng D, Sloan SW & Yu
 * HS
 * Computations Mechanics 26:185-196 (2000)
 * 4 - ShengMohrCoulomb with tension cut-off (not implemented yet)
 * 5 - SloanMohrCoulomb (see Sloan et al... , not implemented yet)
 * 6 - SloanMohrCoulomb with tension cut-off
 * 11- classic Mohr - Coulomb with semi-implicit stress integration
 ***********************************************************************
 */
void
MohrCoulomb::calculateStress(long64 particleID,
                             const Point& pX,
                             const Vector6& dEps,
                             double pShearStrain,
                             double pShearStrainRate,
                             Vector6& stress,
                             Vector6& strain,
                             Vector6& plasticStrain,
                             double& pCohesion,
                             double& pShearModulus,
                             double& pBulkModulus,
                             double& pSuction,
                             double& pSpecificVol)
{
  dbg_doing << "Doing MohrCoulomb::calculateStress\n";

  // Convert to compression +ve form
  Vector6 strainInc = -dEps;
  MohrCoulombState initialState;
  initialState.particleID = particleID;
  initialState.stress = -stress;
  initialState.strain.block<6, 1>(0, 0) = -strain + strainInc;
  initialState.plasticStrain = -plasticStrain;

  // Rate dependence
  if (d_flags.useUndrainedShearTransition) {

    double St = d_params.softeningSt;
    double a1 = d_params.waterInfluenceA1;
    double b1 = d_params.waterInfluenceB1;
    double W  = d_params.waterInfluenceW;

    pCohesion = St * a1 * pow(W, -b1);
    if (pShearStrainRate > d_params.refStrainRate) {
      pCohesion *= std::pow(pShearStrainRate / d_params.refStrainRate,
                            d_params.betaStrainRate);
    }
  }

  // Shear strength linear with depth
  if (d_flags.useLinearlyVaryingCohesion) {
    std::string direction = d_params.depthDirection;
    double yDiff          = 0;
    if (direction == "x-") {
      yDiff = pX.x() - d_params.linearCohesionYRef;
    } else if (direction == "x+") {
      yDiff = d_params.linearCohesionYRef - pX.x();
    } else if (direction == "y-") {
      yDiff = pX.y() - d_params.linearCohesionYRef;
    } else if (direction == "y+") {
      yDiff = d_params.linearCohesionYRef - pX.y();
    } else if (direction == "z-") {
      yDiff = pX.z() - d_params.linearCohesionYRef;
    } else if (direction == "z+") {
      yDiff = d_params.linearCohesionYRef - pX.z();
    } else {
      yDiff = pX.y() - d_params.linearCohesionYRef;
    }
    pCohesion += d_params.linearCohesionA * yDiff;
  }

  // Varying moduli
  if (d_flags.useVariableElasticModulus) {
    pShearModulus = d_params.variableModulusM * pCohesion / 2.0 /
                    (1.0 + d_params.variableModulusNuY);
    pBulkModulus = d_params.variableModulusM * pCohesion / 3.0 /
                   (1.0 - 2 * d_params.variableModulusNuY);
  }

  // Softening
  if (d_flags.useSoftening) {
    if (pShearStrain > pCohesion / pShearModulus) {
      double St = d_params.softeningSt;
      pCohesion *=
        (1.0 / St +
         (1.0 - 1.0 / St) *
           std::pow(2.71, (-3.0 * pShearStrain / d_params.softeningStrain95)));
    }
  }

  // Retention
  if (d_flags.useWaterRetention) {

    // get from current suction value
    double Sr = computeSr(pSuction);

    // constant volume of water - thus we use initial specific volume
    double volWater   = Sr * (pSpecificVol - 1.0) / pSpecificVol;
    double dVolStrain = strainInc(0) + strainInc(1) + strainInc(2);

    // volume of voids - new one, using new specific volume
    double specVol_new = pSpecificVol - dVolStrain * pSpecificVol;
    double volAir      = (specVol_new - 1) / specVol_new - volWater;

    // voids will reduce by dvol, so dSr=VolWater/TotalVolume
    double totalVolume = volAir + volWater;
    double Sr_new      = volWater / totalVolume;
    if (Sr_new > 0.99999) {
      Sr_new = 0.99999;
    } else if (Sr_new < 0.0001) {
      Sr_new = 0.0001;
    }

    double suction_new = computeSuction(Sr_new);
    pCohesion += std::tan(d_params.phi_b * M_PI / 180.0) * suction_new;

    pSuction     = suction_new;
    pSpecificVol = specVol_new;
  }

  // Do the integration
  d_modelP->setModelParameters(pShearModulus,
                               pBulkModulus,
                               pCohesion,
                               d_params.phi,
                               d_params.psi,
                               d_params.pMin);

  MohrCoulombState finalState = initialState;
  Vector7 strainIncExtra      = Vector7::Zero();
  strainIncExtra.block<6, 1>(0, 0) = strainInc;
  if (d_modelType == "classic" || "sheng") {
    finalState = d_modelP->integrate(strainIncExtra, initialState);
  } else if (d_modelType == "classic_semiimplicit") {
    RegionType region;
    finalState = d_modelP->integrate(strainIncExtra, initialState, region);
  } else {
    std::ostringstream err;
    err << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }

  // stress back for output (tension +ve)
  stress        = -finalState.stress;
  strain        = -finalState.strain.block<6, 1>(0, 0);
  plasticStrain = -finalState.plasticStrain;
}

/**
 * Compute saturation and suction
 */
double
MohrCoulomb::computeSr(double suction) const
{
  // VanGenuchten Model
  double m     = d_params.waterRetentionParams[0];
  double n     = d_params.waterRetentionParams[1];
  double alpha = d_params.waterRetentionParams[2];
  double Sr    = alpha * suction;
  Sr           = std::pow(Sr, n);
  Sr           = 1.0 / (1.0 + Sr);
  Sr           = std::pow(Sr, m);
  if (Sr > 1.0) {
    Sr = 1.0;
  } else if (Sr < 0) {
    Sr = 0.0;
  }
  return Sr;
}

double
MohrCoulomb::computeSuction(double Sr) const
{
  // VanGenuchten Model
  double m       = d_params.waterRetentionParams[0];
  double n       = d_params.waterRetentionParams[1];
  double alpha   = d_params.waterRetentionParams[2];
  double suction = std::pow(Sr, 1 / m);
  suction        = (1 - suction) / suction;
  suction        = std::pow(suction, 1 / n);
  suction        = suction / alpha;
  if (suction < 0) {
    suction = 0.0;
  }
  return suction;
}

/**
 * Carry forward model data for RigidMPM
 */
void
MohrCoulomb::carryForward(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  dbg_doing << "Doing MohrCoulomb::carryForward\n";

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch   = patches->get(p);
    int matID            = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    constParticleVariable<Matrix3> pStrain_old;
    ParticleVariable<Matrix3> pStrain_new;
    old_dw->get(pStrain_old, pStrainLabel, pset);
    new_dw->allocateAndPut(pStrain_new, pStrainLabel_preReloc, pset);
    pStrain_new.copyData(pStrain_old);

    constParticleVariable<Matrix3> pPlasticStrain_old;
    ParticleVariable<Matrix3> pPlasticStrain_new;
    old_dw->get(pStrain_old, pPlasticStrainLabel, pset);
    new_dw->allocateAndPut(pStrain_new, pPlasticStrainLabel_preReloc, pset);
    pPlasticStrain_new.copyData(pPlasticStrain_old);

    constParticleVariable<double> pShearModulus_old;
    ParticleVariable<double> pShearModulus_new;
    old_dw->get(pShearModulus_old, pShearModulusLabel, pset);
    new_dw->allocateAndPut(
      pShearModulus_new, pShearModulusLabel_preReloc, pset);
    pShearModulus_new.copyData(pShearModulus_old);

    constParticleVariable<double> pBulkModulus_old;
    ParticleVariable<double> pBulkModulus_new;
    old_dw->get(pBulkModulus_old, pBulkModulusLabel, pset);
    new_dw->allocateAndPut(pBulkModulus_new, pBulkModulusLabel_preReloc, pset);
    pBulkModulus_new.copyData(pBulkModulus_old);

    constParticleVariable<double> pCohesion_old;
    ParticleVariable<double> pCohesion_new;
    old_dw->get(pCohesion_old, pCohesionLabel, pset);
    new_dw->allocateAndPut(pCohesion_new, pCohesionLabel_preReloc, pset);
    pCohesion_new.copyData(pCohesion_old);

    constParticleVariable<double> pSuction_old;
    ParticleVariable<double> pSuction_new;
    old_dw->get(pSuction_old, pSuctionLabel, pset);
    new_dw->allocateAndPut(pSuction_new, pSuctionLabel_preReloc, pset);
    pSuction_new.copyData(pSuction_old);

    constParticleVariable<double> pSpVol_old;
    ParticleVariable<double> pSpVol_new;
    old_dw->get(pSpVol_old, pSpVolLabel, pset);
    new_dw->allocateAndPut(pSpVol_new, pSpVolLabel_preReloc, pset);
    pSpVol_new.copyData(pSpVol_old);

    constParticleVariable<double> pShearStrain_old;
    ParticleVariable<double> pShearStrain_new;
    old_dw->get(pShearStrain_old, pShearStrainLabel, pset);
    new_dw->allocateAndPut(pShearStrain_new, pShearStrainLabel_preReloc, pset);
    pShearStrain_new.copyData(pShearStrain_old);

    constParticleVariable<double> pShearStrainRate_old;
    ParticleVariable<double> pShearStrainRate_new;
    old_dw->get(pShearStrainRate_old, pShearStrainRateLabel, pset);
    new_dw->allocateAndPut(
      pShearStrainRate_new, pShearStrainRateLabel_preReloc, pset);
    pShearStrainRate_new.copyData(pShearStrainRate_old);

    // Don't affect the strain energy or timestep size
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

/**
 *  For particle material model switch mid-simulation
 */
void
MohrCoulomb::allocateCMDataAddRequires(Task* /* task */,
                                       const MPMMaterial* /* matl */,
                                       const PatchSet* /* patch */,
                                       MPMLabel* /* lb */) const
{
  std::ostringstream out;
  out << "Material conversion after failure not implemented for MohrCoulomb.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

void
MohrCoulomb::allocateCMDataAdd(DataWarehouse* /* new_dw */,
                               ParticleSubset* /* subset */,
                               LabelParticleMap* /* newState */,
                               ParticleSubset* /* delset */,
                               DataWarehouse* /* old_dw */)
{
  std::ostringstream out;
  out << "Material conversion after failure not implemented for MohrCoulomb.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
}

/**
 * Hooks for interacting with MPMICE
 */
double
MohrCoulomb::computeRhoMicroCM(double pressure,
                               const double p_ref,
                               const MPMMaterial* matl,
                               double temperature,
                               double rho_guess)
{
  proc0cout
    << "**WARNING** Bulk modulus is modified by MohrCoulomb but assumed"
    << " constant when interacting with ICE.  Results may be erroneous.\n";

  double rho_orig = matl->getInitialDensity();
  double p_gauge  = pressure - p_ref;
  double rho_cur;
  double bulk = d_params.K;

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;
}

void
MohrCoulomb::computePressEOSCM(double rho_cur,
                               double& pressure,
                               double p_ref,
                               double& dp_drho,
                               double& tmp,
                               const MPMMaterial* matl,
                               double temperature)
{
  proc0cout
    << "**WARNING** Bulk modulus is modified by MohrCoulomb but assumed"
    << " constant when interacting with ICE.  Results may be erroneous.\n";

  double bulk     = d_params.K;
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure   = p_ref + p_g;
  dp_drho    = bulk * rho_orig / (rho_cur * rho_cur);
  tmp        = bulk / rho_cur; // speed of sound squared
}

double
MohrCoulomb::getCompressibility()
{
  proc0cout
    << "**WARNING** Bulk modulus is modified by MohrCoulomb but assumed"
    << " constant when interacting with ICE.  Results may be erroneous.\n";

  return 1.0 / d_params.K;

  // return 1; //experimental:1 // found to be irrelevant
}
