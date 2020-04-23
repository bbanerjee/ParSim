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

// sets of external variables for the Sheng Mohr Coulomb algorithm by WTS. Some
// are redundant.

double INCREMENT_TYPE;
double EULER_ITERATIONS;
double STEP_MAX, STEP_MIN, ERROR_DEF, USE_ERROR_STEP, MIN_DIVISION_SIZE;
double LINES_SEC_NO = 50;
int MAX_LOOP,
  ALGORITHM_TYPE; // MAX_LOOP - number of steps to solve, SOLUTION_ALGORITHM -
                  // algorithm to use
double USE_ERROR      = 0.5,
       SAVE_FOR_ERROR = 1; // these are values - 1st - when 'use' the additional
                           // error to speed up the calculations (here 30%)
// SAVEFORERROR says how greater accuracy the program should use; 1 means no
// change. 0.5 means 50% greater accuracy. (tol * 0.5) etc
double CHGEPSINTOL  = 10e-9;
double ADDTOLYIELD  = 0.8;
int USE_NICE_SCHEME = 0;

using namespace Uintah;

static DebugStream cout_doing("MohrCoulomb", false);

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
                                                    d_params.cohesion,
                                                    d_params.phi,
                                                    d_params.psi,
                                                    d_params.pMin);
  } else if (d_modelType == "sheng") {

    d_modelP = std::make_unique<MohrCoulombSheng>(d_params.G,
                                                  d_params.K,
                                                  d_params.cohesion,
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
  d_modelType = cm->d_modelType;
  if (d_modelType == "classic" || d_modelType == "classic_semiimplicit") {

    d_modelP = std::make_unique<MohrCoulombClassic>(cm->d_modelP.get());

  } else if (d_modelType == "sheng") {

    d_modelP = std::make_unique<MohrCoulombSheng>(cm->d_modelP.get());

  } else {
    std::ostringstream err;
    err << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  d_flags  = cm->d_flags;
  d_params = cm->d_params;
  d_int    = cm->d_int;

  // Set up internal variables that can vary with each particle
  initializeLocalMPMLabels();
}

/**
 * Get parameters
 */
void
MohrCoulomb::getInputParameters(ProblemSpecP& ps)
{
  d_modelType = "classic"; // options: classic, sheng
  ps->require("mohr_coulomb_model_type", d_modelType);

  getModelParameters(ps);
  getIntegrationParameters(ps);
}

void
MohrCoulomb::getModelParameters(ProblemSpecP& cm_ps)
{
  ProblemSpecP ps = cm_ps->findBlock("model_parameters");
  if (!child) {
    std::ostringstream out;
    out << "**Error** No  model parameters provided for Mohr-Coulomb model."
        << " Include the <model_parameters> tag in the input "
           "file.\n" throw Uintah::ProblemSetupException(
             out.str(), __FILE__, __LINE__);
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

  d_params.waterRetentionParams = { 0.0, 0.0, 0.0, 0.0 };
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
  d_params.shearStrainRate = 0.0;
  ps->getWithDefault("beta_strain_rate", d_params.betaStrainRate, 0.0);
  ps->getWithDefault("ref_strain_rate", d_params.refStrainRate, 0.0);
  ps->getWithDefault("shearStrainRate", d_params.shearStrainRate, 0.0);

  // use variable elastic modulus
  d_flags.useVariableElasticModulus = false;
  ps->getWithDefault(
    "use_variable_elastic_modulus", d_flags.useVariableElasticModulus, false);

  d_params.variableModulusM           = 0.0;
  d_params.variableModulusNuY         = 0.0;
  d_params.variableModulusShearStrain = 0.0;
  ps->getWithDefault("variable_modulus_m", d_params.variableModulusM, 0.0);
  ps->getWithDefault("variable_modulus_nu_y", d_params.variableModulusNuY, 0.0);
  ps->getWithDefault(
    "variable_modulus_shear_strain", d_params.variableModulusShearStrain, 0.0);

  // use linearly varying cohesion with depth
  d_flags.useLinearlyVaryingCohesion = false;
  ps->getWithDefault(
    "use_linearly_varying_cohesion", d_flags.useLinearlyVaryingCohesion, false);

  d_params.linearCohesionA = 0.0, d_params.linearCohesionYRef = 0.0;
  d_params.depthDirection = "y-"; // "x+", "x-", "y+", "y-", "z+", "z-"
  ps->getWithDefault("linear_cohesion_a", d_params.linearCohesionA, 0.0);
  ps->getWithDefault("linear_cohesion_y_ref", d_params.linearCohesionYRef, 0.0);
  ps->getWithDefault("linear_cohesion_depth_direction", d_params.depthDirection, "y-");

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
  if (d_params.phi < 0.0 || d_params.phi > 90.0)
    std::ostringstream err;
  err << "Friction angle in the Mohr-Coulomb Model is set to: " << d_params.phi
      << " This will cause the code to malfuncion. Any results obtained are "
         "invalid.\n";
  throw ProblemSetupException(err.str(), __FILE__, __LINE__);
}
if (d_params.psi < 0.0 || d_params.psi > 90.0)
  std::ostringstream err;
err << "Dilation angle in the Mohr-Coulomb Model is set to: " << d_params.psi
    << " This will cause the code to malfuncion. Any results obtained are "
       "invalid.\n";
throw ProblemSetupException(err.str(), __FILE__, __LINE__);
}
}

void
MohrCoulomb::getIntegrationParameters(ProblemSpecP& cm_ps)
{
  ProblemSpecP ps = cm_ps->findBlock("integration_parameters");
  if (!child) {
    std::ostringstream out;
    out << "**Error** No integration parameters provided for Mohr-Coulomb "
           "model."
        << " Include the <integration_parameters> tag in the input "
           "file.\n" throw Uintah::ProblemSetupException(
             out.str(), __FILE__, __LINE__);
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
  if (d_int.solutionALgorithm == "modified_Euler") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_SECOND_ORDER_MODIFIED_EULER;
  } else if (d_int.solutionALgorithm == "RK3") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_NYSTROM;
  } else if (d_int.solutionALgorithm == "RK3_Bogacki") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_THIRD_ORDER_BOGACKI;
  } else if (d_int.solutionALgorithm == "RK4") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FOURTH_ORDER;
  } else if (d_int.solutionALgorithm == "RK5_England") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_ENGLAND;
  } else if (d_int.solutionALgorithm == "RK5_Cash") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_CASH;
  } else if (d_int.solutionALgorithm == "RK5_Dormand") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_DORMAND;
  } else if (d_int.solutionALgorithm == "RK5_Bogacki") {
    sol = SolutionAlgorithm::RUNGE_KUTTA_FIFTH_ORDER_BOGACKI;
  } else if (d_int.solutionALgorithm == "extrapolation") {
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
                                     drift,
                                     tol,
                                     sol);
}

/**
 * Set up internal state variables
 */
void
MohrCoulomb::initializeLocalMPMLabels()
{
  pStrainLabel.push_back("p.strainMC",
                         ParticleVariable<Matrix3>::getTypeDescription());
  pStrainLabel_preReloc.push_back(
    "p.strainMC+", ParticleVariable<Matrix3>::getTypeDescription());

  pShearModulusLabel.push_back("p.shearModulusMC",
                         ParticleVariable<double>::getTypeDescription());
  pShearModulusLabel_preReloc.push_back(
    "p.shearModulusMC+", ParticleVariable<double>::getTypeDescription());

  pBulkModulusLabel.push_back("p.bulkModulusMC",
                         ParticleVariable<double>::getTypeDescription());
  pBulkModulusLabel_preReloc.push_back(
    "p.bulkModulusMC+", ParticleVariable<double>::getTypeDescription());

  pCohesionLabel.push_back("p.cohesionMC",
                         ParticleVariable<double>::getTypeDescription());
  pCohesionLabel_preReloc.push_back(
    "p.cohesionMC+", ParticleVariable<double>::getTypeDescription());

  pSuctionLabel.push_back("p.suctionMC",
                         ParticleVariable<double>::getTypeDescription());
  pSuctionLabel_preReloc.push_back(
    "p.suctionMC+", ParticleVariable<double>::getTypeDescription());

  pSpVolLabel.push_back("p.spvolMC",
                         ParticleVariable<double>::getTypeDescription());
  pSpVolLabel_preReloc.push_back(
    "p.spvolMC+", ParticleVariable<double>::getTypeDescription());

  pShearStrainLabel.push_back("p.shearstrainMC",
                         ParticleVariable<double>::getTypeDescription());
  pShearStrainLabel_preReloc.push_back(
    "p.shearstrainMC+", ParticleVariable<double>::getTypeDescription());

  pShearStrainRateLabel.push_back("p.shearstrainrateMC",
                         ParticleVariable<double>::getTypeDescription());
  pShearStrainRateLabel_preReloc.push_back(
    "p.shearstrainrateMC+", ParticleVariable<double>::getTypeDescription());
}

/**
 * Destructor
 */
MohrCoulomb::~MohrCoulomb()
{
  VarLabel::destroy(pStrainLabel);
  VarLabel::destroy(pStrainLabel_preReloc);

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
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type", "mohr_coulomb");
  }

  cm_ps->appendElement("mohr_coulomb_model_type", d_modelType);

  outputModelProblemSpec(ps);
  outputIntegrationProblemSpec(ps);
}

void
MohrCoulomb::outputModelProblemSpec(ProblemSpecP& ps) const
{
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
  cm_ps->appendElement("shearStrainRate", d_params.shearStrainRate);

  cm_ps->appendElement("use_variable_elastic_modulus",
                       d_flags.useVariableElasticModulus);

  cm_ps->appendElement("variable_modulus_m", d_params.variableModulusM);
  cm_ps->appendElement("variable_modulus_nu_y", d_params.variableModulusNuY);
  cm_ps->appendElement("variable_modulus_shear_strain",
                       d_params.variableModulusShearStrain);

  cm_ps->appendElement("use_linearly_varying_cohesion",
                       d_flags.useLinearlyVaryingCohesion);

  cm_ps->appendElement("linear_cohesion_a", d_params.linearCohesionA);
  cm_ps->appendElement("linear_cohesion_y_ref", d_params.linearCohesionYRef);
  cm_ps->appendElement("linear_cohesion_depth_direction", d_params.depthDirection);

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
    case RetentionModel::VAN_GENUCHTEN;
      cm_ps->appendElement("retention_model", "van_genuchten"); break;
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

/*
 * Create clone
 */
MohrCoulomb*
MohrCoulomb::clone()
{
  return scinew MohrCoulomb(*this);
}

/**
 * Add the internal state variables to the particle state
 * for relocation
 */
void
MohrCoulomb::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  from.push_back(pStrainLabel);
  to.push_back(pStrainLabel_preReloc);

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
  cout_doing << "Doing MohrCoulomb::addInitialComputesAndRequires\n";

  task->computes(pStrainLabel, matlset);
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
  cout_doing << "Doing MohrCoulomb::initializeCMData\n";

  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  constParticleVariable<double> pMass, pVolume;
  new_dw->get(pMass, lb->pMassLabel, pset);
  new_dw->get(pVolume, lb->pVolumeLabel, pset);

  ParticleVariable<Matrix3> pStrain;
  new_dw->allocateAndPut(pStrain, pStrainLabel, pset);

  ParticleVariable<double> pShearModulus, pBulkModulus, pCohesion,
                           pSuction, pSpVol,
                           pShearStrain, pShearStrainRate;
  new_dw->allocateAndPut(pShearModulus, pShearModulusLabel, pset);
  new_dw->allocateAndPut(pBulkModulus, pBulkModulusLabel, pset);
  new_dw->allocateAndPut(pCohesion, pCohesionLabel, pset);
  new_dw->allocateAndPut(pSuction, pSuctionLabel, pset);
  new_dw->allocateAndPut(pSpVol, pSpVolLabel, pset);
  new_dw->allocateAndPut(pShearStrain, pShearStrainLabel, pset);
  new_dw->allocateAndPut(pShearStrainRate, pShearStrainRateLabel, pset);

  Matrix3 Zero(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  for (auto idx : *pset) {
    for (auto i = 0u; i < num_ISV, ++ii) {
      ISVs[i][idx] = ISVInitVals[i];
    }
    pStrain[idx] = Zero;
    pShearModulus[idx] = d_params.G;
    pBulkModulus[idx] = d_params.K;
    pCohesion[idx] = d_params.c;
    pSuction[idx] = 0.0;
    pSpVol[idx] = pVolume[idx] / pMass[idx];
    pShearStrain[idx] = 0.0;
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

    double bulkSpeedOfSound =
      std::sqrt((K + 4.0 * G / 3.0) * pVolume[idx] / pMass[idx]);
    Vector velMax = pVelocity[idx].cwiseAbs() + bulkSpeedOfSound;
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
  // Add the computes and requires that are common to all rotated explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addComputesAndRequiresForRotatedExplicit(task, matlset, patches);

  task->computes(lb->p_qLabel_preReloc, matlset);

  // Computes and requires for internal state data
  task->requires(Task::OldDW, pStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearModulusLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pBulkModulusLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pCohesionLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pSuctionLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pSpVolLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pShearStrainRateLabel, matlset, Ghost::None);

  task->computes(pStrainLabel_preReloc, matlset);
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
                                    const bool) const
{
}

void
MohrCoulomb::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl,
                                 DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{
  Matrix3 Identity;
  Identity.Identity();

  double rho_orig      = matl->getInitialDensity();

  delt_vartype delT;
  old_dw->get(delT, lb->delTLabel, getLevel(patches));

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch = patches->get(p);
    Vector dx          = patch->dCell();

    double se    = 0.0;
    double c_dil = 0.0;
    Vector waveSpeed(1.e-12, 1.e-12, 1.e-12);

    int matID            = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    constParticleVariable<Point> pX;
    constParticleVariable<double> pMass, pVolume, pVolume_new, pTemperature;
    constParticleVariable<Vector> pVelocity;

    old_dw->get(pX, lb->pXLabel, pset);
    old_dw->get(pMass, lb->pMassLabel, pset);
    old_dw->get(pVolume, lb->pVolumeLabel, pset);
    old_dw->get(pTemperature, lb->pTemperatureLabel, pset);
    old_dw->get(pVelocity, lb->pVelocityLabel, pset);

    constParticleVariable<Matrix3> pStrain_old;
    constParticleVariable<double> pShearModulus_old, pBulkModulus_old, pCohesion_old,
                                  pSuction_old, pSpVol_old, pShearStrain_old,
                                  pShearStrainRate_old;

    old_dw->get(pStrain_old, pStrainLabel, pset);
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
    new_dw->get(pDefRate_mid, lb->pDefRateMidLabel, pset);
    new_dw->get(pStress_old, lb->pStressUnrotatedLabel, pset);
    new_dw->get(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

    ParticleVariable<double> pdTdt, p_q;
    ParticleVariable<Matrix3> pStress_new;

    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->allocateAndPut(pStress_new, lb->pStressLabel_preReloc, pset);

    ParticleVariable<Matrix3> pStrain_new;
    ParticleVariable<double> pShearModulus_new, pBulkModulus_new, pCohesion_new;
                             pSuction_new, pSpVol_new, pShearStrain_new,
                             pShearStrainRate_new;

    new_dw->allocateAndPut(pStrain_new, pStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pShearModulus_new, pShearModulusLabel_preReloc, pset);
    new_dw->allocateAndPut(pBulkModulus_new, pBulkModulusLabel_preReloc, pset);
    new_dw->allocateAndPut(pCohesion_new, pCohesionLabel_preReloc, pset);
    new_dw->allocateAndPut(pSuction_old, pSuctionLabel_preReloc, pset);
    new_dw->allocateAndPut(pSpVol_old, pSpVolLabel_preReloc, pset);
    new_dw->allocateAndPut(pShearStrain_old, pShearStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pShearStrainRate_old, pShearStrainRateLabel_preReloc, pset);

    for (idx : *pset) {

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // This is unrotated stress and rate of deformation
      // pStress = R^T * pStress * R, D = R^T * D * R
      // Load into 1-D array
      Vector6 stress_vec = Uintah::toVector6(pStress_old[idx]);
      Vector6 D_vec      = Uintah::toVector6(pDefRate_mid[idx]);

      // get the volumetric part of the deformation
      double J = pDefGrad_new[idx].Determinant();

      // Compute the local sound speed
      double rho_cur = rho_orig / J;

      // Shear strain and shear strain rate
      Vector6 strain_vec = Uintah::toVector6(pStrain_old[idx]);
      strain_vec += D_vec * delT;
      pStrain_new[idx] = strain_vec;

      double shearStrain     = computeShearStrain(strain_vec);
      double shearStrainRate = computeShearStrain(D_vec);

      // Non local correction
      bool nonlocal = d_flags.useNonlocalCorrection;
      if (nonlocal) {
        std::tie(shearStrain, shearStrainRate) =   
         computeNonlocalShearStrains(matID, patch, old_dw,
                                     pX[idx], pShearModulus[idx], pCohesion[idx],
                                     shearStrain, shearStrainRate);
      }

      pShearStrain_new[idx] = shearStrain;
      pShearStrainRate_new[idx] = shearStrainRate;

      // Update internal variables and do integration
      Vector6 strainInc = -D_vec * delT;
      double cohesion = pCohesion_old[idx];
      double shearModulus = pShearModulus_old[idx];
      double bulkModulus = pBulkModulus_old[idx];
      double suction = pSuction_old[idx];
      double specificVol = pSpVol_old[idx];
      calculateStress(pX[idx], strainInc, shearStrain, shearStrainRate,
                      stress_vec, cohesion, shearModulus, bulkModulus,
                      suction, specificVol);

      pCohesion_new[idx] = cohesion;
      pShearModulus_new[idx] = shearModulus;
      pBulkModulus_new[idx] = bulkModulus;
      pSuction_new[idx] = suction;
      pSpVol_new[idx] = specificVol;

      // This is the Cauchy stress, still unrotated
      tensorSig(0, 0) = pStress_vec[0];
      tensorSig(1, 1) = pStress_vec[1];
      tensorSig(2, 2) = pStress_vec[2];
      tensorSig(0, 1) = pStress_vec[3];
      tensorSig(1, 0) = pStress_vec[3];
      tensorSig(2, 1) = pStress_vec[4];
      tensorSig(1, 2) = pStress_vec[4];
      tensorSig(2, 0) = pStress_vec[5];
      tensorSig(0, 2) = pStress_vec[5];

      // ROTATE pStress_new: S=R*tensorSig*R^T
      pStress_new[idx] = (tensorR * tensorSig) * (tensorR.Transpose());

      // Compute Ko
      double s_xx = stress_vec[0];
      double s_yy = stress_vec[1];
      double s_xy = stress_vec[3];

      pISV_vec[39] = s_xx;
      pISV_vec[40] = s_yy;
      pISV_vec[45] = s_xy;
      double Ko    = s_xx / s_yy;
      pISV_vec[41] = Ko;

      /*
    // Non local stress
      double n_nonlocal = UI[46];
      double l_nonlocal = UI[47];

      double domain_nonlocal = l_nonlocal * l_nonlocal;
      double Up0 = 0;
      double Up1 = 0;
      double Up2 = 0;
      double Up3 = 0;
      double Up4 = 0;
      double Up5 = 0;
      double Down = 0;
      //double Uptest = 0;

      double weight = 0;
      double rx = 0;
      double ry = 0;
      double rz = 0;
      double r2 = 0;

      double s11 = pStress_new[idx](0,0);
      double s22 = pStress_new[idx](1,1);
      double s33 = pStress_new[idx](2,2);
      double s12 = pStress_new[idx](0,1);
      double s23 = pStress_new[idx](1,2);
      double s13 = pStress_new[idx](2,0);

      for (ParticleSubset::iterator iter1 = pset->begin();
              iter1 != pset->end(); iter1++) {
              particleIndex idx1 = *iter1;

              rx = pX_new[idx1].x() - pX_new[idx].x();
              ry = pX_new[idx1].y() - pX_new[idx].y();
              rz = pX_new[idx1].z() - pX_new[idx].z();

              r2 = rx * rx + ry * ry + rz * rz;

              if (r2 <= (9*domain_nonlocal)) {

                      weight =
  sqrt(r2)*exp(-r2/domain_nonlocal)/domain_nonlocal;
              }
              Up0 += s11 * weight*pVolume_new[idx1];
              Up1 += s22 * weight*pVolume_new[idx1];
              Up2 += s33 * weight*pVolume_new[idx1];
              Up3 += s12 * weight*pVolume_new[idx1];
              Up4 += s23 * weight*pVolume_new[idx1];
              Up5 += s13 * weight*pVolume_new[idx1];
              Down += (weight*pVolume_new[idx1]);

             //Uptest += idx1 * weight*pVolume_new[idx1];

             //cerr << weight << ' ' << r2 << ' ' << domain_nonlocal << endl;

      }

      double Snonlocal[6];
      Snonlocal[0] = (1 - n_nonlocal)*s11 + (n_nonlocal*Up0/Down);
      Snonlocal[1] = (1 - n_nonlocal)*s22 + (n_nonlocal*Up1/Down);
      Snonlocal[2] = (1 - n_nonlocal)*s33 + (n_nonlocal*Up2/Down);
      Snonlocal[3] = (1 - n_nonlocal)*s12 + (n_nonlocal*Up3/Down);
      Snonlocal[4] = (1 - n_nonlocal)*s23 + (n_nonlocal*Up4/Down);
      Snonlocal[5] = (1 - n_nonlocal)*s13 + (n_nonlocal*Up5/Down);

      //double test;
  //test = (1 - n_nonlocal) * idx + (n_nonlocal * Uptest / Down);

  //cerr << n_nonlocal << ' ' << Uptest << ' ' << Down << endl;
      //cerr << test << ' ' << idx <<endl;

      pStress_new[idx](0,0) = Snonlocal[0];
      pStress_new[idx](1,2) = Snonlocal[1];
      pStress_new[idx](2,2) = Snonlocal[2];
      pStress_new[idx](0,1) = Snonlocal[3];
      pStress_new[idx](1,0) = Snonlocal[3];
      pStress_new[idx](1,2) = Snonlocal[4];
      pStress_new[idx](2,1) = Snonlocal[4];
      pStress_new[idx](2,0) = Snonlocal[5];
      pStress_new[idx](0,2) = Snonlocal[5];
      */

      c_dil = sqrt(USM / rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pStress_new[idx] + pStress[idx]) * .5;

      double e = (D(0, 0) * AvgStress(0, 0) + D(1, 1) * AvgStress(1, 1) +
                  D(2, 2) * AvgStress(2, 2) +
                  2. * (D(0, 1) * AvgStress(0, 1) + D(0, 2) * AvgStress(0, 2) +
                        D(1, 2) * AvgStress(1, 2))) *
                 pVolume_new[idx] * delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pVelocity_idx = pVelocity[idx];
      waveSpeed = Vector(Max(c_dil + fabs(pVelocity_idx.x()), waveSpeed.x()),
                         Max(c_dil + fabs(pVelocity_idx.y()), waveSpeed.y()),
                         Max(c_dil + fabs(pVelocity_idx.z()), waveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(UI[1] / rho_cur);
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    waveSpeed       = dx / waveSpeed;
    double delT_new = waveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }
  }
}

double
MohrCoulomb::computeShearStrain(const Vector6& strain) const
{
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
MohrCoulomb::computeNonlocalShearStrains(int matID, const Patch* patch, 
                            DataWarehouse* old_dw,
                            const Point& pX,
                            double pShearModulus,
                            double pCohesion,
                            double shearStrain,
                            double shearStrainRate)
{
  Vector cellSize = patch->dCell();
  double nonlocal_length = d_params.nonlocalL;
  Vector nonlocal_length_vec(nonlocal_length, nonlocal_length, nonlocal_length);

  Point cellIndexLow = patch->getLevel()->positionToIndex(pX - nonlocal_length_vec);
  Point cellIndexHigh = patch->getLevel()->positionToIndex(pX + nonlocal_length_vec);

  ParticleSubset* pset = old_dw->getParticleSubset(matID, patch,
                                                   cellIndexLow,
                                                   cellIndexHigh);

  constParticleVariable<Point> pX_neighbor;
  old_dw->get(pX_neighbor, lb->pXLabel, pset);

  constParticleVariable<double> pVolume_neighbor;
  new_dw->get(pVolume_neighbor, lb->pVolumeLabel_preReloc, pset);


  double domain_nonlocal = nonlocal_length * nonlocal_length;
  double up              = 0;
  double up1             = 0;
  double down            = 0;

  for (auto idx : *pset) {

    double distSq = (pX_neighbor[idx] - pX).lengthSq();
    double dist = std::sqrt(distSq);
    double weight = dist * std::exp(-distSq / domain_nonlocal) / domain_nonlocal;
    double wieghtedVolume = weight * pVolume_neighbor[idx];
    up += (shearStrain * weightedVolume)
    up1 += (shearStrainRate * weightedVolume);
    down += weightedVolume;
  }

  double n_nonlocal = d_params.nonlocalN;
  double shearStrain_nonlocal =
        (1 - n_nonlocal) * shearStrain + (n_nonlocal * up / down);
  double shearStrainRate_nonlocal =
        (1 - n_nonlocal) * shearStrainRate + (n_nonlocal * up1 / down);

  double elastic_strain = pCohesion / pShearModulus;
  bool regularize = d_flags.useRegularizedNonlocalSoftening;
  if (regularize) {
    double tFE = d_params.regularizationTFE;
    double tShear = d_params.regularizationTShear;
    if (shearStrain_nonlocal > elastic_strain) {
      shearStrain_nonlocal     = shearStrain_nonlocal * tFE / tShear;
      shearStrainRate_nonlocal = shearStrainRate_nonlocal * tFE / tShear;
    }
  }

  return std::make_tuple(shearStrain_nonlocal, shearStrainRate_nonlocal);
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
MohrCoulomb::calculateStress(const Point& pX,
                             const Vector6& strainInc,
                             double pShearStrain,
                             double pShearStrainRate,
                             Vector6& stress,
                             double& pCohesion,
                             double& pShearModulus,
                             double& pBulkModulus,
                             double& pSuction,
                             double& pSpecificVol)
{
  // Convert to compression +ve form
  MohrCoulombState initialState;
  initialState.stress = -stress;

  // Rate dependence
  if (d_flags.useUndrainedShearTransition) {

    double St = d_params.softeningSt;
    double a1 = d_params.waterInfluenceA1;
    double b1 = d_params.waterInfluenceB1;
    double W = d_params.waterInfluenceW;

    pCohesion = St * a1 * pow(W, -b1);
    if (pShearStrainRate > d_params.refStrainRate) {
      pCohesion *= 
          std::pow(pShearStrainRate / d_params.refStrainRate, d_params.betaStrainRate);
    } 
  }

  // Shear strength linear with depth
  if (d_flags.useLinearlyVaryingCohesion) {
    std::string direction = d_params.depthDirection;
    double yDiff = 0;
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
    pShearModulus = d_params.variableModulusM * pCohesion / 2.0 / (1.0 + d_params.variableModulusNuY);
    pBulkModulus = d_params.variableModulusM * pCohesion / 3.0 / (1.0 - 2 * d_params.variableModulusNuY);
  }

  // Softening
  if (d_flags.useSoftening) {
    if (pShearStrain > pCohesion / pShearModulus) {
      double St = d_params.softeningSt;
      pCohesion *= (1.0 / St + (1.0 - 1.0 / St) *
                   std::pow(2.71, (-3.0 * pShearStrain / d_params.softeningStrain95)));
    }
  }

  // Retention
  if (d_flags.useWaterRetention) {

    // get from current suction value
    double Sr = computeSr(pSuction);

    // constant volume of water - thus we use initial specific volume
    double volWater = Sr * (pSpecificVol - 1.0) / pSpecificVol; 
    double dVolStrain = strainInc(0) + strainInc(1) + strainInc(2);

    // volume of voids - new one, using new specific volume
    double specVol_new = pSpecificVol - dVolStrain * pSpecificVol; 
    double volAir = (specVol_new - 1) / specVol_new - volWater; 

    // voids will reduce by dvol, so dSr=VolWater/TotalVolume
    double totalVolume = volAir + volWater;
    double Sr_new = volWater / totalVolume; 
    if (Sr_new > 0.99999) {
      Sr_new = 0.99999;
    } else if (Sr_new < 0.0001) {
      Sr_new = 0.0001;
    }

    double suction_new = computeSuction(SrNew);
    pCohesion += std::tan(d_params.phi_b * M_PI / 180.0) * suction_new;

    pSuction = suction_new;
    pSpecificVol = specVol_new;
  }

  // Do the integration
  d_modelP->setModelParameters(pShearModulus, pBulkModulus, pCohesion,
                               d_params.phi, d_params.psi, d_params.pMin);

  MohrCoulombState finalState = initialState;
  if (d_modelType == "classic" || "sheng") {
    finalState = d_modelP->integrate(strainInc, initialState);
  } else if (d_modelType == "classic_semiimplicit") {
    RegionType region;
    finalState = d_modelP->integrate(strainInc, initialState, region);
  } else {
    std::ostringstream err;
    err << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw InvalidValue(err.str(), __FILE__, __LINE__);
  }

  // stress back for output (tension +ve)
  stress = -finalState.stress;
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
    double m     = d_params.waterRetentionParams[0];
    double n     = d_params.waterRetentionParams[1];
    double alpha = d_params.waterRetentionParams[2];
    double suction = std::pow(Sr, 1 / m);
    suction        = (1 - suction) / suction;
    suction        = std::pow(suction, 1 / n);
    suction        = suction / alpha;
    if (suction < 0) {
      suction = 0.0;
    }
    return suction;
}

void
MohrCoulomb::carryForward(const PatchSubset* patches,
                          const MPMMaterial* matl,
                          DataWarehouse* old_dw,
                          DataWarehouse* new_dw)
{
  unsigned int num_ISV = ISVLabels.size();

  for (int p = 0; p < patches->size(); p++) {

    const Patch* patch   = patches->get(p);
    int dwi              = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    std::vector<constParticleVariable<double>> ISVs(num_ISV + 1);
    std::vector<ParticleVariable<double>> ISVs_new(num_ISV + 1);

    for (int i = 0; i < num_ISV; i++) {
      old_dw->get(ISVs[i], ISVLabels[i], pset);
      new_dw->allocateAndPut(ISVs_new[i], ISVLabels_preReloc[i], pset);
      ISVs_new[i].copyData(ISVs[i]);
    }

    // Don't affect the strain energy or timestep size
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.), lb->StrainEnergyLabel);
    }
  }
}

double
MohrCoulomb::computeRhoMicroCM(double pressure,
                               const double p_ref,
                               const MPMMaterial* matl,
                               double temperature,
                               double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge  = pressure - p_ref;
  double rho_cur;
  double bulk = UI[1];

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;

#if 1
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR MohrCoulomb" << endl;
#endif
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

  double bulk     = UI[1];
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure   = p_ref + p_g;
  dp_drho    = bulk * rho_orig / (rho_cur * rho_cur);
  tmp        = bulk / rho_cur; // speed of sound squared

#if 1
  cout << "NO VERSION OF computePressEOSCM EXISTS YET FOR MohrCoulomb" << endl;
#endif
}

double
MohrCoulomb::getCompressibility()
{
  return 1.0 / UI[1];

  // return 1; //experimental:1 // found to be irrelevant
}


