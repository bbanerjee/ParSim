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

#include <CCA/Components/MPM/ConstitutiveModel/MohrCoulomb.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ClassicMohrCoulomb.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ShengMohrCoulomb.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

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

#include <sci_defs/uintah_defs.h>

#include <iostream>
#include <string>

// sets of external variables for the Sheng Mohr Coulomb algorithm by WTS. Some
// are redundant.

double ALFACHECK, ALFACHANGE, ALFARATIO, MAXITER, YIELDTOL, TOL_METHOD,
  INCREMENT_TYPE, BETA_FACT;
double DRIFT_CORRECTION, EULER_ITERATIONS, CRITICAL_STEP_SIZE;
double STEP_MAX, STEP_MIN, ERROR_DEF, USE_ERROR_STEP, MIN_DIVISION_SIZE;
double INTEGRATION_TOL = 0.0001;
double LINES_SEC_NO = 50;
int MAX_LOOP, SOLUTION_ALGORITHM,
  ALGORITHM_TYPE; // MAX_LOOP - number of steps to solve, SOLUTION_ALGORITHM -
                  // algorithm to use
double USE_ERROR = 0.5,
       SAVE_FOR_ERROR = 1; // these are values - 1st - when 'use' the additional
                           // error to speed up the calculations (here 30%)
// SAVEFORERROR says how greater accuracy the program should use; 1 means no
// change. 0.5 means 50% greater accuracy. (tol * 0.5) etc
double CHGEPSINTOL = 10e-9;
double ADDTOLYIELD = 0.8;
double SUCTIONTOL = 0.00000001;
double TINY = 1e-14; // value used for example in checking whether we are not
                     // dividing by zero in CalcStressElast, volumetric strain
                     // must be also larger then tiny
double PMIN =
  0.0001; // value of minimum mean stress to calculate K in CalcStressElast
int USE_NICE_SCHEME = 0;

using namespace Uintah;

MohrCoulomb::MohrCoulomb(ProblemSpecP& ps, MPMFlags* Mflag)
  : ConstitutiveModel(Mflag)
{
  // Read model parameters from the input file
  getInputParameters(ps);

  // Create the model
  if (d_modelType == "classic") {
    
    d_modelP = std::make_unique<MohrCoulombClassic>(d_params.G, d_params.K,
                                                    d_params.cohesion, 
                                                    d_params.phi, d_params.psi);

  } else if (d_modelType == "sheng") {

    d_modelP = std::make_unique<MohrCoulombSheng>(d_params.G, d_params.K,
                                                  d_params.cohesion, 
                                                  d_params.phi, d_params.psi);

  } else {
    std::ostringstream err;
    err << "Version of the Mohr-Coulomb Model is set to: " << d_modelType
        << " This will cause the code to malfuncion. Any results obtained are "
           "invalid.\n";
    throw ProblemSetupException(err.str(), __FILE__, __LINE__);
  }

  // Create VarLabels for GeoModel internal state variables (ISVs)
  int nx;
  nx = d_NBASICINPUTS;

  for (int i = 0; i < nx; i++) {
    rinit[i] = UI[i];
    // cerr<<" UI["<<i<<"]="<<UI[i];
  }
  d_NINSV = nx;
  //  cout << "d_NINSV = " << d_NINSV <<l endl;

  initializeLocalMPMLabels();
}

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
        << " Include the <model_parameters> tag in the input file.\n"
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  
  d_params.G = 10000.0; d_params.K = 20000.0; 
  ps->require("shear_modulus", d_params.G);
  ps->require("bulk_modulus",  d_params.K);

  d_params.c = 0.0; d_params.phi = 30.0; d_params.psi = 30.0;
  ps->require("cohesion", d_params.c);
  ps->require("angle_internal_friction", d_params.phi);
  ps->require("angle_dilation", d_params.psi);

  d_params.initialSuction = 0.0; d_params.phi_b = 0.0;;
  ps->getWithDefault("initial_suction", d_params.initialSuction, 0.0);
  ps->getWithDefault("phi_b", d_params.phi_b, 0.0);

  // whether water retention curve should be used
  d_flags.useWaterRetention = false;
  ps->getWithDefault("use_water_retention", d_flags.useWaterRetention, false); 

  d_params.waterRetentionParams = {0.0, 0.0, 0.0, 0.0};
  ps->getWithDefault("water_retention_param_1", d_params.waterRetentionParams[0], 0.0);
  ps->getWithDefault("water_retention_param_2", d_params.waterRetentionParams[1], 0.0);
  ps->getWithDefault("water_retention_param_3", d_params.waterRetentionParams[2], 0.0);
  ps->getWithDefault("water_retention_param_4", d_params.waterRetentionParams[3], 0.0);

  // whether undrained shear strength transition should be used
  d_flags.useUndrainedShearTransition = false;
  ps->getWithDefault("use_undrained_shear_transition", d_flags.useUndrainedShearTransition, false);

  // water influence parameters 
  d_params.waterInfluenceA1 = 0.0; d_params.waterInfluenceB1 = 0.0; 
  d_params.waterInfluenceW = 0.0;
  ps->getWithDefault("water_influence_A1", d_params.waterInfluenceA1, 0.0); // water influence parameter
  ps->getWithDefault("water_influence_B1", d_params.waterInfluenceB1, 0.0); // water influence parameter
  ps->getWithDefault("water_influence_W", d_params.waterInfluenceW, 0.0);  // water content

  // strain rate influence parameters 
  d_params.betaStrainRate = 0.0; d_params.refStrainRate = 0.0; 
  d_params.shearStrainRate = 0.0;
  ps->getWithDefault("beta_strain_rate", d_params.betaStrainRate, 0.0); 
  ps->getWithDefault("ref_strain_rate", d_params.refStrainRate, 0.0);
  ps->getWithDefault("shear_strain_rate", d_params.shearStrainRate, 0.0);

  // use variable elastic modulus
  d_flags.useVariableElasticModulus = false;
  ps->getWithDefault("use_variable_elastic_modulus", d_flags.useVariableElasticModulus, false);

  d_params.variableModulusM = 0.0; d_params.variableModulusNuY = 0.0; 
  d_params.variableModulusShearStrain = 0.0;
  ps->getWithDefault("variable_modulus_m", d_params.variableModulusM, 0.0);
  ps->getWithDefault("variable_modulus_nu_y", d_params.variableModulusNuY, 0.0);
  ps->getWithDefault("variable_modulus_shear_strain", d_params.variableModulusShearStrain, 0.0);

  // use linearly varying cohesion with depth
  d_flags.useLinearlyVaryingCohesion = false;
  ps->getWithDefault("use_linearly_varying_cohesion", d_flags.useLinearlyVaryingCohesion, false);

  d_params.linearCohesionA = 0.0, d_params.linearCohesionYRef = 0.0;
  ps->getWithDefault("linear_cohesion_a", d_params.linearCohesionA, 0.0);
  ps->getWithDefault("linear_cohesion_y_ref", d_params.linearCohesionYRef, 0.0);

  // use softening model
  d_flags.useSoftening = false;
  ps->getWithDefault("use_softening", d_flags.useSoftening, false);

  d_params.softeningSt = 0.0, d_params.softeningStrain95 = 0.0;
  ps->getWithDefault("softening_St", d_params.softeningSt, 0.0);
  ps->getWithDefault("softening_strain_95", d_params.softeningStrain95, 0.0);

  // use regularized nonlocal softening
  d_flags.useRegularizedNonlocalSoftening = false;
  ps->getWithDefault("use_regularized_nonlocal_softening", d_flags.useRegularizedNonlocalSoftening, false);

  d_params.regularizationTFE = 0.0, d_params.regularizationTShear = 0.0;
  ps->getWithDefault("regularization_t_FE", d_params.regularizationTFE, 0.0);
  ps->getWithDefault("regularization_t_shear", d_params.regularizationTShear, 0.0);

  // use nonlocal correction
  d_flags.useNonlocalCorrection = false;
  ps->getWithDefault("use_nonlocal_correction", d_flags.useNonlocalCorrection, false);

  d_params.nonlocalN = 0.0, d_params.nonlocalL = 0.0;
  ps->getWithDefault("nonlocal_n", d_params.nonlocalN, 0.0);
  ps->getWithDefault("nonlocal_l", d_params.nonlocalL, 0.0);

  // Check that model parameters are valid and allow model to change if needed
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
  if (d_params.c < 0.0)
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
    out << "**Error** No integration parameters provided for Mohr-Coulomb model."
        << " Include the <integration_parameters> tag in the input file.\n"
    throw Uintah::ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  
  d_int.maxIter = 200;
  d_int.alfaCheck = 1; 
  d_int.alfaChange = 0.05;
  d_int.alfaRatio = 10;
  ps->getWithDefault("max_iterations_pegasus", d_int.maxIter, 200);
  ps->getWithDefault("alpha_check_pegasus", d_int.alfaCheck, 1);
  ps->getWithDefault("alpha_change_pegasus", d_int.alfaChange, 0.05);
  ps->getWithDefault("alpha_ratio_pegasus", d_int.alfaRatio, 10);

  d_int.yieldTol = 1.0e-6;
  d_int.integrationTol = 0.01;
  d_int.betaFactor = 0.9; 
  d_int.minMeanStress = -1.0e8;
  d_int.suctionTol = 1.0e-8; 
  ps->getWithDefault("yield_tolerance", d_int.yieldTol, 1.0e-6);
  ps->getWithDefault("integration_tolerance", d_int.integrationTol, 0.01);
  ps->getWithDefault("beta_safety_factor", d_int.betaFactor, 0.9); 
  ps->getWithDefault("minimum_mean_stress", d_int.minMeanStress, -1.0e8);
  ps->getWithDefault("suction_tolerance", d_int.suctionTol, 1.0e-8); 

  d_int.driftCorrection = "at_end";
  d_int.tolMethod = "sloan";
  d_int.solutionAlgorithm = "modified_euler";
  ps->getWithDefault("drift_correction_algorithm", d_int.driftCorrection, "at_end");
  ps->getWithDefault("tolerance_algorithm", d_int.tolMethod, "sloan");
  ps->getWithDefault("solution_algorithm", d_int.solutionAlgorithm, "modified_euler");
}

void
MohrCoulomb::initializeLocalMPMLabels()
{
  vector<string> ISVNames;

  ISVNames.push_back("G");
  ISVNames.push_back("K");
  ISVNames.push_back("c");
  ISVNames.push_back("Phi");
  ISVNames.push_back("Psi");
  ISVNames.push_back("Version");
  ISVNames.push_back("Suction");
  ISVNames.push_back("UseWaterRetention");
  ISVNames.push_back("WR1");
  ISVNames.push_back("WR2");
  ISVNames.push_back("WR3");
  ISVNames.push_back("WR4");
  ISVNames.push_back("SpecificVol");
  ISVNames.push_back("PhiB");
  ISVNames.push_back("Usetransition");
  ISVNames.push_back("A1");
  ISVNames.push_back("B1");
  ISVNames.push_back("W");
  ISVNames.push_back("beta");
  ISVNames.push_back("strain_ref");
  ISVNames.push_back("shear_strain_rate");
  ISVNames.push_back("Usemodul");
  ISVNames.push_back("m_modul");
  ISVNames.push_back("nuy");
  ISVNames.push_back("shear_strain");

  ISVNames.push_back("Use_linear");
  ISVNames.push_back("a");
  ISVNames.push_back("y_ref");

  ISVNames.push_back("strain11");
  ISVNames.push_back("strain22");
  ISVNames.push_back("strain33");
  ISVNames.push_back("strain12");
  ISVNames.push_back("strain23");
  ISVNames.push_back("strain13");

  ISVNames.push_back("Use_softening");
  ISVNames.push_back("St");
  ISVNames.push_back("strain_95");

  ISVNames.push_back("y");
  ISVNames.push_back("n");

  ISVNames.push_back("s_xx");
  ISVNames.push_back("s_yy");
  ISVNames.push_back("Ko");

  ISVNames.push_back("Use_regular");
  ISVNames.push_back("tFE");
  ISVNames.push_back("tShear");
  ISVNames.push_back("s_xy");

  ISVNames.push_back("n_nonlocalMC");
  ISVNames.push_back("l_nonlocal");

  for (int i = 0; i < d_NINSV; i++) {
    ISVLabels.push_back(VarLabel::create(
      ISVNames[i], ParticleVariable<double>::getTypeDescription()));
    ISVLabels_preReloc.push_back(VarLabel::create(
      ISVNames[i] + "+", ParticleVariable<double>::getTypeDescription()));
  }
}
#if 0
MohrCoulomb::MohrCoulomb(const MohrCoulomb* cm) : ConstitutiveModel(cm)
{
  for(int i=0;i<d_NDMMPROP;i++){
    UI[i] = cm->UI[i];
  }

  //Create VarLabels for Diamm internal state variables (ISVs)
  initializeLocalMPMLabels();
}
#endif

MohrCoulomb::~MohrCoulomb()
{
  for (unsigned int i = 0; i < ISVLabels.size(); i++) {
    VarLabel::destroy(ISVLabels[i]);
  }
}

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
  cm_ps->appendElement("bulk_modulus",  d_params.K);

  cm_ps->appendElement("cohesion", d_params.c);
  cm_ps->appendElement("angle_internal_friction", d_params.phi);
  cm_ps->appendElement("angle_dilation", d_params.psi);

  cm_ps->appendElement("initial_suction", d_params.initialSuction);
  cm_ps->appendElement("phi_b", d_params.phi_b);

  cm_ps->appendElement("use_water_retention", d_flags.useWaterRetention);

  cm_ps->appendElement("water_retention_param_1", d_params.waterRetentionParams[0]);
  cm_ps->appendElement("water_retention_param_2", d_params.waterRetentionParams[1]);
  cm_ps->appendElement("water_retention_param_3", d_params.waterRetentionParams[2]);
  cm_ps->appendElement("water_retention_param_4", d_params.waterRetentionParams[3]);

  cm_ps->appendElement("use_undrained_shear_transition", d_flags.useUndrainedShearTransition);

  cm_ps->appendElement("water_influence_A1", d_params.waterInfluenceA1);
  cm_ps->appendElement("water_influence_B1", d_params.waterInfluenceB1);
  cm_ps->appendElement("water_influence_W", d_params.waterInfluenceW);

  cm_ps->appendElement("beta_strain_rate", d_params.betaStrainRate);
  cm_ps->appendElement("ref_strain_rate", d_params.refStrainRate);
  cm_ps->appendElement("shear_strain_rate", d_params.shearStrainRate);

  cm_ps->appendElement("use_variable_elastic_modulus", d_flags.useVariableElasticModulus);

  cm_ps->appendElement("variable_modulus_m", d_params.variableModulusM);
  cm_ps->appendElement("variable_modulus_nu_y", d_params.variableModulusNuY);
  cm_ps->appendElement("variable_modulus_shear_strain", d_params.variableModulusShearStrain);

  cm_ps->appendElement("use_linearly_varying_cohesion", d_flags.useLinearlyVaryingCohesion);

  cm_ps->appendElement("linear_cohesion_a", d_params.linearCohesionA);
  cm_ps->appendElement("linear_cohesion_y_ref", d_params.linearCohesionYRef);

  cm_ps->appendElement("use_softening", d_flags.useSoftening);

  cm_ps->appendElement("softening_St", d_params.softeningSt);
  cm_ps->appendElement("softening_strain_95", d_params.softeningStrain95);

  cm_ps->appendElement("use_regularized_nonlocal_softening", d_flags.useRegularizedNonlocalSoftening);

  cm_ps->appendElement("regularization_t_FE", d_params.regularizationTFE);
  cm_ps->appendElement("regularization_t_shear", d_params.regularizationTShear);

  cm_ps->appendElement("use_nonlocal_correction", d_flags.useNonlocalCorrection);

  cm_ps->appendElement("nonlocal_n", d_params.nonlocalN);
  cm_ps->appendElement("nonlocal_l", d_params.nonlocalL);
}

void
MohrCoulomb::outputIntegrationProblemSpec(ProblemSpecP& ps) const
{
  ProblemSpecP cm_ps = ps->appendChild("integration_parameters");
}

MohrCoulomb*
MohrCoulomb::clone()
{
  return scinew MohrCoulomb(*this);
}

void
MohrCoulomb::initializeCMData(const Patch* patch, const MPMMaterial* matl,
                              DataWarehouse* new_dw)
{
  // Initialize the variables shared by all constitutive models
  // This method is defined in the ConstitutiveModel base class.
  initSharedDataForExplicit(patch, matl, new_dw);

  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  std::vector<ParticleVariable<double>> ISVs(d_NINSV + 1);

  cout << "In initializeCMData" << endl;
  for (int i = 0; i < d_NINSV; i++) {
    new_dw->allocateAndPut(ISVs[i], ISVLabels[i], pset);
    ParticleSubset::iterator iter = pset->begin();
    for (; iter != pset->end(); iter++) {
      ISVs[i][*iter] = rinit[i];
    }
  }

  computeStableTimestep(patch, matl, new_dw);
}

void
MohrCoulomb::addParticleState(std::vector<const VarLabel*>& from,
                              std::vector<const VarLabel*>& to)
{
  // Add the local particle state data for this constitutive model.
  for (int i = 0; i < d_NINSV; i++) {
    from.push_back(ISVLabels[i]);
    to.push_back(ISVLabels_preReloc[i]);
  }
}

void
MohrCoulomb::computeStableTimestep(const Patch* patch, const MPMMaterial* matl,
                                   DataWarehouse* new_dw)
{
  // This is only called for the initial timestep - all other timesteps
  // are computed as a side-effect of computeStressTensor
  Vector dx = patch->dCell();
  int dwi = matl->getDWIndex();
  ParticleSubset* pset = new_dw->getParticleSubset(dwi, patch);
  constParticleVariable<double> pmass, pvolume;
  constParticleVariable<Vector> pvelocity;

  new_dw->get(pmass, lb->pMassLabel, pset);
  new_dw->get(pvolume, lb->pVolumeLabel, pset);
  new_dw->get(pvelocity, lb->pVelocityLabel, pset);

  double c_dil = 0.0;
  Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);

  double bulk = UI[1];
  double G = UI[0]; // modified: K=UI[1], G=UI[0]
  for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
       iter++) {
    particleIndex idx = *iter;

    // Compute wave speed at each particle, store the maximum
    c_dil = sqrt((bulk + 4. * G / 3.) * pvolume[idx] / pmass[idx]);
    WaveSpeed = Vector(Max(c_dil + fabs(pvelocity[idx].x()), WaveSpeed.x()),
                       Max(c_dil + fabs(pvelocity[idx].y()), WaveSpeed.y()),
                       Max(c_dil + fabs(pvelocity[idx].z()), WaveSpeed.z()));
  }
  // UI[14]=matl->getInitialDensity();
  // UI[15]=matl->getRoomTemperature();
  // UI[14]=bulk/matl->getInitialDensity();  ??tim
  // UI[19]=matl->getInitialCv();
  WaveSpeed = dx / WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

void
MohrCoulomb::computeStressTensor(const PatchSubset* patches,
                                 const MPMMaterial* matl, DataWarehouse* old_dw,
                                 DataWarehouse* new_dw)
{

  double rho_orig = matl->getInitialDensity();
  for (int p = 0; p < patches->size(); p++) {
    double se = 0.0;
    const Patch* patch = patches->get(p);

    Matrix3 Identity;
    Identity.Identity();
    double c_dil = 0.0;
    Vector WaveSpeed(1.e-12, 1.e-12, 1.e-12);
    Vector dx = patch->dCell();

    int dwi = matl->getDWIndex();
    // Create array for the particle position
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    constParticleVariable<Matrix3> deformationGradient, pstress;
    ParticleVariable<Matrix3> pstress_new;
    constParticleVariable<Matrix3> deformationGradient_new, velGrad;
    constParticleVariable<double> pmass, pvolume, ptemperature;
    constParticleVariable<double> pvolume_new;
    constParticleVariable<Vector> pvelocity;
    constParticleVariable<Point> px, pxnew;
    delt_vartype delT;

    old_dw->get(delT, lb->delTLabel, getLevel(patches));

    old_dw->get(px, lb->pXLabel, pset);
    old_dw->get(pstress, lb->pStressLabel, pset);
    old_dw->get(pmass, lb->pMassLabel, pset);
    old_dw->get(pvolume, lb->pVolumeLabel, pset);
    old_dw->get(pvelocity, lb->pVelocityLabel, pset);
    old_dw->get(ptemperature, lb->pTemperatureLabel, pset);
    old_dw->get(deformationGradient, lb->pDeformationMeasureLabel, pset);
    new_dw->get(pvolume_new, lb->pVolumeLabel_preReloc, pset);
    new_dw->get(pxnew, lb->pXLabel_preReloc, pset);

    std::vector<constParticleVariable<double>> ISVs(d_NINSV + 1);
    for (int i = 0; i < d_NINSV; i++) {
      old_dw->get(ISVs[i], ISVLabels[i], pset);
    }

    ParticleVariable<double> pdTdt, p_q;

    new_dw->allocateAndPut(pstress_new, lb->pStressLabel_preReloc, pset);
    new_dw->allocateAndPut(pdTdt, lb->pdTdtLabel, pset);
    new_dw->allocateAndPut(p_q, lb->p_qLabel_preReloc, pset);
    new_dw->get(deformationGradient_new, lb->pDeformationMeasureLabel_preReloc,
                pset);
    new_dw->get(velGrad, lb->pVelGradLabel_preReloc, pset);

    std::vector<ParticleVariable<double>> ISVs_new(d_NINSV + 1);
    for (int i = 0; i < d_NINSV; i++) {
      new_dw->allocateAndPut(ISVs_new[i], ISVLabels_preReloc[i], pset);
    }

    for (ParticleSubset::iterator iter = pset->begin(); iter != pset->end();
         iter++) {
      particleIndex idx = *iter;

      // Assign zero internal heating by default - modify if necessary.
      pdTdt[idx] = 0.0;

      // Calculate rate of deformation D, and deviatoric rate DPrime,
      Matrix3 D = (velGrad[idx] + velGrad[idx].Transpose()) * .5;

      // get the volumetric part of the deformation
      double J = deformationGradient_new[idx].Determinant();
      // Check 1: Look at Jacobian
      if (!(J > 0.0)) {
        cerr << getpid();
        constParticleVariable<long64> pParticleID;
        old_dw->get(pParticleID, lb->pParticleIDLabel, pset);
        cerr << "**ERROR** Negative Jacobian of deformation gradient"
             << " in particle " << pParticleID[idx] << endl;
        cerr << "l = " << velGrad[idx] << endl;
        cerr << "F_old = " << deformationGradient[idx] << endl;
        cerr << "F_new = " << deformationGradient_new[idx] << endl;
        cerr << "J = " << J << endl;
        throw InternalError("Negative Jacobian", __FILE__, __LINE__);
      }

      // Compute the local sound speed
      double rho_cur = rho_orig / J;

      // NEED TO FIND R
      Matrix3 tensorR, tensorU;

      // Look into using Rebecca's PD algorithm
      deformationGradient_new[idx].polarDecompositionRMB(tensorU, tensorR);

      // This is the previous timestep Cauchy stress
      // unrotated tensorSig=R^T*pstress*R
      Matrix3 tensorSig = (tensorR.Transpose()) * (pstress[idx] * tensorR);

      // Load into 1-D array for the fortran code
      double sigarg[6];
      sigarg[0] = tensorSig(0, 0);
      sigarg[1] = tensorSig(1, 1);
      sigarg[2] = tensorSig(2, 2);
      sigarg[3] = tensorSig(0, 1);
      sigarg[4] = tensorSig(1, 2);
      sigarg[5] = tensorSig(2, 0);

      // UNROTATE D: S=R^T*D*R
      D = (tensorR.Transpose()) * (D * tensorR);

      // Load into 1-D array for the fortran code
      double Dlocal[6];
      Dlocal[0] = D(0, 0);
      Dlocal[1] = D(1, 1);
      Dlocal[2] = D(2, 2);
      Dlocal[3] = D(0, 1);
      Dlocal[4] = D(1, 2);
      Dlocal[5] = D(2, 0);
      double svarg[d_NINSV];
      double USM = 9e99;
      double dt = delT;
      int nblk = 1;

      // Load ISVs into a 1D array for fortran code
      for (int i = 0; i < d_NINSV; i++) {
        svarg[i] = ISVs[i][idx];
      }

      // Undrained increase linearly with depth
      double n = svarg[38];

      for (int i = 0; i < n; i++) {
        svarg[37] = px[idx](1);
        n = n - 1;
      }
      svarg[38] = n;

      // Compute Ko
      double s_xx = sigarg[0];
      double s_yy = sigarg[1];
      double s_xy = sigarg[3];

      svarg[39] = s_xx;
      svarg[40] = s_yy;
      svarg[45] = s_xy;
      double Ko = s_xx / s_yy;
      svarg[41] = Ko;

      // SHear strain and shear strain rate
      double shear_strain_local = 0;
      double strain11 = svarg[28];
      double strain22 = svarg[29];
      double strain33 = svarg[30];
      double strain12 = svarg[31];
      double strain23 = svarg[32];
      double strain13 = svarg[33];

      double e11 = Dlocal[0];
      double e22 = Dlocal[1];
      double e33 = Dlocal[2];
      double e12 = Dlocal[3];
      double e23 = Dlocal[4];
      double e13 = Dlocal[5];
      double shear_strain_rate = UI[20];

      double Use_regular = UI[42];
      double tFE = UI[43];
      double tShear = UI[44];

      strain11 += e11 * dt;
      strain22 += e22 * dt;
      strain33 += e33 * dt;
      strain12 += e12 * dt;
      strain23 += e23 * dt;
      strain13 += e13 * dt;

      shear_strain_local =
        1.0 / 2.0 * sqrt(2 * (pow((strain11 - strain22), 2.0) +
                              pow((strain11 - strain33), 2.0) +
                              pow((strain22 - strain33), 2.0)) +
                         3.0 * (pow(strain12, 2.0) + pow(strain13, 2.0) +
                                pow(strain23, 2.0)));

      shear_strain_rate =
        1.0 / 2.0 * sqrt(2 * (pow((e11 - e22), 2) + pow((e11 - e33), 2) +
                              pow((e22 - e33), 2)) +
                         3 * (pow(e12, 2) + pow(e13, 2) + pow(e23, 2)));

      svarg[28] = strain11;
      svarg[29] = strain22;
      svarg[30] = strain33;
      svarg[31] = strain12;
      svarg[32] = strain23;
      svarg[33] = strain13;

      /*
// Non local
double n_nonlocal = UI[46];
double l_nonlocal = UI[47];

double domain_nonlocal = l_nonlocal * l_nonlocal;
double Up = 0;
double Up1 = 0;
double Down = 0;
//double Uptest = 0;

double weight = 0;
double rx = 0;
double ry = 0;
double rz = 0;
double r2 = 0;

for (ParticleSubset::iterator iter1 = pset->begin();
      iter1 != pset->end(); iter1++) {
      particleIndex idx1 = *iter1;

      rx = pxnew[idx1].x() - pxnew[idx].x();
      ry = pxnew[idx1].y() - pxnew[idx].y();
      rz = pxnew[idx1].z() - pxnew[idx].z();

      r2 = rx * rx + ry * ry + rz * rz;

      if (r2 <= (9*domain_nonlocal)) {

              weight = sqrt(r2)*exp(-r2/domain_nonlocal)/domain_nonlocal;
      }
      Up += shear_strain_local * weight*pvolume_new[idx1];
      Up1 += shear_strain_rate * weight*pvolume_new[idx1];
      Down += (weight*pvolume_new[idx1]);
}

double shear_strain_nonlocal = 0;
double shear_strain_rate_nonlocal = 0;
shear_strain_nonlocal = (1 - n_nonlocal)*shear_strain_local +
(n_nonlocal*Up/Down);
shear_strain_rate_nonlocal = (1 - n_nonlocal)*shear_strain_rate +
(n_nonlocal*Up1 / Down);


*/

      double elastic_strain = svarg[2] / svarg[0];
      double shear_strain_nonlocal = shear_strain_local;
      double shear_strain_rate_nonlocal = shear_strain_rate;

      if (Use_regular > 0) {
        if (shear_strain_nonlocal > elastic_strain) {
          shear_strain_nonlocal = shear_strain_nonlocal * tFE / tShear;
          shear_strain_rate_nonlocal =
            shear_strain_rate_nonlocal * tFE / tShear;
        }
      }

      svarg[24] = shear_strain_local;
      svarg[20] = shear_strain_rate;

      // Calling the external model here
      CalculateStress(nblk, d_NINSV, dt, UI, sigarg, Dlocal, svarg, USM,
                      shear_strain_nonlocal, shear_strain_rate_nonlocal);

      // Unload ISVs from 1D array into ISVs_new
      for (int i = 0; i < d_NINSV; i++) {
        ISVs_new[i][idx] = svarg[i];
      }

      // This is the Cauchy stress, still unrotated
      tensorSig(0, 0) = sigarg[0];
      tensorSig(1, 1) = sigarg[1];
      tensorSig(2, 2) = sigarg[2];
      tensorSig(0, 1) = sigarg[3];
      tensorSig(1, 0) = sigarg[3];
      tensorSig(2, 1) = sigarg[4];
      tensorSig(1, 2) = sigarg[4];
      tensorSig(2, 0) = sigarg[5];
      tensorSig(0, 2) = sigarg[5];

      // ROTATE pstress_new: S=R*tensorSig*R^T
      pstress_new[idx] = (tensorR * tensorSig) * (tensorR.Transpose());

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

      double s11 = pstress_new[idx](0,0);
      double s22 = pstress_new[idx](1,1);
      double s33 = pstress_new[idx](2,2);
      double s12 = pstress_new[idx](0,1);
      double s23 = pstress_new[idx](1,2);
      double s13 = pstress_new[idx](2,0);

      for (ParticleSubset::iterator iter1 = pset->begin();
              iter1 != pset->end(); iter1++) {
              particleIndex idx1 = *iter1;

              rx = pxnew[idx1].x() - pxnew[idx].x();
              ry = pxnew[idx1].y() - pxnew[idx].y();
              rz = pxnew[idx1].z() - pxnew[idx].z();

              r2 = rx * rx + ry * ry + rz * rz;

              if (r2 <= (9*domain_nonlocal)) {

                      weight =
  sqrt(r2)*exp(-r2/domain_nonlocal)/domain_nonlocal;
              }
              Up0 += s11 * weight*pvolume_new[idx1];
              Up1 += s22 * weight*pvolume_new[idx1];
              Up2 += s33 * weight*pvolume_new[idx1];
              Up3 += s12 * weight*pvolume_new[idx1];
              Up4 += s23 * weight*pvolume_new[idx1];
              Up5 += s13 * weight*pvolume_new[idx1];
              Down += (weight*pvolume_new[idx1]);

             //Uptest += idx1 * weight*pvolume_new[idx1];

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

      pstress_new[idx](0,0) = Snonlocal[0];
      pstress_new[idx](1,2) = Snonlocal[1];
      pstress_new[idx](2,2) = Snonlocal[2];
      pstress_new[idx](0,1) = Snonlocal[3];
      pstress_new[idx](1,0) = Snonlocal[3];
      pstress_new[idx](1,2) = Snonlocal[4];
      pstress_new[idx](2,1) = Snonlocal[4];
      pstress_new[idx](2,0) = Snonlocal[5];
      pstress_new[idx](0,2) = Snonlocal[5];
      */

      c_dil = sqrt(USM / rho_cur);

      // Compute the strain energy for all the particles
      Matrix3 AvgStress = (pstress_new[idx] + pstress[idx]) * .5;

      double e = (D(0, 0) * AvgStress(0, 0) + D(1, 1) * AvgStress(1, 1) +
                  D(2, 2) * AvgStress(2, 2) +
                  2. * (D(0, 1) * AvgStress(0, 1) + D(0, 2) * AvgStress(0, 2) +
                        D(1, 2) * AvgStress(1, 2))) *
                 pvolume_new[idx] * delT;

      se += e;

      // Compute wave speed at each particle, store the maximum
      Vector pvelocity_idx = pvelocity[idx];
      WaveSpeed = Vector(Max(c_dil + fabs(pvelocity_idx.x()), WaveSpeed.x()),
                         Max(c_dil + fabs(pvelocity_idx.y()), WaveSpeed.y()),
                         Max(c_dil + fabs(pvelocity_idx.z()), WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z()) / 3.0;
        double c_bulk = sqrt(UI[1] / rho_cur);
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }
    } // end loop over particles

    WaveSpeed = dx / WaveSpeed;
    double delT_new = WaveSpeed.minComponent();
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se), lb->StrainEnergyLabel);
    }
  }
}

void
MohrCoulomb::carryForward(const PatchSubset* patches, const MPMMaterial* matl,
                          DataWarehouse* old_dw, DataWarehouse* new_dw)
{
  for (int p = 0; p < patches->size(); p++) {
    const Patch* patch = patches->get(p);
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    std::vector<constParticleVariable<double>> ISVs(d_NINSV + 1);
    std::vector<ParticleVariable<double>> ISVs_new(d_NINSV + 1);

    for (int i = 0; i < d_NINSV; i++) {
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

void
MohrCoulomb::addInitialComputesAndRequires(Task* task, const MPMMaterial* matl,
                                           const PatchSet*) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  cout << "In add InitialComputesAnd" << endl;

  // Other constitutive model and input dependent computes and requires
  for (int i = 0; i < d_NINSV; i++) {
    task->computes(ISVLabels[i], matlset);
  }
}

void
MohrCoulomb::addComputesAndRequires(Task* task, const MPMMaterial* matl,
                                    const PatchSet* patches) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);

  // Computes and requires for internal state data
  for (int i = 0; i < d_NINSV; i++) {
    task->requires(Task::OldDW, ISVLabels[i], matlset, Ghost::None);
    task->computes(ISVLabels_preReloc[i], matlset);
  }
}

void
MohrCoulomb::addComputesAndRequires(Task*, const MPMMaterial*, const PatchSet*,
                                    const bool) const
{
}

double
MohrCoulomb::computeRhoMicroCM(double pressure, const double p_ref,
                               const MPMMaterial* matl, double temperature,
                               double rho_guess)
{
  double rho_orig = matl->getInitialDensity();
  double p_gauge = pressure - p_ref;
  double rho_cur;
  double bulk = UI[1];

  rho_cur = rho_orig / (1 - p_gauge / bulk);

  return rho_cur;

#if 1
  cout << "NO VERSION OF computeRhoMicroCM EXISTS YET FOR MohrCoulomb" << endl;
#endif
}

void
MohrCoulomb::computePressEOSCM(double rho_cur, double& pressure, double p_ref,
                               double& dp_drho, double& tmp,
                               const MPMMaterial* matl, double temperature)
{

  double bulk = UI[1];
  double rho_orig = matl->getInitialDensity();

  double p_g = bulk * (1.0 - rho_orig / rho_cur);
  pressure = p_ref + p_g;
  dp_drho = bulk * rho_orig / (rho_cur * rho_cur);
  tmp = bulk / rho_cur; // speed of sound squared

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



/**
 ***********************************************************************
 * The procedure calculates stress for the Mohr-Coulomb model, with several
 * additional options. Those include the variation of the flow rule,
 * several choices of Mohr-Coulomb like yield surfaces and the dependency
 * on the strain rate. Most of the above-described features are 'work in
 * progress'
 *
 *
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
 */
void
MohrCoulomb::CalculateStress(int& nblk, int& ninsv, double& dt, double UI[],
                             double stress[], double D[], double svarg[],
                             double& USM, double shear_strain_nonlocal,
                             double shear_strain_rate_nonlocal)

{

  // this is slightly slow as each model used needs to be declared, but on the
  // other hand allows for keeping things clean

  ShengMohrCoulomb SMCModel;
  ClassicMohrCoulomb CMCModel;

  if (nblk != 1)
    cerr << "Mohr-Coulomb model may only be used with nblk equal to 1. Results "
            "obtained are incorrect."
         << endl;

  double G = UI[0]; // shear modulus [stress units]
  double K = UI[1]; // bulk modulus [stress units]
  double c = UI[2]; // cohesion [stress units]

  double m_modul = UI[22];
  double nuy = UI[23];
  double Phi = UI[3]; // friction angle [degrees]
  double Psi =
    UI[4]; // dilation angle [degrees, at the moment unused, assumed Phi]
  int Flavour = int(UI[5]);
  double a1 = UI[15];
  double b1 = UI[16];
  double W = UI[17];
  double beta = UI[18];
  double strain_ref = UI[19];

  double Use_softening = UI[34];
  double St = UI[35];
  double strain_95 = UI[36];

  double Use_linear = UI[25];
  double a = UI[26];
  double y_ref = UI[27];
  double y = svarg[37];

  /*
  Flavour
  1- classic Mohr - Coulomb,
  2 - classic Mohr-Coulomb with tension cut-off (not implemented yet)
  3 - ShengMohrCoulomb (rounded MC surface see eq 13 in Sheng D, Sloan SW & Yu
  HS
  Computations Mechanics 26:185-196 (2000)
  4 - ShengMohrCoulomb with tension cut-off (not implemented yet)
  5 - SloanMohrCoulomb (see Sloan et al... , not implemented yet)
  6 - SloanMohrCoulomb with tension cut-off
  11- classic Mohr - Coulomb with semi-implicit stress integration
  */

  int UseWaterRetention = UI[7];

  double Temp;

  // cerr<<UI[0]<<UI[1]<<' '<<UI[2]<<' '<<UI[3]<<' '<<UI[4];
  // this 2x3 lines are because the components in the added code are assumed in
  // different sequence.
  // probably this exchange is of no importance, but it is added for peace of
  // mind
  Temp = stress[4];
  stress[4] = stress[5];
  stress[5] = Temp;

  Temp = D[4];
  D[4] = D[5];
  D[5] = Temp;

  BBMPoint InitialPoint;
  double StrainIncrement[6];

  for (int i = 0; i < 6; i++) {
    InitialPoint.stress[i] = -stress[i];
    StrainIncrement[i] = -D[i] * dt;
  }

  double Usetransition = UI[14];
  if (Usetransition > 0) {

    if (shear_strain_rate_nonlocal > strain_ref) {
      c = St * a1 * pow(W, -b1) *
          pow(shear_strain_rate_nonlocal / strain_ref, beta);
    } else {
      c = St * a1 * pow(W, -b1);
    }
  }

  // Shear strength linear with depth
  if (Use_linear > 0) {
    c = c + a * (y - y_ref);
  }

  double Usemodul = UI[21];
  if (Usemodul > 0) {

    G = m_modul * c / 2.0 / (1.0 + nuy);
    K = m_modul * c / 3.0 / (1.0 - 2 * nuy);
  }

  if (Use_softening > 0) {
    if (shear_strain_nonlocal > c / G) {
      c = c * (1.0 / St +
               (1.0 - 1.0 / St) *
                 pow(2.71, (-3.0 * shear_strain_nonlocal / strain_95)));
    }
  }

  if (UseWaterRetention > 0) {
    double WTRParam[5];
    WTRParam[0] = UI[8];
    WTRParam[1] = UI[9];
    WTRParam[2] = UI[10];
    WTRParam[3] = UI[11];
    WTRParam[4] = 0.0;
    double PhiB = UI[13];

    double Suction;
    double SpecVol;
    Suction = svarg[6];
    SpecVol = svarg[12];

    double Sr = GetSr(UseWaterRetention, Suction,
                      WTRParam); // get from current suction value
    double dVolStrain =
      StrainIncrement[0] + StrainIncrement[1] + StrainIncrement[2];
    double SpecVolNew =
      SpecVol - dVolStrain * SpecVol; // 1/(1-TotalVolume);  //adjustment due to
    double VolWater =
      Sr * (SpecVol - 1) /
      SpecVol; // constant volume of water - thus we use initial specific volume
    double VolAir =
      (SpecVolNew - 1) / SpecVolNew -
      VolWater; // volume of voids - new one, using new specific volume
    double TotalVolume = VolAir + VolWater;
    double SrNew =
      VolWater /
      TotalVolume; // voids will reduce by dvol, so dSr=VolWater/TotalVolume
    if (SrNew > 0.99999)
      SrNew = 0.99999;
    if (SrNew < 0.0001)
      SrNew = 0.0001;
    double SuctionNew = GetSuction(
      UseWaterRetention, SrNew, WTRParam); // calculate from specified equation
                                           // (van Genuchten or Gallipoli or...)
    if (SuctionNew < 0)
      SuctionNew = 0.0;
    c = c + tan(PhiB * 3.1415 / 180) * SuctionNew;

    svarg[6] = SuctionNew;
    svarg[12] = SpecVolNew;
  }
  svarg[0] = G;
  svarg[1] = K;
  svarg[2] = c;

  int Region;

  switch (Flavour) {
    case 1: {
      CMCModel.SetModelParameters(G, K, c, Phi, Psi);
      CMCModel.SetDefaultIntegrationParameters();
      CMCModel.Integrate(StrainIncrement, &InitialPoint);

    } break;

    case 11: {
      CMCModel.SetModelParameters(G, K, c, Phi, Psi);
      CMCModel.SetDefaultIntegrationParameters();
      CMCModel.IntegrateMCIClassic(StrainIncrement, &InitialPoint, &Region);

    } break;

    case 12: {
      CMCModel.SetModelParameters(G, K, c, Phi, Psi);
      CMCModel.SetDefaultIntegrationParameters();
      CMCModel.IntegrateMCIClassic(StrainIncrement, &InitialPoint, &Region);

    } break;

    case 3: {
      SMCModel.SetModelParameters(G, K, c, Phi);
      SMCModel.SetDefaultIntegrationParameters();
      SMCModel.Integrate(StrainIncrement, &InitialPoint);
    } break;

    default: {
      cerr << "Error: Mohr Coulomb Model Flavour unspecified or the flavour "
              "demanded not implemented yet"
           << endl;
      cerr << "Flavour of the model is set to:" << Flavour << endl;
      getchar();
    }
  }

  // stress back for output
  for (int i = 0; i < 6; i++)
    stress[i] = -InitialPoint.stress[i];

  // this 2x3 lines are because the components in the added code are assumed in
  // different sequence.
  // probably this exchange is of no importance, but it is added for the peace
  // of mind

  Temp = stress[4];
  stress[4] = stress[5];
  stress[5] = Temp;

  Temp = D[4];
  D[4] = D[5];
  D[5] = Temp;

  double Factor = 5.0; // factor is  a quick fix, as otherwise the USM is too
                       // low and predicted stable time step is way too high

  USM = Factor * (G + 0.3 * K) / 3.0;
  // without Factor it
  // seems that the USM is too low and the analysis fails. Not sure why
  // as the elastic wave should be the quickest and it apparently work in diamm
  // maybe I missed something important there
}

double
MohrCoulomb::GetSr(double UseWaterRetention, double Suction, double* WTRParam)
{
  if (UseWaterRetention == 1.0) {
    // VanGenuchten Model
    double m = WTRParam[0];
    double n = WTRParam[1];
    double alpha = WTRParam[2];
    double Sr = alpha * Suction;
    Sr = pow(Sr, n);
    Sr = 1 / (1 + Sr);
    Sr = pow(Sr, m);
    if (Sr > 1.0)
      Sr = 1.0;
    if (Sr < 0)
      Sr = 0.0;
    return Sr;
  } else
    return 1.0;
}

double
MohrCoulomb::GetSuction(double UseWaterRetention, double Sr, double* WTRParam)
{
  if (UseWaterRetention == 1.0) {
    // VanGenuchten Model
    double m = WTRParam[0];
    double n = WTRParam[1];
    double alpha = WTRParam[2];
    double Suction = pow(Sr, 1 / m);
    Suction = (1 - Suction) / Suction;
    Suction = pow(Suction, 1 / n);
    Suction = Suction / alpha;
    if (Suction < 0)
      Suction = 0.0;
    return Suction;
  } else
    return 0;
}
