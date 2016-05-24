/*
 * The MIT License
 *
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

// Namespace Vaango::
#include <CCA/Components/MPM/ConstitutiveModel/Arenisca3PartiallySaturated.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/YieldConditionFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVariableModelFactory.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/KinematicHardeningModelFactory.h>

// Namespace Uintah::
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Malloc/Allocator.h>

#include <Core/Grid/Box.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Task.h>

#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/VarTypes.h>

#include <Core/Exceptions/ParameterNotFound.h>
#include <Core/Exceptions/InvalidValue.h>

#include <Core/Math/Matrix3.h>
#include <Core/Math/MinMax.h>
#include <Core/Math/MiscMath.h>

#include <sci_values.h>

// Namespace std::
#include <fstream>             
#include <iostream>
#include <limits>

#define WITH_BACKSTRESS

using namespace Vaango;
using SCIRun::VarLabel;
using Uintah::Matrix3;
using std::cout;

// constexpr auto M_PI = std::acos(-1.0);
// constexpr auto M_PIl = std::acosl(-1.0);

const double Arenisca3PartiallySaturated::one_third(1.0/3.0);
const double Arenisca3PartiallySaturated::two_third(2.0/3.0);
const double Arenisca3PartiallySaturated::four_third = 4.0/3.0;
const double Arenisca3PartiallySaturated::sqrt_two = std::sqrt(2.0);
const double Arenisca3PartiallySaturated::one_sqrt_two = 1.0/sqrt_two;
const double Arenisca3PartiallySaturated::sqrt_three = std::sqrt(3.0);
const double Arenisca3PartiallySaturated::one_sqrt_three = 1.0/sqrt_three;
const double Arenisca3PartiallySaturated::one_sixth = 1.0/6.0;
const double Arenisca3PartiallySaturated::one_ninth = 1.0/9.0;
const double Arenisca3PartiallySaturated::pi = M_PI;
const double Arenisca3PartiallySaturated::pi_fourth = 0.25*pi;
const double Arenisca3PartiallySaturated::pi_half = 0.5*pi;
const Matrix3 Arenisca3PartiallySaturated::Identity(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

// Set up constants for rotation 
const int Arenisca3PartiallySaturated::NMAX = 19;  // If this is changed, more entries will 
                                                   // need to be added to sinV cosV.

// Lookup tables for computing the sin() and cos() of the rotation angle.
const std::vector<double> Arenisca3PartiallySaturated::sinV =
                 {0.7071067811865475,-0.5,0.3420201433256687,-0.2306158707424402,0.1545187928078405,
                 -0.1032426220806015,0.06889665647555759,-0.04595133277786571,0.03064021661344469,
                 -0.02042858745187096,0.01361958465478159,-0.009079879062402308,0.006053298918749807,
                 -0.004035546304539714,0.002690368259933135,-0.001793580042002626,0.001195720384163988,
                 -0.0007971470283055577,0.0005314313834717263,-0.00035428759824575,0.0002361917349088998};
const std::vector<double> Arenisca3PartiallySaturated::cosV =
                {0.7071067811865475,0.8660254037844386,0.9396926207859084,0.9730448705798238,
                 0.987989849476809,0.9946562024066014,0.9976238022052647,0.9989436796015769,
                 0.9995304783376449,0.9997913146325693,0.999907249155556,0.9999587770484402,
                 0.9999816786182636,0.999991857149859,0.9999963809527642,0.9999983915340229,
                 0.9999992851261259,0.9999996822782572,0.9999998587903324,0.9999999372401469,
                 0.9999999721067318};

// Requires the necessary input parameters CONSTRUCTORS
Arenisca3PartiallySaturated::Arenisca3PartiallySaturated(Uintah::ProblemSpecP& ps, 
                                                         Uintah::MPMFlags* mpmFlags)
  : Uintah::ConstitutiveModel(mpmFlags)
{
  // Bulk and shear modulus models
  d_elastic = Vaango::ElasticModuliModelFactory::create(ps);
  if(!d_elastic){
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating ElasticModuliModel." << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Internal variable model
  d_intvar = Vaango::InternalVariableModelFactory::create(ps, d_elastic);
  if(!d_intvar){
    ostringstream desc;
    desc << "**ERROR** Internal error while creating InternalVariableModel." << endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Backstress model
  d_backstress = Vaango::KinematicHardeningModelFactory::create(ps, d_intvar);
  if(!d_backstress){
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating KinematicHardeningModel." << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Yield condition model
  d_yield = Vaango::YieldConditionFactory::create(ps, d_intvar);
  if(!d_yield){
    std::ostringstream desc;
    desc << "**ERROR** Internal error while creating YieldConditionModel." << std::endl;
    throw InternalError(desc.str(), __FILE__, __LINE__);
  }

  // Get initial porosity and saturation
  ps->require("initial_porosity",   d_fluidParam.phi0);  // Initial porosity
  ps->require("initial_saturation", d_fluidParam.Sw0);   // Initial water saturation

  ps->getWithDefault("subcycling_characteristic_number",
                     d_cm.subcycling_characteristic_number, 256);    // allowable subcycles
  ps->getWithDefault("use_disaggregation_algorithm",
                     d_cm.use_disaggregation_algorithm, false);

  /*
  d_ev0 = (std::abs(d_Kf) <= std::numeric_limits<double>::epsilon()*std::abs(d_Kf)) ?
    0.0 : d_C1*d_cm.fluid_pressure_initial/(d_Kf*d_Km); // Zero fluid pressure 
  */
  // vol. strain.  
  // (will equal zero if pfi=0)

  // MPM needs three functions to interact with ICE in MPMICE
  // 1) p = f(rho) 2) rho = g(p) 3) C = 1/K(rho)
  // Because the Arenisca3PartiallySaturated bulk modulus model does not have any closed
  // form expressions for these functions, we use a Murnaghan equation of state
  // with parameters K_0 and n = K_0'.  These parameters are read in here.
  // **WARNING** The default values are for Mason sand.
  ps->getWithDefault("K0_Murnaghan_EOS", d_cm.K0_Murnaghan_EOS, 2.5e8);
  ps->getWithDefault("n_Murnaghan_EOS", d_cm.n_Murnaghan_EOS, 13);

  checkInputParameters();

  initializeLocalMPMLabels();
}

void 
Arenisca3PartiallySaturated::checkInputParameters()
{
  
  if (d_cm.subcycling_characteristic_number < 1) {
    ostringstream warn;
    warn << "subcycling characteristic number should be > 1. Default = 256"<<endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }

  if (d_cm.use_disaggregation_algorithm) {
    ostringstream warn;
    warn << "Disaggregation algorithm not currently supported with partial saturation model"<<endl;
    throw ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

Arenisca3PartiallySaturated::Arenisca3PartiallySaturated(const Arenisca3PartiallySaturated* cm)
  : ConstitutiveModel(cm)
{
  d_elastic = Vaango::ElasticModuliModelFactory::createCopy(cm->d_elastic);
  d_yield   = Vaango::YieldConditionFactory::createCopy(cm->d_yield);
  d_intvar  = Vaango::InternalVariableModelFactory::createCopy(cm->d_intvar);
  d_backstress  = Vaango::KinematicHardeningModelFactory::createCopy(cm->d_backstress);

  // Porosity and saturation
  d_fluidParam = cm->d_fluidParam;

  // Subcycling
  d_cm.subcycling_characteristic_number = cm->d_cm.subcycling_characteristic_number;

  // Disaggregation Strain
  d_cm.use_disaggregation_algorithm = cm->d_cm.use_disaggregation_algorithm;

  // For MPMICE Murnaghan EOS
  d_cm.K0_Murnaghan_EOS = cm->d_cm.K0_Murnaghan_EOS;
  d_cm.n_Murnaghan_EOS = cm->d_cm.n_Murnaghan_EOS;

  initializeLocalMPMLabels();
}

// Initialize all labels of the particle variables associated with 
// Arenisca3PartiallySaturated.
void 
Arenisca3PartiallySaturated::initializeLocalMPMLabels()
{
  pPorosityLabel = Uintah::VarLabel::create("p.porosity",
        Uintah::ParticleVariable<double>::getTypeDescription());
  pPorosityLabel_preReloc = Uintah::VarLabel::create("p.porosity+",
        Uintah::ParticleVariable<double>::getTypeDescription());

  pSaturationLabel = Uintah::VarLabel::create("p.saturation",
        Uintah::ParticleVariable<double>::getTypeDescription());
  pSaturationLabel_preReloc = Uintah::VarLabel::create("p.saturation+",
        Uintah::ParticleVariable<double>::getTypeDescription());
      
  //pLocalized
  pLocalizedLabel = VarLabel::create("p.localized",
                                     ParticleVariable<int>::getTypeDescription());
  pLocalizedLabel_preReloc = VarLabel::create("p.localized+",
                                              ParticleVariable<int>::getTypeDescription());

  //pElasticVolStrain
  pElasticVolStrainLabel = VarLabel::create("p.elasticVolStrain",
    ParticleVariable<double>::getTypeDescription());
  pElasticVolStrainLabel_preReloc = VarLabel::create("p.elasticVolStrain+",
    ParticleVariable<double>::getTypeDescription());

  //pStressQS
  pStressQSLabel = VarLabel::create("p.stressQS",
                                    ParticleVariable<Matrix3>::getTypeDescription());
  pStressQSLabel_preReloc = VarLabel::create("p.stressQS+",
                                             ParticleVariable<Matrix3>::getTypeDescription());
}

// DESTRUCTOR
Arenisca3PartiallySaturated::~Arenisca3PartiallySaturated()
{
  VarLabel::destroy(pPorosityLabel);
  VarLabel::destroy(pPorosityLabel_preReloc);

  VarLabel::destroy(pSaturationLabel);
  VarLabel::destroy(pSaturationLabel_preReloc);

  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pLocalizedLabel_preReloc);
  VarLabel::destroy(pElasticVolStrainLabel);              //Elastic Volumetric Strain
  VarLabel::destroy(pElasticVolStrainLabel_preReloc);
  VarLabel::destroy(pStressQSLabel);
  VarLabel::destroy(pStressQSLabel_preReloc);

  delete d_yield;
  delete d_backstress;
  delete d_intvar;
  delete d_elastic;
}

//adds problem specification values to checkpoint data for restart
void 
Arenisca3PartiallySaturated::outputProblemSpec(ProblemSpecP& ps,bool output_cm_tag)
{
  ProblemSpecP cm_ps = ps;
  if (output_cm_tag) {
    cm_ps = ps->appendChild("constitutive_model");
    cm_ps->setAttribute("type","Arenisca3_part_sat");
  }

  d_elastic->outputProblemSpec(cm_ps);
  d_yield->outputProblemSpec(cm_ps);
  d_intvar->outputProblemSpec(cm_ps);
  d_backstress->outputProblemSpec(cm_ps);

  cm_ps->appendElement("initial_porosity",   d_fluidParam.phi0);
  cm_ps->appendElement("initial_saturation", d_fluidParam.Sw0);

  cm_ps->appendElement("subcycling_characteristic_number", d_cm.subcycling_characteristic_number);
  cm_ps->appendElement("use_disaggregation_algorithm",     d_cm.use_disaggregation_algorithm);

  // MPMICE Murnaghan EOS
  cm_ps->appendElement("K0_Murnaghan_EOS", d_cm.K0_Murnaghan_EOS);
  cm_ps->appendElement("n_Murnaghan_EOS",  d_cm.n_Murnaghan_EOS);
}

Arenisca3PartiallySaturated* 
Arenisca3PartiallySaturated::clone()
{
  return scinew Arenisca3PartiallySaturated(*this);
}

//When a particle is pushed from patch to patch, carry information needed for the particle
void 
Arenisca3PartiallySaturated::addParticleState(std::vector<const VarLabel*>& from,
                                              std::vector<const VarLabel*>& to)
{
  // Push back all the particle variables associated with Arenisca.
  // Important to keep from and to lists in same order!
  from.push_back(pPorosityLabel);
  to.push_back(pPorosityLabel_preReloc);

  from.push_back(pSaturationLabel);
  to.push_back(pSaturationLabel_preReloc);

  from.push_back(pLocalizedLabel);
  to.push_back(pLocalizedLabel_preReloc);

  from.push_back(pElasticVolStrainLabel);
  to.push_back(pElasticVolStrainLabel_preReloc);

  from.push_back(pStressQSLabel);
  to.push_back(pStressQSLabel_preReloc);

  // Add the particle state for the internal variable models
  d_intvar->addParticleState(from, to);

  // Add the particle state for the back stress model
  d_backstress->addParticleState(from, to);

  // Add the particle state for the yield condition model
  d_yield->addParticleState(from, to);

}

/*!------------------------------------------------------------------------*/
void 
Arenisca3PartiallySaturated::addInitialComputesAndRequires(Task* task,
                                                           const MPMMaterial* matl,
                                                           const PatchSet* patch) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();

  // Other constitutive model and input dependent computes and requires
  task->computes(pPorosityLabel,         matlset);
  task->computes(pSaturationLabel,       matlset);
  task->computes(pLocalizedLabel,        matlset);
  task->computes(pElasticVolStrainLabel, matlset);
  task->computes(pStressQSLabel,         matlset);

  // Add internal evolution variables computed by internal variable model
  d_intvar->addInitialComputesAndRequires(task, matl, patch);

  // Add internal evolution variables computed by the backstress model
  d_backstress->addInitialComputesAndRequires(task, matl, patch);

  // Add yield function variablity computes
  d_yield->addInitialComputesAndRequires(task, matl, patch);

}

/*!------------------------------------------------------------------------*/
void 
Arenisca3PartiallySaturated::initializeCMData(const Patch* patch,
                                              const MPMMaterial* matl,
                                              DataWarehouse* new_dw)
{
  // Add the initial porosity and saturation to the parameter dictionary
  ParameterDict allParams;
  allParams["phi0"] = d_fluidParam.phi0;
  allParams["Sw0"] = d_fluidParam.Sw0;

  // Get the particles in the current patch
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(),patch);

  // Get the particle volume and mass
  constParticleVariable<double> pVolume, pMass;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pMass,   lb->pMassLabel,   pset);

  // Initialize variables for yield function parameter variability
  d_yield->initializeLocalVariables(patch, pset, new_dw, pVolume);

  // Get yield condition parameters and add to the list of parameters
  ParameterDict yieldParams = d_yield->getParameters();
  allParams.insert(yieldParams.begin(), yieldParams.end());
  std::cout << "Model parameters are: " << std::endl;
  for (auto param : allParams) {
    std::cout << "\t \t" << param.first << " " << param.second << std::endl;
  }

  // Initialize variables for internal variables (needs yield function initialized first)
  d_intvar->initializeInternalVariable(patch, matl, pset, new_dw, lb, allParams);

  // Initialize variables for backstress model
  d_backstress->initializeLocalVariables(patch, pset, new_dw, pVolume);

  // Get backstress parameters and add to the list of parameters
  ParameterDict backstressParams = d_backstress->getParameters();
  allParams.insert(backstressParams.begin(), backstressParams.end());

  // Now initialize the other variables
  ParticleVariable<double>  pdTdt;
  ParticleVariable<Matrix3> pStress;
  ParticleVariable<double>  pPorosity, pSaturation;
  ParticleVariable<int>     pLocalized;
  ParticleVariable<double>  pElasticVolStrain; // Elastic Volumetric Strain
  ParticleVariable<Matrix3> pStressQS;

  new_dw->allocateAndPut(pdTdt,       lb->pdTdtLabel,               pset);
  new_dw->allocateAndPut(pStress,     lb->pStressLabel,             pset);

  new_dw->allocateAndPut(pPorosity,         pPorosityLabel,         pset);
  new_dw->allocateAndPut(pSaturation,       pSaturationLabel,       pset);
  new_dw->allocateAndPut(pLocalized,        pLocalizedLabel,        pset);
  new_dw->allocateAndPut(pElasticVolStrain, pElasticVolStrainLabel, pset);
  new_dw->allocateAndPut(pStressQS,         pStressQSLabel,         pset);

  // To fix : For a material that is initially stressed we need to
  // modify the stress tensors to comply with the initial stress state
  for(auto iter = pset->begin(); iter != pset->end(); iter++){
    pdTdt[*iter]             = 0.0;
    pStress[*iter]           = allParams["Pf0"]*Identity;
    pPorosity[*iter]         = d_fluidParam.phi0;
    pSaturation[*iter]       = d_fluidParam.Sw0;
    pLocalized[*iter]        = 0;
    pElasticVolStrain[*iter] = 0.0;
    pStressQS[*iter]         = pStress[*iter];
  }

  // Compute timestep
  computeStableTimestep(patch, matl, new_dw);
}

// Compute stable timestep based on both the particle velocities
// and wave speed
void 
Arenisca3PartiallySaturated::computeStableTimestep(const Patch* patch,
                                                   const MPMMaterial* matl,
                                                   DataWarehouse* new_dw)
{
  int matID = matl->getDWIndex();

  // Compute initial elastic moduli
  ElasticModuli moduli = d_elastic->getInitialElasticModuli();
  double bulk = moduli.bulkModulus;
  double shear = moduli.shearModulus;

  // Initialize wave speed
  double c_dil = std::numeric_limits<double>::min();
  Vector dx = patch->dCell();
  Vector WaveSpeed(c_dil, c_dil, c_dil);

  // Get the particles in the current patch
  ParticleSubset* pset = new_dw->getParticleSubset(matID, patch);

  // Get particles mass, volume, and velocity
  constParticleVariable<double> pMass, pVolume;
  constParticleVariable<long64> pParticleID;
  constParticleVariable<Vector> pVelocity;

  new_dw->get(pMass,       lb->pMassLabel,       pset);
  new_dw->get(pVolume,     lb->pVolumeLabel,     pset);
  new_dw->get(pParticleID, lb->pParticleIDLabel, pset);
  new_dw->get(pVelocity,   lb->pVelocityLabel,   pset);

  // loop over the particles in the patch
  for(auto iter = pset->begin(); iter != pset->end(); iter++){

    particleIndex idx = *iter;

    // Compute wave speed + particle velocity at each particle,
    // store the maximum
    c_dil = std::sqrt((bulk + four_third*shear)*(pVolume[idx]/pMass[idx]));

    //std::cout << "K = " << bulk << " G = " << shear << " c_dil = " << c_dil << std::endl;
    WaveSpeed = Vector(Max(c_dil+std::abs(pVelocity[idx].x()), WaveSpeed.x()),
                       Max(c_dil+std::abs(pVelocity[idx].y()), WaveSpeed.y()),
                       Max(c_dil+std::abs(pVelocity[idx].z()), WaveSpeed.z()));
  }

  // Compute the stable timestep based on maximum value of
  // "wave speed + particle velocity"
  WaveSpeed = dx/WaveSpeed;
  double delT_new = WaveSpeed.minComponent();
  new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
}

// ------------------------------------- BEGIN COMPUTE STRESS TENSOR FUNCTION
/**
 *  Arenisca3PartiallySaturated::computeStressTensor 
 *  is the core of the Arenisca3PartiallySaturated model which computes
 *  the updated stress at the end of the current timestep along with all other
 *  required data such plastic strain, elastic strain, cap position, etc.
 */
void 
Arenisca3PartiallySaturated::computeStressTensor(const PatchSubset* patches,
                                                 const MPMMaterial* matl,
                                                 DataWarehouse* old_dw,
                                                 DataWarehouse* new_dw)
{
  // Initial variables for internal variables
  ParameterDict yieldParams = d_yield->getParameters();

  // Initial variables for internal variable model
  ParameterDict intvarParam = d_intvar->getParameters();

  // Initial variables for backstress model
  ParameterDict backstressParam = d_backstress->getParameters();

  // Global loop over each patch
  for(int p=0;p<patches->size();p++){

    // Declare and initial value assignment for some variables
    const Patch* patch = patches->get(p);
    Matrix3 D(0.0);

    // Initialize wave speed
    double c_dil = std::numeric_limits<double>::min();
    Vector WaveSpeed(c_dil, c_dil, c_dil);
    Vector dx = patch->dCell();

    // Initialize strain energy
    double se = 0.0;  

    // Get particle subset for the current patch
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Copy yield condition parameter variability to new datawarehouse
    // Each particle has a different set of parameters which remain
    // constant through the simulation
    d_yield->copyLocalVariables(pset, old_dw, new_dw);

    // Get the yield condition parameter variables
    std::vector<constParticleVariable<double> > pYieldParamVars = 
      d_yield->getLocalVariables(pset, old_dw);

    // Get the internal variables
    std::vector<constParticleVariable<double> > pIntVarDouble = 
      d_intvar->getInternalVariables(pset, old_dw, 1.0);
    constParticleVariable<double> pKappa = pIntVarDouble[0];
    constParticleVariable<double> pCapX  = pIntVarDouble[1];
    constParticleVariable<double> pEpv   = pIntVarDouble[2];
    constParticleVariable<double> pP3    = pIntVarDouble[3];

    std::vector<constParticleVariable<Matrix3> > pIntVarMatrix3 = 
      d_intvar->getInternalVariables(pset, old_dw, Identity);
    constParticleVariable<Matrix3> pEp = pIntVarMatrix3[0];

    // Allocate and put internal variables
    // **WARNING** Hardcoded: (Need to know what those labels are right now)
    std::vector<ParticleVariable<double>* > pIntVarDouble_new; 
    ParticleVariable<double>  pKappa_new;
    ParticleVariable<double>  pCapX_new;
    ParticleVariable<double>  pEpv_new;
    ParticleVariable<double>  pP3_new;
    pIntVarDouble_new.push_back(&pKappa_new);
    pIntVarDouble_new.push_back(&pCapX_new);
    pIntVarDouble_new.push_back(&pEpv_new);
    pIntVarDouble_new.push_back(&pP3_new);
    d_intvar->allocateAndPutInternalVariable(pset, new_dw, pIntVarDouble_new);
    
    std::vector<ParticleVariable<Matrix3>* > pIntVarMatrix3_new; 
    ParticleVariable<Matrix3> pEp_new;
    pIntVarMatrix3_new.push_back(&pEp_new);
    d_intvar->allocateAndPutInternalVariable(pset, new_dw, pIntVarMatrix3_new);

    // Get the backstress matrix from the backstress model
    constParticleVariable<Matrix3> pBackStress;
    d_backstress->getBackStress(pset, old_dw, pBackStress);

    // Allocate and put the backstress matrix from the backstress model
    ParticleVariable<Matrix3> pBackStress_new;
    d_backstress->allocateAndPutBackStress(pset, new_dw, pBackStress_new);

    // Get the particle variables
    delt_vartype                   delT;
    constParticleVariable<int>     pLocalized;
    constParticleVariable<double>  pMass,           //used for stable timestep
                                   pElasticVolStrain;
    constParticleVariable<double>  pPorosity_old;
    constParticleVariable<double>  pSaturation_old;
    constParticleVariable<long64>  pParticleID;
    constParticleVariable<Vector>  pVelocity;
    constParticleVariable<Matrix3> pDefGrad,
                                   pStress_old, pStressQS_old;

    old_dw->get(delT,            lb->delTLabel,   getLevel(patches));
    old_dw->get(pMass,           lb->pMassLabel,               pset);
    old_dw->get(pParticleID,     lb->pParticleIDLabel,         pset);
    old_dw->get(pVelocity,       lb->pVelocityLabel,           pset);
    old_dw->get(pDefGrad,        lb->pDefGradLabel,            pset);
    old_dw->get(pStress_old,     lb->pStressLabel,             pset); 

    old_dw->get(pPorosity_old,     pPorosityLabel,               pset); 
    old_dw->get(pSaturation_old,   pSaturationLabel,             pset); 
    old_dw->get(pLocalized,        pLocalizedLabel,              pset); 
    old_dw->get(pElasticVolStrain, pElasticVolStrainLabel,       pset);
    old_dw->get(pStressQS_old,     pStressQSLabel,               pset);

    // Get the particle variables from interpolateToParticlesAndUpdate() in SerialMPM
    constParticleVariable<double>  pVolume;
    constParticleVariable<Matrix3> pVelGrad_new, pDefGrad_new;
    new_dw->get(pVolume,        lb->pVolumeLabel_preReloc,  pset);
    new_dw->get(pVelGrad_new,   lb->pVelGradLabel_preReloc, pset);
    new_dw->get(pDefGrad_new,   lb->pDefGradLabel_preReloc,      pset);

    // Get the particle variables from compute kinematics
    ParticleVariable<double>  p_q, pdTdt; 
    ParticleVariable<Matrix3> pStress_new;
    new_dw->allocateAndPut(p_q,                 lb->p_qLabel_preReloc,         pset);
    new_dw->allocateAndPut(pdTdt,               lb->pdTdtLabel_preReloc,       pset);
    new_dw->allocateAndPut(pStress_new,         lb->pStressLabel_preReloc,     pset);

    ParticleVariable<double>  pPorosity_new, pSaturation_new;
    ParticleVariable<int>     pLocalized_new;
    ParticleVariable<double>  pElasticVolStrain_new;
    ParticleVariable<Matrix3> pStressQS_new;
    new_dw->allocateAndPut(pPorosity_new,         pPorosityLabel_preReloc,         pset);
    new_dw->allocateAndPut(pSaturation_new,       pSaturationLabel_preReloc,       pset);
    new_dw->allocateAndPut(pLocalized_new,        pLocalizedLabel_preReloc,        pset);
    new_dw->allocateAndPut(pElasticVolStrain_new, pElasticVolStrainLabel_preReloc, pset);
    new_dw->allocateAndPut(pStressQS_new,         pStressQSLabel_preReloc,         pset);

    // Loop over the particles of the current patch to update particle
    // stress at the end of the current timestep along with all other
    // required data such plastic strain, elastic strain, cap position, etc.
    for (auto iter = pset->begin(); iter!=pset->end(); iter++) {
      particleIndex idx = *iter;  //patch index
      //cout<<"pID="<<pParticleID[idx]<<endl;

      // A parameter to consider the thermal effects of the plastic work which
      // is not coded in the current source code. Further development of Arenisca
      // may activate this feature.
      pdTdt[idx] = 0.0;

      // Copy particle deletion variable
      pLocalized_new[idx] = pLocalized[idx];

      // Compute the symmetric part of the velocity gradient
      //std::cout << "DefGrad = " << pDefGrad_new[idx] << std::endl;
      //std::cout << "VelGrad = " << pVelGrad_new[idx] << std::endl;
      Matrix3 D = (pVelGrad_new[idx] + pVelGrad_new[idx].Transpose())*.5;

      // Use polar decomposition to compute the rotation and stretch tensors
      Matrix3 FF = pDefGrad[idx];
      Matrix3 tensorR, tensorU;
      FF.polarDecompositionRMB(tensorU, tensorR);

      // Compute the unrotated symmetric part of the velocity gradient
      D = (tensorR.Transpose())*(D*tensorR);

      // To support non-linear elastic properties and to allow for the fluid bulk modulus
      // model to increase elastic stiffness under compression, we allow for the bulk
      // modulus to vary for each substep.  To compute the required number of substeps
      // we use a conservative value for the bulk modulus (the high pressure limit B0+B1)
      // to compute the trial stress and use this to subdivide the strain increment into
      // appropriately sized substeps.  The strain increment is a product of the strain
      // rate and time step, so we pass the strain rate and subdivided time step (rather
      // than a subdivided trial stress) to the substep function.

      // Compute the unrotated stress at the start of the current timestep
      Matrix3 sigma_old = (tensorR.Transpose())*(pStress_old[idx]*tensorR);
      Matrix3 sigmaQS_old = (tensorR.Transpose())*(pStressQS_old[idx]*tensorR);

      //std::cout << "pStress_old = " << pStress_old[idx] << std::endl
      //          << "pStressQS_old = " << pStressQS_old[idx] << std::endl;
      //std::cout << "sigma_old = " << sigma_old << std::endl
      //          << "sigmaQS_old = " << sigmaQS_old << std::endl;

      // initial assignment for the updated values of plastic strains, volumetric
      // part of the plastic strain, volumetric part of the elastic strain, kappa,
      // and the backstress. tentative assumption of elasticity
      ModelState_MasonSand state_old;
      state_old.capX                = pCapX[idx];
      state_old.kappa               = pKappa[idx];
      state_old.pbar_w              = -pBackStress[idx].Trace()/3.0;
      state_old.stressTensor        = sigmaQS_old;
      state_old.plasticStrainTensor = pEp[idx];
      state_old.p3                  = pP3[idx];
      state_old.porosity            = pPorosity_old[idx];
      state_old.saturation          = pSaturation_old[idx];

      //std::cout << "state_old.Stress = " << state_old.stressTensor << std::endl;


      // Get the parameters of the yield surface (for variability)
      for (auto& pYieldParamVar: pYieldParamVars) {
        state_old.yieldParams.push_back(pYieldParamVar[idx]);
      }

      // Compute the elastic moduli at t = t_n
      computeElasticProperties(state_old);
      //std::cout << "State old: " << state_old << std::endl;

      //---------------------------------------------------------
      // Rate-independent plastic step
      // Divides the strain increment into substeps, and calls substep function
      ModelState_MasonSand state_new;
      bool isSuccess = rateIndependentPlasticUpdate(D, delT, yieldParams,
                                                    idx, pParticleID[idx], state_old,
                                                    state_new);

      if (isSuccess) {

        pStressQS_new[idx] = state_new.stressTensor;     // unrotated stress at end of step
        pCapX_new[idx] = state_new.capX;                 // hydrostatic compressive strength at end of step
        pKappa_new[idx] = state_new.kappa;               // branch point
        pBackStress_new[idx] = Identity*(-state_new.pbar_w);  // trace of isotropic backstress at end of step
        pEp_new[idx] = state_new.plasticStrainTensor;    // plastic strain at end of step
        pEpv_new[idx] = pEp_new[idx].Trace();            // Plastic volumetric strain at end of step
        pP3_new[idx] = pP3[idx];

        // Elastic volumetric strain at end of step, compute from updated deformation gradient.
        pElasticVolStrain_new[idx] = log(pDefGrad_new[idx].Determinant()) - pEpv_new[idx];

        pPorosity_new[idx] = state_new.porosity;
        pSaturation_new[idx] = state_new.saturation;
      } else {

        // If the updateStressAndInternalVars function can't converge it will return false.  
        // This indicates substepping has failed, and the particle will be deleted.
        pLocalized_new[idx]=-999;
        cout << "** WARNING ** Bad step, deleting particle"
             << " idx = " << idx 
             << " particleID = " << pParticleID[idx] 
             << ":" << __FILE__ << ":" << __LINE__ << std::endl;
        pStressQS_new[idx] = pStressQS_old[idx];
        pCapX_new[idx] = state_old.capX; 
        pKappa_new[idx] = state_old.kappa;
        pBackStress_new[idx] = Identity*(-state_old.pbar_w);  // trace of isotropic backstress at end of step
        pEp_new[idx] = state_old.plasticStrainTensor;    // plastic strain at end of step
        pEpv_new[idx] = pEp_new[idx].Trace();
        pP3_new[idx] = pP3[idx];
        pElasticVolStrain_new[idx] = pElasticVolStrain[idx];
        pPorosity_new[idx] = pPorosity_old[idx];
        pSaturation_new[idx] = pSaturation_old[idx];
      }

      //---------------------------------------------------------
      // Rate-dependent plastic step
      ModelState_MasonSand stateQS_old(state_old);
      stateQS_old.stressTensor = pStressQS_old[idx];
      ModelState_MasonSand stateQS_new(state_new);
      stateQS_new.stressTensor = pStressQS_new[idx];
      
      //std::cout << "State QS old";
      computeElasticProperties(stateQS_old);
      //std::cout << "State QS new";
      computeElasticProperties(stateQS_new);
 
      rateDependentPlasticUpdate(D, delT, yieldParams, stateQS_old, stateQS_new, state_old,
                                 pStress_new[idx]);

      //---------------------------------------------------------
      // Update the porosity and saturation using the dynamic
      // stress state
      computePorosityAndSaturation(pStress_new[idx], state_old.pbar_w, backstressParam,
                                   pPorosity_new[idx], pSaturation_new[idx]);
      
      //---------------------------------------------------------
      // Use polar decomposition to compute the rotation and stretch tensors.  These checks prevent
      // failure of the polar decomposition algorithm if [F_new] has some extreme values.
      Matrix3 FF_new = pDefGrad_new[idx];
      double Fmax_new = FF_new.MaxAbsElem();
      double JJ_new = FF_new.Determinant();
      if ((Fmax_new > 1.0e16) || (JJ_new < 1.0e-16) || (JJ_new > 1.0e16)) {
        pLocalized_new[idx]=-999;
        std::cout << "Deformation gradient component unphysical: [F] = " << FF << std::endl;
        std::cout << "Resetting [F]=[I] for this step and deleting particle"
                  << " idx = " << idx 
                  << " particleID = " << pParticleID[idx] << std::endl;
        Identity.polarDecompositionRMB(tensorU, tensorR);
      } else {
        FF_new.polarDecompositionRMB(tensorU, tensorR);
      }

      // Compute the rotated dynamic and quasistatic stress at the end of the current timestep
      pStress_new[idx] = (tensorR*pStress_new[idx])*(tensorR.Transpose());
      pStressQS_new[idx] = (tensorR*pStressQS_new[idx])*(tensorR.Transpose());

      //std::cout << "pStress_new = " << pStress_new[idx]
      //          << "pStressQS_new = " << pStressQS_new[idx] << std::endl;

      // Compute wave speed + particle velocity at each particle, store the maximum
      //std::cout << "State QS new rotated";
      computeElasticProperties(stateQS_new); 
      double bulk = stateQS_new.bulkModulus;
      double shear = stateQS_new.shearModulus;
      double rho_cur = pMass[idx]/pVolume[idx];
      c_dil = sqrt((bulk+four_third*shear)/rho_cur);
      //std::cout << "K = " << bulk << " G = " << shear << " c_dil = " << c_dil << std::endl;
      WaveSpeed=Vector(Max(c_dil+std::abs(pVelocity[idx].x()),WaveSpeed.x()),
                       Max(c_dil+std::abs(pVelocity[idx].y()),WaveSpeed.y()),
                       Max(c_dil+std::abs(pVelocity[idx].z()),WaveSpeed.z()));

      // Compute artificial viscosity term
      if (flag->d_artificial_viscosity) {
        double dx_ave = (dx.x() + dx.y() + dx.z())*one_third;
        double c_bulk = sqrt(bulk/rho_cur);
        p_q[idx] = artificialBulkViscosity(D.Trace(), c_bulk, rho_cur, dx_ave);
      } else {
        p_q[idx] = 0.;
      }

      // Compute the averaged stress
      Matrix3 AvgStress = (pStress_new[idx] + pStress_old[idx])*0.5;

      // Compute the strain energy increment associated with the particle
      double e = (D(0,0)*AvgStress(0,0) +
                  D(1,1)*AvgStress(1,1) +
                  D(2,2)*AvgStress(2,2) +
                  2.0*(D(0,1)*AvgStress(0,1) +
                       D(0,2)*AvgStress(0,2) +
                       D(1,2)*AvgStress(1,2))) * pVolume[idx]*delT;

      // Accumulate the total strain energy
      // MH! Note the initialization of se needs to be fixed as it is currently reset to 0
      se += e;
    }

    // Compute the stable timestep based on maximum value of "wave speed + particle velocity"
    WaveSpeed = dx/WaveSpeed; // Variable now holds critical timestep (not speed)

    double delT_new = WaveSpeed.minComponent();

    // Put the stable timestep and total strain enrgy
    new_dw->put(delt_vartype(delT_new), lb->delTLabel, patch->getLevel());
    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(se),        lb->StrainEnergyLabel);
    }
  }
} // -----------------------------------END OF COMPUTE STRESS TENSOR FUNCTION

// ***************************************************************************************
// ***************************************************************************************
// **** HOMEL's FUNCTIONS FOR GENERALIZED RETURN AND NONLINEAR ELASTICITY ****************
// ***************************************************************************************
// ***************************************************************************************
/**
* Function: 
*   rateIndependentPlasticUpdate
*
* Purpose:
*   Divides the strain increment into substeps, and calls substep function
*   All stress values within computeStep are quasistatic.
*/
bool 
Arenisca3PartiallySaturated::rateIndependentPlasticUpdate(const Matrix3& D, 
                                                          const double& delT,
                                                          const ParameterDict& yieldParams,
                                                          particleIndex idx, 
                                                          long64 pParticleID, 
                                                          const ModelState_MasonSand& state_old,
                                                          ModelState_MasonSand& state_new)
{
  // Compute the trial stress
  Matrix3 strain_inc = D*delT;
  Matrix3 stress_trial = computeTrialStress(state_old, strain_inc);

  // Set up a trial state, update the stress invariants, and compute elastic properties
  ModelState_MasonSand state_trial(state_old);
  state_trial.stressTensor = stress_trial;
  computeElasticProperties(state_trial);
  std::cout << "Rate independent update:" << std::endl;
  std::cout << " D = " << D << " delT = " << delT << " strain_inc = " << strain_inc << std::endl;
  std::cout << "\t State old:" << state_old << std::endl;
  std::cout << "\t State trial:" << state_trial << std::endl;
  
  // Determine the number of substeps (nsub) based on the magnitude of
  // the trial stress increment relative to the characteristic dimensions
  // of the yield surface.  Also compare the value of the pressure dependent
  // elastic properties at sigma_old and sigma_trial and adjust nsub if
  // there is a large change to ensure an accurate solution for nonlinear
  // elasticity even with fully elastic loading.
  int nsub = computeStepDivisions(idx, pParticleID, state_old, state_trial, yieldParams);

  // * Upon FAILURE *
  // Delete the particle if the number of substeps is unreasonable
  // Send ParticleDelete Flag to Host Code, Store Inputs to particle data:
  // input values for sigma_new, X_new, Zeta_new, ep_new, along with error flag
  if (nsub < 0) {
    state_new = state_old;
    std::cout << "Step Failed: Particle idx = " << idx << " ID = " << pParticleID << std::endl;
    // bool success  = false;
    return false;
  }

  // Compute a subdivided time step:
  // Loop at least once or until substepping is successful
  const int CHI_MAX = 5;       // max allowed subcycle multiplier

  double dt = delT/nsub;       // substep time increment

  int chi = 1;                 // subcycle multiplier
  double tlocal = 0.0;
  bool isSuccess = false;

  // Set up the initial states for the substeps
  ModelState_MasonSand state_k_old(state_old);
  ModelState_MasonSand state_k_new(state_old);
  do {

    //  Call substep function {sigma_new, ep_new, X_new, Zeta_new}
    //    = computeSubstep(D, dt, sigma_substep, ep_substep, X_substep, Zeta_substep)
    //  Repeat while substeps continue to be successful
    isSuccess = computeSubstep(D, dt, yieldParams, state_k_old, state_k_new);
    if (isSuccess) {

      tlocal += dt;
      std::cout << "K = " << state_k_old.bulkModulus << std::endl;
      std::cout << "capX = " << state_k_new.capX << std::endl;
      std::cout << "pbar_w = " << state_k_new.pbar_w << std::endl;
      state_k_old = state_k_new;
      Matrix3 sig = state_k_new.stressTensor;
      std::cout << "sigma_new = np.array([[" 
                << sig(0,0) << "," << sig(0,1) << "," << sig(0,2) << "],[" 
                << sig(1,0) << "," << sig(1,1) << "," << sig(1,2) << "],[" 
                << sig(2,0) << "," << sig(2,1) << "," << sig(2,2) << "]])"
                << std::endl;
      std::cout << "plot_stress_state(K, G, sigma_trial, sigma_new, 'b')" << std::endl;

    } else {

      // Substepping has failed. Halve the timestep.
      dt /= 2.0;

      // Increase chi to keep track of the number of times the timstep has
      // been halved
      chi *= 2; 
      if (chi > CHI_MAX) {
        state_new = state_k_old;
        return isSuccess; // isSuccess = false;
      }

    }
    std::cout << "tlocal = " << tlocal << " delT = " << delT << " nsub = " << nsub << std::endl;
  } while (tlocal < delT);
    
  state_new = state_k_new;
  return isSuccess;

} 

/** 
 * Method: computeElasticProperties
 *
 * Purpose: 
 *   Compute the bulk and shear modulus at a given state
 */
void 
Arenisca3PartiallySaturated::computeElasticProperties(ModelState_MasonSand& state)
{
  state.updateStressInvariants();
  state.updateVolumetricPlasticStrain();
  ElasticModuli moduli = d_elastic->getCurrentElasticModuli(&state);
  state.bulkModulus = moduli.bulkModulus;
  state.shearModulus = moduli.shearModulus;

  // Modify the moduli if disaggregation is being used
  if (d_cm.use_disaggregation_algorithm) {
    double fac = std::exp(-(state.p3 + state.ep_v));
    double scale = std::max(fac, 0.00001);
    state.bulkModulus *= scale;
    state.shearModulus *= scale;
  }
}

/**
 * Method: computeTrialStress
 * Purpose: 
 *   Compute the trial stress for some increment in strain assuming linear elasticity
 *   over the step.
 */
Matrix3 
Arenisca3PartiallySaturated::computeTrialStress(const ModelState_MasonSand& state_old,
                                                const Matrix3& strain_inc)
{
  // Compute the trial stress
  Matrix3 stress_old = state_old.stressTensor;
  Matrix3 dEps_iso = Identity*(one_third*strain_inc.Trace());
  Matrix3 dEps_dev = strain_inc - dEps_iso;
  Matrix3 stress_trial = stress_old + 
                         dEps_iso*(3.0*state_old.bulkModulus) + 
                         dEps_dev*(2.0*state_old.shearModulus);

  return stress_trial;
} 

/**
 * Method: computeStepDivisions
 * Purpose: 
 *   Compute the number of step divisions (substeps) based on a comparison
 *   of the trial stress relative to the size of the yield surface, as well
 *   as change in elastic properties between sigma_n and sigma_trial.
 * 
 * Caveat:  Uses the mean values of the yield condition parameters.
 */
int 
Arenisca3PartiallySaturated::computeStepDivisions(particleIndex idx,
                                                  long64 particleID, 
                                                  const ModelState_MasonSand& state_old,
                                                  const ModelState_MasonSand& state_trial,
                                                  const ParameterDict& yieldParams)
{
  
  // Get the yield parameters
  double PEAKI1;
  double STREN;
  try {
    PEAKI1 = yieldParams.at("PEAKI1");
    STREN = yieldParams.at("STREN");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameters PEAKI1 and STREN" << std::endl;
    for (auto param : yieldParams) {
      err << param.first << " " << param.second << std::endl;
    }
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  
  // Compute change in bulk modulus:
  double bulk_old = state_old.bulkModulus;
  double bulk_trial = state_trial.bulkModulus;

  int n_bulk = std::ceil(std::abs(bulk_old - bulk_trial)/bulk_old);  
  std::cout << "bulk_old = " << bulk_old 
            << " bulk_trial = " << bulk_trial
            << " n_bulk = " << n_bulk << std::endl;
  
  // Compute trial stress increment relative to yield surface size:
  Matrix3 d_sigma = state_trial.stressTensor - state_old.stressTensor;
  double size = 0.5*(PEAKI1 - state_old.capX);
  if (STREN > 0.0){
    size = std::min(size, STREN);
  }  
  int n_yield = ceil(d_sigma.Norm()/size);

  std::cout << "PEAKI1 = " << PEAKI1 
            << " capX_old = " << state_old.capX
            << " size = " << size 
            << " |dsigma| = " << d_sigma.Norm() 
            << " n_yield = " << n_yield << std::endl;

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = std::max(n_bulk, n_yield);
  int nmax = d_cm.subcycling_characteristic_number;
 
  if (nsub > nmax) {
    std::cout << "\n **WARNING** Too many substeps needed for particle "
              << " idx = " << idx 
              << " particle ID = " << particleID << std::endl;
    std::cout << "\t" << __FILE__ << ":" << __LINE__ << std::endl;
    std::cout << "\t State at t_n: " << state_old;
    
    std::cout << "\t Trial state at t_n+1: " << state_trial;

    std::cout << "\t Ratio of trial bulk modulus to t_n bulk modulus "
              << n_bulk << std::endl;

    std::cout << "\t ||sig_trial - sigma_n|| " << d_sigma.Norm() << std::endl;
    std::cout << "\t Yield surface radius in I1-space: " << size << std::endl;
    std::cout << "\t Ratio of ||sig_trial - sigma_n|| and 10,000*y.s. radius: "
              << n_yield << std::endl;

    std::cout << "** BECAUSE** nsub = " << nsub << " > " 
              << d_cm.subcycling_characteristic_number
              << " : Probably too much tension in the particle."
              << std::endl;
    nsub = -1;
  }
  return nsub;
} 

/** 
 * Method: computeSubstep
 *
 * Purpose: 
 *   Computes the updated stress state for a substep that may be either 
 *   elastic, plastic, or partially elastic.   
 */
bool 
Arenisca3PartiallySaturated::computeSubstep(const Matrix3& D,
                                            const double& dt,
                                            const ParameterDict& yieldParams,
                                            const ModelState_MasonSand& state_k_old,
                                            ModelState_MasonSand& state_k_new)
{
  // Compute the trial stress
  Matrix3 deltaEps = D*dt;
  Matrix3 stress_k_trial = computeTrialStress(state_k_old, deltaEps);
  std::cout << "Inside computeSubstep:" << std::endl;
  std::cout << "K = " << state_k_old.bulkModulus << std::endl;
  std::cout << "capX = " << state_k_old.capX << std::endl;
  std::cout << "pbar_w = " << state_k_old.pbar_w << std::endl;
  Matrix3 sig = stress_k_trial;
  std::cout << "sigma_trial = np.array([[" 
            << sig(0,0) << "," << sig(0,1) << "," << sig(0,2) << "],[" 
            << sig(1,0) << "," << sig(1,1) << "," << sig(1,2) << "],[" 
            << sig(2,0) << "," << sig(2,1) << "," << sig(2,2) << "]])"
            << std::endl;
  std::cout << "plot_stress_state(K, G, sigma_new, sigma_trial, 'r')" << std::endl;
  //std::cout << "\t computeSubstep: sigma_old = " << state_k_old.stressTensor
  //         << " sigma_trial = " << stress_trial 
  //         << " D = " << D << " dt = " << dt
  //         << " deltaEps = " << deltaEps << std::endl;

  // Set up a trial state, update the stress invariants
  ModelState_MasonSand state_k_trial(state_k_old);
  state_k_trial.stressTensor = stress_k_trial;

  // Compute elastic moduli at trial stress state
  // and update stress invariants
  computeElasticProperties(state_k_trial);

  // Evaluate the yield function at the trial stress:
  int isElastic = (int) d_yield->evalYieldCondition(&state_k_trial); 

  // Elastic substep
  if (isElastic == 0 || isElastic == -1) { 
    state_k_new = state_k_trial;
    std::cout << "computeSubstep:Elastic:sigma_new = " << state_k_new.stressTensor
           << __FILE__ << ":" << __LINE__ << std::endl;
    return true; // bool isSuccess = true;
  }

  // Elastic-plastic or fully-plastic substep
  // Compute non-hardening return to initial yield surface:
  // returnFlag would be != 0 if there was an error in the nonHardeningReturn call, but
  // there are currently no tests in that function that could detect such an error.
  Matrix3 sig_fixed(0.0);        // final stress state for non-hardening return
  Matrix3 deltaEps_p_fixed(0.0); // increment in plastic strain for non-hardening return
  std::cout << "\t Doing nonHardeningReturn\n";
  nonHardeningReturn(deltaEps, state_k_old, state_k_trial, yieldParams,
                     sig_fixed, deltaEps_p_fixed);

  // Do "consistency bisection"
  state_k_new = state_k_old;
  std::cout << "\t Doing consistencyBisection\n";
  bool isSuccess = consistencyBisection(deltaEps, state_k_old, state_k_trial,
                                        deltaEps_p_fixed, sig_fixed, yieldParams, 
                                        state_k_new);

  return isSuccess;

} //===================================================================


/**
 * Method: nonHardeningReturn
 * Purpose: 
 *   Computes a non-hardening return to the yield surface in the meridional profile
 *   (constant Lode angle) based on the current values of the internal state variables
 *   and elastic properties.  Returns the updated stress and  the increment in plastic
 *   strain corresponding to this return.
 *
 *   NOTE: all values of r and z in this function are transformed!
 */
void 
Arenisca3PartiallySaturated::nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                                                const ModelState_MasonSand& state_k_old,
                                                const ModelState_MasonSand& state_k_trial,
                                                const ParameterDict& params,
                                                Uintah::Matrix3& sig_fixed,
                                                Uintah::Matrix3& plasticStrain_inc_fixed)
{
  // Get the yield parameters
  double BETA;
  double PEAKI1;
  try {
    BETA = params.at("BETA");
    PEAKI1 = params.at("PEAKI1");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameters BETA and PEAKI1" << std::endl;
    for (auto param : params) {
      err << param.first << " " << param.second << std::endl;
    }
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  // Compute ratio of bulk and shear moduli
  double K_old = state_k_old.bulkModulus;
  double G_old = state_k_old.shearModulus;
  const double sqrt_K_over_G_old = std::sqrt(1.5*K_old/G_old);

  // Save the r and z Lode coordinates for the trial stress state
  double r_trial = BETA*state_k_trial.rr;
  double z_eff_trial = state_k_trial.zz_eff;

  // Compute transformed r coordinates
  double rprime_trial = r_trial*sqrt_K_over_G_old;
  std::cout << " z_trial = " << z_eff_trial 
            << " r_trial = " << rprime_trial/sqrt_K_over_G_old << std::endl;

  // Find closest point
  double z_eff_closest = 0.0, rprime_closest = 0.0;
  bool foundClosestPt = d_yield->getClosestPoint(&state_k_old, z_eff_trial, rprime_trial, 
                                                 z_eff_closest, rprime_closest);
  if (!foundClosestPt) {
    std::ostringstream out;
    out << "**ERROR** Could not find the closest point to the yield surface.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  } 

  std::cout << " z_eff_closest = " << z_eff_closest 
            << " r_closest = " << rprime_closest/sqrt_K_over_G_old << std::endl;

  // Compute updated invariants of total stress
  double I1_closest = std::sqrt(3.0)*z_eff_closest - 3.0*state_k_old.pbar_w;
  double sqrtJ2_closest = 1.0/(sqrt_K_over_G_old*BETA*sqrt_two)*rprime_closest;
  //std::cout << "I1_closest = " << I1_closest 
  //          << " sqrtJ2_closest = " << sqrtJ2_closest << std::endl;
  //std::cout << "Trial state = " << state_trial << std::endl;

  // Compute new stress
  Matrix3 sig_dev = state_k_trial.deviatoricStressTensor;
  if (state_k_trial.sqrt_J2 > 0.0) {
    sig_fixed = one_third*I1_closest*Identity + 
     (sqrtJ2_closest/state_k_trial.sqrt_J2)*sig_dev;
  } else {
    sig_fixed = one_third*I1_closest*Identity + sig_dev;
  }

  // Compute new plastic strain increment
  //  d_ep = d_e - [C]^-1:(sigma_new-sigma_old)
  Matrix3 sig_inc = sig_fixed - state_k_old.stressTensor;
  Matrix3 sig_inc_iso = one_third*sig_inc.Trace()*Identity;
  Matrix3 sig_inc_dev = sig_inc - sig_inc_iso;
  Matrix3 elasticStrain_inc = sig_inc_iso*(one_third/K_old) + sig_inc_dev*(0.5/G_old);
  plasticStrain_inc_fixed = strain_inc - elasticStrain_inc;
  //std::cout << "\t\t\t sig_inc = " << sig_inc << std::endl;
  //std::cout << "\t\t\t strain_inc = " << strain_inc << std::endl;
  //std::cout << "\t\t\t sig_inc_iso = " << sig_inc_iso << std::endl;
  //std::cout << "\t\t\t sig_inc_dev = " << sig_inc_dev << std::endl;
  //std::cout << "\t\t\t plasticStrain_inc_fixed = " << plasticStrain_inc_fixed << std::endl;

} //===================================================================

/**
 * Method: evalYieldCondition
 * Purpose: 
 *   Evaluate the yield condition in transformed Lode coordinates
 *   Returns whether the stress is elastic or not
 */
bool
Arenisca3PartiallySaturated::evalYieldCondition(const double& z_eff_stress, const double& rprime_stress, 
                                                const ModelState_MasonSand& state_old, 
                                                const ParameterDict& params)
{
  // Compute untransformed invariants
  double BETA;
  try {
    BETA = params.at("BETA");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameter BETA" << std::endl;
    for (auto param : params) {
      err << param.first << " " << param.second << std::endl;
    }
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  double G_over_K = std::sqrt(state_old.shearModulus/(1.5*state_old.bulkModulus));
  double I1_eff_stress = std::sqrt(3.0)*z_eff_stress;
  double sqrtJ2_stress = G_over_K*(1.0/(std::sqrt(2.0)*BETA))*rprime_stress;

  // Create a temporary state for evaluation the yield function
  ModelState_MasonSand state_stress(state_old);
  state_stress.I1_eff = I1_eff_stress;
  state_stress.sqrt_J2 = sqrtJ2_stress;

  // Evaluate the yield function
  int isElastic = (int) d_yield->evalYieldCondition(&state_stress);
  if (isElastic == 1) {
    return false;  // Plastic
  }
  
  return true; // Elastic or on yield surface
}

/**
 * Method: consistencyBisection
 * Purpose: 
 *   Find the updated stress for hardening plasticity using the consistency bisection 
 *   algorithm
 *   Returns whether the procedure is sucessful or has failed
 */
bool 
Arenisca3PartiallySaturated::consistencyBisection(const Matrix3& deltaEps_new,
                                                  const ModelState_MasonSand& state_old, 
                                                  const ModelState_MasonSand& state_trial,
                                                  const Matrix3& deltaEps_p_0, 
                                                  const Matrix3& sig_0, 
                                                  const ParameterDict& params, 
                                                  ModelState_MasonSand& state_new)
{
  const double TOLERANCE = 1e-4; // bisection convergence tolerance on eta (if changed, change imax)
  //const double CAPX_TOLERANCE = 1e-6; 
  //const double ZETA_TOLERANCE = 1e-6; 
  const int    IMAX      = 93;   // imax = ceil(-10.0*log(TOL)); // Update this if TOL changes
  const int    JMAX      = 93;   // jmax = ceil(-10.0*log(TOL)); // Update this if TOL changes

  // Get the initial fluid pressure
  ParameterDict backstressParams = d_backstress->getParameters();

  // Local trial stress
  ModelState_MasonSand state_trial_upd(state_trial);
  state_trial_upd.updateStressInvariants();

  // Initialize
  Matrix3 sig_old     = state_old.stressTensor;
  Matrix3 eps_p_old   = state_old.plasticStrainTensor;
  double  eps_p_v_old = eps_p_old.Trace();
  double  pbar_w_old    = state_old.pbar_w;

  Matrix3 sig_new          = sig_0;
  Matrix3 deltaEps_p_new   = deltaEps_p_0;
  double  deltaEps_p_v_0 = deltaEps_p_0.Trace();
  double  deltaEps_p_v_new = deltaEps_p_v_0;

  // Start loop
  int ii = 1;
  double eta_in = 0.0, eta_out = 1.0, eta_mid = 0.5;
  double norm_deltaEps_p_0 = deltaEps_p_0.Norm();
  double norm_deltaEps_p_new = norm_deltaEps_p_0;
  double capX_0 = state_old.capX;
  double capX_new = capX_0;
  double pbar_w_0 = state_old.pbar_w;
  double pbar_w_new = pbar_w_0;

  do { // While (eta_in - eta_out) > TOL)

    int jj = 1;
    bool isElastic = true;

    do { // while isElastic

      eta_mid = 0.5*(eta_in + eta_out); 
      double deltaEps_p_v_mid = eta_mid*deltaEps_p_v_new;
      double eps_p_v_mid = eps_p_v_old + deltaEps_p_v_mid;

      // Update hydrostatic compressive strength
      state_trial_upd.ep_v = eps_p_v_mid;
      state_trial_upd.p3 = state_old.p3;
      capX_0 = capX_new;
      capX_new = computeHydrostaticStrength(state_trial_upd);  
      //std::cout << "\t\t ii = " << ii << " jj = " << jj << std::endl;

      // Update the isotropic backstress
      state_trial_upd.dep_v = deltaEps_p_v_mid;
      state_trial_upd.phi0  = d_fluidParam.phi0;
      state_trial_upd.Sw0   = d_fluidParam.Sw0;
      Matrix3 backStress_new;
      d_backstress->computeBackStress(&state_trial_upd, backStress_new);
      pbar_w_0 = pbar_w_new;

#ifdef WITH_BACKSTRESS
      pbar_w_new = -backStress_new.Trace()/3.0;
#else
      pbar_w_new = 0.0;
#endif

      std::cout << "\t\t " << "eta_in = " << eta_in << " eta_mid = " << eta_mid
                << " eta_out = " << eta_out
                << " capX_0 = " << capX_0 << " capX_new = " << capX_new 
                << " pbar_w_0 = " << pbar_w_0 << " pbar_w_new = " << pbar_w_new 
                << " ||delta eps_p_0|| = " << norm_deltaEps_p_0 
                << " ||delta eps_p_new|| = " << norm_deltaEps_p_new << std::endl;

      // Update the trial stress
      state_trial_upd.capX = capX_new;
      state_trial_upd.pbar_w = pbar_w_new;

      // Test the yield condition
      int yield = (int) d_yield->evalYieldCondition(&state_trial_upd);
      if (yield == 1) {
        isElastic = false;  // Plastic
      } else {
        isElastic = true;   // Elastic or on yield surface
      }

      // If the state is elastic, there is too much plastic strain and
      // the following will be used;
      // otherwise control will break out from the isElastic loop
      eta_out = eta_mid;
      jj++;

      if (std::abs(eta_in - eta_out) < TOLERANCE)  {
        std::cout << "In consistency bisection: The mid_point is equal to the start point" << std::endl;
        break;
      }

      // Too many iterations
      if (jj > JMAX) {
        state_new = state_old;
        // bool isSuccess = false;
        return false;
      }
      
    } while(isElastic);
    //        && (std::abs((capX_0 - capX_new)/capX_new) > CAPX_TOLERANCE) &&
    //        (std::abs((pbar_w_0 - pbar_w_new)/pbar_w_new) > ZETA_TOLERANCE));

    // Update the state and compute elastic properties
    Matrix3 sig_mid = (sig_old + sig_new)*0.5;
    Matrix3 eps_p_mid = eps_p_old + deltaEps_p_0*(eta_mid*0.5);
    double p_bar_w_mid = (pbar_w_old + pbar_w_new)*0.5;

    double phi_mid = 0.0, Sw_mid = 0.0;
    computePorosityAndSaturation(sig_mid, p_bar_w_mid, backstressParams, phi_mid, Sw_mid);

    ModelState_MasonSand state_mid;
    state_mid.stressTensor = sig_mid;
    state_mid.plasticStrainTensor = eps_p_mid;
    state_mid.porosity = phi_mid;
    state_mid.saturation = Sw_mid;
    computeElasticProperties(state_mid);

    // Do non hardening return with updated properties
    state_trial_upd.bulkModulus = state_mid.bulkModulus;
    state_trial_upd.shearModulus = state_mid.shearModulus;
    state_trial_upd.capX = capX_new;
    state_trial_upd.pbar_w = std::max(pbar_w_new, 0.0);
    state_trial_upd.porosity = phi_mid;
    state_trial_upd.saturation = Sw_mid;
    state_trial_upd.updateStressInvariants();
    nonHardeningReturn(deltaEps_new, state_old, state_trial_upd, params,
                       sig_new, deltaEps_p_new);

    // Set up variables for various tests
    Matrix3 sig_trial = state_trial_upd.stressTensor;
    double trial_new = (sig_trial - sig_new).Trace();
    double trial_0 = (sig_trial - sig_0).Trace();
    norm_deltaEps_p_new = deltaEps_p_new.Norm();
    norm_deltaEps_p_0 = eta_mid*deltaEps_p_0.Norm();

    //std::cout << "\t\t K = " << state_trial_upd.bulkModulus
    //          << " G = " << state_trial_upd.shearModulus
    //          << " capX = " << state_trial_upd.capX
    //          << " pbar_w = " << state_trial_upd.pbar_w
    //          << " phi = " << state_trial_upd.porosity
    //          << " Sw = " << state_trial_upd.saturation << std::endl;
    //std::cout << "\t\t sig_0 = " << sig_0 << std::endl;
    //std::cout << "\t\t sig_trial = " << sig_trial << std::endl;
    //std::cout << "\t\t sig_new = " << sig_new << std::endl;
    //std::cout << "\t\t trial_0 = " << trial_0 
    //          << " trial_new = " << trial_new 
    //          << " ||deltaEps_p_0|| = " << norm_deltaEps_p_0
    //          << " ||deltaEps_p_new|| = " << norm_deltaEps_p_new << std::endl;

    // Check whether the isotropic component of the return has changed sign, as this
    // would indicate that the cap apex has moved past the trial stress, indicating
    // too much plastic strain in the return.
    if (std::signbit(trial_new) != std::signbit(trial_0)) {
      eta_out = eta_mid;
      continue;
    }

    // Compare magnitude of plastic strain with prior update
    if (norm_deltaEps_p_new > norm_deltaEps_p_0) {
      eta_in = eta_mid;
    } else {
      eta_out = eta_mid;
    }

    // Increment i and check
    ii++;
    if (ii > IMAX) {
      state_new = state_old;
      // bool isSuccess = false;
      return false;
    }

    //std::cout << "\t\t\t " << "eta_in = " << eta_in << " eta_mid = " << eta_mid
    //          << " eta_out = " << eta_out << std::endl;

  } while (std::abs(eta_out - eta_in) > TOLERANCE);

  // Update the plastic strain
  Matrix3 eps_p_new = eps_p_old + deltaEps_p_new;

  // Update hydrostatic compressive strength
  state_trial_upd.ep_v = eps_p_new.Trace();
  state_trial_upd.p3 = state_old.p3;
  capX_new = computeHydrostaticStrength(state_trial_upd);  

  // Update the isotropic backstress
  state_trial_upd.dep_v = deltaEps_p_new.Trace();
  state_trial_upd.phi0  = d_fluidParam.phi0;
  state_trial_upd.Sw0   = d_fluidParam.Sw0;
  Matrix3 backStress_new;
  d_backstress->computeBackStress(&state_trial_upd, backStress_new);
#ifdef WITH_BACKSTRESS
  pbar_w_new = -backStress_new.Trace()/3.0;
#else
  pbar_w_new = 0.0;
#endif

  // Update the state
  state_new = state_trial_upd;  
  state_new.stressTensor = sig_new;
  state_new.plasticStrainTensor = eps_p_new;
  state_new.capX = capX_new;

  // max() eliminates tensile fluid pressure from explicit integration error
  state_new.pbar_w = std::max(state_new.pbar_w, 0.0);

  // Return success = true  
  // bool isSuccess = true;
  return true;
}

/** 
 * Method: computeHydrostaticStrength
 * Purpose: 
 *   Compute state variable X, the Hydrostatic Compressive strength (cap position)
 *   X is the value of (I1 - Zeta) at which the cap function crosses
 *   the hydrostat. 
 *   In tension, M. Homel's piecewise formulation is used.
 */
double 
Arenisca3PartiallySaturated::computeHydrostaticStrength(const ModelState_MasonSand& state)
{
  double capX = d_intvar->computeInternalVariable(&state);
  return capX;

} //===================================================================

/** 
 * Function: rateDependentPlasticUpdate
 *
 * Purpose:
 *   Rate-dependent plastic step
 *   Compute the new dynamic stress from the old dynamic stress and the new and old QS stress
 *   using Duvaut-Lions rate dependence, as described in "Elements of Phenomenological Plasticity",
 *   by RM Brannon.
 */
bool 
Arenisca3PartiallySaturated::rateDependentPlasticUpdate(const Matrix3& D,
                                                        const double& delT,
                                                        const ParameterDict& yieldParams,
                                                        const ModelState_MasonSand& stateStatic_old,
                                                        const ModelState_MasonSand& stateStatic_new,
                                                        const ModelState_MasonSand& stateDynamic_old,
                                                        Matrix3& pStress_new) 
{
  // Get the T1 & T2 parameters
  double T1 = 0.0, T2 = 0.0;
  try {
    T1 = yieldParams.at("T1");
    T2 = yieldParams.at("T2");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameters T1 and T2" << std::endl;
    for (auto param : yieldParams) {
      err << param.first << " " << param.second << std::endl;
    }
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  // Check if rate-dependent plasticity has been turned on
  if (T1 == 0.0 || T2 == 0.0) {

    // No rate dependence, the dynamic stress equals the static stress.
    pStress_new = stateStatic_new.stressTensor;
    // bool isRateDependent = false;
    return false;

  }

  // This is not straightforward, due to nonlinear elasticity.  The equation requires that we
  // compute the trial stress for the step, but this is not known, since the bulk modulus is
  // evolving through the substeps.  It would be necessary to to loop through the substeps to
  // compute the trial stress assuming nonlinear elasticity, but instead we will approximate
  // the trial stress the average of the elastic moduli at the start and end of the step.

  // Compute midstep bulk and shear modulus
  ModelState_MasonSand stateDynamic(stateDynamic_old);
  stateDynamic.bulkModulus = 0.5*(stateStatic_old.bulkModulus + stateStatic_new.bulkModulus);
  stateDynamic.shearModulus = 0.5*(stateStatic_old.shearModulus + stateStatic_new.shearModulus);

  Matrix3 strain_inc = D*delT;
  Matrix3 sigma_trial = computeTrialStress(stateDynamic, strain_inc);

  // The characteristic time is defined from the rate dependence input parameters and the
  // magnitude of the strain rate.
  // tau = T1*(epsdot)^(-T2) = T1*(1/epsdot)^T2, modified to avoid division by zero.
  double tau = T1*std::pow(1.0/std::max(D.Norm(), 1.0e-15), T2);

  // RH and rh are defined by eq. 6.93 in the RMB book chapter, but there seems to be a sign error
  // in the text, and I've rewritten it to avoid computing the exponential twice.
  double dtbytau = delT/tau;
  double rh  = std::exp(-dtbytau);
  double RH  = (1.0 - rh)/dtbytau;

  // sigma_new = sigmaQS_new + sigma_over_new, as defined by eq. 6.92
  // sigma_over_new = [(sigma_trial_new - sigma_old) - (sigmaQS_new-sigmaQS_old)]*RH + sigma_over_old*rh
  Matrix3 sigmaQS_old = stateStatic_old.stressTensor;
  Matrix3 sigmaQS_new = stateStatic_new.stressTensor;
  Matrix3 sigma_old = stateDynamic_old.stressTensor;
  pStress_new = sigmaQS_new
          + ((sigma_trial - sigma_old) - (sigmaQS_new - sigmaQS_old))*RH
          + (sigma_old - sigmaQS_old)*rh;

  // bool isRateDependent = true;
  return false;
}

/**
 * Method: computePorosityAndSaturation
 *
 * Purpose: 
 *   Compute porosity (phi) and saturation (S_w)
 *
 * TODO:
 *   Don't recompute exponentials of strains
 *
 */
void
Arenisca3PartiallySaturated::computePorosityAndSaturation(const Matrix3& stress,
                                                          const double& pbar_w,
                                                          const ParameterDict& params,
                                                          double& porosity,
                                                          double& saturation)
{
  double pf0 = 0.0;
  try {
    pf0 = params.at("Pf0");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find parameter Pf0" << std::endl;
    for (auto param : params) {
      err << param.first << " " << param.second << std::endl;
    }
    throw InternalError(err.str(), __FILE__, __LINE__);
  }

  double I1_bar = -stress.Trace();
  double phi0 = d_fluidParam.phi0;
  double Sw0 = d_fluidParam.Sw0;

  saturation = computeSaturation(pbar_w, pf0, Sw0);
  porosity = computePorosity(I1_bar, pbar_w, pf0, phi0, Sw0);
}

/**
 * Method: computePorosity
 *
 * Purpose: 
 *   Compute porosity (phi)
 *
 * TODO:
 *   Compute porosity and saturation with one function call
 *
 */
double 
Arenisca3PartiallySaturated::computePorosity(const double& I1_bar,
                                             const double& pbar_w,
                                             const double& pf0,
                                             const double& phi0,
                                             const double& Sw0)
{
  // Compute air and water volume strain
  double exp_ev_a = d_air.computeExpElasticVolumetricStrain(3.0*pbar_w, 0.0);
  double exp_ev_w = d_water.computeExpElasticVolumetricStrain(3.0*pbar_w, pf0);

  // Compute total volumetric strain
  double ev = computeTotalVolStrain(I1_bar, pbar_w, pf0, phi0, Sw0);

  // Compute saturation evolution
  double S_w = computeSaturation(pbar_w, pf0, Sw0);

  // Compute porosity
  double phi = phi0;
  if (Sw0 == 1.0) {
    phi = phi0*exp_ev_w*std::exp(-ev);
  } else {
    phi = phi0*(1-Sw0)/(1-S_w)*exp_ev_a*std::exp(-ev);
  }

  return (phi);
}

/**
 * Method: computeSaturation
 *
 * Purpose: 
 *   Compute water saturation (Sw)
 *
 */
double 
Arenisca3PartiallySaturated::computeSaturation(const double& pbar_w,
                                               const double& pf0,
                                               const double& Sw0)
{
  // Compute air and water volume strain
  double exp_ev_a = d_air.computeExpElasticVolumetricStrain(3.0*pbar_w, 0.0);
  double exp_ev_w = d_water.computeExpElasticVolumetricStrain(3.0*pbar_w, pf0);

  double S_w = Sw0;

  if (Sw0 > 0.0 && Sw0 < 1.0) {

    // Compute C
    double C = Sw0/(1.0 - Sw0)*exp_ev_w/exp_ev_a;

    // Compute saturation
    S_w = C/(1.0 + C);
  }

  return (S_w);

}

/**
 * Method: computeTotalVolStrain
 *
 * Purpose: 
 *   Compute the total volumetric strain (compression positive)
 *
 */
double 
Arenisca3PartiallySaturated::computeTotalVolStrain(const double& I1_bar,
                                                   const double& pbar_w,
                                                   const double& pf0,
                                                   const double& phi0,
                                                   const double& Sw0)
{
  // Compute volume strains in the three components
  double exp_ev_a = d_air.computeExpElasticVolumetricStrain(3.0*pbar_w, 0.0);
  double exp_ev_w = d_water.computeExpElasticVolumetricStrain(3.0*pbar_w, pf0);
  double exp_ev_s = d_solid.computeExpElasticVolumetricStrain(I1_bar - 3.0*pbar_w, 0.0);

  // Compute total vol strain
  double exp_ev = (1.0 - Sw0)*phi0*exp_ev_a + Sw0*phi0*exp_ev_w + (1.0 - phi0)*exp_ev_s;
  double ev = std::log(exp_ev);

  return (-ev);
}


// ****************************************************************************************************
// ****************************************************************************************************
// ************** PUBLIC Uintah MPM constitutive model specific functions *****************************
// ****************************************************************************************************
// ****************************************************************************************************

void Arenisca3PartiallySaturated::addRequiresDamageParameter(Task* task,
                                                             const MPMMaterial* matl,
                                                             const PatchSet* ) const
{
  // Require the damage parameter
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc,matlset,Ghost::None);
}

void Arenisca3PartiallySaturated::getDamageParameter(const Patch* patch,
                                                     ParticleVariable<int>& damage,
                                                     int matID,
                                                     DataWarehouse* old_dw,
                                                     DataWarehouse* new_dw)
{
  // Get the damage parameter
  ParticleSubset* pset = old_dw->getParticleSubset(matID,patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, pLocalizedLabel_preReloc, pset);

  // Loop over the particle in the current patch.
  for (auto iter = pset->begin(); iter != pset->end(); iter++) {
    damage[*iter] = pLocalized[*iter];
  }
}

void Arenisca3PartiallySaturated::carryForward(const PatchSubset* patches,
                                               const MPMMaterial* matl,
                                               DataWarehouse* old_dw,
                                               DataWarehouse* new_dw)
{
  // Carry forward the data.
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    int matID = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(matID, patch);

    // Carry forward the data common to all constitutive models
    // when using RigidMPM.
    // This method is defined in the ConstitutiveModel base class.
    carryForwardSharedData(pset, old_dw, new_dw, matl);

    // Carry forward the data local to this constitutive model
    new_dw->put(delt_vartype(1.e10), lb->delTLabel, patch->getLevel());

    if (flag->d_reductionVars->accStrainEnergy ||
        flag->d_reductionVars->strainEnergy) {
      new_dw->put(sum_vartype(0.0),     lb->StrainEnergyLabel);
    }
  }
}


void Arenisca3PartiallySaturated::addComputesAndRequires(Task* task,
                                                         const MPMMaterial* matl,
                                                         const PatchSet* patches ) const
{
  // Add the computes and requires that are common to all explicit
  // constitutive models.  The method is defined in the ConstitutiveModel
  // base class.
  const MaterialSubset* matlset = matl->thisMaterial();
  addSharedCRForHypoExplicit(task, matlset, patches);
  task->requires(Task::OldDW, lb->pParticleIDLabel,   matlset, Ghost::None);
  task->requires(Task::OldDW, pPorosityLabel,         matlset, Ghost::None);
  task->requires(Task::OldDW, pSaturationLabel,       matlset, Ghost::None);
  task->requires(Task::OldDW, pLocalizedLabel,        matlset, Ghost::None);
  task->requires(Task::OldDW, pElasticVolStrainLabel, matlset, Ghost::None);
  task->requires(Task::OldDW, pStressQSLabel,         matlset, Ghost::None);
  task->computes(pPorosityLabel_preReloc,         matlset);
  task->computes(pSaturationLabel_preReloc,       matlset);
  task->computes(pLocalizedLabel_preReloc,        matlset);
  task->computes(pElasticVolStrainLabel_preReloc, matlset);
  task->computes(pStressQSLabel_preReloc,         matlset);

  // Add Yield Function computes and requires
  d_yield->addComputesAndRequires(task, matl, patches);

  // Add internal variable computes and requires
  d_intvar->addComputesAndRequires(task, matl, patches);

  // Add backstress computes and requires
  d_backstress->addComputesAndRequires(task, matl, patches);
}

//T2D: Throw exception that this is not supported
void Arenisca3PartiallySaturated::addComputesAndRequires(Task* ,
                                                         const MPMMaterial* ,
                                                         const PatchSet* ,
                                                         const bool, 
                                                         const bool ) const
{
  cout << "NO Implicit VERSION OF addComputesAndRequires EXISTS YET FOR Arenisca3PartiallySaturated"<<endl;
}


/*! ---------------------------------------------------------------------------------------
 *  This is needed for converting from one material type to another.  The functionality
 *  has been removed from the main Uintah branch.
 *  ---------------------------------------------------------------------------------------
 */
void 
Arenisca3PartiallySaturated::allocateCMDataAdd(DataWarehouse* new_dw,
                                               ParticleSubset* addset,
                                               ParticleLabelVariableMap* newState,
                                               ParticleSubset* delset,
                                               DataWarehouse* old_dw)
{
  std::ostringstream out;
  out << "Material conversion after failure not implemented for Arenisca.";
  throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  //task->requires(Task::NewDW, pPorosityLabel_preReloc,         matlset, Ghost::None);
  //task->requires(Task::NewDW, pSaturationLabel_preReloc,       matlset, Ghost::None);
}

/*---------------------------------------------------------------------------------------
 * MPMICE Hooks
 *---------------------------------------------------------------------------------------*/
double Arenisca3PartiallySaturated::computeRhoMicroCM(double pressure,
                                                      const double p_ref,
                                                      const MPMMaterial* matl,
                                                      double temperature,
                                                      double rho_guess)
{
  double rho_0 = matl->getInitialDensity();
  double K0 = d_cm.K0_Murnaghan_EOS;
  double n = d_cm.n_Murnaghan_EOS;

  double p_gauge = pressure - p_ref;
  double rho_cur = rho_0*std::pow(((n*p_gauge)/K0 + 1), (1.0/n));

  return rho_cur;
}

void Arenisca3PartiallySaturated::computePressEOSCM(double rho_cur,
                                                    double& pressure, double p_ref,
                                                    double& dp_drho, 
                                                    double& soundSpeedSq,
                                                    const MPMMaterial* matl,
                                                    double temperature)
{
  double rho_0 = matl->getInitialDensity();
  double K0 = d_cm.K0_Murnaghan_EOS;
  double n = d_cm.n_Murnaghan_EOS;

  double eta = rho_cur/rho_0;
  double p_gauge = K0/n*(std::pow(eta, n) - 1.0);

  double bulk = K0 + n*p_gauge;
  // double nu = 0.0;
  double shear = 1.5*bulk;

  pressure = p_ref + p_gauge;
  dp_drho  = K0*std::pow(eta, n-1);
  soundSpeedSq = (bulk + 4.0*shear/3.0)/rho_cur;  // speed of sound squared
}

double Arenisca3PartiallySaturated::getCompressibility()
{
  cout << "NO VERSION OF getCompressibility EXISTS YET FOR Arenisca3PartiallySaturated"
       << endl;
  return 1.0/d_cm.K0_Murnaghan_EOS;
}

