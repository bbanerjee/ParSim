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
    cm_ps->setAttribute("type","Arenisca3PartiallySaturated");
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
  // And the initial porosity and saturation to the parameter dictionary
  ParameterDict allParams;
  allParams["phi0"] = d_fluidParam.phi0;
  allParams["Sw0"] = d_fluidParam.Sw0;

  // Get the particles in the current patch
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(),patch);

  // Get the particle volume and mass
  constParticleVariable<double> pVolume, pMass;
  new_dw->get(pVolume, lb->pVolumeLabel, pset);
  new_dw->get(pMass,   lb->pMassLabel,   pset);

  // Initial variables for yield function parameter variability
  d_yield->initializeLocalVariables(patch, pset, new_dw, pVolume);

  // Initial variables for internal variables (needs yield function initialized first)
  ParameterDict yieldParams = d_yield->getParameters();
  allParams.insert(yieldParams.begin(), yieldParams.end());
  std::cout << "Model parameters are: " << std::endl;
  for (auto param : allParams) {
    std::cout << "\t \t" << param.first << " " << param.second << std::endl;
  }
  d_intvar->initializeInternalVariable(patch, matl, pset, new_dw, lb, allParams);

  // Initial variables for backstress model
  ParameterDict backstressParams = d_backstress->getParameters();
  allParams.insert(backstressParams.begin(), backstressParams.end());
  d_backstress->initializeLocalVariables(patch, pset, new_dw, pVolume);

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
    d_yield->copyLocalVariables(pset, old_dw, new_dw);

    // Get the yield condition parameter variables
    std::vector<constParticleVariable<double> > pYieldParams = 
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
    constParticleVariable<double>  pPhi_old;
    constParticleVariable<double>  pSw_old;
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

    old_dw->get(pPhi_old,          pPorosityLabel,               pset); 
    old_dw->get(pSw_old,           pSaturationLabel,             pset); 
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

    ParticleVariable<double>  pPhi_new, pSw_new;
    ParticleVariable<int>     pLocalized_new;
    ParticleVariable<double>  pElasticVolStrain_new;
    ParticleVariable<Matrix3> pStressQS_new;
    new_dw->allocateAndPut(pPhi_new,              pPorosityLabel_preReloc,         pset);
    new_dw->allocateAndPut(pSw_new,               pSaturationLabel_preReloc,       pset);
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
      // may ativate this feature.
      pdTdt[idx] = 0.0;

      // Copy particle deletion variable
      pLocalized_new[idx] = pLocalized[idx];

      // Compute the symmetric part of the velocity gradient
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

      // Compute the unrotated stress at the first of the current timestep
      Matrix3 sigma_old = (tensorR.Transpose())*(pStress_old[idx]*tensorR);
      Matrix3 sigmaQS_old = (tensorR.Transpose())*(pStressQS_old[idx]*tensorR);

      // initial assignment for the updated values of plastic strains, volumetric
      // part of the plastic strain, volumetric part of the elastic strain, kappa,
      // and the backstress. tentative assumption of elasticity
      ModelState_MasonSand state_old;
      state_old.capX                = pCapX[idx];
      state_old.kappa               = pKappa[idx];
      state_old.zeta                = pBackStress[idx].Trace();
      state_old.stressTensor        = &sigmaQS_old;
      std::cout << "Plastic strain tensor = " << pEp[idx] << std::endl;
      state_old.plasticStrainTensor = &(pEp[idx]);
      state_old.p3                  = pP3[idx];
      state_old.porosity            = pPhi_old[idx];
      state_old.saturation          = pSw_old[idx];


      // Get the parameters of the yield surface (for variability)
      for (auto& pYieldParamVar: pYieldParams) {
        state_old.yieldParams.push_back(pYieldParamVar[idx]);
      }

      // Compute the elastic moduli at t = t_n
      computeElasticProperties(state_old);

      // Rate-independent plastic step
      // Divides the strain increment into substeps, and calls substep function
      ModelState_MasonSand state_new;
      bool isSuccess = rateIndependentPlasticUpdate(D, delT, yieldParams,
                                                    idx, pParticleID[idx], state_old,
                                                    state_new);

      if (isSuccess) {

        pStressQS_new[idx] = *(state_new.stressTensor); // unrotated stress at end of step
        pCapX_new[idx] = state_new.capX;                 // hydrostatic compressive strength at end of step
        pKappa_new[idx] = state_new.kappa;               // branch point
        pBackStress_new[idx] = Identity*state_new.zeta;  // trace of isotropic backstress at end of step
        pEp_new[idx] = *(state_new.plasticStrainTensor); // plastic strain at end of step
        pEpv_new[idx] = pEp_new[idx].Trace();            // Plastic volumetric strain at end of step
        pP3_new[idx] = pP3[idx];

        // Elastic volumetric strain at end of step, compute from updated deformation gradient.
        pElasticVolStrain_new[idx] = log(pDefGrad_new[idx].Determinant()) - pEpv_new[idx];

      } else {

        // If the updateStressAndInternalVars function can't converge it will return false.  
        // This indicates substepping has failed, and the particle will be deleted.
        pLocalized_new[idx]=-999;
        cout << "** WARNING ** Bad step, deleting particle"
             << " idx = " << idx 
             << " particleID = " << pParticleID[idx] << std::endl;
        pStressQS_new[idx] = pStressQS_old[idx];
        pCapX_new[idx] = state_old.capX; 
        pKappa_new[idx] = state_old.kappa;
        pBackStress_new[idx] = Identity*state_old.zeta;  // trace of isotropic backstress at end of step
        pEp_new[idx] = *(state_old.plasticStrainTensor);          // plastic strain at end of step
        pEpv_new[idx] = pEp_new[idx].Trace();
        pP3_new[idx] = pP3[idx];
        pElasticVolStrain_new[idx] = pElasticVolStrain[idx];
      }

      // Rate-dependent plastic step
      ModelState_MasonSand stateQS_old(state_old);
      stateQS_old.stressTensor = &pStressQS_old[idx];
      ModelState_MasonSand stateQS_new(state_new);
      stateQS_new.stressTensor = &pStressQS_new[idx];
      
      computeElasticProperties(stateQS_old);
      computeElasticProperties(stateQS_new);
 
      rateDependentPlasticUpdate(D, delT, yieldParams, stateQS_old, stateQS_new, state_old,
                                 pStress_new[idx]);

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

      // Compute wave speed + particle velocity at each particle, store the maximum
      computeElasticProperties(stateQS_new); 
      double bulk = stateQS_new.bulkModulus;
      double shear = stateQS_new.shearModulus;
      double rho_cur = pMass[idx]/pVolume[idx];
      c_dil = sqrt((bulk+four_third*shear)/rho_cur);
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
                                                          const ModelState_MasonSand& state_n,
                                                          ModelState_MasonSand& state_new)
{
  // Compute the trial stress
  Matrix3 strain_inc = D*delT;
  Matrix3 stress_trial = computeTrialStress(state_n, strain_inc);

  // Set up a trial state, update the stress invariants, and compute elastic properties
  ModelState_MasonSand state_trial;
  state_trial.stressTensor = &stress_trial;
  state_trial.updateStressInvariants();
  computeElasticProperties(state_trial);
  
  // Determine the number of substeps (nsub) based on the magnitude of
  // the trial stress increment relative to the characteristic dimensions
  // of the yield surface.  Also compare the value of the pressure dependent
  // elastic properties at sigma_old and sigma_trial and adjust nsub if
  // there is a large change to ensure an accurate solution for nonlinear
  // elasticity even with fully elastic loading.
  int nsub = computeStepDivisions(idx, pParticleID, state_n, state_trial, yieldParams);

  // * Upon FAILURE *
  // Delete the particle if the number of substeps is unreasonable
  // Send ParticleDelete Flag to Host Code, Store Inputs to particle data:
  // input values for sigma_new, X_new, Zeta_new, ep_new, along with error flag
  if (nsub < 0) {
    state_new = state_n;
    std::cout << "Step Failed: Particle idx = " << idx << " ID = " << pParticleID << std::endl;
    // bool success  = false;
    return false;
  }

  // Compute a subdivided time step:
  // Loop at least once or until substepping is successful
  const int CHI_MAX = 5;       // max allowed subcycle multiplier

  double dt = delT/nsub;       // substep time increment
  ModelState_MasonSand state_old(state_n);

  int chi = 1;                 // subcycle multiplier
  double tlocal = 0.0;
  bool isSuccess = false;

  do {

    // Compute the elastic properties based on the stress and plastic strain at
    // the start of the substep.  These will be constant over the step unless elastic-plastic
    // is used to modify the tangent stiffness in the consistency bisection iteration.
    computeElasticProperties(state_old);

    //  Call substep function {sigma_new, ep_new, X_new, Zeta_new}
    //    = computeSubstep(D, dt, sigma_substep, ep_substep, X_substep, Zeta_substep)
    //  Repeat while substeps continue to be successful
    isSuccess = computeSubstep(D, dt, yieldParams, state_old, state_new);
    if (isSuccess) {

      tlocal += dt;
      state_old = state_new;

    } else {

      // Substepping has failed; increase chi
      chi *= 2; 
      dt /= 2.0;

      if (chi > CHI_MAX) {
        state_new = state_old;
        return isSuccess;
      }

    }
  } while (tlocal < delT);
    
  return isSuccess;

} 

/** 
 * Method: computeElasticProperties
 *
 * Purpose: 
 *   Compute the bulk and shear mdoulus at a given state
 */
void 
Arenisca3PartiallySaturated::computeElasticProperties(ModelState_MasonSand& state)
{
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
  Matrix3 stress_old = *(state_old.stressTensor);
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
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
  }

  int nmax = d_cm.subcycling_characteristic_number;
  
  // Compute change in bulk modulus:
  double bulk_old = state_old.bulkModulus;
  double bulk_trial = state_trial.bulkModulus;

  int n_bulk = std::ceil(std::abs(bulk_old - bulk_trial)/bulk_old);  
  
  // Compute trial stress increment relative to yield surface size:
  Matrix3 d_sigma = *(state_trial.stressTensor) - *(state_old.stressTensor);
  double size = 0.5*(PEAKI1 - state_old.capX);
  if (STREN > 0.0){
    size = std::min(size, STREN);
  }  
  int n_yield = ceil(1.0e-4*d_sigma.Norm()/size);

  // nsub is the maximum of the two values.above.  If this exceeds allowable,
  // throw warning and delete particle.
  int nsub = std::max(n_bulk, n_yield);
 
  if (nsub > d_cm.subcycling_characteristic_number) {
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
  } else {
    nsub = std::min(std::max(nsub,1), nmax);
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
                                            const ModelState_MasonSand& state_old,
                                            ModelState_MasonSand& state_new)
{
  // Compute the trial stress
  Matrix3 deltaEps = D*dt;
  Matrix3 stress_trial = computeTrialStress(state_old, deltaEps);

  // Set up a trial state, update the stress invariants
  ModelState_MasonSand state_trial(state_old);
  state_trial.stressTensor = &stress_trial;
  state_trial.updateStressInvariants();

  // Evaluate the yield function at the trial stress:
  int isElastic = (int) d_yield->evalYieldCondition(&state_trial); 

  // Elastic substep
  if (isElastic == 0 || isElastic == -1) { 
    state_new = state_trial;
    return true; // bool isSuccess = true;
  }

  // Elastic-plastic or fully-plastic substep
  // Compute non-hardening return to initial yield surface:
  // returnFlag would be != 0 if there was an error in the nonHardeningReturn call, but
  // there are currently no tests in that function that could detect such an error.
  Matrix3 sig_0(0.0);               // final stress state for non-hardening return
  Matrix3 deltaEps_p_0(0.0);        // increment in plastic strain for non-hardening return
  int returnFlag = nonHardeningReturn(deltaEps, state_old, state_trial, yieldParams,
                                      sig_0, deltaEps_p_0);
  if (!returnFlag) {
    state_new = state_old;
    return false; // bool isSuccess = false;
  }

  // Do "consistency bisection"
  state_new = state_old;
  bool isSuccess = consistencyBisection(deltaEps, state_old, state_trial,
                                        deltaEps_p_0, sig_0, yieldParams, 
                                        state_new);

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
int 
Arenisca3PartiallySaturated::nonHardeningReturn(const Uintah::Matrix3& strain_inc,
                                                const ModelState_MasonSand& state_old,
                                                const ModelState_MasonSand& state_trial,
                                                const ParameterDict& params,
                                                Uintah::Matrix3& sig_new,
                                                Uintah::Matrix3& plasticStrain_inc_new)
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
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
  }

  // Set up tolerance
  const double TOLERANCE = 1.0e-6;

  // Set up quantities that do not vary
  const double K_over_G = std::sqrt(1.5*state_old.bulkModulus/state_old.shearModulus);

  // Save the r and z Lode coordinates for the trial stress state
  double r_trial = BETA*state_trial.rr;
  double z_trial = state_trial.zz;

  // Compute transformed r coordinates
  double rprime_trial = r_trial*K_over_G;

  // Compute an initial interior point inside the yield surface
  double I1_0 = state_old.zeta + 0.5*(state_old.capX + PEAKI1);
  double sqrtJ2_0 = 0.0;
  double r_0 = BETA*sqrtJ2_0;
  double z_0 = I1_0/std::sqrt(3.0);
  double rprime_0 = r_0*std::sqrt(1.5*state_old.bulkModulus/state_old.shearModulus);

  // Initialize theta
  double theta = 0.0;

  // Initialize new and rotated transformed Lode coordinates
  double z_new = 0.0, z_rot = 0.0;
  double rprime_new = 0.0, rprime_rot = 0.0;
  
  // Loop begin
  do {

    // Apply the bisection algorithm to find new rprime and z
    applyBisectionAlgorithm(z_0, rprime_0, z_trial, rprime_trial, state_old, params,
                            z_new, rprime_new);

    // Apply rotation algorithm to find new internal point
    theta = findNewInternalPoint(z_trial, rprime_trial, z_new, rprime_new, 
                                 theta, state_old, params,
                                 z_rot, rprime_rot);

    // Update transformed Lode coordinates
    z_0 = z_rot;
    rprime_0 = rprime_rot;
    
  } while (std::abs(theta) > TOLERANCE);

  // Compute updated invariants
  double I1_new = std::sqrt(3.0)*z_new;
  double sqrtJ2_new = (1.0/(K_over_G*BETA*std::sqrt(2)))*rprime_new;

  // Compute new stress
  Matrix3 sig_dev = state_trial.deviatoricStressTensor;
  sig_new = one_third*I1_new*Identity + (sqrtJ2_new/state_trial.sqrt_J2)*sig_dev;

  // Compute new plastic strain increment
  Matrix3 sig_inc = sig_new - *(state_old.stressTensor);
  plasticStrain_inc_new = strain_inc - 
     one_third*(1.0/(3.0*state_old.bulkModulus) - 0.5/state_old.shearModulus)*sig_inc.Trace()*Identity - 
     (0.5/state_old.shearModulus)*sig_inc;

  return 0;
} //===================================================================

/**
 * Method: applyBisectionAlgorithm
 * Purpose: 
 *   Uses bisection to find the intersction of a loading path with the yield surface
 *   Returns location of intersection point in transformed stress space
 */
void
Arenisca3PartiallySaturated::applyBisectionAlgorithm(const double& z_0, const double& rprime_0, 
                                                     const double& z_trial, const double& rprime_trial, 
                                                     const ModelState_MasonSand& state_old, 
                                                     const ParameterDict& params,
                                                     double &z_new, double &rprime_new)
{
  const double TOLERANCE = 1.0e-6;

  double eta_in = 0.0, eta_out = 1.0;
  double eta_mid = 0.5, z_mid = 0.0, rprime_mid = 0.0;

  while (!(std::abs(eta_out - eta_in) < TOLERANCE)) {
    eta_mid = 0.5*(eta_in + eta_out);
    z_mid = eta_mid*(z_trial - z_0) + z_0;
    rprime_mid = eta_mid*(rprime_trial - rprime_0) + rprime_0;
    bool isElastic = evalYieldCondition(z_mid, rprime_mid, state_old, params);
    if (isElastic) {
      eta_in = eta_mid;
    } else {
      eta_out = eta_mid;
    }
  } // end while

  z_new = z_mid;
  rprime_new = rprime_mid;
}

/**
 * Method: findNewInternalPoint
 * Purpose: 
 *   Apply rotation algorithm to rotate the stress around the trial state to find
 *   a new internal point
 *   Returns location of the internal point and the angle
 */
double
Arenisca3PartiallySaturated::findNewInternalPoint(const double& z_trial, const double& rprime_trial, 
                                                  const double& z_new, const double& rprime_new, 
                                                  const double& theta_old,
                                                  const ModelState_MasonSand& state_old, 
                                                  const ParameterDict& params,
                                                  double& z_rot, double& rprime_rot)
{
  bool isElastic = true;
  int nn = 0;
  double sinTheta = 0.0, cosTheta = 0.0;
  do {
    
    // To avoid the cost of computing pow() to get theta, and then sin(), cos(),
    // we use the lookup table defined above by sinV and cosV.
    //
    // theta = pi_fourth*Pow(-two_third,n);
    // z_test = z_trial + cos(theta)*(z_0-z_trial) - sin(theta)*(r_0-r_trial);
    // r_test = r_trial + sin(theta)*(z_0-z_trial) + cos(theta)*(r_0-r_trial);
    sinTheta = Arenisca3PartiallySaturated::sinV[nn];
    cosTheta = Arenisca3PartiallySaturated::cosV[nn];
    z_rot = z_trial + cosTheta*(z_new - z_trial) - sinTheta*(rprime_new - rprime_trial);
    rprime_rot = rprime_trial + sinTheta*(z_new - z_trial) + cosTheta*(rprime_new - rprime_trial);

    // Check yield condition
    isElastic = evalYieldCondition(z_rot, rprime_rot, state_old, params);

    // Increment n
    nn += 1;

  } while (isElastic && nn < Arenisca3PartiallySaturated::NMAX);

  double theta_rot = std::asin(sinV[nn-1]);
  return theta_rot;
}

/**
 * Method: evalYieldCondition
 * Purpose: 
 *   Evaluate the yield condition in transformed Lode coordinates
 *   Returns whether the stress is elastic or not
 */
bool
Arenisca3PartiallySaturated::evalYieldCondition(const double& z_stress, const double& rprime_stress, 
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
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
  }

  double G_over_K = std::sqrt(state_old.shearModulus/(1.5*state_old.bulkModulus));
  double I1_stress = std::sqrt(3.0)*z_stress;
  double sqrtJ2_stress = G_over_K*(1.0/(std::sqrt(2.0)*BETA))*rprime_stress;

  // Create a temporary state for evaluation the yield function
  ModelState_MasonSand state_stress(state_old);
  state_stress.I1 = I1_stress;
  state_stress.sqrt_J2 = sqrtJ2_stress;

  // Evaluate the yield function
  bool isElastic = d_yield->evalYieldCondition(&state_stress);
  return isElastic;
}

/**
 * Method: consistencyBisection
 * Purpose: 
 *   Find the updated stress for hardening plasticity using the consistency bisection 
 *   algorithm
 *   Returns whether the procedure is sucessful orhas failed
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
  const int    IMAX      = 93;   // imax = ceil(-10.0*log(TOL)); // Update this if TOL changes
  const int    JMAX      = 93;   // jmax = ceil(-10.0*log(TOL)); // Update this if TOL changes

  // Local trial stress
  ModelState_MasonSand state_trial_upd(state_trial);
  state_trial_upd.updateStressInvariants();

  // Initialize
  Matrix3 sig_old     = *(state_old.stressTensor);
  Matrix3 eps_p_old   = *(state_old.plasticStrainTensor);
  double  eps_p_v_old = eps_p_old.Trace();
  double  zeta_old    = state_old.zeta;

  Matrix3 sig_new          = sig_0;
  Matrix3 deltaEps_p_new   = deltaEps_p_0;
  double  deltaEps_p_v_new = deltaEps_p_0.Trace();
  double  capX_new = 0.0;
  double  zeta_new = 0.0;

  // Start loop
  int ii = 1;
  double eta_in = 0.0, eta_out = 1.0, eta_mid = 0.5;
  double norm_deltaEps_p_new = deltaEps_p_new.Norm();
  double norm_deltaEps_p_0 = eta_mid*norm_deltaEps_p_new;
  do {
    int jj = 1;
    bool isElastic = true;
    while (isElastic) {
      eta_mid = 0.5*(eta_in + eta_out); 
      double deltaEps_p_v_mid = eta_mid*deltaEps_p_v_new;
      double eps_p_v_mid = eps_p_v_old + deltaEps_p_v_mid;

      // Update hydrostatic compressive strength
      state_trial_upd.ep_v = eps_p_v_mid;
      state_trial_upd.p3 = state_old.p3;
      capX_new = computeHydrostaticStrength(state_trial_upd);  

      // Update the isotropic backstress
      state_trial_upd.dep_v = deltaEps_p_v_mid;
      state_trial_upd.phi0  = d_fluidParam.phi0;
      state_trial_upd.Sw0   = d_fluidParam.Sw0;
      Matrix3 backStress_new;
      d_backstress->computeBackStress(&state_trial_upd, backStress_new);
      zeta_new = one_third*backStress_new.Trace();

      // Update the trial stress
      state_trial_upd.capX = capX_new;
      state_trial_upd.zeta = zeta_new;

      // Test the yield condition
      isElastic = d_yield->evalYieldCondition(&state_trial_upd);

      // If the state is elastic, there is too much plastic strain and
      // the following will be used;
      // otherwise control will break out from the isElastic loop
      eta_out = eta_mid;
      jj++;

      // Too many iterations
      if (jj > JMAX) {
        state_new = state_old;
        // bool isSuccess = false;
        return false;
      }
      
    } // end while(isElastic)

    // Update the state andcompute elastic properties
    ModelState_MasonSand state_mid;
    Matrix3 sig_mid = (sig_old + sig_new)*0.5;
    Matrix3 eps_p_mid = eps_p_old + deltaEps_p_0*(eta_mid*0.5);
    state_mid.stressTensor = &sig_mid;
    state_mid.plasticStrainTensor = &eps_p_mid;
    state_mid.porosity = 0.0;  // *TODO*
    state_mid.saturation = 0.0;  // *TODO*
    computeElasticProperties(state_mid);

    // Do non hardening return with updated properties
    state_trial_upd.bulkModulus = state_mid.bulkModulus;
    state_trial_upd.shearModulus = state_mid.shearModulus;
    state_trial_upd.capX = capX_new;
    state_trial_upd.zeta = std::min(zeta_new, 0.0);
    int status = nonHardeningReturn(deltaEps_new, state_old, state_trial_upd, params,
                                    sig_new, deltaEps_p_new);

    // Set up variables for various tests
    Matrix3 sig_trial = *(state_trial_upd.stressTensor);
    double trial_new = (sig_trial - sig_new).Trace();
    double trial_0 = (sig_trial - sig_0).Trace();
    norm_deltaEps_p_new = deltaEps_p_new.Norm();
    norm_deltaEps_p_0 = eta_mid*deltaEps_p_0.Norm();

    if ((std::signbit(trial_new) != std::signbit(trial_0)) ||
        (norm_deltaEps_p_new > norm_deltaEps_p_0)) {
      eta_out = eta_mid;
    } else {
      if (norm_deltaEps_p_new < norm_deltaEps_p_0) {
        eta_in = eta_mid;
      }
    }

    // Increment i and check
    ii++;
    if (ii > IMAX) {
      state_new = state_old;
      // bool isSuccess = false;
      return false;
    }
  } while (std::abs(norm_deltaEps_p_new < norm_deltaEps_p_0) < TOLERANCE);

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
  zeta_new = one_third*backStress_new.Trace();

  // Update the state
  state_new = state_trial_upd;  
  state_new.stressTensor = &sig_new;
  state_new.plasticStrainTensor = &eps_p_new;
  state_new.capX = capX_new;

  // min() eliminates tensile fluid pressure from explicit integration error
  state_new.zeta = std::min(zeta_new, 0.0);

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
  double T1, T2;
  try {
    T1 = yieldParams.at("T1");
    T2 = yieldParams.at("T2");
  } catch (std::out_of_range) {
    std::ostringstream err;
    err << "**ERROR** Could not find yield parameters T1 and T2" << std::endl;
    for (auto param : yieldParams) {
      err << param.first << " " << param.second << std::endl;
      throw InternalError(err.str(), __FILE__, __LINE__);
    }
  }

  // Check if rate-dependent plasticity has been turned on
  if (T1 == 0.0 || T2 == 0.0) {

   // No rate dependence, the dynamic stress equals the static stress.
    pStress_new = *(stateStatic_new.stressTensor);
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
  Matrix3 sigmaQS_old = *(stateStatic_old.stressTensor);
  Matrix3 sigmaQS_new = *(stateStatic_new.stressTensor);
  Matrix3 sigma_old = *(stateDynamic_old.stressTensor);
  pStress_new = sigmaQS_new
          + ((sigma_trial - sigma_old) - (sigmaQS_new - sigmaQS_old))*RH
          + (sigma_old - sigmaQS_old)*rh;

  // bool isRateDependent = true;
  return false;
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

