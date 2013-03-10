#include <CCA/Components/MPM/ConstitutiveModel/BasicDamageModel.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Malloc/Allocator.h>
#include <Core/Math/Gaussian.h>
#include <Core/Math/Weibull.h>
#include <Core/Math/MusilRNG.h>
#include <Core/Util/DebugStream.h>
#include <cmath>
#include <iostream>

using namespace Uintah;
using namespace Vaango;

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

//__________________________________
//  To turn on debug flags
//  csh/tcsh : setenv SCI_DEBUG "BasicDamage:+"
//  bash     : export SCI_DEBUG="BasicDamage:+" )
//  default is OFF

static DebugStream cout_damage("BasicDamage", false);

BasicDamageModel::BasicDamageModel(ProblemSpecP& ps, MPMFlags* mpm_flag)
{
  flag = mpm_flag;
  actuallyCreateDamageModel(ps);
}

BasicDamageModel::BasicDamageModel(const BasicDamageModel* bdm)
{
  flag = bdm->flag;
  actuallyCopyDamageModel(bdm);
}

BasicDamageModel::~BasicDamageModel()
{
  deleteDamageVarLabels();
}

BasicDamageModel* 
BasicDamageModel::clone()
{
  return scinew BasicDamageModel(*this);
}

//-----------------------------------------------------------------------------------
// The objected constructed by the ConstitutiveModel class only contains MPM flags.
// This method actually creates a BasicDamageModel and has to be invoked
// inside an actual constitutive model.
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::actuallyCreateDamageModel(ProblemSpecP& ps)
{
  // read in the failure model data and set the erosion algorithm flags
  getDamageModelData(ps);

  // Initialize local VarLabels
  initializeDamageVarLabels();
}

//-----------------------------------------------------------------------------------
// Get damage flags and damage distribution functions that can be used by 
// all constitutive models
//-----------------------------------------------------------------------------------
void
BasicDamageModel::getDamageModelData(ProblemSpecP& ps)
{
  // Get the brittle damage data
  if (flag->d_erosionAlgorithm == "BrittleDamage") {
    getBrittleDamageData(ps);
  } else {
    ps->require("failure_criteria", d_failure_criteria);

    // Get the failure stress/strain data
    getFailureStressOrStrainData(ps);
  }
    
  // Set the erosion algorithm
  setErosionAlgorithm();
}

void 
BasicDamageModel::getBrittleDamageData(ProblemSpecP& ps)
{
  d_brittle_damage.modulus  = 1.0e9; // Initial young's modulus
  d_brittle_damage.r0b   = 57.0; // Initial energy threshold
  d_brittle_damage.Gf    = 11.2; // Fracture energy
  d_brittle_damage.constant_D = 0.1; // Shape constant in softening function 
  d_brittle_damage.maxDamageInc=0.1; // Maximum damage in a time step
  d_brittle_damage.allowRecovery=false; // Allow recovery
  d_brittle_damage.recoveryCoeff=1.0; // Fraction of recovery if allowed
  d_brittle_damage.printDamage=false;  // Print damage
  ps->get("brittle_damage_modulus",             d_brittle_damage.modulus);
  ps->get("brittle_damage_initial_threshold",   d_brittle_damage.r0b);
  ps->get("brittle_damage_fracture_energy",     d_brittle_damage.Gf);
  ps->get("brittle_damage_constant_D",          d_brittle_damage.constant_D);
  ps->get("brittle_damage_max_damage_increment",d_brittle_damage.maxDamageInc);
  ps->get("brittle_damage_allowRecovery",       d_brittle_damage.allowRecovery);
  ps->get("brittle_damage_recoveryCoeff",       d_brittle_damage.recoveryCoeff);
  ps->get("brittle_damage_printDamage",         d_brittle_damage.printDamage);
  if (d_brittle_damage.recoveryCoeff <0.0 || d_brittle_damage.recoveryCoeff>1.0)
  {
    std::cerr << "brittle_damage_recoveryCoeff must be between 0.0 and 1.0" << endl;
  }     
}

void 
BasicDamageModel::getFailureStressOrStrainData(ProblemSpecP& ps)
{
  d_epsf.mean   = 10.0; // Mean failure stress or strain
  d_epsf.std    = 0.0;  // Std. Dev or Weibull mod. for failure stres or strain
  d_epsf.seed   = 0; // seed for weibull distribution generator
  d_epsf.dist   = "constant";
  d_epsf.scaling = "none";
  // "exponent" is the value of n used in c=(Vbar/V)^(1/n)
  // By setting the default value to DBL_MAX, that makes 1/n=0, which makes c=1
  d_epsf.exponent= DBL_MAX; //Exponent used in vol. scaling of failure criteria
  d_epsf.refVol = 1.0; // Reference volume for scaling failure criteria
  d_epsf.t_char = 1.0e-99; // Characteristic time of damage evolution

  ps->require("failure_criteria", d_failure_criteria);

  if(d_failure_criteria!="MaximumPrincipalStress" &&
     d_failure_criteria!="MaximumPrincipalStrain" &&
     d_failure_criteria!="MohrColoumb"){
     // The above are the only acceptable options.  If not one of them, bail. 
     throw ProblemSetupException("<failure_criteria> must be either MaximumPrincipalStress, MaximumPrincipalStrain or MohrColoumb", __FILE__, __LINE__);

  }

  if(d_failure_criteria=="MohrColoumb"){
    // The cohesion value that MC needs is the "mean" value in the
    // FailureStressOrStrainData struct
    ps->require("friction_angle", d_friction_angle);
    ps->require("tensile_cutoff_fraction_of_cohesion", d_tensile_cutoff);
  }
    
  ps->require("failure_mean",d_epsf.mean); //Mean val. of failure stress/strain
  ps->get("failure_distrib", d_epsf.dist); //"constant", "weibull" or "gauss"

  // Only require std if using a non-constant distribution
  if(d_epsf.dist!="constant"){
    ps->require("failure_std", d_epsf.std); //Std dev (Gauss) or Weibull modulus
  }

  ps->get("scaling", d_epsf.scaling); //"none" or "kayenta"
  if(d_epsf.scaling!="none"){
    // If doing some sort of scaling, require user to provide a reference volume
    ps->require("reference_volume",d_epsf.refVol);
    if(d_epsf.dist=="weibull"){
      d_epsf.exponent=d_epsf.std;// By default, exponent is Weibull modulus, BUT
      ps->get("exponent", d_epsf.exponent); // allow user to choose the exponent
   } else {
      // Force user to choose the exponent
      ps->require("exponent", d_epsf.exponent);
    }
  }
  ps->get("failure_seed",    d_epsf.seed); //Seed for RN generator
  ps->get("char_time",       d_epsf.t_char); //Characteristic time for damage
}

void 
BasicDamageModel::setErosionAlgorithm()
{
  d_setStressToZero = false;
  d_allowNoTension  = false;
  d_allowNoShear    = false;
  d_brittleDamage   = false;
  if (flag->d_doErosion) {
    if (flag->d_erosionAlgorithm == "AllowNoTension") 
      d_allowNoTension  = true;
    else if (flag->d_erosionAlgorithm == "ZeroStress") 
      d_setStressToZero = true;
    else if (flag->d_erosionAlgorithm == "AllowNoShear") 
      d_allowNoShear    = true;
    else if (flag->d_erosionAlgorithm == "BrittleDamage")
      d_brittleDamage   = true;
  }
}

void
BasicDamageModel::initializeDamageVarLabels()
{
  pFailureStressOrStrainLabel = VarLabel::create("p.epsf",
                               ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel             = VarLabel::create("p.localized",
                               ParticleVariable<int>::getTypeDescription());
  pDamageLabel                = VarLabel::create("p.damage",
                               ParticleVariable<double>::getTypeDescription());
  pTimeOfLocLabel             = VarLabel::create("p.timeofloc",
                               ParticleVariable<double>::getTypeDescription());
  pFailureStressOrStrainLabel_preReloc = VarLabel::create("p.epsf+",
                               ParticleVariable<double>::getTypeDescription());
  pLocalizedLabel_preReloc    = VarLabel::create("p.localized+",
                               ParticleVariable<int>::getTypeDescription());
  pDamageLabel_preReloc       = VarLabel::create("p.damage+",
                               ParticleVariable<double>::getTypeDescription());
  pTimeOfLocLabel_preReloc    = VarLabel::create("p.timeofloc+",
                               ParticleVariable<double>::getTypeDescription());
}

//-----------------------------------------------------------------------------------
// The objected constructed by the ConstitutiveModel class does not have a copy constructor.
// This method actually copies a BasicDamageModel into another.
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::actuallyCopyDamageModel(const BasicDamageModel* bdm)
{
  // copy in the failure model data and set the erosion algorithm flags
  setDamageModelData(bdm);

  // Initialize local VarLabels
  initializeDamageVarLabels();
}

//-----------------------------------------------------------------------------------
// Get damage flags and damage distribution functions that can be used by 
// all constitutive models
//-----------------------------------------------------------------------------------
void
BasicDamageModel::setDamageModelData(const BasicDamageModel* bdm)
{
  if (flag->d_erosionAlgorithm == "BrittleDamage") {
    setBrittleDamageData(bdm);
  } else {
    // Set the failure strain data
    setFailureStressOrStrainData(bdm);
    d_failure_criteria = bdm->d_failure_criteria;
    if(d_failure_criteria=="MohrColoumb"){
      d_tensile_cutoff = bdm->d_tensile_cutoff;
      d_friction_angle = bdm->d_friction_angle;
    }
  }

  // Set the erosion algorithm
  setErosionAlgorithm(bdm);
}

void 
BasicDamageModel::setBrittleDamageData(const BasicDamageModel* cm)
{
  d_brittle_damage.modulus = cm->d_brittle_damage.modulus; // Initial modulus
  d_brittle_damage.r0b   = cm->d_brittle_damage.r0b; // Initial energy threshold
  d_brittle_damage.Gf    = cm->d_brittle_damage.Gf; // Fracture energy
  // Shape constant in softening function
  d_brittle_damage.constant_D=cm->d_brittle_damage.constant_D; 
  //maximum damage in a time step 
  d_brittle_damage.maxDamageInc=cm->d_brittle_damage.maxDamageInc; 
  //allow recovery
  d_brittle_damage.allowRecovery=cm->d_brittle_damage.allowRecovery;
  //fraction of recovery if allowed
  d_brittle_damage.recoveryCoeff=cm->d_brittle_damage.recoveryCoeff;
  //print damage
  d_brittle_damage.printDamage = cm->d_brittle_damage.printDamage;
}

void 
BasicDamageModel::setFailureStressOrStrainData(const BasicDamageModel* cm)
{
  d_epsf.mean            = cm->d_epsf.mean;
  d_epsf.std             = cm->d_epsf.std;
  d_epsf.seed            = cm->d_epsf.seed;
  d_epsf.dist            = cm->d_epsf.dist;
  d_epsf.scaling         = cm->d_epsf.scaling;
  d_epsf.exponent        = cm->d_epsf.exponent;
  d_epsf.refVol          = cm->d_epsf.refVol;
  d_epsf.t_char          = cm->d_epsf.t_char;
}

void 
BasicDamageModel::setErosionAlgorithm(const BasicDamageModel* cm)
{
  d_setStressToZero = cm->d_setStressToZero;
  d_allowNoTension  = cm->d_allowNoTension;
  d_allowNoShear    = cm->d_allowNoShear;
  d_brittleDamage   = cm->d_brittleDamage;
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void
BasicDamageModel::deleteDamageVarLabels()
{
  VarLabel::destroy(pFailureStressOrStrainLabel);
  VarLabel::destroy(pFailureStressOrStrainLabel_preReloc);
  VarLabel::destroy(pLocalizedLabel);
  VarLabel::destroy(pLocalizedLabel_preReloc);
  VarLabel::destroy(pDamageLabel);
  VarLabel::destroy(pDamageLabel_preReloc);
  VarLabel::destroy(pTimeOfLocLabel);
  VarLabel::destroy(pTimeOfLocLabel_preReloc);
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void
BasicDamageModel::outputProblemSpecDamage(ProblemSpecP& cm_ps)
{
  if (flag->d_erosionAlgorithm == "BrittleDamage") {
    cm_ps->appendElement("brittle_damage_modulus",
                          d_brittle_damage.modulus);
    cm_ps->appendElement("brittle_damage_initial_threshold",
                          d_brittle_damage.r0b);
    cm_ps->appendElement("brittle_damage_fracture_energy",
                          d_brittle_damage.Gf);
    cm_ps->appendElement("brittle_damage_constant_D",           
                          d_brittle_damage.constant_D);
    cm_ps->appendElement("brittle_damage_max_damage_increment", 
                          d_brittle_damage.maxDamageInc);
    cm_ps->appendElement("brittle_damage_allowRecovery",        
                          d_brittle_damage.allowRecovery);
    cm_ps->appendElement("brittle_damage_recoveryCoeff",        
                          d_brittle_damage.recoveryCoeff);
    cm_ps->appendElement("brittle_damage_printDamage",          
                          d_brittle_damage.printDamage);
  } else {
    cm_ps->appendElement("failure_mean",     d_epsf.mean);
    cm_ps->appendElement("failure_std",      d_epsf.std);
    cm_ps->appendElement("failure_exponent", d_epsf.exponent);
    cm_ps->appendElement("failure_seed" ,    d_epsf.seed);
    cm_ps->appendElement("failure_distrib",  d_epsf.dist);
    cm_ps->appendElement("failure_criteria", d_failure_criteria);
    cm_ps->appendElement("scaling",          d_epsf.scaling);
    cm_ps->appendElement("exponent",         d_epsf.exponent);
    cm_ps->appendElement("reference_volume", d_epsf.refVol);
    cm_ps->appendElement("char_time",        d_epsf.t_char);

    if(d_failure_criteria=="MohrColoumb"){
      cm_ps->appendElement("friction_angle", d_friction_angle);
      cm_ps->appendElement("tensile_cutoff_fraction_of_cohesion",
                                             d_tensile_cutoff);
    }
  } //end if BrittleDamage
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::addInitialComputesAndRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patches,
                                                MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pFailureStressOrStrainLabel, matlset);
  task->computes(pLocalizedLabel,             matlset);
  task->computes(pTimeOfLocLabel,             matlset);
  task->computes(pDamageLabel,                matlset);
  task->computes(lb->TotalLocalizedParticleLabel);
}  

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::addComputesAndRequires(Task* task,
                                         const MPMMaterial* matl,
                                         const PatchSet* patches,
                                         MPMLabel* lb) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  
  task->requires(Task::OldDW, pFailureStressOrStrainLabel,    matlset, Ghost::None);
  task->requires(Task::OldDW, pLocalizedLabel,                matlset, Ghost::None);
  task->requires(Task::OldDW, pTimeOfLocLabel,                matlset, Ghost::None);
  task->requires(Task::OldDW, pDamageLabel,                   matlset, Ghost::None);
  //task->requires(Task::NewDW, lb->pStressLabel_preReloc,      matlset, Ghost::None);
  task->modifies(lb->pStressLabel_preReloc);
    
  task->computes(pFailureStressOrStrainLabel_preReloc,        matlset);
  task->computes(pLocalizedLabel_preReloc,                    matlset);
  task->computes(pTimeOfLocLabel_preReloc,                    matlset);
  task->computes(pDamageLabel_preReloc,                       matlset);
  task->computes(lb->TotalLocalizedParticleLabel);   
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::allocateDamageDataAddRequires(Task* task,
                                                const MPMMaterial* matl,
                                                const PatchSet* patches,
                                                MPMLabel* ) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pFailureStressOrStrainLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pLocalizedLabel_preReloc,             matlset, Ghost::None);
  task->requires(Task::NewDW, pTimeOfLocLabel_preReloc,             matlset, Ghost::None);
  task->requires(Task::NewDW, pDamageLabel_preReloc,                matlset, Ghost::None);
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void
BasicDamageModel::initializeDamageData(const Patch* patch,
                                       const MPMMaterial* matl,
                                       DataWarehouse* new_dw,
                                       MPMLabel* lb)
{

  ParticleVariable<double>      pFailureStrain;
  ParticleVariable<int>         pLocalized;
  ParticleVariable<double>      pTimeOfLoc;
  constParticleVariable<double> pVolume;
  ParticleVariable<double>      pDamage;
    
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  new_dw->get(pVolume,                   lb->pVolumeLabel,            pset);
  new_dw->allocateAndPut(pFailureStrain, pFailureStressOrStrainLabel, pset);
  new_dw->allocateAndPut(pLocalized,     pLocalizedLabel,             pset);
  new_dw->allocateAndPut(pTimeOfLoc,     pTimeOfLocLabel,             pset);
  new_dw->allocateAndPut(pDamage,        pDamageLabel,                pset);
    
  ParticleSubset::iterator iter = pset->begin();

  if (d_brittleDamage) {
    for(;iter != pset->end();iter++){
      pFailureStrain[*iter] = d_brittle_damage.r0b;
      pLocalized[*iter]     = 0;
      pTimeOfLoc[*iter]     = 1.e99;;
      pDamage[*iter]        = 0.0;
    }
  }  else if (d_epsf.dist == "gauss"){
    // Initialize a gaussian random number generator

    // Make the seed differ for each patch, otherwise each patch gets the
    // same set of random #s.
    int patchID = patch->getID();
    int patch_div_32 = patchID/32;
    patchID = patchID%32;
    unsigned int unique_seed = ((d_epsf.seed+patch_div_32+1) << patchID);

    SCIRun::Gaussian gaussGen(d_epsf.mean,d_epsf.std,unique_seed,
                                d_epsf.refVol,d_epsf.exponent);
      
    for(;iter != pset->end();iter++){
      pFailureStrain[*iter] =  fabs(gaussGen.rand(pVolume[*iter]));
      pLocalized[*iter]     = 0;
      pTimeOfLoc[*iter]     = -1.e99;;
      pDamage[*iter]        = 0.0;
    }
  } else if (d_epsf.dist == "weibull"){
    // Initialize a weibull random number generator

    // Make the seed differ for each patch, otherwise each patch gets the
    // same set of random #s.
    int patchID = patch->getID();
    int patch_div_32 = patchID/32;
    patchID = patchID%32;
    unsigned int unique_seed = ((d_epsf.seed+patch_div_32+1) << patchID);

    SCIRun::Weibull weibGen(d_epsf.mean,d_epsf.std,d_epsf.refVol,
                              unique_seed,d_epsf.exponent);
      
    for(;iter != pset->end();iter++){
      pFailureStrain[*iter] = weibGen.rand(pVolume[*iter]);
      pLocalized[*iter]     = 0;
      pTimeOfLoc[*iter]     = -1.e99;;
      pDamage[*iter]        = 0.0;
    }
  } else if (d_epsf.dist == "uniform") {

    // Make the seed differ for each patch, otherwise each patch gets the
    // same set of random #s.
    int patchID = patch->getID();
    int patch_div_32 = patchID/32;
    patchID = patchID%32;
    unsigned int unique_seed = ((d_epsf.seed+patch_div_32+1) << patchID);
    MusilRNG* randGen = scinew MusilRNG(unique_seed);
    for(;iter != pset->end();iter++){
      pLocalized[*iter]     = 0;
      pTimeOfLoc[*iter]     = -1.e99;;

      double rand = (*randGen)(); 
      double range = (2*rand - 1)*d_epsf.std;
      double cc = pow(d_epsf.refVol/pVolume[*iter], 1.0/d_epsf.exponent); 
      double fail_eps = cc*(d_epsf.mean + range);
      pFailureStrain[*iter] = fail_eps;
      pDamage[*iter]        = 0.0;
    }
    delete randGen;

  } else {
    for(;iter != pset->end();iter++){
      pFailureStrain[*iter] = d_epsf.mean;
      pLocalized[*iter]     = 0;
      pTimeOfLoc[*iter]     = -1.e99;;
      pDamage[*iter]        = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::copyDamageDataFromDeletedToAddedParticle(DataWarehouse* new_dw,
                                                            ParticleSubset* addset,
                                                            map<const VarLabel*,
                                                              ParticleVariableBase*>* newState,
                                                            ParticleSubset* delset,
                                                            DataWarehouse* old_dw )
{
  constParticleVariable<double>  o_pFailureStrain;
  constParticleVariable<int>     o_pLocalized;
  constParticleVariable<double>  o_pTimeOfLoc;
  constParticleVariable<double>  o_pDamage;
  new_dw->get(o_pFailureStrain,  pFailureStressOrStrainLabel_preReloc,delset);
  new_dw->get(o_pLocalized,      pLocalizedLabel_preReloc,     delset);
  new_dw->get(o_pDamage,         pDamageLabel_preReloc,        delset);
  new_dw->get(o_pTimeOfLoc,      pTimeOfLocLabel_preReloc,     delset);

  ParticleVariable<double>       pFailureStrain;
  ParticleVariable<int>          pLocalized;
  ParticleVariable<double>       pTimeOfLoc;
  ParticleVariable<double>       pDamage;

  new_dw->allocateTemporary(pFailureStrain, addset);
  new_dw->allocateTemporary(pLocalized,     addset);
  new_dw->allocateTemporary(pTimeOfLoc,     addset);
  new_dw->allocateTemporary(pDamage,        addset);    

  ParticleSubset::iterator o,n     = addset->begin();
  for (o=delset->begin(); o != delset->end(); o++, n++) {
    pFailureStrain[*n]             = o_pFailureStrain[*o];
    pLocalized[*n]                 = o_pLocalized[*o];
    pTimeOfLoc[*n]                 = o_pTimeOfLoc[*o];
    pDamage[*n]                    = o_pDamage[*o];
  }
  (*newState)[pFailureStressOrStrainLabel] = pFailureStrain.clone();
  (*newState)[pLocalizedLabel]     = pLocalized.clone();
  (*newState)[pTimeOfLocLabel]     = pTimeOfLoc.clone();
  (*newState)[pDamageLabel]        = pDamage.clone();
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::carryForwardDamageData(ParticleSubset* pset,
                                         DataWarehouse*  old_dw,
                                         DataWarehouse*  new_dw,
                                         const MPMMaterial* matl)
{
  constParticleVariable<double>  pFailureStrain;
  constParticleVariable<int>     pLocalized;
  constParticleVariable<double>  pTimeOfLoc;
  constParticleVariable<double>  pDamage;
  ParticleVariable<double>       pFailureStrain_new;
  ParticleVariable<int>          pLocalized_new;
  ParticleVariable<double>       pTimeOfLoc_new;
  ParticleVariable<double>       pDamage_new;
      
  old_dw->get(pFailureStrain, pFailureStressOrStrainLabel,     pset);
  old_dw->get(pLocalized,     pLocalizedLabel,                 pset);
  old_dw->get(pTimeOfLoc,     pTimeOfLocLabel,                 pset);
  old_dw->get(pDamage,        pDamageLabel,                    pset);
      
  new_dw->allocateAndPut(pFailureStrain_new,    
                             pFailureStressOrStrainLabel_preReloc, pset);
  new_dw->allocateAndPut(pLocalized_new,      
                             pLocalizedLabel_preReloc,             pset);
  new_dw->allocateAndPut(pTimeOfLoc_new,      
                             pTimeOfLocLabel_preReloc,             pset);
  new_dw->allocateAndPut(pDamage_new,      
                             pDamageLabel_preReloc,                pset);      

  pFailureStrain_new.copyData(pFailureStrain);
  pLocalized_new.copyData(pLocalized);
  pTimeOfLoc_new.copyData(pTimeOfLoc);
  pDamage_new.copyData(pDamage);
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void
BasicDamageModel::addParticleState(std::vector<const VarLabel*>& from,
                                   std::vector<const VarLabel*>& to)
{
  from.push_back(pFailureStressOrStrainLabel);
  from.push_back(pLocalizedLabel);
  from.push_back(pDamageLabel);
  from.push_back(pTimeOfLocLabel);

  to.push_back(pFailureStressOrStrainLabel_preReloc);
  to.push_back(pLocalizedLabel_preReloc);
  to.push_back(pDamageLabel_preReloc);
  to.push_back(pTimeOfLocLabel_preReloc);
}

//-----------------------------------------------------------------------------------
// Add documentation here
//-----------------------------------------------------------------------------------
void 
BasicDamageModel::computeBasicDamage(const PatchSubset* patches,
                                     const MPMMaterial* matl,
                                     DataWarehouse* old_dw,
                                     DataWarehouse* new_dw,
                                     MPMLabel* lb)
{ 
  Ghost::GhostType  gac   = Ghost::AroundCells;
 
  // Normal patch loop
  for(int pp=0;pp<patches->size();pp++){
    const Patch* patch = patches->get(pp);

    long64 totalLocalizedParticle = 0;

    // Get particle info and patch info
    int dwi = matl->getDWIndex();
    ParticleSubset* pset = old_dw->getParticleSubset(dwi, patch);
    double time = d_sharedState->getElapsedTime();

    // Particle and grid data universal to model type
    // Old data containers
    constParticleVariable<int>     pLocalized;
    constParticleVariable<double>  pTimeOfLoc;
    constParticleVariable<double>  pFailureStrain, pDamage;
    constParticleVariable<long64>  pParticleID;

    constParticleVariable<double>  pVolume_new;
    constParticleVariable<Matrix3> pDefGrad_new;

    // New data containers
    ParticleVariable<int>          pLocalized_new;
    ParticleVariable<double>       pTimeOfLoc_new;
    ParticleVariable<double>       pFailureStrain_new, pDamage_new;
    ParticleVariable<Matrix3>      pStress;
    
    
    // Damage gets
    old_dw->get(pParticleID,              lb->pParticleIDLabel,        pset);
    old_dw->get(pLocalized,               pLocalizedLabel,             pset);
    old_dw->get(pTimeOfLoc,               pTimeOfLocLabel,             pset);
    old_dw->get(pFailureStrain,           pFailureStressOrStrainLabel, pset);
    old_dw->get(pDamage,                  pDamageLabel,                pset);

    new_dw->get(pDefGrad_new,             lb->pDefGradLabel_preReloc,  pset);
    new_dw->get(pVolume_new,              lb->pVolumeLabel_preReloc,   pset);

    new_dw->getModifiable(pStress,        lb->pStressLabel_preReloc,   pset);
      
    new_dw->allocateAndPut(pLocalized_new,     pLocalizedLabel_preReloc,              pset);
    new_dw->allocateAndPut(pTimeOfLoc_new,     pTimeOfLocLabel_preReloc,              pset);
    new_dw->allocateAndPut(pFailureStrain_new, pFailureStressOrStrainLabel_preReloc,  pset);
    new_dw->allocateAndPut(pDamage_new,        pDamageLabel_preReloc,                 pset);
      
    // Copy failure strains to new dw
    pFailureStrain_new.copyData(pFailureStrain);
    
    for(ParticleSubset::iterator iter = pset->begin(); iter != pset->end(); iter++){
      particleIndex idx = *iter;

      pDamage_new[idx] = pDamage[idx];

      // Modify the stress if particle has failed/damaged
      if (d_brittleDamage) {
        cout_damage << "Before update: Particle = " << idx << " pDamage = " << pDamage[idx] << " pStress = " << pStress[idx] << endl;
        updateDamageAndModifyStress(pDefGrad_new[idx], pFailureStrain[idx],
                                    pFailureStrain_new[idx], pVolume_new[idx],
                                    pDamage[idx], pDamage_new[idx],
                                    pStress[idx], pParticleID[idx]);
        pLocalized_new[idx] = (pDamage[idx] > 1.0) ? 0 : 1;
        cout_damage << "After update: Particle = " << idx << " pDamage = " << pDamage_new[idx] << " pStress = " << pStress[idx] << endl;
        //pLocalized_new[idx]= pLocalized[idx]; //not really used.
        if (pDamage_new[idx]>0.0) totalLocalizedParticle+=1;
      }
      else {
        updateFailedParticlesAndModifyStress(pDefGrad_new[idx], pFailureStrain[idx], 
                                           pLocalized[idx], pLocalized_new[idx],
                                           pTimeOfLoc[idx], pTimeOfLoc_new[idx],
                                           pStress[idx], pParticleID[idx],time);
        if (pLocalized_new[idx]>0){
          totalLocalizedParticle+=1;
        }
      }
    } // end loop over particles
    
    new_dw->put(sumlong_vartype(totalLocalizedParticle), lb->TotalLocalizedParticleLabel);
 
  }
}

// Modify the stress for brittle damage
// Update pFailureStrain_new (energy threshold)
// pDamage_new (damage; if negative, damage inactive), pStress
void 
BasicDamageModel::updateDamageAndModifyStress(const Matrix3& defGrad, 
                                     const double& pFailureStrain, 
				     double& pFailureStrain_new,
                                     const double& pVolume, 
                                     const double& pDamage,
                                     double& pDamage_new, 
                                     Matrix3& pStress,
				     const long64 particleID)
 {
  Matrix3 Identity, zero(0.0); Identity.Identity();
  double tau_b;  // current 'energy'

  // mean stress
  double pressure = (1.0/3.0)*pStress.Trace();


  // Check for damage (note that pFailureStrain is the energy threshold)
  pFailureStrain_new = pFailureStrain;

  if (pressure <0.0) { 

    //no damage if compressive
    if (pDamage <=0.0) { // previously no damage, do nothing
      return;
    } else { 
      //previously damaged, deactivate damage?
      if (d_brittle_damage.allowRecovery) {  //recovery
        pStress = pStress*d_brittle_damage.recoveryCoeff;
	pDamage_new = -pDamage; //flag damage to be negative
      
       if (d_brittle_damage.printDamage) std::cout << "Particle " << particleID << " damage halted: damage=" << pDamage_new << endl;
      }
      else
	pStress = pStress*(1.0-pDamage); // no recovery (default)
    }
  } //end pDamage <=0.0

  // pressure >0.0; possible damage
  else {

      // Compute Finger tensor (left Cauchy-Green) 
      Matrix3 bb = defGrad*defGrad.Transpose();
      // Compute Eulerian strain tensor
      Matrix3 ee = (Identity - bb.Inverse())*0.5;      
      // Compute the maximum principal strain
      double epsMax=0.,epsMed=0.,epsMin=0.;
      ee.getEigenValues(epsMax,epsMed,epsMin);

      // Young's modulus
      double young = d_brittle_damage.modulus;

      tau_b = sqrt(young*epsMax*epsMax);

      if (tau_b > pFailureStrain) {  
      // further damage
        // equivalent dimension of the particle
        double particleSize = pow(pVolume, 1.0/3.0);
        double r0b = d_brittle_damage.r0b;
	double const_D=d_brittle_damage.constant_D;
	double const_C = r0b*particleSize*(1.0+const_D) \
               /(d_brittle_damage.Gf*const_D)*log(1.0+const_D);
	double d1=1.0+const_D*exp(-const_C*(tau_b-r0b));
	double damage=0.999/const_D*((1.0+const_D)/d1 - 1.0);

	// Restrict the maximum damage in a time step for stability reason.
	if ((damage-pDamage) > d_brittle_damage.maxDamageInc) {
	  damage=pDamage+d_brittle_damage.maxDamageInc;
	}
	// Update threshold and damage
	pFailureStrain_new = tau_b;
	pDamage_new = damage;

	// Update stress
	pStress = pStress*(1.0-damage);
        if (d_brittle_damage.printDamage){
          std::cout << "Particle " << particleID << " damaged: "
               << " damage=" << pDamage_new << " epsMax=" << epsMax 
               << " tau_b=" << tau_b << endl;
        }
      } else {
	if (pDamage==0.0) return; // never damaged

	//current energy less than previous; deactivate damage?
	if (d_brittle_damage.allowRecovery) { //recovery
          pStress = pStress*d_brittle_damage.recoveryCoeff;
	  pDamage_new = -pDamage; //flag it to be negative
          if (d_brittle_damage.printDamage){
            std::cout << "Particle " << particleID << " damage halted: damage=" 
                 << pDamage_new << endl;
          }
	}
	else { //no recovery (default)
	  pStress = pStress*(1.0-pDamage);
          if (d_brittle_damage.printDamage){
            std::cout << "Particle " << particleID << " damaged: " 
                 << " damage=" << pDamage_new << " epsMax=" << epsMax 
                 << " tau_b=" << tau_b << endl;
          }
	}
      } // end if tau_b > pFailureStrain

  } //end if pressure

}

// Modify the stress if particle has failed
void 
BasicDamageModel::updateFailedParticlesAndModifyStress(const Matrix3& defGrad,
                                                 const double& pFailureStr,
                                                 const int& pLocalized,
                                                 int& pLocalized_new,
                                                 const double& pTimeOfLoc,
                                                 double& pTimeOfLoc_new,
                                                 Matrix3& pStress,
                                                 const long64 particleID,
                                                 double time)
{
  Matrix3 Identity, zero(0.0); Identity.Identity();

  // Find if the particle has failed
  pLocalized_new = pLocalized;
  pTimeOfLoc_new = pTimeOfLoc;
  if (pLocalized == 0){
    if(d_failure_criteria=="MaximumPrincipalStress"){
      double maxEigen=0.,medEigen=0.,minEigen=0.;
      pStress.getEigenValues(maxEigen,medEigen,minEigen);
      //The first eigenvalue returned by "eigen" is always the largest 
      if (maxEigen > pFailureStr){
        pLocalized_new = 1;
      }
      if (pLocalized != pLocalized_new) {
        std::cout << "Particle " << particleID << " has failed : MaxPrinStress = "
             << maxEigen << " eps_f = " << pFailureStr << endl;
        pTimeOfLoc_new = time;
      }
    }
    else if(d_failure_criteria=="MaximumPrincipalStrain"){
      // Compute Finger tensor (left Cauchy-Green) 
      Matrix3 bb = defGrad*defGrad.Transpose();
      // Compute Eulerian strain tensor
      Matrix3 ee = (Identity - bb.Inverse())*0.5;

      double maxEigen=0.,medEigen=0.,minEigen=0.;
      ee.getEigenValues(maxEigen,medEigen,minEigen);
      if (maxEigen > pFailureStr){
        pLocalized_new = 1;
      }
      if (pLocalized != pLocalized_new) {
        std::cout << "Particle " << particleID << " has failed : eps = " << maxEigen
             << " eps_f = " << pFailureStr << endl;
        pTimeOfLoc_new = time;
      }
    }
    else if(d_failure_criteria=="MohrColoumb"){
      double maxEigen=0.,medEigen=0.,minEigen=0.;
      pStress.getEigenValues(maxEigen,medEigen,minEigen);
  
      double cohesion = pFailureStr;
  
      double epsMax=0.;
      // Tensile failure criteria (max princ stress > d_tensile_cutoff*cohesion)
      if (maxEigen > d_tensile_cutoff*cohesion){
        pLocalized_new = 1;
        epsMax = maxEigen;
      }

      //  Shear failure criteria (max shear > cohesion + friction)
      double friction_angle = d_friction_angle*(M_PI/180.);

      if ( (maxEigen - minEigen)/2.0 > cohesion*cos(friction_angle)
           - (maxEigen + minEigen)*sin(friction_angle)/2.0){
        pLocalized_new = 2;
        epsMax = (maxEigen - minEigen)/2.0;
      }
      if (pLocalized != pLocalized_new) {
        std::cout << "Particle " << particleID << " has failed : maxPrinStress = "
             << epsMax << " cohesion = " << cohesion << endl;
        pTimeOfLoc_new = time;
      }
    } // Mohr-Coloumb
  } // pLocalized==0

  // If the particle has failed, apply various erosion algorithms
  if (flag->d_doErosion) {
    // Compute pressure
    double pressure = pStress.Trace()/3.0;
    double failTime = time - pTimeOfLoc_new;
    double D = exp(-failTime/d_epsf.t_char);
    if (pLocalized != 0) {
      if (d_allowNoTension) {
        if (pressure > 0.0){
            pStress*=D;
        } else{
            pStress = Identity*pressure;
        }
      } else if (d_allowNoShear){
         pStress = Identity*pressure;
      }
      else if (d_setStressToZero){
        pStress*=D;
      }
    }
  }
}

// This is for the localization flags to be updated
void 
BasicDamageModel::addRequiresLocalizationParameter(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet* ) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
}

// This is for the localization flag to be updated
void 
BasicDamageModel::getLocalizationParameter(const Patch* patch, 
                                           ParticleVariable<int>& isLocalized, int dwi,
                                           DataWarehouse* old_dw,
                                           DataWarehouse* new_dw)
{
  ParticleSubset* pset = old_dw->getParticleSubset(dwi,patch);
  constParticleVariable<int> pLocalized;
  new_dw->get(pLocalized, pLocalizedLabel_preReloc, pset);

  ParticleSubset::iterator iter;
  for (iter = pset->begin(); iter != pset->end(); iter++) {
    isLocalized[*iter] = pLocalized[*iter];
  }
}
