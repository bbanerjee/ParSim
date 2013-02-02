/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
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

#include <CCA/Components/MPM/ConstitutiveModel/ConstitutiveModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>
#include <CCA/Components/MPM/MPMFlags.h>
#include <Core/Math/Matrix3.h>
#include <CCA/Ports/DataWarehouse.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Labels/MPMLabel.h>
#include <Core/Math/FastMatrix.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Malloc/Allocator.h>
#include <cmath>
#include <iostream>

using namespace Uintah;

#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

ConstitutiveModel::ConstitutiveModel(MPMFlags* Mflag)
{
  lb = scinew MPMLabel();
  flag = Mflag;
  if(flag->d_8or27==8){
    NGN=1;
  } else{ 
    NGN=2;
  }
}

ConstitutiveModel::ConstitutiveModel(const ConstitutiveModel* cm)
{
  lb = scinew MPMLabel();
  flag = cm->flag;
  NGN = cm->NGN;
  NGP = cm->NGP;
  d_sharedState = cm->d_sharedState;
}

ConstitutiveModel::~ConstitutiveModel()
{
  delete lb;
}

void 
ConstitutiveModel::addInitialComputesAndRequires(Task* ,
                                                 const MPMMaterial* ,
                                                 const PatchSet*) const
{
}

///////////////////////////////////////////////////////////////////////
/*! Initialize the common quantities that all the explicit constituive
 *  models compute */
///////////////////////////////////////////////////////////////////////
void 
ConstitutiveModel::initSharedDataForExplicit(const Patch* patch,
                                             const MPMMaterial* matl,
                                             DataWarehouse* new_dw)
{
  Matrix3 I; I.Identity();
  Matrix3 zero(0.);
  ParticleSubset* pset = new_dw->getParticleSubset(matl->getDWIndex(), patch);

  ParticleVariable<double>  pdTdt;
  ParticleVariable<Matrix3> pStress;

  ParticleVariable<Matrix3> pDefGrad;
  new_dw->getModifiable(pDefGrad, lb->pDefGradLabel, pset);

  new_dw->allocateAndPut(pdTdt,   lb->pdTdtLabel,    pset);
  new_dw->allocateAndPut(pStress, lb->pStressLabel,  pset);

  // To fix : For a material that is initially stressed we need to
  // modify the stress tensors to comply with the initial stress state
  ParticleSubset::iterator iter = pset->begin();
  for(; iter != pset->end(); iter++){
    particleIndex idx = *iter;
    pdTdt[idx] = 0.0;
    pDefGrad[idx] = I;
    pStress[idx] = zero;
  }
}

void 
ConstitutiveModel::addComputesAndRequires(Task*, 
                                          const MPMMaterial*,
                                          const PatchSet*) const
{
  throw InternalError("Stub Task: ConstitutiveModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
ConstitutiveModel::addComputesAndRequires(Task*, 
                                          const MPMMaterial*,
                                          const PatchSet*,
                                          const bool ,
                                          const bool) const
{
  throw InternalError("Stub Task: ConstitutiveModel::addComputesAndRequires ", __FILE__, __LINE__);  
}

void ConstitutiveModel::scheduleCheckNeedAddMPMMaterial(Task* task, 
                                                        const MPMMaterial*,
                                                        const PatchSet*) const
{
  task->computes(lb->NeedAddMPMMaterialLabel);
}

void 
ConstitutiveModel::addSharedCRForHypoExplicit(Task* task,
                                              const MaterialSubset* matlset,
                                              const PatchSet* p) const
{
  Ghost::GhostType  gnone = Ghost::None;
  addSharedCRForExplicit(task, matlset, p);
  task->requires(Task::OldDW, lb->pStressLabel,             matlset, gnone);
}

void 
ConstitutiveModel::addSharedCRForExplicit(Task* task,
                                          const MaterialSubset* matlset,
                                          const PatchSet* ) const
{
  Ghost::GhostType  gnone = Ghost::None;
  Ghost::GhostType  gac   = Ghost::AroundCells;

  task->requires(Task::OldDW, lb->delTLabel);
  task->requires(Task::OldDW, lb->pXLabel,                  matlset, gnone);
  task->requires(Task::OldDW, lb->pMassLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pVolumeLabel,             matlset, gnone);
  task->requires(Task::OldDW, lb->pTemperatureLabel,        matlset, gnone);
  task->requires(Task::OldDW, lb->pVelocityLabel,           matlset, gnone);
  task->requires(Task::OldDW, lb->pDefGradLabel,            matlset, gnone);
  task->requires(Task::OldDW, lb->pVelGradLabel,            matlset, gnone);
  task->requires(Task::NewDW, lb->pDefGradLabel_preReloc,   matlset, gnone);
  task->requires(Task::NewDW, lb->pVelGradLabel_preReloc,   matlset, gnone);
  //task->requires(Task::OldDW, lb->pDeformationMeasureLabel, matlset, gnone);
  task->requires(Task::NewDW, lb->gVelocityStarLabel,       matlset, gac, NGN);
  if(!flag->d_doGridReset){
    task->requires(Task::NewDW, lb->gDisplacementLabel,     matlset, gac, NGN);
  }
  task->requires(Task::OldDW, lb->pSizeLabel,               matlset, gnone);
  task->requires(Task::OldDW, lb->pTempPreviousLabel,       matlset, gnone);
  if (flag->d_fracture) {
    task->requires(Task::NewDW, lb->pgCodeLabel,            matlset, gnone); 
    task->requires(Task::NewDW, lb->GVelocityStarLabel,     matlset, gac, NGN);
  }

  task->computes(lb->pStressLabel_preReloc,             matlset);
  //task->computes(lb->pDeformationMeasureLabel_preReloc, matlset);
  //task->computes(lb->pVolumeLabel_preReloc,             matlset);
  task->modifies(lb->pVolumeLabel_preReloc,             matlset);
  task->computes(lb->pdTdtLabel_preReloc,               matlset);
  //task->computes(lb->p_qLabel_preReloc,                 matlset);
}

void 
ConstitutiveModel::computeStressTensor(const PatchSubset*,
                                       const MPMMaterial*,
                                       DataWarehouse*,
                                       DataWarehouse*)
{
  throw InternalError("Stub Task: ConstitutiveModel::computeStressTensor ", __FILE__, __LINE__);
}

void 
ConstitutiveModel::computeStressTensorImplicit(const PatchSubset*,
                                               const MPMMaterial*,
                                               DataWarehouse*,
                                               DataWarehouse*)
{
  throw InternalError("Stub Task: ConstitutiveModel::computeStressTensorImplicit ", __FILE__, __LINE__);
}

void ConstitutiveModel::checkNeedAddMPMMaterial(const PatchSubset*,
                                                const MPMMaterial*,
                                                DataWarehouse* new_dw,
                                                DataWarehouse*)
{
  double need_add=0.;
                                                                                
  new_dw->put(sum_vartype(need_add),     lb->NeedAddMPMMaterialLabel);
}


void 
ConstitutiveModel::carryForward(const PatchSubset*,
                                const MPMMaterial*,
                                DataWarehouse*,
                                DataWarehouse*)
{
  throw InternalError("Stub Task: ConstitutiveModel::carryForward ", __FILE__, __LINE__);
}

void
ConstitutiveModel::carryForwardSharedData(ParticleSubset* pset,
                                          DataWarehouse*  old_dw,
                                          DataWarehouse*  new_dw,
                                          const MPMMaterial* matl)
{
  double rho_orig = matl->getInitialDensity();
  Matrix3 Id, Zero(0.0); Id.Identity();

  constParticleVariable<double>  pMass;
  constParticleVariable<Matrix3> pDefGrad_old;
  old_dw->get(pMass,            lb->pMassLabel,               pset);
  old_dw->get(pDefGrad_old,     lb->pDefGradLabel,            pset);

  ParticleVariable<double>  pVol_new;
  ParticleVariable<Matrix3> pDefGrad_new;
  new_dw->getModifiable(pVol_new,     lb->pVolumeLabel_preReloc,  pset);
  new_dw->getModifiable(pDefGrad_new, lb->pDefGradLabel_preReloc, pset);

  ParticleVariable<double>  pIntHeatRate_new,p_q;
  ParticleVariable<Matrix3> pStress_new;
  new_dw->allocateAndPut(pIntHeatRate_new, lb->pdTdtLabel_preReloc,    pset);
  new_dw->allocateAndPut(pStress_new,   lb->pStressLabel_preReloc,     pset);
  new_dw->allocateAndPut(p_q,           lb->p_qLabel_preReloc,         pset);

  ParticleSubset::iterator iter = pset->begin();
  for(; iter != pset->end(); iter++){
    particleIndex idx = *iter;
    pVol_new[idx] = (pMass[idx]/rho_orig);
    pIntHeatRate_new[idx] = 0.0;
    pDefGrad_new[idx] = pDefGrad_old[idx];
    pStress_new[idx] = Zero;
    p_q[idx]=0.;
  }
}

void 
ConstitutiveModel::allocateCMDataAddRequires(Task*, const MPMMaterial*,
                                             const PatchSet*,
                                             MPMLabel*) const
{
  throw InternalError("Stub Task: ConstitutiveModel::allocateCMDataAddRequires ", __FILE__, __LINE__);
}

void 
ConstitutiveModel::addSharedRForConvertExplicit(Task* task,
                                                const MaterialSubset* mset,
                                                const PatchSet*) const
{
  Ghost::GhostType  gnone = Ghost::None;
  task->requires(Task::NewDW,lb->pdTdtLabel_preReloc,              mset,gnone);
  //task->requires(Task::NewDW,lb->pDeformationMeasureLabel_preReloc,mset,gnone);
  task->requires(Task::NewDW,lb->pStressLabel_preReloc,            mset,gnone);
}

void
ConstitutiveModel::copyDelToAddSetForConvertExplicit(DataWarehouse* new_dw,
                                                     ParticleSubset* delset,
                                                     ParticleSubset* addset,
                                                     map<const VarLabel*, ParticleVariableBase*>* newState)
{
  constParticleVariable<double>  pIntHeatRate_del;
  //constParticleVariable<Matrix3> pDefGrad_del;
  constParticleVariable<Matrix3> pStress_del;

  new_dw->get(pIntHeatRate_del, lb->pdTdtLabel_preReloc,               delset);
  //new_dw->get(pDefGrad_del,     lb->pDeformationMeasureLabel_preReloc, delset);
  new_dw->get(pStress_del,      lb->pStressLabel_preReloc,             delset);

  ParticleVariable<double>  pIntHeatRate_add;
  //ParticleVariable<Matrix3> pDefGrad_add;
  ParticleVariable<Matrix3> pStress_add;

  new_dw->allocateTemporary(pIntHeatRate_add, addset);
  //new_dw->allocateTemporary(pDefGrad_add,     addset);
  new_dw->allocateTemporary(pStress_add,      addset);

  ParticleSubset::iterator del = delset->begin();
  ParticleSubset::iterator add = addset->begin();
  for (; del != delset->end(); del++, add++) {
    pIntHeatRate_add[*add] = pIntHeatRate_del[*del];
    //pDefGrad_add[*add] = pDefGrad_del[*del];
    pStress_add[*add]  = pStress_del[*del];
  }

  (*newState)[lb->pdTdtLabel] = pIntHeatRate_add.clone();
  //(*newState)[lb->pDeformationMeasureLabel] = pDefGrad_add.clone();
  (*newState)[lb->pStressLabel] = pStress_add.clone();
}

void 
ConstitutiveModel::addRequiresDamageParameter(Task*, 
                                              const MPMMaterial*,
                                              const PatchSet*) const
{
}

void 
ConstitutiveModel::getDamageParameter(const Patch* ,
                                      ParticleVariable<int>& ,int ,
                                      DataWarehouse* ,
                                      DataWarehouse* )
{
}

Vector 
ConstitutiveModel::getInitialFiberDir()
{
  return Vector(0.,0.,1);
}

//______________________________________________________________________
//______________________________________________________________________
//          HARDWIRE FOR AN IDEAL GAS -Todd 
double 
ConstitutiveModel::computeRhoMicro(double press, double gamma,
                                   double cv, double Temp, double rho_guess)
{
  // Pointwise computation of microscopic density
  return  press/((gamma - 1.0)*cv*Temp);
}

void 
ConstitutiveModel::computePressEOS(double rhoM, double gamma,
                                   double cv, double Temp, double& press, 
                                   double& dp_drho, double& dp_de)
{
  // Pointwise computation of thermodynamic quantities
  press   = (gamma - 1.0)*rhoM*cv*Temp;
  dp_drho = (gamma - 1.0)*cv*Temp;
  dp_de   = (gamma - 1.0)*rhoM;
}
//______________________________________________________________________


// Convert J-integral into stress intensity (for FRACTURE)
void 
ConstitutiveModel::ConvertJToK(const MPMMaterial*,
                               const string&,
                               const Vector&,
                               const double&,
                               const Vector&,
                               Vector& SIF)
{
  SIF=Vector(-9999.,-9999.,-9999.);
}

// Detect if crack propagtes and the propagation direction (for FRACTURE)
short
ConstitutiveModel::CrackPropagates(const double& , const double& , 
                                   const double& , double& theta)
{
  enum {NO=0, YES};
  theta=0.0;
  return NO;
}

double 
ConstitutiveModel::artificialBulkViscosity(double Dkk, 
                                           double c_bulk, 
                                           double rho,
                                           double dx) const 
{
  double q = 0.0;
  if (Dkk < 0.0) {
    double A1 = flag->d_artificialViscCoeff1;
    double A2 = flag->d_artificialViscCoeff2;
    //double c_bulk = sqrt(K/rho);
    q = (A1*fabs(c_bulk*Dkk*dx) + A2*(Dkk*Dkk*dx*dx))*rho;
  }
  return q;
}

void
ConstitutiveModel::computeDeformationGradientFromDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3> &Fold,
                                           Vector dx,
                                           ParticleInterpolator* interp) 
{
  Matrix3 dispGrad,Identity;
  Identity.Identity();
  vector<IntVector> ni(interp->size());
  vector<Vector> d_S(interp->size());
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                            
  for(ParticleSubset::iterator iter = pset->begin();
       iter != pset->end(); iter++){
    particleIndex idx = *iter;
                                                                            
    // Get the node indices that surround the cell
    interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                            
    computeGrad(dispGrad, ni, d_S, oodx, gDisp);

    // Update the deformation gradient tensor to its time n+1 value.
    // Compute the deformation gradient from the displacement gradient
    Fnew[idx] = Identity + dispGrad;

    double J = Fnew[idx].Determinant();
    if (!(J > 0)) {
      ostringstream warn;
      warn << "**ERROR** : ConstitutiveModel::computeDeformationGradientFromDisplacement" << endl << "Negative or zero determinant of Jacobian." << endl;
      warn << "     Particle = " << idx << " J = " << J << " position = " << px[idx] << endl;
      warn << "     Disp Grad = " << dispGrad << endl; 
      warn << "     F_new = " << Fnew[idx] << endl; 
      throw InvalidValue(warn.str(), __FILE__, __LINE__);
    }
  }
}

void 
ConstitutiveModel::computeDeformationGradientFromVelocity(
                                           constNCVariable<Vector> gVel,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> psize,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           ParticleInterpolator* interp,
                                           const double& delT)
{
    Matrix3 velGrad,deformationGradientInc, Identity;
    Identity.Identity();
    vector<IntVector> ni(interp->size());
    vector<Vector> d_S(interp->size());
    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};

    for(ParticleSubset::iterator iter = pset->begin();
      iter != pset->end(); iter++){
      particleIndex idx = *iter;

      // Get the node indices that surround the cell
      interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);

      computeGrad(velGrad, ni, d_S, oodx, gVel);

      // Compute the deformation gradient increment using the time_step
      // velocity gradient
      // F_n^np1 = dudx * dt + Identity
      deformationGradientInc = velGrad * delT + Identity;
                                                                              
      // Update the deformation gradient tensor to its time n+1 value.
      Fnew[idx] = deformationGradientInc * Fold[idx];

      double J = Fnew[idx].Determinant();
      if (!(J > 0)) {
        ostringstream warn;
        warn << "**ERROR** CompNeoHook: Negative or zero determinant of Jacobian."
             << " Particle has inverted." << endl;
        warn << "     Particle = " << idx << ", J = " << J << ", position = " << px[idx]<<endl;
        warn << "          Vel Grad = \n" << velGrad << endl; 
        warn << "          F_inc = \n" << deformationGradientInc << endl; 
        warn << "          F_old = \n" << Fold[idx] << endl; 
        warn << "          F_new = \n" << Fnew[idx] << endl; 
        warn << "          gVelocity:" << endl;
        for(int k = 0; k < flag->d_8or27; k++) {
          warn<< "             node: " << ni[k] << " vel: " << gVel[ni[k]] << endl;
        }
        
        throw InvalidValue(warn.str(), __FILE__, __LINE__);
      }

    }
}

void
ConstitutiveModel::computeDeformationGradientFromTotalDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           ParticleVariable<Matrix3> &Fnew,
                                           constParticleVariable<Matrix3>& Fold,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp)
{
  Matrix3 dispGrad,Identity;
  Identity.Identity();
  vector<IntVector> ni(interp->size());
  vector<double> S(interp->size());
  vector<Vector> d_S(interp->size());
  double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                                
  for(ParticleSubset::iterator iter = pset->begin();
       iter != pset->end(); iter++){
    particleIndex idx = *iter;
                                                                                
    // Get the node indices that surround the cell
    interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                                
    computeGrad(dispGrad, ni, d_S, oodx, gDisp);
                                                                                
    // Update the deformation gradient tensor to its time n+1 value.
    // Compute the deformation gradient from the displacement gradient
    Fnew[idx] = Identity + dispGrad;
  }
}
                                                                                
void
ConstitutiveModel::computeDeformationGradientFromIncrementalDisplacement(
                                           constNCVariable<Vector> gDisp,
                                           ParticleSubset* pset,
                                           constParticleVariable<Point> px,
                                           constParticleVariable<Matrix3> Fold,
                                           ParticleVariable<Matrix3> &Fnew,
                                           Vector dx,
                                           constParticleVariable<Matrix3> psize,
                                           ParticleInterpolator* interp)
{
    Matrix3 IncDispGrad,deformationGradientInc, Identity;
    Identity.Identity();
    vector<IntVector> ni(interp->size());
    vector<double> S(interp->size());
    vector<Vector> d_S(interp->size());

    double oodx[3] = {1./dx.x(), 1./dx.y(), 1./dx.z()};
                                                                                
    for(ParticleSubset::iterator iter = pset->begin();
      iter != pset->end(); iter++){
      particleIndex idx = *iter;
                                                                                
      // Get the node indices that surround the cell
      interp->findCellAndShapeDerivatives(px[idx],ni,d_S,psize[idx],Fold[idx]);
                                                                                
      computeGrad(IncDispGrad, ni, d_S, oodx, gDisp);
                                                                                
      // Compute the deformation gradient increment
      deformationGradientInc = IncDispGrad + Identity;
                                                                                
      // Update the deformation gradient tensor to its time n+1 value.
      Fnew[idx] = deformationGradientInc * Fold[idx];
    }
}

//-----------------------------------------------------------------------------------
// Damage flags and damage distribution functions that can be used by 
// all constitutive models
//-----------------------------------------------------------------------------------
void
ConstitutiveModel::getFailureModelData(ProblemSpecP& ps)
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
ConstitutiveModel::setFailureModelData(ConstitutiveModel* cm)
{
  if (flag->d_erosionAlgorithm == "BrittleDamage") {
    setBrittleDamageData(cm);
  } else {
    // Set the failure strain data
    setFailureStressOrStrainData(cm);
    d_failure_criteria = cm->d_failure_criteria;
    if(d_failure_criteria=="MohrColoumb"){
      d_tensile_cutoff = cm->d_tensile_cutoff;
      d_friction_angle = cm->d_friction_angle;
    }
  }

  // Set the erosion algorithm
  setErosionAlgorithm(cm);
}

void 
ConstitutiveModel::getFailureStressOrStrainData(ProblemSpecP& ps)
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
ConstitutiveModel::outputProblemSpecDamage(ProblemSpecP& cm_ps)
{
  if (flag->d_erosionAlgorithm == "BrittleDamage") {
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

void 
ConstitutiveModel::setFailureStressOrStrainData(const UCNH* cm)
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
ConstitutiveModel::setBrittleDamageData(const ConstitutiveModel* cm)
{
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
ConstitutiveModel::getBrittleDamageData(ProblemSpecP& ps)
{
  d_brittle_damage.r0b   = 57.0; // Initial energy threshold
  d_brittle_damage.Gf    = 11.2; // Fracture energy
  d_brittle_damage.constant_D = 0.1; // Shape constant in softening function 
  d_brittle_damage.maxDamageInc=0.1; // Maximum damage in a time step
  d_brittle_damage.allowRecovery=false; // Allow recovery
  d_brittle_damage.recoveryCoeff=1.0; // Fraction of recovery if allowed
  d_brittle_damage.printDamage=false;  // Print damage
  ps->get("brittle_damage_initial_threshold",   d_brittle_damage.r0b);
  ps->get("brittle_damage_fracture_energy",     d_brittle_damage.Gf);
  ps->get("brittle_damage_constant_D",          d_brittle_damage.constant_D);
  ps->get("brittle_damage_max_damage_increment",d_brittle_damage.maxDamageInc);
  ps->get("brittle_damage_allowRecovery",       d_brittle_damage.allowRecovery);
  ps->get("brittle_damage_recoveryCoeff",       d_brittle_damage.recoveryCoeff);
  ps->get("brittle_damage_printDamage",         d_brittle_damage.printDamage);
  if (d_brittle_damage.recoveryCoeff <0.0 || d_brittle_damage.recoveryCoeff>1.0)
  {
    cerr << "brittle_damage_recoveryCoeff must be between 0.0 and 1.0" << endl;
  }     
}

void 
ConstitutiveModel::setErosionAlgorithm()
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
ConstitutiveModel::setErosionAlgorithm(const ConstitutiveModel* cm)
{
  d_setStressToZero = cm->d_setStressToZero;
  d_allowNoTension  = cm->d_allowNoTension;
  d_allowNoShear    = cm->d_allowNoShear;
  d_brittleDamage   = cm->d_brittleDamage;
}

void
ConstitutiveModel::initializeDamageVarLabels()
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

void
ConstitutiveModel::deleteDamageVarLabels()
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

void
ConstitutiveModel::initializeDamageData(const Patch* patch,
                                        const MPMMaterial* matl,
                                        DataWarehouse* new_dw)
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

void 
ConstitutiveModel::copyDamageDataFromDeletedToAddedParticle(DataWarehouse* new_dw,
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

void 
ConstitutiveModel::allocateDamageDataAddRequires(Task* task,
                                                 const MPMMaterial* matl,
                                                 const PatchSet* patches,
                                                 MPMLabel*lb ) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  Ghost::GhostType  gnone = Ghost::None;
  task->requires(Task::NewDW, pFailureStressOrStrainLabel_preReloc, matlset,
                                                                       gnone);
  task->requires(Task::NewDW, pLocalizedLabel_preReloc,     matlset, gnone);
  task->requires(Task::NewDW, pTimeOfLocLabel_preReloc,     matlset, gnone);
  task->requires(Task::NewDW, pDamageLabel_preReloc,        matlset, gnone);
}

void 
ConstitutiveModel::carryForwardDamageData(ParticleSubset* pset,
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

void 
ConstitutiveModel::addComputesAndRequiresForDamage(Task* task,
                                                   const MPMMaterial* matl,
                                                   const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  
  task->requires(Task::OldDW, pFailureStressOrStrainLabel,    matlset, gnone);
  task->requires(Task::OldDW, pLocalizedLabel,                matlset, gnone);
  task->requires(Task::OldDW, pTimeOfLocLabel,                matlset, gnone);
  task->requires(Task::OldDW, pDamageLabel,                   matlset, gnone);
    
  task->computes(pFailureStressOrStrainLabel_preReloc,        matlset);
  task->computes(pLocalizedLabel_preReloc,                    matlset);
  task->computes(pTimeOfLocLabel_preReloc,                    matlset);
  task->computes(pDamageLabel_preReloc,                       matlset);
  task->computes(lb->TotalLocalizedParticleLabel);   
}

void 
ConstitutiveModel::addInitialComputesAndRequiresForDamage(Task* task,
                                                          const MPMMaterial* matl,
                                                          const PatchSet* patches) const
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pFailureStressOrStrainLabel, matlset);
  task->computes(pLocalizedLabel,             matlset);
  task->computes(pTimeOfLocLabel,             matlset);
  task->computes(pDamageLabel,                matlset);
  task->computes(lb->TotalLocalizedParticleLabel);
}  

void 
ConstitutiveModel::addRequiresDamageParameterDefault(Task* task,
                                                     conts MPMMaterial* matl,
                                                     const PatchSet* patches)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pLocalizedLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pTimeOfLocLabel_preReloc, matlset, Ghost::None);
  task->requires(Task::NewDW, pDamageLabel_preReloc, matlset, Ghost::None);
}

void
ConstitutiveModel::addParticleStateDamage(std::vector<const VarLabel*>& from,
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

// Modify the stress for brittle damage
// Update pFailureStrain_new (energy threshold)
// pDamage_new (damage; if negative, damage inactive), pStress
void 
ConstitutiveModel::updateDamageAndModifyStress(const Matrix3& defGrad, 
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
      
       if (d_brittle_damage.printDamage) cout << "Particle " << particleID << " damage halted: damage=" << pDamage_new << endl;
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
      double young = 9.0*d_initialData.Bulk*d_initialData.tauDev/\
        (3.0*d_initialData.Bulk+d_initialData.tauDev);

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
          cout << "Particle " << particleID << " damaged: "
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
            cout << "Particle " << particleID << " damage halted: damage=" 
                 << pDamage_new << endl;
          }
	}
	else { //no recovery (default)
	  pStress = pStress*(1.0-pDamage);
          if (d_brittle_damage.printDamage){
            cout << "Particle " << particleID << " damaged: " 
                 << " damage=" << pDamage_new << " epsMax=" << epsMax 
                 << " tau_b=" << tau_b << endl;
          }
	}
      } // end if tau_b > pFailureStrain

  } //end if pressure

}

// Modify the stress if particle has failed
void 
ConstitutiveModel::updateFailedParticlesAndModifyStress(const Matrix3& defGrad,
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
        cout << "Particle " << particleID << " has failed : MaxPrinStress = "
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
        cout << "Particle " << particleID << " has failed : eps = " << maxEigen
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
        cout << "Particle " << particleID << " has failed : maxPrinStress = "
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

