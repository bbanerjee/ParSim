/*
 * The MIT License
 *
 * Copyright (c) 1997-2021 The University of Utah
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


#include <CCA/Components/ICE/Materials/ICEMaterial.h>
#include <CCA/Components/ICE/Core/BoundaryCond.h>
#include <CCA/Components/Models/HEChem/JWLpp.h>
#include <CCA/Components/MPM/ConstitutiveModel/MPMMaterial.h>

#include <CCA/Ports/Scheduler.h>

#include <Core/Grid/Level.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/Variables/CellIterator.h>
#include <Core/Grid/Variables/CCVariable.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/ICE/Core/ICELabel.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <Core/Util/DebugStream.h>

#include <iostream>

using namespace Uintah;
using namespace std;
//__________________________________
//  setenv SCI_DEBUG "MODELS_DOING_COUT:+"
//  MODELS_DOING_COUT:   dumps when tasks are scheduled and performed
static DebugStream cout_doing("MODELS_DOING_COUT", false);

JWLpp::JWLpp(const ProcessorGroup* myworld, 
             const MaterialManagerP& materialManager,
             const ProblemSpecP& params,
             const ProblemSpecP& prob_spec)
  : HEChemModel(myworld, materialManager),
    d_params(params), d_prob_spec(prob_spec)
{
  mymatls = 0;
  Ilb = scinew ICELabel();
  d_saveConservedVars = scinew saveConservedVars();
  
  //__________________________________
  //  diagnostic labels
  reactedFractionLabel   = VarLabel::create("F",
                                       CCVariable<double>::getTypeDescription());
  delFLabel              = VarLabel::create("delF",
                                       CCVariable<double>::getTypeDescription());
                     
  totalMassBurnedLabel  = VarLabel::create( "totalMassBurned",
                                            sum_vartype::getTypeDescription() );
  totalHeatReleasedLabel= VarLabel::create( "totalHeatReleased",
                                            sum_vartype::getTypeDescription() );
}

JWLpp::~JWLpp()
{
  delete Ilb;
  delete d_saveConservedVars;

  VarLabel::destroy(reactedFractionLabel);
  VarLabel::destroy(delFLabel);
  VarLabel::destroy(totalMassBurnedLabel);
  VarLabel::destroy(totalHeatReleasedLabel);
    
  if(mymatls && mymatls->removeReference())
    delete mymatls;
}
//______________________________________________________________________
//
void JWLpp::problemSetup(GridP&,
                         [[maybe_unused]]  const bool isRestart)
{
  ProblemSpecP JWL_ps = d_params->findBlock("JWLpp");

  JWL_ps->getWithDefault("ThresholdVolFrac",d_threshold_volFrac, 0.01);
  
  JWL_ps->require("ThresholdPressure", d_threshold_pressure);
  JWL_ps->require("fromMaterial",fromMaterial);
  JWL_ps->require("toMaterial",  toMaterial);
  JWL_ps->require("G",    d_G);
  JWL_ps->require("b",    d_b);
  JWL_ps->require("E0",   d_E0);
  JWL_ps->require("rho0", d_rho0);

  //__________________________________
  //  Are we saving the total burned mass and total burned energy
  ProblemSpecP DA_ps = d_prob_spec->findBlock("DataArchiver");
  for (ProblemSpecP child = DA_ps->findBlock("save"); child != nullptr; child = child->findNextBlock("save") ){
    map<string,string> var_attr;
    child->getAttributes(var_attr);
    
    if (var_attr["label"] == "totalMassBurned"){
      d_saveConservedVars->mass  = true;
    }
    if (var_attr["label"] == "totalHeatReleased"){
      d_saveConservedVars->energy = true;
    }
  }
  
  matl0 = d_materialManager->parseAndLookupMaterial(JWL_ps, "fromMaterial");
  matl1 = d_materialManager->parseAndLookupMaterial(JWL_ps, "toMaterial");

  //__________________________________
  //  define the materialSet
  mymatls = scinew MaterialSet();

  vector<int> m;
  m.push_back(0);                       // needed for the pressure and NC_CCWeight
  m.push_back(matl0->getDWIndex());
  m.push_back(matl1->getDWIndex());

  mymatls->addAll_unique(m);            // elimiate duplicate entries
  mymatls->addReference(); 
}
//______________________________________________________________________
//
void JWLpp::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP model_ps = ps->appendChild("Model");
  model_ps->setAttribute("type","JWLpp");
  ProblemSpecP JWL_ps = model_ps->appendChild("JWLpp");

  JWL_ps->appendElement("ThresholdVolFrac",  d_threshold_volFrac);
  JWL_ps->appendElement("ThresholdPressure", d_threshold_pressure);
  JWL_ps->appendElement("fromMaterial",     fromMaterial);
  JWL_ps->appendElement("toMaterial",       toMaterial);
  JWL_ps->appendElement("G",    d_G);
  JWL_ps->appendElement("b",    d_b);
  JWL_ps->appendElement("E0",   d_E0);
  JWL_ps->appendElement("rho0", d_rho0);
  
}

//______________________________________________________________________
//     
void JWLpp::scheduleInitialize(SchedulerP&,
                               const LevelP&)
{
  // None necessary...
}
//______________________________________________________________________
//      
void JWLpp::scheduleComputeStableTimestep(SchedulerP&,
                                          const LevelP&)
{
  // None necessary...
}

//______________________________________________________________________
//     
void JWLpp::scheduleComputeModelSources(SchedulerP& sched,
                                       const LevelP& level)
{
  Task* t = scinew Task("JWLpp::computeModelSources", this, 
                        &JWLpp::computeModelSources);
  cout_doing << "JWLpp::scheduleComputeModelSources "<<  endl;  

  Ghost::GhostType  gn  = Ghost::None;
  const MaterialSubset* react_matl = matl0->thisMaterial();
  const MaterialSubset* prod_matl  = matl1->thisMaterial();
  MaterialSubset* one_matl     = scinew MaterialSubset();
  one_matl->add(0);
  one_matl->addReference();
  MaterialSubset* press_matl   = one_matl;

  t->needs(Task::OldDW, Ilb->timeStepLabel );
  t->needs(Task::OldDW, Ilb->delTLabel,         level.get_rep());
  //__________________________________
  // Products
  t->needs(Task::NewDW,  Ilb->rho_CCLabel,      prod_matl, gn);

  //__________________________________
  // Reactants
  t->needs(Task::NewDW, Ilb->specificVolume_CCLabel,    react_matl, gn);
  t->needs(Task::OldDW, Ilb->velocity_CCLabel,       react_matl, gn);
  t->needs(Task::OldDW, Ilb->temperature_CCLabel,      react_matl, gn);
  t->needs(Task::NewDW, Ilb->rho_CCLabel,       react_matl, gn);
  t->needs(Task::NewDW, Ilb->vol_frac_CCLabel,  react_matl, gn);

  t->needs(Task::NewDW, Ilb->press_equil_CCLabel, press_matl,gn);
  t->computes(reactedFractionLabel, react_matl);
  t->computes(delFLabel,            react_matl);

  t->modifies(Ilb->modelMass_srcLabel);
  t->modifies(Ilb->modelMom_srcLabel);
  t->modifies(Ilb->modelEng_srcLabel);
  t->modifies(Ilb->modelVol_srcLabel); 

  if(d_saveConservedVars->mass ){
    t->computes(JWLpp::totalMassBurnedLabel);
  }
  if(d_saveConservedVars->energy){
    t->computes(JWLpp::totalHeatReleasedLabel);
  } 
  sched->addTask(t, level->eachPatch(), mymatls);

  if (one_matl->removeReference())
    delete one_matl;
}
//______________________________________________________________________
//
void JWLpp::computeModelSources(const ProcessorGroup*, 
                                const PatchSubset* patches,
                                const MaterialSubset*,
                                DataWarehouse* old_dw,
                                DataWarehouse* new_dw)
{
  const Level* level = getLevel(patches);

  timeStep_vartype timeStep;
  old_dw->get(timeStep, Ilb->timeStepLabel );

  bool isNotInitialTimestep = (timeStep > 0);

  delt_vartype delT;
  old_dw->get(delT, Ilb->delTLabel, level);

  int m0 = matl0->getDWIndex(); /* reactant material */
  int m1 = matl1->getDWIndex(); /* product material */
  double totalBurnedMass = 0;
  double totalHeatReleased = 0;
 
  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);  
    
    cout_doing << "Doing computeModelSources on patch "<< patch->getID()
               <<"\t\t\t\t  JWLpp" << endl;
    CCVariable<double> mass_src_0, mass_src_1, mass_0;
    CCVariable<Vector> momentum_src_0, momentum_src_1;
    CCVariable<double> energy_src_0, energy_src_1;
    CCVariable<double> sp_vol_src_0, sp_vol_src_1;

    new_dw->getModifiable(mass_src_0,    Ilb->modelMass_srcLabel,  m0,patch);
    new_dw->getModifiable(momentum_src_0,Ilb->modelMom_srcLabel,   m0,patch);
    new_dw->getModifiable(energy_src_0,  Ilb->modelEng_srcLabel,   m0,patch);
    new_dw->getModifiable(sp_vol_src_0,  Ilb->modelVol_srcLabel,   m0,patch);

    new_dw->getModifiable(mass_src_1,    Ilb->modelMass_srcLabel,  m1,patch);
    new_dw->getModifiable(momentum_src_1,Ilb->modelMom_srcLabel,   m1,patch);
    new_dw->getModifiable(energy_src_1,  Ilb->modelEng_srcLabel,   m1,patch);
    new_dw->getModifiable(sp_vol_src_1,  Ilb->modelVol_srcLabel,   m1,patch);

    constCCVariable<double> press_CC, cv_reactant,rctVolFrac;
    constCCVariable<double> rctTemp,rctRho,rctSpvol,prodRho;
    constCCVariable<Vector> rctvel_CC;
    CCVariable<double> Fr;
    CCVariable<double> delF;
            
    Vector dx = patch->dCell();
    double cell_vol = dx.x()*dx.y()*dx.z();
    Ghost::GhostType  gn  = Ghost::None;    
   
    //__________________________________
    // Reactant data
    old_dw->get(rctTemp,       Ilb->temperature_CCLabel,      m0,patch,gn, 0); 
    old_dw->get(rctvel_CC,     Ilb->velocity_CCLabel,       m0,patch,gn, 0); 
    new_dw->get(rctRho,        Ilb->rho_CCLabel,       m0,patch,gn, 0); 
    new_dw->get(rctSpvol,      Ilb->specificVolume_CCLabel,    m0,patch,gn, 0); 
    new_dw->get(rctVolFrac,    Ilb->vol_frac_CCLabel,  m0,patch,gn, 0);
    new_dw->allocateAndPut(Fr,   reactedFractionLabel, m0,patch);
    new_dw->allocateAndPut(delF, delFLabel,            m0,patch);
    Fr.initialize(0.);
    delF.initialize(0.);

    //__________________________________
    // Product Data, 
    new_dw->get(prodRho,       Ilb->rho_CCLabel,   m1,patch,gn, 0);

    //__________________________________
    //   Misc.
    new_dw->get(press_CC,         Ilb->press_equil_CCLabel,0,  patch,gn, 0);

    // Get the specific heat, this is the value from the input file
    double cv_rct = -1.0; 
    MPMMaterial* mpm_matl = dynamic_cast<MPMMaterial *>(d_materialManager->getMaterial(m0));
    ICEMaterial* ice_matl = dynamic_cast<ICEMaterial *>(d_materialManager->getMaterial(m0));
    if(mpm_matl) {
      cv_rct = mpm_matl->getSpecificHeat();
    } else if(ice_matl){
      cv_rct = ice_matl->getSpecificHeat();  
    }
    for (CellIterator iter = patch->getCellIterator();!iter.done();iter++){
      IntVector c = *iter;
      if (press_CC[c] > d_threshold_pressure && rctVolFrac[c] > d_threshold_volFrac){          
        //__________________________________
        // Insert Burn Model Here
        double burnedMass;
        double F = prodRho[c]/(rctRho[c]+prodRho[c]);
        if(F >= 0. && F < .99){
          delF[c] = d_G*pow(press_CC[c],d_b)*(1.-F);
        }
        delF[c]*=delT;
        Fr[c] = F;
        double rctMass = rctRho[c]*cell_vol;
        double prdMass = prodRho[c]*cell_vol;
        burnedMass = min(delF[c]*(prdMass+rctMass), rctMass);
        // 20 % burned mass is a hard limit based p. 55 
        //   "JWL++: A Simple Reactive Flow Code Package for Detonation"
        burnedMass = min(burnedMass, .2*d_rho0*cell_vol);

        //__________________________________
        // conservation of mass, momentum and energy                           
        mass_src_0[c]   -= burnedMass;
        mass_src_1[c]   += burnedMass;
        totalBurnedMass += burnedMass;
           
        Vector momX        = rctvel_CC[c] * burnedMass;
        momentum_src_0[c] -= momX;
        momentum_src_1[c] += momX;

        double energyX   = cv_rct*rctTemp[c]*burnedMass; 
        double releasedHeat = burnedMass * d_E0;
        energy_src_0[c]   -= energyX;
        energy_src_1[c]   += energyX + releasedHeat;
        totalHeatReleased += releasedHeat;

        double createdVolx  = burnedMass * rctSpvol[c];
        sp_vol_src_0[c] -= createdVolx;
        sp_vol_src_1[c] += createdVolx;
      }  // if (pressure)
    }  // cell iterator  

    //__________________________________
    //  set symetric BC
    setBC(mass_src_0, "set_if_sym_BC",patch, d_materialManager, m0, new_dw, isNotInitialTimestep);
    setBC(mass_src_1, "set_if_sym_BC",patch, d_materialManager, m1, new_dw, isNotInitialTimestep);
    setBC(delF,       "set_if_sym_BC",patch, d_materialManager, m0, new_dw, isNotInitialTimestep);
    setBC(Fr,         "set_if_sym_BC",patch, d_materialManager, m0, new_dw, isNotInitialTimestep);
  }
  //__________________________________
  //save total quantities
  if(d_saveConservedVars->mass ){
    new_dw->put(sum_vartype(totalBurnedMass),   JWLpp::totalMassBurnedLabel);
  }
  if(d_saveConservedVars->energy){
    new_dw->put(sum_vartype(totalHeatReleased), JWLpp::totalHeatReleasedLabel);
  }
}
