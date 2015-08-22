/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/InternalVar_MasonSand.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_MasonSand.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/InternalError.h>

#include <errno.h>
#include <fenv.h>

using namespace Vaango;
using namespace Uintah;
using namespace std;


InternalVar_MasonSand::InternalVar_MasonSand(ProblemSpecP& ps)
{
  ps->require("p0", d_crushParam.p0);  // Crush Curve Parameter
  ps->require("p1", d_crushParam.p1);  // Crush Curve Parameter
  ps->require("p2", d_crushParam.p2);  // Crush Curve Parameter (not used)
  ps->require("p3", d_crushParam.p3);  // Crush Curve Parameter

  // Initialize internal variable labels for evolution
  pKappaLabel = VarLabel::create("p.kappa",
        ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc = VarLabel::create("p.kappa+",
        ParticleVariable<double>::getTypeDescription());

}
         
InternalVar_MasonSand::InternalVar_MasonSand(const InternalVar_MasonSand* cm)
{
  d_crushParam = cm->d_crushParam;

  // Initialize internal variable labels for evolution
  pKappaLabel = VarLabel::create("p.kappa",
        ParticleVariable<double>::getTypeDescription());
  pKappaLabel_preReloc = VarLabel::create("p.kappa+",
        ParticleVariable<double>::getTypeDescription());

}
         
InternalVar_MasonSand::~InternalVar_MasonSand()
{
  VarLabel::destroy(pKappaLabel);
  VarLabel::destroy(pKappaLabel_preReloc);
}


void InternalVar_MasonSand::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP int_var_ps = ps->appendChild("internal_variable_model");
  int_var_ps->setAttribute("type","mason_sand");

  int_var_ps->appendElement("p0", d_crushParam.p0);
  int_var_ps->appendElement("p1", d_crushParam.p1);
  int_var_ps->appendElement("p2", d_crushParam.p2);
  int_var_ps->appendElement("p3", d_crushParam.p3);
}

         
void 
InternalVar_MasonSand::addInitialComputesAndRequires(Task* task,
                                                     const MPMMaterial* matl ,
                                                     const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->computes(pKappaLabel, matlset);
}

void 
InternalVar_MasonSand::addComputesAndRequires(Task* task,
                                              const MPMMaterial* matl ,
                                              const PatchSet*)
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::OldDW, pKappaLabel, matlset,Ghost::None);
  task->computes(pKappaLabel_preReloc, matlset);
}

void 
InternalVar_MasonSand::addParticleState(std::vector<const VarLabel*>& from,
                                        std::vector<const VarLabel*>& to)
{
  from.push_back(pKappaLabel);
  to.push_back(pKappaLabel_preReloc);
}

void 
InternalVar_MasonSand::allocateCMDataAddRequires(Task* task,
                                                 const MPMMaterial* matl ,
                                                 const PatchSet* ,
                                                 MPMLabel* )
{
  const MaterialSubset* matlset = matl->thisMaterial();
  task->requires(Task::NewDW, pKappaLabel_preReloc, matlset, Ghost::None);
}

void 
InternalVar_MasonSand::allocateCMDataAdd(DataWarehouse* old_dw,
                                         ParticleSubset* addset,
                                             map<const VarLabel*, 
                                               ParticleVariableBase*>* newState,
                                             ParticleSubset* delset,
                                             DataWarehouse* new_dw )
{
  ParticleVariable<double> pKappa;
  constParticleVariable<double> o_kappa;

  new_dw->allocateTemporary(pKappa,addset);

  new_dw->get(o_kappa,pKappaLabel_preReloc,delset);

  ParticleSubset::iterator o,n = addset->begin();
  for(o = delset->begin(); o != delset->end(); o++, n++) {
    pKappa[*n] = o_kappa[*o];
  }

  (*newState)[pKappaLabel]=pKappa.clone();

}

void 
InternalVar_MasonSand::initializeInternalVariable(ParticleSubset* pset,
                                                      DataWarehouse* new_dw)
{
  Uintah::ParticleVariable<double> pKappa;
  new_dw->allocateAndPut(pKappa, pKappaLabel, pset);
  ParticleSubset::iterator iter = pset->begin();
  for(;iter != pset->end(); iter++) {
    pKappa[*iter] = 0.0;
  }
}

void 
InternalVar_MasonSand::getInternalVariable(ParticleSubset* pset ,
                                               DataWarehouse* old_dw,
                                               constParticleVariableBase& pKappa) 
{
  old_dw->get(pKappa, pKappaLabel, pset);
}

void 
InternalVar_MasonSand::allocateAndPutInternalVariable(ParticleSubset* pset,
                                                          DataWarehouse* new_dw,
                                                          ParticleVariableBase& pKappa_new) 
{
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc, pset);
}

void
InternalVar_MasonSand::allocateAndPutRigid(ParticleSubset* pset,
                                               DataWarehouse* new_dw,
                                               constParticleVariableBase& pKappa)
{
  ParticleVariable<double> pKappa_new;
  new_dw->allocateAndPut(pKappa_new, pKappaLabel_preReloc, pset);
  ParticleSubset::iterator iter = pset->begin();
  for(;iter != pset->end(); iter++){
     pKappa_new[*iter] = dynamic_cast<constParticleVariable<double>& >(pKappa)[*iter];
  }
}

//--------------------------------------------------------------------------------------
// Compute kappa_new using Newton's method
//--------------------------------------------------------------------------------------
double 
InternalVar_MasonSand::computeInternalVariable(const ModelStateBase* state_input) const
{
  const ModelState_MasonSand* state = dynamic_cast<const ModelState_MasonSand*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_MasonSand.";
    throw InternalError(out.str(), __FILE__, __LINE__);
  }

  // return the new kappa
  double kappa_new = 0.0;
  return kappa_new;
}

