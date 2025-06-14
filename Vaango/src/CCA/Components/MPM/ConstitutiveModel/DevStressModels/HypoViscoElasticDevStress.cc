/*
 * The MIT License
 *
 * Copyright (c) 1997-2012 The University of Utah
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include "HypoViscoElasticDevStress.h"
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/Constants.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Util/DebugStream.h>


using namespace Uintah;
static DebugStream dbg("HypoViscoElasticDevStress", false);

HypoViscoElasticDevStress::HypoViscoElasticDevStress(ProblemSpecP& ps)
{
  ps->get("tau", d_tau_MW);
  ps->get("mu", d_mu_MW);

  // bulletproofing
  if (d_tau_MW.size() != d_mu_MW.size()) {
    throw ProblemSetupException(
      "ERROR:  The number of maxwell elements for tau != mu", __FILE__,
      __LINE__);
  }

  if (d_tau_MW.size() == 0 || d_mu_MW.size() == 0) {
    throw ProblemSetupException(
      "ERROR:  The number of maxwell elements for tau & mu > 0", __FILE__,
      __LINE__);
  }

  // for speed:
  for (double j : d_tau_MW) {
    d_inv_tau_MW.push_back(1.0 / j);
  }

  // number of Maxwell Elements
  d_MaxwellElements = d_tau_MW.size();

  // create labels for each maxwell element
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
     std::ostringstream name, name2;
    name << "sigmaDev" << j;
    name2 << name.str() << "+";

    // create internal variable labels
    d_sigmaDevLabel.push_back(VarLabel::create(
      name.str(), ParticleVariable<Matrix3>::getTypeDescription()));
    d_sigmaDevLabel_preReloc.push_back(VarLabel::create(
      name2.str(), ParticleVariable<Matrix3>::getTypeDescription()));
  }

  // This is a std::vector of ParticleVariable<Matrix3>
  // Each particle needs d_MaxwellElements matrices
  d_sigmaDev.resize(d_MaxwellElements);
  d_sigmaDev_new.resize(d_MaxwellElements);

  // create the arrays
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    d_sigmaDev[j] = scinew constParticleVariable<Matrix3>();
    d_sigmaDev_new[j] = scinew ParticleVariable<Matrix3>();
  }
}
//______________________________________________________________________
//
HypoViscoElasticDevStress::~HypoViscoElasticDevStress()
{
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    VarLabel::destroy(d_sigmaDevLabel[j]);
    VarLabel::destroy(d_sigmaDevLabel_preReloc[j]);
    delete d_sigmaDev[j];
    delete d_sigmaDev_new[j];
  }
}
//______________________________________________________________________
//
void
HypoViscoElasticDevStress::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP flow_ps = ps->appendChild("deviatoric_stress_model");
  flow_ps->setAttribute("type", "hypoViscoElastic");
  flow_ps->appendElement("tau", d_tau_MW);
  flow_ps->appendElement("mu", d_mu_MW);
}
//______________________________________________________________________
//
void
HypoViscoElasticDevStress::addInitialComputesAndRequires(
  Task* task, const MPMMaterial* matl)
{
  dbg << " hypoViscoElastic::addInitialComputesAndRequires " << std::endl;
  const MaterialSubset* matlset = matl->thisMaterial();

  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    task->computes(d_sigmaDevLabel[j], matlset);
  }
}
//______________________________________________________________________
//    Called by computeStressTensor()
void
HypoViscoElasticDevStress::addComputesAndRequires(Task* task,
                                                  const MPMMaterial* matl)
{
  dbg << " hypoViscoElastic:addComputesAndRequires 1 " << std::endl;
  const MaterialSubset* matlset = matl->thisMaterial();
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {

    task->needs(Task::OldDW, d_sigmaDevLabel[j], matlset, Ghost::None);
    task->computes(d_sigmaDevLabel_preReloc[j], matlset);
  }
}
//______________________________________________________________________
//    Called by computeStressTensorImplicit
void
HypoViscoElasticDevStress::addComputesAndRequires(Task* task,
                                                  const MPMMaterial* matl,
                                                  bool SchedParent)
{

  dbg << " hypoViscoElastic:addComputesAndRequires 2 " << std::endl;
  const MaterialSubset* matlset = matl->thisMaterial();
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    if (SchedParent) {
      task->needs(Task::ParentOldDW, d_sigmaDevLabel[j], matlset,
                     Ghost::None);
    } else {
      task->needs(Task::OldDW, d_sigmaDevLabel[j], matlset, Ghost::None);
    }
  }
}
//______________________________________________________________________
//
void
HypoViscoElasticDevStress::addParticleState(std::vector<const VarLabel*>& from,
                                            std::vector<const VarLabel*>& to)
{
  dbg << " hypoViscoElastic:addParticleState " << std::endl;
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    from.push_back(d_sigmaDevLabel[j]);
    to.push_back(d_sigmaDevLabel_preReloc[j]);
  }
}

//______________________________________________________________________
//
void
HypoViscoElasticDevStress::initializeInternalVars(ParticleSubset* pset,
                                                  DataWarehouse* new_dw)
{
  dbg << " hypoViscoElastic:initializeInternalVars " << std::endl;
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    new_dw->allocateAndPut(*d_sigmaDev_new[j], d_sigmaDevLabel[j], pset);
    for (auto idx : *pset) {
      (*d_sigmaDev_new[j])[idx] = Vaango::Util::Zero;
    }
  }
}
//______________________________________________________________________
//
void
HypoViscoElasticDevStress::getInternalVars(ParticleSubset* pset,
                                           DataWarehouse* old_dw)
{
  dbg << " hypoViscoElastic:getInternalVars " << std::endl;
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    old_dw->get(*d_sigmaDev[j], d_sigmaDevLabel[j], pset);
  }
}
//______________________________________________________________________
//
void
HypoViscoElasticDevStress::allocateAndPutInternalVars(ParticleSubset* pset,
                                                      DataWarehouse* new_dw)
{
  dbg << " hypoViscoElastic:allocateAndPutInternalVars " << std::endl;
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    new_dw->allocateAndPut(*d_sigmaDev_new[j], d_sigmaDevLabel_preReloc[j],
                           pset);
  }
}

//______________________________________________________________________
//  Initializing to zero for the sake of RigidMPM's carryForward
void
HypoViscoElasticDevStress::allocateAndPutRigid(ParticleSubset* pset,
                                               DataWarehouse* new_dw)
{
  dbg << " hypoViscoElastic:allocateAndPutRigid " << std::endl;
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    new_dw->allocateAndPut(*d_sigmaDev_new[j], d_sigmaDevLabel_preReloc[j], pset);
    for (auto idx : *pset) {
      (*d_sigmaDev_new[j])[idx] = Vaango::Util::Zero;
    }
  }
}

//______________________________________________________________________
///
void
HypoViscoElasticDevStress::computeDeviatoricStressInc(
  const particleIndex idx, [[maybe_unused]] const ModelStateBase* plaState,
  DeformationState* defState, const double delT)
{

  dbg << " hypoViscoElastic:computeDevStessInc " << std::endl;

  Matrix3 sigmadot = 0.0;

  // A solution instability was found for constant strain rate, uniaxial
  // compression.
  // For this case, the stress should saturate (sigmadot = 0.0).  It did for a
  // while,
  // then began to grow.  It's not clear what drives this, apparently numerical
  // noise?
  // Explicitly setting sigmadot to zero when it is small relative to the
  // elastic
  // stress increment solved the problem for a range of strain rates
  // investigated.

  for (unsigned int j = 0; j < d_MaxwellElements; j++) {

    Matrix3 sigmadot_elastic = 2.0 * d_mu_MW[j] * defState->devD;
    Matrix3 sigmadot_trial =
      sigmadot_elastic - (*d_sigmaDev[j])[idx] * d_inv_tau_MW[j];

    if (sigmadot_trial.Norm() > 1.e-4 * sigmadot_elastic.Norm()) {
      sigmadot += sigmadot_trial;
    }
  }

  // sigma_dev_trial = sigma_dev_n + sigmadot*delT;    (original Equation.)

  defState->devStressInc = sigmadot * delT;

  // bulletproofing
  for (unsigned int j = 0; j < d_MaxwellElements; j++) {
    if (d_tau_MW[j] < delT) {
       std::ostringstream warn;
      warn << "ERROR: hypoViscoElastic:computeDevStessInc \n"
           << "tau [" << d_tau_MW[j] << "] < delT [" << delT << "] ";
      throw InternalError("warn.str()", __FILE__, __LINE__);
    }
  }
}

//______________________________________________________________________
//
void
HypoViscoElasticDevStress::updateInternalStresses(const particleIndex idx,
                                                  const Matrix3& dp,
                                                  DeformationState* defState,
                                                  const double delT)
{
  dbg << " hypoViscoElastic:updateInternalStresses " << std::endl;
  const Matrix3 devD = defState->devD;

  int j = 0;
  double A = 0.0;
  for (auto& sigmaDev_new : d_sigmaDev_new) {
    (*sigmaDev_new)[idx] += (2.0 * d_mu_MW[j] * (devD - dp) -
                            (*sigmaDev_new)[idx] * d_inv_tau_MW[j]) * delT;
    double B = (*sigmaDev_new)[idx].NormSquared();
    A += B * d_inv_tau_MW[j] / (2.0 * d_mu_MW[j]);
    ++j;
  }

  defState->viscoElasticWorkRate = A;
}

//______________________________________________________________________
//
void
HypoViscoElasticDevStress::rotateInternalStresses(const particleIndex idx,
                                                  const Matrix3& tensorR)
{
  dbg << " hypoViscoElastic:rotateInternalStresses " << std::endl;

  //int j = 0;
  for (auto& sigmaDev_new : d_sigmaDev_new) {
    (*sigmaDev_new)[idx] =
      tensorR.Transpose() * ((*sigmaDev_new)[idx] * tensorR);
    // std::cout << "    d_sigmaDev_new " << ( *sigmaDev_new )[idx] << "
    // d_sigmaDev " << ( *d_sigmaDev[j] )[idx] << std::endl;
    //++j;
  }
}
