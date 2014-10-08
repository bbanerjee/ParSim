/*
 * The MIT License
 *
 * Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

#include <CCA/Components/Peridynamics/MaterialModels/LinearElasticBondModel.h>

using namespace Vaango;

LinearElasticBondModel::LinearElasticBondModel(Uintah::ProblemSpecP& ps,
                                               PeridynamicsFlags* flags)
  : PeridynamicsMaterialModel(flags)
    
{
}

LinearElasticBondModel::LinearElasticBondModel(const LinearElasticBondModel* cm)
  : PeridynamicsMaterialModel(cm)
{
}

LinearElasticBondModel::~LinearElasticBondModel()
{
}

void 
LinearElasticBondModel::outputProblemSpec(Uintah::ProblemSpecP& ps,
                                          bool output_cm_tag)
{
  throw SCIRun::InternalError("Stub Task: LinearElasticBondModel::outputProblemSpec ", __FILE__, __LINE__);
}

void 
LinearElasticBondModel::addInitialComputesAndRequires(Uintah::Task* ,
                                                      const PeridynamicsMaterial* ,
                                                      const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: LinearElasticBondModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

/*! Initialize the variables used in the CM */
void 
LinearElasticBondModel::initialize(const Uintah::Patch* patch,
                                   const PeridynamicsMaterial* matl,
                                   Uintah::DataWarehouse* new_dw)
{
  throw SCIRun::InternalError("Stub Task: LinearElasticBondModel::initialize ", __FILE__, __LINE__);
}

void 
LinearElasticBondModel::addComputesAndRequires(Uintah::Task*, 
                                               const PeridynamicsMaterial*,
                                               const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: LinearElasticBondModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
LinearElasticBondModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                         std::vector<const Uintah::VarLabel*>& to)
{
  throw SCIRun::InternalError("Stub Task: LinearElasticBondModel::addParticleState ", __FILE__, __LINE__);
}

// Make a clone of the constitutive model
LinearElasticBondModel* 
LinearElasticBondModel::clone()
{
  return scinew LinearElasticBondModel(*this);
}

