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

#include <CCA/Components/Peridynamics/FailureModels/BondDamageModel.h>
#include <CCA/Components/Peridynamics/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/PeridynamicsLabel.h>

#include <CCA/Ports/DataWarehouse.h>

#include <Core/Grid/Patch.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/ParticleVariable.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Math/Matrix3.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Malloc/Allocator.h>

#include <cmath>
#include <iostream>

using namespace Vaango;

BondDamageModel::BondDamageModel(Uintah::ProblemSpecP& ps,
                                 PeridynamicsFlags* flags)
  : PeridynamicsFailureModel(flags)
{
}

BondDamageModel::BondDamageModel(const BondDamageModel* cm)
  : PeridynamicsFailureModel(cm)
{
}

BondDamageModel::~BondDamageModel()
{
}

void 
BondDamageModel::outputProblemSpec(Uintah::ProblemSpecP& ps,
                                   bool output_cm_tag)
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::outputProblemSpec", __FILE__, __LINE__);
}
         
void 
BondDamageModel::addInitialComputesAndRequires(Uintah::Task* ,
                                               const PeridynamicsMaterial* ,
                                               const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

void 
BondDamageModel::initialize(const Uintah::Patch* patch,
                            const PeridynamicsMaterial* matl,
                            Uintah::DataWarehouse* new_dw)
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::initialize", __FILE__, __LINE__);
}

void 
BondDamageModel::addComputesAndRequires(Uintah::Task*, 
                                        const PeridynamicsMaterial*,
                                        const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
BondDamageModel::updateDamage(const Uintah::PatchSubset*,
                              const PeridynamicsMaterial*,
                              Uintah::DataWarehouse*,
                              Uintah::DataWarehouse*)
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::updateDamage ", __FILE__, __LINE__);
}

void 
BondDamageModel::addParticleState(std::vector<const Uintah::VarLabel*>& from,
                                  std::vector<const Uintah::VarLabel*>& to)
{
  throw SCIRun::InternalError("Stub Task: BondDamageModel::addParticleState ", __FILE__, __LINE__);
}

BondDamageModel* 
BondDamageModel::clone()
{
  return scinew BondDamageModel(*this);
}
