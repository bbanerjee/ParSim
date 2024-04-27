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

#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsMaterial.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsFlags.h>
#include <CCA/Components/Peridynamics/Core/PeridynamicsLabel.h>

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

PeridynamicsMaterialModel::PeridynamicsMaterialModel(PeridynamicsFlags* flags)
{
  d_label = scinew PeridynamicsLabel();
  d_flag = flags;
  d_numGhostNodes = 2;
}

PeridynamicsMaterialModel::PeridynamicsMaterialModel(const PeridynamicsMaterialModel* cm)
{
  d_label = scinew PeridynamicsLabel();
  d_flag = cm->d_flag;
  d_numGhostNodes = cm->d_numGhostNodes;
  d_mat_manager = cm->d_mat_manager;
}

PeridynamicsMaterialModel::~PeridynamicsMaterialModel()
{
  delete d_label;
}

void 
PeridynamicsMaterialModel::addInitialComputesAndRequires(Uintah::Task* ,
                                                         const PeridynamicsMaterial* ,
                                                         const Uintah::PatchSet*) const
{
  throw Uintah::InternalError("Stub Task: PeridynamicsMaterialModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsMaterialModel::addComputesAndRequires(Uintah::Task*, 
                                                  const PeridynamicsMaterial*,
                                                  const Uintah::PatchSet*) const
{
  throw Uintah::InternalError("Stub Task: PeridynamicsMaterialModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsMaterialModel::computeStressTensor(const Uintah::PatchSubset*,
                                               const PeridynamicsMaterial*,
                                               Uintah::DataWarehouse*,
                                               Uintah::DataWarehouse*)
{
  throw Uintah::InternalError("Stub Task: PeridynamicsMaterialModel::computeStressTensor ", __FILE__, __LINE__);
}

