#include <CCA/Components/Peridynamics/DamageModels/PeridynamicsDamageModel.h>
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

PeridynamicsDamageModel::PeridynamicsDamageModel(PeridynamicsLabel* labels,
                                                 PeridynamicsFlags* flags)
{
  d_label = labels;
  d_flags = flags;
}

PeridynamicsDamageModel::PeridynamicsDamageModel(const PeridynamicsDamageModel* cm)
{
  d_label = cm->d_label;
  d_flags = cm->d_flags;
}

PeridynamicsDamageModel::~PeridynamicsDamageModel()
{
}

void 
PeridynamicsDamageModel::addInitialComputesAndRequires(Uintah::Task* ,
                                                       const PeridynamicsMaterial* ,
                                                       const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: PeridynamicsDamageModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsDamageModel::addComputesAndRequires(Uintah::Task*, 
                                               const PeridynamicsMaterial*,
                                               const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: PeridynamicsDamageModel::addComputesAndRequires ", __FILE__, __LINE__);
}
