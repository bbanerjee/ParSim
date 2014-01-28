#include <CCA/Components/Peridynamics/FailureModels/PeridynamicsFailureModel.h>
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

PeridynamicsFailureModel::PeridynamicsFailureModel(PeridynamicsFlags* flags)
{
  d_varLabel = scinew PeridynamicsLabel();
  d_flag = flags;
}

PeridynamicsFailureModel::PeridynamicsFailureModel(const PeridynamicsFailureModel* cm)
{
  d_varLabel = scinew PeridynamicsLabel();
  d_flag = cm->d_flag;
}

PeridynamicsFailureModel::~PeridynamicsFailureModel()
{
  delete d_varLabel;
}

void 
PeridynamicsFailureModel::addInitialComputesAndRequires(Uintah::Task* ,
                                                         const PeridynamicsMaterial* ,
                                                         const Uintah::PatchSet*) const
{
  throw InternalError("Stub Task: PeridynamicsFailureModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsFailureModel::addComputesAndRequires(Uintah::Task*, 
                                                 const PeridynamicsMaterial*,
                                                 const Uintah::PatchSet*) const
{
  throw InternalError("Stub Task: PeridynamicsFailureModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsFailureModel::updateDamage(const Uintah::PatchSubset*,
                                       const PeridynamicsMaterial*,
                                       Uintah::DataWarehouse*,
                                       Uintah::DataWarehouse*)
{
  throw InternalError("Stub Task: PeridynamicsFailureModel::updateDamage ", __FILE__, __LINE__);
}

