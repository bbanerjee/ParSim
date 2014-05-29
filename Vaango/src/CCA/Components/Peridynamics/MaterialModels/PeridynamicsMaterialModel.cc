#include <CCA/Components/Peridynamics/MaterialModels/PeridynamicsMaterialModel.h>
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
  d_sharedState = cm->d_sharedState;
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
  throw SCIRun::InternalError("Stub Task: PeridynamicsMaterialModel::addInitialComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsMaterialModel::addComputesAndRequires(Uintah::Task*, 
                                                  const PeridynamicsMaterial*,
                                                  const Uintah::PatchSet*) const
{
  throw SCIRun::InternalError("Stub Task: PeridynamicsMaterialModel::addComputesAndRequires ", __FILE__, __LINE__);
}

void 
PeridynamicsMaterialModel::computeStressTensor(const Uintah::PatchSubset*,
                                               const PeridynamicsMaterial*,
                                               Uintah::DataWarehouse*,
                                               Uintah::DataWarehouse*)
{
  throw SCIRun::InternalError("Stub Task: PeridynamicsMaterialModel::computeStressTensor ", __FILE__, __LINE__);
}

