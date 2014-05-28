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
