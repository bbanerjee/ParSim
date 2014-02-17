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

