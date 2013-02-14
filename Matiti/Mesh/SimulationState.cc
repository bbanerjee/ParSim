#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/SimulationState.h>
#include <Core/Grid/Variables/VarLabel.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <Core/Grid/Variables/ReductionVariable.h>
#include <Core/Grid/Variables/Reductions.h>
#include <Core/Grid/Variables/PerPatch.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/Grid.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/SimpleMaterial.h>
#include <CCA/Components/ICE/ICEMaterial.h>
#include <CCA/Components/MPM/ConstitutiveModel/PeriMaterial.h>
#include <CCA/Components/MPM/CohesiveZone/CZMaterial.h>
#include <Core/Containers/StringUtil.h>
#include <Core/Malloc/Allocator.h>

using namespace Matiti;

SimulationState::SimulationState(Uintah::ProblemSpecP &ps)
{
  Uintah::VarLabel* nonconstDelt = 
     Uintah::VarLabel::create("delT", delt_vartype::getTypeDescription() );
  nonconstDelt->allowMultipleComputes();
  delt_label = nonconstDelt;

  all_peri_matls = 0;
  all_matls = 0;
  max_matl_index = 0;
  d_simTime = 0;
}

void 
SimulationState::registerMaterial(Uintah::Material* matl)
{
  matl->registerParticleState(this);
  matl->setDWIndex((int)matls.size());

  matls.push_back(matl);
  if ((int)matls.size() > max_matl_index) {
    max_matl_index = matls.size();
  }

  if(matl->hasName())
    named_matls[matl->getName()] = matl;
}

void 
SimulationState::registerMaterial(Uintah::Material* matl, unsigned int index)
{
  matl->registerParticleState(this);
  matl->setDWIndex(index);

  if (matls.size() <= index)
    matls.resize(index+1);
  matls[index]=matl;

  if ((int)matls.size() > max_matl_index) {
    max_matl_index = matls.size();
  }

  if(matl->hasName())
    named_matls[matl->getName()] = matl;
}

void 
SimulationState::registerPeriMaterial(PeriMaterial* matl)
{
  peri_matls.push_back(matl);
  registerMaterial(matl);
}

void 
SimulationState::registerPeriMaterial(PeriMaterial* matl,unsigned int index)
{
  peri_matls.push_back(matl);
  registerMaterial(matl,index);
}

void 
SimulationState::finalizeMaterials()
{
  if (all_matls && all_matls->removeReference())
    delete all_matls;
  all_matls = new Uintah::MaterialSet();
  all_matls->addReference();
  vector<int> tmp_matls(matls.size());
  for( int i=0; i<(int)matls.size(); i++ ) {
    tmp_matls[i] = matls[i]->getDWIndex();
  }
  all_matls->addAll(tmp_matls);

  if (all_peri_matls && all_peri_matls->removeReference())
    delete all_peri_matls;
  all_peri_matls = new Uintah::MaterialSet();
  all_peri_matls->addReference();
  vector<int> tmp_peri_matls(peri_matls.size());
  for( int i=0; i<(int)peri_matls.size(); i++ ) {
    tmp_peri_matls[i] = peri_matls[i]->getDWIndex();
  }
  all_peri_matls->addAll(tmp_peri_matls);
}

void 
SimulationState::clearMaterials()
{
  for (int i = 0; i < (int)matls.size(); i++)
    old_matls.push_back(matls[i]);

  if(all_matls && all_matls->removeReference())
    delete all_matls;
  
  if(all_peri_matls && all_peri_matls->removeReference())
    delete all_peri_matls;

  matls.clear();
  peri_matls.clear();
  d_meshNodeState.clear();
  d_meshNodeState_preReloc.clear();

  all_matls         = 0;
  all_peri_matls     = 0;
}

SimulationState::~SimulationState()
{
  Uintah::VarLabel::destroy(delt_label);
  clearMaterials();
}

const Uintah::MaterialSet* 
SimulationState::allPeriMaterials() const
{
  //ASSERT(all_peri_matls != 0);
  return all_peri_matls;
}

const Uintah::MaterialSet* 
SimulationState::allMaterials() const
{
  //ASSERT(all_matls != 0);
  return all_matls;
}

Uintah::Material* 
SimulationState::getMaterialByName(const std::string& name) const
{
  map<string, Uintah::Material*>::const_iterator iter = named_matls.find(name);
  if(iter == named_matls.end())
    return 0;
  return iter->second;
}

//__________________________________
//
Uintah::Material* 
SimulationState::parseAndLookupMaterial(Uintah::ProblemSpecP& params,
                                        const std::string& name) const
{
  // for single material problems return matl 0
  Uintah::Material* result = getMaterial(0);

  if( getNumMatls() > 1){
    string matlname;
    if(!params->get(name, matlname)){
      throw Uintah::ProblemSetupException("Cannot find material section", __FILE__, __LINE__);
    }

    result = getMaterialByName(matlname);
    if(!result){ 
      throw Uintah::ProblemSetupException("Cannot find a material named:"+matlname, __FILE__, __LINE__);
    }
  }
  return result;
}

