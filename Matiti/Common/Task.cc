#include <Core/Grid/Task.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Grid.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Containers/StringUtil.h>
#include <set>


using namespace std;
using namespace Matiti;

MaterialSubset* Task::globalMatlSubset = 0;

void Task::initialize()
{
  comp_head = comp_tail = 0;
  req_head = req_tail = 0;
  mod_head = mod_tail = 0;
  matl_set = 0;
  d_hasSubScheduler = false;

  for(int i=0;i<TotalDWs;i++) {
    dwmap[i]=Task::InvalidDW;
  }
  sortedOrder=-1;
  d_phase=-1;
}

Task::ActionBase::~ActionBase()
{
}

Task::ActionGPUBase::~ActionGPUBase()
{
}

Task::~Task()
{
  delete d_action;

  Dependency* dep = req_head;
  while(dep){
    Dependency* next = dep->next;
    delete dep;
    dep=next;
  }

  dep = comp_head;
  while(dep){
    Dependency* next = dep->next;
    delete dep;
    dep=next;
  }

  dep = mod_head;
  while(dep){
    Dependency* next = dep->next;
    delete dep;
    dep=next;
  }
  
  if(matl_set && matl_set->removeReference()) {
    delete matl_set;
  }

  // easier to periodically delete this than to force a call to a cleanup
  // function, and probably not very expensive.
  if (globalMatlSubset && globalMatlSubset->removeReference()) {
    delete globalMatlSubset;
  }

  globalMatlSubset = 0;
}

//__________________________________
void Task::setSets(const PatchSet* ps, const MaterialSet* ms)
{
  //ASSERT(matl_set == 0);
  matl_set=ms;
  if(matl_set) {
    matl_set->addReference();
  }
}

//__________________________________
const MaterialSubset* Task::getGlobalMatlSubset()
{
  if (globalMatlSubset == 0) {
    globalMatlSubset = new MaterialSubset();
    globalMatlSubset->add(-1);
    globalMatlSubset->addReference();
  }
  return globalMatlSubset;
}

void
Task::hasSubScheduler(bool state)
{
  d_hasSubScheduler = state;
}

//__________________________________
void
Task::requires(WhichDW dw, 
               const VarLabel* var,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom)
{
  if (matls != 0 && matls->size() == 0) {
    return; // no materials, no dependency
  }

  Dependency* dep = new Dependency(Requires, this, dw, var, matls, matls_dom);

  dep->next=0;
  if(req_tail)
    req_tail->next=dep;
  else
    req_head=dep;
  req_tail=dep;
  if (dw == OldDW)
    d_requiresOldDW.insert(make_pair(var, dep));
  else
    d_requires.insert(make_pair(var, dep));
}

//__________________________________
void
Task::requires(WhichDW dw, 
               const VarLabel* var,
               const MaterialSubset* matls)
{
  requires(dw, var, matls, NormalDomain);
}

//__________________________________
void
Task::requires(WhichDW dw, 
               const VarLabel* var)
{
  requires(dw, var, 0, NormalDomain);
}

//__________________________________
void
Task::computes(const VarLabel * var,
               const MaterialSubset * matls,
               MaterialDomainSpec matls_dom)
{
  if (var->typeDescription()->isReductionVariable()) {
    if (matls == 0) {
      // default material for a reduction variable is the global material (-1)
      matls = getGlobalMatlSubset();
      matls_dom = OutOfDomain;
    }
  }
  
  Dependency* dep = new Dependency(Computes, this, NewDW, var, false, matls, matls_dom);
  dep->next=0;
  if(comp_tail)
    comp_tail->next=dep;
  else
    comp_head=dep;
  comp_tail=dep;

  d_computes.insert(make_pair(var, dep));
}

//__________________________________
void
Task::computes(const VarLabel * var,
               const MaterialSubset * matls)
{
  computes(var, matls, NormalDomain);
}

//__________________________________
void
Task::computes(const VarLabel * var)
{
  computes(var, 0, NormalDomain);
}


//__________________________________
void
Task::modifies(const VarLabel* var,
               const MaterialSubset* matls,
               MaterialDomainSpec matls_dom)
{
  if (matls == 0 && var->typeDescription()->isReductionVariable()) {
    // default material for a reduction variable is the global material (-1)
    matls = getGlobalMatlSubset();
    matls_dom = OutOfDomain;
  }  

  Dependency* dep = new Dependency(Modifies, this, NewDW, var, matls, matls_dom);
  dep->next=0;
  if (mod_tail)
    mod_tail->next=dep;
  else
    mod_head=dep;
  mod_tail=dep;

  d_requires.insert(make_pair(var, dep));
  d_computes.insert(make_pair(var, dep));
  d_modifies.insert(make_pair(var, dep));
}

//__________________________________
void
Task::modifies(const VarLabel* var,
               const MaterialSubset* matls)
{
  modifies(var, matls, NormalDomain);
}

//__________________________________
void
Task::modifies(const VarLabel* var)
{
  modifies(var, 0, NormalDomain);
}


//__________________________________
bool Task::hasComputes(const VarLabel* var, int matlIndex) const
{
  return isInDepMap(d_computes, var, matlIndex);
}

//__________________________________
bool Task::hasRequires(const VarLabel * var,
                       int matlIndex,
                       WhichDW dw)const
{
  DepMap depMap = d_requires;
  
  if(dw == OldDW){
    depMap = d_requiresOldDW;
  }
  
  Dependency* dep = isInDepMap(depMap, var, matlIndex);  
  
  if (dep) return true;
  return false;
}

//__________________________________
bool Task::hasModifies(const VarLabel* var, int matlIndex) const
{
  return isInDepMap(d_modifies, var, matlIndex);
}

//__________________________________
Task::Dependency* Task::isInDepMap(const DepMap& depMap, 
                                   const VarLabel* var,
                                   int matlIndex) const
{
  DepMap::const_iterator found_iter = depMap.find(var);
  
  // loop over dependency map and search for the right dependency
  
  while (found_iter != depMap.end() && (*found_iter).first->equals(var)) {
  
    Dependency* dep = (*found_iter).second;
    const MaterialSubset* matls = dep->matls;

    bool hasMatls=false;

    if (matls == 0) //if matls==0 then the requierment for matls is satisfied
    {
      hasMatls=true;
    }
    else  //check thta the malts subset contains the matlIndex
    {
      hasMatls=matls->contains(matlIndex);
    }
   
    if(hasMatls)  //if this dependency contains both the matls and patches return the dependency
      return dep;

    found_iter++;
  }
  return 0;
}
//__________________________________
//
Task::Dependency::Dependency(DepType deptype, 
                             Task* task, 
                             WhichDW whichdw,
                             const VarLabel* var,
                             const MaterialSubset* matls,
                             MaterialDomainSpec matls_dom)
                             
  : deptype(deptype), task(task), var(var), matls(matls), matls_dom(matls_dom), whichdw(whichdw)
{
  if (var)
    var->addReference();
  req_head=req_tail=comp_head=comp_tail=0;
  if(matls)
    matls->addReference();
}

Task::Dependency::~Dependency()
{
  VarLabel::destroy(var); // just remove the ref
  if(matls && matls->removeReference())
    delete matls;
}

//__________________________________
constHandle<MaterialSubset>
Task::Dependency::getMaterialsUnderDomain(const MaterialSubset* domainMaterials) const
{
  switch(matls_dom){
  case Task::NormalDomain:
    return MaterialSubset::intersection(matls, domainMaterials);
  case Task::OutOfDomain:
    return matls;
  default:
    throw(Uintah::InternalError(string("Unknown matl domain ") + " type "+SCIRun::to_string(static_cast<int>(matls_dom)),
                            __FILE__, __LINE__));
  }
}

//__________________________________
void
Task::doit(const MaterialSubset* matls,
           vector<DataWarehouseP>& dws)
{
  DataWarehouse* fromDW = mapDataWarehouse(Task::OldDW, dws);
  DataWarehouse* toDW = mapDataWarehouse(Task::NewDW, dws);
  if(d_action) {
    d_action->doit(matls, fromDW, toDW);
  }
}

//__________________________________
void
Task::display( ostream & out ) const
{
  out <<  getName() << " (" << d_tasktype << "): [";
  out << ", ";
  out << *matl_set;
  out << ", DWs: ";
  for(int i=0;i<TotalDWs;i++){
    if(i != 0)
      out << ", ";
    out << dwmap[i];
  }
  out << "]";
}
//__________________________________
namespace Uintah {
  std::ostream &
  operator << ( std::ostream & out, const Uintah::Task::Dependency & dep )
  {
    out << "[";
    out<< left;out.width(20);
    out << *(dep.var) << ", ";

    out << ", MI: ";
    if(dep.matls){
      for(int i=0;i<dep.matls->size();i++){
        if(i>0)
          out << ",";
        out << dep.matls->get(i);
      }
    } else {
      out << "none";
    }
    out << ", ";
    switch(dep.whichdw){
    case Task::OldDW:
      out << "OldDW";
      break;
    case Task::NewDW:
      out << "NewDW";
      break;
    case Task::ParentOldDW:
      out << "ParentOldDW";
      break;
    case Task::ParentNewDW:
      out << "ParentNewDW";
      break;
    default:
      out << "Unknown DW!";
      break;
    }
    out << " (mapped to dw index " << dep.task->mapDataWarehouse(dep.whichdw) << ")";
    out << ", ";
    out << "]";
    return out;
  }
  
//__________________________________
  ostream &
  operator << (ostream& out, const Task& task)
  {
    task.display( out );
    return out;
  }
  
//__________________________________
  ostream&
  operator << (ostream &out, const Task::TaskType & tt)
  {
    switch( tt ) {
    case Task::Normal:
      out << "Normal";
      break;
    case Task::Reduction:
      out << "Reduction";
      break;
    case Task::InitialSend:
      out << "InitialSend";
      break;
    case Task::Output:
      out << "Output";
      break;
    }
    return out;
  }
} // end namespace Uintah

//__________________________________
void
Task::displayAll(ostream& out) const
{
   display(out);
   out << '\n';
   for(Task::Dependency* req = req_head; req != 0; req = req->next)
      out << "  requires: " << *req << '\n';
   for(Task::Dependency* comp = comp_head; comp != 0; comp = comp->next)
      out << "  computes: " << *comp << '\n';
   for(Task::Dependency* mod = mod_head; mod != 0; mod = mod->next)
      out << "  modifies: " << *mod << '\n';
}

//__________________________________
void Task::setMapping(int dwmap[TotalDWs])
{
  for(int i=0;i<TotalDWs;i++)
  {
    this->dwmap[i]=dwmap[i];
  }
}

//__________________________________
int Task::mapDataWarehouse(WhichDW dw) const
{
  //ASSERTRANGE(dw, 0, Task::TotalDWs);
  return dwmap[dw];
}

//__________________________________
DataWarehouse* Task::mapDataWarehouse(WhichDW dw, vector<DataWarehouseP>& dws) const
{
  //ASSERTRANGE(dw, 0, Task::TotalDWs);
  if(dwmap[dw] == Task::NoDW){
    return 0;
  } else {
    //ASSERTRANGE(dwmap[dw], 0, (int)dws.size());
    return dws[dwmap[dw]].get_rep();
  }
}

