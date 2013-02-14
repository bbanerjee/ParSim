#include <Common/VarLabel.h>
#include <map>
#include <iostream>
#include <sstream>

using namespace Matiti;

static map<string, VarLabel*> allLabels;

VarLabel*
VarLabel::create(const string& name,
                 const Matiti::TypeDescription* td,
                 VarType vartype /*= Normal*/)
{
  VarLabel* label = 0;
  map<string, VarLabel*>::iterator iter = allLabels.find(name);
  if(iter != allLabels.end()){
    // two labels with the same name -- make sure they are the same type
    VarLabel* dup = iter->second;
    label = dup;
  }
  else {
    label = new VarLabel(name, td, vartype);
    allLabels[name]=label;
  }
  label->addReference();
  
  return label;
}

bool
VarLabel::destroy(const VarLabel* label)
{
  if (label == 0) return false;
  if (label->removeReference()) {
    map<string, VarLabel*>::iterator iter = allLabels.find(label->d_name);
    if(iter != allLabels.end() && iter->second == label)
      allLabels.erase(iter); 
      
    delete label;
    
    return true;
  }

  return false;
}

VarLabel::VarLabel(const std::string& name, const Matiti::TypeDescription* td,
                   VarType vartype)
  : d_name(name), d_boundaryLayer(boundaryLayer),
    d_vartype(vartype), 
    d_allowMultipleComputes(false)
{
}

VarLabel::~VarLabel()
{
}

void
VarLabel::printAll()
{
  map<string, VarLabel*>::iterator iter = allLabels.begin();
  for (; iter != allLabels.end(); iter++)
    std::cerr << (*iter).second->d_name << std::endl;
}

VarLabel*
VarLabel::find(string name)
{
   map<string, VarLabel*>::iterator found = allLabels.find(name);
   if (found == allLabels.end())
      return NULL;
   else
      return found->second;
}


string
VarLabel::getFullName(int matlIndex) const
{
   ostringstream out;
        out << d_name << "matl=" << matlIndex;
        return out.str();
}                             

namespace Matiti {
  ostream & 
  operator<<( ostream & out, const Matiti::VarLabel & vl )
  {
    out << vl.getName();
    return out;
  }
}

