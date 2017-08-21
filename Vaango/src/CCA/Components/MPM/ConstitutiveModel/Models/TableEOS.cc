#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>

#include <iostream>

using namespace Vaango;
using ProblemSpecP = Uintah::ProblemSpecP;

TableEOS::TableEOS(ProblemSpecP& ps)
{
  std::string indep_vars;
  std::string dep_vars;
  ps->require("independent_variables", indep_vars);
  ps->require("dependent_variables", dep_vars);
  ps->require("filename", d_filename);

  std::vector<std::string> indepVarNames = parseVariableNames(indep_vars);
  std::vector<std::string> depVarNames = parseVariableNames(dep_vars);

  std::cout << "Independent:";
  int index = 0;
  for (const auto& name : indepVarNames) {
    d_indepVars.push_back(std::make_unique<IndependentVar>(name));
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;

  std::cout << "Dependent:";
  index = 0;
  for (const auto& name : depVarNames) {
    d_depVars.push_back(std::make_unique<DependentVar>(name));
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;
}

void 
TableEOS::addIndependentVariable(const std::string& varName)
{
  
}

int  
TableEOS::addDependentVariable(const std::string& varName)
{
  return 0;  
}

void 
TableEOS::setup()
{
  
}

double 
TableEOS::interpolate(int index,
                      std::vector<double>& independents)
{
  return 0; 
}

std::vector<std::string> 
TableEOS::parseVariableNames(const std::string& vars)
{
  std::vector<std::string> varNames = Vaango::Util::split(vars, ',');
  for (auto& name : varNames) {
    name = Vaango::Util::trim(name);
  }
  return varNames;
}

