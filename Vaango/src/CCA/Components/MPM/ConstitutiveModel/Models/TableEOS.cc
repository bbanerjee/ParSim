#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>

#include <Eigen/Sparse>

#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>

using namespace Vaango;
using ProblemSpecP = Uintah::ProblemSpecP;
using ProblemSetupException = Uintah::ProblemSetupException;
using InvalidValue = Uintah::InvalidValue;
using json = nlohmann::json;

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
    addIndependentVariable(name);
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;

  std::cout << "Dependent:";
  index = 0;
  for (const auto& name : depVarNames) {
    addDependentVariable(name);
    std::cout << index << ":" << name << " ";
    index++;
  }
  std::cout << std::endl;
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

void 
TableEOS::addIndependentVariable(const std::string& varName)
{
  d_indepVars.push_back(std::make_unique<IndependentVar>(varName));
}

std::size_t  
TableEOS::addDependentVariable(const std::string& varName)
{
  // Check if dependent variable already exists
  // and just return index if true 
  for (auto ii = 0u; ii < d_depVars.size(); ii++) {
    if (d_depVars[ii]->name == varName) {
      return ii;
    }
  }

  // Name does not exist, add the new variable and return its index
  d_depVars.push_back(std::make_unique<DependentVar>(varName));
  return d_depVars.size()-1;  
}

void 
TableEOS::setup()
{
  if (d_depVars.size() == 4) {
    readJSONTableFromFile<4>(d_filename);
  }
}

template <int dim>
void
TableEOS::readJSONTableFromFile(const std::string& tableFile)
{
  std::ifstream inputFile(tableFile);
  if (!inputFile) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot read tabular input data file "
        << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::stringstream inputStream;
  inputStream << inputFile.rdbuf();
  inputFile.close();

  json doc = loadJSON(inputStream, tableFile);
  readJSONTable<dim>(doc, tableFile);
}

json
TableEOS::loadJSON(std::stringstream& inputStream,
                   const std::string& tableFile)
{
  json doc;
  try {
    doc << inputStream;
  } catch (std::invalid_argument err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot parse tabular input data file "
        << tableFile << "\n"
        << " Please check that the file is valid JSON using a linter";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return doc;
}

namespace Vaango {
template <>
void
TableEOS::readJSONTable<1>(const json& doc,
                           const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  std::vector<double> indepVarData = 
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({IndexKey(0,0,0,0), indepVarData});
  
  for (const auto& depVar : d_depVars) {
    std::vector<double> depVarData = 
      getDoubleArrayJSON(data, depVar->name, tableFile);
    depVar->data.insert({IndexKey(0,0,0,0), depVarData});
  }
}

template <>
void
TableEOS::readJSONTable<4>(const json& doc,
                           const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  std::vector<double> indepVar0 = 
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  json data0 = getDataJSON(data, tableFile);
  
  for (auto ii = 0u; ii < indepVar0.size() ; ii++) {

    std::vector<double> indepVar1 = 
      getDoubleArrayJSON(data0[ii], d_indepVars[1]->name, tableFile);
    json data1 = getDataJSON(data0[ii], tableFile);

    for (auto jj = 0u; jj < indepVar1.size() ; jj++) {

      int index = 0;
      for (const auto& indepVar : d_indepVars) {
        if (index > 1) {
          std::vector<double> indepVar2 = 
            getDoubleArrayJSON(data1[jj], indepVar->name, tableFile);
        }
        index++;
      }

      for (const auto& depVar : d_depVars) {
        std::vector<double> depVar2 = 
          getDoubleArrayJSON(data1[jj], depVar->name, tableFile);
      }
    }
  }
}
}

json
TableEOS::getContentsJSON(const json& doc,
                          const std::string& tableFile)
{
  json contents;
  try {
    contents = doc["Vaango_tabular_data"];
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot find the key \"Vaango_tabular_data\" in "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return contents;
}

std::string
TableEOS::getTitleJSON(const json& contents,
                       const std::string& tableFile)
{
  std::string title;
  try {
    title = contents["Meta"].at("title").get<std::string>();
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Could not find tabular EOS title in "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return title;
}

json
TableEOS::getDataJSON(const json& contents,
                      const std::string& tableFile)
{
  json data;
  try {
    data = contents["Data"];
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Key \"Data\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return data;
}


std::vector<double>
TableEOS::getVectorJSON(const json& object,
                        const std::string key,
                        const std::string& tableFile)
{
  std::string str;
  try {
    str = object.at(key).get<std::string>();
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \""
        << key << "\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  std::vector<std::string> vecStr = Vaango::Util::split(str, ',');
  std::vector<double> vec;
  for (auto& val : vecStr) {
    vec.push_back(Vaango::Util::toDouble(val));
  }
  return vec;
}


std::vector<double>
TableEOS::getDoubleArrayJSON(const json& object,
                             const std::string key,
                             const std::string& tableFile)
{
  // First check if the key exists
  json data;
  try {
    data = object.at(key);
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \""
        << key << "\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::vector<double> vec;
  try {
    vec = data.get<std::vector<double>>();
    /*
    std::cout << key << ":";
    for (const auto& val : vec) {
      std::cout << val << " ";
    }
    std::cout << std::endl;
    */
  } catch (std::exception err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \""
        << key << "\" contains invalid data. "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return vec;
}

double 
TableEOS::interpolate(const std::string& depVarName,
                      std::vector<double>& independents)
{
  return 0; 
}

double
TableEOS::interpolateLinearSpline1D(const double& indepValue,
                                    const std::vector<double>& indepVar,
                                    const std::vector<double>& depVar) const
{
  if (indepValue < *(indepVar.begin()) || indepValue > *(indepVar.end()-1)) {
    std::ostringstream out;
    out << "**ERROR**"
        << " The independent variable value \""
        << indepValue << "\" is outside the range of known data. ";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }
  auto lower = std::find_if(indepVar.begin(), indepVar.end(),
                            [&indepValue](const auto& value) {
                              return (indepValue <= value);
                            });
  double depValue = 0.0;
  if (lower == indepVar.begin()) {
    depValue = *(depVar.begin());
  } else if (lower == indepVar.end()-1) {
    depValue = *(depVar.end()-1);
  } else {
    auto position = lower - indepVar.begin();
    double t = (indepValue - *lower)/(*(lower+1) - *lower);
    depValue = (1 - t)*depVar[position] + t*depVar[position+1];
  }
  return depValue;
}

std::vector<double>
TableEOS::fitCubicSpline1D(const DependentVar depVar,
                          const IndependentVar indepVar) 
{

}

// See: https://cs.uwaterloo.ca/research/tr/1983/CS-83-09.pdf
double
TableEOS::interpolateCubicSpline1D(const double& t) {
  
}

std::vector<double> 
TableEOS::getIndependentVarData(const std::string& name,
                                const IndexKey& index)
{
  auto varIter = std::find_if(d_indepVars.begin(), d_indepVars.end(),
                              [&name](const auto& indepVar) {
                                  return (indepVar->name == name);
                              });
  auto position = std::distance(d_indepVars.begin(), varIter);
  return d_indepVars[position]->data.at(index);
}

std::vector<double> 
TableEOS::getDependentVarData(const std::string& name,
                              const IndexKey& index)
{
  auto varIter = std::find_if(d_depVars.begin(), d_depVars.end(),
                              [&name](const auto& depVar) {
                                  return (depVar->name == name);
                              });
  auto position = std::distance(d_depVars.begin(), varIter);
  return d_depVars[position]->data.at(index);

}