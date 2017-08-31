#include <CCA/Components/MPM/ConstitutiveModel/Models/TableEOS.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TableUtils.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Util/FancyAssert.h>

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

  //std::cout << "Independent:";
  int index = 0;
  for (const auto& name : indepVarNames) {
    addIndependentVariable(name);
    //std::cout << index << ":" << name << " ";
    index++;
  }
  //std::cout << std::endl;

  //std::cout << "Dependent:";
  index = 0;
  for (const auto& name : depVarNames) {
    addDependentVariable(name);
    //std::cout << index << ":" << name << " ";
    index++;
  }
  //std::cout << std::endl;
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
  if (d_indepVars.size() == 1) {
    readJSONTableFromFile<1>(d_filename);
  } else if (d_indepVars.size() == 2) {
    readJSONTableFromFile<2>(d_filename);
  } else if (d_indepVars.size() == 3) {
    readJSONTableFromFile<3>(d_filename);
  } else {
    std::ostringstream out;
    out << "**ERROR**"
        << " More than three independent variables not allowed in "
        << d_filename;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
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

  DoubleVec1D indepVarData = 
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({IndexKey(0,0,0,0), indepVarData});
  
  for (const auto& depVar : d_depVars) {
    DoubleVec1D depVarData = 
      getDoubleArrayJSON(data, depVar->name, tableFile);
    depVar->data.insert({IndexKey(0,0,0,0), depVarData});

    ASSERTEQ(indepVarData.size(), depVarData.size());
  }
}

template <>
void
TableEOS::readJSONTable<2>(const json& doc,
                           const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  DoubleVec1D indepVar0Data = 
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({IndexKey(0,0,0,0), indepVar0Data});

  json data0 = getDataJSON(data, tableFile);
  for (auto ii = 0u; ii < indepVar0Data.size() ; ii++) {

    DoubleVec1D indepVar1Data = 
      getDoubleArrayJSON(data0[ii], d_indepVars[1]->name, tableFile);
    /*
    std::cout << "Read " << d_indepVars[1]->name << " ii = " << ii << " ";
    std::copy(indepVar1Data.begin(), indepVar1Data.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    */
    d_indepVars[1]->data.insert({IndexKey(ii,0,0,0), indepVar1Data});

    for (const auto& depVar : d_depVars) {
      DoubleVec1D depVarData = 
        getDoubleArrayJSON(data0[ii], depVar->name, tableFile);
      depVar->data.insert({IndexKey(ii,0,0,0), depVarData});

      ASSERTEQ(indepVar1Data.size(), depVarData.size());
    }

  }
}

template <>
void
TableEOS::readJSONTable<3>(const json& doc,
                           const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  DoubleVec1D indepVar0Data = 
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({IndexKey(0,0,0,0), indepVar0Data});

  json data0 = getDataJSON(data, tableFile);
  for (auto ii = 0u; ii < indepVar0Data.size() ; ii++) {

    DoubleVec1D indepVar1Data = 
      getDoubleArrayJSON(data0[ii], d_indepVars[1]->name, tableFile);
    d_indepVars[1]->data.insert({IndexKey(ii,0,0,0), indepVar1Data});

    json data1 = getDataJSON(data0[ii], tableFile);
    for (auto jj = 0u; jj < indepVar1Data.size() ; jj++) {

      int index = 0;
      for (const auto& indepVar : d_indepVars) {
        if (index > 1) {
          DoubleVec1D indepVar2Data = 
            getDoubleArrayJSON(data1[jj], indepVar->name, tableFile);
          indepVar->data.insert({IndexKey(ii,jj,0,0), indepVar2Data});
        }
        index++;
      }

      for (const auto& depVar : d_depVars) {
        DoubleVec1D depVar2Data = 
          getDoubleArrayJSON(data1[jj], depVar->name, tableFile);
        depVar->data.insert({IndexKey(ii,jj,0,0), depVar2Data});
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


DoubleVec1D
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
  DoubleVec1D vec;
  for (auto& val : vecStr) {
    vec.push_back(Vaango::Util::toDouble(val));
  }
  return vec;
}


DoubleVec1D
TableEOS::getDoubleArrayJSON(const json& object,
                             const std::string key,
                             const std::string& tableFile)
{
  // First check if the key exists
  json data;
  try {
    data = object.at(key);
    //std::cout << "Data = " << data << std::endl;
  } catch (std::out_of_range err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \""
        << key << "\" not found in tabular EOS input file "
        << tableFile ;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  DoubleVec1D vec;
  try {
    vec = data.get<DoubleVec1D>();
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
                      DoubleVec1D& independents)
{
  return 0; 
}

namespace Vaango {
template <>
DoubleVec1D
TableEOS::interpolateLinearSpline<1>(const std::array<double,1>& indepValues,
                                     const IndepVarPArray& indepVars,
                                     const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  auto numIndepVars = indepVars.size();
  auto numDepVars = depVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 = getIndependentVarData(indepVars[0]->name, 
                                             TableEOS::IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocation(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  //std::cout << "Lo Index 0 = " << segLowIndex0 << " s = " << sval << std::endl;

  DoubleVec1D depVals;
  for (const auto& depVar : depVars) {
    auto depVarData = getDependentVarData(depVar->name, 
                                          TableEOS::IndexKey(0, 0, 0, 0));
    auto depval = computeInterpolated(sval, segLowIndex0, depVarData);
    depVals.push_back(depval);
    //std::cout << "Lo Index 0 = " << segLowIndex0 << " p = " << depval << std::endl;
  }

  return depVals;
}

template <>
DoubleVec1D
TableEOS::interpolateLinearSpline<2>(const std::array<double,2>& indepValues,
                                     const IndepVarPArray& indepVars,
                                     const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  auto numIndepVars = indepVars.size();
  auto numDepVars = depVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 = getIndependentVarData(indepVars[0]->name, 
                                             TableEOS::IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocation(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  //std::cout << "Lo Index 0 = " << segLowIndex0 << " s = " << sval << std::endl;

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  DoubleVec1D depValsT;
  for (auto ii = segLowIndex0; ii <= segLowIndex0+1; ii++) {
    auto indepVarData1 = getIndependentVarData(indepVars[1]->name, 
                                               TableEOS::IndexKey(ii, 0, 0, 0));
    auto segLowIndex1 = findLocation(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);
    //std::cout << "Lo Index 1 = " << segLowIndex1 << " t = " << tval << std::endl;

    for (const auto& depVar : depVars) {
      auto depVarData = getDependentVarData(depVar->name, 
                                            TableEOS::IndexKey(ii, 0, 0, 0));
      auto depvalT = computeInterpolated(tval, segLowIndex1, depVarData);
      depValsT.push_back(depvalT);
      //std::cout << "Lo Index 1 = " << segLowIndex1 << " p = " << depvalT << std::endl;
    }
  }

  DoubleVec1D depVals;
  auto index = 0u;
  for (const auto& depVar : depVars) {
    auto depval = (1 - sval)*depValsT[index] + sval*depValsT[index+numDepVars];
    depVals.push_back(depval);
    //std::cout << "Lo Index 0 = " << segLowIndex0 << " q = " << depval << std::endl;
    index++;
  }

  return depVals;
}

template <>
DoubleVec1D
TableEOS::interpolateLinearSpline<3>(const std::array<double,3>& indepValues,
                                     const IndepVarPArray& indepVars,
                                     const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  auto numIndepVars = indepVars.size();
  auto numDepVars = depVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 = getIndependentVarData(indepVars[0]->name, 
                                             TableEOS::IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocation(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  DoubleVec1D depValsT;
  for (auto ii = segLowIndex0; ii <= segLowIndex0+1; ii++) {
    auto indepVarData1 = getIndependentVarData(indepVars[1]->name, 
                                               TableEOS::IndexKey(ii, 0, 0, 0));
    auto segLowIndex1 = findLocation(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);

    // Choose the last two vectors containing the relevant independent variable data
    // and find the segments containing the data
    DoubleVec1D depValsU;
    for (auto jj = segLowIndex1; jj <= segLowIndex1+1; jj++) {
      auto indepVarData2 = getIndependentVarData(indepVars[2]->name, 
                                                 TableEOS::IndexKey(ii, jj, 0, 0));
      auto segLowIndex2 = findLocation(indepValues[2], indepVarData2);
      auto uval = computeParameter(indepValues[2], segLowIndex2, indepVarData2);

      for (const auto& depVar : depVars) {
        auto depVarData = getDependentVarData(depVar->name, 
                                              TableEOS::IndexKey(ii, jj, 0, 0));
        auto depvalU = computeInterpolated(uval, segLowIndex2, depVarData);
        depValsU.push_back(depvalU);
      }
    }

    auto index = 0u;
    for (const auto& depVar : depVars) {
      auto depvalT = (1 - tval)*depValsU[index] + tval*depValsU[index+numDepVars];
      depValsT.push_back(depvalT);
      index++;
    }
  }

  DoubleVec1D depVals;
  auto index = 0u;
  for (const auto& depVar : depVars) {
    auto depval = (1 - sval)*depValsT[index] + sval*depValsT[index+numDepVars];
    depVals.push_back(depval);
    index++;
  }

  return depVals;
}
} // end namespace for template specialization (needed by gcc)

std::size_t
TableEOS::findLocation(const double& value,
                       const DoubleVec1D& varData) const
{
  if (value < varData.front() || value > varData.back()) {
    std::ostringstream out;
    out << "**ERROR**"
        << " The independent variable value \""
        << value << "\" is outside the range of known data. ";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  auto lower = std::lower_bound(varData.begin(), varData.end(), value);
  return (lower - varData.begin() - 1);
}

double
TableEOS::computeParameter(const double& input,
                           const std::size_t& startIndex,
                           const DoubleVec1D& data) const
{
  if (data.size() < 2) return 0.0;
  auto t = (input - data[startIndex]) / (data[startIndex+1] - data[startIndex]);
  return t;
}

double
TableEOS::computeInterpolated(const double& tval,
                              const std::size_t& startIndex,
                              const DoubleVec1D& data) const
{
  if (data.empty()) {
    std::ostringstream out;
    out << "**ERROR**"
        << " No data available for interpolation.";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  } 

  return (data.size() == 1) ? data.front() : 
           (1 - tval)*data[startIndex] + tval*data[startIndex+1];
}

DoubleVec1D
TableEOS::fitCubicSpline1D(const DependentVar depVar,
                           const IndependentVar indepVar) 
{

}

// See: https://cs.uwaterloo.ca/research/tr/1983/CS-83-09.pdf
double
TableEOS::interpolateCubicSpline1D(const double& t) {
  
}

DoubleVec1D 
TableEOS::getIndependentVarData(const std::string& name,
                                const IndexKey& index) const
{
  auto varIter = std::find_if(d_indepVars.begin(), d_indepVars.end(),
                              [&name](const auto& indepVar) {
                                  return (indepVar->name == name);
                              });
  auto position = std::distance(d_indepVars.begin(), varIter);
  return d_indepVars[position]->data.at(index);
}

DoubleVec1D 
TableEOS::getDependentVarData(const std::string& name,
                              const IndexKey& index) const
{
  auto varIter = std::find_if(d_depVars.begin(), d_depVars.end(),
                              [&name](const auto& depVar) {
                                  return (depVar->name == name);
                              });
  auto position = std::distance(d_depVars.begin(), varIter);
  return d_depVars[position]->data.at(index);
}