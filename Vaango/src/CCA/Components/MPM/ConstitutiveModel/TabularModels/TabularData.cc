/*
 * The MIT License
 *
 * Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */
#include <CCA/Components/MPM/ConstitutiveModel/Utilities/TableUtils.h>
#include <CCA/Components/MPM/ConstitutiveModel/TabularModels/TabularData.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Util/FancyAssert.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace Vaango;
using IndependentVar = TableContainers::IndependentVar;
using DependentVar = TableContainers::DependentVar;
using ProblemSpecP = Uintah::ProblemSpecP;
using ProblemSetupException = Uintah::ProblemSetupException;
using InvalidValue = Uintah::InvalidValue;
using Vector = Uintah::Vector;
using json = nlohmann::json;

//#define DEBUG_TABLE_INTERPOLATION

TabularData::TabularData(ProblemSpecP& ps)
{
  ps->require("independent_variables", d_indepVarNames);
  ps->require("dependent_variables", d_depVarNames);
  ps->require("filename", d_filename);
  //ps->getWithDefault("interpolation", d_interpType, "linear");
  ProblemSpecP interp = ps->findBlock("interpolation");
  if (!interp) {
    d_interpType = "linear";
  } else {
    if (!interp->getAttribute("type", d_interpType)) {
      throw ProblemSetupException("**ERROR** Interpolation \
        tag needs type=linear/cubic", __FILE__, __LINE__);
    }
  }
  initialize();
}

TabularData::TabularData(const TabularData& table) {
  *this = table;
}

TabularData& 
TabularData::operator=(const TabularData& table) {

  if (this != &table) {
    d_indepVarNames = table.d_indepVarNames;
    d_depVarNames = table.d_depVarNames;
    d_filename = table.d_filename;
    for (auto ii = 0u; ii < table.d_indepVars.size(); ii++) {
      addIndependentVariable(table.d_indepVars[ii]->name);
      d_indepVars[ii]->data = table.d_indepVars[ii]->data;
    }
    for (auto ii = 0u; ii < table.d_depVars.size(); ii++) {
      addDependentVariable(table.d_depVars[ii]->name);
      d_depVars[ii]->data = table.d_depVars[ii]->data;
    }
  }
  return *this;
}

void
TabularData::initialize() 
{
  std::vector<std::string> indepVarNames = parseVariableNames(d_indepVarNames);
  std::vector<std::string> depVarNames = parseVariableNames(d_depVarNames);

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

void
TabularData::outputProblemSpec(ProblemSpecP& ps)
{
  ps->appendElement("independent_variables", d_indepVarNames);
  ps->appendElement("dependent_variables", d_depVarNames);
  ps->appendElement("filename", d_filename);
  auto interp_ps = ps->appendChild("interpolation");
  interp_ps->setAttribute("type", d_interpType);
 
}

std::vector<std::string>
TabularData::parseVariableNames(const std::string& vars)
{
  std::vector<std::string> varNames = Vaango::Util::split(vars, ',');
  for (auto& name : varNames) {
    name = Vaango::Util::trim(name);
  }
  return varNames;
}

void
TabularData::addIndependentVariable(const std::string& varName)
{
  d_indepVars.push_back(std::make_unique<IndependentVar>(varName));
}

std::size_t
TabularData::addDependentVariable(const std::string& varName)
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
  return d_depVars.size() - 1;
}

void
TabularData::setup()
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
TabularData::readJSONTableFromFile(const std::string& tableFile)
{
  std::ifstream inputFile(tableFile);
  if (!inputFile) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot read tabular input data file " << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }

  std::stringstream inputStream;
  inputStream << inputFile.rdbuf();
  inputFile.close();

  json doc = loadJSON(inputStream, tableFile);
  readJSONTable<dim>(doc, tableFile);
}

json
TabularData::loadJSON(std::stringstream& inputStream,
                      const std::string& tableFile)
{
  json doc;
  try {
    inputStream >> doc;
  } catch (std::invalid_argument const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot parse tabular input data file " << tableFile << "\n"
        << " Please check that the file is valid JSON using a linter";
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return doc;
}

namespace Vaango {
template <>
void
TabularData::readJSONTable<1>(const json& doc, const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  DoubleVec1D indepVarData =
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({ IndexKey(0, 0, 0, 0), indepVarData });
  //std::cout << "Read " << d_indepVars[0]->name << " ";
  //std::copy(indepVarData.begin(), indepVarData.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  for (const auto& depVar : d_depVars) {
    DoubleVec1D depVarData = getDoubleArrayJSON(data, depVar->name, tableFile);
    depVar->data.insert({ IndexKey(0, 0, 0, 0), depVarData });
    //std::cout << "Read " << depVar->name << " ";
    //std::copy(depVarData.begin(), depVarData.end(),
    //          std::ostream_iterator<double>(std::cout, " "));
    //std::cout << std::endl;

    ASSERTEQ(indepVarData.size(), depVarData.size());
  }
}

template <>
void
TabularData::readJSONTable<2>(const json& doc, const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  DoubleVec1D indepVar0Data =
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({ IndexKey(0, 0, 0, 0), indepVar0Data });

  json data0 = getDataJSON(data, tableFile);
  for (auto ii = 0u; ii < indepVar0Data.size(); ii++) {

    DoubleVec1D indepVar1Data =
      getDoubleArrayJSON(data0[ii], d_indepVars[1]->name, tableFile);
    //std::cout << "Read " << d_indepVars[1]->name << " ii = " << ii << " ";
    //std::copy(indepVar1Data.begin(), indepVar1Data.end(),
    //          std::ostream_iterator<double>(std::cout, " "));
    //std::cout << std::endl;
    d_indepVars[1]->data.insert({ IndexKey(ii, 0, 0, 0), indepVar1Data });

    for (const auto& depVar : d_depVars) {
      DoubleVec1D depVarData =
        getDoubleArrayJSON(data0[ii], depVar->name, tableFile);
      depVar->data.insert({ IndexKey(ii, 0, 0, 0), depVarData });

      ASSERTEQ(indepVar1Data.size(), depVarData.size());
    }
  }
}

template <>
void
TabularData::readJSONTable<3>(const json& doc, const std::string& tableFile)
{
  json contents = getContentsJSON(doc, tableFile);
  std::string title = getTitleJSON(contents, tableFile);
  json data = getDataJSON(contents, tableFile);

  DoubleVec1D indepVar0Data =
    getDoubleArrayJSON(data, d_indepVars[0]->name, tableFile);
  d_indepVars[0]->data.insert({ IndexKey(0, 0, 0, 0), indepVar0Data });

  json data0 = getDataJSON(data, tableFile);
  for (auto ii = 0u; ii < indepVar0Data.size(); ii++) {

    DoubleVec1D indepVar1Data =
      getDoubleArrayJSON(data0[ii], d_indepVars[1]->name, tableFile);
    d_indepVars[1]->data.insert({ IndexKey(ii, 0, 0, 0), indepVar1Data });

    json data1 = getDataJSON(data0[ii], tableFile);
    for (auto jj = 0u; jj < indepVar1Data.size(); jj++) {

      int index = 0;
      for (const auto& indepVar : d_indepVars) {
        if (index > 1) {
          DoubleVec1D indepVar2Data =
            getDoubleArrayJSON(data1[jj], indepVar->name, tableFile);
          indepVar->data.insert({ IndexKey(ii, jj, 0, 0), indepVar2Data });
        }
        index++;
      }

      for (const auto& depVar : d_depVars) {
        DoubleVec1D depVar2Data =
          getDoubleArrayJSON(data1[jj], depVar->name, tableFile);
        depVar->data.insert({ IndexKey(ii, jj, 0, 0), depVar2Data });
      }
    }
  }
}
}

json
TabularData::getContentsJSON(const json& doc, const std::string& tableFile)
{
  json contents;
  try {
    contents = doc["Vaango_tabular_data"];
  } catch (std::out_of_range const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Cannot find the key \"Vaango_tabular_data\" in " << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return contents;
}

std::string
TabularData::getTitleJSON(const json& contents, const std::string& tableFile)
{
  std::string title;
  try {
    title = contents["Meta"].at("title").get<std::string>();
  } catch (std::out_of_range const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Could not find tabular EOS title in " << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return title;
}

json
TabularData::getDataJSON(const json& contents, const std::string& tableFile)
{
  json data;
  try {
    data = contents["Data"];
  } catch (std::out_of_range const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Key \"Data\" not found in tabular EOS input file " << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return data;
}

DoubleVec1D
TabularData::getVectorJSON(const json& object, const std::string key,
                           const std::string& tableFile)
{
  std::string str;
  try {
    str = object.at(key).get<std::string>();
  } catch (std::out_of_range const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" not found in tabular EOS input file "
        << tableFile;
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
TabularData::getDoubleArrayJSON(const json& object, const std::string key,
                                const std::string& tableFile)
{
  // First check if the key exists
  json data;
  try {
    data = object.at(key);
    // std::cout << "Data = " << data << std::endl;
  } catch (std::out_of_range const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" not found in tabular EOS input file "
        << tableFile;
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
  } catch (std::exception const& err) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Variable \"" << key << "\" contains invalid data. " << tableFile;
    throw ProblemSetupException(out.str(), __FILE__, __LINE__);
  }
  return vec;
}

namespace Vaango {
template <>
void
TabularData::translateIndepVar1ByIndepVar0<2>() 
{
  auto indepVarData0 =
    getIndependentVarData(d_indepVars[0]->name, IndexKey(0, 0, 0, 0));
  //std::cout << "Read " << d_indepVars[0]->name << " ";
  //std::copy(indepVarData0.begin(), indepVarData0.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  for (auto ii = 0u; ii < indepVarData0.size() ; ii++) {
    auto indepVarData1 =
      getIndependentVarData(d_indepVars[1]->name, IndexKey(ii, 0, 0, 0));
    //std::cout << "Read " << d_indepVars[1]->name << " ";
    //std::copy(indepVarData1.begin(), indepVarData1.end(),
    //          std::ostream_iterator<double>(std::cout, " "));
    //std::cout << std::endl;
    for (auto& val : indepVarData1) {
      val -= indepVarData0[ii];
    }
    d_indepVars[1]->data[IndexKey(ii, 0, 0, 0)] =  indepVarData1;
  }
}

template <>
void
TabularData::translateAlongNormals<1>(const std::vector<Vector>& normals,
                                      const double& shift) 
{
  auto indepVarData0 =
    getIndependentVarData(d_indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto depVarData0 =
    getDependentVarData(d_depVars[0]->name, IndexKey(0, 0, 0, 0));
  //std::cout << "Read " << d_indepVars[0]->name << " ";
  //std::copy(indepVarData0.begin(), indepVarData0.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  for (auto ii = 0u; ii < indepVarData0.size(); ii++) {
    indepVarData0[ii] += normals[ii].x()*shift;
    depVarData0[ii] += normals[ii].y()*shift;
  }
  d_indepVars[0]->data[IndexKey(0, 0, 0, 0)] = indepVarData0;
  d_depVars[0]->data[IndexKey(0, 0, 0, 0)] = depVarData0;

  //auto data =
  //  getIndependentVarData(d_indepVars[0]->name, IndexKey(0, 0, 0, 0));
  //std::cout << "Translated " << d_indepVars[0]->name << " ";
  //std::copy(data.begin(), data.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;
}
}

template <int dim>
DoubleVec1D
TabularData::findLimits(const std::array<double, dim>& indepValues) const
{
  return findLimitValues<dim>(indepValues, d_indepVars, d_depVars);
}

namespace Vaango {
template <>
DoubleVec1D
TabularData::findLimitValues<1>(
  const std::array<double, 1>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  sval = (sval < 0) ? 0 : sval;
  sval = (sval > 1) ? 1 : sval;

  DoubleVec1D indepVals;
  auto indepval = computeInterpolated(sval, segLowIndex0, indepVarData0);
  indepVals.push_back(indepval);

  return indepVals;
}

template <>
DoubleVec1D
TabularData::findLimitValues<2>(
  const std::array<double, 2>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  sval = (sval < 0) ? 0 : sval;
  sval = (sval > 1) ? 1 : sval;

  DoubleVec1D indepVals;
  auto indepval = computeInterpolated(sval, segLowIndex0, indepVarData0);
  indepVals.push_back(indepval);

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  DoubleVec1D depValsT;
  for (auto ii = segLowIndex0; ii <= segLowIndex0 + 1; ii++) {
    auto indepVarData1 =
      getIndependentVarData(indepVars[1]->name, IndexKey(ii, 0, 0, 0));
    indepVals.push_back(indepVarData1.front());
    indepVals.push_back(indepVarData1.back());
    auto segLowIndex1 = findLocationNoThrow(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);
    tval = (tval < 0) ? 0 : tval;
    tval = (tval > 1) ? 1 : tval;
    auto indepval = computeInterpolated(tval, segLowIndex1, indepVarData1);
    indepVals.push_back(indepval);
  }

  return indepVals;
}

template <>
DoubleVec1D
TabularData::findLimitValues<3>(
  const std::array<double, 3>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  sval = (sval < 0) ? 0 : sval;
  sval = (sval > 1) ? 1 : sval;

  DoubleVec1D indepVals;
  auto indepval = computeInterpolated(sval, segLowIndex0, indepVarData0);
  indepVals.push_back(indepval);

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  for (auto ii = segLowIndex0; ii <= segLowIndex0 + 1; ii++) {
    auto indepVarData1 =
      getIndependentVarData(indepVars[1]->name, IndexKey(ii, 0, 0, 0));
    auto segLowIndex1 = findLocationNoThrow(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);
    tval = (tval < 0) ? 0 : tval;
    tval = (tval > 1) ? 1 : tval;

    auto indepval = computeInterpolated(tval, segLowIndex1, indepVarData1);
    indepVals.push_back(indepval);

    // Choose the last two vectors containing the relevant independent variable
    // data
    // and find the segments containing the data
    for (auto jj = segLowIndex1; jj <= segLowIndex1 + 1; jj++) {
      auto indepVarData2 =
        getIndependentVarData(indepVars[2]->name, IndexKey(ii, jj, 0, 0));
      indepVals.push_back(indepVarData2.front());
      indepVals.push_back(indepVarData2.back());
    }
  }

  return indepVals;
}
} // end namespace for template specialization (needed by gcc)

template <int dim>
DoubleVec1D
TabularData::interpolate(const std::array<double, dim>& indepValues) const
{
  return interpolateLinearSpline<dim>(indepValues, d_indepVars, d_depVars);
}

namespace Vaango {
template <>
DoubleVec1D
TabularData::interpolateLinearSpline<1>(
  const std::array<double, 1>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  // auto numIndepVars = indepVars.size();
  // auto numDepVars = depVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);
  // std::cout << "Lo Index 0 = " << segLowIndex0 << " s = " << sval <<
  // std::endl;

  DoubleVec1D depVals;
  for (const auto& depVar : depVars) {
    auto depVarData = getDependentVarData(depVar->name, IndexKey(0, 0, 0, 0));
    auto depval = computeInterpolated(sval, segLowIndex0, depVarData);
    depVals.push_back(depval);
    // std::cout << "Lo Index 0 = " << segLowIndex0 << " p = " << depval <<
    // std::endl;
  }

  return depVals;
}

/*!
 * Assumes that data are in ascending order of the first independent variable
 * followed by ascending order of the second independent variable.
 * y1 = y(x11, x21)
 * y2 = y(x12, x22) 
 * ...
 * yp = y(x1p, x2p)
 *
 * x11 < x12 < ... < x1p, 
 *
 * x21 = {a1, a2, .., an}, a1 < a2 < ... < an
 * x22 = {b1, b2, .., bm}, b1 < b2 < ... < bm
 * ...
 * x2p = {q1, q2, .., qr}, q1 < q2 < ... < qr
 *
 * y1 = {y1_1, y1_2, ..., y1_n}
 * y2 = {y2_1, y2_2, ..., y2_m}
 * ...
 * yp = {yp_1, yp_2, ..., yp_r}
 */
template <>
DoubleVec1D
TabularData::interpolateLinearSpline<2>(
  const std::array<double, 2>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  // auto numIndepVars = indepVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));

  #ifdef DEBUG_TABLE_INTERPOLATION
  std::cout << "Read " << indepVars[0]->name << " ";
  //std::copy(indepVarData0.begin(), indepVarData0.end(),
  //          std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  #endif

  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);

  #ifdef DEBUG_TABLE_INTERPOLATION
  std::cout << "t = " << sval << "\n";
  if (sval < 0 || sval > 1) {
    std::cout << "Lo Index 0 = [" << segLowIndex0 << "," << (segLowIndex0+1) << "]"
              << " s = " << sval << " value = " << indepValues[0] 
              << " data = [" << indepVarData0[segLowIndex0]
              << "," << indepVarData0[segLowIndex0+1] << "]" << std::endl;
  }
  #endif

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  DoubleVec1D depValsT;
  for (auto ii = segLowIndex0; ii <= segLowIndex0 + 1; ii++) {
    auto indepVarData1 =
      getIndependentVarData(indepVars[1]->name, IndexKey(ii, 0, 0, 0));

    #ifdef DEBUG_TABLE_INTERPOLATION
    std::cout << "Read " << indepVars[1]->name << " ";
    //std::copy(indepVarData1.begin(), indepVarData1.end(),
    //          std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    #endif

    auto segLowIndex1 = findLocationNoThrow(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);

    #ifdef DEBUG_TABLE_INTERPOLATION
    std::cout << "s [" << ii << "]= " << tval << " start = " << segLowIndex1
              << " end = " << segLowIndex1 + 1 << "\n";
    if (sval < 0 || sval > 1) {
      if (tval < 0 || tval > 1) {
        std::cout << "Lo Index 1 = [" << segLowIndex1 << "," << (segLowIndex1+1) << "]"
                  << " t = " << tval << " value = " << indepValues[1] 
                  << " data = [" << indepVarData1[segLowIndex1]
                  << "," << indepVarData1[segLowIndex1+1] << "]" << std::endl;
      }
    }
    #endif

    for (const auto& depVar : depVars) {
      auto depVarData =
        getDependentVarData(depVar->name, IndexKey(ii, 0, 0, 0));
      auto depvalT = computeInterpolated(tval, segLowIndex1, depVarData);
      depValsT.push_back(depvalT);

      #ifdef DEBUG_TABLE_INTERPOLATION
      std::cout << "p = " << depvalT << "\n";
      if (sval < 0 || sval > 1) {
        if (tval < 0 || tval > 1) {
          std::cout << "Lo Index 1 = " << segLowIndex1 << " p = " << depvalT 
                    << " data = [" << depVarData[segLowIndex1]
                    << "," << depVarData[segLowIndex1+1] << "]" << std::endl;
        }
      }
      #endif
    }
  }

  auto numDepVars = depVars.size();
  DoubleVec1D depVals;
  for (auto index = 0u; index < numDepVars; index++) {
    auto depval =
      (1 - sval) * depValsT[index] + sval * depValsT[index + numDepVars];
    depVals.push_back(depval);

    #ifdef DEBUG_TABLE_INTERPOLATION
    if (sval < 0 || sval > 1) {
      std::cout << "Lo Index 0 = " << segLowIndex0 << " q = " << depval 
                << " data = [" << depValsT[index]
                << "," << depValsT[index+numDepVars] << "]" << std::endl;
    }
    #endif

  }

  return depVals;
}


template <>
DoubleVec1D
TabularData::interpolateLinearSpline<3>(
  const std::array<double, 3>& indepValues, const IndepVarPArray& indepVars,
  const DepVarPArray& depVars) const
{
  // Get the numbers of vars
  // auto numIndepVars = indepVars.size();
  auto numDepVars = depVars.size();

  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData0 =
    getIndependentVarData(indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto segLowIndex0 = findLocationNoThrow(indepValues[0], indepVarData0);
  auto sval = computeParameter(indepValues[0], segLowIndex0, indepVarData0);

  // Choose the two vectors containing the relevant independent variable data
  // and find the segments containing the data
  DoubleVec1D depValsT;
  for (auto ii = segLowIndex0; ii <= segLowIndex0 + 1; ii++) {
    auto indepVarData1 =
      getIndependentVarData(indepVars[1]->name, IndexKey(ii, 0, 0, 0));
    auto segLowIndex1 = findLocationNoThrow(indepValues[1], indepVarData1);
    auto tval = computeParameter(indepValues[1], segLowIndex1, indepVarData1);

    // Choose the last two vectors containing the relevant independent variable
    // data
    // and find the segments containing the data
    DoubleVec1D depValsU;
    for (auto jj = segLowIndex1; jj <= segLowIndex1 + 1; jj++) {
      auto indepVarData2 =
        getIndependentVarData(indepVars[2]->name, IndexKey(ii, jj, 0, 0));
      auto segLowIndex2 = findLocationNoThrow(indepValues[2], indepVarData2);
      auto uval = computeParameter(indepValues[2], segLowIndex2, indepVarData2);

      for (const auto& depVar : depVars) {
        auto depVarData =
          getDependentVarData(depVar->name, IndexKey(ii, jj, 0, 0));
        auto depvalU = computeInterpolated(uval, segLowIndex2, depVarData);
        depValsU.push_back(depvalU);
      }
    }

    for (auto index = 0u; index < numDepVars; index++) {
      auto depvalT =
        (1 - tval) * depValsU[index] + tval * depValsU[index + numDepVars];
      depValsT.push_back(depvalT);
    }
  }

  DoubleVec1D depVals;
  for (auto index = 0u; index < numDepVars; index++) {
    auto depval =
      (1 - sval) * depValsT[index] + sval * depValsT[index + numDepVars];
    depVals.push_back(depval);
  }

  return depVals;
}
} // end namespace for template specialization (needed by gcc)


double
TabularData::computeDerivative(double indepValue, double depValue) const
{
  // First find the segment containing the first independent variable value
  // and the value of parameter s
  auto indepVarData =
    getIndependentVarData(d_indepVars[0]->name, IndexKey(0, 0, 0, 0));
  auto index = findLocationNoThrow(indepValue, indepVarData);
  auto sval = computeParameter(indepValue, index, indepVarData);

  // Choose the two dependent variable vectors that will be used for interpolation
  auto depVarData0 = getDependentVarData(d_depVars[0]->name, IndexKey(index, 0, 0, 0));
  auto depVarData1 = getDependentVarData(d_depVars[0]->name, IndexKey(index+1, 0, 0, 0));

  // Choose the two independent variable vectors that will be used for interpolation
  auto indepVarData0 = getIndependentVarData(d_indepVars[1]->name, IndexKey(index, 0, 0, 0));
  auto indepVarData1 = getIndependentVarData(d_indepVars[1]->name, IndexKey(index+1, 0, 0, 0));
  auto minmax_indep0 = std::minmax_element(indepVarData0.begin(), indepVarData0.end());
  auto minmax_indep1 = std::minmax_element(indepVarData1.begin(), indepVarData1.end());

  // Set up line segments for intersection
  auto indepLo = std::min(*(minmax_indep0.first), *(minmax_indep1.first));
  auto indepHi = std::max(*(minmax_indep0.second), *(minmax_indep1.second));
  auto depLo = depValue - 1.0e-6*depValue;
  auto depHi = depValue + 1.0e-6*depValue;
  depLo = (depLo == 0) ? -1.0e-6 : depLo;
  depHi = (depHi == 0) ?  1.0e-6 : depHi;

  // Find the intersection points
  bool status;
  double t_p, t_q;

  double indepValLo0, depValLo0;
  std::tie(status, index, t_p, t_q, indepValLo0, depValLo0) = 
   Vaango::Util::findIntersectionTableLinearSearch(indepVarData0, depVarData0,
                                                   indepLo, depLo, indepHi, depLo);
  //std::cout << "status = " << status << " index = " << index 
  //          << " t_p = " << t_p << " t_q = " << t_q
  //          << " x0 = " << indepValLo0 << " y0 = " << depValLo0 
  //          << " seg = [[" << indepLo << "," << depLo << "],["
  //                         << indepHi << "," << depLo << "]\n";
  
  double indepValLo1, depValLo1;
  std::tie(status, index, t_p, t_q, indepValLo1, depValLo1) = 
   Vaango::Util::findIntersectionTableLinearSearch(indepVarData1, depVarData1,
                                                   indepLo, depLo, indepHi, depLo);
  //std::cout << "status = " << status << " index = " << index 
  //          << " t_p = " << t_p << " t_q = " << t_q
  //          << " x0 = " << indepValLo1 << " y0 = " << depValLo1 
  //          << " seg = [[" << indepLo << "," << depLo << "],["
  //                         << indepHi << "," << depLo << "]\n";

  double indepValHi0, depValHi0;
  std::tie(status, index, t_p, t_q, indepValHi0, depValHi0) = 
   Vaango::Util::findIntersectionTableLinearSearch(indepVarData0, depVarData0,
                                                   indepLo, depHi, indepHi, depHi);
  //std::cout << "status = " << status << " index = " << index 
  //          << " t_p = " << t_p << " t_q = " << t_q
  //          << " x0 = " << indepValHi0 << " y0 = " << depValHi0 
  //          << " seg = [[" << indepLo << "," << depHi << "],["
  //                         << indepHi << "," << depHi << "]\n";

  double indepValHi1, depValHi1;
  std::tie(status, index, t_p, t_q, indepValHi1, depValHi1) = 
   Vaango::Util::findIntersectionTableLinearSearch(indepVarData1, depVarData1,
                                                   indepLo, depHi, indepHi, depHi);
  //std::cout << "status = " << status << " index = " << index 
  //          << " t_p = " << t_p << " t_q = " << t_q
  //          << " x0 = " << indepValHi1 << " y0 = " << depValHi1 
  //          << " seg = [[" << indepLo << "," << depHi << "],["
  //                         << indepHi << "," << depHi << "]\n";


  // Interpolate between the high and low values
  //std::cout << "ep_v(0,0) = " << indepValLo0 << " ep_v(1,0) = " << indepValLo1
  //          << "ep_v(0,1) = " << indepValHi0 << " ep_v(1,1) = " << indepValHi1 << "\n";
  auto indepValLo = (1 - sval)*indepValLo0 + sval*indepValLo1;
  auto indepValHi = (1 - sval)*indepValHi0 + sval*indepValHi1;

  // Compute derivative
  auto d_dep = depHi - depLo;
  auto d_indep = indepValHi - indepValLo;
  auto deriv = d_dep/d_indep;

  return deriv;
}

std::size_t
TabularData::findLocation(const double& value, const DoubleVec1D& varData) const
{
  if (value < varData.front() || value > varData.back()) {
    std::ostringstream out;
    out << "**WARNING**"
        << " The independent variable value \"" 
        << std::setprecision(10) << std::scientific << value
        << "\" is outside the range of known data ["
        << varData.front() << "," << varData.back() << "]";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  auto lower = std::lower_bound(varData.begin(), varData.end(), value);
  auto index = (lower == varData.begin()) 
               ? lower - varData.begin() 
               : lower - varData.begin() - 1;
  return index;
}

std::size_t
TabularData::findLocationNoThrow(const double& value, const DoubleVec1D& varData) const
{
  auto second = std::lower_bound(varData.begin(), varData.end(), value);
  auto first = second-1;
  if (second == varData.end()) {
    first = varData.end() - 2;
    second = varData.end() - 1;
  } else if (second == varData.begin()) {
    first = second;
    second = first+1;
  }

  auto index = first - varData.begin();
  return index;
}

double
TabularData::computeParameter(const double& input,
                              const std::size_t& startIndex,
                              const DoubleVec1D& data) const
{
  if (data.size() < 2)
    return 0.0;
  auto t =
    (input - data[startIndex]) / (data[startIndex + 1] - data[startIndex]);
  //std::cout << "input = " << input << " start_val = " << data[startIndex]
  //          << " end_val = " << data[startIndex+1] << "\n";
  return t;
}

double
TabularData::computeInterpolated(const double& tval,
                                 const std::size_t& startIndex,
                                 const DoubleVec1D& data) const
{
  if (data.empty()) {
    std::ostringstream out;
    out << "**ERROR**"
        << " No data available for interpolation.";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  return (data.size() == 1)
           ? data.front()
           : (1 - tval) * data[startIndex] + tval * data[startIndex + 1];
}

// See: https://cs.uwaterloo.ca/research/tr/1983/CS-83-09.pdf
/*
double
TabularData::interpolateCubicSpline1D(const double& t) {
}
*/

/*
 * Find the intersection of 1D tabular data with a line segment
 *  assuming x = independent variable
 *           y = dependent variable
 *  indepLo = x_lo of line segment
 *  indepHi = x_hi of line sgement
 *  depLo = y_lo of line segment
 *  depHi = y_hi of line sgement
 *
 * Returns:
 *  bool : status indicating if there is an intersection or not
 *  double: x - location of intersection
 *  double: y - location of intersection
 */
std::tuple<bool, double, double>
TabularData::intersect1D(double indepLo, double indepHi,
                         double depHi, double depLo)
{
  // Check that the table contains 1D data
  if (getNumIndependents() + getNumDependents() > 2) {
    std::ostringstream out;
    out << "**ERROR**"
        << " Intersection of a line with tabular data no allowed"
        << " if there is more than one independent and one dependent variable";
    throw InvalidValue(out.str(), __FILE__, __LINE__);
  }

  // Get the data 
  auto x_coords = (*d_indepVars.begin())->data.begin()->second;
  auto y_coords = (*d_depVars.begin())->data.begin()->second;

  // Find the intersection point
  bool status;
  std::size_t index;
  double t_p, t_q, indepVal, depVal;
  std::tie(status, index, t_p, t_q, indepVal, depVal) = 
   Vaango::Util::findIntersectionTableLinearSearch(x_coords, y_coords,
                                                   indepLo, depLo,
                                                   indepHi, depHi);
  return std::make_tuple(status, indepVal, depVal);
}


DoubleVec1D
TabularData::getIndependentVarData(const std::string& name,
                                   const IndexKey& index) const
{
  auto varIter = std::find_if(
    d_indepVars.begin(), d_indepVars.end(),
    [&name](const IndependentVarP& indepVar) { return (indepVar->name == name); });
  if (varIter == d_indepVars.end()) {
    return DoubleVec1D();
  }

  return (*varIter)->data.at(index);
}

DoubleVec1D
TabularData::getDependentVarData(const std::string& name,
                                 const IndexKey& index) const
{
  auto varIter = std::find_if(
    d_depVars.begin(), d_depVars.end(),
    [&name](const DependentVarP& depVar) { return (depVar->name == name); });
  if (varIter == d_indepVars.end()) {
    return DoubleVec1D();
  }

  return (*varIter)->data.at(index);
}

void 
TabularData::setIndependentVarData(const std::string& name,
                                   const IndexKey& index,
                                   const DoubleVec1D& data)
{
  auto varIter = std::find_if(
    d_indepVars.begin(), d_indepVars.end(),
    [&name](const IndependentVarP& indepVar) { return (indepVar->name == name); });
  (*varIter)->data[index] = data;
}

void 
TabularData::setDependentVarData(const std::string& name,
                                 const IndexKey& index,
                                 const DoubleVec1D& data)
{
  auto varIter = std::find_if(
    d_depVars.begin(), d_depVars.end(),
    [&name](const DependentVarP& depVar) { return (depVar->name == name); });
  (*varIter)->data[index] = data;
}

namespace Vaango {
template DoubleVec1D
TabularData::interpolate<1>(const std::array<double, 1>& indepValues) const;
template DoubleVec1D
TabularData::interpolate<2>(const std::array<double, 2>& indepValues) const;
template DoubleVec1D
TabularData::interpolate<3>(const std::array<double, 3>& indepValues) const;
template DoubleVec1D
TabularData::findLimits<1>(const std::array<double, 1>& indepValues) const;
template DoubleVec1D
TabularData::findLimits<2>(const std::array<double, 2>& indepValues) const;
template DoubleVec1D
TabularData::findLimits<3>(const std::array<double, 3>& indepValues) const;
}
