#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABULAR_DATA_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABULAR_DATA_H

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableContainers.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <submodules/json/src/json.hpp>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace Vaango {

using IndexKey = TableContainers::IndexKey;

using DoubleVec1D = std::vector<double>;
using DoubleVec2D = std::vector<DoubleVec1D>;

using IndependentVarP = std::unique_ptr<TableContainers::IndependentVar>;
using DependentVarP = std::unique_ptr<TableContainers::DependentVar>;

using IndepVarPArray = std::vector<IndependentVarP>;
using DepVarPArray = std::vector<DependentVarP>;

class TabularData
{

public:
  TabularData() {};
  TabularData(Uintah::ProblemSpecP& ps);
  TabularData(const TabularData& table);
  ~TabularData() = default;
  TabularData& operator=(const TabularData& table);

  void initialize();
  void outputProblemSpec(Uintah::ProblemSpecP& ps);

  void addIndependentVariable(const std::string& varName);
  std::size_t addDependentVariable(const std::string& varName);

  void setup();

  template <int dim>
  DoubleVec1D interpolate(const std::array<double, dim>& indepValues) const;

  template <int dim>
  void readJSONTableFromFile(const std::string& tableFile);

  template <int dim>
  void readJSONTable(const nlohmann::json& doc, const std::string& tableFile);

  template <int dim>
  DoubleVec1D interpolateLinearSpline(
    const std::array<double, dim>& indepValues, const IndepVarPArray& indepVars,
    const DepVarPArray& depVars) const;

  // double interpolateCubicSpline1D(const double& t);

  std::size_t getNumIndependents() const { return d_indepVars.size(); }
  std::size_t getNumDependents() const { return d_depVars.size(); }

  DoubleVec1D getIndependentVarData(const std::string& name,
                                    const IndexKey& index) const;
  DoubleVec1D getDependentVarData(const std::string& name,
                                  const IndexKey& index) const;

  const IndepVarPArray& getIndependentVars() const { return d_indepVars; }
  const DepVarPArray& getDependentVars() const { return d_depVars; }

private:
  std::string d_filename;
  std::string d_indepVarNames;
  std::string d_depVarNames;
  std::string d_interpType;

  IndepVarPArray d_indepVars;
  DepVarPArray d_depVars;

  std::vector<std::string> parseVariableNames(const std::string& vars);

  nlohmann::json loadJSON(std::stringstream& inputStream,
                          const std::string& fileName);
  nlohmann::json getContentsJSON(const nlohmann::json& doc,
                                 const std::string& fileName);
  std::string getTitleJSON(const nlohmann::json& contents,
                           const std::string& fileName);
  nlohmann::json getDataJSON(const nlohmann::json& contents,
                             const std::string& fileName);
  DoubleVec1D getVectorJSON(const nlohmann::json& object, const std::string key,
                            const std::string& tableFile);
  DoubleVec1D getDoubleArrayJSON(const nlohmann::json& object,
                                 const std::string key,
                                 const std::string& tableFile);
  std::size_t findLocation(const double& value,
                           const DoubleVec1D& varData) const;
  double computeParameter(const double& input, const std::size_t& startIndex,
                          const DoubleVec1D& data) const;
  double computeInterpolated(const double& tval, const std::size_t& startIndex,
                             const DoubleVec1D& data) const;
};
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABULAR_DATA_H
