#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableBase.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <submodules/json/src/json.hpp>

#include <string>
#include <vector>
#include <array>
#include <memory>

namespace Vaango {

  using DoubleVec1D = std::vector<double>;
  using DoubleVec2D = std::vector<DoubleVec1D>;

  using IndependentVarP = std::unique_ptr<TableBase::IndependentVar>;
  using DependentVarP = std::unique_ptr<TableBase::DependentVar>;

  using IndepVarPArray = std::vector<IndependentVarP>;
  using DepVarPArray = std::vector<DependentVarP>;

  class TableEOS : public TableBase {

  public:
    TableEOS() = delete;
    TableEOS(Uintah::ProblemSpecP& ps);
    ~TableEOS() override = default;

    void addIndependentVariable(const std::string& varName) override;
    std::size_t addDependentVariable(const std::string& varName) override;

    void setup() override;

    double interpolate(const std::string& depvarName,
                       DoubleVec1D& independents) override;

    template <int dim>
    void readJSONTableFromFile(const std::string& tableFile);

    template <int dim>
    void readJSONTable(const nlohmann::json& doc,
                       const std::string& tableFile);

    template <int dim>
    DoubleVec1D interpolateLinearSpline(const std::array<double,dim>& indepValues,
                                        const IndepVarPArray& indepVars,
                                        const DepVarPArray& depVars) const;

    DoubleVec1D fitCubicSpline1D(const DependentVar depVar,
                                         const IndependentVar indepVar);

    double interpolateCubicSpline1D(const double& t);

    std::size_t getNumIndependents() const {return d_indepVars.size();}
    std::size_t getNumDependents() const {return d_depVars.size();}

    DoubleVec1D getIndependentVarData(const std::string& name,
                                      const IndexKey& index) const;
    DoubleVec1D getDependentVarData(const std::string& name,
                                    const IndexKey& index) const;

    const IndepVarPArray& getIndependentVars() const {return d_indepVars;} 
    const DepVarPArray& getDependentVars() const {return d_depVars;} 

    private:
      std::string d_filename;
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
      DoubleVec1D getVectorJSON(const nlohmann::json& object,
                                        const std::string key,
                                        const std::string& tableFile);
      DoubleVec1D getDoubleArrayJSON(const nlohmann::json& object,
                                             const std::string key,
                                             const std::string& tableFile);
      std::size_t findLocation(const double& value,
                               const DoubleVec1D& varData) const;
      double computeParameter(const double& input,
                              const std::size_t& startIndex,
                              const DoubleVec1D& data) const;
      double computeInterpolated(const double& tval,
                                 const std::size_t& startIndex,
                                 const DoubleVec1D& data) const;
  };
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
