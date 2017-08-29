#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableBase.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <submodules/json/src/json.hpp>

#include <string>
#include <vector>
#include <memory>

namespace Vaango {

  using IndependentVarP = std::unique_ptr<TableBase::IndependentVar>;
  using DependentVarP = std::unique_ptr<TableBase::DependentVar>;

  class TableEOS : public TableBase {

  public:
    TableEOS() = delete;
    TableEOS(Uintah::ProblemSpecP& ps);
    ~TableEOS() override = default;

    void addIndependentVariable(const std::string& varName) override;
    std::size_t addDependentVariable(const std::string& varName) override;

    void setup() override;

    double interpolate(const std::string& depvarName,
                       std::vector<double>& independents) override;

    template <int dim>
    void readJSONTableFromFile(const std::string& tableFile);

    template <int dim>
    void readJSONTable(const nlohmann::json& doc,
                       const std::string& tableFile);

    double interpolateLinearSpline1D(const double& indepValue,
                                     const std::vector<double>& indepVar,
                                     const std::vector<double>& depVar) const;

    std::vector<double> fitCubicSpline1D(const DependentVar depVar,
                                         const IndependentVar indepVar);

    double interpolateCubicSpline1D(const double& t);

    std::size_t getNumIndependents() const {return d_indepVars.size();}
    std::size_t getNumDependents() const {return d_depVars.size();}

    std::vector<double> getIndependentVarData(const std::string& name,
                                              const IndexKey& index);
    std::vector<double> getDependentVarData(const std::string& name,
                                            const IndexKey& index);

    private:
      std::string d_filename;
      std::vector<IndependentVarP> d_indepVars;
      std::vector<DependentVarP> d_depVars;

      std::vector<std::string> parseVariableNames(const std::string& vars);

      nlohmann::json loadJSON(std::stringstream& inputStream,
                              const std::string& fileName);
      nlohmann::json getContentsJSON(const nlohmann::json& doc,
                                     const std::string& fileName);
      std::string getTitleJSON(const nlohmann::json& contents,
                               const std::string& fileName);
      nlohmann::json getDataJSON(const nlohmann::json& contents,
                                 const std::string& fileName);
      std::vector<double> getVectorJSON(const nlohmann::json& object,
                                        const std::string key,
                                        const std::string& tableFile);
      std::vector<double> getDoubleArrayJSON(const nlohmann::json& object,
                                             const std::string key,
                                             const std::string& tableFile);
  };
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
