#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableBase.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <submodules/json/src/json.hpp>

#include <string>
#include <vector>
#include <memory>

namespace Vaango {

  class TableEOS : public TableBase {

  public:
    TableEOS() = delete;
    TableEOS(Uintah::ProblemSpecP& ps);
    ~TableEOS() override = default;

    void addIndependentVariable(const std::string& varName) override;
    std::size_t addDependentVariable(const std::string& varName) override;

    void setup() override;

    double interpolate(int index,
                       std::vector<double>& independents) override;

    void readJSONTableFromFile(const std::string& tableFile);

    void readJSONTable(const nlohmann::json& doc,
                       const std::string& tableFile);

    private:
      std::string d_filename;

      struct InterpolationAxis {
        std::vector<double> weights;
        std::vector<long> offset;
      };
      using InterpolationAxisP = std::unique_ptr<InterpolationAxis>;

      struct IndependentVar {
        std::string name;
        IndependentVar() = delete;
        IndependentVar(const std::string name) {this->name = name;}
      };
      using IndependentVarP = std::unique_ptr<IndependentVar>;
      std::vector<IndependentVarP> d_indepVars;

      struct DependentVar {
        std::string name;
        std::vector<IndependentVarP> independents;
        std::vector<InterpolationAxisP> axes;
        double* data;
        DependentVar() = delete;
        DependentVar(const std::string name) {this->name = name;}
      };
      using DependentVarP = std::unique_ptr<DependentVar>;
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
