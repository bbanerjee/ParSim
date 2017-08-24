#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_EOS_H

#include <CCA/Components/MPM/ConstitutiveModel/Models/TableBase.h>
#include <Core/ProblemSpec/ProblemSpec.h>

#include <submodules/json/src/json.hpp>

#include <string>
#include <vector>
#include <unordered_map>
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

    template <int dim>
    void readJSONTableFromFile(const std::string& tableFile);

    template <int dim>
    void readJSONTable(const nlohmann::json& doc,
                       const std::string& tableFile);

    private:
      std::string d_filename;

      struct IndexKey {
        std::uint8_t _ii;
        std::uint8_t _jj;
        std::uint8_t _kk;
        std::uint8_t _ll;
        IndexKey(std::uint8_t ii, std::uint8_t jj, 
                 std::uint8_t kk, std::uint8_t ll)
                 : _ii(ii), _jj(jj), _kk(kk), _ll(ll)
        {
        }
      };

      struct IndexEqual {
        bool operator()(const IndexKey& lhs, const IndexKey& rhs) const {
          return (lhs._ii == rhs._ii && lhs._jj == rhs._jj && 
                  lhs._kk == rhs._kk && lhs._ll == rhs._ll);
        }
      };

      struct IndexHash {
        std::size_t operator()(const IndexKey& key) const {
          std::size_t hashval = key._ii;
          hashval *= 37;
          hashval += key._jj;
          hashval *= 37;
          hashval += key._kk;
          hashval *= 37;
          hashval += key._ll;
          return hashval;
        }
      };

      struct IndependentVar {
        std::string name;
        std::unordered_map<IndexKey, std::vector<double>, IndexHash, IndexEqual> data;
        IndependentVar() = delete;
        IndependentVar(const std::string name) {this->name = name;}
      };
      using IndependentVarP = std::unique_ptr<IndependentVar>;
      std::vector<IndependentVarP> d_indepVars;

      struct DependentVar {
        std::string name;
        std::unordered_map<IndexKey, std::vector<double>, IndexHash, IndexEqual> data;
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
