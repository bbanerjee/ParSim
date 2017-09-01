#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_CONTAINERS_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_CONTAINERS_H

#include <string>
#include <unordered_map>
#include <vector>

namespace Vaango {

class TableContainers
{

public:
  TableContainers() = default;
  ~TableContainers() = default;

  struct IndexKey
  {
    std::uint8_t _ii;
    std::uint8_t _jj;
    std::uint8_t _kk;
    std::uint8_t _ll;
    IndexKey(std::uint8_t ii, std::uint8_t jj, std::uint8_t kk, std::uint8_t ll)
      : _ii(ii)
      , _jj(jj)
      , _kk(kk)
      , _ll(ll)
    {
    }
  };

  struct IndexEqual
  {
    bool operator()(const IndexKey& lhs, const IndexKey& rhs) const
    {
      return (lhs._ii == rhs._ii && lhs._jj == rhs._jj && lhs._kk == rhs._kk &&
              lhs._ll == rhs._ll);
    }
  };

  struct IndexHash
  {
    std::size_t operator()(const IndexKey& key) const
    {
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

  struct IndependentVar
  {
    std::string name;
    std::unordered_map<IndexKey, std::vector<double>, IndexHash, IndexEqual>
      data;
    IndependentVar() = delete;
    IndependentVar(const std::string name) { this->name = name; }
  };

  struct DependentVar
  {
    std::string name;
    std::unordered_map<IndexKey, std::vector<double>, IndexHash, IndexEqual>
      data;
    DependentVar() = delete;
    DependentVar(const std::string name) { this->name = name; }
  };
};
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_CONTAINERS_H
