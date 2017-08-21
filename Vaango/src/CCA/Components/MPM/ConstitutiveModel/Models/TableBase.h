#ifndef VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_BASE_H
#define VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_BASE_H

#include <string>
#include <vector>

namespace Vaango {

  class TableBase {

  public:
    TableBase() = default;
    virtual ~TableBase() = default;

    virtual void addIndependentVariable(const std::string& varName) = 0;
    virtual int  addDependentVariable(const std::string& varName) = 0;

    virtual void setup() = 0;

    virtual double interpolate(int index,
                               std::vector<double>& independents) = 0;

  };
}

#endif // VAANGO_MPM_CONSTITUTIVE_MODEL_TABLE_BASE_H
