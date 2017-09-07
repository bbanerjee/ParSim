/*
 * The MIT License
 *
 * Copyright (c) 2015-2107 Parresia Research Limited, New Zealand
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

#ifndef ___ELASTIC_MODULI_TABULAR_MODEL_H__
#define ___ELASTIC_MODULI_TABULAR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/TabularData.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelStateBase.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Air.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Granite.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/Pressure_Water.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <limits>

namespace Vaango {

/*! \class ElasticModuli_Tabular
 *  \brief A tabular elasticity model
 *  \author Biswajit Banerjee,
 */
class ElasticModuli_Tabular : public ElasticModuliModel
{

public:

  /*! Construct a constant elasticity model. */
  ElasticModuli_Tabular(Uintah::ProblemSpecP& ps);

  ElasticModuli_Tabular(const ElasticModuli_Tabular& smm) = delete;

  /*! Construct a copy of constant elasticity model. */
  ElasticModuli_Tabular(const ElasticModuli_Tabular* smm);

  /*! Destructor of constant elasticity model.   */
  ~ElasticModuli_Tabular() = default;

  ElasticModuli_Tabular& operator=(const ElasticModuli_Tabular& smm) = delete;

  void outputProblemSpec(Uintah::ProblemSpecP& ps) override;

  /*! Get parameters */
  std::map<std::string, double> getParameters() const override
  {
    std::map<std::string, double> params;
    params["G0"] = d_shear.G0;
    params["nu"] = d_shear.nu;
    return params;
  }

  /*! Compute the elasticity */
  ElasticModuli getInitialElasticModuli() const override;
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) override;

  ElasticModuli getElasticModuliLowerBound() const override
  {
    return getInitialElasticModuli();
  }
  ElasticModuli getElasticModuliUpperBound() const override
  {
    return ElasticModuli(std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
  }

private:

  /* Tangent bulk modulus parameters */
  struct BulkModulusParameters
  {
    TabularData table;
    BulkModulusParameters() = default;
    BulkModulusParameters(Uintah::ProblemSpecP& ps) : table(ps) {
      table.setup();
    }
    BulkModulusParameters(const BulkModulusParameters& bulk) {
      table = bulk.table;
    }
    BulkModulusParameters&
    operator=(const BulkModulusParameters& bulk) {
      if (this != &bulk) {
        table = bulk.table;
      }
      return *this;
    }
  };

  /* Tangent shear modulus parameters */
  struct ShearModulusParameters
  {
    double G0;
    double nu;
  };

  BulkModulusParameters d_bulk;
  ShearModulusParameters d_shear;

  void checkInputParameters();

  double computeBulkModulus(const double& totalStrain,
                            const double& plasticStrain) const;
  double computeShearModulus(const double& bulkModulus) const;

};
} // End namespace Uintah

#endif // __ELASTIC_MODULI_TABULAR_MODEL_H__
