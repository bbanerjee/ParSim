/*
 * The MIT License
 *
 * Copyright (c) 2015-2023 Biswajit Banerjee
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

#ifndef ___ELASTIC_MODULI_SUPPORT_VECTOR_MODEL_H__
#define ___ELASTIC_MODULI_SUPPORT_VECTOR_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <submodules/json/single_include/nlohmann/json.hpp>

#include <limits>
#include <iomanip>

namespace Vaango {

using MatrixN1d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using MatrixN2d = Eigen::Matrix<double, Eigen::Dynamic, 2>;

/** 
 * ElasticModuli from a model fitted using support vector regression
 */
class ElasticModuli_SupportVector : public ElasticModuliModel
{

public:

  ElasticModuli_SupportVector() = delete;
  ElasticModuli_SupportVector(const ElasticModuli_SupportVector& smm) = delete;
  ~ElasticModuli_SupportVector() = default;

  ElasticModuli_SupportVector(Uintah::ProblemSpecP& ps);
  ElasticModuli_SupportVector(const ElasticModuli_SupportVector* smm);
  ElasticModuli_SupportVector& operator=(const ElasticModuli_SupportVector& smm) = delete;

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
  ElasticModuli getCurrentElasticModuli(const ModelStateBase* state) const override;

  ElasticModuli getElasticModuliLowerBound() const override
  {
    return getInitialElasticModuli();
  }
  ElasticModuli getElasticModuliUpperBound() const override
  {
    return ElasticModuli(std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::max());
  }

  /* Get the elastic moduli and their derivatives with respect to a single
     plastic internal variable */
  std::pair<ElasticModuli, ElasticModuli>
    getElasticModuliAndDerivatives(const ModelStateBase* state_input) const override;

  /*! Compute derivatives of moduli with respect to internal variables */
  std::vector<ElasticModuli> computeDModuliDIntVar(const ModelStateBase* state) const override;

  /*! Compute moduli and derivatives of moduli with respect to internal variables */
  std::pair<ElasticModuli, std::vector<ElasticModuli>>
  computeModuliAndDModuliDIntVar(const ModelStateBase* state) const override;

private:

  /* Tangent bulk modulus parameters */
  struct BulkModulusParameters
  {
    std::string filename;
    Eigen::Vector2d strain_conversion_factor;
    Eigen::Vector2d strain_scale_factor;
    Eigen::Vector2d strain_min;
    Eigen::Vector2d strain_max;
    double pressure_conversion_factor;
    double pressure_scale_factor;
    double pressure_min;
    double pressure_max;
    double rbf_kernel_gamma;
    MatrixN2d svr_support_vectors;
    MatrixN1d svr_dual_coeffs;
    //Eigen::Matrix<double, Eigen::Dynamic, 2> svr_support_vectors;
    //Eigen::Matrix<double, Eigen::Dynamic, 1> svr_dual_coeffs;
    double svr_intercept;
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

  double computeBulkModulus(const double& elasticVolStrain,
                            const double& plasticVolStrain) const;
  double computeBulkModulusPressure(double pressure,
                                    double plasticVolStrain) const;
  double computeShearModulus(const double& bulkModulus) const;

  /* Functions for handling JSON files */
  void
  readJSONFile(const std::string& filename);

  nlohmann::json 
  loadJSON(std::stringstream& inputStream,
           const std::string& filename) const;

  double 
  getDoubleJSON(const nlohmann::json& object, 
                const std::string key,
                const std::string& filename) const;

  double 
  getVector1dJSON(const nlohmann::json& object, 
                  const std::string key,
                  const std::string& filename) const;

  Eigen::Vector2d 
  getVector2dJSON(const nlohmann::json& object, 
                  const std::string key,
                  const std::string& filename) const;

  //Eigen::Matrix<double, Eigen::Dynamic, 1>
  MatrixN1d
  getMatrixN1dJSON(const nlohmann::json& object,
                   const std::string key,
                   const std::string& filename) const;

  //Eigen::Matrix<double, Eigen::Dynamic, 2>
  MatrixN2d
  getMatrixN2dJSON(const nlohmann::json& object,
                   const std::string key,
                   const std::string& filename) const;

};
} // End namespace Vaango

#endif // __ELASTIC_MODULI_SUPPORT_VECTOR_MODEL_H__
