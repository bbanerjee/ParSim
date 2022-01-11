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

#ifndef ___ELASTIC_MODULI_NEURAL_NET_BULK_MODULUS_MODEL_H__
#define ___ELASTIC_MODULI_NEURAL_NET_BULK_MODULUS_MODEL_H__

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuliModel.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelStateBase.h>
#include <Core/ProblemSpec/ProblemSpecP.h>

#include <Eigen/Core>
#include <submodules/json/single_include/nlohmann/json.hpp>
#include "H5Cpp.h"

#include <limits>
#include <iomanip>

namespace Vaango {

/*! \class ElasticModuli_NeuralNet_Bulk
 *  \brief A fully connected multilayer perceptron neural net bulk modulus model
 *         with shear modulus either constant or related by a Poisson's ratio
 *         to the bulk modulus
 *  \author Biswajit Banerjee,
 */
class ElasticModuli_NeuralNet_Bulk : public ElasticModuliModel
{

public:

  ElasticModuli_NeuralNet_Bulk() = delete;
  ElasticModuli_NeuralNet_Bulk(const ElasticModuli_NeuralNet_Bulk& smm) = delete;
  ~ElasticModuli_NeuralNet_Bulk() = default;

  ElasticModuli_NeuralNet_Bulk(Uintah::ProblemSpecP& ps);
  ElasticModuli_NeuralNet_Bulk(const ElasticModuli_NeuralNet_Bulk* smm);
  ElasticModuli_NeuralNet_Bulk& operator=(const ElasticModuli_NeuralNet_Bulk& smm) = delete;

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

  template <typename T>
  struct NeuralNetworkLayer {
    using EigenMatrixRowMajor = 
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    std::string name;
    std::string activation;
    int units;
    int input_size;
    EigenMatrixRowMajor weights;
    EigenMatrixRowMajor bias; 
  };

  template <typename T>
  struct NeuralNetworkModel {
    std::vector<NeuralNetworkLayer<T>> layers;
    std::vector<double> meanStrain;
    std::vector<double> stdStrain;
    double meanBulkModulus;
    double stdBulkModulus;

    int inputLayerSize() const { return layers[0].input_size; } 
    void readNeuralNetworkHDF5(const std::string& filename);
    double predict(double elasticVolStrain, double plasticVolStrain) const;
  };

  /* Tangent bulk modulus parameters */
  template <typename T>
  struct BulkModulusParameters
  {
    std::string filename;
    NeuralNetworkModel<T> model;

    BulkModulusParameters() = default;
    BulkModulusParameters(Uintah::ProblemSpecP& ps) {
      ps->require("filename", filename);
      double strain_scaling;
      ps->require("mean_elastic_strain", strain_scaling);
      model.meanStrain.push_back(strain_scaling);
      ps->require("mean_plastic_strain", strain_scaling);
      model.meanStrain.push_back(strain_scaling);
      ps->require("std_dev_elastic_strain", strain_scaling);
      model.stdStrain.push_back(strain_scaling);
      ps->require("std_dev_plastic_strain", strain_scaling);
      model.stdStrain.push_back(strain_scaling);
      ps->require("mean_bulk_modulus", model.meanBulkModulus);
      ps->require("std_dev_bulk_modulus", model.stdBulkModulus);

      /*
      std::cout << "Inputs: \n"
                << "\t filename = " << filename << "\n"
                << "\t eps_e_mean = " << model.meanStrain[0] << "\n"
                << "\t eps_p_mean = " << model.meanStrain[1] << "\n"
                << "\t eps_e_std = " << model.stdStrain[0] << "\n"
                << "\t eps_p_std = " << model.stdStrain[1] << "\n"
                << "\t K_mean = " << model.meanBulkModulus << "\n"
                << "\t K_std = " << model.stdBulkModulus << "\n";
      */

      model.readNeuralNetworkHDF5(filename);
      //std::cout << "Completed reading model\n";
    }

    BulkModulusParameters(const BulkModulusParameters& bulk) {
      filename = bulk.filename;
      model = bulk.model;
    }

    BulkModulusParameters&
    operator=(const BulkModulusParameters& bulk) {
      if (this != &bulk) {
        filename = bulk.filename;
        model = bulk.model;
      }
      return *this;
    }

    void
    outputProblemSpec(Uintah::ProblemSpecP& ps) {
      ps->appendElement("filename", filename);
      ps->appendElement("mean_elastic_strain", model.meanStrain[0]);
      ps->appendElement("mean_plastic_strain", model.meanStrain[1]);
      ps->appendElement("std_dev_elastic_strain", model.stdStrain[0]);
      ps->appendElement("std_dev_plastic_strain", model.stdStrain[1]);
      ps->appendElement("mean_bulk_modulus", model.meanBulkModulus);
      ps->appendElement("std_dev_bulk_modulus", model.stdBulkModulus);
    }
  };

  /* Tangent shear modulus parameters */
  struct ShearModulusParameters
  {
    double G0;
    double nu;
  };

  BulkModulusParameters<double> d_bulk;
  ShearModulusParameters d_shear;

  void checkInputParameters();

  double computeBulkModulus(const double& elasticStrain,
                            const double& plasticStrain) const;
  double computeShearModulus(const double& bulkModulus) const;

};
} // End namespace Vaango

#endif // __ELASTIC_MODULI_NEURAL_NET_BULK_MODULUS_MODEL_H__
