/*
 * The MIT License
 *
 * Copyright (c) 2015-2018 Parresia Research Limited, New Zealand
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

#include <CCA/Components/MPM/ConstitutiveModel/Models/ElasticModuli_NeuralNet.h>
#include <CCA/Components/MPM/ConstitutiveModel/Models/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

using namespace Vaango;

// Construct a default elasticity model.
ElasticModuli_NeuralNet::ElasticModuli_NeuralNet(Uintah::ProblemSpecP& ps)
  : d_bulk(ps)
{
  ps->require("G0", d_shear.G0);
  ps->require("nu", d_shear.nu);

  checkInputParameters();
}

//--------------------------------------------------------------
// Check that the input parameters are reasonable
//--------------------------------------------------------------
void
ElasticModuli_NeuralNet::checkInputParameters()
{
  std::ostringstream warn;

  /* TODO : Add checks for neural network parameters */
  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_NeuralNet::ElasticModuli_NeuralNet(const ElasticModuli_NeuralNet* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

void
ElasticModuli_NeuralNet::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "neural_net");

  d_bulk.outputProblemSpec(elasticModuli_ps);

  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("nu", d_shear.nu);
}

// Compute the elastic moduli
ElasticModuli
ElasticModuli_NeuralNet::getInitialElasticModuli() const
{
  double K = computeBulkModulus(0, 0);
  double G = computeShearModulus(K);
  return ElasticModuli(K, G);
}

ElasticModuli
ElasticModuli_NeuralNet::getCurrentElasticModuli(const ModelStateBase* state_input)
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Make sure the quantities are positive in compression
  double ev_e_bar = -(state->elasticStrainTensor).Trace();
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  // Compute the elastic moduli
  //if (ev_p_bar != 0) {
  //  std::cout << "ev_e = " << ev_e_bar << " ev_p = " << ev_p_bar;
  //}
  double K = computeBulkModulus(ev_e_bar, ev_p_bar);
  double G = computeShearModulus(K);
  //if (ev_p_bar != 0) {
  //  std::cout << " K = " << K << " G = " << G << std::endl;
  //}

  return ElasticModuli(K, G);
}

double 
ElasticModuli_NeuralNet::computeBulkModulus(const double& eps_v_e,
                                            const double& eps_v_p) const
{
  double epsilon = 1.0e-6;
  double eps_v = eps_v_e + eps_v_p;

  
  double pressure_lo = d_bulk.d_model.predict(eps_v - epsilon, std::max(eps_v_p, 0.0));
  double pressure_hi = d_bulk.d_model.predict(eps_v + epsilon, std::max(eps_v_p, 0.0));

  double K = (pressure_hi - pressure_lo)/(2*epsilon);
  if (K < 1.0e-6) {
    std::cout << std::setprecision(16) << "ev_e = " << eps_v_e
              << " ev_p = " << eps_v_p
              << " ev- = " << eps_v - epsilon
              << " ev+ = " << eps_v + epsilon
              << " p_lo = " << pressure_lo << " p_hi = " << pressure_hi
              << " K = " << K << std::endl;
  }
  return K;
}

double 
ElasticModuli_NeuralNet::computeShearModulus(const double& K) const
{
  double nu = d_shear.nu;
  double G = (nu > -1.0 && nu < 0.5) 
             ? 1.5*K*(1.0 - 2.0*nu)/(1.0 + nu) 
             : d_shear.G0;
  return G;
}

/* Get the elastic moduli and their derivatives with respect to a single
   plastic internal variable */
std::pair<ElasticModuli, ElasticModuli>
ElasticModuli_NeuralNet::getElasticModuliAndDerivatives(
  const ModelStateBase* state_input) const
{
  const ModelState_Tabular* state =
    dynamic_cast<const ModelState_Tabular*>(state_input);
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }

  // Make sure the quantities are positive in compression
  double ev_e_bar = -(state->elasticStrainTensor).Trace();
  double ev_p_bar = -(state->plasticStrainTensor).Trace();

  // Compute the elastic moduli
  double epsilon = 1.0e-8;
  double K = computeBulkModulus(ev_e_bar, ev_p_bar);
  double K_min = computeBulkModulus(ev_e_bar, ev_p_bar-epsilon);
  double K_max = computeBulkModulus(ev_e_bar, ev_p_bar+epsilon);
  double G = computeShearModulus(K);
  double G_min = computeShearModulus(K_min);
  double G_max = computeShearModulus(K_max);

  // Compute derivatives
  double dK_deps_p = -(K_max - K_min)/(2*epsilon);
  double dG_deps_p = -(G_max - G_min)/(2*epsilon);

  return std::make_pair(ElasticModuli(K, G),
                        ElasticModuli(dK_deps_p, dG_deps_p));
}

template<typename T>
void 
ElasticModuli_NeuralNet::NeuralNetworkModel<T>::readNeuralNetworkHDF5(const std::string& filename)
{
  using EigenMatrixRowMajor = 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  const H5std_string FILE_NAME(filename);
  const H5std_string DATASET_NAME("");
  try {

    H5::Exception::dontPrint();

    H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

    // Read the configuration from the HDF5 file 
    H5std_string model_config_json;
    auto group = file.openGroup("/");
    auto attribute = group.openAttribute("model_config");
    auto datatype = attribute.getDataType();
    attribute.read(datatype, model_config_json);

    // Parse the JSON configuration and store layer information
    std::stringstream ss;
    ss.str(model_config_json);
    nlohmann::json doc;
    doc << ss;

    if (doc["class_name"] != "Sequential") {
      std::ostringstream out;
      out << "**ERROR** The input Keras neural network model is required to be Sequential";
      throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
    }

    NeuralNetworkLayer<T> prev_layer, current_layer;
    auto config = doc["config"];
    for (auto it = config.begin(); it != config.end(); ++it) {
      auto layer_class_name = (*it)["class_name"];
      auto layer_config = (*it)["config"];
      for (auto l_it = layer_config.begin(); l_it != layer_config.end(); ++l_it) {
        if (it == config.begin()) {
          if (l_it.key() == "batch_input_shape") {
            const std::size_t offset = l_it.value().front().is_null() ? 1 : 0;
            current_layer.input_size = l_it.value()[0 + offset];
          }
        } else {
          current_layer.input_size = prev_layer.units;
        }
        if (l_it.key() == "name") {
          current_layer.name = l_it.value();
        }
        if (l_it.key() == "activation") {
          current_layer.activation = l_it.value();
        }
        if (l_it.key() == "units") {
          current_layer.units = l_it.value();
        }
      }
      d_layers.push_back(current_layer);
      prev_layer = current_layer;
    }

    // Reads the weights and biases (HDF5)
    group = file.openGroup("/model_weights");
    attribute = group.openAttribute("layer_names");
    datatype = attribute.getDataType();
    auto dataspace = attribute.getSpace();
    auto size = datatype.getSize();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims;
    dataspace.getSimpleExtentDims(&dims, nullptr);

    char layer_names[dims][size];
    attribute.read(datatype, (void *)layer_names);

    for (auto ii = 0u; ii < dims; ++ii) {
      std::string layer_name(layer_names[ii], size);
      group = file.openGroup("/model_weights/"+layer_name);
      attribute = group.openAttribute("weight_names");
      datatype = attribute.getDataType();
      dataspace = attribute.getSpace();
      auto size_weights_bias = datatype.getSize();
      hsize_t dims_weights_bias;
      dataspace.getSimpleExtentDims(&dims_weights_bias, nullptr);

      char data_labels[dims_weights_bias][size_weights_bias];
      attribute.read(datatype, (void *)data_labels);

      for (auto jj=0u; jj < dims_weights_bias; ++jj) {
        std::string weights_name(data_labels[jj], size_weights_bias);
        auto dataset = group.openDataSet(weights_name);
        datatype = dataset.getDataType();
        dataspace = dataset.getSpace();
        rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_data[rank];
        dataspace.getSimpleExtentDims(dims_data, nullptr);
        
        if (rank == 2) {
          EigenMatrixRowMajor mat(dims_data[0], dims_data[1]);
          dataset.read((void *)mat.data(), datatype);
          d_layers[ii].weights = mat.transpose();
          //std::cout << "weights = \n" << mat << std::endl;
        } else {
          EigenMatrixRowMajor mat(dims_data[0], 1);
          dataset.read((void *)mat.data(), datatype);
          d_layers[ii].bias = mat;
          //std::cout << "biases = \n" << mat << std::endl;
        }
      }
    }
    
  } catch (H5::FileIException error) {
    error.printError();
  } catch (H5::DataSetIException error) {
    error.printError();
  } catch (H5::DataSpaceIException error) {
    error.printError();
  } catch (H5::DataTypeIException error) {
    error.printError();
  } catch (H5::AttributeIException error) {
    error.printError();
  }
}

template<typename T>
double 
ElasticModuli_NeuralNet::NeuralNetworkModel<T>::predict(double totalVolStrain, double plasticVolStrain) const
{
  using EigenMatrixRowMajor = 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  T minStrain = static_cast<T>(d_minStrain);
  T maxStrain = static_cast<T>(d_maxStrain);
  T minPressure = static_cast<T>(d_minPressure);
  T maxPressure = static_cast<T>(d_maxPressure);

  EigenMatrixRowMajor input(2, 1);
  input(0, 0) = static_cast<T>(totalVolStrain);
  input(1, 0) = static_cast<T>(plasticVolStrain);
  
  input = input.unaryExpr([&minStrain, &maxStrain](T x) -> T 
    {
      return (x - minStrain)/(maxStrain - minStrain);
    });

  for (const auto& layer : d_layers) {
    EigenMatrixRowMajor output = layer.weights * input + layer.bias;
    if (layer.activation == "sigmoid") { 
      EigenMatrixRowMajor layer_output = output.unaryExpr([](T x) -> T 
        {
          T divisor = 1 + std::exp(-x);
          if (divisor == 0) {
            divisor = std::numeric_limits<T>::min();
          }
          return 1 / divisor;
        });
      input = layer_output;
    } else if (layer.activation == "relu") {
      EigenMatrixRowMajor layer_output = output.unaryExpr([](T x) -> T 
        {
          return std::max<T>(x, 0);
        });
      input = layer_output;
    } else if (layer.activation == "linear") {
      EigenMatrixRowMajor layer_output = output;
      input = layer_output;
    }
  }

  input = input.unaryExpr([&minPressure, &maxPressure](T x) -> T 
    {
      return minPressure + x * (maxPressure - minPressure);
    });
  return input(0, 0);
}
