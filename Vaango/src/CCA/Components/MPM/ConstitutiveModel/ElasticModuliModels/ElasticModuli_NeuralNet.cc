/*
 * The MIT License
 *
 * Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_NeuralNet.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

#include <format>
#include <span>
#include <ranges>

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
ElasticModuli_NeuralNet::getCurrentElasticModuli(const ModelStateBase* state_input) const
{
  const ModelState_Tabular* state =
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

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
    static_cast<const ModelState_Tabular*>(state_input);
  /*
  if (!state) {
    std::ostringstream out;
    out << "**ERROR** The correct ModelState object has not been passed."
        << " Need ModelState_Tabular.";
    throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
  }
  */

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
  //std::cout << "Reading hdf5\n";
  using EigenMatrixRowMajor = 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  const H5std_string FILE_NAME(filename);
  const H5std_string DATASET_NAME("");
  try {

    H5::Exception::dontPrint();

    H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

    //std::cout << "Reading configuration from hdf5\n";
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
    ss >> doc;

    //std::cout << doc.dump() << "\n";

    //std::cout << "Storing layer info from hdf5\n";
    if (doc["class_name"] != "Sequential") {
      std::ostringstream out;
      out << "**ERROR** The input Keras neural network model is required to be Sequential";
      throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
    }
    //std::cout << "Read class name info from hdf5\n";

    NeuralNetworkLayer<T> prev_layer, current_layer;
#ifndef TENSORFLOW_2_0
    auto config = doc["config"];
#else
    auto config = doc["config"]["layers"];
#endif
    //std::cout << config.dump(2) << "\n";
    //std::cout << "Read config info from hdf5\n";
    for (auto it = config.begin(); it != config.end(); ++it) {
      auto layer_class_name = (*it)["class_name"];
      //std::cout << "Read config class name info from hdf5\n";
      auto layer_config = (*it)["config"];
      //std::cout << "Config info from hdf5 " << layer_class_name << "\n";
      for (auto l_it = layer_config.begin(); l_it != layer_config.end(); ++l_it) {
        //std::cout << "Config layer info from hdf5 " << l_it.value() << "\n";
        if (it == config.begin()) {
          if (l_it.key() == "batch_input_shape") {
            const std::size_t offset = l_it.value().front().is_null() ? 1 : 0;
            current_layer.input_size = l_it.value()[0 + offset];
          }
        } else {
          current_layer.input_size = prev_layer.units;
        }
        if (l_it.key() == "name") {
          current_layer.name = l_it.value().get<std::string>();
          //std::cout << "layer name from hdf5 " << l_it.value() << "\n";
        }
        if (l_it.key() == "activation") {
          current_layer.activation = l_it.value().get<std::string>();
          //std::cout << "layer activation from hdf5 " << l_it.value() << "\n";
        }
        if (l_it.key() == "units") {
          current_layer.units = l_it.value().get<int>();
          //std::cout << "layer units from hdf5 " << l_it.value() << "\n";
        }
      }
      d_layers.push_back(current_layer);
      prev_layer = current_layer;
    }

    //std::cout << "Read weights and biases from hdf5 \n";
    // Reads the weights and biases (HDF5)
    auto mw_group = file.openGroup("/model_weights");
    auto mw_attribute = mw_group.openAttribute("layer_names");
    auto mw_datatype = mw_attribute.getDataType();
    auto mw_dataspace = mw_attribute.getSpace();
    auto size = mw_datatype.getSize();
    int rank = mw_dataspace.getSimpleExtentNdims();
    hsize_t dims;
    mw_dataspace.getSimpleExtentDims(&dims, nullptr);

    // Create a 2D buffer
    std::vector<std::vector<char>> layer_names_buffer(dims);
    for (auto& name_buffer : layer_names_buffer) {
      name_buffer.resize(size);
    }

    // Create a flat buffer for HDF5 reading
    std::vector<char> flat_layers_buffer(dims * size);
    mw_attribute.read(mw_datatype, flat_layers_buffer.data());

    // Convert flat buffer to individual strings
    std::vector<std::string> layer_names;
    layer_names.reserve(dims);

    for (const auto i : std::views::iota(0uz, dims)) {
      const auto start_pos = i * size;
      std::span<const char> name_span{ flat_layers_buffer.data() + start_pos,
                                       size };

      // Find the actual string length (stop at null terminator)
      const auto actual_length =
        std::ranges::find(name_span, '\0') - name_span.begin();
      layer_names.emplace_back(name_span.data(), actual_length);
    }

    //std::cout << "Dims from hdf5 = " << dims << "\n";
    for (const auto& [ii, layer_name]: std::views::enumerate(layer_names)) {
      auto ln_group = file.openGroup(std::format("/model_weights/{}", layer_name));
      auto ln_attribute = ln_group.openAttribute("weight_names");
      auto ln_datatype = ln_attribute.getDataType();
      auto ln_dataspace = ln_attribute.getSpace();
      auto size_weights_bias = ln_datatype.getSize();
      hsize_t dims_weights_bias;
      ln_dataspace.getSimpleExtentDims(&dims_weights_bias, nullptr);

      std::vector<std::vector<char>> weight_names_buffer(dims_weights_bias);
      for (auto& name_buffer : weight_names_buffer) {
        name_buffer.resize(size_weights_bias);
      }

      // Create a flat buffer for HDF5 reading
      std::vector<char> flat_buffer(dims_weights_bias * size_weights_bias);
      ln_attribute.read(ln_datatype, flat_buffer.data());

      // Convert flat buffer to individual strings
      std::vector<std::string> weight_names;
      weight_names.reserve(dims_weights_bias);

      for (const auto i : std::views::iota(0uz, dims_weights_bias)) {
        const auto start_pos = i * size_weights_bias;
        std::span<const char> name_span{ flat_buffer.data() + start_pos,
                                         size_weights_bias };

        // Find the actual string length (stop at null terminator)
        const auto actual_length =
          std::ranges::find(name_span, '\0') - name_span.begin();
        weight_names.emplace_back(name_span.data(), actual_length);
      }

      //std::cout << "Read attribute from hdf5 = " << ii << "\n";
      //std::cout << "Dims weights/bias from hdf5 = " << dims_weights_bias << "\n";

      for (const auto& [jj, weights_name]: std::views::enumerate(weight_names)) {

        auto wn_dataset = ln_group.openDataSet(weights_name);
        auto wn_datatype = wn_dataset.getDataType();
        auto wn_dataspace = wn_dataset.getSpace();
        rank = wn_dataspace.getSimpleExtentNdims();
        std::vector<hsize_t> dims_data(rank);
        wn_dataspace.getSimpleExtentDims(dims_data.data(), nullptr);
        
        //std::cout << "Read dims/weights from hdf5 = " << jj << "\n";
        if (rank == 2) {
          EigenMatrixRowMajor mat(dims_data[0], dims_data[1]);
          wn_dataset.read((void *)mat.data(), wn_datatype);
          d_layers[ii].weights = mat.transpose();
          //std::cout << "weights = \n" << mat << std::endl;
        } else {
          EigenMatrixRowMajor mat(dims_data[0], 1);
          wn_dataset.read((void *)mat.data(), wn_datatype);
          d_layers[ii].bias = mat;
          //std::cout << "biases = \n" << mat << std::endl;
        }
      }
    }
    
  } catch (H5::FileIException const& error) {
    std::cout << "File Input Exception Reading hdf5\n";
    error.printErrorStack();
  } catch (H5::DataSetIException const& error) {
    std::cout << "Data set input Exception Reading hdf5\n";
    error.printErrorStack();
  } catch (H5::DataSpaceIException const& error) {
    std::cout << "Data space input Exception Reading hdf5\n";
    error.printErrorStack();
  } catch (H5::DataTypeIException const& error) {
    std::cout << "Data type input Exception Reading hdf5\n";
    error.printErrorStack();
  } catch (H5::AttributeIException const& error) {
    std::cout << "Attribute input Exception Reading hdf5\n";
    error.printErrorStack();
  }

  //std::cout << "Done reading hdf5\n";
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

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_NeuralNet::computeDModuliDIntVar([[maybe_unused]] const ModelStateBase* state) const
{
  std::vector<ElasticModuli> derivs;
  return derivs;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_NeuralNet::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  ElasticModuli moduli = getCurrentElasticModuli(state);
  std::vector<ElasticModuli> derivs;
  return std::make_pair(moduli, derivs);
}
