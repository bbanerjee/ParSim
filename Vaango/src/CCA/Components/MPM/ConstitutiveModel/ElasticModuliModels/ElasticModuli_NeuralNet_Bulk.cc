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

#include <CCA/Components/MPM/ConstitutiveModel/ElasticModuliModels/ElasticModuli_NeuralNet_Bulk.h>
#include <CCA/Components/MPM/ConstitutiveModel/ModelState/ModelState_Tabular.h>
#include <Core/Exceptions/InternalError.h>
#include <Core/Exceptions/InvalidValue.h>
#include <Core/Exceptions/ProblemSetupException.h>

#define TENSORFLOW_2_0
//#define DEBUG_HDF5_READ

using namespace Vaango;

// Construct a default elasticity model.
ElasticModuli_NeuralNet_Bulk::ElasticModuli_NeuralNet_Bulk(Uintah::ProblemSpecP& ps)
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
ElasticModuli_NeuralNet_Bulk::checkInputParameters()
{
  std::ostringstream warn;

  if (d_shear.G0 <= 0.0) {
    warn << "G0 must be positive. G0 = " << d_shear.G0 << std::endl;
    throw Uintah::ProblemSetupException(warn.str(), __FILE__, __LINE__);
  }
}

// Construct a copy of a elasticity model.
ElasticModuli_NeuralNet_Bulk::ElasticModuli_NeuralNet_Bulk(const ElasticModuli_NeuralNet_Bulk* model)
{
  d_bulk = model->d_bulk;
  d_shear = model->d_shear;
}

void
ElasticModuli_NeuralNet_Bulk::outputProblemSpec(Uintah::ProblemSpecP& ps)
{
  Uintah::ProblemSpecP elasticModuli_ps = ps->appendChild("elastic_moduli_model");
  elasticModuli_ps->setAttribute("type", "neural_net_bulk");

  d_bulk.outputProblemSpec(elasticModuli_ps);

  elasticModuli_ps->appendElement("G0", d_shear.G0);
  elasticModuli_ps->appendElement("nu", d_shear.nu);
}

// Compute the elastic moduli
ElasticModuli
ElasticModuli_NeuralNet_Bulk::getInitialElasticModuli() const
{
  double K = computeBulkModulus(0, 0);
  double G = computeShearModulus(K);
  return ElasticModuli(K, G);
}

ElasticModuli
ElasticModuli_NeuralNet_Bulk::getCurrentElasticModuli(const ModelStateBase* state_input) const
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
ElasticModuli_NeuralNet_Bulk::computeBulkModulus(const double& eps_v_e,
                                                 const double& eps_v_p) const
{
  double K = d_bulk.model.predict(std::max(eps_v_e, 0.0), std::max(eps_v_p, 0.0));
  if (K < 1.0e-6) {
    std::cout << std::setprecision(16) << "ev_e = " << eps_v_e
              << " ev_p = " << eps_v_p
              << " K = " << K << std::endl;
  }
  return K;
}

double 
ElasticModuli_NeuralNet_Bulk::computeShearModulus(const double& K) const
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
ElasticModuli_NeuralNet_Bulk::getElasticModuliAndDerivatives(
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
ElasticModuli_NeuralNet_Bulk::NeuralNetworkModel<T>::readNeuralNetworkHDF5(const std::string& filename)
{
  //std::cout << "Reading hdf5\n";
  using EigenMatrixRowMajor = 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  const H5std_string FILE_NAME(filename);
  const H5std_string DATASET_NAME("");
  try {

    H5::Exception::dontPrint();

    H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

    #ifdef DEBUG_HDF5_READ
    std::cout << "Reading configuration from hdf5\n";
    #endif

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

    #ifdef DEBUG_HDF5_READ
    std::cout << doc.dump() << "\n";
    #endif

    #ifdef DEBUG_HDF5_READ
    std::cout << "Storing layer info from hdf5\n";
    #endif
    if (doc["class_name"] != "Sequential") {
      std::ostringstream out;
      out << "**ERROR** The input Keras neural network model is required to be Sequential";
      throw Uintah::InternalError(out.str(), __FILE__, __LINE__);
    }
    #ifdef DEBUG_HDF5_READ
    std::cout << "Read class name info from hdf5\n";
    #endif

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
          current_layer.name = l_it.value();
          //std::cout << "layer name from hdf5 " << l_it.value() << "\n";
        }
        if (l_it.key() == "activation") {
          current_layer.activation = l_it.value();
          //std::cout << "layer activation from hdf5 " << l_it.value() << "\n";
        }
        if (l_it.key() == "units") {
          current_layer.units = l_it.value();
          //std::cout << "layer units from hdf5 " << l_it.value() << "\n";
        }
      }
      layers.push_back(current_layer);
      prev_layer = current_layer;
    }

    #ifdef DEBUG_HDF5_READ
    std::cout << "Read weights and biases from hdf5 \n";
    #endif
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

    #ifdef DEBUG_HDF5_READ
    std::cout << "Dims from hdf5 = " << dims << "\n";
    #endif
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
      #ifdef DEBUG_HDF5_READ
      std::cout << "Read attribute from hdf5 = " << ii << "\n";
      #endif

      #ifdef DEBUG_HDF5_READ
      std::cout << "Dims weights/bias from hdf5 = " << dims_weights_bias << "\n";
      #endif
      for (auto jj=0u; jj < dims_weights_bias; ++jj) {
        std::string weights_name(data_labels[jj], size_weights_bias);
        auto dataset = group.openDataSet(weights_name);
        datatype = dataset.getDataType();
        dataspace = dataset.getSpace();
        rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_data[rank];
        dataspace.getSimpleExtentDims(dims_data, nullptr);
        
        #ifdef DEBUG_HDF5_READ
        std::cout << "Read dims/weights from hdf5 = " << jj << "\n";
        #endif
        if (rank == 2) {
          EigenMatrixRowMajor mat(dims_data[0], dims_data[1]);
          dataset.read((void *)mat.data(), datatype);
          layers[ii].weights = mat.transpose();
          #ifdef DEBUG_HDF5_READ
          std::cout << "weights = \n" << mat << std::endl;
          #endif
        } else {
          EigenMatrixRowMajor mat(dims_data[0], 1);
          dataset.read((void *)mat.data(), datatype);
          layers[ii].bias = mat;
          #ifdef DEBUG_HDF5_READ
          std::cout << "biases = \n" << mat << std::endl;
          #endif
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

  #ifdef DEBUG_HDF5_READ
  std::cout << "Done reading hdf5\n";
  #endif
}

template<typename T>
double 
ElasticModuli_NeuralNet_Bulk::NeuralNetworkModel<T>::predict(double elasticVolStrain, double plasticVolStrain) const
{
  using EigenMatrixRowMajor = 
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  // Scale input data
  double eps_e_v_translate = elasticVolStrain - meanStrain[0];
  double eps_p_v_translate = plasticVolStrain - meanStrain[1];
  double eps_e_v = (stdStrain[0] == 0.0) ? eps_e_v_translate : eps_e_v_translate/stdStrain[0];
  double eps_p_v = (stdStrain[1] == 0.0) ? eps_p_v_translate : eps_p_v_translate/stdStrain[1];

  T K_mean = static_cast<T>(meanBulkModulus);
  T K_std = static_cast<T>(stdBulkModulus);

  EigenMatrixRowMajor input(2, 1);
  input(0, 0) = static_cast<T>(eps_e_v);
  input(1, 0) = static_cast<T>(eps_p_v);
  
  for (const auto& layer : layers) {
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

  input = input.unaryExpr([&K_mean, &K_std](T x) -> T 
    {
      return K_mean + x * K_std;
    });
  return input(0, 0);
}

/*! Compute derivatives of moduli with respect to internal variables */
std::vector<ElasticModuli> 
ElasticModuli_NeuralNet_Bulk::computeDModuliDIntVar(const ModelStateBase* state) const
{
  std::pair<ElasticModuli, ElasticModuli> K_dK = getElasticModuliAndDerivatives(state);
  std::vector<ElasticModuli> derivs;
  derivs.push_back(K_dK.second);
  return derivs;
}

/*! Compute moduli and derivatives of moduli with respect to internal variables */
std::pair<ElasticModuli, std::vector<ElasticModuli>>
ElasticModuli_NeuralNet_Bulk::computeModuliAndDModuliDIntVar(const ModelStateBase* state) const
{
  auto K_dK = getElasticModuliAndDerivatives(state);
  std::vector<ElasticModuli> derivs;
  derivs.push_back(K_dK.second);
  return std::make_pair(K_dK.first, derivs);
}
