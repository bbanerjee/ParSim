#include <iostream>
#include <string>

#include "H5Cpp.h"
#include <submodules/json/src/json.hpp>
#include <Eigen/Core>

#include <gtest/gtest.h>

using EigenMatrixF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct Layer {
  std::string name;
  std::string activation;
  int units;
  int input_size;
  EigenMatrixF weights;
  EigenMatrixF bias;
};

TEST(HDF5Tests, readTest)
{
  std::vector<Layer> layers;

  const H5std_string FILE_NAME("mlp_regression_keras_total_scaled.h5");
  const H5std_string DATASET_NAME("");
  try {

    H5::Exception::dontPrint();

    H5::H5File file(FILE_NAME, H5F_ACC_RDONLY);

    H5std_string test;
    auto group = file.openGroup("/");
    auto attribute = group.openAttribute("model_config");
    auto datatype = attribute.getDataType();
    attribute.read(datatype, test);

    //std::cout << test << std::endl;
    std::stringstream ss;
    ss.str(test);
    nlohmann::json doc;
    doc << ss;
    //std::cout << doc;

    if (doc["class_name"] != "Sequential") {
      std::cout << "Not sequential!\n";
    }

    auto config = doc["config"];
    //std::cout << config << "\n";
    //std::cout << config.size() << "\n";
    
    Layer prev_layer, current_layer;
    for (auto it = config.begin(); it != config.end(); ++it) {
      //std::cout << *it << "\n";
      //std::cout << (*it).size() << "\n";
      auto layer_class_name = (*it)["class_name"];
      //std::cout << "layer type = " << layer_class_name << "\n";
      auto layer_config = (*it)["config"];
      //std::cout << "layer_config = " << layer_config.size() << "\n";
      for (auto l_it = layer_config.begin(); l_it != layer_config.end(); ++l_it) {
        if (it == config.begin()) {
          std::cout << "\t" << l_it.key() << " = " << l_it.value() << "\n";
          if (l_it.key() == "batch_input_shape") {
            EXPECT_TRUE(l_it.value().is_array());
            EXPECT_TRUE(l_it.value().size() > 0);
            const std::size_t offset = l_it.value().front().is_null() ? 1 : 0;
            std::cout << "offset = " << offset << "\n";
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
      layers.push_back(current_layer);
      prev_layer = current_layer;
    }

    group = file.openGroup("/model_weights");
    attribute = group.openAttribute("layer_names");
    datatype = attribute.getDataType();
    auto dataspace = attribute.getSpace();
    auto size = datatype.getSize();
    std::cout << "size = " << datatype.getSize() << std::endl;
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims;
    int ndims = dataspace.getSimpleExtentDims(&dims, nullptr);
    std::cout << "rank = " << rank << " ndims = " << ndims << " " << dims << std::endl;
    char test1[dims][size];
    attribute.read(datatype, (void *)test1);
    //for (auto ii=0u; ii < dims; ++ii) {
    //  std::cout << std::string(test1[ii], size) << std::endl;
    //}

    for (auto ii = 0u; ii < dims; ++ii) {
      std::string layer_name(test1[ii], size);
      std::cout << "layer_name = " << layer_name << std::endl;
      group = file.openGroup("/model_weights/"+layer_name);
      attribute = group.openAttribute("weight_names");
      datatype = attribute.getDataType();
      dataspace = attribute.getSpace();
      auto size_wn = datatype.getSize();
      std::cout << "size = " << size_wn << std::endl;
      hsize_t dims_wn;
      dataspace.getSimpleExtentDims(&dims_wn, nullptr);
      std::cout << "dims = " << dims_wn << std::endl;
      char test2[dims_wn][size_wn];
      attribute.read(datatype, (void *)test2);
      for (auto jj=0u; jj < dims_wn; ++jj) {
        std::string weights_name(test2[jj], size_wn);
        std::cout << "weights_name = " << weights_name << std::endl;

        auto dataset = group.openDataSet(weights_name);
        datatype = dataset.getDataType();
        dataspace = dataset.getSpace();
        rank = dataspace.getSimpleExtentNdims();
        std::cout << "rank = " << rank << std::endl;
        hsize_t dims_data[rank];
        dataspace.getSimpleExtentDims(dims_data, nullptr);
        for (auto kk=0; kk < rank; ++kk) {
          std::cout << "dims = " << dims_data[kk] << std::endl;
        }
        
        if (rank == 2) {
          //float test3[dims_data[0]][dims_data[1]];
          //dataset.read((void *)test3, datatype);
          //std::cout << "test3 = " ;
          //for (auto kk=0u; kk < dims_data[0]; ++kk) {
          //  for (auto ll=0u; ll < dims_data[1]; ++ll) {
          //    std::cout << test3[kk][ll] << " ";
          //  }
          //}
          //std::cout << "\n";
          EigenMatrixF mat(dims_data[0], dims_data[1]);
          dataset.read((void *)mat.data(), datatype);
          std::cout << "mat = " ;
          for (auto kk=0u; kk < dims_data[0]; ++kk) {
            for (auto ll=0u; ll < dims_data[1]; ++ll) {
              std::cout << mat(kk, ll) << " ";
            }
          }
          std::cout << "\n";
          layers[ii].weights = mat.transpose();
        } else {
          //float test3[dims_data[0]];
          //dataset.read((void *)test3, datatype);
          //std::cout << "test3 = " ;
          //for (auto kk=0u; kk < dims_data[0]; ++kk) {
          //  std::cout << test3[kk] << " ";
          //}
          //std::cout << "\n";
          EigenMatrixF mat(dims_data[0], 1);
          dataset.read((void *)mat.data(), datatype);
          std::cout << "mat = " ;
          for (auto kk=0u; kk < dims_data[0]; ++kk) {
            std::cout << mat(kk, 0) << " ";
          }
          std::cout << "\n";
          layers[ii].bias = mat;
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


  float total_strain_min = 0.0;
  float total_strain_max = 0.452;
  float pressure_min = 0.0;
  float pressure_max = 1.0e6;
  
  EigenMatrixF input_orig(layers[0].input_size, 1);
  input_orig(0, 0) = -10;
  input_orig(1, 0) = 0;
  EigenMatrixF input = input_orig.unaryExpr([&total_strain_min, &total_strain_max](float x) -> float 
    {
      return (x - total_strain_min)/(total_strain_max - total_strain_min);
    });

  for (const auto& layer : layers) {
    std::cout << "name = " << layer.name << " activation = " << layer.activation
              << " input_size = " << layer.input_size << " units = " << layer.units << std::endl;
    //std::cout << "weights = " << layer.weights << "\n";
    //std::cout << "bias = " << layer.bias << "\n";

    EigenMatrixF output = layer.weights * input + layer.bias;
    std::cout << "output =" << output << "\n";
    if (layer.activation == "sigmoid") { 
      EigenMatrixF transformed = output.unaryExpr([](float x) -> float {
        float divisor = 1 + std::exp(-x);
        if (divisor == 0)
        {
            divisor = std::numeric_limits<float>::min();
        }
        return 1 / divisor;
      });
      input = transformed;
    } else if (layer.activation == "relu") {
      EigenMatrixF transformed = output.unaryExpr([](float x) -> float {
        return std::max<float>(x, 0);
      });
      input = transformed;
    } else if (layer.activation == "linear") {
      input = output;
    }
    std::cout << "input =" << input << "\n";
  }

  EigenMatrixF output = input.unaryExpr([&pressure_min, &pressure_max](float x) -> float 
    {
      return pressure_min + x * (pressure_max - pressure_min);
    });
  std::cout << "Prediction = " << output << std::endl;

}
