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

    ASSERT_EQ(doc["class_name"], "Sequential");

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
          //std::cout << "\t" << l_it.key() << " = " << l_it.value() << "\n";
          if (l_it.key() == "batch_input_shape") {
            EXPECT_TRUE(l_it.value().is_array());
            EXPECT_TRUE(l_it.value().size() > 0);
            const std::size_t offset = l_it.value().front().is_null() ? 1 : 0;
            //std::cout << "offset = " << offset << "\n";
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
    //std::cout << "size = " << datatype.getSize() << std::endl;
    ASSERT_EQ(datatype.getSize(), 7);

    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims;
    int ndims = dataspace.getSimpleExtentDims(&dims, nullptr);
    //std::cout << "rank = " << rank << " ndims = " << ndims << " " << dims << std::endl;
    ASSERT_EQ(rank, 1);
    ASSERT_EQ(ndims, 1);
    ASSERT_EQ(dims, 4);

    char test1[dims][size];
    attribute.read(datatype, (void *)test1);
    //for (auto ii=0u; ii < dims; ++ii) {
    //  std::cout << std::string(test1[ii], size) << std::endl;
    //}

    for (auto ii = 0u; ii < dims; ++ii) {
      std::string layer_name(test1[ii], size);
      //std::cout << "layer_name = " << layer_name << std::endl;
      switch(ii) {
        case 0:  ASSERT_EQ(layer_name, "dense_1"); break;
        case 1:  ASSERT_EQ(layer_name, "dense_2"); break;
        case 2:  ASSERT_EQ(layer_name, "dense_3"); break;
        case 3:  ASSERT_EQ(layer_name, "dense_4"); break;
      };

      group = file.openGroup("/model_weights/"+layer_name);
      attribute = group.openAttribute("weight_names");
      datatype = attribute.getDataType();
      dataspace = attribute.getSpace();
      auto size_wn = datatype.getSize();
      hsize_t dims_wn;
      dataspace.getSimpleExtentDims(&dims_wn, nullptr);

      //std::cout << "size = " << size_wn << std::endl;
      //std::cout << "dims = " << dims_wn << std::endl;
      ASSERT_EQ(size_wn, 16); 
      ASSERT_EQ(dims_wn, 2); 

      char test2[dims_wn][size_wn];
      attribute.read(datatype, (void *)test2);
      for (auto jj=0u; jj < dims_wn; ++jj) {
        std::string weights_name(test2[jj], size_wn);

        auto dataset = group.openDataSet(weights_name);
        datatype = dataset.getDataType();
        dataspace = dataset.getSpace();
        rank = dataspace.getSimpleExtentNdims();
        hsize_t dims_data[rank];
        dataspace.getSimpleExtentDims(dims_data, nullptr);

        //std::cout << "weights_name = " << weights_name << std::endl;
        //std::cout << "rank = " << rank << std::endl;
        //for (auto kk=0; kk < rank; ++kk) {
        //  std::cout << "dims = " << dims_data[kk] << std::endl;
        //}
        
        // Remove extra \0 characters in the bias name
        auto pos = weights_name.find_last_not_of('\0');
        if (pos == std::string::npos) weights_name.clear();
        else weights_name.resize(pos + 1);

        switch(ii) {
          case 0:  switch(jj) {
                     case 0: ASSERT_EQ(weights_name, "dense_1/kernel:0"); 
                             ASSERT_EQ(rank, 2);
                             ASSERT_EQ(dims_data[0], 2);
                             ASSERT_EQ(dims_data[1], 64);
                             break;
                     case 1: ASSERT_EQ(weights_name, "dense_1/bias:0"); 
                             ASSERT_EQ(rank, 1);
                             ASSERT_EQ(dims_data[0], 64);
                             break;
                   };
                   break;
          case 1:  switch(jj) {
                     case 0: ASSERT_EQ(weights_name, "dense_2/kernel:0"); 
                             ASSERT_EQ(rank, 2);
                             ASSERT_EQ(dims_data[0], 64);
                             ASSERT_EQ(dims_data[1], 32);
                             break;
                     case 1: ASSERT_EQ(weights_name, "dense_2/bias:0"); 
                             ASSERT_EQ(rank, 1);
                             ASSERT_EQ(dims_data[0], 32);
                             break;
                   };
                   break;
          case 2:  switch(jj) {
                     case 0: ASSERT_EQ(weights_name, "dense_3/kernel:0"); 
                             ASSERT_EQ(rank, 2);
                             ASSERT_EQ(dims_data[0], 32);
                             ASSERT_EQ(dims_data[1], 32);
                             break;
                     case 1: ASSERT_EQ(weights_name, "dense_3/bias:0"); 
                             ASSERT_EQ(rank, 1);
                             ASSERT_EQ(dims_data[0], 32);
                             break;
                   };
                   break;
          case 3:  switch(jj) {
                     case 0: ASSERT_EQ(weights_name, "dense_4/kernel:0"); 
                             ASSERT_EQ(rank, 2);
                             ASSERT_EQ(dims_data[0], 32);
                             ASSERT_EQ(dims_data[1], 1);
                             break;
                     case 1: ASSERT_EQ(weights_name, "dense_4/bias:0"); 
                             ASSERT_EQ(rank, 1);
                             ASSERT_EQ(dims_data[0], 1);
                             break;
                   };
                   break;
        };

        if (rank == 2) {
          EigenMatrixF weights(dims_data[0], dims_data[1]);
          dataset.read((void *)weights.data(), datatype);
          layers[ii].weights = weights.transpose();

          //std::cout << "weights = \n" ;
          //for (auto kk=0u; kk < dims_data[0]; ++kk) {
          //  for (auto ll=0u; ll < dims_data[1]; ++ll) {
          //    std::cout << "[" << kk << "," << ll << "]=" << weights(kk, ll) << " ";
          //  }
          //  std::cout << "\n";
          //}
          //std::cout << "\n";

          switch(ii) {
            case 0: ASSERT_NEAR(weights(0, 0), 2.93526, 1.0e-5);
                    ASSERT_NEAR(weights(1, 31), 2.58096, 1.0e-5);
                    ASSERT_NEAR(weights(1, 63), -0.357476, 1.0e-5);
                    break;
            case 1: ASSERT_NEAR(weights(0, 0), 0.915178, 1.0e-5);
                    ASSERT_NEAR(weights(31, 15), -1.12838, 1.0e-5);
                    ASSERT_NEAR(weights(63, 31), 1.08655, 1.0e-5);
                    break;
            case 2: EXPECT_NEAR(weights(0, 0), -8.95587, 1.0e-5);
                    EXPECT_NEAR(weights(15, 15), -10.352586, 1.0e-5);
                    ASSERT_NEAR(weights(31, 31), -0.056907, 1.0e-5);
                    break;
            case 3: EXPECT_NEAR(weights(0, 0), -4.41087, 1.0e-5);
                    EXPECT_NEAR(weights(15, 0), -7.43057, 1.0e-5);
                    ASSERT_NEAR(weights(31, 0), -0.01294, 1.0e-5);
                    break;
          };

        } else {
          EigenMatrixF bias(dims_data[0], 1);
          dataset.read((void *)bias.data(), datatype);
          layers[ii].bias = bias;

          //std::cout << "bias = \n" ;
          //for (auto kk=0u; kk < dims_data[0]; ++kk) {
          //  std::cout << "[" << kk << "]=" << bias(kk, 0) << " ";
          //}
          //std::cout << "\n\n";

          switch(ii) {
            case 0: ASSERT_NEAR(bias(0, 0), -3.04455, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -1.833757, 1.0e-5);
                    ASSERT_NEAR(bias(63, 0), -3.00222, 1.0e-5);
                    break;
            case 1: EXPECT_NEAR(bias(0, 0), -2.38528, 1.0e-5);
                    EXPECT_NEAR(bias(15, 0), -2.47038, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -2.41679, 1.0e-5);
                    break;
            case 2: EXPECT_NEAR(bias(0, 0), 2.69647, 1.0e-5);
                    EXPECT_NEAR(bias(15, 0), 3.03744, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -0.015684, 1.0e-5);
                    break;
            case 3: ASSERT_NEAR(bias(0, 0), 2.44017, 1.0e-5);
                    break;
          };
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
  //input_orig(0, 0) = -10;
  //input_orig(1, 0) = 0;
  input_orig(0, 0) = 1.0e-6;
  input_orig(1, 0) = 0;
  EigenMatrixF input = input_orig.unaryExpr([&total_strain_min, &total_strain_max](float x) -> float 
    {
      return (x - total_strain_min)/(total_strain_max - total_strain_min);
    });

  int layer_num = 0;
  for (const auto& layer : layers) {
    //std::cout << "name = " << layer.name << " activation = " << layer.activation
    //          << " input_size = " << layer.input_size << " units = " << layer.units << std::endl;
    //std::cout << "weights = " << layer.weights << "\n";
    //std::cout << "bias = " << layer.bias << "\n";

    switch(layer_num) {
      case 0:  EXPECT_EQ(layer.name, "dense_1");
               EXPECT_EQ(layer.activation, "sigmoid");
               EXPECT_EQ(layer.input_size, 2);
               ASSERT_EQ(layer.units, 64);
               EXPECT_NEAR(layer.weights(0, 0), 2.93526, 1.0e-5);
               EXPECT_NEAR(layer.weights(31, 1), 2.58096, 1.0e-5);
               ASSERT_NEAR(layer.weights(63, 1), -0.357476, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), -3.04455, 1.0e-5);
               EXPECT_NEAR(layer.bias(31, 0), -1.833757, 1.0e-5);
               ASSERT_NEAR(layer.bias(63, 0), -3.00222, 1.0e-5);
               break;
      case 1:  EXPECT_EQ(layer.name, "dense_2");
               EXPECT_EQ(layer.activation, "sigmoid");
               EXPECT_EQ(layer.input_size, 64);
               ASSERT_EQ(layer.units, 32);
               EXPECT_NEAR(layer.weights(0, 0), 0.915178, 1.0e-5);
               EXPECT_NEAR(layer.weights(15, 31), -1.12838, 1.0e-5);
               ASSERT_NEAR(layer.weights(31, 63), 1.08655, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), -2.38528, 1.0e-5);
               EXPECT_NEAR(layer.bias(15, 0), -2.47038, 1.0e-5);
               ASSERT_NEAR(layer.bias(31, 0), -2.41679, 1.0e-5);
               break;
      case 2:  EXPECT_EQ(layer.name, "dense_3");
               EXPECT_EQ(layer.activation, "relu");
               EXPECT_EQ(layer.input_size, 32);
               ASSERT_EQ(layer.units, 32);
               EXPECT_NEAR(layer.weights(0, 0), -8.95587, 1.0e-5);
               EXPECT_NEAR(layer.weights(15, 15), -10.352586, 1.0e-5);
               ASSERT_NEAR(layer.weights(31, 31), -0.056907, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), 2.69647, 1.0e-5);
               EXPECT_NEAR(layer.bias(15, 0), 3.03744, 1.0e-5);
               ASSERT_NEAR(layer.bias(31, 0), -0.015684, 1.0e-5);
               break;
      case 3:  EXPECT_EQ(layer.name, "dense_4");
               EXPECT_EQ(layer.activation, "linear");
               EXPECT_EQ(layer.input_size, 32);
               ASSERT_EQ(layer.units, 1);
               EXPECT_NEAR(layer.weights(0, 0), -4.41087, 1.0e-5);
               EXPECT_NEAR(layer.weights(0, 15), -7.43057, 1.0e-5);
               ASSERT_NEAR(layer.weights(0, 31), -0.01294, 1.0e-5);
               ASSERT_NEAR(layer.bias(0, 0), 2.44017, 1.0e-5);
               break;
    };

    EigenMatrixF output = layer.weights * input + layer.bias;

    //std::cout << "input =\n" << input << "\n";
    //std::cout << "output =\n" << output << "\n";

    switch(layer_num) {
      case 0: EXPECT_EQ(input.rows(), 2);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 2.21239e-6, 1.0e-8);
              ASSERT_EQ(input(1, 0), 0);
              EXPECT_EQ(output.rows(), 64);
              EXPECT_EQ(output.cols(), 1);
              EXPECT_NEAR(output(0, 0), -3.04454, 1.0e-5);
              EXPECT_NEAR(output(31, 0), -1.83375, 1.0e-5);
              ASSERT_NEAR(output(63, 0), -3.00222, 1.0e-5);
              break;
      case 1: EXPECT_EQ(input.rows(), 64);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 0.0454537, 1.0e-6);
              ASSERT_NEAR(input(63, 0), 0.0473257, 1.0e-6);
              EXPECT_EQ(output.rows(), 32);
              ASSERT_EQ(output.cols(), 1);
              EXPECT_NEAR(output(0, 0), -9.95303, 1.0e-5);
              ASSERT_NEAR(output(31, 0), -9.94701, 1.0e-5);
              break;
      case 2: EXPECT_EQ(input.rows(), 32);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 4.75808e-5, 1.0e-8);
              ASSERT_NEAR(input(31, 0), 4.78684e-5, 1.0e-8);
              EXPECT_EQ(output.rows(), 32);
              ASSERT_EQ(output.cols(), 1);
              EXPECT_NEAR(output(0, 0), 2.69267, 1.0e-5);
              ASSERT_NEAR(output(31, 0), -0.0171661, 1.0e-5);
              break;
      case 3: EXPECT_EQ(input.rows(), 32);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 2.69267, 1.0e-5);
              ASSERT_NEAR(input(30, 0), 4.42327, 1.0e-5);
              EXPECT_EQ(output.rows(), 1);
              ASSERT_EQ(output.cols(), 1);
              ASSERT_NEAR(output(0, 0), -1.45476, 1.0e-5);
              break;
    };

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
    //std::cout << " transformed input =\n" << input << "\n";
    ++layer_num;
  }

  EigenMatrixF output = input.unaryExpr([&pressure_min, &pressure_max](float x) -> float 
    {
      return pressure_min + x * (pressure_max - pressure_min);
    });
  //std::cout << "Prediction = " << output << std::endl;
  ASSERT_NEAR(output(0, 0), -1.45476e6, 1.0);

}
