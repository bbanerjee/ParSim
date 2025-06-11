#include <iostream>
#include <string>
#include <format>
#include <ranges>

#include "H5Cpp.h"
#include <submodules/json/single_include/nlohmann/json.hpp>
#include <Eigen/Core>

#include <gtest/gtest.h>

using EigenMatrixF = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

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
    ss >> doc;
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
          current_layer.name = l_it.value().get<std::string>();
        }
        if (l_it.key() == "activation") {
          current_layer.activation = l_it.value().get<std::string>();
        }
        if (l_it.key() == "units") {
          current_layer.units = l_it.value().get<int>();
        }
      }
      layers.push_back(current_layer);
      prev_layer = current_layer;
    }

    auto layer_group = file.openGroup("/model_weights");
    auto layer_attribute = layer_group.openAttribute("layer_names");
    auto layer_datatype = layer_attribute.getDataType();
    auto layer_dataspace = layer_attribute.getSpace();
    auto size = layer_datatype.getSize();
    //std::cout << "size = " << datatype.getSize() << std::endl;
    ASSERT_EQ(layer_datatype.getSize(), 7);

    int rank = layer_dataspace.getSimpleExtentNdims();
    hsize_t dims;
    int ndims = layer_dataspace.getSimpleExtentDims(&dims, nullptr);
    //std::cout << "rank = " << rank << " ndims = " << ndims << " " << dims << std::endl;
    ASSERT_EQ(rank, 1);
    ASSERT_EQ(ndims, 1);
    ASSERT_EQ(dims, 4);


    // Create a 2D buffer
    std::vector<std::vector<char>> layer_names_buffer(dims);
    for (auto& name_buffer : layer_names_buffer) {
      name_buffer.resize(size);
    }

    // Create a flat buffer for HDF5 reading
    std::vector<char> flat_layers_buffer(dims * size);
    layer_attribute.read(layer_datatype, flat_layers_buffer.data());

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

    //for (auto layer_name: layer_names) {
    //  std::cout << layer_name << std::endl;
    //}

    for (const auto& [ii, layer_name]: std::views::enumerate(layer_names)) {
      //std::cout << "layer_name = " << layer_name << std::endl;
      switch(ii) {
        case 0:  ASSERT_EQ(layer_name, "dense_1"); break;
        case 1:  ASSERT_EQ(layer_name, "dense_2"); break;
        case 2:  ASSERT_EQ(layer_name, "dense_3"); break;
        case 3:  ASSERT_EQ(layer_name, "dense_4"); break;
      };

      auto weight_group = file.openGroup(std::format("/model_weights/{}", layer_name));
      auto weight_attribute = weight_group.openAttribute("weight_names");
      auto weight_datatype = weight_attribute.getDataType();
      auto weight_dataspace = weight_attribute.getSpace();
      auto size_wn = weight_datatype.getSize();
      hsize_t dims_wn;
      weight_dataspace.getSimpleExtentDims(&dims_wn, nullptr);

      //std::cout << "size = " << size_wn << std::endl;
      //std::cout << "dims = " << dims_wn << std::endl;
      ASSERT_EQ(size_wn, 16); 
      ASSERT_EQ(dims_wn, 2);

      std::vector<std::vector<char>> weight_names_buffer(dims_wn);
      for (auto& name_buffer : weight_names_buffer) {
        name_buffer.resize(size_wn);
      }

      // Create a flat buffer for HDF5 reading
      std::vector<char> flat_buffer(dims_wn * size_wn);
      weight_attribute.read(weight_datatype, flat_buffer.data());

      // Convert flat buffer to individual strings
      std::vector<std::string> weight_names;
      weight_names.reserve(dims_wn);

      for (const auto i : std::views::iota(0uz, dims_wn)) {
        const auto start_pos = i * size_wn;
        std::span<const char> name_span{ flat_buffer.data() + start_pos,
                                         size_wn };

        // Find the actual string length (stop at null terminator)
        const auto actual_length =
          std::ranges::find(name_span, '\0') - name_span.begin();
        weight_names.emplace_back(name_span.data(), actual_length);
      }

      for (const auto& [jj, weights_name]: std::views::enumerate(weight_names)) {

        auto wn_dataset = weight_group.openDataSet(weights_name);
        auto wn_datatype = wn_dataset.getDataType();
        auto wn_dataspace = wn_dataset.getSpace();
        rank = wn_dataspace.getSimpleExtentNdims();
        std::vector<hsize_t> dims_data(rank);
        wn_dataspace.getSimpleExtentDims(dims_data.data(), nullptr);

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
          wn_dataset.read((void *)weights.data(), wn_datatype);
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
            case 0: ASSERT_NEAR(weights(0, 0), 2.66871, 1.0e-5);
                    ASSERT_NEAR(weights(1, 31), 2.26223, 1.0e-5);
                    ASSERT_NEAR(weights(1, 63), -0.207627, 1.0e-5);
                    break;
            case 1: ASSERT_NEAR(weights(0, 0), 1.24319, 1.0e-5);
                    ASSERT_NEAR(weights(31, 15), -1.25566, 1.0e-5);
                    ASSERT_NEAR(weights(63, 31), 0.718625, 1.0e-5);
                    break;
            case 2: EXPECT_NEAR(weights(0, 0), 5.0301, 1.0e-5);
                    EXPECT_NEAR(weights(15, 15), 8.63985, 1.0e-5);
                    ASSERT_NEAR(weights(31, 31), 0.0212148, 1.0e-5);
                    break;
            case 3: EXPECT_NEAR(weights(0, 0), 5.0117, 1.0e-5);
                    EXPECT_NEAR(weights(15, 0), -29.7871, 1.0e-4);
                    ASSERT_NEAR(weights(31, 0), 0.0211935, 1.0e-5);
                    break;
          };

        } else {
          EigenMatrixF bias(dims_data[0], 1);
          wn_dataset.read((void *)bias.data(), wn_datatype);
          layers[ii].bias = bias;

          //std::cout << "bias = \n" ;
          //for (auto kk=0u; kk < dims_data[0]; ++kk) {
          //  std::cout << "[" << kk << "]=" << bias(kk, 0) << " ";
          //}
          //std::cout << "\n\n";

          switch(ii) {
            case 0: ASSERT_NEAR(bias(0, 0), -2.86732, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -0.493266, 1.0e-5);
                    ASSERT_NEAR(bias(63, 0), -2.86736, 1.0e-5);
                    break;
            case 1: EXPECT_NEAR(bias(0, 0), -2.97205, 1.0e-5);
                    EXPECT_NEAR(bias(15, 0), -3.25083, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -2.20559, 1.0e-5);
                    break;
            case 2: EXPECT_NEAR(bias(0, 0), 4.68216, 1.0e-5);
                    EXPECT_NEAR(bias(15, 0), 0.633401, 1.0e-5);
                    ASSERT_NEAR(bias(31, 0), -0.0314119, 1.0e-5);
                    break;
            case 3: ASSERT_NEAR(bias(0, 0), 3.10022, 1.0e-5);
                    break;
          };
        }
      }

    }
    
  } catch (H5::FileIException& error) {
    error.printErrorStack();
  } catch (H5::DataSetIException& error) {
    error.printErrorStack();
  } catch (H5::DataSpaceIException& error) {
    error.printErrorStack();
  } catch (H5::DataTypeIException& error) {
    error.printErrorStack();
  } catch (H5::AttributeIException& error) {
    error.printErrorStack();
  }


  double total_strain_min = 0.0;
  double total_strain_max = 0.452;
  double pressure_min = 0.0;
  double pressure_max = 1.0e6;
  
  EigenMatrixF input_orig(layers[0].input_size, 1);
  //input_orig(0, 0) = -10;
  //input_orig(1, 0) = 0;
  input_orig(0, 0) = 1.0e-6;
  input_orig(1, 0) = 0;
  EigenMatrixF input = input_orig.unaryExpr([&total_strain_min, &total_strain_max](double x) -> double 
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
               EXPECT_NEAR(layer.weights(0, 0), 2.66871, 1.0e-5);
               EXPECT_NEAR(layer.weights(31, 1), 2.26223, 1.0e-5);
               ASSERT_NEAR(layer.weights(63, 1), -0.207627, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), -2.86732, 1.0e-5);
               EXPECT_NEAR(layer.bias(31, 0), -0.493266, 1.0e-5);
               ASSERT_NEAR(layer.bias(63, 0), -2.86736, 1.0e-5);
               break;
      case 1:  EXPECT_EQ(layer.name, "dense_2");
               EXPECT_EQ(layer.activation, "sigmoid");
               EXPECT_EQ(layer.input_size, 64);
               ASSERT_EQ(layer.units, 32);
               EXPECT_NEAR(layer.weights(0, 0), 1.24319, 1.0e-5);
               EXPECT_NEAR(layer.weights(15, 31), -1.25566, 1.0e-5);
               ASSERT_NEAR(layer.weights(31, 63), 0.718625, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), -2.97205, 1.0e-5);
               EXPECT_NEAR(layer.bias(15, 0), -3.25083, 1.0e-5);
               ASSERT_NEAR(layer.bias(31, 0), -2.20559, 1.0e-5);
               break;
      case 2:  EXPECT_EQ(layer.name, "dense_3");
               EXPECT_EQ(layer.activation, "relu");
               EXPECT_EQ(layer.input_size, 32);
               ASSERT_EQ(layer.units, 32);
               EXPECT_NEAR(layer.weights(0, 0), 5.0301, 1.0e-5);
               EXPECT_NEAR(layer.weights(15, 15), 8.63985, 1.0e-5);
               ASSERT_NEAR(layer.weights(31, 31), 0.0212148, 1.0e-5);
               EXPECT_NEAR(layer.bias(0, 0), 4.68216, 1.0e-5);
               EXPECT_NEAR(layer.bias(15, 0), 0.633401, 1.0e-5);
               ASSERT_NEAR(layer.bias(31, 0), -0.0314119, 1.0e-5);
               break;
      case 3:  EXPECT_EQ(layer.name, "dense_4");
               EXPECT_EQ(layer.activation, "linear");
               EXPECT_EQ(layer.input_size, 32);
               ASSERT_EQ(layer.units, 1);
               EXPECT_NEAR(layer.weights(0, 0), 5.0117, 1.0e-5);
               EXPECT_NEAR(layer.weights(0, 15), -29.7871, 1.0e-4);
               ASSERT_NEAR(layer.weights(0, 31), 0.0211935, 1.0e-5);
               ASSERT_NEAR(layer.bias(0, 0), 3.10022, 1.0e-5);
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
              EXPECT_NEAR(output(0, 0), -2.867314, 1.0e-5);
              EXPECT_NEAR(output(31, 0), -0.493263, 1.0e-5);
              ASSERT_NEAR(output(63, 0), -2.867358, 1.0e-5);
              break;
      case 1: EXPECT_EQ(input.rows(), 64);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 0.0537932, 1.0e-6);
              ASSERT_NEAR(input(63, 0), 0.05379095, 1.0e-6);
              EXPECT_EQ(output.rows(), 32);
              ASSERT_EQ(output.cols(), 1);
              EXPECT_NEAR(output(0, 0), -11.326995, 1.0e-5);
              ASSERT_NEAR(output(31, 0), -8.612443, 1.0e-5);
              break;
      case 2: EXPECT_EQ(input.rows(), 32);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 1.2043238e-5, 1.0e-8);
              ASSERT_NEAR(input(31, 0), 0.000181796, 1.0e-8);
              EXPECT_EQ(output.rows(), 32);
              ASSERT_EQ(output.cols(), 1);
              EXPECT_NEAR(output(0, 0), 4.776572, 1.0e-5);
              ASSERT_NEAR(output(31, 0), -0.0325998, 1.0e-5);
              break;
      case 3: EXPECT_EQ(input.rows(), 32);
              EXPECT_EQ(input.cols(), 1);
              EXPECT_NEAR(input(0, 0), 4.7765725, 1.0e-5);
              ASSERT_NEAR(input(30, 0), 0, 1.0e-5);
              EXPECT_EQ(output.rows(), 1);
              ASSERT_EQ(output.cols(), 1);
              ASSERT_NEAR(output(0, 0), -2.1679472, 1.0e-5);
              break;
    };

    if (layer.activation == "sigmoid") { 
      EigenMatrixF transformed = output.unaryExpr([](double x) -> double {
        double divisor = 1 + std::exp(-x);
        if (divisor == 0)
        {
            divisor = std::numeric_limits<double>::min();
        }
        return 1 / divisor;
      });
      input = transformed;
    } else if (layer.activation == "relu") {
      EigenMatrixF transformed = output.unaryExpr([](double x) -> double {
        return std::max<double>(x, 0);
      });
      input = transformed;
    } else if (layer.activation == "linear") {
      input = output;
    }
    //std::cout << " transformed input =\n" << input << "\n";
    ++layer_num;
  }

  EigenMatrixF output = input.unaryExpr([&pressure_min, &pressure_max](double x) -> double 
    {
      return pressure_min + x * (pressure_max - pressure_min);
    });
  //std::cout << "Prediction = " << output << std::endl;
  ASSERT_NEAR(output(0, 0), -2167947.2, 1.0);

}
