#include <iostream>
#include <string>

#include "H5Cpp.h"
#include <submodules/json/src/json.hpp>

#include <gtest/gtest.h>

TEST(HDF5Tests, readTest)
{
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
    std::cout << config.size() << "\n";
    for (auto it = config.begin(); it != config.end(); ++it) {
      //std::cout << *it << "\n";
      //std::cout << (*it).size() << "\n";
      auto layer_class_name = (*it)["class_name"];
      std::cout << "layer type = " << layer_class_name << "\n";
      auto layer_config = (*it)["config"];
      std::cout << "layer_config = " << layer_config.size() << "\n";
      for (auto l_it = layer_config.begin(); l_it != layer_config.end(); ++l_it) {
        std::cout << "\t" << l_it.key() << " = " << l_it.value() << "\n";
      }
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
          float test3[dims_data[0]][dims_data[1]];
          dataset.read((void *)test3, datatype);
          for (auto kk=0u; kk < dims_data[0]; ++kk) {
            for (auto ll=0u; ll < dims_data[1]; ++ll) {
              std::cout << test3[kk][ll] << " ";
            }
          }
          std::cout << "\n";
        } else {
          float test3[dims_data[0]];
          dataset.read((void *)test3, datatype);
          for (auto kk=0u; kk < dims_data[0]; ++kk) {
            std::cout << test3[kk] << " ";
          }
          std::cout << "\n";
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
