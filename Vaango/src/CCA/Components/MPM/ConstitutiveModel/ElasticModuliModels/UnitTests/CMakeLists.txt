#
# The MIT License
#
# Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

#message(STATUS "Binary dir = ${CMAKE_CURRENT_BINARY_DIR}")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "HDF5 library dir = ${HDF5_LIBRARY_DIRS}")
#message(STATUS "HDF5 library = ${HDF5_LIBRARY}")

#------------------------
# Copy input files
#------------------------
configure_file(table_elastic.json table_elastic.json COPYONLY)
configure_file(mlp_regression_keras_total_scaled.json mlp_regression_keras_total_scaled.json COPYONLY)
configure_file(mlp_regression_keras_total_scaled.h5 mlp_regression_keras_total_scaled.h5 COPYONLY)
configure_file(ARL_Sand_SVR_fit_10_001.json ARL_Sand_SVR_fit_10_001.json COPYONLY)
configure_file(mlp_regression_keras_bulk_elastic_ten_com_scaled_relu_32_400.h5 mlp_regression_keras_bulk_elastic_ten_com_scaled_relu_32_400.h5 COPYONLY)
configure_file(DrySand_ElasticData_Ten_Com.json DrySand_ElasticData_Ten_Com.json COPYONLY)
configure_file(DrySand_BulkModulusData_Ten_Com.json DrySand_BulkModulusData_Ten_Com.json COPYONLY)

#------------------------
# Vaango libs
#------------------------
SET(VAANGO_LIBS 
        Vaango_Core_Thread       
        Vaango_Core_DataArchive 
        Vaango_Core_Grid        
        Vaango_Core_Parallel    
        Vaango_Core_Labels      
        Vaango_Core_Util        
        Vaango_Core_Math        
        Vaango_Core_Exceptions  
        Vaango_Core_OS  
        Vaango_CCA_Ports        
        Vaango_CCA_Components_Parent 
        Vaango_CCA_Components_DataArchiver  
        Vaango_CCA_Components_LoadBalancers 
        Vaango_CCA_Components_Regridder     
        Vaango_Core_ProblemSpec             
        Vaango_CCA_Components_SimulationController 
        Vaango_CCA_Components_Schedulers           
        Vaango_CCA_Components_ProblemSpecification 
        Vaango_CCA_Components_Solvers
        ${HDF5_LIBRARY}
)

#------------------------
# Google tests
#------------------------
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/submodules/googletest/googletest/include")
set(GTEST_LIB gtest_main gtest)
include_directories(${GTEST_INCLUDE_DIR})

add_executable(testElasticModuli_Arena testElasticModuli_Arena.cc)
target_link_libraries(testElasticModuli_Arena
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testElasticModuli_Arena)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testElasticModuli_MetalIso testElasticModuli_MetalIso.cc)
target_link_libraries(testElasticModuli_MetalIso
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testElasticModuli_MetalIso)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testElasticModuli_Tabular testElasticModuli_Tabular.cc)
target_link_libraries(testElasticModuli_Tabular
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testElasticModuli_Tabular)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testElasticModuli_NeuralNet testElasticModuli_NeuralNet.cc)
target_link_libraries(testElasticModuli_NeuralNet
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testElasticModuli_NeuralNet)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST} > ${UNIT_TEST}.out
  VERBATIM
)

set(UNIT_TEST testElasticModuli_NeuralNet_Bulk)
add_executable(${UNIT_TEST} ${UNIT_TEST}.cc)
target_link_libraries(${UNIT_TEST}
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

set(UNIT_TEST testElasticModuli_SupportVector)
add_executable(${UNIT_TEST} ${UNIT_TEST}.cc)
target_link_libraries(${UNIT_TEST}
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

