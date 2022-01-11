#
# The MIT License
#
# Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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
configure_file(tabular_eos.json tabular_eos.json COPYONLY)
configure_file(HydrostaticCompression_PrescribedDeformation.inp 
               HydrostaticCompression_PrescribedDeformation.inp COPYONLY)

configure_file(tabular_linear_elastic.json tabular_linear_elastic.json COPYONLY)
configure_file(tabular_von_mises.json tabular_von_mises.json COPYONLY)
configure_file(tabular_drucker_prager.json tabular_drucker_prager.json COPYONLY)
configure_file(tabular_cap.json tabular_cap.json COPYONLY)
configure_file(DrySand_HydrostatData.json DrySand_HydrostatData.json COPYONLY)
configure_file(UniaxialStrainRotate_PrescribedDeformation.inp 
               UniaxialStrainRotate_PrescribedDeformation.inp COPYONLY)

configure_file(mlp_regression_keras_total_scaled.h5 
               mlp_regression_keras_total_scaled.h5 COPYONLY)

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

add_executable(testTabularData testTabularData.cc)
target_link_libraries(testTabularData
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testTabularData)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testTabularEOS testTabularEOS.cc)
target_link_libraries(testTabularEOS
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testTabularEOS)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testTabularPlasticity testTabularPlasticity.cc)
target_link_libraries(testTabularPlasticity
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testTabularPlasticity)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testNeuralNetTabularPlasticity testNeuralNetTabularPlasticity.cc)
target_link_libraries(testNeuralNetTabularPlasticity
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testNeuralNetTabularPlasticity)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testTabularPlasticityCap testTabularPlasticityCap.cc)
target_link_libraries(testTabularPlasticityCap
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testTabularPlasticityCap)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

