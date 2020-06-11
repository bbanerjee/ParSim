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
configure_file(polyline_sim_data.csv polyline_sim_data.csv COPYONLY)
configure_file(mlp_regression_keras_total_scaled.json mlp_regression_keras_total_scaled.json COPYONLY)
configure_file(mlp_regression_keras_total_scaled.h5 mlp_regression_keras_total_scaled.h5 COPYONLY)

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

add_executable(testHDF5Reader testHDF5Reader.cc)
target_link_libraries(testHDF5Reader
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testHDF5Reader)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testTensorUtils testTensorUtils.cc)
target_link_libraries(testTensorUtils
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testTensorUtils)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testYieldCondUtils testYieldCondUtils.cc)
target_link_libraries(testYieldCondUtils
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testYieldCondUtils)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testEigen testEigen.cc)
target_link_libraries(testEigen
  ${GTEST_LIB}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testEigen)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testBoostAutodiff testBoostAutodiff.cc)
target_link_libraries(testBoostAutodiff
  ${GTEST_LIB}
  ${Boost_LIBRARIES}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testBoostAutodiff)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

add_executable(testBoostClosestPoint testBoostClosestPoint.cc)
target_link_libraries(testBoostClosestPoint
  ${GTEST_LIB}
  ${Boost_LIBRARIES}
  ${VAANGO_LIBS}
)
set(UNIT_TEST testBoostClosestPoint)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)
