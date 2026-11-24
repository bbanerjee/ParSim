#
# The MIT License
#
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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
# Vaango libs
#------------------------
SET(VAANGO_LIBS 
  Vaango_Core_DataArchive 
  Vaango_Core_Grid        
  Vaango_Core_Parallel    
  Vaango_Core_Util        
  Vaango_Core_Math        
  Vaango_Core_Malloc        
  Vaango_Core_Exceptions  
  Vaango_Core_OS  
  Vaango_Core_ProblemSpec             
  Vaango_CCA_Ports        
  Vaango_CCA_Components_Parent 
  Vaango_CCA_Components_DataArchiver  
  Vaango_CCA_Components_LoadBalancers 
  Vaango_CCA_Components_Regridder     
  Vaango_CCA_Components_Schedulers           
  Vaango_CCA_Components_SimulationController 
  Vaango_CCA_Components_ProblemSpecification 
  Vaango_CCA_Components_Solvers
  ${HDF5_TARGET}
  $<$<BOOL:${HAVE_HYPRE}>:${HYPRE_TARGET}>
  $<$<BOOL:${HAVE_KOKKOS}>:Kokkos::kokkos>
)

#------------------------
# Google tests
#------------------------

set(GTEST_LIB gtest_main gtest)

if (HAVE_KOKKOS)
  set(GPU_SRC ${CMAKE_SOURCE_DIR}/CCA/Components/Schedulers/GPUDataWarehouse.cc)
endif()

set(UNIT_TEST testPressure_Air)
set(SOURCES ${UNIT_TEST}.cc)
add_executable(${UNIT_TEST} ${SOURCES} ${GPU_SRC})
target_link_libraries(${UNIT_TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testPressure_Water testPressure_Water.cc ${GPU_SRC})
target_link_libraries(testPressure_Water PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testPressure_Water)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testPressure_Granite testPressure_Granite.cc ${GPU_SRC})
target_link_libraries(testPressure_Granite PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testPressure_Granite)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

set(TEST testEOS_BorjaT)
set(SOURCES ${TEST}.cc)
add_executable(${TEST} ${SOURCES} ${GPU_SRC})
target_link_libraries(${TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST ${TEST})
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testMieGruneisenEOS testMieGruneisenEOS.cc ${GPU_SRC})
target_link_libraries(testMieGruneisenEOS PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testMieGruneisenEOS)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testHyperElasticEOS testHyperElasticEOS.cc ${GPU_SRC})
target_link_libraries(testHyperElasticEOS PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testHyperElasticEOS)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testDefaultHypoElasticEOS testDefaultHypoElasticEOS.cc ${GPU_SRC})
target_link_libraries(testDefaultHypoElasticEOS PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testDefaultHypoElasticEOS)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()

add_executable(testBorjaEOS testBorjaEOS.cc ${GPU_SRC})
target_link_libraries(testBorjaEOS PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ConstitutiveModel
  MPM_Diffusion
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testBorjaEOS)
if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()
