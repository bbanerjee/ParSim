#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2014-2023 Biswajit Banerjee
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
# CMakeLists.txt for Vaango/src/CCA/Components/MPM/PhysicalBC

SET(Vaango_CCA_Components_MPM_PhysicalBC_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/MPMPhysicalBCFactory.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/ForceBC.cc                
  ${CMAKE_CURRENT_SOURCE_DIR}/MomentBC.cc            
  ${CMAKE_CURRENT_SOURCE_DIR}/VelocityBC.cc            
  ${CMAKE_CURRENT_SOURCE_DIR}/PressureBC.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CrackBC.cc      
  ${CMAKE_CURRENT_SOURCE_DIR}/PolynomialData.cc    
  ${CMAKE_CURRENT_SOURCE_DIR}/FluxBCModelFactory.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/FluxBCModel.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/AutoCycleFluxBC.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/ScalarFluxBC.cc    
  ${CMAKE_CURRENT_SOURCE_DIR}/ArchesHeatFluxBC.cc    
  ${CMAKE_CURRENT_SOURCE_DIR}/HeatFluxBC.cc      
)

add_library(MPM_PhysicalBC ${Vaango_CCA_Components_MPM_PhysicalBC_SRCS})

target_link_libraries(MPM_PhysicalBC 
  ${LIBS}
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_MPM_PhysicalBC)
ENDIF(BUILD_SHARED_LIBS)
