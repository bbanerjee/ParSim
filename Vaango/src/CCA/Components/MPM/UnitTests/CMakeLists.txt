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

#------------------------
# Vaango libs
#------------------------
SET(VAANGO_LIBS 
        Vaango_Core_Grid 
        Vaango_Core_DataArchive 
        Vaango_Core_Parallel    
        Vaango_Core_Util        
        Vaango_Core_Math        
        Vaango_Core_Exceptions  
        Vaango_Core_OS  
        Vaango_Core_GeometryPiece
        Vaango_CCA_Ports        
        Vaango_CCA_Components_Parent 
        Vaango_CCA_Components_DataArchiver  
        Vaango_CCA_Components_LoadBalancers 
        Vaango_CCA_Components_Regridder     
        Vaango_CCA_Components_Schedulers     
        Vaango_Core_ProblemSpec             
        Vaango_CCA_Components_SimulationController 
        Vaango_CCA_Components_ProblemSpecification 
        Vaango_CCA_Components_Solvers
        Vaango_CCA_Components_Models
        Vaango_CCA_Components_MPM
)

#------------------------
# Google tests
#------------------------
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/submodules/googletest/googletest/include")
set(GTEST_LIB gtest_main gtest)
include_directories(${GTEST_INCLUDE_DIR})

add_executable(testHybridDEM testHybridDEM.cc)
target_link_libraries(testHybridDEM PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  MPM_ParticleCreator
  MPM_GradientComputer
  MPM_ConstitutiveModel
  MPM_PhysicalBC
  MPM_MMS
  MPM_Diffusion
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)
