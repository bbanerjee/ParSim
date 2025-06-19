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

# CMakeLists.txt for Vaango/src/CCA/Components/MPM/ReactionDiffusion

SET(MPM_Diffusion_ConductivityModels_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/ConductivityModels/BinaryEquation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ConductivityModels/ConductivityEquation.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ConductivityModels/FixedEquation.cc
   )

SET(MPM_Diffusion_DiffusionInterfaces_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionInterfaces/CommonIFConcDiff.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionInterfaces/SDInterfaceModel.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionInterfaces/SimpleDiffusionContact.cc
   )

SET(MPM_Diffusion_DiffusionModels_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/ConstantRate.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/JGConcentrationDiffusion.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/NonLinearDiff1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/NonLinearDiff2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/RFConcDiffusion1MPM.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/DiffusionModels/ScalarDiffusionModel.cc
   )

SET(Vaango_CCA_Components_MPM_Diffusion_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/ScalarDiffusionModelFactory.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/SDInterfaceModelFactory.cc
  ${MPM_Diffusion_ConductivityModels_SRCS}
  ${MPM_Diffusion_DiffusionInterfaces_SRCS}
  ${MPM_Diffusion_DiffusionModels_SRCS}
)

set(LIBS 
    ${LIBS} 
    Vaango_Core_OS
    MPM_Core 
    MPM_Contact
   )

add_library(MPM_Diffusion ${Vaango_CCA_Components_MPM_Diffusion_SRCS})

target_link_libraries(MPM_Diffusion 
  ${LIBS}
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_MPM_Diffusion)
ENDIF(BUILD_SHARED_LIBS)

install(TARGETS
  MPM_Diffusion
)
