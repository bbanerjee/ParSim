# 
#  The MIT License
# 
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#  Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

SET(Vaango_CCA_Components_ICE_SRCS
  ICE.cc 
  AMRICE.cc 
  impICE.cc 
  impAMRICE.cc 
  customInitialize.cc   
)

add_subdirectory(Core)
add_subdirectory(Advection)
add_subdirectory(CustomBCs)
add_subdirectory(EOS)
add_subdirectory(Materials)
add_subdirectory(SpecificHeatModel)
add_subdirectory(TurbulenceModel)
add_subdirectory(WallShearStressModel)

add_library(Vaango_CCA_Components_ICE ${Vaango_CCA_Components_ICE_SRCS})

SET(LIBS
  ICE_Core
  Vaango_CCA_Components_Models
  Vaango_CCA_Components_OnTheFlyAnalysis
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_Math
  Vaango_Core_GeometryPiece
  Vaango_Core_ProblemSpec
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
)

IF(APPLE)
  SET(LIBS ${LIBS} Vaango_CCA_Components_MPM)
ENDIF(APPLE)

IF(HYPRE_TARGET)
  SET(LIBS ${LIBS} ${HYPRE_TARGET})
ENDIF(HYPRE_TARGET)

target_link_libraries(Vaango_CCA_Components_ICE PRIVATE
  ${LIBS}
  LibXml2::LibXml2
  MPI::MPI_CXX
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_CCA_Components_ICE)
ENDIF(BUILD_SHARED_LIBS)

install(TARGETS
  Vaango_CCA_Components_ICE
)