#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Wasatch/Expressions

SET(CCA_Components_Wasatch_Expressions_SRCS
  BasicExprBuilder.cc
  ConvectiveFlux.cc    
  DiffusiveFlux.cc    
  DiffusiveVelocity.cc
  Dilatation.cc
  MomentumPartialRHS.cc 
  MomentumRHS.cc 
  MonolithicRHS.cc
  PrimVar.cc
  ScalarRHS.cc
  ScalabilityTestSrc.cc
  SetCurrentTime.cc
  Stress.cc 
  VelocityMagnitude.cc 
  Vorticity.cc 
  Pressure.cc     
  PoissonExpression.cc
)

ADD_LIBRARY(CCA_Components_Wasatch_Expressions ${CCA_Components_Wasatch_Expressions_SRCS})

TARGET_LINK_LIBRARIES(CCA_Components_Wasatch_Expressions
  Packages_Uintah_CCA_Components_Schedulers 
  Packages_Uintah_CCA_Ports                 
  Packages_Uintah_Core_Exceptions           
  Packages_Uintah_Core_IO                   
  Packages_Uintah_Core_GeometryPiece        
  Packages_Uintah_Core_Grid                 
  Packages_Uintah_Core_Math                 
  Packages_Uintah_Core_Parallel             
  Packages_Uintah_Core_ProblemSpec          
  Packages_Uintah_Core_Util                 
  ${XML2_LIBRARY} 
  ${MPI_LIBRARY} 
  ${M_LIBRARY}        
  ${EXPRLIB_LIBRARY} 
  ${SPATIALOPS_LIBRARY}           
  ${TABPROPS_LIBRARY} 
  ${HDF5_LIBRARY}                
  ${BOOST_LIBRARY} 
  ${BLAS_LIBRARY} 
  ${LAPACK_LIBRARY} 
)

include_directories(
  ${MPI_INCLUDE_PATH}
  ${SPATIALOPS_INCLUDE} 
  ${EXPRLIB_INCLUDE} 
  ${HDF5_INCLUDE} 
  ${TABPROPS_INCLUDE} 
  ${BOOST_INCLUDE} 
  ${LAPACK_INCLUDE}
)

add_subdirectory(MMS)
add_subdirectory(PBE)
add_subdirectory(PostProcessing)
add_subdirectory(Turbulence)


IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_CCA_Components_Wasatch_Expressions)
ENDIF(BUILD_SHARED_LIBS)

