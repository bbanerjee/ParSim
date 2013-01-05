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

# CMakeLists.txt for Packages/Uintah/Core/Math

SET(Packages_Uintah_Core_Math_SRCS
  CubicPWI.cc          
  Gaussian.cc          
  Weibull.cc           
  LinAlg.c             
  LinearPWI.cc         
  Mat.c                
  MiscMath.cc          
  MusilRNG.cc          
  PiecewiseInterp.cc   
  TrigTable.cc         
  sci_lapack.cc        
  fft.c                
  ssmult.c             
  FastMatrix.cc        
  Primes.cc            
  Matrix3.cc           
  SymmMatrix3.cc       
  CubeRoot.cc          
  Sparse.cc            
  Short27.cc           
  TangentModulusTensor.cc  
  ${CMAKE_SOURCE_DIR}/Core/Geometry/BBox.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/CompGeom.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/IntVector.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Plane.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Point.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Polygon.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Quaternion.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Ray.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Tensor.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Transform.cc 
  ${CMAKE_SOURCE_DIR}/Core/Geometry/Vector.cc 
  ${CMAKE_SOURCE_DIR}/Core/Disclosure/TypeDescription.cc 
  ${CMAKE_SOURCE_DIR}/Core/Disclosure/TypeUtils.cc 
)

IF(WIN32)
  SET(Packages_Uintah_Core_Math_SRCS ${Packages_Uintah_Core_Math_SRCS} Rand48.cc)
ENDIF(WIN32)

ADD_LIBRARY(Packages_Uintah_Core_Math ${Packages_Uintah_Core_Math_SRCS})

set(LIBS
 ${M_LIBRARY} 
 ${DL_LIBRARY} 
 ${LAPACK_LIBRARY} 
 ${BLAS_LIBRARY} 
 ${F_LIBRARY} 
 ${XML2_LIBRARY} 
 ${MPI_LIBRARY} 
 ${TEEM_LIBRARY} 
 ${PNG_LIBRARY} 
 ${Z_LIBRARY}
)

TARGET_LINK_LIBRARIES(Packages_Uintah_Core_Math
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_Thread
  Packages_Uintah_Core_Persistent
  ${LIBS}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_Core_Math)
ENDIF(BUILD_SHARED_LIBS)
