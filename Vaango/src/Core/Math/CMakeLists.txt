#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2023 Biswajit Banerjee
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

# CMakeLists.txt for Vaango/src/Core/Math

#set(CMAKE_VERBOSE_MAKEFILE 1)

SET(Vaango_Core_Math_SRCS
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
  Matrix3Rotation.cc           
  SymmMatrix3.cc       
  SymmMatrix6.cc       
  CubeRoot.cc          
  Sparse.cc            
  Short27.cc           
  Int130.cc           
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
  SET(Vaango_Core_Math_SRCS ${Vaango_Core_Math_SRCS} Rand48.cc)
ENDIF(WIN32)

ADD_LIBRARY(Vaango_Core_Math ${Vaango_Core_Math_SRCS})

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

TARGET_LINK_LIBRARIES(Vaango_Core_Math
  Vaango_Core_Exceptions
  Vaango_Core_Util
  Vaango_Core_Thread
  Vaango_Core_Persistent
  ${LIBS}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_Core_Math)
ENDIF(BUILD_SHARED_LIBS)
