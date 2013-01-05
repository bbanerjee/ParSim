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

# CMakeLists.txt for Packages/Uintah/Core/Datatypes

SET(Packages_Uintah_Core_Datatypes_SRCS
  Clipper.cc                    
  Color.cc                      
  ColumnMatrix.cc               
  Datatype.cc                   
  DenseColMajMatrix.cc          
  DenseMatrix.cc                
  Field.cc                      
  Image.cc                      
  Matrix.cc                     
  MatrixOperations.cc           
  Mesh.cc                       
  PropertyManager.cc            
  SearchGrid.cc                 
  SparseRowMatrix.cc            
  String.cc                     
  FieldInterfaceAux.cc          
  TypeName.cc                   
  Archive.cc                  
  ScalarParticles.cc     
  VectorParticles.cc     
  TensorParticles.cc     
  PSet.cc       
  VariableCache.cc     
)

if (HAVE_TEEM)
  SET(Packages_Uintah_Core_Datatypes_SRCS
    ${Packages_Uintah_Core_Datatypes_SRCS}
    NrrdData.cc                   
  )
  include_directories(${TEEM_INCLUDE})
else ()
  include_directories(${BLAS_INCLUDE})
endif ()

ADD_LIBRARY(Packages_Uintah_Core_Datatypes ${Packages_Uintah_Core_Datatypes_SRCS})

set(LIBS  
  ${M_LIBRARY} 
  ${BLAS_LIBRARY} 
  ${F_LIBRARY}           
  ${TEEM_LIBRARY} 
  ${PNG_LIBRARY} 
  ${Z_LIBRARY}         
  ${BLAS_LIBRARY} 
  ${LAPACK_LIBRARY} 
  ${THREAD_LIBRARY} 
  ${MPI_LIBRARY}
)

TARGET_LINK_LIBRARIES(Packages_Uintah_Core_Datatypes
  Packages_Uintah_Core_Basis        
  Packages_Uintah_Core_Exceptions   
  Packages_Uintah_Core_Grid   
  Packages_Uintah_Core_Math         
  Packages_Uintah_Core_Persistent   
  Packages_Uintah_Core_Thread       
  Packages_Uintah_Core_Util         
  Packages_Uintah_Core_ProblemSpec
  ${LIBS}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Packages_Uintah_Core_Datatypes)
ENDIF(BUILD_SHARED_LIBS)
