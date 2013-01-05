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

# CMakeLists.txt for Packages/Uintah/CCA/Components/MPM/ConstitutiveModel/fortran

set(FORTRAN_CM_FLAGS "${CMAKE_Fortran_FLAGS} -DI1SOFTLIM -DNEWXSV -DDILAT_LIM")

set(FORTRAN_CM_SRC 
  ${CMAKE_CURRENT_SOURCE_DIR}/Hooke.F  
  ${CMAKE_CURRENT_SOURCE_DIR}/Kayenta_calcs.F  
  ${CMAKE_CURRENT_SOURCE_DIR}/MIGUtilsF.F  
  ${CMAKE_CURRENT_SOURCE_DIR}/MIGUtils.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/diamm.F
)

set_source_files_properties(
  ${FORTRAN_CM_SRC}
  PROPERTIES
    COMPILE_FLAGS ${FORTRAN_CM_FLAGS} 
)

SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${FORTRAN_CM_SRC}
  PARENT_SCOPE
)

