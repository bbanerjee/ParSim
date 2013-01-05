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

# CMakeLists.txt for Packages/Uintah/Core/Grid/BoundaryConditions

SET(Packages_Uintah_Core_Grid_SRCS
  ${Packages_Uintah_Core_Grid_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/AnnulusBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BCDataArray.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BCGeomBase.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BCUtils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BoundCondFactory.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/BoundCondReader.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CircleBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/DifferenceBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/EllipseBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RectangleBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/SideBCData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/UnionBCData.cc
  PARENT_SCOPE
)

