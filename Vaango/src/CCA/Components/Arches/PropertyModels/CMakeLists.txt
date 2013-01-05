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

# CMakeLists.txt for Packages/Uintah/CCA/Components/Arches/PropertyModels

SET(Packages_Uintah_CCA_Components_Arches_SRCS
  ${Packages_Uintah_CCA_Components_Arches_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/PropertyModelBase.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/PropertyModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ABSKP.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ExtentRxn.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TabStripFactor.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/EmpSoot.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/AlgebraicScalarDiss.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ScalarVarianceScaleSim.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/HeatLoss.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ConstProperty.cc
  PARENT_SCOPE
)

