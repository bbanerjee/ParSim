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

# CMakeLists.txt for Packages/Uintah/Core/Grid/Variables

SET(Packages_Uintah_Core_Grid_SRCS
  ${Packages_Uintah_Core_Grid_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/Iterator.cc                   
  ${CMAKE_CURRENT_SOURCE_DIR}/CellIterator.cc               
  ${CMAKE_CURRENT_SOURCE_DIR}/NodeIterator.cc               
  ${CMAKE_CURRENT_SOURCE_DIR}/GridIterator.cc               
  ${CMAKE_CURRENT_SOURCE_DIR}/GridSurfaceIterator.cc        
  ${CMAKE_CURRENT_SOURCE_DIR}/ListOfCellsIterator.cc        
  ${CMAKE_CURRENT_SOURCE_DIR}/DifferenceIterator.cc         
  ${CMAKE_CURRENT_SOURCE_DIR}/UnionIterator.cc              
  ${CMAKE_CURRENT_SOURCE_DIR}/ComputeSet.cc                 
  ${CMAKE_CURRENT_SOURCE_DIR}/ComputeSet_special.cc         
  ${CMAKE_CURRENT_SOURCE_DIR}/GridVariableBase.cc           
  ${CMAKE_CURRENT_SOURCE_DIR}/LocallyComputedPatchVarMap.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ParticleSubset.cc             
  ${CMAKE_CURRENT_SOURCE_DIR}/ParticleVariableBase.cc       
  ${CMAKE_CURRENT_SOURCE_DIR}/ParticleVariable_special.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/PerPatchBase.cc               
  ${CMAKE_CURRENT_SOURCE_DIR}/PSPatchMatlGhost.cc           
  ${CMAKE_CURRENT_SOURCE_DIR}/PSPatchMatlGhostRange.cc      
  ${CMAKE_CURRENT_SOURCE_DIR}/ReductionVariableBase.cc      
  ${CMAKE_CURRENT_SOURCE_DIR}/ReductionVariable_special.cc  
  ${CMAKE_CURRENT_SOURCE_DIR}/ScrubItem.cc                  
  ${CMAKE_CURRENT_SOURCE_DIR}/SoleVariableBase.cc           
  ${CMAKE_CURRENT_SOURCE_DIR}/Stencil7.cc                   
  ${CMAKE_CURRENT_SOURCE_DIR}/Stencil4.cc                   
  ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cc                      
  ${CMAKE_CURRENT_SOURCE_DIR}/ugc_templates.cc              
  ${CMAKE_CURRENT_SOURCE_DIR}/VarLabel.cc                   
  ${CMAKE_CURRENT_SOURCE_DIR}/Variable.cc                   
  PARENT_SCOPE
)

