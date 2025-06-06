#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

# CMakeLists.txt for Vaango/src/CCA/Components/ICE/EOS

SET(Vaango_CCA_Components_ICE_SRCS
  ${Vaango_CCA_Components_ICE_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/EquationOfState.cc    
  ${CMAKE_CURRENT_SOURCE_DIR}/EquationOfStateFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/IdealGas.cc               
  ${CMAKE_CURRENT_SOURCE_DIR}/HardSphereGas.cc          
  ${CMAKE_CURRENT_SOURCE_DIR}/TST.cc                    
  ${CMAKE_CURRENT_SOURCE_DIR}/JWL.cc                    
  ${CMAKE_CURRENT_SOURCE_DIR}/JWLC.cc                   
  ${CMAKE_CURRENT_SOURCE_DIR}/Gruneisen.cc              
  ${CMAKE_CURRENT_SOURCE_DIR}/Tillotson.cc              
  ${CMAKE_CURRENT_SOURCE_DIR}/Thomsen_Hartka_water.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/KnaussSeaWater.cc         
  ${CMAKE_CURRENT_SOURCE_DIR}/KumariDass.cc             
  ${CMAKE_CURRENT_SOURCE_DIR}/Murnaghan.cc              
  ${CMAKE_CURRENT_SOURCE_DIR}/BirchMurnaghan.cc
  PARENT_SCOPE
)
