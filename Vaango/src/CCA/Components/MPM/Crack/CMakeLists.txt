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

# CMakeLists.txt for Vaango/src/CCA/Components/MPM/Crack

SET(Vaango_CCA_Components_MPM_SRCS
  ${Vaango_CCA_Components_MPM_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/ReadAndDiscretizeCracks.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/CrackGeometryFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/CrackGeometry.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/NullCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/QuadCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/CurvedQuadCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/TriangularCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ArcCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/EllipticCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/PartialEllipticCrack.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ParticleNodePairVelocityField.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/CrackSurfaceContact.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/FractureParametersCalculation.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/CrackPropagation.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MoveCracks.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/UpdateCrackFront.cc 
  PARENT_SCOPE
)

