#
# The MIT License
#
# Copyright (c) 2017- Parresia Research Limited, New Zealand
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

SET(MATH_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/IntVec.cpp          
  ${CMAKE_CURRENT_SOURCE_DIR}/Vec.cpp          
  ${CMAKE_CURRENT_SOURCE_DIR}/Matrix2.cpp          
  ${CMAKE_CURRENT_SOURCE_DIR}/Matrix3.cpp          
  ${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cpp          
  ${CMAKE_CURRENT_SOURCE_DIR}/balanc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hqr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/nrutil.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ran.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CubeRoot.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/root6.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/root6_old.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/zrhqr.cpp
)

SET(ELLIP3D_SRCS
  ${ELLIP3D_SRCS}
  ${MATH_SRCS}
  PARENT_SCOPE
)
