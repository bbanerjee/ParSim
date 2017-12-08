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
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/googletest/googletest/include")
set(TEST_LIB gtest_main gtest)
#set(CMAKE_VERBOSE_MAKEFILE 1)
include_directories(${GTEST_INCLUDE_DIR})

#--------------------------------------
# Input files
#--------------------------------------
file(COPY ./4002/axisymmetric_strain_pb_two.000 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#--------------------------------------
# DEMComputeStressStrain
#--------------------------------------
include_directories(${QHULL_INCLUDE_DIRS})
add_executable(test4002 4002/testComputeStressStrain.cpp)

target_link_libraries(test4002
  ${TEST_LIB}
  ELLIP3D_LIB
  ${MPI_LIBRARY}
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${QHULL_LIBRARIES}
  stdc++fs
)

set(UNIT_TEST test4002)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

