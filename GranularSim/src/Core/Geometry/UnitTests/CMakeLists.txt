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
# FaceEdge
#--------------------------------------
add_executable(testFaceEdge testFaceEdge.cpp)

target_link_libraries(testFaceEdge
  ${TEST_LIB}
  ELLIP3D_LIB
  ${MPI_LIBRARY}
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ZLIB_LIBRARIES}
)

set(UNIT_TEST testFaceEdge)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

#--------------------------------------
# OrientedBox
#--------------------------------------
add_executable(testOrientedBox testOrientedBox.cpp)

target_link_libraries(testOrientedBox
  ${TEST_LIB}
  ELLIP3D_LIB
  ${MPI_LIBRARY}
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ZLIB_LIBRARIES}
)

set(UNIT_TEST testOrientedBox)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)

#--------------------------------------
# Ellipsoid
#--------------------------------------
add_executable(testEllipsoid testEllipsoid.cpp)

target_link_libraries(testEllipsoid
  ${TEST_LIB}
  ELLIP3D_LIB
  ${MPI_LIBRARY}
  ${Boost_LIBRARIES}
  ${VTK_LIBRARIES}
  ${ZLIB_LIBRARIES}
)

set(UNIT_TEST testEllipsoid)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)
