#
# The MIT License
#
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

enable_testing()
include(GoogleTest)

SET(VAANGO_LIBS
    Vaango_Core_OS
    Vaango_Core_Exceptions
    Vaango_Core_Parallel
    )

#------------------------
# Google tests
#------------------------
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/submodules/googletest/googletest/include")
set(GTEST_LIB GTest::gtest_main)
include_directories(${GTEST_INCLUDE_DIR})

add_executable(testDir testDir.cc)
target_include_directories(testDir PRIVATE
                           ${CMAKE_SOURCE_DIR}/Core/OS) 
target_link_libraries(testDir
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  ${CMAKE_DL_LIBS}
  Vaango::Options
  Vaango::Warnings
)
set(UNIT_TEST testDir)
add_custom_command(
  TARGET ${UNIT_TEST}
  POST_BUILD
  COMMAND ${UNIT_TEST}
)
gtest_discover_tests(${UNIT_TEST})
