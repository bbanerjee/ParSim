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

set(VAANGO_LIBS
    Vaango_Core_ProblemSpec
    Vaango_Core_Util
    Vaango_Core_Exceptions
    Vaango_Core_Math
    Vaango_Core_Malloc
)

set(GTEST_LIB GTest::gtest_main)

#---------------------------------------------------------
# testArrhenius
#---------------------------------------------------------
set(UNIT_TEST testArrhenius)
# Compile Arrhenius.cc directly to avoid linking against the entire Models library
add_executable(${UNIT_TEST} testArrhenius.cc ../Arrhenius.cc)
target_link_libraries(${UNIT_TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)

if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()
gtest_discover_tests(${UNIT_TEST})

#---------------------------------------------------------
# testModifiedArrhenius
#---------------------------------------------------------
set(UNIT_TEST testModifiedArrhenius)
add_executable(${UNIT_TEST} testModifiedArrhenius.cc ../ModifiedArrhenius.cc)
target_link_libraries(${UNIT_TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)

if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()
gtest_discover_tests(${UNIT_TEST})

#---------------------------------------------------------
# testNthOrderModel
#---------------------------------------------------------
set(UNIT_TEST testNthOrderModel)
add_executable(${UNIT_TEST} testNthOrderModel.cc ../NthOrderModel.cc)
target_link_libraries(${UNIT_TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)

if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()
gtest_discover_tests(${UNIT_TEST})

#---------------------------------------------------------
# testPowerModel
#---------------------------------------------------------
set(UNIT_TEST testPowerModel)
add_executable(${UNIT_TEST} testPowerModel.cc ../PowerModel.cc)
target_link_libraries(${UNIT_TEST} PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)

if (RUN_UNIT_TESTS_ON_BUILD)
  add_custom_command(
    TARGET ${UNIT_TEST}
    POST_BUILD
    COMMAND ${UNIT_TEST}
  )
endif()
gtest_discover_tests(${UNIT_TEST})
