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
# CMakeLists.txt for Vaango/src/R_Tester
# Based on the deal.II CMakeLists.txt

# Set up the regression tests
#   create_goldstandards - creates the benchmark simulations and saves data
#   clean_tests          - cleans up all test data

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

if (DEFINED HAVE_REGRESSION_TESTS)

  MESSAGE(STATUS "Setting up regression tests")

  # Write minimalistic CTestTestfile.cmake files to CMAKE_BINARY_DIR and
  # CMAKE_BINARY_DIR/tests:
  #
  FILE(WRITE ${CMAKE_BINARY_DIR}/CTestTestfile.cmake "SUBDIRS(R_Tester)")

endif()

FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake "")

#
# Find all regression tests, i.e., every directory that contains a
# CMakeLists.txt file
#
SET(_categories)
FILE(GLOB _dirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*
  )
FOREACH(_dir ${_dirs})
  IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${_dir}/CMakeLists.txt)
    LIST(APPEND _categories ${_dir})
  ENDIF()
ENDFOREACH()

#
# Custom targets for the regression tests:
#

# Setup tests:
ADD_CUSTOM_TARGET(create_goldstandard)

# Remove all tests:
ADD_CUSTOM_TARGET(clean_tests)

# Loop through the directories
FOREACH(_category ${_categories})
  SET(_category_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_category})

  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category})

  SET(_command ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_category_dir} > /dev/null)

  ADD_CUSTOM_TARGET(create_goldstandard_${_category}
    COMMAND ${_command}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(create_goldstandard create_goldstandard_${_category})

  ADD_CUSTOM_TARGET(clean_tests_${_category}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_category}
    COMMENT "Processing tests/${_category}"
    )
  ADD_DEPENDENCIES(clean_tests clean_tests_${_category})

  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
    "SUBDIRS(${_category})\n"
    )
ENDFOREACH()

MESSAGE(STATUS "Setting up regression tests - Done")

  


