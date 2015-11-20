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
#   goldstandards - creates the benchmark simulations and saves data
#   run_tests     - runs regression tests
#   clean_tests   - cleans up all test data

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
# Find all regression tests, i.e., all the Python scripts in this directory
#
SET(_test_files)
FILE(GLOB _components RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py
  )
FOREACH(_component ${_components})
  string(REPLACE ".py" "" _component_name ${_component})
  if (NOT ${_component_name} MATCHES "__init__") 
    LIST(APPEND _test_files ${_component_name})
    message(STATUS "component_name = ${_component_name}")
  endif()
ENDFOREACH()

#
# Custom targets for the regression tests:
#

# Setup goldstandards:
ADD_CUSTOM_TARGET(goldstandard)

# Run all tests
ADD_CUSTOM_TARGET(run_tests)

# Remove all tests:
ADD_CUSTOM_TARGET(clean_tests)

# Set the python path and python command
set(PYTHON_COMMAND python)
set(_pythonpath "${CMAKE_CURRENT_SOURCE_DIR}/toplevel:${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "python path = ${_pythonpath}")
set(ENV{PYTHONPATH} ${_pythonpath}:\$ENV{PYTHONPATH})

# Set up variable names to match Uintah
set(IS_DEBUG "no")
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(IS_DEBUG "yes")
endif()
message(STATUS "build type = ${CMAKE_BUILD_TYPE}; IS_DEBUG = ${IS_DEBUG}")

set(OBJTOP_ABS ${CMAKE_BINARY_DIR})
message(STATUS "OBJTOP_ABS = ${OBJTOP_ABS}")

set(SRCTOP_ABS ${CMAKE_SOURCE_DIR})
message(STATUS "SRCTOP_ABS = ${SRCTOP_ABS}")

set(SCI_MALLOC_ON "no")  # Default to no for now

# Loop through the test files
FOREACH(_test_file ${_test_files})

  SET(_test_file_dir ${CMAKE_CURRENT_SOURCE_DIR}/${_test_file})

  FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_test_file})

  SET(_command ${CMAKE_COMMAND} -G${CMAKE_GENERATOR} ${_options} ${_test_file_dir} > /dev/null)

  # Create gold standards
  ADD_CUSTOM_TARGET(goldstandard_${_test_file}

    # Create gold standard directory
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/gold_standards
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/gold_standards

    # Run the python command
    COMMAND export PYTHONPATH=${_pythonpath} &&
            ${PYTHON_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/toplevel/generateGoldStandards.py 
            -d ${IS_DEBUG} -b ${OBJTOP_ABS} -s ${SRCTOP_ABS} -m ${SCI_MALLOC_ON} 
            -t ${_test_file}

    COMMENT "Creating gold standards/${_test_file}"
    )
  ADD_DEPENDENCIES(goldstandard goldstandard_${_test_file})

  ADD_CUSTOM_TARGET(run_tests_${_test_file}
    COMMAND ${_command}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_test_file}
    COMMENT "Processing tests/${_test_file}"
    )
  ADD_DEPENDENCIES(run_tests run_tests_${_test_file})

  ADD_CUSTOM_TARGET(clean_tests_${_test_file}
    COMMAND ${CMAKE_COMMAND} -E remove_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_test_file}
    COMMAND ${CMAKE_COMMAND} -E make_directory
      ${CMAKE_CURRENT_BINARY_DIR}/${_test_file}
    COMMENT "Processing tests/${_test_file}"
    )
  ADD_DEPENDENCIES(clean_tests clean_tests_${_test_file})

  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/CTestTestfile.cmake
    "SUBDIRS(${_test_file})\n"
    )
ENDFOREACH()

MESSAGE(STATUS "Setting up regression tests - Done")

  


