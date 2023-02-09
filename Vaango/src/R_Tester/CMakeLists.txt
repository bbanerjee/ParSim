#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2023 Biswajit Banerjee
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
#   gold_standard - creates the benchmark simulations and saves data
#   localRT       - runs regression tests
#   altLocalRT       - runs regression tests

# *WARNING* Uses symbolic links and "export" and will not work with non linux systems.

CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

#######################################################################
# Set up python path/command and variable names to match Uintah
#######################################################################
set(PYTHON_COMMAND python)
set(_pythonpath "${CMAKE_CURRENT_SOURCE_DIR}/toplevel:${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "python path = ${_pythonpath}")

set(IS_DEBUG "no")
set(DBG_OPT "opt")
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(IS_DEBUG "yes")
  set(DBG_OPT "dbg")
endif()
message(STATUS "build type = ${CMAKE_BUILD_TYPE}; IS_DEBUG = ${IS_DEBUG}")

set(SCI_MALLOC_ON "no")  #**TODO** Find out how to check memory usage

set(MAX_PARALLELISM 8)   # Default to 8 (identical to Uintah)

set(OBJTOP_ABS ${CMAKE_BINARY_DIR})
message(STATUS "OBJTOP_ABS = ${OBJTOP_ABS}")

set(SRCTOP_ABS ${CMAKE_SOURCE_DIR})
message(STATUS "SRCTOP_ABS = ${SRCTOP_ABS}")

# Set the execuatble directory
set(VAANGO_DIR ${CMAKE_BINARY_DIR}/StandAlone)

# Set the inputs directory
set(INPUTS_DIR ${CMAKE_SOURCE_DIR}/StandAlone/inputs)
message(STATUS "INPUTS_ABS = ${INPUTS_DIR}")

# Set the gold standard directory
set(GOLDSTANDARD_DIR ${CMAKE_CURRENT_BINARY_DIR}/gold_standards)

# Set the RT output directory
set(LOCALRT_DIR ${CMAKE_CURRENT_BINARY_DIR}/localRTData)

#######################################################################
# Find all regression tests, i.e., all the Python scripts in the
# R_Tester directory
#######################################################################
SET(_test_names)
SET(_test_names_py)
FILE(GLOB _components_py RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py
  )
FOREACH(_component_py ${_components_py})
  string(REPLACE ".py" "" _component_name ${_component_py})
  if (NOT ${_component_name} MATCHES "__init__") 
    LIST(APPEND _test_names ${_component_name})
    LIST(APPEND _test_names_py ${_component_py})
    message(STATUS "component_name = ${_component_name}")
  endif()
ENDFOREACH()

#######################################################################
# Custom targets for the gold standards:
#######################################################################
# Setup goldstandards:
ADD_CUSTOM_TARGET(gold_standard)

# Loop through the test files
FOREACH(_test ${_test_names})

  # Create gold standards
  ADD_CUSTOM_TARGET(gold_standard_${_test}

    # Create gold standard directory
    COMMAND ${CMAKE_COMMAND} -E make_directory ${GOLDSTANDARD_DIR}

    # Run the python command
    COMMAND export PYTHONPATH=${_pythonpath} &&
            ${PYTHON_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/toplevel/generateGoldStandards.py 
            -d ${IS_DEBUG} -b ${OBJTOP_ABS} -s ${SRCTOP_ABS} -m ${SCI_MALLOC_ON} 
            -t ${_test} -o ${GOLDSTANDARD_DIR}

    COMMENT "Creating gold standards/${_test}"
    )
  ADD_DEPENDENCIES(gold_standard gold_standard_${_test})

ENDFOREACH()

MESSAGE(STATUS "Setting up gold standards - Done")

#######################################################################
# Custom target for the local regression tests:
#######################################################################
# Run all tests
ADD_CUSTOM_TARGET(localRT)

#message(STATUS "test_names_py = ${_test_names_py}")

# Loop through the test files
FOREACH(_test_file ${_test_names_py})

  # Extract the test name if necessary
  string(REPLACE ".py" "" _test ${_test_file})
  #message(STATUS "test_file = ${_test_file}; test = ${_test}")

  # Get absolute path for the python script
  set(TEST_ABS ${CMAKE_CURRENT_SOURCE_DIR}/${_test_file})

  # Run local regression tests
  ADD_CUSTOM_TARGET(localRT${_test}

    # Create regression test directory
    COMMAND ${CMAKE_COMMAND} -E make_directory ${LOCALRT_DIR}

    # Run the python command
    COMMAND export PYTHONPATH=${_pythonpath} &&
            export OS=Linux &&
            export SCI_DEBUG="" &&
            export LOCAL_OR_NIGHTLY_TEST=local &&
            export WHICH_TESTS=local &&
            export SCI_MALLOC_ENABLED=no &&
            export LOCALRT_DIR=${LOCALRT_DIR} &&
            ${PYTHON_COMMAND} -u ${TEST_ABS}
            ${VAANGO_DIR} ${INPUTS_DIR} ${GOLDSTANDARD_DIR} 
            ${DBG_OPT} ${MAX_PARALLELISM}

    COMMENT "Running regression tests: ${_test_file}"
    )
  ADD_DEPENDENCIES(localRT localRT${_test})

ENDFOREACH()

#######################################################################
# Custom target for the local regression tests gainst alternative gold_standard:
#######################################################################
# Run all tests
ADD_CUSTOM_TARGET(altLocalRT)

# Set the gold standard directory
set(ALT_GOLDSTANDARD_DIR ${CMAKE_BINARY_DIR}/../gold_standards/OptPreAMR)

# Set the RT output directory
set(ALT_LOCALRT_DIR ${CMAKE_CURRENT_BINARY_DIR}/altLocalRTData)


# Loop through the test files
FOREACH(_test_file ${_test_names_py})

  # Extract the test name if necessary
  string(REPLACE ".py" "" _test ${_test_file})
  #message(STATUS "test_file = ${_test_file}; test = ${_test}")

  # Get absolute path for the python script
  set(TEST_ABS ${CMAKE_CURRENT_SOURCE_DIR}/${_test_file})

  # Run local regression tests
  ADD_CUSTOM_TARGET(altLocalRT${_test}

    # Run a shell command
    #COMMAND echo "Starting with ${_test_file}"

    # Create regression test directory
    COMMAND ${CMAKE_COMMAND} -E make_directory ${ALT_LOCALRT_DIR}

    # Run the python command
    COMMAND export PYTHONPATH=${_pythonpath} &&
            export OS=Linux &&
            export SCI_DEBUG="" &&
            export LOCAL_OR_NIGHTLY_TEST=local &&
            export WHICH_TESTS=local &&
            export SCI_MALLOC_ENABLED=no &&
            export LOCALRT_DIR=${ALT_LOCALRT_DIR} &&
            ${PYTHON_COMMAND} -u ${TEST_ABS}
            ${VAANGO_DIR} ${INPUTS_DIR} ${ALT_GOLDSTANDARD_DIR} 
            ${DBG_OPT} ${MAX_PARALLELISM}

    # Run a shell command
    #COMMAND echo "Done with ${_test_file}"

    COMMENT "Running regression tests: ${_test_file}"
    )
  ADD_DEPENDENCIES(altLocalRT altLocalRT${_test})

ENDFOREACH()

MESSAGE(STATUS "Setting up regression tests - Done")

  


