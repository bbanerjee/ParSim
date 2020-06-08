#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2020 Parresia Research Limited, New Zealand
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
#-----------------------------------------------------------------------------
# For reading JSON input files
#-----------------------------------------------------------------------------
#set(JSON_INCLUDES ${Vaango_SOURCE_DIR}/submodules/json/src)
#message(STATUS "JSON include dir = ${JSON_INCLUDES}")
#include_directories(${JSON_INCLUDES})
#set(CMAKE_VERBOSE_MAKEFILE 1)

#-----------------------------------------------------------------------------
# The models
#-----------------------------------------------------------------------------
SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVariableModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_SoilModelBrannonKappa.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_BorjaPressure.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_Arena.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/InternalVar_TabularCap.cc 
  PARENT_SCOPE
)

#message(STATUS "Compiler is ${CMAKE_CXX_COMPILER_ID}")
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_COMPILER_IS_INTEL 1)
endif()
IF(NOT CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_COMPILER_IS_INTEL) 
  add_subdirectory(UnitTests)
ENDIF()
