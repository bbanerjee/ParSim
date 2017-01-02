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

# Override this from configure if you want to install -public (to VisIt 
# install dir instead of user's .visit directory).
set(VISIT_PLUGIN_INSTALL_TYPE "-private")
set(VISIT_PLUGIN_SRC_DIR ${Vaango_SOURCE_DIR}/VisIt/udaReaderMTMD) 
set(VISIT_PLUGIN_BIN_DIR ${Vaango_BINARY_DIR}/VisIt/udaReaderMTMD) 
set(UDA2VIS_SRC_DIR ${Vaango_SOURCE_DIR}/StandAlone/tools/uda2vis) 
set(UDA2VIS_BIN_DIR ${Vaango_BINARY_DIR}/StandAlone/tools/uda2vis) 

if (WIN32)
  # Windows untested - BB
  set(LN_CMD mklink)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDFileFormat.h ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDFileFormat.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.C ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.C WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.h ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  execute_process(COMMAND $LN_CMD$ ${UDA2VIS_BIN_DIR}/udaData.h ${UDA2VIS_SRC_DIR}/udaData.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
else ()
  set(LN_CMD ln -fs)
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDFileFormat.h ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDFileFormat.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  message(STATUS "-- Created symbolic link to avtudaReaderMTMDFileFormat.h")
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.C ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.C WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  message(STATUS "-- Created symbolic link to avtudaReaderMTMDOptions.C")
  execute_process(COMMAND ln -fs ${VISIT_PLUGIN_SRC_DIR}/avtudaReaderMTMDOptions.h ${VISIT_PLUGIN_BIN_DIR}/avtudaReaderMTMDOptions.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  message(STATUS "-- Created symbolic link to avtudaReaderMTMDOptions.h")
  execute_process(COMMAND ln -fs ${UDA2VIS_SRC_DIR}/udaData.h ${UDA2VIS_BIN_DIR}/udaData.h WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
  message(STATUS "-- Created symbolic link to udaData.h")
endif ()


# Execute the xml2cmake command to create CMake files
# execute_process(COMMAND ${VISIT_DIR}/bin/xml2cmake ${VISIT_PLUGIN_INSTALL_TYPE} -clobber ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
message(STATUS "-- Running ${VISIT_DIR}/bin/xml2cmake on ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml")
execute_process(COMMAND ${VISIT_DIR}/bin/xml2cmake ${VISIT_PLUGIN_INSTALL_TYPE} -clobber ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml 
  WORKING_DIRECTORY ${VISIT_PLUGIN_BIN_DIR} 
  RESULT_VARIABLE visit_xml2cmake
  OUTPUT_VARIABLE visit_xml2cmake_output)
message(STATUS "xml2cmake: ${visit_xml2cmake} : ${visit_xml2cmake_output}")

execute_process(COMMAND ${VISIT_DIR}/bin/xml2info -clobber ${VISIT_PLUGIN_SRC_DIR}/udaReaderMTMD.xml WORKING_DIRECTORY "${VISIT_PLUGIN_BIN_DIR}" OUTPUT_QUIET)
#execute_process(COMMAND cmake -DVISIT_DISABLE_SETTING_COMPILER:BOOL=TRUE -DCMAKE_CXX_COMPILER:FILEPATH=${CXX} -DCMAKE_CXX_FLAGS:STRING="-I${Vaango_BINARY_DIR} -I${Vaango_SOURCE_DIR} ${CXXFLAGS}" ${VISIT_PLUGIN_BIN_DIR})

include(ExternalProject)
include(${VISIT_DIR}/CMake/PluginMacros.cmake)
set(VISIT_CXX_FLAGS_C -DCMAKE_CXX_FLAGS:STRING=-I${Vaango_BINARY_DIR})
set(VISIT_INSTALL_FLAGS -DVISIT_SOURCE_DIR=${VISIT_DIR} -DCMAKE_INSTALL_PREFIX=${VISIT_PLUGIN_BIN_DIR})
set(VISIT_CXX_FLAGS "${VISIT_CXX_FLAGS_C} -rdynamic")
message("VISIT CXX Flags = ${VISIT_CXX_FLAGS}")
ExternalProject_Add (
  udaReaderMTMD 
  PREFIX ${VISIT_PLUGIN_BIN_DIR}
  DOWNLOAD_DIR ${VISIT_PLUGIN_SRC_DIR}
  SOURCE_DIR ${VISIT_PLUGIN_BIN_DIR}
  BINARY_DIR ${VISIT_PLUGIN_BIN_DIR}
  INSTALL_DIR ${VISIT_PLUGIN_BIN_DIR}
  INSTALL_COMMAND ""
  CMAKE_ARGS ${VISIT_CXX_FLAGS} ${VISIT_INSTALL_FLAGS}
)

