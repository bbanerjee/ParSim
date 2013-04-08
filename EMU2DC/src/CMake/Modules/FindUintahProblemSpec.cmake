# Find where the Uintah ProblemSpec libraries and headers are located
# These are needed if you wish to access the Uintah machinery for reading XML input files etc.

include(CheckFunctionExists)

set(INCLUDEPATH ${CMAKE_SOURCE_DIR}/../../Vaango/src)
message(STATUS ${INCLUDEPATH})
set(LIBRARYPATH ${CMAKE_SOURCE_DIR}/../../Vaango/opt/lib)
message(STATUS ${LIBRARYPATH})

find_path(PROBLEMSPEC_INCLUDE_DIR 
          NAMES CCA
          PATHS "${INCLUDEPATH}")

message(STATUS "${PROBLEMSPEC_INCLUDE_DIR}")

if (PROBLEMSPEC_INCLUDE_DIR)
  message(STATUS "Uintah::ProblemSpec.h found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Uintah::ProblemSpec.h not found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR}")
endif()


find_library(PROBLEMSPEC_LIBRARY 
             NAMES Vaango_Core_ProblemSpec 
             PATHS "${LIBRARYPATH}")

if (PROBLEMSPEC_LIBRARY)
  set(PROBLEMSPEC_FOUND)
  message(STATUS "Uintah::ProblemSpec.so found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR} LIBRARY = ${PROBLEMSPEC_LIBRARY}")
else()
  set(PROBLEMSPEC_LIBRARY "")
  message(FATAL_ERROR "Uintah::ProblemSpec.so not found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR} LIBRARY = ${PROBLEMSPEC_LIBRARY}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Vaango_Core_ProblemSpec DEFAULT_MSG PROBLEMSPEC_LIBRARY PROBLEMSPEC_INCLUDE_DIR)

mark_as_advanced(PROBLEMSPEC_INCLUDE_DIR PROBLEMSPEC_LIBRARY)

