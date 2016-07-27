# Find where the Uintah ProblemSpec libraries and headers are located
# These are needed if you wish to access the Uintah machinery for reading XML input files etc.

include(CheckFunctionExists)

# ProblemSpec includes
set(INCLUDEPATH ${CMAKE_SOURCE_DIR}/../../Vaango/src)
message(STATUS ${INCLUDEPATH})

find_path(PROBLEMSPEC_INCLUDE_DIR 
          NAMES CCA
          PATHS "${INCLUDEPATH}")

message(STATUS "${PROBLEMSPEC_INCLUDE_DIR}")

if (PROBLEMSPEC_INCLUDE_DIR)
  message(STATUS "Uintah::ProblemSpec.h found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Uintah::ProblemSpec.h not found: INCLUDE Path = ${PROBLEMSPEC_INCLUDE_DIR}")
endif()

# SCI definitions 
set(SCIDEFS_PATH ${CMAKE_SOURCE_DIR}/../../Vaango/opt/include)
message(STATUS ${SCIDEFS_PATH})

find_path(SCIDEFS_INCLUDE_DIR 
          NAMES sci_defs
          PATHS "${SCIDEFS_PATH}")

message(STATUS "${SCIDEFS_INCLUDE_DIR}")

if (SCIDEFS_INCLUDE_DIR)
  message(STATUS "Uintah::sci_defs found: INCLUDE Path = ${SCIDEFS_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "Uintah::sci_defs not found: INCLUDE Path = ${SCIDEFS_INCLUDE_DIR}")
endif()

set(PROBLEMSPEC_INCLUDE_DIR
    ${PROBLEMSPEC_INCLUDE_DIR} 
    ${SCIDEFS_INCLUDE_DIR})

# ProblemSpec library
set(LIBRARYPATH ${CMAKE_SOURCE_DIR}/../../Vaango/opt/lib)
message(STATUS ${LIBRARYPATH})
find_library(PROBLEMSPEC_LIBRARY 
             NAMES Vaango_Core_ProblemSpec 
             PATHS "${LIBRARYPATH}")

if (PROBLEMSPEC_LIBRARY)
  set(PROBLEMSPEC_FOUND)
  message(STATUS "Uintah::ProblemSpec.so found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${PROBLEMSPEC_LIBRARY}")
else()
  set(PROBLEMSPEC_LIBRARY "")
  message(FATAL_ERROR "Uintah::ProblemSpec.so not found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${PROBLEMSPEC_LIBRARY}")
endif()

# Containers library
find_library(CONTAINERS_LIBRARY 
             NAMES Vaango_Core_Containers 
             PATHS "${LIBRARYPATH}")

if (CONTAINERS_LIBRARY)
  set(CONTAINERS_FOUND)
  message(STATUS "Uintah::Containers.so found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${CONTAINERS_LIBRARY}")
else()
  set(CONTAINERS_LIBRARY "")
  message(FATAL_ERROR "Uintah::Containers.so not found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${CONTAINERS_LIBRARY}")
endif()

# Util library
find_library(UTIL_LIBRARY 
             NAMES Vaango_Core_Util 
             PATHS "${LIBRARYPATH}")

if (UTIL_LIBRARY)
  set(UTIL_FOUND)
  message(STATUS "Uintah::Util.so found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${UTIL_LIBRARY}")
else()
  set(UTIL_LIBRARY "")
  message(FATAL_ERROR "Uintah::Util.so not found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${UTIL_LIBRARY}")
endif()

# Math library
find_library(MATH_LIBRARY 
             NAMES Vaango_Core_Math 
             PATHS "${LIBRARYPATH}")

if (MATH_LIBRARY)
  set(MATH_FOUND)
  message(STATUS "Uintah::Math.so found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${MATH_LIBRARY}")
else()
  set(MATH_LIBRARY "")
  message(FATAL_ERROR "Uintah::Math.so not found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${MATH_LIBRARY}")
endif()

#set(PROBLEMSPEC_LIBRARY
#    ${CONTAINERS_LIBRARY}
#    ${MATH_LIBRARY}
#    ${UTIL_LIBRARY}
#    ${PROBLEMSPEC_LIBRARY})

set(PROBLEMSPEC_LIBRARY
    ${UTIL_LIBRARY}
    ${MATH_LIBRARY}
    ${CONTAINERS_LIBRARY}
    ${PROBLEMSPEC_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Vaango_Core_ProblemSpec DEFAULT_MSG PROBLEMSPEC_LIBRARY PROBLEMSPEC_INCLUDE_DIR)

mark_as_advanced(PROBLEMSPEC_INCLUDE_DIR PROBLEMSPEC_LIBRARY)

