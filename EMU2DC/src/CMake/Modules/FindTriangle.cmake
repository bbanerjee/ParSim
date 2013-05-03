# Find where the Triangle libraries and headers are located
# These are needed if you wish to access the Delaunay triangulation machinery

include(CheckFunctionExists)

# Triangle includes
set(INCLUDEPATH ${CMAKE_SOURCE_DIR}/../Triangle/include)
message(STATUS ${INCLUDEPATH})

find_path(TRIANGLE_INCLUDE_DIR 
          NAMES del_interface.hpp
          PATHS "${INCLUDEPATH}")

message(STATUS "${TRIANGLE_INCLUDE_DIR}")

if (TRIANGLE_INCLUDE_DIR)
  message(STATUS "tpp::triangle.h found: INCLUDE Path = ${TRIANGLE_INCLUDE_DIR}")
else()
  message(FATAL_ERROR "tpp::triangle.h not found: INCLUDE Path = ${TRIANGLE_INCLUDE_DIR}")
endif()

# Triangle library
set(LIBRARYPATH ${CMAKE_SOURCE_DIR}/../Triangle/opt)
message(STATUS ${LIBRARYPATH})
find_library(TRIANGLE_LIBRARY 
             NAMES Triangle++ 
             PATHS "${LIBRARYPATH}")

if (TRIANGLE_LIBRARY)
  set(TRIANGLE_FOUND)
  message(STATUS "tpp::Triangle.so found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${TRIANGLE_LIBRARY}")
else()
  set(TRIANGLE_LIBRARY "")
  message(FATAL_ERROR "tpp::Triangle.so not found: INCLUDE Path = ${LIBRARYPATH} LIBRARY = ${TRIANGLE_LIBRARY}")
endif()


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Triangle++ DEFAULT_MSG TRIANGLE_LIBRARY TRIANGLE_INCLUDE_DIR)

mark_as_advanced(TRIANGLE_INCLUDE_DIR TRIANGLE_LIBRARY)

