# Additional target to perform clang-format/clang-tidy run
# Requires clang-format and clang-tidy

# Get all project files
file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.hpp *.cc *.h)

#message(STATUS "Binary dir = ${CMAKE_BINARY_DIR}")
#message(STATUS "Include dir = ${BASEPATH}")
#message(STATUS "MPI lib =   ${MPI_INCLUDE_PATH_CONV}")
#message(STATUS "Boost lib =   ${Boost_INCLUDE_DIR}")
#message(STATUS "VTK lib =   ${VTK_INCLUDE_DIRS}")

add_custom_target(
        clang-format
        COMMAND /usr/bin/clang-format
        -style=Mozilla
        -i
        ${ALL_SOURCE_FILES}
)

add_custom_target(
        clang-tidy
        COMMAND /usr/bin/clang-tidy
        ${ALL_SOURCE_FILES}
        -config=''
        -export-fixes='clang-tidy-fixes.dat'
        --
        -std=c++17
        -I${BASEPATH}
        -I${MPI_INCLUDE_PATH_CONV}
        -I${Boost_INCLUDE_DIR}
        -I${VTK_INCLUDE_DIRS}
)
