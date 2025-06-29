# cmake/ConfigureHypre.cmake

include(ExternalProject)

function(configure_hypre)
    # Set hypre configuration
    set(HYPRE_VERSION "2.33.0")
    set(HYPRE_TARGET_NAME "hypre_external")
    set(HYPRE_INTERFACE_TARGET_NAME "hypre_lib") 

    # Define directories
    set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/hypre/src)
    set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/hypre-install)

    file(MAKE_DIRECTORY "${HYPRE_INSTALL_DIR}")

    # Check if hypre submodule exists
    if(NOT EXISTS ${HYPRE_SOURCE_DIR}/CMakeLists.txt AND NOT EXISTS ${HYPRE_SOURCE_DIR}/configure)
        message(FATAL_ERROR 
            "Hypre submodule not found at ${HYPRE_SOURCE_DIR}. "
            "Please run: git submodule update --init --recursive")
    endif()

    # Find required tools
    find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)
    if(NOT MAKE_EXECUTABLE)
        message(FATAL_ERROR "Could not find 'make' executable required for building hypre.")
    endif()

    # Set hypre configure options
    set(HYPRE_CONFIGURE_CMD
        ${HYPRE_SOURCE_DIR}/configure
        --prefix=${HYPRE_INSTALL_DIR}
        --enable-shared
        --enable-bigint=no
        --with-MPI
        --enable-fortran=no
    )

    # Convert list to space-separated string for shell script
    string(REPLACE ";" " " HYPRE_CONFIGURE_ARGS "${HYPRE_CONFIGURE_CMD}")

    # Set up configure script
    set(HYPRE_CONFIGURE_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/configure_hypre.sh)
    file(WRITE ${HYPRE_CONFIGURE_SCRIPT}
        "#!/bin/bash\n"
        "cd ${HYPRE_SOURCE_DIR}\n"
        "${HYPRE_CONFIGURE_ARGS}\n"
    )
    file(CHMOD ${HYPRE_CONFIGURE_SCRIPT} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
    
    # Configure hypre as external project
    message(STATUS "HYPRE_TARGET_NAME: ${HYPRE_TARGET_NAME}")
    message(STATUS "HYPRE_SOURCE_DIR: ${HYPRE_SOURCE_DIR}")
    message(STATUS "HYPRE_CONFIGURE_OPTIONS: ${HYPRE_CONFIGURE_ARGS}")
    message(STATUS "HYPRE_INSTALL_DIR: ${HYPRE_INSTALL_DIR}")
    message(STATUS "MAKE_EXECUTABLE: ${MAKE_EXECUTABLE}")
    message(STATUS "CMAKE_BUILD_PARALLEL_LEVEL: ${CMAKE_BUILD_PARALLEL_LEVEL}")
    ExternalProject_Add(${HYPRE_TARGET_NAME}
        SOURCE_DIR ${HYPRE_SOURCE_DIR}
        CONFIGURE_COMMAND ${HYPRE_CONFIGURE_SCRIPT}
        BUILD_COMMAND ${MAKE_EXECUTABLE} -j${CMAKE_BUILD_PARALLEL_LEVEL}
        INSTALL_COMMAND ${MAKE_EXECUTABLE} install
        BUILD_IN_SOURCE 1
        INSTALL_DIR ${HYPRE_INSTALL_DIR}
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
        LOG_OUTPUT_ON_FAILURE ON
        BUILD_BYPRODUCTS
            "${HYPRE_INSTALL_DIR}/lib/libHYPRE.so"
            "${HYPRE_INSTALL_DIR}/include/HYPRE.h"
    )

    # Create imported target for hypre
    add_library(${HYPRE_INTERFACE_TARGET_NAME} INTERFACE)
    target_include_directories(${HYPRE_INTERFACE_TARGET_NAME} INTERFACE ${HYPRE_INSTALL_DIR}/include)
    target_link_directories(${HYPRE_INTERFACE_TARGET_NAME} INTERFACE ${HYPRE_INSTALL_DIR}/lib) # Add link directories
    target_link_libraries(${HYPRE_INTERFACE_TARGET_NAME} INTERFACE HYPRE) # Link against the library name

    # Make sure the imported target depends on the external project
    add_dependencies(${HYPRE_INTERFACE_TARGET_NAME} ${HYPRE_TARGET_NAME})

    # Export hypre variables to parent scope
    set(HYPRE_INCLUDE_DIR ${HYPRE_INSTALL_DIR}/include PARENT_SCOPE)
    set(HYPRE_LIBRARY ${HYPRE_INSTALL_DIR}/lib/libHYPRE.so PARENT_SCOPE)
    set(HYPRE_LIBRARIES ${HYPRE_INSTALL_DIR}/lib/libHYPRE.so PARENT_SCOPE)
    set(HYPRE_FOUND TRUE PARENT_SCOPE)
    set(HYPRE_TARGET ${HYPRE_INTERFACE_TARGET_NAME} PARENT_SCOPE)

    execute_process(COMMAND git -C ${HYPRE_SOURCE_DIR} describe --match v* --abbrev=0
                    OUTPUT_VARIABLE HYPRE_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE git_result)
    message(STATUS "Hyper git tag: ${HYPRE_VERSION}")
    string(REGEX MATCH "^v?([0-9]+)\\.([0-9]+)\\.([0-9]+)$" _ "${HYPRE_VERSION}")
    set(HYPRE_VERSION_MAJOR "${CMAKE_MATCH_1}" PARENT_SCOPE)
    set(HYPRE_VERSION_MINOR "${CMAKE_MATCH_2}" PARENT_SCOPE)
    set(HYPRE_VERSION_PATCH "${CMAKE_MATCH_3}" PARENT_SCOPE)

    # Print the variables to verify
    #message(STATUS "Original String: ${develop_string}")
    #message(STATUS "Inside fn: HYPRE_MAJOR_VERSION: ${HYPRE_VERSION_MAJOR}")
    #message(STATUS "Inside fn: HYPRE_MINOR_VERSION: ${HYPRE_VERSION_MINOR}")
    #message(STATUS "Inside fn: HYPRE_VERSION_PATCH: ${HYPRE_VERSION_PATCH}")
endfunction()
