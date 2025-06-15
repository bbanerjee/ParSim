# cmake/ConfigureHypre.cmake

include(ExternalProject)

function(configure_hypre)
    # Set hypre configuration
    set(HYPRE_VERSION "2.33.0")
    set(HYPRE_TARGET_NAME "hypre")

    # Define directories
    set(HYPRE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/hypre/src)
    set(HYPRE_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/hypre-install)

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
    set(HYPRE_CONFIGURE_OPTIONS
        --prefix=${HYPRE_INSTALL_DIR}
        --enable-shared
        --enable-bigint=no
        --with-MPI
        --enable-fortran=no
    )

    # Convert list to space-separated string for shell script
    string(REPLACE ";" " " HYPRE_CONFIGURE_ARGS "${HYPRE_CONFIGURE_OPTIONS}")
    
    # Configure hypre as external project
    message(STATUS "HYPRE_TARGET_NAME: ${HYPRE_TARGET_NAME}")
    message(STATUS "HYPRE_SOURCE_DIR: ${HYPRE_SOURCE_DIR}")
    message(STATUS "HYPRE_CONFIGURE_OPTIONS: ${HYPRE_CONFIGURE_ARGS}")
    message(STATUS "HYPRE_INSTALL_DIR: ${HYPRE_INSTALL_DIR}")
    message(STATUS "MAKE_EXECUTABLE: ${MAKE_EXECUTABLE}")
    message(STATUS "CMAKE_BUILD_PARALLEL_LEVEL: ${CMAKE_BUILD_PARALLEL_LEVEL}")
    ExternalProject_Add(${HYPRE_TARGET_NAME}
        SOURCE_DIR ${HYPRE_SOURCE_DIR}
        CONFIGURE_COMMAND 
            ${HYPRE_SOURCE_DIR}/configure ${HYPRE_CONFIGURE_ARGS}
        BUILD_COMMAND 
            ${MAKE_EXECUTABLE} -j${CMAKE_BUILD_PARALLEL_LEVEL}
        INSTALL_COMMAND 
            ${MAKE_EXECUTABLE} install
        BUILD_IN_SOURCE 1
        INSTALL_DIR ${HYPRE_INSTALL_DIR}
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
        LOG_OUTPUT_ON_FAILURE ON
    )

    # Create imported target for hypre
    add_library(hypre::hypre STATIC IMPORTED)
    set_target_properties(hypre::hypre PROPERTIES
        IMPORTED_LOCATION ${HYPRE_INSTALL_DIR}/lib/libHYPRE.so
        INTERFACE_INCLUDE_DIRECTORIES ${HYPRE_INSTALL_DIR}/include
    )

    # Make sure the imported target depends on the external project
    add_dependencies(hypre::hypre ${HYPRE_TARGET_NAME})

    # Export hypre variables to parent scope
    set(HYPRE_INCLUDE_DIR ${HYPRE_INSTALL_DIR}/include PARENT_SCOPE)
    set(HYPRE_LIBRARY ${HYPRE_INSTALL_DIR}/lib/libHYPRE.so PARENT_SCOPE)
    set(HYPRE_LIBRARIES ${HYPRE_INSTALL_DIR}/lib/libHYPRE.so PARENT_SCOPE)
    set(HYPRE_FOUND TRUE PARENT_SCOPE)
endfunction()
