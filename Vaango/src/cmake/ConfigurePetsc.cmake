# ConfigurePETSc.cmake
# CMake module for configuring, building, and installing PETSc as an external project
#
# Usage:
#   include(cmake/ConfigurePETSc.cmake)
#   configure_petsc(
#     SOURCE_DIR path/to/petsc/source
#     INSTALL_DIR path/to/install/location
#     [DEBUGGING ON|OFF]
#     [SHARED_LIBS ON|OFF]
#     [WITH_MPI ON|OFF]
#     [DOWNLOAD_DEPS ON|OFF]
#     [EXTRA_CONFIGURE_ARGS arg1 arg2 ...]
#   )
#
# Creates targets:
#   - petsc: ExternalProject target
#   - PETSc::petsc: Imported library target
#
# Creates function:
#   - target_link_petsc(target_name): Link PETSc to your target

include(ExternalProject)

# Default PETSc configuration options
set(PETSC_DEFAULT_DEBUGGING OFF)
set(PETSC_DEFAULT_SHARED_LIBS ON)
set(PETSC_DEFAULT_WITH_MPI ON)
set(PETSC_DEFAULT_DOWNLOAD_DEPS ON)

function(set_petsc_arch)
    # Default fallback
    set(PETSC_ARCH_DEFAULT "arch-unknown")
    
    # Get OS name
    string(TOLOWER ${CMAKE_SYSTEM_NAME} OS_NAME)
    
    # Get architecture
    string(TOLOWER ${CMAKE_SYSTEM_PROCESSOR} ARCH_NAME)
    if(ARCH_NAME MATCHES "x86_64|amd64")
        set(ARCH_NAME "x86_64")
    elseif(ARCH_NAME MATCHES "aarch64|arm64")
        set(ARCH_NAME "arm64")
    endif()
    
    # Get compiler
    string(TOLOWER ${CMAKE_CXX_COMPILER_ID} COMPILER_NAME)
    if(COMPILER_NAME STREQUAL "gnu")
        set(COMPILER_NAME "gcc")
    elseif(COMPILER_NAME STREQUAL "appleclang")
        set(COMPILER_NAME "clang")
    endif()
    
    # Get build type
    string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    #if(BUILD_TYPE STREQUAL "debug")
    if(PETSC_DEBUGGING)
        set(BUILD_SUFFIX "debug")
    else()
        set(BUILD_SUFFIX "opt")
    endif()
    
    # Construct PETSC_ARCH
    set(PETSC_ARCH "${OS_NAME}-${COMPILER_NAME}-${BUILD_SUFFIX}" PARENT_SCOPE)
endfunction()


# Function to configure PETSc
function(configure_petsc)
    # Parse arguments
    set(options)
    set(oneValueArgs SOURCE_DIR INSTALL_DIR DEBUGGING SHARED_LIBS WITH_MPI DOWNLOAD_DEPS TARGET_NAME)
    set(multiValueArgs EXTRA_CONFIGURE_ARGS)
    
    cmake_parse_arguments(PETSC "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Validate required arguments
    if(NOT PETSC_SOURCE_DIR)
        message(FATAL_ERROR "configure_petsc: SOURCE_DIR is required")
    endif()
    
    if(NOT PETSC_INSTALL_DIR)
        message(FATAL_ERROR "configure_petsc: INSTALL_DIR is required")
    endif()
    
    # Set defaults
    if(NOT DEFINED PETSC_DEBUGGING)
        set(PETSC_DEBUGGING ${PETSC_DEFAULT_DEBUGGING})
    endif()
    
    if(NOT DEFINED PETSC_SHARED_LIBS)
        set(PETSC_SHARED_LIBS ${PETSC_DEFAULT_SHARED_LIBS})
    endif()
    
    if(NOT DEFINED PETSC_WITH_MPI)
        set(PETSC_WITH_MPI ${PETSC_DEFAULT_WITH_MPI})
    endif()
    
    if(NOT DEFINED PETSC_DOWNLOAD_DEPS)
        set(PETSC_DOWNLOAD_DEPS ${PETSC_DEFAULT_DOWNLOAD_DEPS})
    endif()
    
    if(NOT PETSC_TARGET_NAME)
        set(PETSC_TARGET_NAME "petsc_external")
    endif()

    # Check if source directory exists
    if(NOT EXISTS ${PETSC_SOURCE_DIR})
         message(FATAL_ERROR "PETSc source directory does not exist: ${PETSC_SOURCE_DIR}")
    endif()
    
    # Check PETSc version if possible
    if(EXISTS ${PETSC_SOURCE_DIR}/include/petscversion.h)
        file(READ ${PETSC_SOURCE_DIR}/include/petscversion.h PETSC_VERSION_FILE)
        string(REGEX MATCH "#define PETSC_VERSION_MAJOR[ ]+([0-9]+)" _ ${PETSC_VERSION_FILE})
        set(PETSC_VERSION_MAJOR ${CMAKE_MATCH_1})
        string(REGEX MATCH "#define PETSC_VERSION_MINOR[ ]+([0-9]+)" _ ${PETSC_VERSION_FILE})
        set(PETSC_VERSION_MINOR ${CMAKE_MATCH_1})
        string(REGEX MATCH "#define PETSC_VERSION_SUBMINOR[ ]+([0-9]+)" _ ${PETSC_VERSION_FILE})
        set(PETSC_VERSION_SUBMINOR ${CMAKE_MATCH_1})
        message(STATUS "PETSc: Detected version ${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")
    endif()
    
    # Build configure command
    set(PETSC_CONFIGURE_CMD "")
    #list(APPEND PETSC_CONFIGURE_CMD --prefix=${PETSC_INSTALL_DIR})
    
    # Debugging options
    if(PETSC_DEBUGGING)
        list(APPEND PETSC_CONFIGURE_CMD --with-debugging=1)
        message(STATUS "PETSc: Enabling debugging")
    else()
        list(APPEND PETSC_CONFIGURE_CMD --with-debugging=0)
        list(APPEND PETSC_CONFIGURE_CMD COPTFLAGS=-O3)
        list(APPEND PETSC_CONFIGURE_CMD CXXOPTFLAGS=-O3)
        list(APPEND PETSC_CONFIGURE_CMD FOPTFLAGS=-O3)
        message(STATUS "PETSc: Optimized build (no debugging)")
    endif()
    
    # Shared libraries
    if(PETSC_SHARED_LIBS)
        list(APPEND PETSC_CONFIGURE_CMD --with-shared-libraries=1)
        message(STATUS "PETSc: Building shared libraries")
    else()
        list(APPEND PETSC_CONFIGURE_CMD --with-shared-libraries=0)
        message(STATUS "PETSc: Building static libraries")
    endif()
    
    # MPI configuration
    if(PETSC_WITH_MPI)
        list(APPEND PETSC_CONFIGURE_CMD --with-mpi=1)
        message(STATUS "PETSc: Enabling MPI support")
        
        # Try to find system MPI
        find_package(MPI QUIET)
        if(MPI_FOUND)
            message(STATUS "PETSc: Found system MPI")
            if(MPI_C_COMPILER)
                set(ENV{MPICC} ${MPI_C_COMPILER})
            endif()
            if(MPI_CXX_COMPILER)
                set(ENV{MPICXX} ${MPI_CXX_COMPILER})
            endif()
            if(MPI_Fortran_COMPILER)
                set(ENV{MPIFC} ${MPI_Fortran_COMPILER})
            endif()
        else()
            message(STATUS "PETSc: System MPI not found, PETSc will use detected MPI")
        endif()
    else()
        list(APPEND PETSC_CONFIGURE_CMD --with-mpi=0)
        message(STATUS "PETSc: Disabling MPI support")
    endif()

    
    # Download dependencies
    if(PETSC_DOWNLOAD_DEPS)
        #list(APPEND PETSC_CONFIGURE_CMD --download-openblas=yes)
        list(APPEND PETSC_CONFIGURE_CMD --download-hypre=yes)
        #list(APPEND PETSC_CONFIGURE_CMD --download-fblaslapack=1)
        #list(APPEND PETSC_CONFIGURE_CMD --download-cmake=1)
        message(STATUS "PETSc: Will download dependencies (BLAS/LAPACK, Hypre)")
    endif()
    
    # Add extra configure arguments
    if(PETSC_EXTRA_CONFIGURE_ARGS)
        list(APPEND PETSC_CONFIGURE_CMD ${PETSC_EXTRA_CONFIGURE_ARGS})
        message(STATUS "PETSc: Extra configure args: ${PETSC_EXTRA_CONFIGURE_ARGS}")
    endif()
    
    # Set compiler environment variables
    if(CMAKE_C_COMPILER)
        set(ENV{CC} ${CMAKE_C_COMPILER})
    endif()
    if(CMAKE_CXX_COMPILER)
        set(ENV{CXX} ${CMAKE_CXX_COMPILER})
    endif()
    if(CMAKE_Fortran_COMPILER)
        set(ENV{FC} ${CMAKE_Fortran_COMPILER})
    endif()
    
    # Create the external project
    message(STATUS "PETSc: Configuring external project")
    message(STATUS "PETSc: Source dir: ${PETSC_SOURCE_DIR}")
    message(STATUS "PETSc: Install dir: ${PETSC_INSTALL_DIR}")

    # Convert list to space-separated string for shell script
    string(REPLACE ";" " " PETSC_CONFIGURE_ARGS "${PETSC_CONFIGURE_CMD}")

    # Set PETSC_ARCH
    set_petsc_arch()
    message(STATUS "PETSC_ARCH set to: ${PETSC_ARCH}")

    # Set PETSC_DIR
    set(PETSC_DIR "${PETSC_SOURCE_DIR}")

    # For cross-platform compatibility, create configure wrapper
    if(WIN32)
        set(PETSC_CONFIGURE_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/configure_petsc.bat)
        file(WRITE ${PETSC_CONFIGURE_SCRIPT}
            "@echo off\n"
            "cd /d ${PETSC_SOURCE_DIR}\n"
            "python configure ${PETSC_CONFIGURE_ARGS}\n"
        )
    else()
        set(PETSC_CONFIGURE_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/configure_petsc.sh)
        file(WRITE ${PETSC_CONFIGURE_SCRIPT}
            "#!/bin/bash\n"
            "cd ${PETSC_SOURCE_DIR}\n"
            "./configure PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} ${PETSC_CONFIGURE_ARGS}\n"
        )
        file(CHMOD ${PETSC_CONFIGURE_SCRIPT} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
    endif()
    
    # Find the 'make' executable (crucial for Ninja)
    find_program(MAKE_EXECUTABLE NAMES make REQUIRED)
    if(NOT MAKE_EXECUTABLE)
        message(FATAL_ERROR "Could not find 'make' executable. It is required for building PETSc.")
    endif()

    message(STATUS "PETSC_CONFIGURE_ARGS = ${PETSC_CONFIGURE_ARGS}")
    ExternalProject_Add(${PETSC_TARGET_NAME}
        SOURCE_DIR ${PETSC_SOURCE_DIR}
        CONFIGURE_COMMAND 
            ${PETSC_SOURCE_DIR}/configure PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} 
            --prefix=${PETSC_INSTALL_DIR}
            "${PETSC_CONFIGURE_ARGS}"
        BUILD_COMMAND 
            ${MAKE_EXECUTABLE} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
        INSTALL_COMMAND 
            ${MAKE_EXECUTABLE} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} install
        BUILD_IN_SOURCE 1
        INSTALL_DIR ${PETSC_INSTALL_DIR}
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
        BUILD_BYPRODUCTS
            "${PETSC_INSTALL_DIR}/lib/libpetsc.so" # Linux/Unix shared library
            "${PETSC_INSTALL_DIR}/include/petscsys.h"
            "${PETSC_INSTALL_DIR}/lib/pkgconfig/PETSc.pc"
    )
    
    # Determine library extension
    if(PETSC_SHARED_LIBS)
        set(PETSC_LIB_EXT ${CMAKE_SHARED_LIBRARY_SUFFIX})
        set(PETSC_LIB_TYPE SHARED)
    else()
        set(PETSC_LIB_EXT ${CMAKE_STATIC_LIBRARY_SUFFIX})
        set(PETSC_LIB_TYPE STATIC)
    endif()
    
    # Create imported target
    set(PETSC_INTERFACE_TARGET "PETSc_lib")
    # Create an INTERFACE library
    add_library(${PETSC_INTERFACE_TARGET} INTERFACE)

    # Set properties for the INTERFACE library
    target_include_directories(${PETSC_INTERFACE_TARGET} INTERFACE ${PETSC_INSTALL_DIR}/include)
    target_link_directories(${PETSC_INTERFACE_TARGET} INTERFACE ${PETSC_INSTALL_DIR}/lib)
    target_link_libraries(${PETSC_INTERFACE_TARGET} INTERFACE petsc) # Link against the actual library name, not the target name

    add_dependencies(${PETSC_INTERFACE_TARGET} ${PETSC_TARGET_NAME})
    
    message(STATUS "PETSc: Created interface target ${PETSC_INTERFACE_TARGET}")
    
    # Export variables to parent scope
    set(PETSC_FOUND TRUE PARENT_SCOPE)
    set(PETSC_INCLUDE_DIRS ${PETSC_INSTALL_DIR}/include PARENT_SCOPE)
    set(PETSC_LIBRARIES ${PETSC_INSTALL_DIR}/lib/libpetsc${PETSC_LIB_EXT} PARENT_SCOPE)
    set(PETSC_DIR ${PETSC_SOURCE_DIR} PARENT_SCOPE)
    set(PETSC_ARCH ${PETSC_ARCH} PARENT_SCOPE)
    set(PETSC_TARGET ${PETSC_INTERFACE_TARGET} PARENT_SCOPE)
    
    # Create utility targets
    add_custom_target(clean-petsc
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/${PETSC_TARGET_NAME}-prefix
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${PETSC_INSTALL_DIR}
        COMMENT "Cleaning PETSc build and installation"
    )
    
    add_custom_target(check-petsc
        COMMAND ${PETSC_INSTALL_DIR}/bin/petsc-config --version
        DEPENDS ${PETSC_TARGET_NAME}
        COMMENT "Checking PETSc installation"
    )
    
    message(STATUS "PETSc: Configuration complete")
endfunction()

# Convenience function to link PETSc to a target
function(target_link_petsc target_name access_type)
    if(NOT TARGET ${target_name})
        message(FATAL_ERROR "target_link_petsc: Target ${target_name} does not exist")
    endif()
    
    if(NOT TARGET PETSc_lib)
        message(FATAL_ERROR "target_link_petsc: PETSc has not been configured. Call configure_petsc() first.")
    endif()
    
    add_dependencies(${target_name} petsc_external)
    target_link_libraries(${target_name} ${access_type} ${PETSC_TARGET})
    
    message(STATUS "PETSc: Linked to target ${target_name}")
endfunction()

# Function to get PETSc configuration info
function(petsc_get_info)
    if(TARGET ${PETSC_TARGET})
        get_target_property(PETSC_LOCATION ${PETSC_TARGET} IMPORTED_LOCATION)
        get_target_property(PETSC_INCLUDES ${PETSC_TARGET} INTERFACE_INCLUDE_DIRECTORIES)
        
        message(STATUS "PETSc Library: ${PETSC_LOCATION}")
        message(STATUS "PETSc Includes: ${PETSC_INCLUDES}")
        message(STATUS "PETSc Dir: ${PETSC_DIR}")
    else()
        message(STATUS "PETSc has not been configured")
    endif()
endfunction()
