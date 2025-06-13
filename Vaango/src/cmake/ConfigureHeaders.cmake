# cmake/ConfigureHeaders.cmake
# Modern CMake 3.28 function for configuring header files

#[=======================================================================[.rst:
ConfigureHeaders
----------------

This module provides functions for configuring header files from templates.

.. command:: configure_sci_headers

  Configure multiple science definition header files from templates::

    configure_sci_headers(
        SOURCE_DIR <source_directory>
        BINARY_DIR <binary_directory>
        [INCLUDE_DIRECTORIES <dir1> <dir2> ...]
        [HEADERS <header1> <header2> ...]
    )

  ``SOURCE_DIR``
    Source directory containing the .in template files
  
  ``BINARY_DIR``
    Binary directory where configured files will be generated
  
  ``INCLUDE_DIRECTORIES``
    Optional list of include directories to add
  
  ``HEADERS``
    Optional list of specific headers to configure. If not provided,
    all standard sci_defs headers will be configured.

#]=======================================================================]

include_guard(GLOBAL)

# Default list of all sci_defs headers
set(_SCI_DEFS_HEADERS
    audio_defs
    babel_defs
    bits_defs
    blas_defs
    boost_defs
    chromium_defs
    collab_vis_defs
    compile_defs
    config_defs
    crypto_defs
    cuda_defs
    dataflow_defs
    dynamic_cast_defs
    environment_defs
    error_defs
    exe_defs
    framework_defs
    git_info
    globus_defs
    gperftools_defs
    hashmap_defs
    hdf5_defs
    hypre_defs
    ieeefp_defs
    image_defs
    kepler_defs
    kokkos_defs
    lapack_defs
    loki_defs
    malloc_defs
    mdsplus_defs
    mpi_defs
    osx_defs
    papi_defs
    petsc_defs
    pidx_defs
    ptolemy_defs
    ruby_defs
    scisock_defs
    ssl_defs
    stat64_defs
    template_defs
    tena_defs
    thread_defs
    uintah_defs
    vdt_defs
    visit_defs
    vtk_defs
    wx_defs
    z_defs
)

# Use a macro instead of function to inherit all parent scope variables automatically
macro(configure_sci_headers)
    # Parse arguments
    cmake_parse_arguments(
        CSH                          # prefix
        ""                          # options
        "SOURCE_DIR;BINARY_DIR"     # one_value_keywords
        "INCLUDE_DIRECTORIES;HEADERS" # multi_value_keywords
        ${ARGN}
    )
    
    # Validate required arguments
    if(NOT CSH_SOURCE_DIR)
        message(FATAL_ERROR "configure_sci_headers: SOURCE_DIR is required")
    endif()
    
    if(NOT CSH_BINARY_DIR)
        message(FATAL_ERROR "configure_sci_headers: BINARY_DIR is required")
    endif()
    
    # Convert to absolute paths
    if(NOT IS_ABSOLUTE "${CSH_SOURCE_DIR}")
        set(CSH_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${CSH_SOURCE_DIR}")
    endif()
    
    if(NOT IS_ABSOLUTE "${CSH_BINARY_DIR}")
        set(CSH_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CSH_BINARY_DIR}")
    endif()
    
    # Use provided headers or default list
    if(CSH_HEADERS)
        set(headers_to_configure ${CSH_HEADERS})
    else()
        set(headers_to_configure ${_SCI_DEFS_HEADERS})
    endif()
    
    # Add include directories if specified
    if(CSH_INCLUDE_DIRECTORIES)
        foreach(inc_dir ${CSH_INCLUDE_DIRECTORIES})
            if(IS_ABSOLUTE "${inc_dir}")
                target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE "${inc_dir}")
            else()
                target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE 
                    "${CMAKE_CURRENT_SOURCE_DIR}/${inc_dir}")
            endif()
        endforeach()
    endif()
    
    # Create output directory if it doesn't exist
    file(MAKE_DIRECTORY "${CSH_BINARY_DIR}/include/sci_defs")
    
    # Configure each header file
    foreach(header ${headers_to_configure})
        # Convert output name to input name pattern
        string(REPLACE "_defs" "_testdefs" input_name "${header}")
        
        # Determine paths
        if(header STREQUAL "git_info")
            set(input_file "${CSH_SOURCE_DIR}/include/sci_defs/${header}.h.in")
            set(output_file "${CSH_BINARY_DIR}/include/sci_defs/${header}.h")
        else()
            set(input_file "${CSH_SOURCE_DIR}/include/sci_defs/${input_name}.h.in")
            set(output_file "${CSH_BINARY_DIR}/include/sci_defs/${header}.h")
        endif()
        
        # Check if input file exists
        if(NOT EXISTS "${input_file}")
            message(WARNING "configure_sci_headers: Input file not found: ${input_file}")
            continue()
        endif()
        
        # Configure the file
        configure_file("${input_file}" "${output_file}" @ONLY)
        
        message(STATUS "Configured: ${header}.h")
    endforeach()
    
    # Add configured include directories to the project
    #target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE
    #    "${CSH_SOURCE_DIR}/include"
    #    "${CSH_BINARY_DIR}/include"
    #    "${CSH_BINARY_DIR}/include/sci_defs"
    #)
    
    # Get the count of configured headers
    list(LENGTH headers_to_configure header_count)
    message(STATUS "Configured ${header_count} header files")
endmacro()

# Convenience macro for the standard sci_defs configuration
macro(configure_standard_sci_headers)
    configure_sci_headers(
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
        BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}"
        #INCLUDE_DIRECTORIES 
        #    "include"
        #    "${CMAKE_CURRENT_BINARY_DIR}/include"
        #    "${CMAKE_CURRENT_BINARY_DIR}/include/sci_defs"
        ${ARGN}  # Pass through any additional arguments
    )
endmacro()