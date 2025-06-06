#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2014-2025 Biswajit Banerjee, Parresia Research Limited, NZ
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

### Modern usage
## add_library(CCA_Components_Arches source1.cpp source2.cpp)
## add_fortran_dependency(
##     FORTRAN_FILE fortran/bcscalar_fort.F
##     TARGET CCA_Components_Arches
## )

### With custom FSPEC tool
## add_fortran_dependency(
##     FORTRAN_FILE fortran/another_fort.F
##     TARGET CCA_Components_Arches
##     FSPEC_TOOL /custom/path/to/fspec
## )

### Legacy compatibility (deprecated but still works)
## FORTRAN_DEPENDENCY(fortran/bcscalar_fort.F CCA_Components_Arches)

#[[
Modern Fortran dependency function that generates header files from .fspec files
and manages Fortran source lists for targets.

Usage:
  add_fortran_dependency(
    FORTRAN_FILE fortran/bcscalar_fort.F
    TARGET my_target
    [FSPEC_TOOL /path/to/fspec]
  )

This replaces the old FORTRAN_DEPENDENCY macro with better practices:
- Uses functions instead of macros for better scoping
- Uses target-based approach instead of global variables
- Proper dependency tracking with custom commands
- Better error handling and validation
]]

function(add_fortran_dependency)
    # Parse arguments
    set(options "")
    set(oneValueArgs FORTRAN_FILE TARGET FSPEC_TOOL)
    set(multiValueArgs "")
    
    cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Validate required arguments
    if(NOT ARG_FORTRAN_FILE)
        message(FATAL_ERROR "add_fortran_dependency: FORTRAN_FILE is required")
    endif()
    
    if(NOT ARG_TARGET)
        message(FATAL_ERROR "add_fortran_dependency: TARGET is required")
    endif()
    
    # Set default FSPEC tool if not provided
    if(NOT ARG_FSPEC_TOOL)
        if(DEFINED FSPEC)
            set(ARG_FSPEC_TOOL ${FSPEC})
        else()
            find_program(PERL_EXECUTABLE perl REQUIRED)
            find_file(FSPEC_SCRIPT fspec 
                PATHS ${CMAKE_SOURCE_DIR}/scripts ${CMAKE_SOURCE_DIR}/tools
                DOC "FSPEC tool for generating Fortran headers")
            if(NOT FSPEC_SCRIPT)
                message(FATAL_ERROR "add_fortran_dependency: FSPEC tool not found. Please set FSPEC_TOOL argument.")
            endif()
            set(ARG_FSPEC_TOOL ${PERL_EXECUTABLE} ${FSPEC_SCRIPT})
        endif()
    endif()
    
    # Generate file paths
    get_filename_component(fortran_name ${ARG_FORTRAN_FILE} NAME)
    string(REPLACE ".F" "_fort.h" header_name ${fortran_name})
    string(REPLACE ".F" ".fspec" fspec_name ${fortran_name})
    
    set(source_fspec "${CMAKE_CURRENT_SOURCE_DIR}/${ARG_FORTRAN_FILE}")
    string(REPLACE ".F" ".fspec" source_fspec ${source_fspec})
    
    set(output_header "${CMAKE_CURRENT_BINARY_DIR}/fortran/${header_name}")
    
    # Ensure output directory exists
    get_filename_component(output_dir ${output_header} DIRECTORY)
    file(MAKE_DIRECTORY ${output_dir})
    
    # Validate input files exist
    if(NOT EXISTS ${source_fspec})
        message(FATAL_ERROR "add_fortran_dependency: fspec file does not exist: ${source_fspec}")
    endif()
    
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${ARG_FORTRAN_FILE}")
        message(FATAL_ERROR "add_fortran_dependency: Fortran file does not exist: ${ARG_FORTRAN_FILE}")
    endif()
    
    # Create custom command to generate header from fspec
    add_custom_command(
        OUTPUT ${output_header}
        COMMAND ${ARG_FSPEC_TOOL} ${source_fspec} ${output_header}
        DEPENDS ${source_fspec}
        COMMENT "Generating Fortran header: ${header_name}"
        VERBATIM
    )
    
    # Create a custom target for the generated header (optional, for visibility)
    get_filename_component(header_target_name ${header_name} NAME_WE)
    set(header_target_name "${ARG_TARGET}_${header_target_name}_header")
    
    add_custom_target(${header_target_name}
        DEPENDS ${output_header}
    )
    
    # Add sources to the target
    target_sources(${ARG_TARGET} PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_FORTRAN_FILE}
        ${output_header}  # Include generated header as source for dependency tracking
    )
    
    # Ensure the target depends on the header generation
    add_dependencies(${ARG_TARGET} ${header_target_name})
    
    # Add binary directory to include path so generated headers can be found
    target_include_directories(${ARG_TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    
    # Store information for debugging/status
    set_property(TARGET ${ARG_TARGET} APPEND PROPERTY FORTRAN_GENERATED_HEADERS ${output_header})
    
    message(STATUS "Added Fortran dependency: ${ARG_FORTRAN_FILE} -> ${header_name} for target ${ARG_TARGET}")
    
endfunction()

#[[
Legacy compatibility macro (deprecated - use add_fortran_dependency function instead)
]]
macro(FORTRAN_DEPENDENCY fortran_file lib)
    message(DEPRECATION "FORTRAN_DEPENDENCY macro is deprecated. Use add_fortran_dependency function instead.")
    
    add_fortran_dependency(
        FORTRAN_FILE ${fortran_file}
        TARGET ${lib}
    )
endmacro()

#[[
Utility function to get all Fortran generated headers for a target
]]
function(get_target_fortran_headers target_name output_var)
    get_property(headers TARGET ${target_name} PROPERTY FORTRAN_GENERATED_HEADERS)
    set(${output_var} ${headers} PARENT_SCOPE)
endfunction()


