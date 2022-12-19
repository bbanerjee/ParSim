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

# - find where HYPRE and friends are located.
# HYPRE_FOUND - system has dynamic linking interface available
# HYPRE_INCLUDE_DIR - where dlfcn.h is located.
# HYPRE_LIBRARIES - libraries needed to use dlopen

include(CheckFunctionExists)

find_path(HYPRE_INCLUDE NAMES HYPRE.h)
find_library(HYPRE_LIBRARY NAMES HYPRE)
if(HYPRE_LIBRARY)
  set(HYPRE_FOUND)
  message(STATUS "HYPRE found: INCLUDE Path = ${HYPRE_INCLUDE} LIBRARY Path: ${HYPRE_LIBRARY}")
else()
  set(HYPRE_LIBRARY "")
  message(FATAL_ERROR "HYPRE not found: INCLUDE Path = ${HYPRE_INCLUDE} LIBRARY Path: ${HYPRE_LIBRARY}")
endif()
#find_library(HYPRE_LIB_DIST_MATRIX NAMES HYPRE_DistributedMatrix)
#find_library(HYPRE_LIB_DIST_MATRIX_PILUT NAMES HYPRE_DistributedMatrixPilutSolver)
#find_library(HYPRE_LIB_EUCLID NAMES HYPRE_Euclid)
#find_library(HYPRE_LIB_FEI NAMES HYPRE_FEI)
#find_library(HYPRE_LIB_FEI_FGMRES NAMES HYPRE_FEI_fgmres)
#find_library(HYPRE_LIB_IJ_MV NAMES HYPRE_IJ_mv)
#find_library(HYPRE_LIB_KRYLOV NAMES HYPRE_krylov)
#find_library(HYPRE_LIB_LSI NAMES HYPRE_LSI)
#find_library(HYPRE_LIB_MATRIX NAMES HYPRE_MatrixMatrix)
#find_library(HYPRE_LIB_MLI NAMES HYPRE_mli)
#find_library(HYPRE_LIB_MULTIVECTOR NAMES HYPRE_multivector)
#find_library(HYPRE_LIB_PARASAILS NAMES HYPRE_ParaSails)
#find_library(HYPRE_LIB_PARSCR_BLOCK NAMES HYPRE_parcsr_block_mv)
#find_library(HYPRE_LIB_PARSCR_LS NAMES HYPRE_parcsr_ls)
#find_library(HYPRE_LIB_PARSCR_MV NAMES HYPRE_parcsr_mv)
#find_library(HYPRE_LIB_SEQ_MV NAMES HYPRE_seq_mv)
#find_library(HYPRE_LIB_SSTRUCT_LS NAMES HYPRE_sstruct_ls)
#find_library(HYPRE_LIB_SSTRUCT_MV NAMES HYPRE_sstruct_mv)
#find_library(HYPRE_LIB_STRUCT_LS NAMES HYPRE_struct_ls)
#find_library(HYPRE_LIB_STRUCT_MV NAMES HYPRE_struct_mv)
#find_library(HYPRE_LIB_UTIL NAMES HYPRE_utilities)
#set(HYPRE_LIBRARY
#  ${HYPRE_LIBRARY}
#  ${HYPRE_LIB_DIST_MATRIX}
#  ${HYPRE_LIB_DIST_MATRIX_PILUT}
#  ${HYPRE_LIB_EUCLID}
#  ${HYPRE_LIB_FEI}
#  ${HYPRE_LIB_FEI_FGMRES}
#  ${HYPRE_LIB_IJ_MV}
#  ${HYPRE_LIB_KRYLOV}
#  ${HYPRE_LIB_LSI}
#  ${HYPRE_LIB_MATRIX}
#  ${HYPRE_LIB_MLI}
#  ${HYPRE_LIB_MULTIVECTOR}
#  ${HYPRE_LIB_PARASAILS}
#  ${HYPRE_LIB_PARSCR_BLOCK}
#  ${HYPRE_LIB_PARSCR_LS}
#  ${HYPRE_LIB_PARSCR_MV}
#  ${HYPRE_LIB_SEQ_MV}
#  ${HYPRE_LIB_SSTRUCT_LS}
#  ${HYPRE_LIB_SSTRUCT_MV}
#  ${HYPRE_LIB_STRUCT_LS}
#  ${HYPRE_LIB_STRUCT_MV}
#  ${HYPRE_LIB_UTIL}
#)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG HYPRE_LIBRARY HYPRE_INCLUDE)

mark_as_advanced(HYPRE_INCLUDE HYPRE_LIBRARY)



