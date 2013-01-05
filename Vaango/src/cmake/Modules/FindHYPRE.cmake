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
find_library(HYPRE_LIB_DIST_MATRIX NAMES HYPRE_DistributedMatrix)
find_library(HYPRE_LIB_DIST_MATRIX_PILUT NAMES HYPRE_DistributedMatrixPilutSolver)
find_library(HYPRE_LIB_EUCLID NAMES HYPRE_Euclid)
find_library(HYPRE_LIB_FEI NAMES HYPRE_FEI)
find_library(HYPRE_LIB_FEI_FGMRES NAMES HYPRE_FEI_fgmres)
find_library(HYPRE_LIB_IJ_MV NAMES HYPRE_IJ_mv)
find_library(HYPRE_LIB_KRYLOV NAMES HYPRE_krylov)
find_library(HYPRE_LIB_LSI NAMES HYPRE_LSI)
find_library(HYPRE_LIB_MATRIX NAMES HYPRE_MatrixMatrix)
find_library(HYPRE_LIB_MLI NAMES HYPRE_mli)
find_library(HYPRE_LIB_MULTIVECTOR NAMES HYPRE_multivector)
find_library(HYPRE_LIB_PARASAILS NAMES HYPRE_ParaSails)
find_library(HYPRE_LIB_PARSCR_BLOCK NAMES HYPRE_parcsr_block_mv)
find_library(HYPRE_LIB_PARSCR_LS NAMES HYPRE_parcsr_ls)
find_library(HYPRE_LIB_PARSCR_MV NAMES HYPRE_parcsr_mv)
find_library(HYPRE_LIB_SEQ_MV NAMES HYPRE_seq_mv)
find_library(HYPRE_LIB_SSTRUCT_LS NAMES HYPRE_sstruct_ls)
find_library(HYPRE_LIB_SSTRUCT_MV NAMES HYPRE_sstruct_mv)
find_library(HYPRE_LIB_STRUCT_LS NAMES HYPRE_struct_ls)
find_library(HYPRE_LIB_STRUCT_MV NAMES HYPRE_struct_mv)
find_library(HYPRE_LIB_UTIL NAMES HYPRE_utilities)
set(HYPRE_LIBRARY
  ${HYPRE_LIBRARY}
  ${HYPRE_LIB_DIST_MATRIX}
  ${HYPRE_LIB_DIST_MATRIX_PILUT}
  ${HYPRE_LIB_EUCLID}
  ${HYPRE_LIB_FEI}
  ${HYPRE_LIB_FEI_FGMRES}
  ${HYPRE_LIB_IJ_MV}
  ${HYPRE_LIB_KRYLOV}
  ${HYPRE_LIB_LSI}
  ${HYPRE_LIB_MATRIX}
  ${HYPRE_LIB_MLI}
  ${HYPRE_LIB_MULTIVECTOR}
  ${HYPRE_LIB_PARASAILS}
  ${HYPRE_LIB_PARSCR_BLOCK}
  ${HYPRE_LIB_PARSCR_LS}
  ${HYPRE_LIB_PARSCR_MV}
  ${HYPRE_LIB_SEQ_MV}
  ${HYPRE_LIB_SSTRUCT_LS}
  ${HYPRE_LIB_SSTRUCT_MV}
  ${HYPRE_LIB_STRUCT_LS}
  ${HYPRE_LIB_STRUCT_MV}
  ${HYPRE_LIB_UTIL}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG HYPRE_LIBRARY HYPRE_INCLUDE)

mark_as_advanced(HYPRE_INCLUDE HYPRE_LIBRARY)



