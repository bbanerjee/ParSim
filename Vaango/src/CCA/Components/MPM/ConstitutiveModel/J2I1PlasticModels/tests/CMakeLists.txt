#------------------------
# Libraries needed
#------------------------
SET(VAANGO_LIBS 
        Vaango_Core_Parallel    
        MPM_ConstitutiveModel
        MPM_ConstitutiveModel_FORTRAN
        MPM_ConstitutiveModel_FORTRAN_Utils
        ${HDF5_TARGET}
        $<$<BOOL:${HAVE_HYPRE}>:${HYPRE_TARGET}>
)
#------------------------
# Google tests
#------------------------

set(GTEST_LIB gtest_main gtest)

set(TEST_NAME test_cpp_diamm)

set(TEST_SOURCES
    test_DruckerPragerInducedAniso.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../DruckerPragerInducedAniso.cpp
)

add_executable(${TEST_NAME} ${TEST_SOURCES})

target_link_libraries(${TEST_NAME}
    ${VAANGO_LIBS}
    ${GTEST_LIB}
    Vaango::Options
    Vaango::Warnings
)

target_include_directories(${TEST_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/CCA/Components/MPM/ConstitutiveModel/J2I1PlasticModels
    ${PROJECT_SOURCE_DIR}/CCA/Components/MPM/ConstitutiveModel/fortran
)

# add_custom_command(
#   TARGET ${TEST_NAME}
#   POST_BUILD
#   COMMAND ${TEST_NAME}
# )
# gtest_discover_tests(${TEST_NAME})
