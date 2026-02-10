#------------------------
# Libraries needed
#------------------------
SET(VAANGO_LIBS 
        Vaango_Core_Parallel    
        MPM_ConstitutiveModel_FORTRAN
)
#------------------------
# Google tests
#------------------------
set(GTEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/submodules/googletest/googletest/include")
set(GTEST_LIB gtest_main gtest)
include_directories(${GTEST_INCLUDE_DIR})

set(FORTRAN_CM_CPP_SRC 
    ${CMAKE_CURRENT_SOURCE_DIR}/../MIGUtils.cc   
)

set(TEST_NAME test_fortran_diamm)

set(TEST_SOURCES
    test_diamm.c
    ${FORTRAN_CM_CPP_SRC}
)

add_executable(${TEST_NAME} ${TEST_SOURCES})

target_link_libraries(${TEST_NAME}
    ${VAANGO_LIBS}
    ${GTEST_LIB}
)

target_include_directories(${TEST_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/CCA/Components/MPM/ConstitutiveModel/fortran
)

# add_custom_command(
#   TARGET ${TEST_NAME}
#   POST_BUILD
#   COMMAND ${TEST_NAME}
# )
# gtest_discover_tests(${TEST_NAME})