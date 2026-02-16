#------------------------
# Vaango libs
#------------------------
SET(VAANGO_LIBS 
        Vaango_Core_Grid 
        Vaango_Core_Parallel    
        Vaango_Core_Util        
        Vaango_Core_Math        
        Vaango_Core_Exceptions  
        Vaango_Core_OS  
        Vaango_Core_GeometryPiece
        Vaango_Core_ProblemSpec
)

#------------------------
# Google tests
#------------------------

set(GTEST_LIB gtest_main gtest)


add_executable(testLocalSDF testLocalSDF.cc)
target_link_libraries(testLocalSDF PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)

add_executable(testDynamicSDFGeometry testDynamicSDFGeometry.cc)
target_link_libraries(testDynamicSDFGeometry PRIVATE
  ${GTEST_LIB}
  ${VAANGO_LIBS}
  LibXml2::LibXml2
  Vaango::Options
  Vaango::Warnings
)
