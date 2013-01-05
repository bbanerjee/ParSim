# CMakeLists.txt for Packages/Uintah/CCA/Components/MPM/GradientComputer

SET(Packages_Uintah_CCA_Components_MPM_SRCS
  ${Packages_Uintah_CCA_Components_MPM_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/DeformationGradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/DisplacementGradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/GradientComputer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/VelocityGradientComputer.cc
  PARENT_SCOPE
)

