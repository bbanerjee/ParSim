# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics/InternalForceComputer

SET(Vaango_CCA_Components_Peridynamics_SRCS
  ${Vaango_CCA_Components_Peridynamics_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/BondInternalForceComputer.cc  
  ${CMAKE_CURRENT_SOURCE_DIR}/ParticleInternalForceComputer.cc  
  PARENT_SCOPE
)

