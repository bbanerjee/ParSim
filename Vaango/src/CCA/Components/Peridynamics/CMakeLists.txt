# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics

SET(Vaango_CCA_Components_Peridynamics_SRCS
  Peridynamics.cc
  PeridynamicsFlags.cc
  PeridynamicsLabel.cc
  PeridynamicsMaterial.cc
  PeridynamicsDomainBoundCond.cc
)

add_subdirectory(ParticleCreator)
add_subdirectory(MaterialModels)
add_subdirectory(DamageModels)
add_subdirectory(FamilyComputer)
add_subdirectory(GradientComputer)
#add_subdirectory(BondBasedFailureModels)
if (BUILD_UNIT_TESTS)
  add_subdirectory(unit_test)
endif ()

SET(LIBS
  Vaango_CCA_Ports
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_Disclosure
  Vaango_Core_Labels
  Vaango_Core_ProblemSpec
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_Core_GeometryPiece
  Vaango_Core_Math
  Peridynamics_MaterialModels
  Peridynamics_DamageModels
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

ADD_LIBRARY(Vaango_CCA_Components_Peridynamics ${Vaango_CCA_Components_Peridynamics_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_Peridynamics ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_CCA_Components_Peridynamics)
ENDIF(BUILD_SHARED_LIBS)
