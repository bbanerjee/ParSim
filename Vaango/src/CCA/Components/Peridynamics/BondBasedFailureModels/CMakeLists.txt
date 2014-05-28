# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics/FailureModels

SET(Peridynamics_FailureModel_SRCS
  PeridynamicsFailureModel.cc
  PeridynamicsFailureModelFactory.cc
  BondDamageModel.cc
)

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
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

ADD_LIBRARY(Peridynamics_FailureModels ${Peridynamics_FailureModel_SRCS})

TARGET_LINK_LIBRARIES(Peridynamics_FailureModels ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Peridynamics_FailureModels)
ENDIF(BUILD_SHARED_LIBS)




