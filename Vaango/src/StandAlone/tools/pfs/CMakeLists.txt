ADD_EXECUTABLE(pfs pfs.cc)
ADD_EXECUTABLE(pfs2 pfs2.cc)

SET(PFS_LIBS
  Vaango_Core_Exceptions
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_Parallel
  Vaango_Core_ProblemSpec
  Vaango_Core_Util
  Vaango_Core_GeometryPiece
  Vaango_CCA_Components_ProblemSpecification
)

TARGET_LINK_LIBRARIES(pfs ${PFS_LIBS})
TARGET_LINK_LIBRARIES(pfs2 ${PFS_LIBS})
