IF(NOT WIN32)
  # doesn't build on windows....
  ADD_EXECUTABLE(graphview graphview.cc DaVinci.cc GV_TaskGraph.cc)

  TARGET_LINK_LIBRARIES(graphview
    Vaango_Core_DataArchive
    Vaango_Core_ProblemSpec
    Vaango_Core_Util
    Vaango_CCA_Components_ProblemSpecification
    Vaango_Core_Exceptions)

ENDIF(NOT WIN32)
