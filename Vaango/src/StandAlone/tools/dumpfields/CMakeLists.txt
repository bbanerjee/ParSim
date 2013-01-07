IF(NOT WIN32)
  ADD_EXECUTABLE(dumpfields dumpfields.cc
    utils.cc
    Args.cc
    FieldSelection.cc
    FieldDiags.cc
    ScalarDiags.cc
    VectorDiags.cc
    TensorDiags.cc
    FieldDumper.cc
    TextDumper.cc
    EnsightDumper.cc
    InfoDumper.cc
    HistogramDumper.cc
  )

  TARGET_LINK_LIBRARIES(dumpfields
     Vaango_Core_DataArchive   
     Vaango_Core_OS
  )
ENDIF(NOT WIN32)
