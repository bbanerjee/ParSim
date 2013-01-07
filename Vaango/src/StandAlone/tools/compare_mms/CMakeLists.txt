ADD_EXECUTABLE(compare_mms compare_mms.cc ExpMMS.cc LinearMMS.cc SineMMS.cc)
TARGET_LINK_LIBRARIES(compare_mms 
        Vaango_Core_DataArchive 
)
