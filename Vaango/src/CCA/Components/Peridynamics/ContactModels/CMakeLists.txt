SET(Vaango_CCA_Components_Peridynamics_SRCS
  ${Vaango_CCA_Components_Peridynamics_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/SingleVelocityContact.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ContactModelFactory.cc   
  ${CMAKE_CURRENT_SOURCE_DIR}/ContactModelList.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/NullContact.cc      
  ${CMAKE_CURRENT_SOURCE_DIR}/ContactModelBase.cc
  PARENT_SCOPE
)

