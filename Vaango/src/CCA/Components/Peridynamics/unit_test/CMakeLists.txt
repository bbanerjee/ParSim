# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics/ParticleCreator

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(UT_SRCS utParticleTest.cc)
set(VAANGO_LIBS 
        Vaango_Core_Thread       
        Vaango_Core_DataArchive 
        Vaango_Core_Grid        
        Vaango_Core_Parallel    
        Vaango_Core_Labels      
        Vaango_Core_Util        
        Vaango_Core_Math        
        Vaango_Core_Disclosure  
        Vaango_Core_Exceptions  
        Vaango_Core_OS  
        Vaango_CCA_Ports        
        Vaango_CCA_Components_Parent 
        Vaango_CCA_Components_DataArchiver  
        Vaango_CCA_Components_LoadBalancers 
        Vaango_CCA_Components_Regridder     
        Vaango_Core_ProblemSpec             
        Vaango_CCA_Components_SimulationController 
        Vaango_CCA_Components_Schedulers           
        Vaango_CCA_Components_ProblemSpecification 
        Vaango_CCA_Components_Solvers
        Vaango_CCA_Components_MPM
        Vaango_CCA_Components_Peridynamics
        Peridynamics_DamageModels
        Peridynamics_MaterialModels
)

add_library(ut_lib ${UT_SRCS})
target_link_libraries(ut_lib ${VAANGO_LIBS} ${Boost_LIBRARIES}) 

add_executable(boostTest boostTest.cc)
target_link_libraries(boostTest ut_lib ${VAANGO_LIBS} ${Boost_LIBRARIES}) 

#add_executable(utTest1 utTest1.cc)
#target_link_libraries(utTest1 ut_lib ${VAANGO_LIBS} ${Boost_LIBRARIES})  

#add_executable(utTest0 utTest0.cc)
#target_link_libraries(utTest0 ${Boost_LIBRARIES} ${MPI_LIBRARIES})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_ut_lib)
ENDIF(BUILD_SHARED_LIBS)
