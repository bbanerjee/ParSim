# CMakeLists.txt for Vaango/src/CCA/Components/Peridynamics/ParticleCreator

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
add_executable(utVaango utVaango.cc)
add_executable(utTest0 utTest0.cc)
add_executable(utTest1 utTest1.cc)

set(UT_INCS utBlank.h)
set(UT_SRCS utBlank.cc)
add_library(UT_LIBS ${UT_SRCS})

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
)

target_link_libraries(utVaango ${VAANGO_LIBS})
target_link_libraries(utTest1 ${VAANGO_LIBS} ${UT_LIBS}) 
target_link_libraries(utTest0 ${Boost_LIBRARIES} ${MPI_LIBRARIES})
#add_test(utTest0 ./utTest0)
