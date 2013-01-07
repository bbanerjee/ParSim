ADD_EXECUTABLE(SimpleMath SimpleMath.cc)

TARGET_LINK_LIBRARIES(SimpleMath
        Vaango_Core_Exceptions    
        Vaango_Core_Grid          
        Vaango_Core_Util          
        Vaango_Core_Math          
        Vaango_Core_ProblemSpec   
        Vaango_Core_Parallel      
        Vaango_Core_Disclosure    
        Vaango_Core_Persistent  
        Vaango_Core_Thread      
        Vaango_Core_Containers
)
