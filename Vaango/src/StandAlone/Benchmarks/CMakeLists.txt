ADD_EXECUTABLE(SimpleMath SimpleMath.cc)

TARGET_LINK_LIBRARIES(SimpleMath
        Packages_Uintah_Core_Exceptions    
        Packages_Uintah_Core_Grid          
        Packages_Uintah_Core_Util          
        Packages_Uintah_Core_Math          
        Packages_Uintah_Core_ProblemSpec   
        Packages_Uintah_Core_Parallel      
        Packages_Uintah_Core_Disclosure    
        Packages_Uintah_Core_Persistent  
        Packages_Uintah_Core_Thread      
        Packages_Uintah_Core_Containers
)
