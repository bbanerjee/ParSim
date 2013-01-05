#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#

# CMakeLists.txt for Packages/Uintah/CCA/Components/MPM/ConstitutiveModel

SET(MPM_ConstitutiveModel_SRCS
  RigidMaterial.cc              
  CompMooneyRivlin.cc           
  ConstitutiveModelFactory.cc   
  ConstitutiveModel.cc          
  ImplicitCM.cc                 
  MPMMaterial.cc                
  CNH_MMS.cc                    
  TransIsoHyper.cc              
  TransIsoHyperImplicit.cc      
  ViscoTransIsoHyper.cc         
  ViscoTransIsoHyperImplicit.cc 
  ViscoScram.cc                 
  ViscoSCRAMHotSpot.cc          
  HypoElastic.cc                
  HypoElasticImplicit.cc        
  ViscoScramImplicit.cc         
  MWViscoElastic.cc             
  IdealGasMP.cc                 
  Membrane.cc                   
  ShellMaterial.cc              
  ElasticPlasticHP.cc           
  Water.cc                      
  ViscoPlastic.cc               
  MurnaghanMPM.cc               
  ProgramBurn.cc                
  JWLppMPM.cc                   
  UCNH.cc                       
  P_Alpha.cc                    
  SoilFoam.cc	             
  NonLocalDruckerPrager.cc      
  Arenisca.cc
)

if (FORTRAN)
  SET(MPM_ConstitutiveModel_SRCS
    ${MPM_ConstitutiveModel_SRCS}
    Diamm.cc                      
    HypoElasticFortran.cc         
    Kayenta.cc                    
  )
  add_subdirectory(fortran)
endif(FORTRAN)

add_subdirectory(PlasticityModels)
add_subdirectory(Biswajit)

SET(LIBS
  Packages_Uintah_CCA_Ports
  Packages_Uintah_Core_Grid
  Packages_Uintah_Core_Util
  Packages_Uintah_Core_Disclosure
  Packages_Uintah_Core_Labels
  Packages_Uintah_Core_ProblemSpec
  Packages_Uintah_Core_Parallel
  Packages_Uintah_Core_Exceptions
  Packages_Uintah_Core_GeometryPiece
  Packages_Uintah_Core_Math
  MPM_ConstitutiveModel_PlasticityModels
  MPM_ConstitutiveModel_Biswajit
  ${XML2_LIBRARY}
  ${VT_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
)

ADD_LIBRARY(MPM_ConstitutiveModel ${MPM_ConstitutiveModel_SRCS})

TARGET_LINK_LIBRARIES(MPM_ConstitutiveModel ${LIBS})

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_MPM_ConstitutiveModel)
ENDIF(BUILD_SHARED_LIBS)




