#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
#  Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# 
#  
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and_or sell copies of the Software, and to permit persons to whom the
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

# CMakeLists.txt for Vaango_StandAlone_tools_uda2vtk

set(UDA2VTK_SRC
  uda2vtk.cc
)

set(UDA2VTK_LIBS
  Vaango_Core_Exceptions    
  Vaango_Core_Grid          
  Vaango_Core_Util          
  Vaango_Core_Math          
  Vaango_Core_ProblemSpec   
  Vaango_Core_DataArchive   
  Vaango_Core_Parallel      
  Vaango_CCA_Ports          
  Vaango_CCA_Components_ProblemSpecification 
  Vaango_Core_Basis        
  Vaango_Core_Containers   
  Vaango_Core_Datatypes    
  ${XML2_LIBRARY}
  ${Z_LIBRARY}
  ${M_LIBRARY}
  ${MPI_LIBRARY}
  ${F_LIBRARY}
)

add_executable(uda2vtk ${UDA2VTK_SRC})
target_link_libraries(uda2vtk ${UDA2VTK_LIBS})
	
