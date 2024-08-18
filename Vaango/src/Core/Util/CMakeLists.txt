#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2024 Biswajit Banerjee
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#
# CMakeLists.txt for Vaango/src/Core/Util

set(Vaango_Core_Util_SRCS
  DebugStream.cc        
  DOUT.cpp
  Endian.cc             
  Environment.cc        
  FileUtils.cc          
  ProgressiveWarning.cc 
  RWS.cc                
  SizeTypeConvert.cc    
  soloader.cc           
  StringUtil.cc
  SysUtils.cc
  Util.cc
  XMLUtils.cc
  ${REFCOUNT_IMPL}
)

if (HAVE_CUDA)
  set(Vaango_Core_Util_SRCS
      ${Vaango_Core_Util_SRCS}
      GPU.cu
  )
endif (HAVE_CUDA)

add_library(Vaango_Core_Util ${Vaango_Core_Util_SRCS})

set(LIBS 
  Vaango_Core_Containers
  Vaango_Core_Exceptions
  Vaango_Core_Malloc
)
if (WIN32)
  set(LIBS
    ${LIBS}
    Vaango_Core_OS
  )
endif()

if (HAVE_CUDA)
  set(LIBS
    ${LIBS}
    GPU.o
  )
endif (HAVE_CUDA)

# message(STATUS "Core/Util Linking to libraries ${LIBS} ${DL_LIBRARY} ${THREAD_LIBRARY} ${SOCKET_LIBRARY}")

target_link_libraries(Vaango_Core_Util ${LIBS}
   ${DL_LIBRARY}
   ${Z_LIBRARY} 
   ${CUDA_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_Core_Util)
ENDIF(BUILD_SHARED_LIBS)
