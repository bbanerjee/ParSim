#
# The MIT License
#
# Copyright (c) 2017- Parresia Research Limited, New Zealand
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

SET(IO_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/InputParameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OutputTecplot.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OutputVTK.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMBoundaryContactFileWriterCSV.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMBoundaryContactFileWriterXML.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMParticleContactFileWriterCSV.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMParticleContactFileWriterXML.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMParticleContactFileReaderXML.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMParticleFileReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DEMParticleFileWriter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PeriParticleFileReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/IOUtils.cpp
)

SET(ELLIP3D_SRCS
  ${ELLIP3D_SRCS}
  ${IO_SRCS}
  PARENT_SCOPE
)

add_subdirectory(UnitTests)
