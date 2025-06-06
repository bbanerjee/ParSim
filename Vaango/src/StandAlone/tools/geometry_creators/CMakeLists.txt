#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2025 Biswajit Banerjee, Parresia Research Ltd., NZ
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

add_executable(createFileGeomPieceFromUda createFileGeomPieceFromUda.cc)
add_executable(gambitFileReader gambitFileReader.cc)
add_executable(imageFromGeom ImageFromGeom.cc)
add_executable(objFileReader OBJFileReader.cc)
add_executable(particleFileSplitter particleFileSplitter.cc)
add_executable(particleImageFileSplitter particleImageFileSplitter.cc)
add_executable(plyFileReader PLYFileReader.cc)
add_executable(rawToUniqueGrains rawToUniqueGrains.cc)

SET(GEOM_CREATE_LIBS
  Vaango_Core_Exceptions
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_Parallel
  Vaango_Core_ProblemSpec
  Vaango_Core_Util
  Vaango_Core_GeometryPiece
  Vaango_CCA_Components_ProblemSpecification
  Vaango::Options
  Vaango::Warnings
)

target_link_libraries(createFileGeomPieceFromUda 
                      ${VAANGO_LIBS})
target_link_libraries(gambitFileReader ${GEOM_CREATE_LIBS})
target_link_libraries(imageFromGeom ${GEOM_CREATE_LIBS})
target_link_libraries(objFileReader ${GEOM_CREATE_LIBS})
target_link_libraries(particleFileSplitter ${GEOM_CREATE_LIBS})
target_link_libraries(particleImageFileSplitter ${GEOM_CREATE_LIBS})
target_link_libraries(plyFileReader ${GEOM_CREATE_LIBS})
target_link_libraries(rawToUniqueGrains ${GEOM_CREATE_LIBS})
