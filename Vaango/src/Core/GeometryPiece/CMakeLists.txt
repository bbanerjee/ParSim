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
# CMakeLists.txt for Vaango/src/Core/GeometryPiece

SET(Vaango_Core_GeometryPiece_SRCS
  GeometryObject.cc 
  GeometryPiece.cc 
  SphereGeometryPiece.cc 
  SphereMembraneGeometryPiece.cc 
  NaaBoxGeometryPiece.cc 
  BoxGeometryPiece.cc 
  CylinderGeometryPiece.cc 
  TorusGeometryPiece.cc 
  ConeGeometryPiece.cc 
  TriGeometryPiece.cc 
  UniformGrid.cc 
  UnionGeometryPiece.cc 
  DifferenceGeometryPiece.cc 
  IntersectionGeometryPiece.cc 
  FileGeometryPiece.cc 
  NullGeometryPiece.cc 
  GeometryPieceFactory.cc 
  ShellGeometryPiece.cc 
  ShellGeometryFactory.cc 
  CylinderShellPiece.cc 
  PlaneShellPiece.cc 
  SphereShellPiece.cc 
  SpecialGeomPiece.cc 
  SmoothCylGeomPiece.cc 
  SmoothSphereGeomPiece.cc 
  CorrugEdgeGeomPiece.cc 
  EllipsoidGeometryPiece.cc 
  AbaqusMeshGeometryPiece.cc
)

add_library(Vaango_Core_GeometryPiece ${Vaango_Core_GeometryPiece_SRCS})

target_link_libraries(Vaango_Core_GeometryPiece
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_Math
  Vaango_Core_Parallel
  Vaango_Core_ProblemSpec
  Vaango_Core_Util
  Vaango::Options
  Vaango::Warnings
)

IF(BUILD_SHARED_LIBS)
  add_definitions(-DBUILD_Vaango_Core_GeometryPiece)
ENDIF(BUILD_SHARED_LIBS)

install(TARGETS
  Vaango_Core_GeometryPiece
)