#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
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

cmake_minimum_required(VERSION 2.8.6)
find_package(Java)
include(UseJava)
project(UintahGui Java)
set(CMAKE_VERBOSE_MAKEFILE 1)
set(JAVA_INCLUDE_PATH "/usr/share/java")
set(J3DCore ${UintahGui_SOURCE_DIR}/ThirdParty/lib/ext/j3dcore.jar)
set(J3DUtils ${UintahGui_SOURCE_DIR}/ThirdParty/lib/ext/j3dutils.jar)
set(VecMath ${UintahGui_SOURCE_DIR}/ThirdParty/lib/ext/vecmath.jar)
set(CMAKE_JAVA_INCLUDE_PATH ${J3DCore})
set(CMAKE_JAVA_INCLUDE_PATH ${J3DUtils})
set(CMAKE_JAVA_INCLUDE_PATH ${VecMath})
set(JARS
    ${J3DCore} 
    ${J3DUtils} 
    ${VecMath} 
   )
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "${UintahGUI_SOURCE_DIR}/ThirdParty/lib/amd64")
#SET(CMAKE_JAVA_COMPILE_FLAGS "-g -Xlint:unchecked")
#include(${Vaango_SOURCE_DIR}/cmake/Macros/fortranDependency.cmake)
set(SRC_DIR vaango_ui)
set(BUILD_DIR "${UintahGui_SOURCE_DIR}/build")
message(STATUS ${UintahGui_SOURCE_DIR})
message(STATUS ${SRC_DIR})
set(SRC
  ${SRC_DIR}/UintahGui.java
  ${SRC_DIR}/BoxGeomPiece.java
  ${SRC_DIR}/BoxGeomPiecePanel.java
  ${SRC_DIR}/CompNeoHookMaterialPanel.java
  ${SRC_DIR}/ComputeParticleLocPanel.java
  ${SRC_DIR}/CreateGeomObjectPanel.java
  ${SRC_DIR}/CreateGeomPiecePanel.java
  ${SRC_DIR}/CylinderGeomPiece.java
  ${SRC_DIR}/DecimalField.java
  ${SRC_DIR}/DecimalVectorField.java
  ${SRC_DIR}/DifferenceGeomPiece.java
  ${SRC_DIR}/DisplayGeometryFrame.java
  ${SRC_DIR}/DisplayGeometryPanel.java
  ${SRC_DIR}/DisplayPartDistFrame.java
  ${SRC_DIR}/DisplayPartDistPanel.java
  ${SRC_DIR}/DisplayParticleLocFrame.java
  ${SRC_DIR}/DisplayParticleLocPanel.java
  ${SRC_DIR}/ElasticPlasticMaterialPanel.java
  ${SRC_DIR}/GeneralInputsPanel.java
  ${SRC_DIR}/GeometryPanel.java
  ${SRC_DIR}/GeomObject.java
  ${SRC_DIR}/GeomObjectPanel.java
  ${SRC_DIR}/GeomPiece.java
  ${SRC_DIR}/GeomPiecePanel.java
  ${SRC_DIR}/GridBCPanel.java
  ${SRC_DIR}/HelpAboutFrame.java
  ${SRC_DIR}/HypoElasticMaterialPanel.java
  ${SRC_DIR}/ICEInputsPanel.java
  ${SRC_DIR}/ICEMaterialInputPanel.java
  ${SRC_DIR}/ICEMaterialsPanel.java
  ${SRC_DIR}/InputGeometryPanel.java
  ${SRC_DIR}/InputPartDistPanel.java
  ${SRC_DIR}/IntegerField.java
  ${SRC_DIR}/IntegerVectorField.java
  ${SRC_DIR}/LightWeightCanvas.java
  ${SRC_DIR}/MPMContactInputPanel.java
  ${SRC_DIR}/MPMICEExchangePanel.java
  ${SRC_DIR}/MPMInputsPanel.java
  ${SRC_DIR}/MPMMaterialInputPanel.java
  ${SRC_DIR}/MPMMaterialsPanel.java
  ${SRC_DIR}/Packing.java
  ${SRC_DIR}/ParticleGeneratePanel.java
  ${SRC_DIR}/Particle.java
  ${SRC_DIR}/ParticleList.java
  ${SRC_DIR}/ParticleLocGeneratePanel.java
  ${SRC_DIR}/ParticleSizeDistInputPanel.java
  ${SRC_DIR}/ParticleSize.java
  ${SRC_DIR}/Point.java
  ${SRC_DIR}/PolygonDouble.java
  ${SRC_DIR}/RigidMaterialPanel.java
  ${SRC_DIR}/SmoothCylGeomPiece.java
  ${SRC_DIR}/TimeInputPanel.java
  ${SRC_DIR}/UintahInputPanel.java
  ${SRC_DIR}/UnionGeomPiece.java
  ${SRC_DIR}/VariableSaveInputPanel.java
  ${SRC_DIR}/Vector3D.java
  ${SRC_DIR}/ViscoSCRAMMaterialPanel.java
  ${SRC_DIR}/Voronoi.java
  ${SRC_DIR}/WholeNumberField.java
)
message(STATUS ${SRC})
add_jar(UintahGui 
        SOURCES ${SRC} 
        INCLUDE_JARS ${JARS} 
        ENTRY_POINT ${BUILD_DIR}/UintahGui.class
        OUTPUT_DIR ${BUILD_DIR})
