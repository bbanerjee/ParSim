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
#SET(CMAKE_JAVA_COMPILE_FLAGS "-g -Xlint:unchecked")
set(SRC
  UintahGui.java
  BoxGeomPiece.java
  BoxGeomPiecePanel.java
  CompNeoHookMaterialPanel.java
  ComputeParticleLocPanel.java
  CreateGeomObjectPanel.java
  CreateGeomPiecePanel.java
  CylinderGeomPiece.java
  DecimalField.java
  DecimalVectorField.java
  DifferenceGeomPiece.java
  DisplayGeometryFrame.java
  DisplayGeometryPanel.java
  DisplayPartDistFrame.java
  DisplayPartDistPanel.java
  DisplayParticleLocFrame.java
  DisplayParticleLocPanel.java
  ElasticPlasticMaterialPanel.java
  GeneralInputsPanel.java
  GeometryPanel.java
  GeomObject.java
  GeomObjectPanel.java
  GeomPiece.java
  GeomPiecePanel.java
  GridBCPanel.java
  HelpAboutFrame.java
  HypoElasticMaterialPanel.java
  ICEInputsPanel.java
  ICEMaterialInputPanel.java
  ICEMaterialsPanel.java
  InputGeometryPanel.java
  InputPartDistPanel.java
  IntegerField.java
  IntegerVectorField.java
  LightWeightCanvas.java
  MPMContactInputPanel.java
  MPMICEExchangePanel.java
  MPMInputsPanel.java
  MPMMaterialInputPanel.java
  MPMMaterialsPanel.java
  Packing.java
  ParticleGeneratePanel.java
  Particle.java
  ParticleList.java
  ParticleLocGeneratePanel.java
  ParticleSizeDistInputPanel.java
  ParticleSize.java
  Point.java
  PolygonDouble.java
  RigidMaterialPanel.java
  SmoothCylGeomPiece.java
  TimeInputPanel.java
  UintahInputPanel.java
  UnionGeomPiece.java
  VariableSaveInputPanel.java
  Vector3D.java
  ViscoSCRAMMaterialPanel.java
  Voronoi.java
  WholeNumberField.java
)
add_jar(UintahGui ${SRC} UintahGui.class)
