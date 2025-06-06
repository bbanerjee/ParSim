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

add_subdirectory(SimulationCommon)
add_subdirectory(Schedulers)
add_subdirectory(Regridder)
add_subdirectory(SimulationController)
add_subdirectory(SwitchingCriteria)
add_subdirectory(DataArchiver)
add_subdirectory(Examples)
add_subdirectory(LoadBalancers)
add_subdirectory(Parent)
add_subdirectory(ProblemSpecification)
add_subdirectory(Solvers)
add_subdirectory(PostProcessUda)

IF(MPM)
  add_subdirectory(MPM)
ENDIF(MPM)
IF(ICE)
  add_subdirectory(ICE)
  # only do Models with ICE
  add_subdirectory(Models)
ENDIF(ICE)
IF(MPM AND ICE)
  add_subdirectory(MPMICE)
ENDIF(MPM AND ICE)

add_subdirectory(Peridynamics)

add_subdirectory(OnTheFlyAnalysis)
