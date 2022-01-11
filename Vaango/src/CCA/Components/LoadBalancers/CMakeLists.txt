#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2022 Parresia Research Limited, New Zealand
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
# CMakeLists.txt for Vaango/src/CCA/Components/LoadBalancers

SET(Vaango_CCA_Components_LoadBalancers_SRCS
  LoadBalancerCommon.cc 
  LoadBalancerFactory.cc 
  RoundRobinLoadBalancer.cc 
  DynamicLoadBalancer.cc 
  SimpleLoadBalancer.cc 
  SingleProcessorLoadBalancer.cc 
  CostProfiler.cc 
  ProfileDriver.cc 
  CostModelForecaster.cc 
  ParticleLoadBalancer.cc
)

ADD_LIBRARY(Vaango_CCA_Components_LoadBalancers ${Vaango_CCA_Components_LoadBalancers_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_LoadBalancers
  Vaango_CCA_Components_ProblemSpecification
  Vaango_Core_Grid
  Vaango_Core_Util
  Vaango_Core_DataArchive
  Vaango_Core_ProblemSpec
  Vaango_CCA_Ports
  Vaango_Core_Parallel
  Vaango_Core_Exceptions
  Vaango_CCA_Components_ProblemSpecification
  Vaango_Core_Thread
  ${XML2_LIBRARY}
  ${MPI_LIBS}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_CCA_Components_LoadBalancers)
ENDIF(BUILD_SHARED_LIBS)
