#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
# Copyright (c) 2015-2023 Biswajit Banerjee
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

SET(Vaango_CCA_Components_SwitchingCriteria_SRCS
  SwitchingCriteriaFactory.cc
  DDT1Criterion.cc
  None.cc
  TimestepNumber.cc
  SimpleBurn.cc
  SteadyBurn.cc
  SteadyState.cc
)

ADD_LIBRARY(Vaango_CCA_Components_SwitchingCriteria ${Vaango_CCA_Components_SwitchingCriteria_SRCS})

TARGET_LINK_LIBRARIES(Vaango_CCA_Components_SwitchingCriteria
  Vaango_CCA_Ports
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_Parallel
  
  Vaango_Core_Parallel
  Vaango_Core_Util
  Vaango_Core_ProblemSpec
  ${XML2_LIBRARY}
)

IF(BUILD_SHARED_LIBS)
  ADD_DEFINITIONS(-DBUILD_Vaango_CCA_Components_SwitchingCriteria)
ENDIF(BUILD_SHARED_LIBS)
