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
SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/FlowStressModel.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/FlowStressModelFactory.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/LinearHardeningFlow.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/JohnsonCookFlow.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/ZAFlow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ZAPolymerFlow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MTSFlow.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/SCGFlow.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/PTWFlow.cc 
  PARENT_SCOPE
)
