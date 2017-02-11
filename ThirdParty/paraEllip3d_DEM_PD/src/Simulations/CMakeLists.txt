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

SET(COMMAND_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/CavityExpansionCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CavityExpansionResumeCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CommandHandler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CoupledFluidFlowCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DepositIntoContainerCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DepositIntoContainerResumeCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/IsotropicLoadingCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OedometerLoadingCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PeridynamicsPullOutCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PeridynamicsRigidInclusionCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PlaneStrainLoadingCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ProceedFromPresetCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TriaxialLoadingCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TrimParticlesCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TrueTriaxialLoadingCommand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/TuneMassPercentageCommand.cpp
)

SET(ELLIP3D_SRCS
  ${ELLIP3D_SRCS}
  ${COMMAND_SRCS}
  PARENT_SCOPE
)
