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
#-----------------------------------------------------------------------------
# The models
#-----------------------------------------------------------------------------
#message(STATUS "Compiling FORTRAN code? ${FORTRAN}")
if (FORTRAN)
  SET(MPM_Fortran_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/Diamm.cc                      
    ${CMAKE_CURRENT_SOURCE_DIR}/Kayenta.cc                    
  )
endif(FORTRAN)

SET(MPM_ConstitutiveModel_SRCS
  ${MPM_ConstitutiveModel_SRCS}
  ${CMAKE_CURRENT_SOURCE_DIR}/Arena.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ArenaMixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/Arenisca.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/Arenisca3.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/Arenisca4.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CamClay.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulomb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/NonLocalDruckerPrager.cc      
  ${CMAKE_CURRENT_SOURCE_DIR}/SoilModelBrannon.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombState.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombTypes.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombBase.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombSheng.cc 
  ${CMAKE_CURRENT_SOURCE_DIR}/MohrCoulombClassic.cc 
  ${MPM_Fortran_SRCS}
  PARENT_SCOPE
)
