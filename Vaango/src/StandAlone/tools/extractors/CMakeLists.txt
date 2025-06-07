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

add_executable(selectpart selectpart.cc)
add_executable(partvarRange partvarRange.cc)
add_executable(timeextract timeextract.cc)
add_executable(lineextract lineextract.cc)
add_executable(partextract partextract.cc)
add_executable(faceextract faceextract.cc)
add_executable(extractF extractF.cc)
add_executable(extractV extractV.cc)
add_executable(extractPvec extractPvector.cc)
add_executable(extractPmat extractPmatrix.cc)
add_executable(extractPscalar extractPscalar.cc)
add_executable(extractS extractS.cc)
add_executable(extractPos extractPos.cc)
add_executable(extractPosVelMasVol extractPosVelMasVol.cc)

SET(EXTRACTOR_LIBS
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_DataArchive
  Vaango_Core_Util
  Vaango::Options
  Vaango::Warnings
)

target_link_libraries(selectpart ${EXTRACTOR_LIBS})
target_link_libraries(partvarRange ${EXTRACTOR_LIBS})
target_link_libraries(timeextract ${EXTRACTOR_LIBS})
target_link_libraries(lineextract ${EXTRACTOR_LIBS})
target_link_libraries(partextract ${EXTRACTOR_LIBS})
target_link_libraries(faceextract ${EXTRACTOR_LIBS})
target_link_libraries(extractF ${EXTRACTOR_LIBS})
target_link_libraries(extractV ${EXTRACTOR_LIBS})
target_link_libraries(extractPvec ${EXTRACTOR_LIBS})
target_link_libraries(extractPmat ${EXTRACTOR_LIBS})
target_link_libraries(extractPscalar ${EXTRACTOR_LIBS})
target_link_libraries(extractS ${EXTRACTOR_LIBS})
target_link_libraries(extractPos ${EXTRACTOR_LIBS})
target_link_libraries(extractPosVelMasVol ${EXTRACTOR_LIBS})
