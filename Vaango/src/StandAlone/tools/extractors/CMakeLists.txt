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

ADD_EXECUTABLE(selectpart selectpart.cc)
ADD_EXECUTABLE(partvarRange partvarRange.cc)
ADD_EXECUTABLE(timeextract timeextract.cc)
ADD_EXECUTABLE(lineextract lineextract.cc)
ADD_EXECUTABLE(partextract partextract.cc)
ADD_EXECUTABLE(faceextract faceextract.cc)
ADD_EXECUTABLE(extractF extractF.cc)
ADD_EXECUTABLE(extractV extractV.cc)
ADD_EXECUTABLE(extractPvec extractPvector.cc)
ADD_EXECUTABLE(extractPmat extractPmatrix.cc)
ADD_EXECUTABLE(extractPscalar extractPscalar.cc)
ADD_EXECUTABLE(extractS extractS.cc)
ADD_EXECUTABLE(extractPos extractPos.cc)
ADD_EXECUTABLE(extractPosVelMasVol extractPosVelMasVol.cc)

SET(EXTRACTOR_LIBS
  Vaango_Core_Exceptions
  Vaango_Core_Grid
  Vaango_Core_DataArchive
  Vaango_Core_Util
)

TARGET_LINK_LIBRARIES(selectpart ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(partvarRange ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(timeextract ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(lineextract ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(partextract ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(faceextract ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractF ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractV ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractPvec ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractPmat ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractPscalar ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractS ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractPos ${EXTRACTOR_LIBS})
TARGET_LINK_LIBRARIES(extractPosVelMasVol ${EXTRACTOR_LIBS})
