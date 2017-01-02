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

import numpy as np

def getMomentum( data, dwi ):
    t = data.keys()
    t.sort()
    dt = t[1]-t[0]
    pw = []
    pf = []
    for t0 in t:
        pw0 = data[t0].dw['pw',dwi]
        pm0 = data[t0].dw['pm',dwi]
        pa0 = data[t0].dw['pvI',dwi]
        pf0 = pa0*pm0*dt
        mpw0 = np.sqrt((pw0*pw0).sum(axis=1)[:,np.newaxis])
        mpf0 = np.sqrt((pf0*pf0).sum(axis=1)[:,np.newaxis])
        pw.append(sum(pw0))
        pf.append(sum(mpf0))
        
    t = np.array(t)
    pw = np.array(pw)
    pf = np.array(pf)
    return t,pw,pf