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
import operator
import warnings

def ellipseLvl( r, c, n=2, x=0 ):
    pwr = np.power
    tol = 1.e-14
    n = 1.*n
    def lvl(x):
        return 1. - pwr( sum(((x-c)/r)**n), 1/n )
    return lvl

def rectLvl( x0, x1, x=0 ):
    c = (x0+x1)/2.
    r = (x1-x0)/2.
    return ellipseLvl( r, c, 200. )


def fillLvl( lvl, patch ):
    nn = pCeil((patch.X1 - patch.X0) / (patch.dX/patch.ppe))
    ps = (patch.X1-patch.X0)/nn
    vol = patch.thick * ps[0] * ps[1]
    parts = []
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = patch.X0 + ps*ns
            if (patch.inPatch(pt) and lvl(pt) >= 0):
                parts.append(pt)

    parts = np.array(parts)
    return (parts, vol)    


def pCeil( x ):
    tol = 1.e-14
    return np.ceil(x-tol)