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
from numba import njit
from itertools import count
from dateutil.relativedelta import relativedelta as reldelta

# Utils for moving data between particles and grid

@njit(cache=True)
def integrate( cIdx, cW, pp, gg ):
    # Integrate particle values to grid (p->g)
    nParts = pp.shape[0]
    nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj]
            w = cW[ii,jj]
            gg[ixc,0] += pp[ii,0] * w
            gg[ixc,1] += pp[ii,1] * w
    return gg

@njit(cache=True)
def interpolate( cIdx, cW, pp, gg ):
    # Interpolate grid values to particles pp
    nParts = pp.shape[0]
    nContrib = cIdx.shape[1]
    for ii in range(nParts):
        pp[ii,0] = 0.0
        pp[ii,1] = 0.0
        for jj in range(nContrib):
            ixc = cIdx[ii,jj]
            w = cW[ii,jj]
            pp[ii,0] += gg[ixc,0] * w
            pp[ii,1] += gg[ixc,1] * w
    return pp

@njit(cache=True)
def gradient( cIdx, cGrad, pp, gg ):
    # Interpolate a gradient
    nParts = pp.shape[0]
    nContrib = cIdx.shape[1]
    for ii in range(nParts):
        pp[ii,0,0] = 0.0
        pp[ii,0,1] = 0.0
        pp[ii,1,0] = 0.0
        pp[ii,1,1] = 0.0
        for jj in range(nContrib):
            ixc = cIdx[ii,jj]
            cg0 = cGrad[ii,jj,0]
            cg1 = cGrad[ii,jj,1]
            pp[ii,0,0] += gg[ixc,0] * cg0
            pp[ii,0,1] += gg[ixc,0] * cg1
            pp[ii,1,0] += gg[ixc,1] * cg0
            pp[ii,1,1] += gg[ixc,1] * cg1
    return pp

@njit(cache=True)
def divergence( cIdx, cGrad, pp, gg ):
    # Send divergence of particle field to the grid
    nParts = pp.shape[0]
    nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj]
            cg0 = cGrad[ii,jj,0]
            cg1 = cGrad[ii,jj,1]
            gg[ixc,0] -= pp[ii,0,0] * cg0 + pp[ii,1,0] * cg1
            gg[ixc,1] -= pp[ii,1,0] * cg0 + pp[ii,1,1] * cg1
    return gg   

@njit(cache=True)
def gradscalar( cIdx, cGrad, pp, gg ):
    # Send gradient of particle field to the grid
    nParts = pp.shape[0]
    nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj]
            gg[ixc,0] += pp[ii,0] * cGrad[ii,jj,0]
            gg[ixc,1] += pp[ii,0] * cGrad[ii,jj,1]
    return gg   

@njit(cache=True)
def dotAdd( pp, qq ):
    # return pp += qq dot pp
    nParts = pp.shape[0]
    dot = np.zeros((2,2))
    for ii in range(nParts):
        for jj in range(2):
            for kk in range(2):
                dot[jj,kk] = qq[ii,jj,0]*pp[ii,0,kk]+qq[ii,jj,1]*pp[ii,1,kk]
        for jj in range(2):
            for kk in range(2):
                pp[ii,jj,kk] += dot[jj,kk]

def readableTime( tt ):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']    
    h_read = lambda delta: ['%d %s' % (getattr(delta, attr), 
                   attr if getattr(delta, attr) > 1 else attr[:-1]) 
    for attr in attrs if getattr(delta, attr)]
    tm = h_read(reldelta(seconds=tt))
    
    st = ''
    for t in tm:
        st += t + ', '
    return st[:-2]