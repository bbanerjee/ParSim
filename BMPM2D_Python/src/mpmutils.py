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
from itertools import izip, count
from dateutil.relativedelta import relativedelta as reldelta
# Utils for moving data between particles and grid

def integrate( cIdx, cW, pp, gg ):
    # Integrate particle values to grid (p->g)
    for (ppi,idx,w) in izip(pp,cIdx,cW):
        gg[idx] += ppi * w[:,np.newaxis]
    return gg

def interpolate( cIdx, cW, pp, gg ):
    # Interpolate grid values to particles pp
    for (ii,idx,w) in izip(count(),cIdx,cW):
        pp[ii] = np.sum( gg[idx] * w[:,np.newaxis], 0 )
    return pp

def gradient( cIdx, cGrad, pp, gg ):
    # Interpolate a gradient
    for (ii,idxi,gradi) in izip(count(),cIdx,cGrad):
        pp[ii] = [0]
        for (idx,grad) in izip( idxi, gradi ):
            gR = np.reshape( gg[idx], (2,1) )
            cg = np.reshape( grad, (1,2) )
            pp[ii] += np.dot( gR, cg )
    return pp

def divergence( cIdx, cGrad, pp, gg ):
    # Send divergence of particle field to the grid
    for (ppi,idxi,gradi) in izip(pp,cIdx,cGrad):
        for idx,grad in izip( idxi, gradi ):
            cg = np.reshape( grad, (2,1) )            
            gg[idx] -= np.reshape( np.dot( ppi, cg ), 2 )
    return gg   

def gradscalar( cIdx, cGrad, pp, gg ):
    # Send gradient of particle field to the grid
    for (ppi,idxi,gradi) in izip(pp,cIdx,cW):
        for idx,grad in izip( idxi, gradi ):
            gg[idx] += np.reshape( ppi * grad, 2 )
    return gg   


def dotAdd( pp, qq ):
    # return pp += qq dot pp
    for (ii,ppi,qqi) in izip(count(),pp,qq):
        pp[ii] += np.dot( qqi, ppi )


def readableTime( tt ):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']    
    h_read = lambda delta: ['%d %s' % (getattr(delta, attr), 
                   getattr(delta, attr) > 1 and attr or attr[:-1]) 
    for attr in attrs if getattr(delta, attr)]
    tm = h_read(reldelta(seconds=tt))
    
    st = ''
    for t in tm:
        st += t + ', '
    return st[:-2]