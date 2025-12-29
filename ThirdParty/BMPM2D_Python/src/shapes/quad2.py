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

@njit(cache=True)
def updateContribs(inpf, idxs, ng, th, px, gx, cIdx, cW, cGrad, gDist):
    nParts = px.shape[0]
    hh = np.array([inpf[0,0], inpf[0,1]])
    hm = np.sqrt(hh[0]*hh[0]+hh[1]*hh[1])
    
    for ii in range(nParts):
        pp = np.zeros(2)
        cix = np.zeros(2, dtype=np.int64)
        
        # Get Cell
        for kk in range(2):
            pp[kk] = px[ii,kk]
            val = (pp[kk] - inpf[1,kk])/inpf[2,kk] + ng
            cix[kk] = int(np.floor(val))
            if (val - 1.0*cix[kk]) < 0.5:
                cix[kk] -= 1
        
        cc = cix[1]*idxs[3] + cix[0]

        for jj in range(9):
            idx = cc + idxs[jj]
            dx = pp[0]-gx[idx,0]
            dy = pp[1]-gx[idx,1]
            d = np.sqrt( dx*dx + dy*dy )
            
            S = np.zeros(2)
            G = np.zeros(2)
        
            for kk in range(2):
                x = pp[kk] - gx[idx,kk]
                r = np.abs(x)
                h = hh[kk]
                if x >= 0: sgn = 1.0
                else: sgn = -1.0
        
                if ( r < 0.5*h ):
                    S[kk] = -r*r/(h*h) + 0.75
                    G[kk] = -2.0*x/(h*h)
                elif ( r < 1.5*h ): 
                    S[kk] = r*r/(2.0*h*h) - 1.5*r/h + 1.125
                    G[kk] = x/(h*h) - sgn*1.5/h        
                else: 
                    S[kk] = 0.0
                    G[kk] = 0.0

            cIdx[ii,jj] = idx
            cW[ii,jj] = S[0]*S[1]
            cGrad[ii,jj,0] = S[1]*G[0]
            cGrad[ii,jj,1] = S[0]*G[1]
            
            val = 1.0 - d/hm
            old = gDist[idx]
            if val > old: old = val
            if 0.0 > old: old = 0.0
            gDist[idx] = old

def updateContribList( dw, patch, dwi ):
    # Update node contribution list
    nx = patch.Nc[0]
    th = patch.thick
    h = patch.dX
    ng = patch.nGhost
    inpf = np.zeros((3,2))
    inpf[0,:] = h
    inpf[1,:] = patch.X0
    inpf[2,:] = patch.dX
    
    idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2], dtype=np.int64)
    labels = ['px','gx','cIdx','cW','cGrad','gDist']
    px,gx,cIdx,cW,cGrad, gDist = dw.getMult(labels,dwi)
    updateContribs( inpf, idxs, int(ng), float(th), px, gx, cIdx, cW, cGrad, gDist )