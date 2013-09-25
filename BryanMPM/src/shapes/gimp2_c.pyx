# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
cimport cython
import numpy as np
cimport numpy as np
ctypedef np.int_t ITYPE_t
ctypedef np.double_t FTYPE_t
cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double x)
    double floor(double x)
    double copysign(double x, double y)
    
ITYPE = np.int
FTYPE = np.float


#===============================================================================    
def updateContribList( dw, patch, dwi ):
    # Update node contribution list
    nx = patch.Nc[0]
    th = patch.thick
    h = patch.dX
    dxdy = h[::-1]/h
    ng = patch.nGhost
    inpf = np.array([h,dxdy,patch.X0, patch.dX])
    idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2])  
    labels = ['px','pVol','pF','gx','cIdx','cW','cGrad','gDist']
    px,pVol,pF,gx,cIdx,cW,cGrad,gDist = dw.getMult(labels,dwi)
    updateContribs( inpf, idxs, ng, th, px, pVol, pF, 
		     gx, cIdx, cW, cGrad, gDist )


#===============================================================================        
def updateContribs( np.ndarray[FTYPE_t, ndim=2] inpf, 
                    np.ndarray[ITYPE_t, ndim=1] idxs,
                    ITYPE_t ng, FTYPE_t th, np.ndarray[FTYPE_t, ndim=2] px, 
                    np.ndarray[FTYPE_t, ndim=1] pVol, 
                    np.ndarray[FTYPE_t, ndim=3] pF, 
                    np.ndarray[FTYPE_t, ndim=2] gx, 
                    np.ndarray[ITYPE_t, ndim=2] cIdx, 
                    np.ndarray[FTYPE_t, ndim=2] cW, 
                    np.ndarray[FTYPE_t, ndim=3] cGrad,
                    np.ndarray[FTYPE_t, ndim=1] gDist ):
    # inpf - float input vector - h, dxdy, patch.X0, patch.dX
    cdef int nParts = px.shape[0]
    cdef int cc, idx 
    cdef int ii, jj, kk
    cdef double x, r, h, l, sgn, hm, xx, xy
    cdef double* hh = [inpf[0,0],inpf[0,1]]
    cdef double* pp = [0.,0.]
    cdef double* ll = [0.,0.]
    cdef double* S = [0.,0.]
    cdef double* G = [0.,0.]
    cdef int* cix = [0,0]
    
    hm = min(hh[0],hh[1])    
    hm = sqrt(hh[0]*hh[0]+hh[1]*hh[1])
        
    for ii in range(nParts):
	    
	# Get Cell
        for kk in range(2):
            pp[kk] = px[ii,kk]
            ll[kk] = sqrt(pVol[ii] / (4.*th*inpf[1,kk])) * pF[ii,kk,kk]
	    
            x = (pp[kk] - inpf[2,kk])/inpf[3,kk] + ng;
            cix[kk] = int(floor(x))
            if (x - 1.*cix[kk]) < 0.5:    cix[kk] -= 1
                
        cc = int(cix[1]*idxs[3] + cix[0]);

        for jj in range(9):
            idx = cc + idxs[jj];
            xx = pp[0]-gx[idx,0]
            xy = pp[1]-gx[idx,1]
            d = sqrt( xx*xx + xy*xy )
		
            for kk in range(2):
                x = pp[kk] - gx[idx,kk]
                r = fabs(x)
                l = ll[kk]
                h = hh[kk]
                sgn = copysign(1.,x)		    
                if (r < l):
                    S[kk] = 1. - (r*r+l*l)/(2.*h*l)
                    G[kk] = -x/(h*l)
                elif (r < h-l):
                    S[kk] = 1. - r/h
                    G[kk] = -sgn/h
                elif (r < h+l):
                    S[kk] = (h+l-r)*(h+l-r) / (4.*h*l)
                    G[kk] = (h+l-r) / (-2.*sgn*h*l)
                else:
                    S[kk] = 0.
                    G[kk] = 0.
		    		
            cIdx[ii,jj] = idx
            cW[ii,jj] = S[0]*S[1]
            cGrad[ii,jj,0] = S[1]*G[0]
            cGrad[ii,jj,1] = S[0]*G[1]
            gDist[idx] = max(0, max(gDist[idx], (1. - d/hm)))
	
    return 0