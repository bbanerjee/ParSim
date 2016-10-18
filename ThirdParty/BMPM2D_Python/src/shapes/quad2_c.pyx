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
    inpf = np.array([h,patch.X0, patch.dX])
    idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2])
    labels = ['px','gx','cIdx','cW','cGrad','gDist']
    px,gx,cIdx,cW,cGrad, gDist = dw.getMult(labels,dwi)
    updateContribs( inpf, idxs, ng, th, px, gx, cIdx, cW, cGrad, gDist )


#===============================================================================        
def updateContribs( np.ndarray[FTYPE_t, ndim=2] inpf, 
                    np.ndarray[ITYPE_t, ndim=1] idxs,
                    ITYPE_t ng, FTYPE_t th, np.ndarray[FTYPE_t, ndim=2] px, 
                    np.ndarray[FTYPE_t, ndim=2] gx, 
                    np.ndarray[ITYPE_t, ndim=2] cIdx, 
                    np.ndarray[FTYPE_t, ndim=2] cW, 
                    np.ndarray[FTYPE_t, ndim=3] cGrad,
		    np.ndarray[FTYPE_t, ndim=1] gDist):
    # inpf - float input vector - h, dxdy, patch.X0, patch.dX
    cdef int nParts = px.shape[0]
    cdef int cc, idx 
    cdef int ii, jj, kk
    cdef double r, h, sgn, hm
    cdef double* hh = [inpf[0,0],inpf[0,1]]
    cdef double* x = [0.,0.]
    cdef double* pp = [0.,0.]
    cdef double* S = [0.,0.]
    cdef double* G = [0.,0.]
    cdef int* cix = [0,0]

    hm = min(hh[0],hh[1])    
    hm = sqrt(hh[0]*hh[0]+hh[1]*hh[1])
    
    for ii in range(nParts):
	    
	# Get Cell
        for kk in range(2):
            pp[kk] = px[ii,kk]
	    
            x[kk] = (pp[kk] - inpf[1,kk])/inpf[2,kk] + ng;
            cix[kk] = int(floor(x[kk]))
            if (x[kk] - 1.*cix[kk]) < 0.5:    cix[kk] -= 1
                
        cc = int(cix[1]*idxs[3] + cix[0]);

        for jj in range(9):
            idx = cc + idxs[jj];
            x[0] = pp[0]-gx[idx,0]
            x[1] = pp[1]-gx[idx,1]
            d = sqrt( x[0]*x[0] + x[1]*x[1] )
		
            for kk in range(2):
                r = fabs(x[kk])
                h = hh[kk]
                sgn = copysign(1.,x[kk])
		
                if ( r < 0.5*h ):
                    S[kk] = -r*r/(h*h) + 3./4
                    G[kk] = -2.*x[kk]/(h*h)
                elif ( r < 1.5*h ): 
                    S[kk] = r*r/(2.*h*h) - 3.*r/(2.*h) + 9./8
                    G[kk] = x[kk]/(h*h) - sgn*3./(2*h)        
                else: 
                    S[kk] = 0.
                    G[kk] = 0.	

            cIdx[ii,jj] = idx
            cW[ii,jj] = S[0]*S[1]
            cGrad[ii,jj,0] = S[1]*G[0]
            cGrad[ii,jj,1] = S[0]*G[1]
            gDist[idx] = max(0, max(gDist[idx], (1. - d/hm)))
	
    return 0