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
    
ITYPE = np.int
FTYPE = np.float

#===============================================================================
class Shape:
    #  Shape functions - compute nodal contributions to particle values
    def __init__(self):
        self.dim = 2;
        self.S = np.zeros([self.dim,1])    # Value of Shape function
        self.G = np.zeros([self.dim,1])    # Value of Shape function derivative	


#===============================================================================
class GIMP(Shape):
    def __init__(self):
        self.nSupport = 9
        self.nGhost = 2
        Shape.__init__(self)
	    

    def updateContribList( self, dw, patch, mIdx ):
        # Update node contribution list
        nx = patch.Nc[0]
        th = patch.thick
        h = patch.dX
        dxdy = h[::-1]/h
        ng = patch.nGhost
        inpf = np.array([h,dxdy,patch.X0, patch.dX])
        idxs = np.array([0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2])
        GIMP.updateContribs( inpf, idxs, ng, th, dw.px, dw.pVol, dw.pF, 
			     dw.gx, dw.cIdx, dw.cW, dw.cGrad, mIdx )

        
    @staticmethod	
    def updateContribs( np.ndarray[FTYPE_t, ndim=2] inpf, 
                        np.ndarray[ITYPE_t, ndim=1] idxs,
                        ITYPE_t ng, FTYPE_t th, np.ndarray[FTYPE_t, ndim=2] px, 
                        np.ndarray[FTYPE_t, ndim=1] pVol, 
                        np.ndarray[FTYPE_t, ndim=3] pF, 
                        np.ndarray[FTYPE_t, ndim=2] gx, 
                        np.ndarray[ITYPE_t, ndim=2] cIdx, 
                        np.ndarray[FTYPE_t, ndim=2] cW, 
                        np.ndarray[FTYPE_t, ndim=3] cGrad,
			np.ndarray[ITYPE_t, ndim=1] mIdx ):
        # inpf - float input vector - h, dxdy, patch.X0, patch.dX
        cdef int nParts = mIdx.shape[0]
        cdef int cc, idx 
        cdef int ii, jj, kk, mm
        cdef double x, r, h, l, sgn
        cdef double* hh = [inpf[0,0],inpf[0,1]]
        cdef double* pp = [0.,0.]
        cdef double* ll = [0.,0.]
        cdef double* S = [0.,0.]
        cdef double* G = [0.,0.]
        cdef int* cix = [0,0]
        
        for mm in range(nParts):
            ii = mIdx[mm]
	    
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
		
                for kk in range(2):
                    x = pp[kk] - gx[idx,kk]
                    r = fabs(x)
                    l = ll[kk]
                    h = hh[kk]
                    sgn = 1. if x>0 else -1.
		    
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
	
        return 0