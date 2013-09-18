# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
cimport cython
import numpy as np
cimport numpy as np
ctypedef np.int_t ITYPE_t
ctypedef np.double_t FTYPE_t
ITYPE = np.int
FTYPE = np.float


# Utils for moving data between particles and grid
def integrate( np.ndarray[ITYPE_t, ndim=2] cIdx, 
               np.ndarray[FTYPE_t, ndim=2] cW, 
               np.ndarray[FTYPE_t, ndim=2] pp,
               np.ndarray[FTYPE_t, ndim=2] gg ):
    # Integrate particle values to grid (p->g)
    cdef int ii, jj, ixc
    cdef double w
    cdef int nParts = pp.shape[0]
    cdef int nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj];  
            w = cW[ii,jj]
            gg[ixc,0] += pp[ii,0] * w
            gg[ixc,1] += pp[ii,1] * w
            

def interpolate( np.ndarray[ITYPE_t, ndim=2] cIdx, 
                 np.ndarray[FTYPE_t, ndim=2] cW, 
                 np.ndarray[FTYPE_t, ndim=2] pp,
                 np.ndarray[FTYPE_t, ndim=2] gg ):
    # Interpolate grid values to particles pp
    cdef int ii, jj, ixc
    cdef double w
    cdef int nParts = pp.shape[0]
    cdef int nContrib = cIdx.shape[1]
    for ii in range(nParts):
        pp[ii,0] = 0.0
        pp[ii,1] = 0.0
        for jj in range(nContrib):
            ixc = cIdx[ii,jj];  
            w = cW[ii,jj]
            pp[ii,0] += gg[ixc,0] * w; 
            pp[ii,1] += gg[ixc,1] * w; 
    return 0       


def gradient( np.ndarray[ITYPE_t, ndim=2] cIdx, 
              np.ndarray[FTYPE_t, ndim=3] cGrad, 
              np.ndarray[FTYPE_t, ndim=3] pp,
              np.ndarray[FTYPE_t, ndim=2] gg ):
    # Interpolate grid values to particles pp
    cdef int ii, jj, ixc
    cdef double cg0, cg1
    cdef int nParts = pp.shape[0]
    cdef int nContrib = cIdx.shape[1]
    for ii in range(nParts):
        pp[ii,0,0] = 0.0
        pp[ii,0,1] = 0.0
        pp[ii,1,0] = 0.0
        pp[ii,1,1] = 0.0
        for jj in range(nContrib):
            ixc = cIdx[ii,jj];  
            cg0 = cGrad[ii,jj,0]
            cg1 = cGrad[ii,jj,1]
            pp[ii,0,0] += gg[ixc,0] * cg0; 
            pp[ii,0,1] += gg[ixc,0] * cg1; 
            pp[ii,1,0] += gg[ixc,1] * cg0; 
            pp[ii,1,1] += gg[ixc,1] * cg1; 
    return 0       


def divergence( np.ndarray[ITYPE_t, ndim=2] cIdx, 
              np.ndarray[FTYPE_t, ndim=3] cGrad, 
              np.ndarray[FTYPE_t, ndim=3] pp,
              np.ndarray[FTYPE_t, ndim=2] gg ):
    # Send divergence of particle field to the grid
    cdef int ii, jj, ixc
    cdef double cg0, cg1
    cdef int nParts = pp.shape[0]
    cdef int nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj];  
            cg0 = cGrad[ii,jj,0]
            cg1 = cGrad[ii,jj,1]
            gg[ixc,0] -= pp[ii,0,0] * cg0 + pp[ii,1,0] * cg1; 
            gg[ixc,1] -= pp[ii,1,0] * cg0 + pp[ii,1,1] * cg1; 
    return 0      

def gradscalar( np.ndarray[ITYPE_t, ndim=2] cIdx, 
              np.ndarray[FTYPE_t, ndim=3] cGrad, 
              np.ndarray[FTYPE_t, ndim=2] pp,
              np.ndarray[FTYPE_t, ndim=2] gg ):
    # Send gradient of scalar particle field to grid
    cdef int ii, jj, ixc
    cdef int nParts = pp.shape[0]
    cdef int nContrib = cIdx.shape[1]
    for ii in range(nParts):
        for jj in range(nContrib):
            ixc = cIdx[ii,jj];  
            gg[ixc,0] += pp[ii,0] * cGrad[ii,jj,0];
            gg[ixc,1] += pp[ii,0] * cGrad[ii,jj,1];
    return 0      

def dotAdd( np.ndarray[FTYPE_t, ndim=3] pp,
            np.ndarray[FTYPE_t, ndim=3] qq ):
    # return pp += qq dot pp
    cdef int ii, jj, kk
    cdef int nParts = pp.shape[0]
    cdef np.ndarray[FTYPE_t,ndim=2] dot = np.zeros([2,2], dtype=FTYPE)

    for ii in range(nParts):
        for jj in range(2):
            for kk in range(2):
                dot[jj,kk] = qq[ii,jj,0]*pp[ii,0,kk]+qq[ii,jj,1]*pp[ii,1,kk]
                
                
        for jj in range(2):
            for kk in range(2):
                pp[ii,jj,kk] += dot[jj,kk]
    return 0      


def readableTime( tt ):
    from dateutil.relativedelta import relativedelta as reldelta
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']    
    h_read = lambda delta: ['%d %s' % (getattr(delta, attr), 
                         getattr(delta, attr) > 1 and attr or attr[:-1]) 
    for attr in attrs if getattr(delta, attr)]
    tm = h_read(reldelta(seconds=tt))
    
    st = ''
    for t in tm:
        st += t + ', '
    return st[:-2]