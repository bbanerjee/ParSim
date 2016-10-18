# cython: profile=True
# cython: cdivision=True
# cython: boundscheck = False
cimport cython
import numpy as np
cimport numpy as np
ctypedef np.double_t DTYPE_t
DTYPE = np.double
cdef extern from "math.h":
    double log(double x)
    double sqrt(double x)

def makeArray( props, modelName ):
    if( modelName == 'planeStrainNeoHookean' ):
        arr = np.zeros(3)
        arr[0] = props['modulus']
        arr[1] = props['poisson']
        arr[2] = props['density']
    elif( modelName == 'planeStrainNeoHookeanMaxStress' ):        
        arr = np.zeros(4)
        arr[0] = props['modulus']
        arr[1] = props['poisson']
        arr[2] = props['density']
        arr[3] = props['maxStress']
    else:
        arr = np.array([])
    
    return arr
    

#===============================================================================
class MaterialModel:
    # Defines material models - accessed using getStress 
    #  - actual computation done in static methods   
    # Returns stress tensor and jacobian of deformation
    def __init__(self, modelName, props):
        self.modelName = modelName               # Selects Material Model
        self.props = makeArray(props, modelName)
        
    def getStress( self, F ):
        model = getattr( self, self.modelName )
        S,Ja = model(self.props, F);    
        return (S,Ja)
        
    def changeProps( self, props ):
        self.props = makeArray(props)
    

    @staticmethod
    def planeStrainNeoHookean( np.ndarray[DTYPE_t, ndim=1] props, 
                               np.ndarray[DTYPE_t, ndim=2] F ):
        cdef np.ndarray[DTYPE_t,ndim=2] S = np.zeros([2,2], dtype=DTYPE)
        cdef double v = props[1]
        cdef double E = props[0]
        cdef double l = E * v / ((1.+v)*(1.-2.*v))
        cdef double m = 0.5 * E / (1.+v)
        cdef double Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
           
        S[0,0] = l*log(Ja)/Ja + m/Ja * (F[0,0]*F[0,0]+F[0,1]*F[0,1] - 1.)
        S[0,1] = m/Ja * (F[0,0]*F[1,0]+F[0,1]*F[1,1])
        S[1,0] = m/Ja * (F[1,0]*F[0,0]+F[1,1]*F[0,1])
        S[1,1] = l*log(Ja)/Ja + m/Ja * (F[1,0]*F[1,0]+F[1,1]*F[1,1] - 1.)
     
        return (S,Ja)
    
    
    @staticmethod
    def planeStrainNeoHookeanMaxStress( np.ndarray[DTYPE_t, ndim=1] props, 
                                        np.ndarray[DTYPE_t, ndim=2] F ):
        cdef np.ndarray[DTYPE_t,ndim=2] S = np.zeros([2,2], dtype=DTYPE)
        cdef double v = props[1]
        cdef double E = props[0]
        cdef double sMax = props[3]
        cdef double l = E * v / ((1.+v)*(1.-2.*v))
        cdef double m = 0.5 * E / (1.+v)
        cdef double Ja = F[0,0]*F[1,1] - F[1,0]*F[0,1]
        cdef double vm
            
        S[0,0] = l*log(Ja)/Ja + m/Ja * (F[0,0]*F[0,0]+F[0,1]*F[0,1] - 1.)
        S[0,1] = m/Ja * (F[0,0]*F[1,0]+F[0,1]*F[1,1])
        S[1,0] = m/Ja * (F[1,0]*F[0,0]+F[1,1]*F[0,1])
        S[1,1] = l*log(Ja)/Ja + m/Ja * (F[1,0]*F[1,0]+F[1,1]*F[1,1] - 1.)
        
        vm = sqrt( S[0,0]*S[0,0]-S[0,0]*S[1,1]+S[1,1]*S[1,1]+3.*S[1,0]*S[0,1] )     
        if vm > sMax:
            S[0,0] = 0.
            S[0,1] = 0.
            S[1,0] = 0.
            S[1,1] = 0.  
            Ja = 1.
            
        return (S,Ja)