import numpy as np
import Material
import cProfile
from itertools import izip, count
from inspect import isfunction
import mpmutils as util

try:
    import materialmodel2d_c as mmodel_c
    import mpmutils_c as util_c
except Exception:
    mmodel_c = mmodel
    util_c = util


#===============================================================================
class Material_Stationary(Material):
    # Stationary Rigid Material - 
    def __init__(self, props, model, dwi, shape, useCython=True):
        self.props = props
        self.dwi = dwi
        self.shape = shape
                       
        if useCython:
            self.util = util_c
            self.mmodel = mmodel_c
        else:
            self.util = util
            self.mmodel = mmodel
    
    def getdwi( self ):
        return self.dwi


    def updateContributions( self, dw, patch ):
        dw.zeroGrid( self.dwi )
        self.shape.updateContribList( dw, patch, self.dwi )                
    
    
    def setVelocity( self, dw, v ):
        pass
            
    def setExternalLoad( self, dw, fe ):
        pass
                
                
    def setExternalAcceleration( self, dw, acc ):
        pass


    def applyExternalLoads( self, dw, patch ):
        # Apply external loads to each material
        pass

            
    def interpolateParticlesToGrid( self, dw, patch ):
        # Interpolate particle mass and momentum to the grid
        cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )     

        pp = dw.get( 'pm', self.dwi )                          # Mass
        gg = dw.get( 'gm', self.dwi)
        self.util.integrate( cIdx, cW, pp, gg )

     
    def computeStressTensor( self, dw, patch ):    
        pass 
                        
                        
    def computeInternalForce( self, dw, patch ):  
        # Compute internal body forces - integrate divergence of stress to grid
        pass


    def computeAndIntegrateAcceleration( self, dw, patch, tol ):
        # Integrate grid acceleration
        pass
            
            
    def interpolateToParticlesAndUpdate( self, dw, patch ):         
        pass
    