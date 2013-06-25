import numpy as np
import mpmutils as util

try:
    import mpmutils_c as util_c
except Exception:
    util_c = util

def vnorm( x ):
    return np.sqrt((x*x).sum(axis=1)[:,np.newaxis])

class SimpleContact:
    def __init__( self, dwis ):
        self.dwis = dwis
        self.nodes = []

    def findIntersection( self, dw ):
        self.findIntersectionSimple( dw )
    
    def findIntersectionSimple( self, dw ):
        # Assumes all materials share a common grid
        gm0 = dw.get('gm',self.dwis[0])
        gm1 = dw.get('gm',self.dwis[1])
        self.nodes = np.where( (gm0>0)*(gm1>0) == True )[0]
    
    def exchMomentumInterpolated( self, dw ):
        pass
                 
    def exchForceInterpolated( self, dw ):
        pass
    
    def exchMomentumIntegrated( self, dw ):
        pass        
    

class FreeContact(SimpleContact):
    def __init__( self, dwis ):
        SimpleContact.__init__(self, dwis)
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersection( dw )
        if self.nodes.any():
            self.exchVals( 'gm', dw )
            self.exchVals( 'gw', dw )
         
    def exchForceInterpolated( self, dw ):
        if self.nodes.any():
            self.exchVals( 'gfi', dw )

    
    def exchVals( self, lbl, dw ):
        g0 = dw.get(lbl,self.dwis[0])
        g1 = dw.get(lbl,self.dwis[1])
       
        g0[self.nodes] += g1[self.nodes]
        g1[self.nodes] = g0[self.nodes]

        
class FrictionlessContact(SimpleContact):
    def __init__(self, dwis, bCython=True ):
        SimpleContact.__init__(self, dwis)
        if useCython:  self.util = util_c
        else:          self.util = util
        
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersection( dw )       
        
        for dwi in self.dwis:
            pm = dw.get( 'pm', dwi )
            pVol = dw.get( 'pVol', dwi )
            gGm = dw.get( 'gGm', dwi )
            self.util.gradscalar( cIdx, cGrad, pm, gGm )
            
        gGm0 = dw.get('gGm',dwis[0])
        gGm1 = dw.get('gGm',dwis[1])
        
        nGm0 = vnorm( gGm0[self.nodes] )
        nGm1 = vnorm( gGm1[self.nodes] )
        bGm0 = nGm0 > nGm1
        
        gGm0[self.nodes] = bGm0 * gGm0[self.nodes] + (bGm0-1)*gGm1[self.nodes]
        gGm1[self.nodes] = -gGm0[self.nodes]
        