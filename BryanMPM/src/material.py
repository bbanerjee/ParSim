import numpy as np
import materialmodel2d as mmodel
import cProfile
import mpmutils as util

try:
    import materialmodel2d_c as mmodel_c
    import mpmutils_c as util_c
except Exception:
    mmodel_c = mmodel
    util_c = util


#===============================================================================
class Error(Exception):
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg    

class JacobianError(Error):
    def __init__(self, expr, msg):
        Error.__init__(self,expr,msg)

#===============================================================================
class Material:
    # Material - holds update functions - default is deformable
    # overridden by RigidMaterial for rigid materials
    def __init__(self, matid, props, model, useCython=True, ignoreNegJ=False):
        self.matid = matid
        self.props = props
        self.hasParts = False
        self.pIdx = []   
        self.ignoreNegJ = ignoreNegJ
        if useCython:
            self.util = util_c
            self.mmodel = mmodel_c
        else:
            self.util = util
            self.mmodel = mmodel
            
        self.mm = self.mmodel.MaterialModel( model, props )
        

    def getParticles( self, dw, patch, sh ):
        self.pIdx = dw.getMatIndex( self.matid )
        self.hasParts = ( len(self.pIdx) > 0 )
        if self.hasParts:
            sh.updateContribList( dw, patch, self.pIdx )
            
            
    def setVelocity( self, dw, v ):
        pw = dw.getData( 'pw' )
        pm = dw.getData( 'pm' )
        for ii in self.pIdx:
            pw[ii] = v * pm[ii]

        
    def setExternalLoad( self, dw, fe ):
        pfe = dw.getData( 'pfe' )
        for ii in self.pIdx: 
            pfe[ii] = fe
                
                
    def setExternalAcceleration( self, dw, acc ):
        pfe = dw.getData( 'pfe' )
        pm = dw.getData( 'pm' )        
        for ii in self.pIdx: 
            pfe[ii] = acc * pm[ii]


    def applyExternalLoads( self, dw, patch ):
        # Apply external loads to each material
        pp = dw.getData( 'pfe' )                         # External force
        gg = dw.getData( 'gfe')
        self.util.integrate( dw.cIdx, dw.cW, pp, gg, self.pIdx )

            
    def interpolateParticlesToGrid( self, dw, patch ):
        # Interpolate particle mass and momentum to the grid
        pp = dw.getData( 'pm' )                          # Mass
        gg = dw.getData( 'gm')
        self.util.integrate( dw.cIdx, dw.cW, pp, gg, self.pIdx )
        
        pp = dw.getData( 'pw' )                          # Momentum
        gg = dw.getData( 'gw')
        self.util.integrate( dw.cIdx, dw.cW, pp, gg, self.pIdx )        
     
    def computeStressTensor( self, dw, patch ):    
        pf  = dw.getData( 'pF' )                      # Deformation Gradient
        pvs = dw.getData( 'pVS' )                     # Volume * Stress
        pv  = dw.getData( 'pVol' )                    # Volume
        for ii in self.pIdx:
            S,Ja = self.mm.getStress( pf[ii] )        # Get stress and det(pf)
            pvs[ii] = S * pv[ii] * Ja                 # Stress * deformed volume     
            if not self.ignoreNegJ:
                if Ja < 0: 
                    raise JacobianError('computeStressTensor', 
                                        'Negative Jacobian')            
                        
    def computeInternalForce( self, dw, patch ):  
        # Compute internal body forces - integrate divergence of stress to grid
        pp = dw.getData( 'pVS' )                          # Stress*Volume
        gg = dw.getData( 'gfi')
        self.util.divergence( dw.cIdx, dw.cGrad, pp, gg, self.pIdx )   
            
    def interpolateToParticlesAndUpdate( self, dw, patch ):         
        pvI = dw.getData( 'pvI' )
        pxI = dw.getData( 'pxI' )
        pGv = dw.getData( 'pGv' )
        ga  = dw.getData( 'ga' )
        gv  = dw.getData( 'gv' )
        
        self.util.interpolate( dw.cIdx, dw.cW, pvI, ga, self.pIdx )
        self.util.interpolate( dw.cIdx, dw.cW, pxI, gv, self.pIdx )
        self.util.gradient( dw.cIdx, dw.cGrad, pGv, gv, self.pIdx )
        
        px = dw.getData( 'px' )
        pw = dw.getData( 'pw' )
        pm = dw.getData( 'pm'  )        
        pF = dw.getData( 'pF' )
        
        pw[self.pIdx] += pvI[self.pIdx] * pm[self.pIdx] * patch.dt
        px[self.pIdx] += pxI[self.pIdx] * patch.dt
        
        self.util.dotAdd( pF, pGv*patch.dt, self.pIdx )  # pF += (pGv*dt).pF