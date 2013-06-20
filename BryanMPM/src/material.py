import numpy as np
import materialmodel2d as mmodel
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
    def __init__(self, props, model, dwi, shape, useCython=True):
        self.props = props
        self.dwi = dwi
        self.shape = shape
        
        try:
            self.ignoreNegJ = props['ignoreNegJ']
        except Exception:
            self.ignoreNegJ = False
            
        if useCython:
            self.util = util_c
            self.mmodel = mmodel_c
        else:
            self.util = util
            self.mmodel = mmodel
            
        self.mm = self.mmodel.MaterialModel( model, props )

    
    def getdwi( self ):
        return self.dwi


    def updateContributions( self, dw, patch ):
        dw.zeroGrid( self.dwi )
        self.shape.updateContribList( dw, patch, self.dwi )                
    
    
    def setVelocity( self, dw, v ):
        pw,pm,px = dw.getMult( ['pw','pm','px'], self.dwi )
        
        for (ii,pxi,pmi) in izip(count(),px,pm):
            if isfunction(v):
                pw[ii] = v(pxi) * pmi
            else:
                pw[ii] = v * pmi
    
            
    def setExternalLoad( self, dw, fe ):
        pfe = dw.get( 'pfe', self.dwi )
        for pfei in pfe:
            pfei = fe
                
                
    def setExternalAcceleration( self, dw, acc ):
        pfe,pm = dw.getMult( ['pfe','pm'], self.dwi )
        pfe = acc * pm


    def applyExternalLoads( self, dw, patch ):
        # Apply external loads to each material
        cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )
        
        pp = dw.get( 'pfe', self.dwi )                         # External force
        gg = dw.get( 'gfe', self.dwi )        
        self.util.integrate( cIdx, cW, pp, gg )

            
    def interpolateParticlesToGrid( self, dw, patch ):
        # Interpolate particle mass and momentum to the grid
        cIdx,cW = dw.getMult( ['cIdx','cW'], self.dwi )     

        pp = dw.get( 'pm', self.dwi )                          # Mass
        gg = dw.get( 'gm', self.dwi)
        self.util.integrate( cIdx, cW, pp, gg )
        
        pp = dw.get( 'pw', self.dwi )                          # Momentum
        gg = dw.get( 'gw', self.dwi )
        self.util.integrate( cIdx, cW, pp, gg )            

     
    def computeStressTensor( self, dw, patch ):    
        pf  = dw.get( 'pF', self.dwi )                # Deformation Gradient
        pvs = dw.get( 'pVS', self.dwi )               # Volume * Stress
        pv  = dw.get( 'pVol', self.dwi )              # Volume
        
        for (ii,pfi,pvi) in izip(count(),pf,pv):
            S,Ja = self.mm.getStress( pfi )     # Get stress and det(pf)
            pvs[ii] = S * pvi * Ja              # Stress * deformed volume     
            if not self.ignoreNegJ:
                if Ja < 0:  raise JacobianError('computeStressTensor','Neg J')            
                        
                        
    def computeInternalForce( self, dw, patch ):  
        # Compute internal body forces - integrate divergence of stress to grid
        cIdx,cGrad = dw.getMult( ['cIdx','cGrad'], self.dwi )

        pp = dw.get( 'pVS', self.dwi )                          # Stress*Volume
        gg = dw.get( 'gfi', self.dwi)
        self.util.divergence( cIdx, cGrad, pp, gg )   


    def computeAndIntegrateAcceleration( self, dw, patch, tol ):
        # Integrate grid acceleration
        dwi = self.dwi
        a_leap = 1. - (patch.it==0) * 0.5             # Initializes leap-frog
        
        gm = dw.get( 'gm', dwi )                      # Mass
        gw = dw.get( 'gw', dwi )                      # Momentum
        gfi = dw.get( 'gfi', dwi )                    # Internal Force
        gfe = dw.get( 'gfe', dwi )                    # External Force
        gv = dw.get( 'gv', dwi )                      # Velocity
        ga = dw.get( 'ga', dwi )
        
        gm[:] += tol
        gv[:] = gw/gm
        ga[:] = a_leap * (gfe+gfi)/gm
        gv[:] += ga*patch.dt
            
            
    def interpolateToParticlesAndUpdate( self, dw, patch ):         
        dwi = self.dwi
        cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], self.dwi )

        pvI = dw.get( 'pvI', dwi )
        pxI = dw.get( 'pxI', dwi )
        pGv = dw.get( 'pGv', dwi )
        ga  = dw.get( 'ga', dwi )
        gv  = dw.get( 'gv', dwi )
        
        self.util.interpolate( cIdx, cW, pvI, ga )
        self.util.interpolate( cIdx, cW, pxI, gv )
        self.util.gradient( cIdx, cGrad, pGv, gv )
        
        px = dw.get( 'px', dwi )
        pw = dw.get( 'pw', dwi )
        pm = dw.get( 'pm', dwi )        
        pF = dw.get( 'pF', dwi )
        
        pw += pvI * pm * patch.dt
        px += pxI * patch.dt
        
        self.util.dotAdd( pF, pGv*patch.dt )                # pF += (pGv*dt).pF