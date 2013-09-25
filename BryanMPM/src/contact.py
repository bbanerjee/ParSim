import numpy as np
import mpmutils as util
from itertools import izip, count

try:
    import mpmutils_c as util_c
except Exception:
    util_c = util

def vnorm( x ):
    return np.sqrt((x*x).sum(axis=1)[:,np.newaxis])

def vdot( x, y ):
    return (x*y).sum(axis=1)[:,np.newaxis]

def vmin( x, y ):
    return x*((x-y)<0) + y*((x-y)>=0)

def vnormalize( x ):
    xx = vnorm(x)
    xnorm = (xx>0)*xx + (xx==0)*1.
    return x/xnorm


#===============================================================================
class Contact:
    def __init__( self, dwis, ppe=1, useCython=True ):
        self.dwis = dwis
        self.ppe = ppe
        self.nodes = []
        self.mtol = 1.e-16;

        if useCython:  self.util = util_c
        else:          self.util = util        
    
    def findIntersectionSimple( self, dw ):
        # Assumes all materials share a common grid
        gm0 = dw.get('gm',self.dwis[0])
        gm1 = dw.get('gm',self.dwis[1])
        self.nodes = np.where( (gm0>self.mtol)*(gm1>self.mtol) == True )[0]
    
    def exchMomentumInterpolated( self, dw ):
        pass
                 
    def exchForceInterpolated( self, dw ):
        pass
    
    def exchMomentumIntegrated( self, dw ):
        pass        
    
#===============================================================================
class FreeContact(Contact):
    def __init__( self, dwis, nx=1, ppe=1 ):
        Contact.__init__(self, dwis)
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersectionSimple( dw )
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


#===============================================================================
class FrictionlessContact(Contact):
    # See Pan et al - 3D Multi-Mesh MPM for Solving Collision Problems
    def __init__(self, dwis, nx, ppe=1, bCython=True ):
        Contact.__init__(self, dwis, ppe)
        self.nx = nx
    
    
    def findIntersection( self, dw ):
        self.nodes = np.array([],int)
        nodes = []
        #idxs = [-self.nx-1,-self.nx,-self.nx+1,-1,0,1,self.nx-1,self.nx,self.nx+1]      
        idxs = [0]

        lvl0 = 1. - 1./self.ppe
        gd0 = dw.get('gDist', self.dwis[0]) - lvl0
        gd1 = dw.get('gDist', self.dwis[1]) - lvl0
        gmask = (gd0>-lvl0)*(gd1>-lvl0)*((gd0+gd1)>(-0.5/self.ppe))
        nodes0 = np.where( gmask == True )[0]

        if nodes0.any():
            for idx in idxs:
                nodes += list( nodes0 + idx )  
            self.nodes = np.array(list(set(nodes)))        
                
        for dwi in self.dwis:
            cIdx,cW = dw.getMult( ['cIdx','cW'], dwi )            
            pn = dw.get( 'pn', dwi )
            gn = dw.get( 'gn', dwi )
            self.util.integrate( cIdx, cW, pn, gn )  
                    
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersection( dw )       
        ii = self.nodes
        
        mr = dw.get('gm',self.dwis[0])
        ms = dw.get('gm',self.dwis[1])
        Pr = dw.get('gw',self.dwis[0])
        Ps = dw.get('gw',self.dwis[1])
        Prc = dw.get('gwc',self.dwis[0])
        Psc = dw.get('gwc',self.dwis[1])
        nr = dw.get('gn', self.dwis[0])
        ns = dw.get('gn', self.dwis[1])
        nn = (vnormalize(nr[ii])- vnormalize(ns[ii]))/2.    
        #nn = vnormalize(nr[ii])    
        dp0 = 1/(mr[ii]+ms[ii]+self.mtol)*(ms[ii]*Pr[ii]-mr[ii]*Ps[ii])
        dp = vdot(dp0,nn)
        dp = dp * (dp>=0.)
        
        Prc[ii] = -dp * nn
        Psc[ii] = dp * nn

        
    def exchForceInterpolated( self, dw ):
        ii = self.nodes
        mr = dw.get('gm',self.dwis[0])
        ms = dw.get('gm',self.dwis[1])  
        fr = dw.get('gfe', self.dwis[0])
        fs = dw.get('gfe', self.dwis[1])
        fir = dw.get('gfi', self.dwis[0])
        fis = dw.get('gfi', self.dwis[1])
        
        nr = dw.get('gn', self.dwis[0])
        ns = dw.get('gn', self.dwis[1])
        nn = (vnormalize(nr[ii])- vnormalize(ns[ii]))/2.    
        #nn = vnormalize(nr[ii])    
        psi = vdot( ms[ii]*fir[ii]-mr[ii]*fis[ii], nn )
        psi = psi * (psi>=0.)
        
        fnor = (1./(mr[ii]+ms[ii]+self.mtol)*psi) * nn
        fr[ii] = -fnor
        fs[ii] = +fnor
        
        
class FrictionContact(FrictionlessContact):
    def __init__(self, dwis, mu, dt, nx, ppe=1, bCython=True ):
        FrictionlessContact.__init__(self, dwis, nx, ppe)
        self.mu = mu
        self.dt = dt
        
    def exchForceInterpolated( self, dw ):
        ii = self.nodes
        dt = self.dt
        mu = self.mu
        
        mr = dw.get('gm',self.dwis[0])[ii]             # Mass
        ms = dw.get('gm',self.dwis[1])[ii]
        Pr = dw.get('gw',self.dwis[0])[ii]             # Momentum
        Ps = dw.get('gw',self.dwis[1])[ii]        
        fir = dw.get('gfi', self.dwis[0])[ii]          # Internal Force
        fis = dw.get('gfi', self.dwis[1])[ii]        
        fr = dw.get('gfe', self.dwis[0])               # External Force
        fs = dw.get('gfe', self.dwis[1])
        gm = dw.get('gGm', self.dwis[0])[ii]           # Mass Gradient
        vr = Pr/mr                                     # Velocity
        vs = Ps/ms
                
        nn = gm/vnorm(gm)                              # Surface Normal
        tt0 = (vr-vs)-vdot(vr-vs,nn)*nn
        tt = tt0/vnorm(tt0)                            # Velocity Tangent
                
        psi = vdot( ms*fir-mr*fis, nn )
        psi = psi * (psi>0)
                
        fnor = (1/(mr+ms+self.mtol)*psi) * nn          # Normal Force
        fr[ii] -= fnor
        fs[ii] += fnor        
        
        ftan = vdot((ms*Pr-mr*Ps)+(ms*fir-mr*fis)*dt,tt)/((mr+ms+self.mtol)*dt)
        ffric = vmin(mu*vnorm(fnor), vnorm(ftan)) * tt
        
        fr[ii] -= ffric
        fs[ii] += ffric