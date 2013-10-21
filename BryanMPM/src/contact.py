import numpy as np
import mpmutils as util
import skfmm

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

# Git test comment

#===============================================================================
class Contact:
    def __init__( self, dwis, patch, useCython=True ):
        self.dwis = dwis
        self.patch = patch
        self.nodes = []
        self.mtol = 1.e-10;
        if useCython:  self.util = util_c
        else:          self.util = util        

        
    def findIntersection( self, dw ):
        lvl0 = 1. - 1./self.patch.ppe
        tol0 = max(self.patch.dX)/2.
        tol = 0. #max(self.patch.dX)*0.1/self.patch.ppe
        gdi0 = dw.get('gDist',self.dwis[0])
        gdi1 = dw.get('gDist',self.dwis[1])
        gd0 = lvl0 - dw.get('gDist', self.dwis[0])
        gd1 = lvl0 - dw.get('gDist', self.dwis[1])
        sh = gd0.shape
        phi0 = gd0.reshape(self.patch.Nc)
        phi1 = gd1.reshape(self.patch.Nc)
        dist0 = skfmm.distance( phi0, self.patch.dX )
        dist1 = skfmm.distance( phi1, self.patch.dX )
        gdi0[:] = dist0.reshape(sh)
        gdi1[:] = dist1.reshape(sh)
        gmask = (dist0.reshape(sh)<tol0)*(dist1.reshape(sh)<tol0)
        gmask = gmask*((dist0.reshape(sh)+dist1.reshape(sh))<tol)
        self.nodes = np.where( gmask == True )[0]                
        
        
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
    def __init__( self, dwis, patch ):
        Contact.__init__(self, dwis, patch)
        
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
    def __init__(self, dwis, patch, bCython=True ):
        Contact.__init__(self, dwis, patch, bCython)
    
    
    def findIntersection( self, dw ):
        for dwi in self.dwis:
            cIdx,cGrad = dw.getMult( ['cIdx','cGrad'], dwi )            
            pm = dw.get( 'pm', dwi )
            pVol = dw.get( 'pVol', dwi )
            gGm = dw.get( 'gGm', dwi )
            self.util.gradscalar( cIdx, cGrad, pm, gGm )        
        
        Contact.findIntersection( self, dw )       
        
        
    def exchMomentumInterpolated( self, dw ):
        self.findIntersection( dw )       
        ii = self.nodes
        
        mr = dw.get('gm',self.dwis[0])
        ms = dw.get('gm',self.dwis[1])
        Pr = dw.get('gw',self.dwis[0])
        Ps = dw.get('gw',self.dwis[1])
        gfc = dw.get('gfc', self.dwis[1])
        
        gmr = dw.get('gGm', self.dwis[0])
        gms = dw.get('gGm', self.dwis[1])        
        nn = (vnormalize(gmr[ii])- vnormalize(gms[ii]))/2.            
        
        dp0 = 1/(mr[ii]+ms[ii])*(ms[ii]*Pr[ii]-mr[ii]*Ps[ii])
        dp = vdot(dp0,nn)
        dp = dp * (dp>0)
        
        gfc[ii] = 1.
        Pr[ii] -= dp * nn
        Ps[ii] += dp * nn

        
    def exchForceInterpolated( self, dw ):
        ii = self.nodes
        mr = dw.get('gm',self.dwis[0])
        ms = dw.get('gm',self.dwis[1])  
        
        fr = dw.get('gfe', self.dwis[0])
        fs = dw.get('gfe', self.dwis[1])
        
        fir = dw.get('gfi', self.dwis[0])
        fis = dw.get('gfi', self.dwis[1])
        
        gmr = dw.get('gGm', self.dwis[0])
        gms = dw.get('gGm', self.dwis[1])        
        nn = (vnormalize(gmr[ii])- vnormalize(gms[ii]))/2.            
        
        psi = vdot( ms[ii]*fir[ii]-mr[ii]*fis[ii], nn )
        psi = psi * (psi>0)
        
        fnor = (1/(mr[ii]+ms[ii])*psi) * nn
        fr[ii] -= vdot( fr[ii], nn ) * nn
        fs[ii] -= vdot( fs[ii], nn ) * nn
        fr[ii] -= fnor
        fs[ii] += fnor

        
#===============================================================================        
class FrictionContact(FrictionlessContact):
    def __init__(self, dwis, patch, mu, bCython=True ):
        FrictionlessContact.__init__(self, dwis, patch, bCython)
        self.mu = mu
        self.dt = patch.dt
        
    def exchForceInterpolated( self, dw ):
        ii = self.nodes
        dt = self.dt
        mu = self.mu
        
        mr = dw.get('gm',self.dwis[0])[ii] + self.mtol            # Mass
        ms = dw.get('gm',self.dwis[1])[ii] + self.mtol
        Pr = dw.get('gw',self.dwis[0])[ii]             # Momentum
        Ps = dw.get('gw',self.dwis[1])[ii]        
        fir = dw.get('gfi', self.dwis[0])[ii]          # Internal Force
        fis = dw.get('gfi', self.dwis[1])[ii]        
        fr = dw.get('gfe', self.dwis[0])               # External Force
        fs = dw.get('gfe', self.dwis[1])
        vr = Pr/mr                                     # Velocity
        vs = Ps/ms
        gmr = dw.get('gGm', self.dwis[0])
        gms = dw.get('gGm', self.dwis[1])
        nn = (vnormalize(gmr[ii])- vnormalize(gms[ii]))/2.               
                
        tt0 = (vr-vs)-vdot(vr-vs,nn)*nn
        tt = tt0/(vnorm(tt0)+self.mtol)                # Velocity Tangent
        psi = vdot( ms*fir-mr*fis, nn )
        psi = psi * (psi>0)
                
        fnor = (1/(mr+ms)*psi) * nn                    # Normal Force
        fr[ii] -= vdot( fr[ii], nn ) * nn
        fs[ii] -= vdot( fs[ii], nn ) * nn       
        fr[ii] -= fnor
        fs[ii] += fnor        
        ftan = vdot((ms*Pr-mr*Ps)+(ms*fir-mr*fis)*dt,tt) / ((mr+ms)*dt)
        ffricr = vmin(mu*vnorm(fnor)+vdot(fr[ii],tt), vnorm(ftan)) * tt
        ffrics = vmin(mu*vnorm(fnor)+vdot(fs[ii],tt), vnorm(ftan)) * tt
        fr[ii] -= ffricr
        fs[ii] += ffrics


class FrictionContactTest(FrictionContact):
    def __init__(self, dwis, patch, mu, bCython=True ):
            FrictionContact.__init__(self, dwis, patch, mu, bCython)    

    def findIntersection( self, dw ):
            # Assumes all materials share a common grid
            gm0 = dw.get('gm',self.dwis[0])
            gm1 = dw.get('gm',self.dwis[1])
            self.nodes = np.where( (gm0>self.mtol)*(gm1>self.mtol) == True )[0]

#===============================================================================
class VelocityContact(FrictionlessContact):
    def __init__(self, dwis, patch, bCython=True ):
        FrictionlessContact.__init__(self, dwis, patch, bCython)
    
    
    def findIntersection( self, dw ):
        for dwi in self.dwis:
            cIdx,cGrad = dw.getMult( ['cIdx','cGrad'], dwi )            
            pm = dw.get( 'pm', dwi )
            pVol = dw.get( 'pVol', dwi )
            gGm = dw.get( 'gGm', dwi )
            self.util.gradscalar( cIdx, cGrad, pm, gGm )        
        
        Contact.findIntersectionSimple( self, dw )       
        
        
    def exchMomentumInterpolated( self, dw ):
        FrictionlessContact.exchMomentumInterpolated( self, dw )

    def exchForceInterpolated( self, dw ):
        pass
        
    def exchMomentumIntegrated( self, dw ):
        self.findIntersection( dw )       
        ii = self.nodes
        
        mr = dw.get('gm',self.dwis[0])
        ms = dw.get('gm',self.dwis[1])
        Pr = dw.get('gw',self.dwis[0])
        Ps = dw.get('gw',self.dwis[1])
        Vr = dw.get('gv',self.dwis[0])
        Vs = dw.get('gv',self.dwis[1])
        gfc = dw.get('gfc', self.dwis[1])
        
        gmr = dw.get('gGm', self.dwis[0])
        gms = dw.get('gGm', self.dwis[1])
        nn = (vnormalize(gmr[ii])- vnormalize(gms[ii]))/2.            
        
        dp0 = 1/(mr[ii]+ms[ii])*(ms[ii]*Pr[ii]-mr[ii]*Ps[ii])
        dp = vdot(dp0,nn)
        dp = dp * (dp>0)
        
        gfc[ii] = 1.
        Vr[ii] -= dp * nn / mr[ii]
        Vs[ii] += dp * nn / ms[ii]