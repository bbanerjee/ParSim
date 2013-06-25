import numpy as np


#===============================================================================
class Patch:
    # The computational domain - called Patch to match with Vaango/Uintah
    def __init__(self,X0,X1,Nc,nGhost,t0,tf,dt,th):
        dim = 2
        self.X0 = X0                 # Bottom corner of patch domain
        self.X1 = X1                 # Top corner of patch domain
        self.Nc = Nc+1+2*nGhost      # Vector of node counts
        self.thick = th              # Thickness
        self.nGhost = nGhost         # Number of Ghost nodes
        self.dX = (X1-X0)/(Nc)       # Cell size
        self.t = t0                  # Time
        self.tf = tf                 # Final time
        self.dt = dt                 # Time increment
        self.it = 0                  # Timestep
        self.tol = 1.e-15            # Global tolerance
        self.bcs = []

            
    def initGrid(self):
        dg = self.nGhost*self.dX
        x = np.linspace( self.X0[0]-dg[0], self.X1[0]+dg[0], self.Nc[0] )
        y = np.linspace( self.X0[1]-dg[1], self.X1[1]+dg[1], self.Nc[1] )
        xx, yy = np.meshgrid( x, y )
        gx = np.append(xx.reshape(xx.size,1), yy.reshape(yy.size,1), axis=1)
        
        return gx
        
                
    def inPatch( self, pt ):
        if (pt[0] < self.X0[0]) or (pt[1] <self.X0[1]):
            return False
        if (pt[0] > self.X1[0]) or (pt[1] >self.X1[1]):
            return False
        return True
    
    def allInPatch( self, pts ):
        for pt in pts:
            if not self.inPatch( pt ):
                return False
        return True
    
    def stepTime( self ):
        self.t += self.dt
        self.it += 1