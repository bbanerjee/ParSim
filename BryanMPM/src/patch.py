import numpy as np


#===============================================================================
class Patch:
    # The computational domain - called Patch to match with Vaango/Uintah
    def __init__(self,X0,X1,Nc,nGhost,t0,tf,dt,th,dw=0):
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
        
        if not (dw==0):
            self.initGrid(dw)        # If specified, initialize nodes in data
                                     # warehouse 
            
    def initGrid(self, dw):
        for jj in range(self.Nc[1]):
            yy = (jj-self.nGhost)*self.dX[1] + self.X0[1]
            for ii in range(self.Nc[0]):
                xx = (ii-self.nGhost)*self.dX[0] + self.X0[0]
                XX = np.array( [xx, yy] )
                dw.addNode( XX )
                
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