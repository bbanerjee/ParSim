#
# The MIT License
#
# Copyright (c) 2013-2014 Callaghan Innovation, New Zealand
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#

import numpy as np


#===============================================================================
class Patch:
    # The computational domain - called Patch to match with Vaango/Uintah
    def __init__(self,X0,X1,Nc,nGhost,t0,tf,dt,th,ppe):
        dim = 2
        self.ppe = ppe               # Particles per element
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