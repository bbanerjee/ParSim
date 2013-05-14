import numpy as np
import collections
import os
from evtk.hl import pointsToVTK
from copy import deepcopy as copy

def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )


#===============================================================================	
class DataWarehouse:
    # Holds all the data particle and node for an individual 
    # timestep used for data access
    def __init__(self, savedt, fname, ddir='.', 
                 saveidx=0, t=0., idx=0):
	self.t = t                               # Current time
	self.idx = idx                           # Time index
	self.saveidx = saveidx                   # Index for output file
	self.savedt = savedt                     # Output interval
	self.ddir = ddir                         # Output directory
	self.fname = fname                       # Output file name
	self.nzeros = 4                          # Number of digits in filename

	try:  os.mkdir( ddir )	    
	except Exception:  pass
	    
	# Particle variable lists
	self.pX   = []     # Initial Position
	self.px   = []     # Position	
	self.pMat = []     # Material ID
	self.pm   = []     # Mass
	self.pVol = []     # Initial Volume
	        
	# Node variable lists
	self.nNodes = 0      # Number of nodes
        self.gx = []         # Position	
        
        
    def initArrays( self, shSize ):
	npt = len(self.pX)    
	ng  = len(self.gx)
	
	self.pX = np.array(self.pX)
	self.px = np.array(self.px)
	self.pm = np.array(self.pm)
	self.pMat = np.array(self.pMat)
	self.pVol = np.array(self.pVol)
	self.gx = np.array(self.gx)
	
	# Particle Variables
	self.pw  = np.zeros((npt,2))           # Momentum
	self.pvI = np.zeros((npt,2))           # Velocity Increment
	self.pxI = np.zeros((npt,2))           # Position Increment
	self.pfe = np.zeros((npt,2))           # External Force
	self.pGv = np.zeros((npt,2,2))         # Velocity Gradient
	self.pVS = np.zeros((npt,2,2))         # Stress * Volume
	self.pF  = np.zeros((npt,2,2))         # Deformation gradient
	for ii in range(len(self.pF)): 
	    self.pF[ii] = np.diag(np.ones(2))
	    
	# Node Variables
	self.gm  = np.zeros((ng,2))           # Mass
	self.gv  = np.zeros((ng,2))           # Velocity
	self.gw  = np.zeros((ng,2))           # Momentum
	self.ga  = np.zeros((ng,2))           # Acceleration
	self.gfe = np.zeros((ng,2))           # External Force
	self.gfi = np.zeros((ng,2))           # Internal Force	  
	
	# Contribution Variables
	self.cIdx  = np.zeros((npt,shSize),dtype=np.int)
	self.cW    = np.zeros((npt,shSize))
	self.cGrad = np.zeros((npt,shSize,2))
	
	
    def getData( self, name ):
	source = getattr(self, name)
	return source

    
    def getMatIndex( self, matid ):
	matIdx = np.where( self.pMat==matid )[0]
	return matIdx	
 
    
    def addParticle( self, mat, pos, mass, vol ):
	self.pX.append( pos.copy() )              # Initial Position
	self.px.append( pos.copy() )              # Position
	self.pMat.append( mat )                   # Material ID
	self.pm.append( np.array([mass,mass]) )   # Mass
	self.pVol.append( vol )                   # Initial Volume	
	
	
    def addNode( self, pos ):
	self.gx.append( pos )
	

    def resetNodes( self ):
	#  Reset nodal variables (other than position) to zero
	self.gm  *= 0.0;    self.gv  *= 0.0
	self.gw  *= 0.0;    self.gfe *= 0.0
	self.gfi *= 0.0;    self.ga  *= 0.0   		


    def saveDataAndAdvance( self, dt ):
	rem = self.t % self.savedt / self.savedt
	tol = dt/2./self.savedt
	if (rem < tol) or ((1-rem) < tol):
	    self.saveData()
	    self.saveidx += 1
	
	self.t += dt
	self.idx += 1	

	
    def saveData( self ):	
	strIdx = str(self.saveidx).zfill(self.nzeros)
	fNamePoint = ( self.ddir + '/' + self.fname + 
	               str(self.saveidx).zfill(self.nzeros) )
	fNameNode = ( self.ddir + '/' + self.fname + "_n" + 
	              str(self.saveidx).zfill(self.nzeros) )
	
	self.savePointData(fNamePoint)
	self.saveNodeData(fNameNode)		

    
    def savePointData( self, fName ):
	px = copy(self.px[:,0])
	py = copy(self.px[:,1])
	pz = np.zeros(px.shape)
	vx = self.pxI[:][:,0]
	vy = self.pxI[:][:,1]
	vv = np.sqrt( vx*vx + vy*vy ) * np.sign(vx)
	pS = [self.pVS[ii]/self.pVol[ii] for ii in range(len(self.pVol))]
	mises = np.array([vonMises(ii) for ii in pS])
	
	vsdat = {"vonMises":mises, "v":vv, "mat":self.pMat}
	pointsToVTK(fName, px, py, pz, data = vsdat)
	
	
    def saveNodeData( self, fName ):
	gx = copy(self.gx[:][:,0])
	gy = copy(self.gx[:][:,1])
	gz = np.zeros(gx.shape)
	ax = copy(self.ga[:][:,0])
	ay = copy(self.ga[:][:,1])
	aa = np.sqrt( ax*ax + ay*ay )    
		
	vsdat = {"ax":ax, "ay":ay, "aa":aa}
	pointsToVTK(fName, gx, gy, gz, data = vsdat)	