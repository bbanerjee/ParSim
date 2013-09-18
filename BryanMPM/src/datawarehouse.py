import numpy as np
import collections
from itertools import izip
import os
from saveutil import SaveUtil
from evtk.hl import pointsToVTK
from copy import deepcopy as copy

def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )
def toArray( val ):
    if type(val) is np.ndarray:
	return val
    else:
	return np.array(val)

#===============================================================================	
class DataWarehouse:
    # Holds all the data particle and node for an individual 
    # timestep used for data access
    def __init__(self, ddir, fname, dt, idx=0, t=0. ):
	self.dw = dict()
	self.out_idx = idx                       # Index for output file
	self.idx = idx                           # Iteration index
	self.t = t                               # Initial time
	self.dt = dt
	self.fname = ddir + '/' + fname          # Output file name
	self.nzeros = 4                          # Number of digits in filename
	self.saveUtil = SaveUtil(dt,self.fname)  # Saving Utility Class
	
	try:  os.mkdir( ddir )	    
	except Exception:  pass	

    def saveData( self, dt, matlist ):
	if self.checkSave( dt ):
	    self.out_idx = self.saveUtil.saveData( self.out_idx, matlist, self )
	self.t += dt
	self.idx += 1
	
    def dumpData( self, dt, matlist ):
	if self.checkSave( dt ):
	    self.out_idx = self.saveUtil.dumpData( self.out_idx, matlist, self )
	self.t += dt
	self.idx += 1	
	
    def checkSave( self, dt ):
	dr = self.t/self.dt
	dt0 = self.dt * min( dr-np.floor(dr), np.ceil(dr)-dr )
	return dt0 < dt/2.

    def init( self, label, dwi, val ):
	self.dw[label,dwi] = toArray(val)
	
    def append( self, label, dwi, val ):
	self.dw[label,dwi] = np.append( self.dw[label,dwi], toArray(val), axis=0 )
	
    def add( self, label, dwi, val ):
	if len(self.get(label,dwi)): self.append( label, dwi, val )
	else: self.init( label, dwi, val )
        
    def zero( self, label, dwi ):
	self.dw[label,dwi] *= 0.    
	
    def get( self, label, dwi ):
	try:
	    return self.dw[label,dwi]
	except Exception:
	    return []
	    
    def getMult( self, labels, dwi ):
	output = []
	for label in labels:
	    output.append( self.get( label, dwi ) )
	return output    
      
    def addParticles( self, dwi, pX, pVol, pN, density, shSize ):
	npt = len(pX)
	labels = ['pw','pvI','pxI','pfe','pGv','pVS','pfi','pfc','pwc']
	sh0 = (npt,2)
	shapes = [sh0,sh0,sh0,sh0,(npt,2,2),(npt,2,2),sh0,sh0,sh0,sh0]

	# Add initial position, position, volume, mass, and node contributions
	self.add( 'pX',    dwi, pX )
	self.add( 'px',    dwi, pX )
	self.add( 'pN',    dwi, pN )
	self.add( 'pn',    dwi, pN )
	self.add( 'pVol',  dwi, pVol*np.ones(npt) )
	self.add( 'pm',    dwi, pVol*density*np.ones((npt,2)) )
	self.add( 'cIdx',  dwi, np.zeros((npt,shSize),dtype=np.int) )
	self.add( 'cW',    dwi, np.zeros((npt,shSize)) )
	self.add( 'cGrad', dwi, np.zeros((npt,shSize,2)) )		

	# Add variables contained in "labels"
	for (lbl,shp) in izip( labels, shapes ):
	    self.add( lbl, dwi, np.zeros(shp) )
	
	# Create and add initial deformation (identity matrix)
	pF = np.zeros((npt,2,2))
	for pFi in pF:
	    pFi += np.eye(2)
	self.add( 'pF', dwi, pF )

    def createGrid( self, dwi, patch ):
	gx = patch.initGrid()
	self.dw['gx',dwi] = toArray(gx)
	self.zeroGrid(dwi)
	
    def zeroGrid( self, dwi ):
	gx = self.get('gx',dwi)
	labels = ['gm','gv','gw','ga','gfe','gfi','gn','gfc','gwc','gGm']
	for label in labels:
	    self.init( label, dwi, np.zeros(gx.shape) )
	    
	self.init('gDist', dwi, np.zeros(len(gx)) )