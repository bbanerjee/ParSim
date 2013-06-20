import numpy as np
import os
from evtk.hl import pointsToVTK
from copy import deepcopy as copy

def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )


#===============================================================================	
class SaveUtil:
    def __init__(self, dt, fname):
	self.dt = dt                             # Output interval
	self.fname = fname                       # Output file name
	self.nzeros = 4                          # Number of digits in filename

	try:  os.mkdir( ddir )	    
	except Exception:  pass
    

    def saveData( self, idx, matlist, dw ):
	self.save( matlist, idx, dw )
	return idx + 1
	
	
    def save( self, matlist, idx, dw ):	
	fName = self.fname + str(idx).zfill(self.nzeros)
	
	x = [];	   y = [];    v = [];    ms = [];
	matid = []
	
	for mat in matlist:
	    dwi = mat.dwi
	    px,pv,pVS,pVol = dw.getMult( ['px','pxI','pVS','pVol'], dwi )
	    nn = len(pVol)
	    matid += [dwi]*nn
	    
	    x += list(px[:,0])
            y += list(px[:,1]) 
	
	    vx = pv[:,0]
	    vy = pv[:,1]
	    v += list( np.sqrt( vx*vx + vy*vy ) * np.sign(vx) )
	    
	    pS = [pVS[ii]/pVol[ii] for ii in range(nn)]
	    ms += [vonMises(pp) for pp in pS]
	    
	x = np.array(x)
	y = np.array(y)
	z = np.zeros(x.shape)
	v = np.array(v)
	ms = np.array(ms)
	matid = np.array(matid)
	
	vsdat = {"vonMises":ms, "v":v, "mat":matid}
	pointsToVTK(fName, x, y, z, data = vsdat)