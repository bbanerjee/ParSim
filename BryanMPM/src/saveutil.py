import numpy as np
import os
import time
from evtk.hl import pointsToVTK
from copy import deepcopy as copy

def vonMises( S ):
    return np.sqrt( S[0,0]*S[0,0] - S[0,0]*S[1,1] + S[1,1]*S[1,1] +
                    3.*S[1,0]*S[0,1] )

def vnorm( x ):
    return np.sqrt((x*x).sum(axis=1))


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
	

    def dumpData( self, idx, matlist, dw ):
	self.dump( matlist, idx, dw )
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
	
    def dump( self, matlist, idx, dw ):	
	fName = self.fname + str(idx).zfill(self.nzeros)
	fNodeName = self.fname + 'nodes_' + str(idx).zfill(self.nzeros)

        # Standard Save	
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
	pdat = {"vonMises":ms, "v":v, "mat":matid}
	
	# Create Node Grid
	gX = dw.get('gx',matlist[0].dwi)
	gx = np.array(gX[:,0])
	gy = np.array(gX[:,1])
	gz = np.zeros(gx.shape)
	gdat = {"g_test":gz}
	
	partList = ['pX','pxI','pvI','pw','pfi','pfe','pfc','pwc','pn','pm']
	nodeList = ['gv','gw','ga','gfe','gfi','gn','gfc','gwc','gDist','gm']
	
	for var in partList:
	    tmp1 = []
	    tmp2 = []
	    tmp3 = []
	    for mat in matlist:
		dwi = mat.dwi
		pvar = dw.get( var, dwi )
		tmp1 += list(pvar[:,0])
		tmp2 += list(pvar[:,1])
		tmp3 += list(vnorm(pvar))
	    pdat[var+'_x'] = np.array(tmp1)
	    pdat[var+'_y'] = np.array(tmp2)
	    pdat[var] = np.array(tmp3)
	    
	for var in nodeList:
	    for mat in matlist:
		dwi = mat.dwi
		gvar = dw.get( var, dwi )
		try:
		    gdat[var+str(dwi)+'_x'] = np.array(list(gvar[:,0]))
		    gdat[var+str(dwi)+'_y'] = np.array(list(gvar[:,1]))
		    gdat[var+str(dwi)] = vnorm(gvar)	
		except Exception:
		    gdat[var+str(dwi)] = gvar		    
		    
	pointsToVTK(fName, x, y, z, data = pdat)
	pointsToVTK(fNodeName, gx, gy, gz, data = gdat)	