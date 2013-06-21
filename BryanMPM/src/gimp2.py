import numpy as np
    
def uSG( x, h, l ):
    r = abs(x)
    sgnx = cmp(x,-x)
	
    if (r<l):      S = 1. - (r*r+l*l) / (2.*h*l);  G = -x/(h*l)
    elif(r<h-l):   S = 1. - r/h;                   G = -sgnx/h
    elif(r<h+l):   S = (h+l-r)*(h+l-r) / (4.*h*l); G = (h+l-r) / (-2.*sgnx*h*l)
    else:          S = G = 0.
    return (S,G)
	

def getCell( patch, pos ):    
    # Gets lower left node of 4-cell block
    x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
    idx = np.floor(x_sc)
    rem = (x_sc - 1.*idx) >= 0.5
    ii = idx[0] if rem[0] else idx[0]-1
    jj = idx[1] if rem[1] else idx[1]-1
	
    return int(jj * patch.Nc[0] + ii)
    

def updateContribList( dw, patch, dwi ):
    # Update node contribution list
    nx = patch.Nc[0]
    h = patch.dX
    dxdy = h[::-1]/h
    idxs = [0,1,2,nx,nx+1,nx+2,2*nx,2*nx+1,2*nx+2]
    S = np.zeros(h.size)
    G = np.zeros(h.size)	
    
    cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], dwi )
    px,gx,pVol,pF = dw.getMult( ['px','gx','pVol','pF'], dwi )

    for ii in range(len(pVol)):
        cc = getCell( patch, px[ii] )	           
        l = np.sqrt(pVol[ii]/(4.*patch.thick*dxdy)) * np.diag(pF[ii])

        for jj in range(9):	
            idx = idxs[jj] + cc 
            r = px[ii] - gx[idx]	
		
            for kk in range(len(r)):
                S[kk],G[kk] = uSG( r[kk], h[kk], l[kk] )
		
            cIdx[ii][jj] = idx
            cW[ii][jj] = S[0]*S[1]
            cGrad[ii][jj] = G * S[::-1] 