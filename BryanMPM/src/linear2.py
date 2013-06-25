import numpy as np
    
def uSG( x, h ):
    r = abs(x)
    sgnx = cmp(x,-x)
    
    if ( r < h ):
        S = 1. - r/h
        G = -sgnx/h
    else: 
        S = G = 0.
    return( S,G )
	

def getCell( patch, pos ):    
    # Gets lower left node of 4-cell block
    x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
    idx = np.floor(x_sc)
    ii = idx[0] 
    jj = idx[1] 
	
    return int(jj * patch.Nc[0] + ii)
    

def updateContribList( dw, patch, dwi ):
    # Update node contribution list
    nx = patch.Nc[0]
    h = patch.dX
    idxs = [0,1,nx,nx+1]
    S = np.zeros(h.size)
    G = np.zeros(h.size)	
    
    cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], dwi )
    px,gx = dw.getMult( ['px','gx'], dwi )

    for ii in range(len(pVol)):
        cc = getCell( patch, px[ii] )	           

        for jj in range(4):	
            idx = idxs[jj] + cc 
            r = px[ii] - gx[idx]	
		
            for kk in range(len(r)):
                S[kk],G[kk] = uSG( r[kk], h[kk] )
		
            cIdx[ii][jj] = idx
            cW[ii][jj] = S[0]*S[1]
            cGrad[ii][jj] = G * S[::-1] 