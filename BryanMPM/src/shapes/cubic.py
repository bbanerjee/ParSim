import numpy as np
    
def uSG( x, h ):
    r = abs(x)
    sgnx = cmp(x,-x)
    
    if ( r < 0.5*h ):
        S = -r*r/(h*h) + 3./4
        G = -2.*x/(h*h)
    elif ( r < 1.5*h ): 
        S = r*r/(2.*h*h) - 3.*r/(2.*h) + 9./8
        G = x/(h*h) - sgnx*3./(2*h)        
    else: 
        S = G = 0.
    return( S,G )
	

def getCell( patch, pos ):    
    # Gets lower left node of 1-cell block
    x_sc = (pos - patch.X0)/patch.dX + patch.nGhost
    idx = np.floor(x_sc)
    ii = idx[0]-1 
    jj = idx[1]-1 
	
    return int(jj * patch.Nc[0] + ii)
    

def updateContribList( dw, patch, dwi ):
    # Update node contribution list
    nx = patch.Nc[0]
    h = patch.dX
    idxs = [0,1,2,3,nx,nx+1,nx+2,nx+3,2*nx,2*nx+1,2*nx+2,2*nx+3]
    idxs = idxs + [3*nx, 3*nx+1, 3*nx+2, 3*nx+3]
    S = np.zeros(h.size)
    G = np.zeros(h.size)	
    
    cIdx,cW,cGrad = dw.getMult( ['cIdx','cW','cGrad'], dwi )
    px,gx = dw.getMult( ['px','gx'], dwi )

    for ii in range(len(pVol)):
        cc = getCell( patch, px[ii] )	           

        for jj in range(9):	
            idx = idxs[jj] + cc 
            r = px[ii] - gx[idx]	
		
            for kk in range(len(r)):
                S[kk],G[kk] = uSG( r[kk], h[kk] )
		
            cIdx[ii][jj] = idx
            cW[ii][jj] = S[0]*S[1]
            cGrad[ii][jj] = G * S[::-1] 