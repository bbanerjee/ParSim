import numpy as np

def fillRectangle( pt1, pt2, ppe, patch, dw, matid, density ):
    nn = pCeil( (pt2-pt1) / (patch.dX/ppe) )
    ps = (pt2-pt1)/nn
    vol = patch.thick * ps[0] * ps[1]
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = pt1 + ps*ns
            if patch.inPatch( pt ):
                dw.addParticle( matid, pt, vol*density, vol )


def fillAnnulus( pt1, r, ppe, patch, dw, matid, density ):
    nn = pCeil( 2*r[1] / (patch.dX/ppe) )
    ps = 2.0*r[1]/nn
    vol = patch.thick * ps[0] * ps[1]
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = pt1 - r[1] + ps*ns
            if patch.inPatch( pt ):
                if ( r[0] <= np.linalg.norm( pt1 - pt ) <= r[1] ):
                    dw.addParticle( matid, pt, vol*density, vol )    


def pCeil( x ):
    tol = 1.e-14
    return np.ceil(x-tol)