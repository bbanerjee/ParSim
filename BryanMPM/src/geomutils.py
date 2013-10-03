import numpy as np
import operator
import warnings

def ellipseLvl( r, c, n=2, x=0 ):
    pwr = np.power
    tol = 1.e-14
    n = 1.*n
    def lvl(x):
        return 1. - pwr( sum(((x-c)/r)**n), 1/n )
    def normal(x):
        warnings.filterwarnings("ignore")
        nml = 1./r * pwr( sum(((x-c)/r)**n), 1/n-1 ) * ((x-c)/r)**(n-1)
        nml = nml/np.linalg.norm(nml)
        if np.isnan(nml[0]): return x*0.
        else: return nml
    return lvl, normal

def rectLvl( x0, x1, x=0 ):
    c = (x0+x1)/2.
    r = (x1-x0)/2.
    return ellipseLvl( r, c, 200. )


def fillLvl( lvl, normal, patch ):
    nn = pCeil(patch.X1 - patch.X0) / (patch.dX/patch.ppe)
    ps = (patch.X1-patch.X0)/nn
    vol = patch.thick * ps[0] * ps[1]
    parts = []
    normals = []
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = patch.X0 + ps*ns
            if (patch.inPatch(pt) and lvl(pt) >= 0):
                parts.append(pt)
                normals.append( normal(pt) )

    parts = np.array(parts)
    normals = np.array( normals )
    return (parts, vol, normals)    


def pCeil( x ):
    tol = 1.e-14
    return np.ceil(x-tol)