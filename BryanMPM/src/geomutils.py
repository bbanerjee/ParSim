import numpy as np
import operator
import warnings

def ellipseLvl( r, c, n=2, x=0 ):
    pwr = np.power
    tol = 1.e-14
    n = 1.*n
    def lvl(x):
        return 1. - pwr( sum(((x-c)/r)**n), 1/n )
    return lvl

def rectLvl( x0, x1, x=0 ):
    c = (x0+x1)/2.
    r = (x1-x0)/2.
    return ellipseLvl( r, c, 200. )


def fillLvl( lvl, patch ):
    nn = pCeil((patch.X1 - patch.X0) / (patch.dX/patch.ppe))
    ps = (patch.X1-patch.X0)/nn
    vol = patch.thick * ps[0] * ps[1]
    parts = []
    for jj in range(int(nn[1])):
        for ii in range(int(nn[0])):
            ns = np.array([ii+0.5,jj+0.5])
            pt = patch.X0 + ps*ns
            if (patch.inPatch(pt) and lvl(pt) >= 0):
                parts.append(pt)

    parts = np.array(parts)
    return (parts, vol)    


def pCeil( x ):
    tol = 1.e-14
    return np.ceil(x-tol)