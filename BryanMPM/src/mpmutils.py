import numpy as np
from itertools import izip, count
from dateutil.relativedelta import relativedelta as reldelta
# Utils for moving data between particles and grid

def integrate( cIdx, cW, pp, gg ):
    # Integrate particle values to grid (p->g)
    for (ppi,idx,w) in izip(pp,cIdx,cW):
        gg[idx] += ppi * w[:,np.newaxis]
    return gg

def interpolate( cIdx, cW, pp, gg ):
    # Interpolate grid values to particles pp
    for (ii,idx,w) in izip(count(),cIdx,cW):
        pp[ii] = np.sum( gg[idx] * w[:,np.newaxis], 0 )
    return pp

def gradient( cIdx, cGrad, pp, gg ):
    # Interpolate a gradient
    for (ii,idxi,gradi) in izip(count(),cIdx,cGrad):
        pp[ii] = [0]
        for (idx,grad) in izip( idxi, gradi ):
            gR = np.reshape( gg[idx], (2,1) )
            cg = np.reshape( grad, (1,2) )
            pp[ii] += np.dot( gR, cg )
    return pp

def divergence( cIdx, cGrad, pp, gg ):
    # Send divergence of particle field to the grid
    for (ppi,idxi,gradi) in izip(pp,cIdx,cW):
        for idx,grad in izip( idxi, gradi ):
            cg = np.reshape( grad, (2,1) )            
            gg[idx] -= np.reshape( np.dot( ppi, cg ), 2 )
    return gg   

def gradscalar( cIdx, cGrad, pp, gg ):
    # Send gradient of particle field to the grid
    for (ppi,idxi,gradi) in izip(pp,cIdx,cW):
        for idx,grad in izip( idxi, gradi ):
            gg[idx] += np.reshape( ppi * grad, 2 )
    return gg   


def dotAdd( pp, qq ):
    # return pp += qq dot pp
    for (ii,ppi,qqi) in izip(count(),pp,qq):
        pp[ii] += np.dot( qqi, ppi )


def readableTime( tt ):
    attrs = ['years', 'months', 'days', 'hours', 'minutes', 'seconds']    
    h_read = lambda delta: ['%d %s' % (getattr(delta, attr), 
                   getattr(delta, attr) > 1 and attr or attr[:-1]) 
    for attr in attrs if getattr(delta, attr)]
    tm = h_read(reldelta(seconds=tt))
    
    st = ''
    for t in tm:
        st += t + ', '
    return st[:-2]