import numpy as np
from dateutil.relativedelta import relativedelta as reldelta
# Utils for moving data between particles and grid

def integrate( cIdx, cW, pp, gg, pIdx ):
    # Integrate particle values to grid (p->g)
    for ii in pIdx:
        gg[cIdx[ii]] += pp[ii] * cW[ii][:,np.newaxis]
    return gg

def interpolate( cIdx, cW, pp, gg, pIdx ):
    # Interpolate grid values to particles pp
    for ii in pIdx:
        pp[ii] = np.sum( gg[cIdx[ii]] * cW[ii][:,np.newaxis], 0 )
    return pp

def gradient( cIdx, cGrad, pp, gg, pIdx ):
    # Interpolate a gradient
    for ii in pIdx:
        pp[ii] = [0]
        for idx,grad in zip( cIdx[ii], cGrad[ii] ):
            gR = np.reshape( gg[idx], (2,1) )
            cg = np.reshape( grad, (1,2) )
            pp[ii] += np.dot( gR, cg )
    return pp

def divergence( cIdx, cGrad, pp, gg, pIdx ):
    # Send divergence of particle field to the grid
    for ii in pIdx:
        for idx,grad in zip( cIdx[ii], cGrad[ii] ):
            cg = np.reshape( grad, (2,1) )            
            gg[idx] -= np.reshape( np.dot( pp[ii], cg ), 2 )
    return gg   

def dotAdd( pp, qq, pIdx ):
    # return pp += qq dot pp
    for ii in pIdx:
        pp[ii] += np.dot( qq[ii], pp[ii] )


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