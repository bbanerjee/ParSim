import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy as copy
from matplotlib.widgets import Slider, Button, RadioButtons

def vnorm( x ):
    xx = np.copysign(np.sqrt((x*x).sum(axis=1)[:,np.newaxis]), x[:,0])
    return np.reshape( xx, xx.size )

def pPosition( data, t ):
    px = copy( data[t].get('px',1) )
    px = np.append(px,data[t].get('px',2), axis=0 )
    return px

def pVelocity( data, t ):
    pv = copy( data[t].get('pxI',1) )
    pv = np.append(pv,data[t].get('pxI',2), axis=0 )
    return vnorm(pv)

def pVonMises( data, t ):
    pv = copy( data[t].get('pxI',1) )
    pv = np.append(pv,data[t].get('pxI',2), axis=0 )
    return vnorm(pv)


def plotInteractive( data ):
    ax = plt.subplot(111)
    plt.subplots_adjust(bottom=0.25)
    t = sorted(data.keys())
    
    pl = plt.scatter([0],[0],s=10,linewidths=(0,0,0))
    plt.axis([0, 2, 0, 2])        
    
    axTime = plt.axes([0.25, 0.1, 0.65, 0.03]) 
    sTime = Slider(axTime, 'Time', min(t), max(t), valinit=min(t))
    
    def updateTime(val):
        t0 = sTime.val
        dt = abs(np.array(t)-t0)
        t_in = t[dt.argmin()]
        px = getpx(data,t_in)
        pv = getpv(data,t_in)
        pl.set_offsets(px)  
        pl.set_array(pv)
        plt.draw()
        
    updateTime(min(t))
    sTime.on_changed(updateTime)
        
    plt.show(block=False)